/*
Pablo Requeijo, November 26 2026.
HPC coursework: PIC.
Python module pic_prequeijo (pybind11): a second front end to the same
simulation code as the C++ driver (src/simulation.hpp), returning NumPy
arrays named like the driver's HDF5 datasets.

  import pic_prequeijo as pic
  sim = pic.HybridSimulation.from_file("examples/alfven_R.txt")
  sim.run()                        # to the file's steps= (or sim.run(100))
  history = sim.history()          # t, x_primal, x_dual, Ex, ..., Ti
  sim.save("build/results.h5", "/alfven_R")

  trajectories = pic.run_test_particles("examples/cyclotron.txt")

Every array returned is a copy: modifying it does not change the run.
*/
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <H5Cpp.h>

#include <array>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
#include "grid.hpp"
#include "hdf5_writer.hpp"
#include "hybrid_settings.hpp"
#include "hybrid_step.hpp"
#include "input_file.hpp"
#include "loader.hpp"
#include "parse.hpp"
#include "simulation.hpp"
#include "snapshot.hpp"
#include "trajectory.hpp"
#include "vector3.hpp"

namespace py = pybind11;

// A NumPy argument seen as contiguous doubles; other dtypes and Python
// sequences are converted (copied only when needed).
using DoubleArray =
    py::array_t<double, py::array::c_style | py::array::forcecast>;

// ---------------------------------------------------------------------------
// C++ <-> NumPy conversions
// ---------------------------------------------------------------------------

// Copies a vector of doubles into a 1D NumPy array.
py::array_t<double> toArray(const std::vector<double> &values) {
    return py::array_t<double>(static_cast<py::ssize_t>(values.size()),
                               values.data());
}

// Copies row-major doubles into a (rows x columns) NumPy array.
py::array_t<double> toArray(const std::vector<double> &values, std::size_t rows,
                            std::size_t columns) {
    const std::vector<py::ssize_t> shape = {static_cast<py::ssize_t>(rows),
                                            static_cast<py::ssize_t>(columns)};
    return py::array_t<double>(shape, values.data());
}

// Copies 3D vectors into an (N x 3) NumPy array.
py::array_t<double> toArray(const std::vector<Vector3> &vectors) {
    py::array_t<double> array(
        std::vector<py::ssize_t>{static_cast<py::ssize_t>(vectors.size()), 3});
    auto out = array.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < out.shape(0); ++i) {
        const Vector3 &vector = vectors[static_cast<std::size_t>(i)];
        out(i, 0) = vector.x;
        out(i, 1) = vector.y;
        out(i, 2) = vector.z;
    }
    return array;
}

// Copies a 1D array argument 'name' into a vector of doubles.
std::vector<double> toVector(const py::handle &value, const std::string &name) {
    const DoubleArray array = DoubleArray::ensure(value);
    require(array && array.ndim() == 1,
            "'" + name + "' must be a 1D array of numbers");
    return std::vector<double>(array.data(), array.data() + array.size());
}

// Copies an (N x 3) array argument 'name' into 3D vectors.
std::vector<Vector3> toVector3List(const py::handle &value,
                                   const std::string &name) {
    const DoubleArray array = DoubleArray::ensure(value);
    require(array && array.ndim() == 2 && array.shape(1) == 3,
            "'" + name + "' must be an (N, 3) array of numbers");
    const auto in = array.unchecked<2>();
    std::vector<Vector3> vectors;
    for (py::ssize_t i = 0; i < in.shape(0); ++i) {
        vectors.push_back(Vector3{in(i, 0), in(i, 1), in(i, 2)});
    }
    return vectors;
}

// Reads a 3-component argument 'name' (sequence or array) into a Vector3.
Vector3 toVector3(const py::handle &value, const std::string &name) {
    const std::vector<double> components = toVector(value, name);
    require(components.size() == 3, "'" + name + "' must have 3 components");
    return Vector3{components[0], components[1], components[2]};
}

// True if 'value' is a single number rather than an array.
bool isScalar(const py::handle &value) {
    const DoubleArray array = DoubleArray::ensure(value);
    return array && array.ndim() == 0;
}

// ---------------------------------------------------------------------------
// Hybrid input from keyword arguments
// ---------------------------------------------------------------------------

// Applies one keyword argument to the Input being built.
using KeySetter = std::function<void(Input &, const py::handle &)>;

// Keyword arguments without a default, as in the .txt hybrid mode (steps
// may be left out: run() takes a step count).
constexpr std::array<const char *, 5> kRequiredHybridKeys = {
    "dt", "B", "length", "nx", "ppc"};

// A scalar moment given as a number (stored in 'uniform') or as one value
// per primal node (stored in 'profile'), like density= and vth=.
void setScalarMoment(double &uniform, std::vector<double> &profile,
                     const py::handle &value, const std::string &name) {
    if (isScalar(value)) {
        uniform = value.cast<double>();
    } else {
        profile = toVector(value, name);
    }
}

// density=: a number (uniform) or one value per primal node (profile).
void setDensity(Input &input, const py::handle &value) {
    PopulationSettings &population = input.hybrid->population;
    setScalarMoment(population.density, population.profiles.density, value,
                    "density");
}

// vth=: a number (uniform) or one value per primal node (profile).
void setThermalSpeed(Input &input, const py::handle &value) {
    PopulationSettings &population = input.hybrid->population;
    setScalarMoment(population.vth, population.profiles.thermalSpeed, value,
                    "vth");
}

// bulk_velocity=: an (nx, 3) array, one ion bulk velocity per primal node.
void setBulkVelocity(Input &input, const py::handle &value) {
    input.hybrid->population.profiles.bulkVelocity =
        toVector3List(value, "bulk_velocity");
}

// seed=: a non-negative integer, checked like the .txt key.
void setSeed(Input &input, const py::handle &value) {
    input.hybrid->population.seed =
        parseSeed(std::to_string(value.cast<int>()), "seed");
}

// The setters of the hybrid keyword arguments, named like the .txt keys
// (bulk_velocity replaces the ux, uy, uz profile columns).
const std::map<std::string, KeySetter> &hybridKeySetters() {
    using H = const py::handle &;
    static const std::map<std::string, KeySetter> setters = {
        {"dt", [](Input &in, H v) { in.dt = v.cast<double>(); }},
        {"steps", [](Input &in, H v) { in.numSteps = v.cast<int>(); }},
        {"B", [](Input &in, H v) { in.B = toVector3(v, "B"); }},
        {"length", [](Input &in, H v) { in.length = v.cast<double>(); }},
        {"nx", [](Input &in, H v) { in.hybrid->nx = v.cast<int>(); }},
        {"iterations",
         [](Input &in, H v) { in.hybrid->iterations = v.cast<int>(); }},
        {"output_every",
         [](Input &in, H v) { in.hybrid->outputEvery = v.cast<int>(); }},
        {"Te",
         [](Input &in, H v) {
             in.hybrid->electrons.temperature = v.cast<double>();
         }},
        {"eta",
         [](Input &in, H v) {
             in.hybrid->electrons.resistivity = v.cast<double>();
         }},
        {"nu",
         [](Input &in, H v) {
             in.hybrid->electrons.hyperResistivity = v.cast<double>();
         }},
        {"ppc",
         [](Input &in, H v) { in.hybrid->population.ppc = v.cast<int>(); }},
        {"density", setDensity},
        {"vth", setThermalSpeed},
        {"bulk_velocity", setBulkVelocity},
        {"seed", setSeed},
        {"wave_mode",
         [](Input &in, H v) { in.hybrid->wave.mode = v.cast<int>(); }},
        {"wave_amplitude",
         [](Input &in, H v) { in.hybrid->wave.amplitude = v.cast<double>(); }},
        {"wave_branch", [](Input &in, H v) {
             in.hybrid->wave.branch = parseWaveBranch(v.cast<std::string>());
         }}};
    return setters;
}

// Applies one keyword argument; an unknown name raises TypeError, like a
// wrong keyword in a Python call.
void applyKeyword(Input &input, const std::string &key,
                  const py::handle &value) {
    const auto &setters = hybridKeySetters();
    const auto setter = setters.find(key);
    if (setter == setters.end()) {
        throw py::type_error("unexpected keyword argument '" + key + "'");
    }
    setter->second(input, value);
}

// Builds a hybrid Input from keyword arguments named like the .txt keys.
// Values are range-checked later, by startHybridRun.
Input hybridInputFromKwargs(const py::kwargs &kwargs) {
    Input input;
    input.hybrid = HybridSettings{};
    for (const auto &[key, value] : kwargs) {
        applyKeyword(input, key.cast<std::string>(), value);
    }
    for (const char *key : kRequiredHybridKeys) {
        require(kwargs.contains(key),
                std::string("missing required key '") + key + "'");
    }
    return input;
}

// Parses a .txt input and checks it is in the mode the caller handles.
Input readInputFile(const std::filesystem::path &path, bool hybrid) {
    Input input = parseInputFile(path.string());
    require(input.hybrid.has_value() == hybrid,
            path.string() +
                (hybrid ? " is a test-particle input (no nx=): use "
                          "run_test_particles"
                        : " is a hybrid input (nx=): use HybridSimulation"));
    return input;
}

// ---------------------------------------------------------------------------
// Hybrid simulation
// ---------------------------------------------------------------------------

// Adds the primal and dual node positions x_primal, x_dual to 'arrays'.
void addNodePositions(py::dict &arrays, const Grid1D &grid) {
    arrays["x_primal"] = toArray(nodePositions(grid, Centering::Primal));
    arrays["x_dual"] = toArray(nodePositions(grid, Centering::Dual));
}

// The snapshots of 'run' as NumPy arrays, named like the HDF5 datasets:
// t, x_primal, x_dual and one (snapshots x nx) array per saved quantity.
py::dict historyArrays(const HybridRun &run) {
    py::dict arrays;
    arrays["t"] = toArray(snapshotTimes(run.history));
    addNodePositions(arrays, run.state.grid);
    const auto nx = static_cast<std::size_t>(run.state.grid.nx);
    for (std::size_t q = 0; q < kSavedFieldCount; ++q) {
        arrays[kSavedFieldNames[q]] =
            toArray(stackSnapshots(run.history, q), run.history.size(), nx);
    }
    return arrays;
}

// The saved grid quantities at the current time, one 1D array each, plus
// t (a number), x_primal and x_dual.
py::dict currentFieldArrays(const HybridRun &run) {
    const FieldSnapshot snapshot = takeSnapshot(run.state, currentTime(run));
    const auto fields = snapshotFields(snapshot);
    py::dict arrays;
    arrays["t"] = snapshot.time;
    addNodePositions(arrays, run.state.grid);
    for (std::size_t q = 0; q < kSavedFieldCount; ++q) {
        arrays[kSavedFieldNames[q]] = toArray(fields[q]->values);
    }
    return arrays;
}

// The macroparticles at the current time: x (N), v (N x 3), weight (N).
py::dict particleArrays(const HybridState &state) {
    std::vector<double> position, weight;
    std::vector<Vector3> velocity;
    for (const Particle1D &particle : state.particles) {
        position.push_back(particle.position);
        velocity.push_back(particle.velocity);
        weight.push_back(particle.weight);
    }
    py::dict arrays;
    arrays["x"] = toArray(position);
    arrays["v"] = toArray(velocity);
    arrays["weight"] = toArray(weight);
    return arrays;
}

// Magnetic, kinetic and total energy of the current state, as a dict.
py::dict energyDict(const HybridState &state) {
    const Energies w = energies(state);
    py::dict values;
    values["magnetic"] = w.magnetic;
    values["kinetic"] = w.kinetic;
    values["total"] = w.total;
    return values;
}

// Writes the history and settings of 'run' to group 'groupName' of the
// HDF5 file 'path', exactly as the C++ driver does.
void saveRun(const HybridRun &run, const std::filesystem::path &path,
             const std::string &groupName) {
    H5::H5File file = openOutputFile(path.string());
    H5::Group group = recreateGroup(file, groupName);
    Input settings = run.input;
    settings.numSteps = run.step;
    writeRunAttributes(group, settings);
    writeFieldHistory(group, run.history, run.state.grid);
}

// Emits each stability warning of 'run' as a Python RuntimeWarning.
void warnStability(const HybridRun &run) {
    for (const std::string &warning :
         stabilityWarnings(run.input, run.state.grid)) {
        if (PyErr_WarnEx(PyExc_RuntimeWarning, warning.c_str(), 1) != 0) {
            throw py::error_already_set();
        }
    }
}

// A hybrid run owned by Python, advanced on demand.
class HybridSimulation {
  public:
    // Starts the run of 'input' and emits its stability warnings.
    explicit HybridSimulation(const Input &input)
        : run_(startHybridRun(input)) {
        warnStability(run_);
    }

    // Advances 'steps' steps; by default, up to the input's steps=.
    void run(std::optional<int> steps) {
        const int count = steps.value_or(run_.input.numSteps - run_.step);
        require(count >= 0, "steps must be >= 0");
        // The loop touches no Python object: let other threads run.
        py::gil_scoped_release release;
        advanceHybridRun(run_, count);
    }

    // Read access to the run, for the module functions above.
    const HybridRun &state() const {
        return run_;
    }

    // Short description: grid size, particle count, step and time.
    std::string repr() const {
        return "HybridSimulation(nx=" + std::to_string(run_.state.grid.nx) +
               ", particles=" + std::to_string(run_.state.particles.size()) +
               ", step=" + std::to_string(run_.step) +
               ", t=" + std::to_string(currentTime(run_)) + ")";
    }

  private:
    HybridRun run_;
};

// ---------------------------------------------------------------------------
// Test particles
// ---------------------------------------------------------------------------

// One particle's trajectory as NumPy arrays t, x, y, z, vx, vy, vz (the
// HDF5 dataset names) plus its q and m.
py::dict trajectoryArrays(const Trajectory &trajectory,
                          const Species &species) {
    py::dict arrays;
    for (const TrajectoryColumn &column : trajectoryColumns(trajectory)) {
        arrays[column.name] = toArray(column.values);
    }
    arrays["q"] = species.charge;
    arrays["m"] = species.mass;
    return arrays;
}

// Runs every particle of a test-particle input; one dict per particle.
py::list runTestParticles(const Input &input) {
    validateTestParticles(input);
    std::vector<Trajectory> trajectories;
    {
        py::gil_scoped_release release;
        trajectories = simulateTestParticles(input);
    }
    py::list result;
    for (std::size_t p = 0; p < trajectories.size(); ++p) {
        result.append(
            trajectoryArrays(trajectories[p], input.particles[p].species));
    }
    return result;
}

// Builds test particles from per-particle arrays: q, m (N), x, v (N x 3).
std::vector<Particle3D> particlesFromArrays(const py::handle &q,
                                            const py::handle &m,
                                            const py::handle &x,
                                            const py::handle &v) {
    const std::vector<double> charge = toVector(q, "q");
    const std::vector<double> mass = toVector(m, "m");
    const std::vector<Vector3> position = toVector3List(x, "x");
    const std::vector<Vector3> velocity = toVector3List(v, "v");
    const std::size_t count = charge.size();
    require(mass.size() == count && position.size() == count &&
                velocity.size() == count,
            "q, m, x and v must describe the same number of particles");
    std::vector<Particle3D> particles;
    for (std::size_t p = 0; p < count; ++p) {
        particles.push_back(
            Particle3D{position[p], velocity[p], Species{charge[p], mass[p]}});
    }
    return particles;
}

// ---------------------------------------------------------------------------
// Analytic references
// ---------------------------------------------------------------------------

// Applies 'frequency' to every wavenumber of k: a number for a scalar k,
// an array of the same shape otherwise.
py::object mapWavenumbers(const DoubleArray &k,
                          const std::function<double(double)> &frequency) {
    py::array_t<double> omega(
        std::vector<py::ssize_t>(k.shape(), k.shape() + k.ndim()));
    for (py::ssize_t i = 0; i < k.size(); ++i) {
        omega.mutable_data()[i] = frequency(k.data()[i]);
    }
    if (k.ndim() == 0) {
        return py::float_(*omega.data());
    }
    return std::move(omega);
}

// Cold-plasma eigenfrequency of branch "R" or "L" at wavenumber k: a
// number for a scalar k, an array of the same shape otherwise.
py::object waveFrequencies(const DoubleArray &k, const std::string &name) {
    const WaveBranch branch = parseWaveBranch(name);
    return mapWavenumbers(
        k, [branch](double kk) { return waveFrequency(kk, branch); });
}

// Eigenfrequency of branch "R" or "L" on the discrete scheme (Yee curl,
// linear shape, 'ppc' evenly spaced cold ions per cell of size dx), with
// the shape of k.
py::object discreteWaveFrequencies(const DoubleArray &k, double dx, int ppc,
                                   const std::string &name) {
    const WaveBranch branch = parseWaveBranch(name);
    require(dx > 0.0, "dx must be > 0");
    require(ppc > 0, "ppc must be > 0");
    return mapWavenumbers(k, [=](double kk) {
        return discreteWaveFrequency(kk, dx, ppc, branch);
    });
}

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------

constexpr const char *kHybridSimulationDoc = R"doc(
Self-consistent 1D3V hybrid PIC run (kinetic ions, massless fluid
electrons, periodic Yee mesh), advanced on demand.

Build it from a .txt input with ``HybridSimulation.from_file(path)``, or
from keyword arguments named like the .txt hybrid keys (README.md):

Parameters
----------
dt, length : float
    Time step and periodic box length (required).
B : sequence of 3 floats
    Initial uniform magnetic field B0 (required).
nx, ppc : int
    Number of cells and macroparticles per cell (required).
steps : int, optional
    Default step count of ``run()`` (default 0).
iterations, output_every : int, optional
    ICN corrections per step (default 2), snapshot period in steps
    (default 1).
Te, eta, nu : float, optional
    Electron temperature, resistivity, hyper-resistivity (default 0).
density, vth : float or array of nx floats, optional
    Uniform value (defaults 1 and 0) or one value per primal node.
bulk_velocity : array of shape (nx, 3), optional
    Ion bulk velocity per primal node (default 0).
seed : int, optional
    Seed of the thermal velocity draws (default 1).
wave_mode, wave_amplitude, wave_branch : optional
    R/L eigenmode: wavelengths in the box, amplitude, "R" or "L".

Raises
------
ValueError
    If a value is out of range, as for the .txt input.
TypeError
    For an unknown keyword argument.

Warns
-----
RuntimeWarning
    If dt exceeds the approximate stability limit, or iterations=1.
)doc";

constexpr const char *kRunTestParticlesFileDoc = R"doc(
Runs a test-particle .txt input (Boris pusher in uniform E, B).

Parameters
----------
path : str or os.PathLike
    Input file in test-particle mode (no nx= line).

Returns
-------
list of dict
    One dict per particle with arrays t, x, y, z, vx, vy, vz (the HDF5
    dataset names) and the numbers q, m.
)doc";

constexpr const char *kRunTestParticlesDoc = R"doc(
Runs test particles in uniform static fields (Boris pusher).

Parameters
----------
q, m : array of N floats
    Charge and mass (> 0) of each particle.
x, v : array of shape (N, 3)
    Initial position and velocity of each particle.
E, B : sequence of 3 floats
    Uniform electric and magnetic fields.
dt : float
    Time step (> 0).
steps : int
    Number of steps (>= 0).
length : float, optional
    Periodic box length along x; None (default) for an unbounded box.

Returns
-------
list of dict
    One dict per particle with arrays t, x, y, z, vx, vy, vz and the
    numbers q, m.
)doc";

constexpr const char *kWaveFrequencyDoc = R"doc(
Analytic eigenfrequency of the parallel R or L wave in a cold plasma
with massless electrons, for exp(i(kx - omega t)): the roots of
omega^2 + k^2 omega - k^2 = 0 (normalized units, B0 = n0 = 1).

Parameters
----------
k : float or array
    Wavenumber, k > 0, in 1/d_i.
branch : {"R", "L"}
    "L" gives omega > 0, "R" (whistler) gives omega < 0.

Returns
-------
float or ndarray
    omega in units of Omega_ci, with the shape of k.
)doc";

constexpr const char *kDiscreteWaveFrequencyDoc = R"doc(
Eigenfrequency of the parallel R or L wave on the discrete scheme of this
code: the roots of omega^2 + (K^2 - sigma) omega - K^2 = 0, with the Yee
wavenumber K = (2/dx) sin(k dx/2) and the linear-shape ion response
deficit sigma = 4 sin^2(k dx/2) (1/6 + 1/(12 ppc^2)) (src/loader.hpp).
It is the exact linear dispersion of the code's spatial discretization
for cold ions loaded by the driver (evenly spaced), and tends to
wave_frequency as dx -> 0.

Parameters
----------
k : float or array
    Wavenumber, k > 0, in 1/d_i.
dx : float
    Cell size, in d_i.
ppc : int
    Macroparticles per cell.
branch : {"R", "L"}
    "L" gives omega > 0, "R" (whistler) gives omega < 0.

Returns
-------
float or ndarray
    omega in units of Omega_ci, with the shape of k.
)doc";

// Binds HybridSimulation and its methods.
void bindHybridSimulation(py::module_ &module) {
    py::class_<HybridSimulation>(module, "HybridSimulation",
                                 kHybridSimulationDoc)
        .def(py::init([](const py::kwargs &kwargs) {
            return HybridSimulation(hybridInputFromKwargs(kwargs));
        }))
        .def_static(
            "from_file",
            [](const std::filesystem::path &path) {
                return HybridSimulation(readInputFile(path, true));
            },
            py::arg("path"),
            "Builds the run of a hybrid .txt input (nx= line).")
        .def("run", &HybridSimulation::run, py::arg("steps") = py::none(),
             "Advances `steps` steps (default: up to the input's steps=), "
             "taking a snapshot every output_every steps.")
        .def(
            "history",
            [](const HybridSimulation &sim) {
                return historyArrays(sim.state());
            },
            "Snapshots so far as a dict of arrays named like the HDF5 "
            "datasets: t, x_primal, x_dual, and Ex..Ez, Bx..Bz, n, ux..uz, "
            "Ti of shape (snapshots, nx).")
        .def(
            "fields",
            [](const HybridSimulation &sim) {
                return currentFieldArrays(sim.state());
            },
            "Grid quantities now: t, x_primal, x_dual and one (nx,) array "
            "per saved quantity.")
        .def(
            "particles",
            [](const HybridSimulation &sim) {
                return particleArrays(sim.state().state);
            },
            "Macroparticles now: x (N,), v (N, 3), weight (N,).")
        .def(
            "energies",
            [](const HybridSimulation &sim) {
                return energyDict(sim.state().state);
            },
            "Magnetic, kinetic (from the particles) and total energy now.")
        .def(
            "max_div_b",
            [](const HybridSimulation &sim) {
                return maxDivergenceB(sim.state().state);
            },
            "Largest |dBx/dx| on the mesh (0 up to round-off).")
        .def(
            "save",
            [](const HybridSimulation &sim, const std::filesystem::path &path,
               const std::string &group) { saveRun(sim.state(), path, group); },
            py::arg("path"), py::arg("group"),
            "Writes the snapshots and settings to an HDF5 group, like the "
            "C++ driver (other groups of the file are kept).")
        .def_property_readonly(
            "step",
            [](const HybridSimulation &sim) { return sim.state().step; },
            "Number of steps done.")
        .def_property_readonly(
            "time",
            [](const HybridSimulation &sim) {
                return currentTime(sim.state());
            },
            "Current time, step * dt.")
        .def_property_readonly(
            "timestep_limit",
            [](const HybridSimulation &sim) {
                const HybridRun &run = sim.state();
                return hybridTimestepLimit(run.input, run.state.grid);
            },
            "Approximate stability limit on dt (grid-scale whistler and "
            "hyper-resistivity).")
        .def("__repr__", &HybridSimulation::repr);
}

// Binds both forms of run_test_particles: from a file, from arrays.
void bindTestParticles(py::module_ &module) {
    module.def(
        "run_test_particles",
        [](const std::filesystem::path &path) {
            return runTestParticles(readInputFile(path, false));
        },
        py::arg("path"), kRunTestParticlesFileDoc);
    module.def(
        "run_test_particles",
        [](const py::handle &q, const py::handle &m, const py::handle &x,
           const py::handle &v, const py::handle &E, const py::handle &B,
           double dt, int steps, std::optional<double> length) {
            Input input;
            input.particles = particlesFromArrays(q, m, x, v);
            input.E = toVector3(E, "E");
            input.B = toVector3(B, "B");
            input.dt = dt;
            input.numSteps = steps;
            input.length = length;
            return runTestParticles(input);
        },
        py::kw_only(), py::arg("q"), py::arg("m"), py::arg("x"), py::arg("v"),
        py::arg("E"), py::arg("B"), py::arg("dt"), py::arg("steps"),
        py::arg("length") = py::none(), kRunTestParticlesDoc);
}

// Turns HDF5 errors into Python RuntimeError with the library's message.
void translateHdf5Errors(std::exception_ptr error) {
    try {
        if (error) {
            std::rethrow_exception(error);
        }
    } catch (const H5::Exception &hdf5Error) {
        PyErr_SetString(PyExc_RuntimeError,
                        ("HDF5 error: " + hdf5Error.getDetailMsg()).c_str());
    }
}

PYBIND11_MODULE(pic_prequeijo, module) {
    module.doc() = "1D3V hybrid PIC and Boris test particles (C++ core).";
    // HDF5 errors are reported by the translator, not printed by HDF5.
    H5::Exception::dontPrint();
    py::register_exception_translator(translateHdf5Errors);
    bindHybridSimulation(module);
    bindTestParticles(module);
    module.def("wave_frequency", &waveFrequencies, py::arg("k"),
               py::arg("branch"), kWaveFrequencyDoc);
    module.def("discrete_wave_frequency", &discreteWaveFrequencies,
               py::arg("k"), py::arg("dx"), py::arg("ppc"), py::arg("branch"),
               kDiscreteWaveFrequencyDoc);
}
