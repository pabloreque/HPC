/*
Pablo Requeijo, September 24 2026 (Session 2); October 1 2026 (Session 3).
HPC coursework: PIC.
Generic driver: reads .txt initial conditions, writes HDF5 results.

Reads arbitrary initial conditions from a .txt input file (see examples/)
and runs one of two modes (see input.hpp):
- test-particle mode: Boris pusher in uniform static E/B, optionally in a
  periodic box; writes one trajectory per particle;
- hybrid mode: self-consistent 1D hybrid PIC (hybrid_step.hpp); writes
  field snapshots (t, x_primal, x_dual, Ex..Ez, Bx..Bz, n, ux..uz).
Results go to an HDF5 file under the requested group.

Usage:
  pic_prequeijo <input.txt> [output.h5] [group]
  pic_prequeijo <input.txt> [output.h5] --group <group>
Defaults: output.h5 = "results.h5", group = "/trajectory".
If output.h5 exists it is opened read-write and the group is replaced;
otherwise it is created. Each notebook dataset is one example run into
its own group (full list in README.md), e.g.:
  ./pic_prequeijo ../examples/cyclotron.txt results.h5 /cyclotron
  ./pic_prequeijo ../examples/drift.txt results.h5 /drift
  ./pic_prequeijo ../examples/periodic_drift.txt results.h5 /periodic_drift
  ./pic_prequeijo ../examples/alfven_R.txt results.h5 /alfven_R
  ./pic_prequeijo ../examples/alfven_L.txt results.h5 /alfven_L
  ./pic_prequeijo ../examples/thermal_noise.txt results.h5 /thermal_noise
*/

#include <H5Cpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
#include "cli.hpp"
#include "field_solver.hpp"
#include "grid.hpp"
#include "hdf5_writer.hpp"
#include "hybrid_step.hpp"
#include "input.hpp"
#include "loader.hpp"
#include "vector3.hpp"

// Averages the velocity over all recorded samples.
Vector3 meanVelocity(const Trajectory &trajectory) {
    Vector3 sum;
    for (const Vector3 &velocity : trajectory.velocity) {
        sum = sum + velocity;
    }
    return (1.0 / static_cast<double>(trajectory.velocity.size())) * sum;
}

// Prints |v0|, |v_final| and mean velocity of one trajectory.
void printSummary(const Trajectory &trajectory) {
    const double speed0 = norm(trajectory.velocity.front());
    const double speedN = norm(trajectory.velocity.back());
    const Vector3 mean = meanVelocity(trajectory);
    std::printf("  steps=%zu |v0|=%.12f |v_final|=%.12f drift=%.3e "
                "mean(v)=(%.6f, %.6f, %.6f)\n",
                trajectory.velocity.size(), speed0, speedN,
                std::fabs(speedN - speed0), mean.x, mean.y, mean.z);
}

// Writes one trajectory: flat into 'group' for a single-particle run,
// into a particle_<index> subgroup otherwise.
void writeParticle(H5::Group &group, const Trajectory &trajectory,
                   std::size_t index, std::size_t particleCount) {
    if (particleCount == 1) {
        writeTrajectory(group, trajectory);
        return;
    }
    H5::Group sub = group.createGroup("particle_" + std::to_string(index));
    writeTrajectory(sub, trajectory);
}

// Runs every particle in input, writing each trajectory under 'group'
// (flat for a single particle, particle_N subgroups otherwise).
void runTestParticles(const Input &input, H5::Group &group,
                      const RunConfig &config) {
    const EMFields fields{input.E, input.B};
    const std::size_t count = input.particles.size();
    for (std::size_t p = 0; p < count; ++p) {
        const Trajectory trajectory = simulateBorisPusher(
            input.particles[p], fields, input.dt, input.numSteps, input.length);
        writeParticle(group, trajectory, p, count);
        std::cout << "particle " << p << " -> " << config.groupName << ":";
        printSummary(trajectory);
    }
    std::cout << "Wrote " << count << " trajector" << (count == 1 ? "y" : "ies")
              << " from '" << config.inputPath << "' to '" << config.outputPath
              << config.groupName << "'\n";
}

// Magnetic energy sum(B^2/2) dx over the grid.
double magneticEnergy(const HybridState &state) {
    double sum = 0.0;
    for (const Field *component : {&state.B.x, &state.B.y, &state.B.z}) {
        for (double b : component->values) {
            sum += 0.5 * b * b;
        }
    }
    return sum * state.grid.dx;
}

// Ion kinetic energy sum(w m v^2/2) dx over all macroparticles.
double kineticEnergy(const HybridState &state) {
    double sum = 0.0;
    for (const Particle1D &particle : state.particles) {
        sum += 0.5 * particle.weight * particle.species.mass *
               dot(particle.velocity, particle.velocity);
    }
    return sum * state.grid.dx;
}

// Largest |dBx/dx|: the discrete div B of the 1D Yee mesh.
double maxDivergenceB(const HybridState &state) {
    const Field divB = derivativeX(state.B.x, state.grid);
    double largest = 0.0;
    for (double value : divB.values) {
        largest = std::max(largest, std::fabs(value));
    }
    return largest;
}

// Prints magnetic, kinetic and total energy plus max |div B| of a state.
// W = W_B + W_K leaves out the electron thermal energy and the ohmic
// losses, so it is only expected to be conserved for Te = eta = nu = 0.
void printEnergies(const char *label, const HybridState &state) {
    const double wB = magneticEnergy(state);
    const double wK = kineticEnergy(state);
    std::printf("  %s: W_B=%.10f W_K=%.10f W=%.10f max|divB|=%.2e\n", label, wB,
                wK, wB + wK, maxDivergenceB(state));
}

// Warns on stderr if dt exceeds the approximate stability limit of the
// grid-scale whistler and hyper-resistive damping.
void warnIfUnstable(const Input &input, const Grid1D &grid) {
    const HybridSettings &settings = *input.hybrid;
    const double limit =
        timestepLimit(grid, norm(input.B), settings.population.density,
                      settings.electrons.hyperResistivity);
    if (input.dt > limit) {
        std::cerr << "warning: dt=" << input.dt << " exceeds the stability "
                  << "limit ~" << limit << " (grid whistler"
                  << (settings.electrons.hyperResistivity > 0.0
                          ? " + hyper-resistivity)\n"
                          : ", n dx^2 / 2|B0|)\n");
    }
}

// Warns on stderr if a single ICN correction is requested: it is first
// order and lets whistlers grow at every dt.
void warnIfSingleCorrection(const Input &input) {
    if (input.hybrid->iterations == 1) {
        std::cerr << "warning: iterations=1 is unstable for whistlers "
                  << "(|g|^2 = 1 + (omega dt)^4/4); use 2\n";
    }
}

// Builds the initial hybrid state from the input settings.
HybridState initialHybridState(const Input &input, const Grid1D &grid) {
    const HybridSettings &settings = *input.hybrid;
    return initializeState(
        grid, initialMagneticField(input.B, settings.wave, grid),
        loadUniformPlasma(settings.population, settings.wave, grid),
        settings.electrons);
}

// Advances 'state' input.numSteps hybrid steps, taking a snapshot at
// t = 0 and every output_every steps.
std::vector<FieldSnapshot> runHybridLoop(HybridState &state,
                                         const Input &input) {
    const HybridSettings &settings = *input.hybrid;
    const HybridParameters params{input.dt, settings.iterations,
                                  settings.electrons};
    std::vector<FieldSnapshot> history{takeSnapshot(state, 0.0)};
    for (int step = 1; step <= input.numSteps; ++step) {
        hybridStep(state, params);
        if (step % settings.outputEvery == 0) {
            history.push_back(takeSnapshot(state, step * input.dt));
        }
    }
    return history;
}

// Runs the self-consistent hybrid PIC simulation and writes the field
// history under 'group'.
void runHybrid(const Input &input, H5::Group &group, const RunConfig &config) {
    const Grid1D grid = makeGrid(input.hybrid->nx, *input.length);
    warnIfUnstable(input, grid);
    warnIfSingleCorrection(input);
    HybridState state = initialHybridState(input, grid);
    std::cout << "hybrid: nx=" << grid.nx
              << " particles=" << state.particles.size()
              << " steps=" << input.numSteps << "\n";
    printEnergies("initial", state);
    const std::vector<FieldSnapshot> history = runHybridLoop(state, input);
    printEnergies("final  ", state);
    writeFieldHistory(group, history, grid);
    std::cout << "Wrote " << history.size() << " snapshot"
              << (history.size() == 1 ? "" : "s") << " from '"
              << config.inputPath << "' to '" << config.outputPath
              << config.groupName << "'\n";
}

// Opens (or creates) the output file, replaces the target group and runs
// the mode selected by the input file.
void runSimulation(const Input &input, const RunConfig &config) {
    H5::Exception::dontPrint();
    // Reopen an existing .h5 read-write to keep its other groups;
    // otherwise create it from scratch.
    H5::H5File file(config.outputPath,
                    std::filesystem::exists(config.outputPath) ? H5F_ACC_RDWR
                                                               : H5F_ACC_TRUNC);
    replaceGroup(file, config.groupName);
    H5::Group group = file.createGroup(config.groupName);
    if (input.hybrid) {
        runHybrid(input, group, config);
    } else {
        runTestParticles(input, group, config);
    }
}

// Parses the CLI and input file, runs the simulation, and returns
// 1 with a message on stderr if any of these steps fails.
int main(int argc, char *argv[]) {
    RunConfig config;
    try {
        config = parseArgs(argc, argv);
    } catch (const std::exception &error) {
        std::cerr << error.what() << "\n";
        return 1;
    }

    Input input;
    try {
        input = parseInputFile(config.inputPath);
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }

    try {
        runSimulation(input, config);
    } catch (const H5::Exception &) {
        std::cerr << "HDF5 error while writing '" << config.outputPath
                  << config.groupName << "'\n";
        return 1;
    }
    return 0;
}
