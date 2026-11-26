/*
Pablo Requeijo, September 24 2026 (Session 2); October 1 2026 (Session 3).
HPC coursework: PIC.
Generic driver: reads .txt initial conditions, writes HDF5 results.

Reads arbitrary initial conditions from a .txt input file (see examples/)
and runs one of two modes (see input.hpp):
- test-particle mode: Boris pusher in uniform static E/B, optionally in a
  periodic box; writes one trajectory per particle;
- hybrid mode: self-consistent 1D hybrid PIC (hybrid_step.hpp); writes
  field snapshots (t, x_primal, x_dual, Ex..Ez, Bx..Bz, n, ux..uz, Ti).
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
  ./pic_prequeijo ../examples/density_sine.txt results.h5 /density_sine
*/

#include <H5Cpp.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
#include "cli.hpp"
#include "hdf5_writer.hpp"
#include "hybrid_step.hpp"
#include "input.hpp"
#include "simulation.hpp"
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

// Writes one trajectory and its particle's q, m: flat into 'group' for a
// single-particle run, into a particle_<index> subgroup otherwise.
void writeParticle(H5::Group &group, const Trajectory &trajectory,
                   const Species &species, std::size_t index,
                   std::size_t particleCount) {
    H5::Group target =
        particleCount == 1
            ? group
            : group.createGroup("particle_" + std::to_string(index));
    writeTrajectory(target, trajectory);
    writeSpeciesAttributes(target, species);
}

// Runs every particle in input, writing each trajectory under 'group'
// (flat for a single particle, particle_N subgroups otherwise).
void runTestParticles(const Input &input, H5::Group &group,
                      const RunConfig &config) {
    const std::vector<Trajectory> trajectories = simulateTestParticles(input);
    const std::size_t count = trajectories.size();
    for (std::size_t p = 0; p < count; ++p) {
        writeParticle(group, trajectories[p], input.particles[p].species, p,
                      count);
        std::cout << "particle " << p << " -> " << config.groupName << ":";
        printSummary(trajectories[p]);
    }
    std::cout << "Wrote " << count << " trajector" << (count == 1 ? "y" : "ies")
              << " from '" << config.inputPath << "' to '" << config.outputPath
              << config.groupName << "'\n";
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

// Prints every stability warning of a hybrid run on stderr.
void printStabilityWarnings(const HybridRun &run) {
    for (const std::string &warning :
         stabilityWarnings(run.input, run.state.grid)) {
        std::cerr << "warning: " << warning << "\n";
    }
}

// Runs the self-consistent hybrid PIC simulation and writes the field
// history under 'group'.
void runHybrid(const Input &input, H5::Group &group, const RunConfig &config) {
    HybridRun run = startHybridRun(input);
    printStabilityWarnings(run);
    std::cout << "hybrid: nx=" << run.state.grid.nx
              << " particles=" << run.state.particles.size()
              << " steps=" << input.numSteps << "\n";
    printEnergies("initial", run.state);
    advanceHybridRun(run, input.numSteps);
    printEnergies("final  ", run.state);
    writeFieldHistory(group, run.history, run.state.grid);
    std::cout << "Wrote " << run.history.size() << " snapshot"
              << (run.history.size() == 1 ? "" : "s") << " from '"
              << config.inputPath << "' to '" << config.outputPath
              << config.groupName << "'\n";
}

// Opens (or creates) the output file, replaces the target group and runs
// the mode selected by the input file.
void runSimulation(const Input &input, const RunConfig &config) {
    H5::Exception::dontPrint();
    H5::H5File file = openOutputFile(config.outputPath);
    H5::Group group = recreateGroup(file, config.groupName);
    writeRunAttributes(group, input);
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
