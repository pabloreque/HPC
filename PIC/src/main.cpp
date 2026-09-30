/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
Generic driver: reads .txt initial conditions, writes HDF5 trajectory.

Reads arbitrary initial conditions from a .txt input file (see examples/),
runs the pusher, and appends the trajectory to an HDF5 file under the
requested group. Assumes uniform static E/B shared by all particles.

Usage:
  session_2_prequeijo <input.txt> [output.h5] [group]
  session_2_prequeijo <input.txt> [output.h5] --group <group>
Defaults: output.h5 = "results.h5", group = "/trajectory".
If output.h5 exists it is opened read-write and the group is replaced;
otherwise it is created. To regenerate both notebook datasets:
  ./session_2_prequeijo ../examples/cyclotron.txt results.h5 /cyclotron
  ./session_2_prequeijo ../examples/drift.txt results.h5 /drift
*/

#include <H5Cpp.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <iostream>
#include <string>

#include "boris_pusher.hpp"
#include "cli.hpp"
#include "hdf5_writer.hpp"
#include "input.hpp"
#include "vector3.hpp"

// Prints |v0|, |v_final| and mean velocity of one trajectory.
void printSummary(const Trajectory &trajectory) {
    const double speed0 = norm(trajectory.velocity.front());
    const double speedN = norm(trajectory.velocity.back());
    double meanVx = 0.0, meanVy = 0.0, meanVz = 0.0;
    for (const Vector3 &v : trajectory.velocity) {
        meanVx += v.x;
        meanVy += v.y;
        meanVz += v.z;
    }
    const double n = static_cast<double>(trajectory.velocity.size());
    std::printf("  steps=%zu |v0|=%.12f |v_final|=%.12f drift=%.3e mean(v)=(%.6f, %.6f, %.6f)\n",
                trajectory.velocity.size(), speed0, speedN, std::fabs(speedN - speed0), meanVx / n, meanVy / n,
                meanVz / n);
}

// Runs every particle in input, writing each trajectory under groupName
// (flat for a single particle, particle_N subgroups otherwise).
void runSimulation(const Input &input, const RunConfig &config) {
    H5::Exception::dontPrint();
    H5::H5File file(config.outputPath,
                    fileExists(config.outputPath) ? H5F_ACC_RDWR : H5F_ACC_TRUNC);
    replaceGroup(file, config.groupName);
    H5::Group group = file.createGroup(config.groupName);

    const EMFields fields{input.E, input.B};
    for (std::size_t p = 0; p < input.particles.size(); ++p) {
        const Trajectory trajectory =
            simulateBorisPusher(input.particles[p], fields, input.dt, input.numSteps);
        if (input.particles.size() == 1) {
            writeTrajectory(group, trajectory);
        } else {
            H5::Group sub = group.createGroup("particle_" + std::to_string(p));
            writeTrajectory(sub, trajectory);
        }
        std::cout << "particle " << p << " -> " << config.groupName << ":";
        printSummary(trajectory);
    }
    std::cout << "Wrote " << input.particles.size() << " trajector"
              << (input.particles.size() == 1 ? "y" : "ies") << " from '" << config.inputPath << "' to '"
              << config.outputPath << config.groupName << "'\n";
}

int main(int argc, char *argv[]) {
    RunConfig config;
    try {
        config = parseArgs(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    Input input;
    try {
        input = parseInputFile(config.inputPath);
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    try {
        runSimulation(input, config);
    } catch (const H5::Exception &e) {
        std::cerr << "HDF5 error while writing '" << config.outputPath << config.groupName << "'\n";
        return 1;
    }
    return 0;
}
