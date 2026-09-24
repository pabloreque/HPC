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

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.txt> [output.h5] [group]\n"
                  << "       " << argv[0] << " <input.txt> [output.h5] --group <group>\n";
        return 1;
    }
    const std::string inputPath = argv[1];
    std::string outputPath = "results.h5";
    std::string groupName = "/trajectory";
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--group") {
            if (++i >= argc) {
                std::cerr << "error: --group needs a value\n";
                return 1;
            }
            groupName = argv[i];
        } else if (arg.rfind("--", 0) == 0) {
            std::cerr << "error: unknown option '" << arg << "'\n";
            return 1;
        } else if (outputPath == "results.h5" && (i == 2)) {
            // First positional after input is the output file... unless the
            // caller passed only a group; disambiguate by extension.
            if (arg.size() >= 3 && arg.substr(arg.size() - 3) == ".h5") {
                outputPath = arg;
            } else {
                groupName = arg;
            }
        } else if (groupName == "/trajectory") {
            groupName = arg;
        } else {
            std::cerr << "error: unexpected argument '" << arg << "'\n";
            return 1;
        }
    }
    // Edge case: "./prog input.txt results.h5" leaves group default -- correct.
    // "./prog input.txt /cyclotron" is detected above by missing .h5 suffix.

    Input input;
    try {
        input = parseInputFile(inputPath);
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    try {
        H5::Exception::dontPrint();
        H5::H5File file(outputPath,
                        fileExists(outputPath) ? H5F_ACC_RDWR : H5F_ACC_TRUNC);
        replaceGroup(file, groupName);
        H5::Group group = file.createGroup(groupName);

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
            std::cout << "particle " << p << " -> " << groupName << ":";
            printSummary(trajectory);
        }
        std::cout << "Wrote " << input.particles.size() << " trajector"
                  << (input.particles.size() == 1 ? "y" : "ies") << " from '" << inputPath << "' to '" << outputPath
                  << groupName << "'\n";
    } catch (const H5::Exception &e) {
        std::cerr << "HDF5 error while writing '" << outputPath << groupName << "'\n";
        return 1;
    }
    return 0;
}
