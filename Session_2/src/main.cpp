// Runs the Boris pusher on two scenarios in normalized units (q = m = 1,
// |B| = 1) and exports both trajectories to results.h5 for analysis.ipynb:
//   - cyclotron: pure magnetic field -> circular gyration
//   - exb_drift: uniform E perpendicular to B -> gyration plus E x B drift
// See README.md for the physics and how the results are validated.

#include <H5Cpp.h>
#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
#include "scenarios.hpp"

std::array<std::vector<double>, 3> splitComponents(const std::vector<Vector3> &vectors) {
    std::array<std::vector<double>, 3> components;
    components[0].reserve(vectors.size());
    components[1].reserve(vectors.size());
    components[2].reserve(vectors.size());
    for (const Vector3 &v : vectors) {
        components[0].push_back(v.x);
        components[1].push_back(v.y);
        components[2].push_back(v.z);
    }
    return components;
}

void writeDataset(H5::Group &group, const std::string &name, const std::vector<double> &data,
                   const H5::DataSpace &dataspace) {
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

void writeTrajectory(H5::H5File &file, const std::string &groupName, const Trajectory &trajectory) {
    H5::Group group = file.createGroup(groupName);
    hsize_t dims[1] = {static_cast<hsize_t>(trajectory.time.size())};
    H5::DataSpace dataspace(1, dims);

    writeDataset(group, "t", trajectory.time, dataspace);

    const auto position = splitComponents(trajectory.position);
    writeDataset(group, "x", position[0], dataspace);
    writeDataset(group, "y", position[1], dataspace);
    writeDataset(group, "z", position[2], dataspace);

    const auto velocity = splitComponents(trajectory.velocity);
    writeDataset(group, "vx", velocity[0], dataspace);
    writeDataset(group, "vy", velocity[1], dataspace);
    writeDataset(group, "vz", velocity[2], dataspace);
}

Trajectory runScenario(const Scenario &scenario) {
    return simulateBorisPusher(scenario.initialState, scenario.charge, scenario.mass, scenario.E, scenario.B,
                                scenario.dt, scenario.numSteps);
}

void printSpeedConservation(const Trajectory &trajectory) {
    const double speed0 = norm(trajectory.velocity.front());
    const double speedN = norm(trajectory.velocity.back());
    std::printf("  |v0| = %.12f, |v_final| = %.12f, drift = %.3e\n", speed0, speedN, std::fabs(speedN - speed0));
}

void printDriftEstimate(const Trajectory &trajectory, const Vector3 &analyticDrift) {
    double meanVy = 0.0;
    for (const Vector3 &v : trajectory.velocity) {
        meanVy += v.y;
    }
    meanVy /= static_cast<double>(trajectory.velocity.size());
    std::printf("  mean(vy) = %.6f, analytic v_drift.y = %.6f\n", meanVy, analyticDrift.y);
}

int main() {
    const std::string outFile = "results.h5";
    H5::H5File file(outFile, H5F_ACC_TRUNC);

    const Scenario cyclotron = cyclotronScenario();
    const Trajectory cyclotronTrajectory = runScenario(cyclotron);
    writeTrajectory(file, "/cyclotron", cyclotronTrajectory);

    const Scenario drift = exbDriftScenario();
    const Trajectory driftTrajectory = runScenario(drift);
    writeTrajectory(file, "/drift", driftTrajectory);

    std::cout << "cyclotron scenario (speed conservation check):\n";
    printSpeedConservation(cyclotronTrajectory);

    std::cout << "exb_drift scenario (drift velocity check):\n";
    printDriftEstimate(driftTrajectory, analyticDriftVelocity(drift.E, drift.B));

    std::cout << "\nResults exported to " << outFile << "\n";

    return 0;
}
