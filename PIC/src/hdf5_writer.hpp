/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
HDF5 output: trajectory datasets and group handling.
*/
#pragma once

#include <H5Cpp.h>

#include <array>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
#include "vector3.hpp"

// Splits vector<Vector3> into one double array per axis (x, y, z).
inline std::array<std::vector<double>, 3>
splitComponents(const std::vector<Vector3> &vectors) {
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

// Writes one 1D array of doubles to disk inside 'group'.
inline void writeDataset(H5::Group &group, const std::string &name,
                         const std::vector<double> &data,
                         const H5::DataSpace &dataspace) {
    H5::DataSet dataset =
        group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

// Writes a full trajectory as 7 datasets (t, x, y, z, vx, vy, vz), all
// sharing one shape: one value per time step.
inline void writeTrajectory(H5::Group &group, const Trajectory &trajectory) {
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

// Deletes 'groupName' if it already exists, so re-running the driver
// replaces the old trajectory instead of failing.
inline void replaceGroup(H5::H5File &file, const std::string &groupName) {
    if (H5Lexists(file.getId(), groupName.c_str(), H5P_DEFAULT) > 0) {
        file.unlink(groupName);
    }
}
