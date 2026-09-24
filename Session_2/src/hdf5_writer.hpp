/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
HDF5 output: trajectory datasets and group handling.
*/
#pragma once

#include <H5Cpp.h>

#include <array>
#include <cstdio>
#include <string>
#include <vector>

#include "boris_pusher.hpp"

inline std::array<std::vector<double>, 3> splitComponents(const std::vector<Vector3> &vectors) {
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

inline void writeDataset(H5::Group &group, const std::string &name, const std::vector<double> &data,
                         const H5::DataSpace &dataspace) {
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

inline void writeTrajectory(H5::H5Location &parent, const Trajectory &trajectory) {
    hsize_t dims[1] = {static_cast<hsize_t>(trajectory.time.size())};
    H5::DataSpace dataspace(1, dims);

    writeDataset(static_cast<H5::Group &>(parent), "t", trajectory.time, dataspace);

    const auto position = splitComponents(trajectory.position);
    writeDataset(static_cast<H5::Group &>(parent), "x", position[0], dataspace);
    writeDataset(static_cast<H5::Group &>(parent), "y", position[1], dataspace);
    writeDataset(static_cast<H5::Group &>(parent), "z", position[2], dataspace);

    const auto velocity = splitComponents(trajectory.velocity);
    writeDataset(static_cast<H5::Group &>(parent), "vx", velocity[0], dataspace);
    writeDataset(static_cast<H5::Group &>(parent), "vy", velocity[1], dataspace);
    writeDataset(static_cast<H5::Group &>(parent), "vz", velocity[2], dataspace);
}

inline bool fileExists(const std::string &path) {
    // Avoid <filesystem> link requirements on older toolchains.
    if (FILE *f = std::fopen(path.c_str(), "rb")) {
        std::fclose(f);
        return true;
    }
    return false;
}

inline void replaceGroup(H5::H5File &file, const std::string &groupName) {
    // groupName like "/cyclotron" or "cyclotron".
    const std::string name = groupName.front() == '/' ? groupName : "/" + groupName;
    if (H5Lexists(file.getId(), name.c_str(), H5P_DEFAULT) > 0) {
        file.unlink(name);
    }
}
