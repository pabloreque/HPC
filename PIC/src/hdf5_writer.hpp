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
#include "vector3.hpp"

// Splits vector<Vector3> into one double array per axis (x, y, z).
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

// Writes one 1D array of doubles to disk inside 'group'.
//   H5::Group         : the HDF5 "folder" the dataset will live in.
//   H5::DataSpace(1, dims): describes the shape: 1 dimension with N points.
//   createDataSet     : creates the dataset with a name, a type
//                       (NATIVE_DOUBLE = double in memory layout) and the shape.
//   dataset.write     : copies the contiguous buffer (data.data()) into the file.
inline void writeDataset(H5::Group &group, const std::string &name, const std::vector<double> &data,
                         const H5::DataSpace &dataspace) {
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

// Writes a full trajectory as 7 datasets (t, x, y, z, vx, vy, vz).
//   hsize_t dims[1] : HDF5's own size type; here, number of time steps.
//   H5::DataSpace   : one shared shape reused for all 7 datasets
//                     (all have exactly one value per time step).
//   H5::H5Location  : 'parent' arrives as this generic type so the function
//                     works for a top-level group (/cyclotron) and for
//                     subgroups (particle_N) alike.
//   static_cast<H5::Group&> : creating datasets is only possible on a Group,
//                     not on a generic location, so we tell the compiler to
//                     treat 'parent' as a Group.
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

// True if 'path' already exists on disk. main() uses it to open an
// existing .h5 read-write (keeping the other groups) or create a new
// one from scratch. Plain fopen instead of <filesystem> to avoid extra
// link requirements on older toolchains.
inline bool fileExists(const std::string &path) {
    if (FILE *f = std::fopen(path.c_str(), "rb")) {
        std::fclose(f);
        return true;
    }
    return false;
}

// Deletes 'groupName' if it already exists, so re-running the driver
// replaces the old trajectory instead of failing. Accepts both
// "/cyclotron" and "cyclotron" spellings.
//   H5Lexists : C-API check: does this link (group) already exist in the file?
//   unlink    : removes the link, freeing the name for createGroup.
inline void replaceGroup(H5::H5File &file, const std::string &groupName) {
    const std::string name = groupName.front() == '/' ? groupName : "/" + groupName;
    if (H5Lexists(file.getId(), name.c_str(), H5P_DEFAULT) > 0) {
        file.unlink(name);
    }
}
