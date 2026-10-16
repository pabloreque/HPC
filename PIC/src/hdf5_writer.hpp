/*
Pablo Requeijo, September 24 2026 (Session 2); October 1 2026 (Session 3).
HPC coursework: PIC.
HDF5 output: trajectory datasets, field histories and group handling.
*/
#pragma once

#include <H5Cpp.h>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
#include "grid.hpp"
#include "hybrid_step.hpp"
#include "vector3.hpp"

// Splits vector<Vector3> into one double array per axis (x, y, z).
inline std::array<std::vector<double>, 3>
splitComponents(const std::vector<Vector3> &vectors) {
    std::array<std::vector<double>, 3> components;
    components[0].reserve(vectors.size());
    components[1].reserve(vectors.size());
    components[2].reserve(vectors.size());
    for (const Vector3 &vector : vectors) {
        components[0].push_back(vector.x);
        components[1].push_back(vector.y);
        components[2].push_back(vector.z);
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

// Attaches a string attribute "centering" = "primal"/"dual" to a dataset.
inline void writeCenteringAttribute(H5::DataSet &dataset, Centering centering) {
    const H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
    H5::Attribute attribute =
        dataset.createAttribute("centering", type, H5::DataSpace(H5S_SCALAR));
    const std::string value =
        centering == Centering::Primal ? "primal" : "dual";
    attribute.write(type, value);
}

constexpr std::size_t kSavedFieldCount = 10;

// Dataset names of the saved grid quantities, in snapshotFields() order.
constexpr std::array<const char *, kSavedFieldCount> kSavedFieldNames = {
    "Ex", "Ey", "Ez", "Bx", "By", "Bz", "n", "ux", "uy", "uz"};

// The saved grid quantities of one snapshot, in kSavedFieldNames order.
inline std::array<const Field *, kSavedFieldCount>
snapshotFields(const FieldSnapshot &snapshot) {
    return {&snapshot.E.x,
            &snapshot.E.y,
            &snapshot.E.z,
            &snapshot.B.x,
            &snapshot.B.y,
            &snapshot.B.z,
            &snapshot.moments.density,
            &snapshot.moments.velocity.x,
            &snapshot.moments.velocity.y,
            &snapshot.moments.velocity.z};
}

// Stacks grid quantity number quantityIndex of every snapshot into a row-major
// (snapshots x nx) array.
inline std::vector<double>
stackSnapshots(const std::vector<FieldSnapshot> &history,
               std::size_t quantityIndex) {
    std::vector<double> data;
    for (const FieldSnapshot &snapshot : history) {
        const std::vector<double> &row =
            snapshotFields(snapshot)[quantityIndex]->values;
        data.insert(data.end(), row.begin(), row.end());
    }
    return data;
}

// Writes grid quantity number quantityIndex over time as a (snapshots x nx)
// dataset, tagged with its centering.
inline void writeFieldDataset(H5::Group &group,
                              const std::vector<FieldSnapshot> &history,
                              std::size_t quantityIndex) {
    const Field &first = *snapshotFields(history.front())[quantityIndex];
    const hsize_t dims[2] = {static_cast<hsize_t>(history.size()),
                             static_cast<hsize_t>(first.values.size())};
    H5::DataSet dataset = group.createDataSet(kSavedFieldNames[quantityIndex],
                                              H5::PredType::NATIVE_DOUBLE,
                                              H5::DataSpace(2, dims));
    dataset.write(stackSnapshots(history, quantityIndex).data(),
                  H5::PredType::NATIVE_DOUBLE);
    writeCenteringAttribute(dataset, first.centering);
}

// Writes the primal and dual node positions, x_primal and x_dual.
inline void writeNodePositions(H5::Group &group, const Grid1D &grid) {
    std::vector<double> primal, dual;
    for (int i = 0; i < grid.nx; ++i) {
        primal.push_back(nodePosition(i, Centering::Primal, grid));
        dual.push_back(nodePosition(i, Centering::Dual, grid));
    }
    const hsize_t dims[1] = {static_cast<hsize_t>(grid.nx)};
    const H5::DataSpace dataspace(1, dims);
    writeDataset(group, "x_primal", primal, dataspace);
    writeDataset(group, "x_dual", dual, dataspace);
}

// Writes the snapshot times t (one value per snapshot).
inline void writeSnapshotTimes(H5::Group &group,
                               const std::vector<FieldSnapshot> &history) {
    std::vector<double> time;
    for (const FieldSnapshot &snapshot : history) {
        time.push_back(snapshot.time);
    }
    const hsize_t dims[1] = {static_cast<hsize_t>(time.size())};
    writeDataset(group, "t", time, H5::DataSpace(1, dims));
}

// Writes a hybrid run: t, x_primal, x_dual and one (snapshots x nx)
// dataset per grid quantity (Ex..Ez, Bx..Bz, n, ux..uz).
inline void writeFieldHistory(H5::Group &group,
                              const std::vector<FieldSnapshot> &history,
                              const Grid1D &grid) {
    writeSnapshotTimes(group, history);
    writeNodePositions(group, grid);
    for (std::size_t q = 0; q < kSavedFieldCount; ++q) {
        writeFieldDataset(group, history, q);
    }
}
