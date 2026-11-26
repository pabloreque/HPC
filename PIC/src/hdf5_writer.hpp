/*
Pablo Requeijo, September 24 2026 (Session 2); October 1 2026 (Session 3).
HPC coursework: PIC.
HDF5 output: trajectory datasets, field histories, run settings as
attributes, and file/group handling.
*/
#pragma once

#include <H5Cpp.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
#include "grid.hpp"
#include "hybrid_input.hpp"
#include "hybrid_step.hpp"
#include "input.hpp"
#include "loader.hpp"
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

// Opens 'path' read-write if it exists, to keep its other groups;
// otherwise creates it.
inline H5::H5File openOutputFile(const std::string &path) {
    return H5::H5File(path, std::filesystem::exists(path) ? H5F_ACC_RDWR
                                                          : H5F_ACC_TRUNC);
}

// Replaces group 'groupName' of 'file' by a new empty one.
inline H5::Group recreateGroup(H5::H5File &file, const std::string &groupName) {
    replaceGroup(file, groupName);
    return file.createGroup(groupName);
}

// Attaches a scalar double attribute to a group or dataset.
inline void writeAttribute(H5::H5Object &object, const std::string &name,
                           double value) {
    H5::Attribute attribute = object.createAttribute(
        name, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR));
    attribute.write(H5::PredType::NATIVE_DOUBLE, &value);
}

// Attaches a scalar integer attribute to a group or dataset.
inline void writeAttribute(H5::H5Object &object, const std::string &name,
                           int value) {
    H5::Attribute attribute = object.createAttribute(
        name, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
    attribute.write(H5::PredType::NATIVE_INT, &value);
}

// Attaches a 3-component double attribute to a group or dataset.
inline void writeAttribute(H5::H5Object &object, const std::string &name,
                           const Vector3 &value) {
    const double components[3] = {value.x, value.y, value.z};
    const hsize_t dims[1] = {3};
    H5::Attribute attribute = object.createAttribute(
        name, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dims));
    attribute.write(H5::PredType::NATIVE_DOUBLE, components);
}

// Attaches a string attribute to a group or dataset.
inline void writeAttribute(H5::H5Object &object, const std::string &name,
                           const std::string &value) {
    const H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
    H5::Attribute attribute =
        object.createAttribute(name, type, H5::DataSpace(H5S_SCALAR));
    attribute.write(type, value);
}

// Attaches a string attribute "centering" = "primal"/"dual" to a dataset.
inline void writeCenteringAttribute(H5::DataSet &dataset, Centering centering) {
    writeAttribute(
        dataset, "centering",
        std::string(centering == Centering::Primal ? "primal" : "dual"));
}

// Writes the charge q and mass m of a trajectory's particle.
inline void writeSpeciesAttributes(H5::Group &group, const Species &species) {
    writeAttribute(group, "q", species.charge);
    writeAttribute(group, "m", species.mass);
}

// Writes dt, steps, B, and length if the box is periodic.
inline void writeCommonAttributes(H5::Group &group, const Input &input) {
    writeAttribute(group, "dt", input.dt);
    writeAttribute(group, "steps", input.numSteps);
    writeAttribute(group, "B", input.B);
    if (input.length) {
        writeAttribute(group, "length", *input.length);
    }
}

// Writes nx, iterations, output_every and the electron model.
inline void writeSolverAttributes(H5::Group &group,
                                  const HybridSettings &settings) {
    writeAttribute(group, "nx", settings.nx);
    writeAttribute(group, "iterations", settings.iterations);
    writeAttribute(group, "output_every", settings.outputEvery);
    writeAttribute(group, "Te", settings.electrons.temperature);
    writeAttribute(group, "eta", settings.electrons.resistivity);
    writeAttribute(group, "nu", settings.electrons.hyperResistivity);
}

// Writes ppc, seed, and the uniform density and vth when no profile
// replaces them (a profile is in the t = 0 snapshot instead).
inline void writePopulationAttributes(H5::Group &group,
                                      const HybridSettings &settings) {
    const PopulationSettings &population = settings.population;
    writeAttribute(group, "ppc", population.ppc);
    writeAttribute(group, "seed", static_cast<int>(population.seed));
    if (population.profiles.density.empty()) {
        writeAttribute(group, "density", population.density);
    }
    if (population.profiles.thermalSpeed.empty()) {
        writeAttribute(group, "vth", population.vth);
    }
    if (!settings.profileFile.empty()) {
        writeAttribute(group, "profile_file", settings.profileFile);
    }
}

// Writes wave_mode, wave_amplitude and wave_branch.
inline void writeWaveAttributes(H5::Group &group, const WaveSettings &wave) {
    writeAttribute(group, "wave_mode", wave.mode);
    writeAttribute(group, "wave_amplitude", wave.amplitude);
    writeAttribute(group, "wave_branch",
                   std::string(wave.branch == WaveBranch::R ? "R" : "L"));
}

// Writes the input settings as attributes of 'group', named like the .txt
// keys, so the analysis reads them instead of copying them by hand.
inline void writeRunAttributes(H5::Group &group, const Input &input) {
    writeCommonAttributes(group, input);
    if (!input.hybrid) {
        writeAttribute(group, "E", input.E);
        return;
    }
    writeSolverAttributes(group, *input.hybrid);
    writePopulationAttributes(group, *input.hybrid);
    writeWaveAttributes(group, input.hybrid->wave);
}

constexpr std::size_t kSavedFieldCount = 11;

// Dataset names of the saved grid quantities, in snapshotFields() order.
constexpr std::array<const char *, kSavedFieldCount> kSavedFieldNames = {
    "Ex", "Ey", "Ez", "Bx", "By", "Bz", "n", "ux", "uy", "uz", "Ti"};

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
            &snapshot.moments.velocity.z,
            &snapshot.ionTemperature};
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

// Positions of the nx nodes of one centering.
inline std::vector<double> nodePositions(const Grid1D &grid,
                                         Centering centering) {
    std::vector<double> positions;
    for (int i = 0; i < grid.nx; ++i) {
        positions.push_back(nodePosition(i, centering, grid));
    }
    return positions;
}

// Writes the primal and dual node positions, x_primal and x_dual.
inline void writeNodePositions(H5::Group &group, const Grid1D &grid) {
    const hsize_t dims[1] = {static_cast<hsize_t>(grid.nx)};
    const H5::DataSpace dataspace(1, dims);
    writeDataset(group, "x_primal", nodePositions(grid, Centering::Primal),
                 dataspace);
    writeDataset(group, "x_dual", nodePositions(grid, Centering::Dual),
                 dataspace);
}

// Times of the snapshots, one value per snapshot.
inline std::vector<double>
snapshotTimes(const std::vector<FieldSnapshot> &history) {
    std::vector<double> time;
    for (const FieldSnapshot &snapshot : history) {
        time.push_back(snapshot.time);
    }
    return time;
}

// Writes the snapshot times t (one value per snapshot).
inline void writeSnapshotTimes(H5::Group &group,
                               const std::vector<FieldSnapshot> &history) {
    const hsize_t dims[1] = {static_cast<hsize_t>(history.size())};
    writeDataset(group, "t", snapshotTimes(history), H5::DataSpace(1, dims));
}

// Writes a hybrid run: t, x_primal, x_dual and one (snapshots x nx)
// dataset per grid quantity (Ex..Ez, Bx..Bz, n, ux..uz, Ti).
inline void writeFieldHistory(H5::Group &group,
                              const std::vector<FieldSnapshot> &history,
                              const Grid1D &grid) {
    writeSnapshotTimes(group, history);
    writeNodePositions(group, grid);
    for (std::size_t q = 0; q < kSavedFieldCount; ++q) {
        writeFieldDataset(group, history, q);
    }
}
