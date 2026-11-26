/*
Pablo Requeijo, November 26 2026.
HPC coursework: PIC.
Field snapshots of a hybrid run: which grid quantities are saved at an
output time, in which order and under which names. The HDF5 writer and
the Python module both read them from here, so a file and a NumPy dict
always hold the same quantities under the same names.
*/
#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "grid.hpp"
#include "hybrid_step.hpp"
#include "shape.hpp"

// Grid quantities saved at one output time.
struct FieldSnapshot {
    double time = 0.0;
    VecField E;
    VecField B;
    Moments moments;
    // Ion temperature on primal nodes (diagnostic, see depositTemperature).
    Field ionTemperature;
};

// Copies the grid quantities of 'state' at 'time' and computes the ion
// temperature.
inline FieldSnapshot takeSnapshot(const HybridState &state, double time) {
    return FieldSnapshot{time, state.E, state.B, state.moments,
                         depositTemperature(state.particles, state.grid)};
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

// Times of the snapshots, one value per snapshot.
inline std::vector<double>
snapshotTimes(const std::vector<FieldSnapshot> &history) {
    std::vector<double> time;
    for (const FieldSnapshot &snapshot : history) {
        time.push_back(snapshot.time);
    }
    return time;
}
