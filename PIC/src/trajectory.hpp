/*
Pablo Requeijo, September 24 2026 (Session 2); November 26 2026.
HPC coursework: PIC.
Test-particle trajectories: recording a Boris run step by step, and the
columns it is saved as (t, x, y, z, vx, vy, vz). The HDF5 writer and the
Python module both take the columns from here, so a file and a NumPy dict
always use the same names.
*/
#pragma once

#include <array>
#include <cstddef> // std::size_t: sample count for the trajectory vectors
#include <optional>
#include <utility>
#include <vector>

#include "boris_pusher.hpp"
#include "vector3.hpp"

struct Trajectory {
    std::vector<double> time;
    // vector<Vector3> is a vector of 3D vectors, representing the
    // position (or velocity) of the particle at each time step
    std::vector<Vector3> position;
    std::vector<Vector3> velocity;
};

// Appends the time, position and velocity of one step to the trajectory.
inline void recordSample(Trajectory &trajectory, double time,
                         const Particle3D &state) {
    trajectory.time.push_back(time);
    trajectory.position.push_back(state.position);
    trajectory.velocity.push_back(state.velocity);
}

// Integrates numSteps Boris steps from initialState, recording every
// step plus the initial state. With periodicLength set, x is wrapped into
// [0, periodicLength) after every step.
inline Trajectory
simulateBorisPusher(const Particle3D &initialState, const EMFields &fields,
                    double dt, int numSteps,
                    const std::optional<double> &periodicLength = {}) {
    Trajectory trajectory;
    // We want to reserve space for the number of steps in our
    // simulation and the initial state, hence numSteps + 1.
    const std::size_t samples = static_cast<std::size_t>(numSteps) + 1;
    trajectory.time.reserve(samples);
    trajectory.position.reserve(samples);
    trajectory.velocity.reserve(samples);

    Particle3D state = initialState;
    recordSample(trajectory, 0.0, state);
    for (int step = 1; step <= numSteps; ++step) {
        state = borisStep(state, fields, dt);
        applyPeriodicX(state, periodicLength);
        recordSample(trajectory, step * dt, state);
    }
    return trajectory;
}

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

// One saved column of a trajectory: its dataset name and one value per
// recorded step.
struct TrajectoryColumn {
    const char *name;
    std::vector<double> values;
};

// The 7 columns t, x, y, z, vx, vy, vz of a trajectory, in saved order.
inline std::array<TrajectoryColumn, 7>
trajectoryColumns(const Trajectory &trajectory) {
    auto position = splitComponents(trajectory.position);
    auto velocity = splitComponents(trajectory.velocity);
    return {TrajectoryColumn{"t", trajectory.time},
            TrajectoryColumn{"x", std::move(position[0])},
            TrajectoryColumn{"y", std::move(position[1])},
            TrajectoryColumn{"z", std::move(position[2])},
            TrajectoryColumn{"vx", std::move(velocity[0])},
            TrajectoryColumn{"vy", std::move(velocity[1])},
            TrajectoryColumn{"vz", std::move(velocity[2])}};
}
