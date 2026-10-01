/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
Boris pusher step and trajectory simulation (generic, uniform fields).

Time-centered integration of dv/dt = (q/m)(E + v x B), dx/dt = v,
one step (r^n, v^n) -> (r^(n+1), v^(n+1)):

1. Half position drift: r^(n+1/2) = r^n + v^n*dt/2.

2. Half electric kick: v^- = v^n + (q*dt/2m)*E.

3. Magnetic rotation vector: t = (q*dt/2m)*B.

4. v' = v^- + v^- x t.

5. s = 2t/(1 + t.t).

6. v^+ = v^- + v' x s.

7. Second half electric kick: v^(n+1) = v^+ + (q*dt/2m)*E.

8. Closing half drift: r^(n+1) = r^(n+1/2) + v^(n+1)*dt/2.

Rotation never changes a vector's length, so |v| is preserved exactly
(to floating-point precision) for any dt.
*/
#pragma once

#include <cstddef> // std::size_t: sample count for the trajectory vectors
#include <vector>

#include "vector3.hpp"

struct Species {
    double charge;
    double mass;
};

struct Particle {
    Vector3 position;
    Vector3 velocity;
    Species species;
};

struct EMFields {
    Vector3 E;
    Vector3 B;
};

struct Trajectory {
    std::vector<double> time;
    // vector<Vector3> is a vector of 3D vectors, representing the
    // position (or velocity) of the particle at each time step
    std::vector<Vector3> position;
    std::vector<Vector3> velocity;
};

// Advances position by half a step at constant velocity: r + v*dt/2.
inline Vector3 halfPositionDrift(const Vector3 &position,
                                 const Vector3 &velocity, double dt) {
    return position + (0.5 * dt) * velocity;
}

// Accelerates velocity by half a step in E: v + (q*dt/2m)*E.
inline Vector3 halfElectricKick(const Vector3 &velocity, const Species &species,
                                const Vector3 &E, double dt) {
    return velocity + (species.charge * (0.5 * dt) / species.mass) * E;
}

// Rotates v^- around B into v^+ (steps 3-6 above); preserves |v|.
inline Vector3 magneticRotation(const Vector3 &vMinus, const Species &species,
                                const Vector3 &B, double dt) {
    const Vector3 t = (species.charge * (0.5 * dt) / species.mass) * B;
    const Vector3 vPrime = vMinus + cross(vMinus, t);
    const Vector3 s = (2.0 / (1.0 + dot(t, t))) * t;
    return vMinus + cross(vPrime, s);
}

// Advances one particle by one full Boris step (steps 1-8 above).
inline Particle borisStep(const Particle &state, const EMFields &fields,
                          double dt) {
    const Vector3 xHalf = halfPositionDrift(state.position, state.velocity, dt);
    const Vector3 vMinus =
        halfElectricKick(state.velocity, state.species, fields.E, dt);
    const Vector3 vPlus = magneticRotation(vMinus, state.species, fields.B, dt);
    const Vector3 vNew = halfElectricKick(vPlus, state.species, fields.E, dt);
    const Vector3 xNew = halfPositionDrift(xHalf, vNew, dt);
    return {xNew, vNew, state.species};
}

// Appends the time, position and velocity of one step to the trajectory.
inline void recordSample(Trajectory &trajectory, double time,
                         const Particle &state) {
    trajectory.time.push_back(time);
    trajectory.position.push_back(state.position);
    trajectory.velocity.push_back(state.velocity);
}

// Integrates numSteps Boris steps from initialState, recording every
// step plus the initial state.
inline Trajectory simulateBorisPusher(const Particle &initialState,
                                      const EMFields &fields, double dt,
                                      int numSteps) {
    Trajectory trajectory;
    // We want to reserve space for the number of steps in our
    // simulation and the initial state, hence numSteps + 1.
    const std::size_t samples = static_cast<std::size_t>(numSteps) + 1;
    trajectory.time.reserve(samples);
    trajectory.position.reserve(samples);
    trajectory.velocity.reserve(samples);

    Particle state = initialState;
    recordSample(trajectory, 0.0, state);
    for (int step = 1; step <= numSteps; ++step) {
        state = borisStep(state, fields, dt);
        recordSample(trajectory, step * dt, state);
    }
    return trajectory;
}
