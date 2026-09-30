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
    std::vector<Vector3> position;    // vector<Vector3> is a vector of 3D vectors, 
    std::vector<Vector3> velocity;    // representing the position of the particle at each time step
};

inline Vector3 halfPositionDrift(const Vector3 &position, const Vector3 &velocity, double dt) {
    return position + (0.5 * dt) * velocity;
}

inline Vector3 halfElectricKick(const Vector3 &velocity, const Species &species, const Vector3 &E,
                                double dt) {
    return velocity + (species.charge * dt / (2.0 * species.mass)) * E;
}

inline Vector3 magneticRotation(const Vector3 &velocity, const Species &species, const Vector3 &B,
                                double dt) {
    const Vector3 t = (species.charge * dt / (2.0 * species.mass)) * B;
    const Vector3 vPrime = velocity + cross(velocity, t);
    const Vector3 s = (2.0 / (1.0 + dot(t, t))) * t;
    return velocity + cross(vPrime, s);
}

inline Particle borisStep(const Particle &state, const EMFields &fields, double dt) {
    const Vector3 xHalf = halfPositionDrift(state.position, state.velocity, dt);
    const Vector3 vMinus = halfElectricKick(state.velocity, state.species, fields.E, dt);
    const Vector3 vPlus = magneticRotation(vMinus, state.species, fields.B, dt);
    const Vector3 vNew = halfElectricKick(vPlus, state.species, fields.E, dt);
    const Vector3 xNew = halfPositionDrift(xHalf, vNew, dt);
    return {xNew, vNew, state.species};
}

inline Trajectory simulateBorisPusher(const Particle &initialState, const EMFields &fields, double dt,
                                      int numSteps) {
    Trajectory trajectory;
    trajectory.time.reserve(numSteps + 1);      // We want to reserve space for the number of steps in our
    trajectory.position.reserve(numSteps + 1);  // simulation and the initial state, hence numSteps + 1.
    trajectory.velocity.reserve(numSteps + 1);

    Particle state = initialState;
    trajectory.time.push_back(0.0);
    trajectory.position.push_back(state.position);
    trajectory.velocity.push_back(state.velocity);

    for (int step = 1; step <= numSteps; ++step) {
        state = borisStep(state, fields, dt);
        trajectory.time.push_back(step * dt);
        trajectory.position.push_back(state.position);
        trajectory.velocity.push_back(state.velocity);
    }

    return trajectory;
}
