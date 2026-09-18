// The Boris pusher: the standard time-centered algorithm for integrating the
// motion of a charged particle in uniform electric and magnetic fields,
//   dv/dt = (q/m) (E + v x B),   dx/dt = v,
// via a symmetric drift-kick-drift (Strang-split) leapfrog step: a half
// position drift, a velocity kick (two electric half-kicks around an exact
// magnetic rotation), then a second half position drift. This keeps
// position and velocity synchronized at integer time levels while the
// velocity kick effectively acts at the midpoint of the position step,
// which is what makes the position update second-order accurate. The
// rotation sub-step is what gives the method its other defining property:
// it preserves |v| exactly (to floating-point precision) for any time step,
// because rotating a vector never changes its length.
#pragma once

#include <vector>

#include "vector3.hpp"

struct ParticleState {
    Vector3 position;
    Vector3 velocity;
};

struct Trajectory {
    std::vector<double> time;
    std::vector<Vector3> position;
    std::vector<Vector3> velocity;
};

inline Vector3 halfElectricKick(const Vector3 &velocity, double charge, double mass, const Vector3 &E, double dt) {
    return velocity + (charge * dt / (2.0 * mass)) * E;
}

inline Vector3 magneticRotation(const Vector3 &velocity, double charge, double mass, const Vector3 &B, double dt) {
    const Vector3 t = (charge * dt / (2.0 * mass)) * B;
    const Vector3 vPrime = velocity + cross(velocity, t);
    const Vector3 s = (2.0 / (1.0 + dot(t, t))) * t;
    return velocity + cross(vPrime, s);
}

inline ParticleState borisStep(const ParticleState &state, double charge, double mass, const Vector3 &E,
                                const Vector3 &B, double dt) {
    const Vector3 xHalf = state.position + (0.5 * dt) * state.velocity;
    const Vector3 vMinus = halfElectricKick(state.velocity, charge, mass, E, dt);
    const Vector3 vPlus = magneticRotation(vMinus, charge, mass, B, dt);
    const Vector3 vNew = halfElectricKick(vPlus, charge, mass, E, dt);
    const Vector3 xNew = xHalf + (0.5 * dt) * vNew;
    return {xNew, vNew};
}

inline Trajectory simulateBorisPusher(const ParticleState &initialState, double charge, double mass,
                                       const Vector3 &E, const Vector3 &B, double dt, int numSteps) {
    Trajectory trajectory;
    trajectory.time.reserve(numSteps + 1);
    trajectory.position.reserve(numSteps + 1);
    trajectory.velocity.reserve(numSteps + 1);

    ParticleState state = initialState;
    trajectory.time.push_back(0.0);
    trajectory.position.push_back(state.position);
    trajectory.velocity.push_back(state.velocity);

    for (int step = 1; step <= numSteps; ++step) {
        state = borisStep(state, charge, mass, E, B, dt);
        trajectory.time.push_back(step * dt);
        trajectory.position.push_back(state.position);
        trajectory.velocity.push_back(state.velocity);
    }

    return trajectory;
}
