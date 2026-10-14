/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
One hybrid PIC time step: iterated Crank-Nicolson predictor-corrector.

The PIC loop (move particles, accumulate moments, solve the
field equations) has a time-centering problem: the Boris pusher needs E
and B at n+1/2, but E comes from Ohm's law, which needs the
ion moments at the same time level as E. The iterated Crank-Nicolson
scheme breaks the loop by guessing the n+1 fields and
refining the guess:

  predict:  B* = Faraday(B^n, E^n),  E* = Ohm(moments^n, B*)
  repeat 'iterations' times (correct):
      Ebar = (E^n + E*)/2,  Bbar = (B^n + B*)/2     fields at n+1/2
      particles* = Boris(particles^n, Ebar, Bbar), x wrapped periodically
      moments*   = deposit(particles*)               moments at n+1
      B* = Faraday(B^n, Ebar),  E* = Ohm(moments*, B*)
  commit: state^(n+1) = (E*, B*, particles*, moments*)

Every correction restarts from the state at n, so the guesses converge to
the time-centered (Crank-Nicolson) solution; two corrections give the
second-order accuracy of the scheme.
*/
#pragma once

#include <cmath>
#include <utility>
#include <vector>

#include "boris_pusher.hpp"
#include "field_solver.hpp"
#include "grid.hpp"
#include "shape.hpp"

// Everything that evolves in time, at one time level: the fields, the
// particles and their moments on the grid.
struct HybridState {
    Grid1D grid;
    VecField E;
    VecField B;
    std::vector<Particle> particles;
    Moments moments;
};

// Fixed settings of the time integration.
struct HybridParameters {
    double dt = 0.0;
    // Number of ICN corrections per step: 1 is unstable for whistlers
    // (|g|^2 = 1 + (omega dt)^4/4), 2 (default) is stable and second order.
    int iterations = 2;
    ElectronModel electrons;
};

// Builds the initial state: deposits the moments and computes E^0 from
// Ohm's law, so the run starts with self-consistent fields.
inline HybridState initializeState(const Grid1D &grid, VecField B,
                                   std::vector<Particle> particles,
                                   const ElectronModel &electrons) {
    Moments moments = depositMoments(particles, grid);
    VecField E = ohm(moments, B, electrons, grid);
    return HybridState{grid, std::move(E), std::move(B), std::move(particles),
                       std::move(moments)};
}

// Field sampler interpolating Yee-mesh E and B at a particle position.
struct GridFieldSampler {
    const VecField &E;
    const VecField &B;
    const Grid1D &grid;

    // Interpolates E and B at 'position' with the linear shape function.
    EMFields operator()(const Vector3 &position) const {
        return gatherFields(E, B, position, grid);
    }
};

// Pushes every particle one Boris step in the grid fields E, B and wraps
// the new positions into the periodic box.
inline std::vector<Particle> pushAll(const std::vector<Particle> &particles,
                                     const VecField &E, const VecField &B,
                                     const Grid1D &grid, double dt) {
    const GridFieldSampler sampler{E, B, grid};
    std::vector<Particle> pushed;
    pushed.reserve(particles.size());
    for (const Particle &p : particles) {
        Particle next = borisStep(p, sampler, dt);
        next.position.x = wrapPosition(next.position.x, grid.length);
        pushed.push_back(next);
    }
    return pushed;
}

// Predictor: first guess of B^(n+1) and E^(n+1) from the fields at n,
// keeping the moments at n. Particles are filled in by correct().
inline HybridState predict(const HybridState &state,
                           const HybridParameters &params) {
    VecField B = faraday(state.B, state.E, params.dt, state.grid);
    VecField E = ohm(state.moments, B, params.electrons, state.grid);
    return HybridState{
        state.grid, std::move(E), std::move(B), {}, state.moments};
}

// Corrector: re-advances from 'state' (time n) using fields averaged with
// the current guess, and overwrites 'guess' with the improved n+1 state.
inline void correct(const HybridState &state, HybridState &guess,
                    const HybridParameters &params) {
    const Grid1D &grid = state.grid;
    const VecField Ebar = timeAverage(state.E, guess.E);
    const VecField Bbar = timeAverage(state.B, guess.B);
    guess.particles = pushAll(state.particles, Ebar, Bbar, grid, params.dt);
    guess.moments = depositMoments(guess.particles, grid);
    guess.B = faraday(state.B, Ebar, params.dt, grid);
    guess.E = ohm(guess.moments, guess.B, params.electrons, grid);
}

// Advances the hybrid state by one time step (predict, then correct
// params.iterations times, then commit).
inline void hybridStep(HybridState &state, const HybridParameters &params) {
    HybridState guess = predict(state, params);
    for (int k = 0; k < params.iterations; ++k) {
        correct(state, guess, params);
    }
    state = std::move(guess);
}

// Grid quantities saved at one output time.
struct FieldSnapshot {
    double time = 0.0;
    VecField E;
    VecField B;
    Moments moments;
};

// Copies the grid quantities of 'state' at 'time'.
inline FieldSnapshot takeSnapshot(const HybridState &state, double time) {
    return FieldSnapshot{time, state.E, state.B, state.moments};
}

// Approximate ICN (2 corrections) stability limit. The fastest grid mode,
// k = 2/dx, oscillates as a whistler (omega ~ k^2 |B0|/n) and decays by
// hyper-resistivity (rate nu k^4); dt*|nu k^4 + i omega| must stay inside
// the ICN stability region, of radius 2 on the imaginary axis (nu = 0:
// dt < n dx^2/(2|B0|)) and of radius >= 1.75 elsewhere.
inline double timestepLimit(const Grid1D &grid, double B0, double density,
                            double nu) {
    const double k2 = 4.0 / (grid.dx * grid.dx);
    const double omega = k2 * B0 / density;
    const double decay = nu * k2 * k2;
    const double radius = nu > 0.0 ? 1.75 : 2.0;
    return radius / std::hypot(omega, decay);
}
