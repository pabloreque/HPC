/*
Pablo Requeijo, September 24 2026 (Session 2); October 1 2026 (Session 3).
HPC coursework: PIC.
Input of a run: initial conditions and settings of either mode, and their
range checks.

Scope: what an input is, independent of where it comes from. The .txt
parser (input_file.hpp) builds it from a file, the Python module
(python/bindings.cpp) from keyword arguments; both end in the same checks.
The hybrid-model settings and their checks are in hybrid_settings.hpp.
*/
#pragma once

#include <optional>
#include <vector>

#include "boris_pusher.hpp"
#include "hybrid_settings.hpp"
#include "parse.hpp"
#include "vector3.hpp"

struct Input {
    double dt = 0.0;
    int numSteps = 0;
    Vector3 E;
    Vector3 B;
    std::vector<Particle3D> particles;
    // Periodic box length along x; unset = unbounded (test-particle mode).
    std::optional<double> length;
    // Hybrid PIC settings, present only in hybrid mode (nx= line).
    std::optional<HybridSettings> hybrid;
};

// Range checks of dt, steps and length, shared by the .txt parser and
// inputs built in memory (Python). Throws invalid_argument on failure.
inline void validateRunSettings(const Input &input) {
    require(input.dt > 0.0, "dt must be > 0");
    require(input.numSteps >= 0, "steps must be >= 0");
    require(!input.length || *input.length > 0.0, "length must be > 0");
}

// Checks of a test-particle input built in memory: the parser checks the
// same line by line. Throws invalid_argument on failure.
inline void validateTestParticles(const Input &input) {
    validateRunSettings(input);
    require(!input.particles.empty(), "need at least one particle");
    for (const Particle3D &particle : input.particles) {
        require(particle.species.mass > 0.0, "m must be > 0");
    }
}
