/*
Pablo Requeijo, November 26 2026.
HPC coursework: PIC.
Simulation runs shared by the C++ driver (main.cpp) and the Python module
(python/bindings.cpp), so both front ends simulate exactly the same thing:

- test-particle mode: one Boris trajectory per particle;
- hybrid mode: initial state, time loop with field snapshots, energy and
  div B diagnostics, stability warnings.
*/
#pragma once

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
#include "field_solver.hpp"
#include "grid.hpp"
#include "hybrid_input.hpp"
#include "hybrid_step.hpp"
#include "input.hpp"
#include "loader.hpp"
#include "parse.hpp"
#include "vector3.hpp"

// Integrates every particle of a test-particle input in its uniform E, B.
inline std::vector<Trajectory> simulateTestParticles(const Input &input) {
    const EMFields fields{input.E, input.B};
    std::vector<Trajectory> trajectories;
    for (const Particle3D &particle : input.particles) {
        trajectories.push_back(simulateBorisPusher(
            particle, fields, input.dt, input.numSteps, input.length));
    }
    return trajectories;
}

// A hybrid simulation in progress: its input, the current state, the
// number of steps done so far and the snapshots taken so far.
struct HybridRun {
    Input input;
    HybridParameters params;
    HybridState state;
    int step = 0;
    std::vector<FieldSnapshot> history;
};

// Range checks of a hybrid input, whether parsed from a .txt file or built
// in memory (Python). Throws invalid_argument on failure.
inline void validateHybridInput(const Input &input) {
    validateRunSettings(input);
    require(input.hybrid.has_value(), "not a hybrid input (no nx=)");
    require(input.length.has_value(), "missing required key 'length'");
    validateHybridSettings(*input.hybrid, input.B);
}

// Builds the initial hybrid state from the input settings.
inline HybridState initialHybridState(const Input &input) {
    const HybridSettings &settings = *input.hybrid;
    const Grid1D grid = makeGrid(settings.nx, *input.length);
    return initializeState(grid,
                           initialMagneticField(input.B, settings.wave, grid),
                           loadPlasma(settings.population, settings.wave, grid),
                           settings.electrons);
}

// Validates 'input', builds its initial state and takes the t = 0
// snapshot.
inline HybridRun startHybridRun(const Input &input) {
    validateHybridInput(input);
    const HybridSettings &settings = *input.hybrid;
    HybridRun run{
        input,
        HybridParameters{input.dt, settings.iterations, settings.electrons},
        initialHybridState(input)};
    run.history.push_back(takeSnapshot(run.state, 0.0));
    return run;
}

// Time reached by 'run'.
inline double currentTime(const HybridRun &run) {
    return run.step * run.input.dt;
}

// Advances 'run' by 'steps' hybrid steps, taking a snapshot every
// output_every steps counted from t = 0.
inline void advanceHybridRun(HybridRun &run, int steps) {
    const int outputEvery = run.input.hybrid->outputEvery;
    for (int k = 0; k < steps; ++k) {
        hybridStep(run.state, run.params);
        ++run.step;
        if (run.step % outputEvery == 0) {
            run.history.push_back(takeSnapshot(run.state, currentTime(run)));
        }
    }
}

// Magnetic energy sum(B^2/2) dx over the grid.
inline double magneticEnergy(const HybridState &state) {
    double sum = 0.0;
    for (const Field *component : {&state.B.x, &state.B.y, &state.B.z}) {
        for (double b : component->values) {
            sum += 0.5 * b * b;
        }
    }
    return sum * state.grid.dx;
}

// Ion kinetic energy sum(w m v^2/2) dx over all macroparticles.
inline double kineticEnergy(const HybridState &state) {
    double sum = 0.0;
    for (const Particle1D &particle : state.particles) {
        sum += 0.5 * particle.weight * particle.species.mass *
               dot(particle.velocity, particle.velocity);
    }
    return sum * state.grid.dx;
}

// Largest |dBx/dx|: the discrete div B of the 1D Yee mesh.
inline double maxDivergenceB(const HybridState &state) {
    const Field divB = derivativeX(state.B.x, state.grid);
    double largest = 0.0;
    for (double value : divB.values) {
        largest = std::max(largest, std::fabs(value));
    }
    return largest;
}

// Approximate stability limit on dt of a hybrid input (see timestepLimit).
inline double hybridTimestepLimit(const Input &input, const Grid1D &grid) {
    const HybridSettings &settings = *input.hybrid;
    return timestepLimit(grid, norm(input.B),
                         minimumDensity(settings.population),
                         settings.electrons.hyperResistivity);
}

// Warning if dt exceeds the approximate stability limit of the grid-scale
// whistler and hyper-resistive damping; empty otherwise.
inline std::string timestepWarning(const Input &input, const Grid1D &grid) {
    const double limit = hybridTimestepLimit(input, grid);
    if (input.dt <= limit) {
        return "";
    }
    std::ostringstream message;
    message << "dt=" << input.dt << " exceeds the stability limit ~" << limit
            << " (grid whistler"
            << (input.hybrid->electrons.hyperResistivity > 0.0
                    ? " + hyper-resistivity)"
                    : ", n dx^2 / 2|B0|)");
    return message.str();
}

// Warning if a single ICN correction is requested (first order, lets
// whistlers grow at every dt); empty otherwise.
inline std::string correctionWarning(const Input &input) {
    if (input.hybrid->iterations != 1) {
        return "";
    }
    return "iterations=1 is unstable for whistlers "
           "(|g|^2 = 1 + (omega dt)^4/4); use 2";
}

// Every stability warning that applies to a hybrid input on 'grid'.
inline std::vector<std::string> stabilityWarnings(const Input &input,
                                                  const Grid1D &grid) {
    std::vector<std::string> warnings;
    for (const std::string &warning :
         {timestepWarning(input, grid), correctionWarning(input)}) {
        if (!warning.empty()) {
            warnings.push_back(warning);
        }
    }
    return warnings;
}
