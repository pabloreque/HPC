/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Hybrid-mode input keys: grid, electron model, population and wave.

An nx= line switches the input file to hybrid mode (see input.hpp):

  nx            number of cells (the box length comes from length=)
  iterations    ICN corrections per step (default 2)
  output_every  field snapshot period in steps (default 1)
  Te, eta       electron temperature and resistivity (default 0)
  nu            hyper-resistivity (default 0)
  ppc           particles per cell
  density       ion density (default 1)
  vth           thermal speed, per-component std dev (default 0: cold)
  seed          random seed of the thermal velocities (default 1)
  wave_mode     wavelengths in the box of the eigenmode (default 0: none)
  wave_amplitude, wave_branch (R or L)
*/
#pragma once

#include <stdexcept>
#include <string>

#include "field_solver.hpp"
#include "loader.hpp"
#include "parse.hpp"
#include "vector3.hpp"

struct HybridSettings {
    int nx = 0;
    int iterations = 2;
    int outputEvery = 1;
    ElectronModel electrons;
    PopulationSettings population;
    WaveSettings wave;
};

// Throws invalid_argument with 'message' if 'condition' is false.
inline void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::invalid_argument(message);
    }
}

// Parses "R" or "L". Throws invalid_argument otherwise.
inline WaveBranch parseWaveBranch(const std::string &value) {
    if (value == "R") {
        return WaveBranch::R;
    }
    if (value == "L") {
        return WaveBranch::L;
    }
    throw std::invalid_argument("bad value for 'wave_branch', expected R or L");
}

// Parses a non-negative integer seed. Throws invalid_argument otherwise.
inline unsigned parseSeed(const std::string &value, const std::string &key) {
    const int seed = parseInt(value, key);
    // Read as int first: a negative value would silently wrap to a huge
    // unsigned instead of being rejected.
    require(seed >= 0, "seed must be >= 0");
    return static_cast<unsigned>(seed);
}

// Handles nx, iterations, output_every, Te, eta, nu. Returns false
// otherwise.
inline bool applySolverKey(HybridSettings &settings, const std::string &key,
                           const std::string &value) {
    if (key == "nx") {
        settings.nx = parseInt(value, key);
    } else if (key == "iterations") {
        settings.iterations = parseInt(value, key);
    } else if (key == "output_every") {
        settings.outputEvery = parseInt(value, key);
    } else if (key == "Te") {
        settings.electrons.temperature = parseDouble(value, key);
    } else if (key == "eta") {
        settings.electrons.resistivity = parseDouble(value, key);
    } else if (key == "nu") {
        settings.electrons.hyperResistivity = parseDouble(value, key);
    } else {
        return false;
    }
    return true;
}

// Handles ppc, density, vth, seed. Returns false otherwise.
inline bool applyPopulationKey(PopulationSettings &population,
                               const std::string &key,
                               const std::string &value) {
    if (key == "ppc") {
        population.ppc = parseInt(value, key);
    } else if (key == "density") {
        population.density = parseDouble(value, key);
    } else if (key == "vth") {
        population.vth = parseDouble(value, key);
    } else if (key == "seed") {
        population.seed = parseSeed(value, key);
    } else {
        return false;
    }
    return true;
}

// Handles wave_mode, wave_amplitude, wave_branch. Returns false otherwise.
inline bool applyWaveKey(WaveSettings &wave, const std::string &key,
                         const std::string &value) {
    if (key == "wave_mode") {
        wave.mode = parseInt(value, key);
    } else if (key == "wave_amplitude") {
        wave.amplitude = parseDouble(value, key);
    } else if (key == "wave_branch") {
        wave.branch = parseWaveBranch(value);
    } else {
        return false;
    }
    return true;
}

// Handles every hybrid-mode key. Returns false if 'key' is not one.
inline bool applyHybridKey(HybridSettings &settings, const std::string &key,
                           const std::string &value) {
    return applySolverKey(settings, key, value) ||
           applyPopulationKey(settings.population, key, value) ||
           applyWaveKey(settings.wave, key, value);
}

// Range checks of the grid, time scheme and electron model settings.
inline void validateSolverSettings(const HybridSettings &settings) {
    require(settings.nx > 0, "nx must be > 0");
    require(settings.iterations >= 1, "iterations must be >= 1");
    require(settings.outputEvery >= 1, "output_every must be >= 1");
    require(settings.electrons.temperature >= 0.0, "Te must be >= 0");
    require(settings.electrons.resistivity >= 0.0, "eta must be >= 0");
    require(settings.electrons.hyperResistivity >= 0.0, "nu must be >= 0");
}

// Range checks of the ion population settings.
inline void validatePopulationSettings(const PopulationSettings &population) {
    require(population.ppc > 0, "ppc must be > 0");
    require(population.density > 0.0, "density must be > 0");
    require(population.vth >= 0.0, "vth must be >= 0");
}

// The eigenmode relation in loader.hpp assumes B0 = (1, 0, 0) and n0 = 1.
inline void validateWaveSettings(const HybridSettings &settings,
                                 const Vector3 &B0) {
    require(settings.wave.mode >= 0, "wave_mode must be >= 0");
    if (settings.wave.mode == 0) {
        require(settings.wave.amplitude == 0.0,
                "wave_amplitude needs wave_mode > 0");
        return;
    }
    require(B0.x == 1.0 && B0.y == 0.0 && B0.z == 0.0,
            "a wave needs B=1 0 0 (normalized eigenmode)");
    require(settings.population.density == 1.0, "a wave needs density=1");
}

// Validates all hybrid settings. Throws invalid_argument on failure.
inline void validateHybridSettings(const HybridSettings &settings,
                                   const Vector3 &B0) {
    validateSolverSettings(settings);
    validatePopulationSettings(settings.population);
    validateWaveSettings(settings, B0);
}
