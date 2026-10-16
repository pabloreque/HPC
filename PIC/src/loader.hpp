/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Initial conditions: proton population from initial moments on the mesh,
optional wave.

Population, loaded in three steps:
1. Initial moments: the ion density n, bulk velocity u and thermal speed
   vth, each either uniform (density=, u = 0, vth=) or given on the
   primal nodes x_i = i*dx (NodeProfiles, e.g. read from a profile file,
   or later passed from Python) and linearly interpolated in between,
   with the same shape function as the deposit.
2. 'ppc' protons (q = m = 1) per cell, evenly spaced in x.
3. Each macroparticle at x_p gets the weight w_p = n(x_p)/ppc and the
   velocity v_p = u(x_p) + vth(x_p) * N(0, 1) per component (fixed seed,
   reproducible), plus the wave velocity if any. Depositing the loaded
   particles gives back the mesh moments: exactly for uniform ones, to
   O(dx^2) otherwise, plus the statistical noise of the thermal draws.
   The weights never change during the run: a macroparticle always
   stands for the same number of physical ions.

Wave: a circularly polarized parallel eigenmode of the cold-ion hybrid
equations with B0 = x (|B0| = 1, n0 = 1). With b = By + i*Bz and
w = uy + i*uz, linearizing the ion momentum, Ohm's (Te = 0) and Faraday's
laws for b = A exp(i(kx - omega t)) gives

  w = -(k/omega) b,     omega^2 + k^2 omega - k^2 = 0,

whose roots are the two branches (frequencies in units of Omega_ci):
  L (ion cyclotron, rotates with the ions, propagates +x):
      omega = +(k^2/2) (sqrt(1 + 4/k^2) - 1)
  R (whistler, rotates with the electrons, propagates -x):
      omega = -(k^2/2) (sqrt(1 + 4/k^2) + 1)

Why a wave: the uniform plasma is an equilibrium where nothing happens.
A small eigenmode on top of it has a known analytic frequency, so
measuring it validates the whole hybrid loop at once (pusher, deposit,
interpolation, Ampere, Faraday, Ohm, ICN). Loading the matching ion
velocity, not only the B perturbation, excites a single clean branch
instead of a mix of R and L.
*/
#pragma once

#include <algorithm>
#include <cmath>
#include <complex> // std::complex: wave as b = By + i*Bz
#include <cstddef>
#include <optional>
#include <random>
#include <vector>

#include "boris_pusher.hpp"
#include "grid.hpp"
#include "shape.hpp"
#include "vector3.hpp"

constexpr double kTwoPi = 6.28318530717958647692;

// Initial ion moments on the nx primal nodes x_i = i*dx. An empty vector
// means that moment is uniform (see PopulationSettings).
struct NodeProfiles {
    std::vector<double> density;
    std::vector<Vector3> bulkVelocity;
    std::vector<double> thermalSpeed;
};

struct PopulationSettings {
    // Macroparticles per cell.
    int ppc = 0;
    // Uniform ion density n0, used when profiles.density is empty.
    double density = 1.0;
    // Uniform thermal speed, used when profiles.thermalSpeed is empty:
    // per-component std dev of the Maxwellian, in v_A (T_i = vth^2,
    // beta_i = 2 n vth^2/|B0|^2); 0 for a cold population.
    double vth = 0.0;
    // Seed of the thermal velocity draws (same seed, same run).
    unsigned seed = 1;
    // Non-uniform initial moments; the bulk velocity is 0 when not given.
    NodeProfiles profiles;
};

enum class WaveBranch { R, L };

struct WaveSettings {
    // Number of wavelengths in the box (k = 2*pi*mode/L); 0 = no wave.
    int mode = 0;
    double amplitude = 0.0;
    WaveBranch branch = WaveBranch::R;
};

// Wavenumber k = 2*pi*mode/L of the perturbation.
inline double waveNumber(const WaveSettings &wave, const Grid1D &grid) {
    return kTwoPi * wave.mode / grid.length;
}

// Signed eigenfrequency of the branch for exp(i(kx - omega t)), k > 0:
// the L (+) and R (-) roots of omega^2 + k^2 omega - k^2 = 0, derived in
// the header comment of this file.
inline double waveFrequency(double k, WaveBranch branch) {
    const double root = std::sqrt(1.0 + 4.0 / (k * k));
    return branch == WaveBranch::L ? 0.5 * k * k * (root - 1.0)
                                   : -0.5 * k * k * (root + 1.0);
}

// Complex transverse magnetic perturbation b = By + i*Bz at x, t = 0.
inline std::complex<double> waveB(double x, const WaveSettings &wave,
                                  const Grid1D &grid) {
    const double k = waveNumber(wave, grid);
    return wave.amplitude * std::polar(1.0, k * x);
}

// Ion velocity perturbation (0, Re w, Im w) with w = -(k/omega) b at x.
inline Vector3 waveVelocity(double x, const WaveSettings &wave,
                            const Grid1D &grid) {
    if (wave.mode == 0) {
        return Vector3{};
    }
    const double k = waveNumber(wave, grid);
    const std::complex<double> w =
        -(k / waveFrequency(k, wave.branch)) * waveB(x, wave, grid);
    return Vector3{0.0, w.real(), w.imag()};
}

// Uniform B0 on the Yee mesh plus the wave perturbation on By, Bz (dual).
inline VecField initialMagneticField(const Vector3 &B0,
                                     const WaveSettings &wave,
                                     const Grid1D &grid) {
    VecField B = makeMagneticField(grid);
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, B.y.centering, grid);
        const std::complex<double> b =
            wave.mode == 0 ? 0.0 : waveB(x, wave, grid);
        B.x.values[i] = B0.x;
        B.y.values[i] = B0.y + b.real();
        B.z.values[i] = B0.z + b.imag();
    }
    return B;
}

// Draws a Maxwellian thermal velocity (std dev vth per component); zero
// for a cold population.
inline Vector3 thermalVelocity(std::mt19937 &rng, double vth) {
    if (vth == 0.0) {
        return Vector3{};
    }
    std::normal_distribution<double> normal(0.0, vth);
    const double vx = normal(rng);
    const double vy = normal(rng);
    const double vz = normal(rng);
    return Vector3{vx, vy, vz};
}

// A scalar initial moment: uniform, or given on the primal nodes and
// linearly interpolated in between. The uniform value is returned as is,
// so a uniform load stays exact.
struct ScalarProfile {
    double uniform = 0.0;
    std::optional<Field> nodes;

    // Value at position x.
    double at(double x, const Grid1D &grid) const {
        return nodes ? interpolate(*nodes, x, grid) : uniform;
    }
};

// A vector initial moment: zero, or given on the primal nodes and
// linearly interpolated in between.
struct VectorProfile {
    std::optional<VecField> nodes;

    // Value at position x.
    Vector3 at(double x, const Grid1D &grid) const {
        return nodes ? interpolate(*nodes, x, grid) : Vector3{};
    }
};

// Scalar profile from node values, or uniform if 'values' is empty.
inline ScalarProfile scalarProfile(const std::vector<double> &values,
                                   double uniform) {
    if (values.empty()) {
        return ScalarProfile{uniform, std::nullopt};
    }
    return ScalarProfile{uniform, Field{Centering::Primal, values}};
}

// Vector profile from node values, or zero if 'values' is empty.
inline VectorProfile vectorProfile(const std::vector<Vector3> &values,
                                   const Grid1D &grid) {
    if (values.empty()) {
        return VectorProfile{};
    }
    VecField field = makePrimalVecField(grid);
    for (std::size_t i = 0; i < values.size(); ++i) {
        field.x.values[i] = values[i].x;
        field.y.values[i] = values[i].y;
        field.z.values[i] = values[i].z;
    }
    return VectorProfile{field};
}

// Initial ion density on the primal nodes (uniform or from the profile).
inline Field initialDensity(const PopulationSettings &population,
                            const Grid1D &grid) {
    if (!population.profiles.density.empty()) {
        return Field{Centering::Primal, population.profiles.density};
    }
    Field density = makeField(Centering::Primal, grid);
    std::fill(density.values.begin(), density.values.end(), population.density);
    return density;
}

// Smallest initial ion density: where the whistler is fastest.
inline double minimumDensity(const PopulationSettings &population) {
    const std::vector<double> &profile = population.profiles.density;
    if (profile.empty()) {
        return population.density;
    }
    return *std::min_element(profile.begin(), profile.end());
}

// Loads protons one at a time from the initial moments: weight n/ppc,
// velocity u + vth * N(0, 1) + wave velocity, all evaluated at x.
struct ParticleLoader {
    ScalarProfile density;
    VectorProfile bulkVelocity;
    ScalarProfile thermalSpeed;
    WaveSettings wave;
    Grid1D grid;
    int ppc;
    std::mt19937 rng;

    // Builds the loader of 'population' (and 'wave') on 'grid'.
    ParticleLoader(const PopulationSettings &population,
                   const WaveSettings &wave, const Grid1D &grid)
        : density(
              scalarProfile(population.profiles.density, population.density)),
          bulkVelocity(vectorProfile(population.profiles.bulkVelocity, grid)),
          thermalSpeed(
              scalarProfile(population.profiles.thermalSpeed, population.vth)),
          wave(wave), grid(grid), ppc(population.ppc), rng(population.seed) {
    }

    // Loads one proton at x.
    Particle1D particleAt(double x) {
        const Species proton{1.0, 1.0};
        const Vector3 v = bulkVelocity.at(x, grid) +
                          thermalVelocity(rng, thermalSpeed.at(x, grid)) +
                          waveVelocity(x, wave, grid);
        return Particle1D{x, v, proton, density.at(x, grid) / ppc};
    }
};

// Loads ppc protons per cell, evenly spaced, from the initial moments.
inline std::vector<Particle1D> loadPlasma(const PopulationSettings &population,
                                          const WaveSettings &wave,
                                          const Grid1D &grid) {
    ParticleLoader loader(population, wave, grid);
    std::vector<Particle1D> particles;
    particles.reserve(static_cast<std::size_t>(grid.nx) * population.ppc);
    for (int cell = 0; cell < grid.nx; ++cell) {
        for (int k = 0; k < population.ppc; ++k) {
            const double x = (cell + (k + 0.5) / population.ppc) * grid.dx;
            particles.push_back(loader.particleAt(x));
        }
    }
    return particles;
}
