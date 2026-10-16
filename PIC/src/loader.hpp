/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Initial conditions: proton population on a density mesh, optional wave.

Population, loaded in three steps:
1. Density mesh: the initial ion density n is given on the primal nodes
   x_i = i*dx, either uniform (density=) or as one value per node
   (densityProfile, e.g. read from a file, or later passed from Python).
2. 'ppc' protons (q = m = 1) per cell, evenly spaced in x, with Maxwellian
   velocities of per-component standard deviation vth (fixed seed,
   reproducible).
3. Each macroparticle gets the weight w_p = n(x_p)/ppc, with n
   interpolated at x_p by the same linear shape function as the deposit,
   so depositing the loaded particles gives back the density mesh (exactly
   for a uniform n, to O(dx^2) otherwise). The weights never change
   during the run: a macroparticle always stands for the same number of
   physical ions.

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
#include <random>
#include <vector>

#include "boris_pusher.hpp"
#include "grid.hpp"
#include "shape.hpp"
#include "vector3.hpp"

constexpr double kTwoPi = 6.28318530717958647692;

struct PopulationSettings {
    // Macroparticles per cell.
    int ppc = 0;
    // Uniform ion density n0, used when densityProfile is empty.
    double density = 1.0;
    // Thermal speed: per-component std dev of the Maxwellian, in v_A
    // (T_i = vth^2, beta_i = 2 n vth^2/|B0|^2); 0 for a cold population.
    double vth = 0.0;
    // Seed of the thermal velocity draws (same seed, same run).
    unsigned seed = 1;
    // Initial ion density on the nx primal nodes x_i = i*dx; empty for a
    // uniform density.
    std::vector<double> densityProfile;
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

// Initial ion density on the primal nodes: the profile if given,
// otherwise the uniform density.
inline Field initialDensity(const PopulationSettings &population,
                            const Grid1D &grid) {
    if (!population.densityProfile.empty()) {
        return Field{Centering::Primal, population.densityProfile};
    }
    Field density = makeField(Centering::Primal, grid);
    std::fill(density.values.begin(), density.values.end(), population.density);
    return density;
}

// Smallest initial ion density: where the whistler is fastest.
inline double minimumDensity(const PopulationSettings &population) {
    if (population.densityProfile.empty()) {
        return population.density;
    }
    return *std::min_element(population.densityProfile.begin(),
                             population.densityProfile.end());
}

// Loads ppc protons per cell, evenly spaced, thermal plus wave velocity,
// each weighted by the density mesh at its position.
inline std::vector<Particle1D> loadPlasma(const PopulationSettings &population,
                                          const WaveSettings &wave,
                                          const Grid1D &grid) {
    std::mt19937 rng(population.seed);
    const Species proton{1.0, 1.0};
    const Field density = initialDensity(population, grid);
    std::vector<Particle1D> particles;
    particles.reserve(static_cast<std::size_t>(grid.nx) * population.ppc);
    for (int cell = 0; cell < grid.nx; ++cell) {
        for (int k = 0; k < population.ppc; ++k) {
            const double x = (cell + (k + 0.5) / population.ppc) * grid.dx;
            const Vector3 v = thermalVelocity(rng, population.vth) +
                              waveVelocity(x, wave, grid);
            const double weight =
                interpolate(density, x, grid) / population.ppc;
            particles.push_back(Particle1D{x, v, proton, weight});
        }
    }
    return particles;
}
