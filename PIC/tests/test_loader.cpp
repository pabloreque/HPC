/*
Pablo Requeijo, October 16 2026.
Session 3, HPC coursework: hybrid PIC.
Catch2 tests of the initial ion load (src/loader.hpp). Analytic moment
profiles n(x), u(x), vth(x) are sampled on the primal nodes, particles
are loaded from them, and the deposited moments are compared with the
analytic functions node by node.

A cold load is deterministic: its only error is the O((k dx)^2)
smoothing of interpolating (to load) and then depositing with the
linear shape function, which must converge at second order. A warm load
adds the statistical noise of the thermal draws, so its tolerances are
a few standard deviations of the sample mean and variance.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "grid.hpp"
#include "loader.hpp"
#include "shape.hpp"

namespace {

constexpr double kBoxLength = 4.0;
// One wavelength of every profile in the box.
constexpr double kProfileK = kTwoPi / kBoxLength;

// Analytic initial density.
double analyticDensity(double x) {
    return 1.0 + 0.5 * std::sin(kProfileK * x);
}

// Analytic initial bulk velocity.
Vector3 analyticBulkVelocity(double x) {
    return Vector3{0.1 * std::cos(kProfileK * x), 0.3 * std::sin(kProfileK * x),
                   -0.2 * std::cos(kProfileK * x)};
}

// Analytic initial thermal speed.
double analyticThermalSpeed(double x) {
    return 0.1 * (1.0 + 0.5 * std::cos(kProfileK * x));
}

// Population of ppc particles per cell whose node profiles sample the
// analytic n and u, and also vth for a warm one (cold otherwise).
PopulationSettings sampledPopulation(const Grid1D &grid, int ppc, bool warm) {
    PopulationSettings population{ppc, 1.0, 0.0, 1};
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, Centering::Primal, grid);
        population.profiles.density.push_back(analyticDensity(x));
        population.profiles.bulkVelocity.push_back(analyticBulkVelocity(x));
        if (warm) {
            population.profiles.thermalSpeed.push_back(analyticThermalSpeed(x));
        }
    }
    return population;
}

// Largest node error of each deposited moment against its analytic value.
struct LoadErrors {
    double density = 0.0;
    double bulkVelocity = 0.0;
    double thermalSpeed = 0.0;
};

// Largest node errors of the deposited n and u.
void addMomentErrors(LoadErrors &errors, const Moments &moments,
                     const Grid1D &grid) {
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, Centering::Primal, grid);
        const Vector3 u{moments.velocity.x.values[i],
                        moments.velocity.y.values[i],
                        moments.velocity.z.values[i]};
        errors.density =
            std::max(errors.density,
                     std::fabs(moments.density.values[i] - analyticDensity(x)));
        errors.bulkVelocity =
            std::max(errors.bulkVelocity, norm(u - analyticBulkVelocity(x)));
    }
}

// Largest node error of sqrt(T_i), which equals vth for protons.
void addThermalSpeedError(LoadErrors &errors, const Field &temperature,
                          const Grid1D &grid) {
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, Centering::Primal, grid);
        errors.thermalSpeed = std::max(
            errors.thermalSpeed, std::fabs(std::sqrt(temperature.values[i]) -
                                           analyticThermalSpeed(x)));
    }
}

// Loads a population sampled on nx cells, deposits it and returns the
// largest node errors of n, u and, for a warm load, vth.
LoadErrors loadErrors(int nx, int ppc, bool warm) {
    const Grid1D grid = makeGrid(nx, kBoxLength);
    const std::vector<Particle1D> particles =
        loadPlasma(sampledPopulation(grid, ppc, warm), WaveSettings{}, grid);
    LoadErrors errors;
    addMomentErrors(errors, depositMoments(particles, grid), grid);
    if (warm) {
        addThermalSpeedError(errors, depositTemperature(particles, grid), grid);
    }
    return errors;
}

} // namespace

TEST_CASE("A cold load deposits the analytic density and bulk velocity",
          "[loader]") {
    // 32 cells: k dx = 0.39, smoothing error ~ (k dx)^2/6 of the amplitude.
    const LoadErrors errors = loadErrors(32, 10, false);
    REQUIRE(errors.density < 5e-3);
    REQUIRE(errors.bulkVelocity < 5e-3);
}

TEST_CASE("Cold-load errors converge at second order in dx", "[loader]") {
    const LoadErrors coarse = loadErrors(16, 10, false);
    const LoadErrors fine = loadErrors(32, 10, false);
    REQUIRE(coarse.density / fine.density == Catch::Approx(4.0).margin(0.3));
    REQUIRE(coarse.bulkVelocity / fine.bulkVelocity ==
            Catch::Approx(4.0).margin(0.3));
}

TEST_CASE("Load weights follow the density mesh and conserve the charge",
          "[loader]") {
    const Grid1D grid = makeGrid(16, kBoxLength);
    const PopulationSettings population = sampledPopulation(grid, 20, false);
    const std::vector<Particle1D> particles =
        loadPlasma(population, WaveSettings{}, grid);
    const Moments moments = depositMoments(particles, grid);
    // Every cell holds ppc particles at the same offsets, so the total
    // deposited charge is exactly the sum of the mesh density.
    double deposited = 0.0, mesh = 0.0;
    for (int i = 0; i < grid.nx; ++i) {
        deposited += moments.density.values[i];
        mesh += population.profiles.density[i];
    }
    REQUIRE(deposited == Catch::Approx(mesh).epsilon(1e-12));
    // Heavier at the density peak (x = L/4) than at the trough (x = 3L/4).
    const std::size_t ppc = population.ppc;
    const std::size_t peakCell = grid.nx / 4, troughCell = 3 * grid.nx / 4;
    REQUIRE(particles[peakCell * ppc].weight >
            2.0 * particles[troughCell * ppc].weight);
}

TEST_CASE("A warm load deposits the analytic moments within the noise",
          "[loader]") {
    // 16 cells, 4000 particles per cell. Deterministic smoothing: ~0.013
    // (n), ~0.009 (u), ~0.007 (vth, including the bulk velocity spread
    // within one shape function). Noise of the sample mean and variance
    // over ~6000 effective particles per node: sigma ~0.003 (|u|) and
    // ~0.001 (vth); the bounds allow about 5 sigma on top.
    const LoadErrors errors = loadErrors(16, 4000, true);
    REQUIRE(errors.density < 0.015);
    REQUIRE(errors.bulkVelocity < 0.025);
    REQUIRE(errors.thermalSpeed < 0.015);
}
