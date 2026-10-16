/*
Pablo Requeijo, October 16 2026.
Session 3, HPC coursework: hybrid PIC.
Catch2 tests of the initial ion load (src/loader.hpp). Analytic moment
profiles n(x), u(x), vth(x) are sampled on the primal nodes, particles
are loaded from them, and the deposited moments are compared with the
analytic functions node by node.

A cold load is deterministic: its only error is the smoothing of
interpolating the node values to the particles and depositing them back
with the linear shape function, which is O(dx^2) with a coefficient
derived below. A warm load adds the statistical noise of the thermal
draws, whose standard deviation is also derived below.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "grid.hpp"
#include "loader.hpp"
#include "shape.hpp"
#include "tolerances.hpp"

namespace {

constexpr double kBoxLength = 4.0;
// One wavelength of every profile in the box.
constexpr double kProfileK = kTwoPi / kBoxLength;

// Profile mean + cosine cos(kx) + sine sin(kx), with k = kProfileK.
struct SinusoidalProfile {
    double mean;
    double cosine;
    double sine;

    // Profile value at x.
    double value(double x) const {
        return mean + cosine * std::cos(kProfileK * x) +
               sine * std::sin(kProfileK * x);
    }

    // First derivative at x.
    double slope(double x) const {
        return kProfileK * (sine * std::cos(kProfileK * x) -
                            cosine * std::sin(kProfileK * x));
    }

    // Second derivative at x.
    double curvature(double x) const {
        return -kProfileK * kProfileK * (value(x) - mean);
    }
};

const SinusoidalProfile kDensity{1.0, 0.0, 0.5};
const std::array<SinusoidalProfile, 3> kBulkVelocity = {
    SinusoidalProfile{0.0, 0.1, 0.0}, SinusoidalProfile{0.0, 0.0, 0.3},
    SinusoidalProfile{0.0, -0.2, 0.0}};
const SinusoidalProfile kThermalSpeed{0.1, 0.05, 0.0};

// Analytic initial bulk velocity.
Vector3 analyticBulkVelocity(double x) {
    return Vector3{kBulkVelocity[0].value(x), kBulkVelocity[1].value(x),
                   kBulkVelocity[2].value(x)};
}

// First derivative of the bulk velocity.
Vector3 bulkVelocitySlope(double x) {
    return Vector3{kBulkVelocity[0].slope(x), kBulkVelocity[1].slope(x),
                   kBulkVelocity[2].slope(x)};
}

// Second derivative of the bulk velocity.
Vector3 bulkVelocityCurvature(double x) {
    return Vector3{kBulkVelocity[0].curvature(x), kBulkVelocity[1].curvature(x),
                   kBulkVelocity[2].curvature(x)};
}

// Second moment of the smoothing kernel. Loading interpolates the node
// values linearly to each particle and the deposit spreads it back with
// the same hat: a profile f comes back as f + h f'', and a product of
// interpolated profiles f g loses h f' g'.
double smoothingScale(double dx) {
    return dx * dx / 6.0;
}

// Leading smoothing error of the deposited density, h n''.
double densitySmoothing(double x, double dx) {
    return smoothingScale(dx) * std::fabs(kDensity.curvature(x));
}

// Leading smoothing error of the deposited bulk velocity: the ratio of
// the deposits of n u and n gives h (u'' + n' u' / n).
double bulkVelocitySmoothing(double x, double dx) {
    const double densityGradient = kDensity.slope(x) / kDensity.value(x);
    return smoothingScale(dx) * norm(bulkVelocityCurvature(x) +
                                     densityGradient * bulkVelocitySlope(x));
}

// Leading smoothing error of sqrt(T_i): the deposited temperature is off
// by h (|u'|^2/3 + vth'^2 + 2 vth vth'' + 2 n' vth vth' / n), the first
// term being the spread of the bulk velocity across the kernel.
double thermalSpeedSmoothing(double x, double dx) {
    const double vth = kThermalSpeed.value(x);
    const double slope = kThermalSpeed.slope(x);
    const double temperatureError =
        smoothingScale(dx) *
        (dot(bulkVelocitySlope(x), bulkVelocitySlope(x)) / 3.0 + slope * slope +
         2.0 * vth * kThermalSpeed.curvature(x) +
         2.0 * kDensity.slope(x) * vth * slope / kDensity.value(x));
    return std::fabs(temperatureError) / (2.0 * vth);
}

// Effective number of independent particles behind one node moment: the
// node collects ppc particles from each neighbouring cell with hat
// weights S, so (sum S)^2 / sum S^2 = ppc^2 / (2 ppc / 3).
double effectiveParticles(int ppc) {
    return 1.5 * ppc;
}

// Standard deviation of the deposited |u - u_exact|: vth / sqrt(N) per
// component, combined over the three components.
double bulkVelocityNoise(double x, int ppc) {
    return std::sqrt(3.0) * kThermalSpeed.value(x) /
           std::sqrt(effectiveParticles(ppc));
}

// Standard deviation of the deposited sqrt(T_i): the sample variance of
// 3N Gaussian draws has relative spread sqrt(2/(3N)), halved by the root.
double thermalSpeedNoise(double x, int ppc) {
    return kThermalSpeed.value(x) / std::sqrt(6.0 * effectiveParticles(ppc));
}

// Population of ppc particles per cell whose node profiles sample the
// analytic n and u, and also vth for a warm one (cold otherwise).
PopulationSettings sampledPopulation(const Grid1D &grid, int ppc, bool warm) {
    PopulationSettings population{ppc, 1.0, 0.0, 1};
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, Centering::Primal, grid);
        population.profiles.density.push_back(kDensity.value(x));
        population.profiles.bulkVelocity.push_back(analyticBulkVelocity(x));
        if (warm) {
            population.profiles.thermalSpeed.push_back(kThermalSpeed.value(x));
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
                     std::fabs(moments.density.values[i] - kDensity.value(x)));
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
                                           kThermalSpeed.value(x)));
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

// Allowed node errors of a load on nx cells with ppc particles per cell:
// kSafetyFactor times the largest smoothing error, plus kNoiseSigmas
// standard deviations of the thermal noise for a warm load.
LoadErrors loadTolerances(int nx, int ppc, bool warm) {
    const Grid1D grid = makeGrid(nx, kBoxLength);
    LoadErrors tolerances;
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, Centering::Primal, grid);
        const double noise = warm ? kNoiseSigmas : 0.0;
        tolerances.density = std::max(
            tolerances.density, kSafetyFactor * densitySmoothing(x, grid.dx));
        tolerances.bulkVelocity =
            std::max(tolerances.bulkVelocity,
                     kSafetyFactor * bulkVelocitySmoothing(x, grid.dx) +
                         noise * bulkVelocityNoise(x, ppc));
        tolerances.thermalSpeed =
            std::max(tolerances.thermalSpeed,
                     kSafetyFactor * thermalSpeedSmoothing(x, grid.dx) +
                         noise * thermalSpeedNoise(x, ppc));
    }
    return tolerances;
}

} // namespace

TEST_CASE("A cold load deposits the analytic density and bulk velocity",
          "[loader]") {
    const LoadErrors errors = loadErrors(32, 10, false);
    const LoadErrors tolerances = loadTolerances(32, 10, false);
    REQUIRE(errors.density < tolerances.density);
    REQUIRE(errors.bulkVelocity < tolerances.bulkVelocity);
}

TEST_CASE("Cold-load errors converge at second order in dx", "[loader]") {
    const LoadErrors coarse = loadErrors(16, 10, false);
    const LoadErrors fine = loadErrors(32, 10, false);
    // Analytic reference: no reference error to subtract.
    const double expected =
        expectedConvergenceRatio(2, std::numeric_limits<double>::infinity());
    const double tolerance = convergenceRatioTolerance(
        2, spatialTruncation(kProfileK, kBoxLength / 16));
    REQUIRE(coarse.density / fine.density ==
            Catch::Approx(expected).margin(tolerance));
    REQUIRE(coarse.bulkVelocity / fine.bulkVelocity ==
            Catch::Approx(expected).margin(tolerance));
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
    REQUIRE(deposited ==
            Catch::Approx(mesh).margin(roundoffTolerance(
                particles.size() * kShapeOperations + grid.nx, mesh)));
    // Each weight is n(x_p)/ppc, with n linearly interpolated between
    // nodes: off by at most dx^2 |n''| / 8.
    const double interpolationError =
        grid.dx * grid.dx * kProfileK * kProfileK * kDensity.sine / 8.0;
    for (const Particle1D &p : particles) {
        CAPTURE(p.position);
        REQUIRE(p.weight * population.ppc ==
                Catch::Approx(kDensity.value(p.position))
                    .margin(kSafetyFactor * interpolationError));
    }
}

TEST_CASE("A warm load deposits the analytic moments within the noise",
          "[loader]") {
    const LoadErrors errors = loadErrors(16, 4000, true);
    const LoadErrors tolerances = loadTolerances(16, 4000, true);
    // The positions and weights do not depend on the random draws, so the
    // density has no noise.
    REQUIRE(errors.density < tolerances.density);
    REQUIRE(errors.bulkVelocity < tolerances.bulkVelocity);
    REQUIRE(errors.thermalSpeed < tolerances.thermalSpeed);
}
