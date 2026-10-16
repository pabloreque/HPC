/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Catch2 tests of the hybrid field equations (src/field_solver.hpp):
second-order Yee derivatives, Ampere/Faraday signs, div B, Ohm's law terms
(ideal, pressure, hyper-resistive) and a uniform equilibrium that must stay
exactly at rest.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "field_solver.hpp"
#include "grid.hpp"
#include "hybrid_step.hpp"
#include "loader.hpp"
#include "scenarios.hpp"
#include "tolerances.hpp"

namespace {

// Floating-point operations of Ohm's law at one node (interpolations,
// cross product, pressure gradient, divisions), rounded up.
constexpr double kOhmOperations = 100.0;

// Allowed error of a second-order spatial operator on a mode of
// wavenumber k and amplitude 'amplitude', with node spacing dx.
double spatialTolerance(double k, double dx, double amplitude) {
    return kSafetyFactor * spatialTruncation(k, dx) * amplitude;
}

// Field f(x_i) on the nodes of the given centering.
Field sampledField(const std::function<double(double)> &f, Centering centering,
                   const Grid1D &grid) {
    Field out = makeField(centering, grid);
    for (int i = 0; i < grid.nx; ++i) {
        out.values[i] = f(nodePosition(i, centering, grid));
    }
    return out;
}

// Largest |field(x_i) - exact(x_i)| over the nodes of 'field'.
double maxError(const Field &field, const std::function<double(double)> &exact,
                const Grid1D &grid) {
    double largest = 0.0;
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, field.centering, grid);
        largest = std::max(largest, std::fabs(field.values[i] - exact(x)));
    }
    return largest;
}

// Max error of d/dx sin(kx) from 'centering' on a grid of nx cells.
double derivativeError(int nx, Centering centering) {
    const Grid1D grid = makeGrid(nx, 2.0 * kPi);
    const Field f =
        sampledField([](double x) { return std::sin(x); }, centering, grid);
    return maxError(
        derivativeX(f, grid), [](double x) { return std::cos(x); }, grid);
}

// Uniform moments: density n0 and bulk velocity u0 on every node.
Moments uniformMoments(double n0, const Vector3 &u0, const Grid1D &grid) {
    Moments moments{makeField(Centering::Primal, grid),
                    makePrimalVecField(grid)};
    for (int i = 0; i < grid.nx; ++i) {
        moments.density.values[i] = n0;
        moments.velocity.x.values[i] = u0.x;
        moments.velocity.y.values[i] = u0.y;
        moments.velocity.z.values[i] = u0.z;
    }
    return moments;
}

// Uniform B0 on the Yee mesh.
VecField uniformB(const Vector3 &B0, const Grid1D &grid) {
    return initialMagneticField(B0, WaveSettings{}, grid);
}

// Largest |value| over the three components of a vector field.
double maxAbs(const VecField &field) {
    double largest = 0.0;
    for (const Field *f : {&field.x, &field.y, &field.z}) {
        for (double value : f->values) {
            largest = std::max(largest, std::fabs(value));
        }
    }
    return largest;
}

} // namespace

TEST_CASE("Yee derivatives are second-order accurate in both directions",
          "[fields][derivative]") {
    // Analytic reference (cos x): no reference error to subtract.
    const double expected =
        expectedConvergenceRatio(2, std::numeric_limits<double>::infinity());
    const double tolerance =
        convergenceRatioTolerance(2, spatialTruncation(1.0, 2.0 * kPi / 32));
    for (Centering centering : {Centering::Primal, Centering::Dual}) {
        const double ratio =
            derivativeError(32, centering) / derivativeError(64, centering);
        REQUIRE(ratio == Catch::Approx(expected).margin(tolerance));
    }
}

TEST_CASE("Ampere gives j = (0, -dBz/dx, dBy/dx) on the E centering",
          "[fields][ampere]") {
    const Grid1D grid = makeGrid(256, 2.0 * kPi);
    VecField B = uniformB({1.0, 0.0, 0.0}, grid);
    B.y = sampledField([](double x) { return std::sin(x); }, Centering::Dual,
                       grid);
    B.z = sampledField([](double x) { return std::cos(x); }, Centering::Dual,
                       grid);
    const VecField j = ampere(B, grid);
    REQUIRE(j.y.centering == Centering::Primal);
    REQUIRE(j.z.centering == Centering::Primal);
    REQUIRE(std::all_of(j.x.values.begin(), j.x.values.end(),
                        [](double v) { return v == 0.0; }));
    const double tolerance = spatialTolerance(1.0, grid.dx, 1.0);
    REQUIRE(
        maxError(j.y, [](double x) { return std::sin(x); }, grid) < tolerance);
    REQUIRE(
        maxError(j.z, [](double x) { return std::cos(x); }, grid) < tolerance);
}

TEST_CASE("Faraday gives dB/dt = -curl E and leaves Bx (div B) untouched",
          "[fields][faraday]") {
    const Grid1D grid = makeGrid(256, 2.0 * kPi);
    const VecField Bn = uniformB({1.0, 0.0, 0.0}, grid);
    VecField E = makeElectricField(grid);
    E.x = sampledField([](double x) { return std::sin(3.0 * x); },
                       Centering::Dual, grid);
    E.y = sampledField([](double x) { return std::sin(x); }, Centering::Primal,
                       grid);
    E.z = sampledField([](double x) { return std::cos(x); }, Centering::Primal,
                       grid);
    const double dt = 1e-3;
    const VecField B = faraday(Bn, E, dt, grid);
    REQUIRE(B.x.values == Bn.x.values);
    // dBy/dt = dEz/dx = -sin(x), dBz/dt = -dEy/dx = -cos(x).
    const Field dBy = addScaled(B.y, Bn.y, -1.0);
    const Field dBz = addScaled(B.z, Bn.z, -1.0);
    const double tolerance = spatialTolerance(1.0, grid.dx, dt);
    REQUIRE(maxError(
                dBy, [dt](double x) { return -dt * std::sin(x); }, grid) <
            tolerance);
    REQUIRE(maxError(
                dBz, [dt](double x) { return -dt * std::cos(x); }, grid) <
            tolerance);
}

TEST_CASE("Ohm's ideal term gives E = -u x B for a uniform flow",
          "[fields][ohm]") {
    const Grid1D grid = makeGrid(16, 4.0);
    const Moments moments = uniformMoments(1.0, {0.0, 0.2, 0.0}, grid);
    const VecField E =
        ohm(moments, uniformB({1.0, 0.0, 0.0}, grid), ElectronModel{}, grid);
    // -u x B = -(0.2 y) x (x) = 0.2 z, exact up to round-off.
    const double tolerance = roundoffTolerance(kOhmOperations, 0.2);
    for (int i = 0; i < grid.nx; ++i) {
        REQUIRE(E.x.values[i] == Catch::Approx(0.0).margin(tolerance));
        REQUIRE(E.y.values[i] == Catch::Approx(0.0).margin(tolerance));
        REQUIRE(E.z.values[i] == Catch::Approx(0.2).margin(tolerance));
    }
}

TEST_CASE("Ohm's pressure term gives Ex = -Te (dn/dx) / n", "[fields][ohm]") {
    const Grid1D grid = makeGrid(256, 2.0 * kPi);
    Moments moments = uniformMoments(1.0, {}, grid);
    moments.density =
        sampledField([](double x) { return 1.0 + 0.1 * std::sin(x); },
                     Centering::Primal, grid);
    const ElectronModel electrons{0.5, 0.0};
    const VecField E =
        ohm(moments, uniformB({1.0, 0.0, 0.0}, grid), electrons, grid);
    const auto exact = [](double x) {
        return -0.5 * 0.1 * std::cos(x) / (1.0 + 0.1 * std::sin(x));
    };
    // Largest |Ex|, at the steepest density drop.
    const double amplitude = 0.5 * 0.1 / (1.0 - 0.1);
    REQUIRE(E.x.centering == Centering::Dual);
    REQUIRE(maxError(E.x, exact, grid) <
            spatialTolerance(1.0, grid.dx, amplitude));
}

TEST_CASE("Ohm's hyper-resistive term gives E = -nu d^2j/dx^2",
          "[fields][ohm]") {
    const Grid1D grid = makeGrid(256, 2.0 * kPi);
    // No Bx and no flow: the ideal term vanishes from Ey.
    VecField B = uniformB({0.0, 0.0, 0.0}, grid);
    B.z = sampledField([](double x) { return std::cos(2.0 * x); },
                       Centering::Dual, grid);
    ElectronModel electrons;
    electrons.hyperResistivity = 0.01;
    const VecField E = ohm(uniformMoments(1.0, {}, grid), B, electrons, grid);
    // jy = -dBz/dx = 2 sin 2x, so Ey = -nu d^2jy/dx^2 = 8 nu sin 2x.
    const auto exact = [](double x) { return 0.08 * std::sin(2.0 * x); };
    REQUIRE(maxError(E.y, exact, grid) < spatialTolerance(2.0, grid.dx, 0.08));
}

TEST_CASE("A uniform equilibrium has E = 0 and stays exactly at rest",
          "[fields][equilibrium]") {
    const Grid1D grid = makeGrid(32, 8.0);
    const PopulationSettings population{5, 1.0, 0.0, 1};
    const ElectronModel electrons{0.3, 0.01, 0.001};
    const VecField B0 = uniformB({1.0, 0.3, -0.2}, grid);
    HybridState state = initializeState(
        grid, B0, loadPlasma(population, WaveSettings{}, grid), electrons);
    // E vanishes exactly in exact arithmetic. The largest round-off term
    // is the pressure gradient Te dn/dx / n of the deposit's round-off
    // in n, over one cell.
    const double densityRoundoff =
        roundoffTolerance(2.0 * population.ppc * kShapeOperations, 1.0);
    const double tolerance =
        electrons.temperature * densityRoundoff / grid.dx +
        roundoffTolerance(kOhmOperations, norm(Vector3{1.0, 0.3, -0.2}));
    REQUIRE(maxAbs(state.E) < tolerance);
    const HybridParameters params{0.01, 2, electrons};
    for (int step = 0; step < 50; ++step) {
        hybridStep(state, params);
    }
    REQUIRE(maxAbs(state.E) < tolerance);
    REQUIRE(state.B.y.values == B0.y.values);
    REQUIRE(state.B.z.values == B0.z.values);
}
