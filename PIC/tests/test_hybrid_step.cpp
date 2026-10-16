/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Catch2 validation of the full hybrid PIC step (src/hybrid_step.hpp)
against the analytic dispersion relation of parallel circularly
polarized waves (ideal and hyper-resistive), plus the convergence order
of the ICN scheme.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>

#include "hybrid_step.hpp"
#include "loader.hpp"
#include "scenarios.hpp"

namespace {

// Wavenumber of a grid-scale wave: 4 cells per wavelength on 64 cells
// spanning 2 pi.
constexpr double kGridScaleK = 16.0;

// Fourier amplitude of b = By + i*Bz at wavenumber k (dual nodes).
std::complex<double> modeAmplitude(const VecField &B, const Grid1D &grid,
                                   double k) {
    std::complex<double> sum = 0.0;
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, B.y.centering, grid);
        sum += std::complex<double>(B.y.values[i], B.z.values[i]) *
               std::polar(1.0, -k * x);
    }
    return sum / static_cast<double>(grid.nx);
}

struct WaveMeasurement {
    double omega;
    double amplitudeRatio;
};

// Runs a wave scenario for 'periods' analytic periods and measures its
// frequency from the accumulated phase of the mode, b ~ exp(-i omega t).
WaveMeasurement measureWave(HybridScenario &scenario, double periods) {
    HybridState &state = scenario.state;
    const double k = waveNumber(scenario.wave, state.grid);
    const double omega = waveFrequency(k, scenario.wave.branch);
    const int steps = static_cast<int>(std::lround(
        periods * 2.0 * kPi / std::fabs(omega) / scenario.params.dt));
    const std::complex<double> initial = modeAmplitude(state.B, state.grid, k);
    std::complex<double> previous = initial;
    double phase = 0.0;
    for (int step = 0; step < steps; ++step) {
        hybridStep(state, scenario.params);
        const std::complex<double> current =
            modeAmplitude(state.B, state.grid, k);
        phase += std::arg(current / previous);
        previous = current;
    }
    return WaveMeasurement{-phase / (steps * scenario.params.dt),
                           std::abs(previous) / std::abs(initial)};
}

// Growth rate of the mode amplitude over 'duration': least-squares slope
// of ln|b(t)/b(0)|, which averages out the small beating with the other
// branch excited by the ideal (nu = 0) initial eigenvector.
double measureGrowthRate(HybridScenario &scenario, double duration) {
    HybridState &state = scenario.state;
    const double dt = scenario.params.dt;
    const double k = waveNumber(scenario.wave, state.grid);
    const std::complex<double> initial = modeAmplitude(state.B, state.grid, k);
    const int steps = static_cast<int>(std::lround(duration / dt));
    double sumT = 0.0, sumY = 0.0, sumTT = 0.0, sumTY = 0.0;
    for (int step = 1; step <= steps; ++step) {
        hybridStep(state, scenario.params);
        const double t = step * dt;
        const double y =
            std::log(std::abs(modeAmplitude(state.B, state.grid, k) / initial));
        sumT += t;
        sumY += y;
        sumTT += t * t;
        sumTY += t * y;
    }
    return (steps * sumTY - sumT * sumY) / (steps * sumTT - sumT * sumT);
}

// Exact complex frequency of the cold parallel modes with an effective
// resistivity eta' = nu k^2: omega^2 + a omega - a = 0, a = k^2 (1 + i eta').
std::complex<double> dissipativeFrequency(double k, double nu,
                                          WaveBranch branch) {
    const std::complex<double> a =
        k * k * std::complex<double>(1.0, nu * k * k);
    const std::complex<double> root = std::sqrt(a * a + 4.0 * a);
    return branch == WaveBranch::R ? (-a - root) / 2.0 : (-a + root) / 2.0;
}

// Measured growth rate of the given branch with hyper-resistivity 'nu'.
double hyperResistiveGrowthRate(WaveBranch branch, double nu) {
    ElectronModel electrons;
    electrons.hyperResistivity = nu;
    HybridScenario scenario =
        waveScenario(branch, 64, 0.002, 2, 0.01, electrons);
    return measureGrowthRate(scenario, 20.0);
}

// Largest |By|, |Bz| difference between two magnetic fields.
double maxDifference(const VecField &a, const VecField &b) {
    double largest = 0.0;
    for (std::size_t i = 0; i < a.y.values.size(); ++i) {
        largest = std::max(largest, std::fabs(a.y.values[i] - b.y.values[i]));
        largest = std::max(largest, std::fabs(a.z.values[i] - b.z.values[i]));
    }
    return largest;
}

// B after a fixed time 2.0 of the R-mode run with the given dt and number
// of ICN corrections (coarse grid, larger amplitude: time error dominates).
VecField finalField(double dt, int iterations) {
    HybridScenario scenario =
        waveScenario(WaveBranch::R, 32, dt, iterations, 0.05);
    const int steps = static_cast<int>(std::lround(2.0 / dt));
    for (int step = 0; step < steps; ++step) {
        hybridStep(scenario.state, scenario.params);
    }
    return scenario.state.B;
}

// Error ratio e(dt)/e(dt/2) against a dt/16 reference: ~2^order.
double convergenceRatio(int iterations) {
    const double dt = 0.004;
    const VecField reference = finalField(dt / 16.0, iterations);
    return maxDifference(finalField(dt, iterations), reference) /
           maxDifference(finalField(dt / 2.0, iterations), reference);
}

// Adds a circularly polarized perturbation b = amplitude exp(i k x) to
// By + i*Bz.
void addCircularPerturbation(VecField &B, const Grid1D &grid, double k,
                             double amplitude) {
    for (int i = 0; i < grid.nx; ++i) {
        const double x = nodePosition(i, B.y.centering, grid);
        const std::complex<double> b = std::polar(amplitude, k * x);
        B.y.values[i] += b.real();
        B.z.values[i] += b.imag();
    }
}

// Cold uniform plasma in B0 = x whose only perturbation is a grid-scale
// magnetic wave, 4 cells per wavelength (k = kGridScaleK on 64 cells).
HybridScenario gridScaleScenario(double nu) {
    ElectronModel electrons;
    electrons.hyperResistivity = nu;
    HybridScenario scenario =
        waveScenario(WaveBranch::R, 64, 0.002, 2, 0.0, electrons);
    HybridState &state = scenario.state;
    VecField B = state.B;
    addCircularPerturbation(B, state.grid, kGridScaleK, 1e-3);
    state = initializeState(state.grid, std::move(B),
                            std::move(state.particles), electrons);
    return scenario;
}

// Ratio |b_k(t)| / |b_k(0)| of the mode k after running for 'duration'.
double amplitudeRatioAfter(HybridScenario &scenario, double k,
                           double duration) {
    HybridState &state = scenario.state;
    const std::complex<double> initial = modeAmplitude(state.B, state.grid, k);
    const int steps =
        static_cast<int>(std::lround(duration / scenario.params.dt));
    for (int step = 0; step < steps; ++step) {
        hybridStep(state, scenario.params);
    }
    return std::abs(modeAmplitude(state.B, state.grid, k)) / std::abs(initial);
}

} // namespace

TEST_CASE("The R-mode (whistler) eigenmode has the analytic frequency",
          "[hybrid][wave]") {
    HybridScenario scenario = waveScenario(WaveBranch::R, 64, 0.002);
    const double k = waveNumber(scenario.wave, scenario.state.grid);
    const WaveMeasurement measured = measureWave(scenario, 3.0);
    // k d_i = 1: omega = -(1 + sqrt 5)/2.
    REQUIRE(waveFrequency(k, WaveBranch::R) ==
            Catch::Approx(-(1.0 + std::sqrt(5.0)) / 2.0));
    REQUIRE(measured.omega ==
            Catch::Approx(waveFrequency(k, WaveBranch::R)).epsilon(5e-3));
    REQUIRE(measured.amplitudeRatio == Catch::Approx(1.0).epsilon(1e-3));
    for (const Particle &p : scenario.state.particles) {
        REQUIRE(p.position.x >= 0.0);
        REQUIRE(p.position.x < scenario.state.grid.length);
    }
}

TEST_CASE("The L-mode (ion cyclotron) eigenmode has the analytic frequency",
          "[hybrid][wave]") {
    HybridScenario scenario = waveScenario(WaveBranch::L, 64, 0.002);
    const double k = waveNumber(scenario.wave, scenario.state.grid);
    const WaveMeasurement measured = measureWave(scenario, 2.0);
    // k d_i = 1: omega = (sqrt 5 - 1)/2.
    REQUIRE(waveFrequency(k, WaveBranch::L) ==
            Catch::Approx((std::sqrt(5.0) - 1.0) / 2.0));
    REQUIRE(measured.omega ==
            Catch::Approx(waveFrequency(k, WaveBranch::L)).epsilon(5e-3));
    REQUIRE(measured.amplitudeRatio == Catch::Approx(1.0).epsilon(1e-3));
}

TEST_CASE("ICN is second order with 2 corrections, first order with 1",
          "[hybrid][convergence]") {
    REQUIRE(convergenceRatio(2) == Catch::Approx(4.0).margin(0.3));
    REQUIRE(convergenceRatio(1) == Catch::Approx(2.0).margin(0.3));
}

TEST_CASE("Hyper-resistivity damps the R mode and drives the L mode at the "
          "analytic rates",
          "[hybrid][wave][hyper]") {
    const double nu = 0.002;
    // k d_i = 1, so eta' = nu: Im(omega) = -0.00234 (R), +0.00034 (L).
    const double gammaR = dissipativeFrequency(1.0, nu, WaveBranch::R).imag();
    const double gammaL = dissipativeFrequency(1.0, nu, WaveBranch::L).imag();
    REQUIRE(gammaR < 0.0);
    REQUIRE(gammaL > 0.0);
    REQUIRE(hyperResistiveGrowthRate(WaveBranch::R, nu) ==
            Catch::Approx(gammaR).epsilon(0.02));
    REQUIRE(hyperResistiveGrowthRate(WaveBranch::L, nu) ==
            Catch::Approx(gammaL).epsilon(0.15));
}

TEST_CASE("The timestep limit matches the ICN stability boundary",
          "[hybrid][stability]") {
    const Grid1D grid = makeGrid(64, 2.0 * kPi);
    // nu = 0: the classic whistler limit n dx^2 / (2|B0|).
    REQUIRE(timestepLimit(grid, 1.0, 1.0, 0.0) ==
            Catch::Approx(grid.dx * grid.dx / 2.0));
    // Exact ICN boundaries for this grid, from |g(dt (nu k^4 + i k^2))| = 1:
    // the estimate stays below them by less than 20 %.
    for (const auto &[nu, exact] :
         {std::pair{0.002, 0.00339}, std::pair{0.005, 0.00193},
          std::pair{0.01, 0.00110}}) {
        CAPTURE(nu);
        const double limit = timestepLimit(grid, 1.0, 1.0, nu);
        REQUIRE(limit < exact);
        REQUIRE(limit > 0.8 * exact);
    }
}

TEST_CASE("Small hyper-resistivity leaves the long waves unchanged",
          "[hybrid][wave][hyper]") {
    for (const double nu : {1e-4, 1e-3}) {
        for (const WaveBranch branch : {WaveBranch::R, WaveBranch::L}) {
            CAPTURE(nu, branch == WaveBranch::R ? "R" : "L");
            ElectronModel electrons;
            electrons.hyperResistivity = nu;
            HybridScenario scenario =
                waveScenario(branch, 64, 0.002, 2, 0.01, electrons);
            const double k = waveNumber(scenario.wave, scenario.state.grid);
            const WaveMeasurement measured = measureWave(scenario, 2.0);
            REQUIRE(measured.omega ==
                    Catch::Approx(waveFrequency(k, branch)).epsilon(5e-3));
            REQUIRE(measured.amplitudeRatio == Catch::Approx(1.0).margin(0.02));
        }
    }
}

TEST_CASE("Small hyper-resistivity removes grid-scale whistlers",
          "[hybrid][wave][hyper]") {
    // Ideal electrons: the grid-scale wave survives. Its ratio ~0.83 is the
    // ICN numerical damping, |g|^2 = 1 - y^4/4 + y^6/16 per step at
    // y = omega dt ~ 0.42, i.e. ~0.4 % per step over 50 steps.
    HybridScenario ideal = gridScaleScenario(0.0);
    REQUIRE(amplitudeRatioAfter(ideal, kGridScaleK, 0.1) > 0.5);
    // Same nu as the long-wave test: rate ~nu k^4 ~ 65 wipes it out.
    HybridScenario damped = gridScaleScenario(1e-3);
    REQUIRE(amplitudeRatioAfter(damped, kGridScaleK, 0.1) < 0.05);
}
