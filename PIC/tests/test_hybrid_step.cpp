/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Catch2 validation of the full hybrid PIC step (src/hybrid_step.hpp)
against the analytic dispersion relation of parallel circularly
polarized waves (ideal and hyper-resistive), plus the convergence order
and the stability of the ICN scheme: a grid-scale perturbation must decay
below the timestep limit and grow above it (or at any dt with a single
correction), at the analytic ICN rate. TESTS.md gives the physics behind
each check and the origin of each tolerance.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <utility>

#include "hybrid_step.hpp"
#include "loader.hpp"
#include "scenarios.hpp"
#include "tolerances.hpp"

namespace {

// Grid of every wave test: 64 cells spanning 2 pi (one wavelength of the
// k d_i = 1 eigenmodes).
constexpr int kWaveCells = 64;
constexpr double kWaveDx = 2.0 * kPi / kWaveCells;
// Timestep of the wave tests, well below the timestep limit (0.0048).
constexpr double kWaveDt = 0.002;
// Wavenumber of a grid-scale wave: 4 cells per wavelength.
constexpr double kGridScaleK = 16.0;
// Shortest wave of the same grid, 2 cells per wavelength: the fastest
// whistler, which sets the timestep limit.
constexpr double kNyquistK = 32.0;
// Steps of the Nyquist stability runs: enough to separate growth from
// decay by orders of magnitude, few enough to stay linear (no overflow).
constexpr int kStabilitySteps = 10;
// Duration of the grid-scale damping runs.
constexpr double kGridScaleDuration = 0.1;
// Duration of the growth-rate fits, about 3 periods of the L mode.
constexpr double kGrowthDuration = 20.0;
// Coarse timestep and grid of the convergence runs, and refinement of
// their reference run.
constexpr double kConvergenceDt = 0.004;
constexpr int kConvergenceCells = 32;
constexpr double kReferenceRefinement = 16.0;
// Operations of the closed-form waveFrequency, rounded up.
constexpr double kDispersionOperations = 10.0;
// Largest relative velocity error of linear interpolation of exp(ikx)
// inside a cell, f (1 - f)(k dx)^2/2 at most, divided by (k dx)^2.
constexpr double kInterpolationMismatch = 1.0 / 8.0;

// Macroparticles per cell of a hybrid state (the loader puts the same
// number in every cell).
int particlesPerCell(const HybridState &state) {
    return static_cast<int>(state.particles.size()) / state.grid.nx;
}

// Relative mismatch between the loaded eigenmode and the discrete one of
// the same branch. The loader gives each ion the continuous velocity
// v = -(k/omega) b at its own position, while the discrete eigenmode has
// v/b = -omega_d/(K (1 - omega_d)) (ion equation with the Yee wavenumber
// K), linearly interpolated between nodes. The difference excites the
// other branch, which beats with the mode; its relative amplitude stays
// below this mismatch (measured: 2.6e-4 for R, 6.8e-4 for L, against
// bounds of 1.3e-3 and 3.0e-3; free ion gyration stays near 1e-6).
double loadMismatch(double k, double dx, int ppc, WaveBranch branch) {
    const double omega = waveFrequency(k, branch);
    const double omegaDiscrete = discreteWaveFrequency(k, dx, ppc, branch);
    const double K = yeeWavenumber(k, dx);
    const double loaded = -k / omega;
    const double discrete = -omegaDiscrete / (K * (1.0 - omegaDiscrete));
    return std::fabs(discrete / loaded - 1.0) +
           kInterpolationMismatch * (k * dx) * (k * dx);
}

// Relative error of a frequency measured from the phase accumulated over
// 'periods' periods, when other modes of relative amplitude 'mismatch'
// beat with the measured one: the phase is off by at most 'mismatch' at
// each end.
double beatBias(double mismatch, double periods) {
    return 2.0 * mismatch / (2.0 * kPi * periods);
}

// Relative frequency error of ICN with 2 corrections on an oscillation:
// arg g(i omega dt) = -omega dt (1 - (omega dt)^2/12) + O(dt^5).
double icnFrequencyError(double omega, double dt) {
    return omega * dt * omega * dt / 12.0;
}

// Allowed relative error of a frequency measured over 'periods' periods
// against discreteWaveFrequency, which is exact for the spatial
// discretization: only the time step and the beats remain.
double discreteFrequencyTolerance(double omega, double dt, double mismatch,
                                  double periods) {
    return kSafetyFactor *
           (icnFrequencyError(omega, dt) + beatBias(mismatch, periods));
}

// ICN amplification factor per step of db/dt = -lambda b, z = lambda dt:
// 1 - z + z^2/2 with 1 correction, minus z^3/4 with 2.
std::complex<double> icnFactor(std::complex<double> z, int iterations) {
    const std::complex<double> g = 1.0 - z + 0.5 * z * z;
    return iterations == 1 ? g : g - 0.25 * z * z * z;
}

// Analytic amplitude ratio of a pure whistler (B0 = n = 1) of Yee
// wavenumber kEff after 'steps' ICN steps: it oscillates at kEff^2 and
// decays by hyper-resistivity at nu kEff^4.
double whistlerAmplification(double kEff, double nu, double dt, int iterations,
                             int steps) {
    const double k2 = kEff * kEff;
    const std::complex<double> lambda{nu * k2 * k2, k2};
    return std::pow(std::abs(icnFactor(lambda * dt, iterations)), steps);
}

// Allowed |ln(measured/expected)| of a whistler amplitude after 'steps'
// steps: the analytic factor treats the pure whistler, while the ions
// carry a fraction 1/(k d_i)^2 of its response, which can change ln|b|
// by that much per step.
double ionResponseTolerance(int steps, double kEff) {
    return kSafetyFactor * steps / (kEff * kEff);
}

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
    double duration;
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
    const double duration = steps * scenario.params.dt;
    return WaveMeasurement{-phase / duration,
                           std::abs(previous) / std::abs(initial), duration};
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
        waveScenario(branch, kWaveCells, kWaveDt, 2, 0.01, electrons);
    return measureGrowthRate(scenario, kGrowthDuration);
}

// Relative discretization error of the k d_i = 1 frequency of 'branch' on
// the wave-test grid (4 ions per cell): discrete against continuous.
double discretizationError(WaveBranch branch) {
    const double omega = waveFrequency(1.0, branch);
    return std::fabs(discreteWaveFrequency(1.0, kWaveDx, 4, branch) / omega -
                     1.0);
}

// Allowed error of a growth rate 'gamma' of 'branch' fitted over
// kGrowthDuration at k d_i = 1 with effective resistivity etaPrime. The
// loaded ideal eigenmode misses the discrete dissipative one by a
// relative loadMismatch + eta', which excites other modes with that
// relative amplitude; a least-squares line through a beat of amplitude A
// and frequency W = omega_L - omega_R over a time T has a slope of at
// most 12 A / (W T^2). The discretization changes gamma itself by its
// relative error.
double growthRateTolerance(double gamma, double etaPrime, WaveBranch branch) {
    const double mismatch = loadMismatch(1.0, kWaveDx, 4, branch) + etaPrime;
    const double beatFrequency =
        waveFrequency(1.0, WaveBranch::L) - waveFrequency(1.0, WaveBranch::R);
    const double beatingBias =
        12.0 * mismatch / (beatFrequency * kGrowthDuration * kGrowthDuration);
    return kSafetyFactor *
           (beatingBias + discretizationError(branch) * std::fabs(gamma));
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
        waveScenario(WaveBranch::R, kConvergenceCells, dt, iterations, 0.05);
    const int steps = static_cast<int>(std::lround(2.0 / dt));
    for (int step = 0; step < steps; ++step) {
        hybridStep(scenario.state, scenario.params);
    }
    return scenario.state.B;
}

// Error ratio e(dt)/e(dt/2) against a dt/kReferenceRefinement reference.
double convergenceRatio(int iterations) {
    const VecField reference =
        finalField(kConvergenceDt / kReferenceRefinement, iterations);
    return maxDifference(finalField(kConvergenceDt, iterations), reference) /
           maxDifference(finalField(kConvergenceDt / 2.0, iterations),
                         reference);
}

// Checks that ICN with 'iterations' corrections converges at 'order':
// the next-order terms are smaller by omega dt.
void checkConvergenceOrder(int iterations, int order) {
    const double omegaDt =
        std::fabs(waveFrequency(1.0, WaveBranch::R)) * kConvergenceDt;
    REQUIRE(convergenceRatio(iterations) ==
            Catch::Approx(expectedConvergenceRatio(order, kReferenceRefinement))
                .margin(convergenceRatioTolerance(order, omegaDt)));
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

// Ratio |b_k(t)| / |b_k(0)| of the mode k after 'steps' steps.
double amplitudeRatioAfter(HybridScenario &scenario, double k, int steps) {
    HybridState &state = scenario.state;
    const std::complex<double> initial = modeAmplitude(state.B, state.grid, k);
    for (int step = 0; step < steps; ++step) {
        hybridStep(state, scenario.params);
    }
    return std::abs(modeAmplitude(state.B, state.grid, k)) / std::abs(initial);
}

// Cold uniform plasma in B0 = x on 64 cells whose only perturbation is a
// small magnetic wave of wavenumber k, advanced with 'params'.
HybridScenario perturbedScenario(const HybridParameters &params, double k) {
    HybridScenario scenario =
        waveScenario(WaveBranch::R, kWaveCells, params.dt, params.iterations,
                     0.0, params.electrons);
    HybridState &state = scenario.state;
    VecField B = state.B;
    addCircularPerturbation(B, state.grid, k, 1e-3);
    state = initializeState(state.grid, std::move(B),
                            std::move(state.particles), params.electrons);
    return scenario;
}

// Checks that a grid-scale perturbation (4 cells per wavelength) with
// hyper-resistivity 'nu' keeps the analytic amplitude of the pure
// whistler after kGridScaleDuration.
void checkGridScaleAmplification(double nu) {
    ElectronModel electrons;
    electrons.hyperResistivity = nu;
    HybridScenario scenario =
        perturbedScenario(HybridParameters{kWaveDt, 2, electrons}, kGridScaleK);
    const int steps =
        static_cast<int>(std::lround(kGridScaleDuration / kWaveDt));
    const double kEff = yeeWavenumber(kGridScaleK, kWaveDx);
    const double measured = amplitudeRatioAfter(scenario, kGridScaleK, steps);
    const double expected = whistlerAmplification(kEff, nu, kWaveDt, 2, steps);
    REQUIRE(std::log(measured / expected) ==
            Catch::Approx(0.0).margin(ionResponseTolerance(steps, kEff)));
}

// Checks the amplitude ratio of a Nyquist-scale perturbation after
// kStabilitySteps steps of dt = fraction * timestepLimit (ideal
// electrons) against the analytic ICN amplification, and returns it.
double checkNyquistAmplification(double fraction, int iterations) {
    const Grid1D grid = makeGrid(kWaveCells, 2.0 * kPi);
    const double dt = fraction * timestepLimit(grid, kB0, 1.0, 0.0);
    HybridScenario scenario =
        perturbedScenario(HybridParameters{dt, iterations, {}}, kNyquistK);
    // At 2 cells per wavelength the Yee wavenumber is 2/dx.
    const double kEff = yeeWavenumber(kNyquistK, grid.dx);
    const double measured =
        amplitudeRatioAfter(scenario, kNyquistK, kStabilitySteps);
    const double expected =
        whistlerAmplification(kEff, 0.0, dt, iterations, kStabilitySteps);
    REQUIRE(
        std::log(measured / expected) ==
        Catch::Approx(0.0).margin(ionResponseTolerance(kStabilitySteps, kEff)));
    return measured;
}

// Exact ICN (2 corrections) stability limit of the Nyquist whistler on
// 'grid' with hyper-resistivity nu: the dt at which |g| reaches 1, found
// by bisection down to double precision.
double exactStabilityLimit(const Grid1D &grid, double nu) {
    const double k2 = 4.0 / (grid.dx * grid.dx);
    const std::complex<double> lambda{nu * k2 * k2, k2};
    // |g| <= 1 for small dt, and |g| >= 1 once |z| reaches the largest
    // radius of the stability region.
    double stable = 0.0;
    double unstable = kIcnOscillationRadius / std::abs(lambda);
    for (int i = 0; i < std::numeric_limits<double>::digits; ++i) {
        const double middle = 0.5 * (stable + unstable);
        if (std::abs(icnFactor(lambda * middle, 2)) <= 1.0) {
            stable = middle;
        } else {
            unstable = middle;
        }
    }
    return stable;
}

// Checks that the eigenmode of 'branch' at k d_i = 1, run for 'periods'
// periods, oscillates at the discrete frequency of the scheme and keeps
// its amplitude, and returns the final scenario.
HybridScenario checkEigenmode(WaveBranch branch, double periods) {
    HybridScenario scenario = waveScenario(branch, kWaveCells, kWaveDt);
    const double k = waveNumber(scenario.wave, scenario.state.grid);
    const int ppc = particlesPerCell(scenario.state);
    const double omega = discreteWaveFrequency(k, kWaveDx, ppc, branch);
    const double mismatch = loadMismatch(k, kWaveDx, ppc, branch);
    const WaveMeasurement measured = measureWave(scenario, periods);
    REQUIRE(measured.omega ==
            Catch::Approx(omega).epsilon(
                discreteFrequencyTolerance(omega, kWaveDt, mismatch, periods)));
    // Ideal electrons: no damping. The amplitude only beats with the modes
    // excited by the load mismatch, by at most 2 mismatch.
    REQUIRE(measured.amplitudeRatio ==
            Catch::Approx(1.0).margin(kSafetyFactor * 2.0 * mismatch));
    return scenario;
}

} // namespace

TEST_CASE("The R-mode (whistler) eigenmode has the frequency of the "
          "discrete scheme",
          "[hybrid][wave]") {
    // k d_i = 1: the continuous frequency is -(1 + sqrt 5)/2.
    const double omega = waveFrequency(1.0, WaveBranch::R);
    REQUIRE(omega == Catch::Approx(-(1.0 + std::sqrt(5.0)) / 2.0)
                         .margin(roundoffTolerance(kDispersionOperations,
                                                   std::fabs(omega))));
    const HybridScenario scenario = checkEigenmode(WaveBranch::R, 3.0);
    for (const Particle1D &p : scenario.state.particles) {
        REQUIRE(p.position >= 0.0);
        REQUIRE(p.position < scenario.state.grid.length);
    }
}

TEST_CASE("The L-mode (ion cyclotron) eigenmode has the frequency of the "
          "discrete scheme",
          "[hybrid][wave]") {
    // k d_i = 1: the continuous frequency is (sqrt 5 - 1)/2.
    const double omega = waveFrequency(1.0, WaveBranch::L);
    REQUIRE(omega == Catch::Approx((std::sqrt(5.0) - 1.0) / 2.0)
                         .margin(roundoffTolerance(kDispersionOperations,
                                                   std::fabs(omega))));
    checkEigenmode(WaveBranch::L, 2.0);
}

TEST_CASE("The discrete frequency tends to the analytic one as dx^2",
          "[hybrid][wave][dispersion]") {
    // Pure arithmetic of discreteWaveFrequency: halving dx divides its
    // distance to waveFrequency by 4, up to the next order (k dx)^2.
    for (const WaveBranch branch : {WaveBranch::R, WaveBranch::L}) {
        for (const double k : {0.25, 1.0, 3.0}) {
            CAPTURE(branch == WaveBranch::R ? "R" : "L", k);
            const double omega = waveFrequency(k, branch);
            const double coarse =
                discreteWaveFrequency(k, kWaveDx, 4, branch) - omega;
            const double fine =
                discreteWaveFrequency(k, kWaveDx / 2.0, 4, branch) - omega;
            REQUIRE(
                coarse / fine ==
                Catch::Approx(expectedConvergenceRatio(
                                  2, std::numeric_limits<double>::infinity()))
                    .margin(convergenceRatioTolerance(
                        2, spatialTruncation(k, kWaveDx))));
        }
    }
}

TEST_CASE("ICN is second order with 2 corrections, first order with 1",
          "[hybrid][convergence]") {
    checkConvergenceOrder(2, 2);
    checkConvergenceOrder(1, 1);
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
            Catch::Approx(gammaR).margin(
                growthRateTolerance(gammaR, nu, WaveBranch::R)));
    REQUIRE(hyperResistiveGrowthRate(WaveBranch::L, nu) ==
            Catch::Approx(gammaL).margin(
                growthRateTolerance(gammaL, nu, WaveBranch::L)));
}

TEST_CASE("The timestep limit matches the ICN stability boundary",
          "[hybrid][stability]") {
    const Grid1D grid = makeGrid(kWaveCells, 2.0 * kPi);
    // nu = 0: the classic whistler limit n dx^2 / (2|B0|), which is the
    // exact boundary.
    const double whistlerLimit = grid.dx * grid.dx / 2.0;
    const double roundoff =
        roundoffTolerance(kDispersionOperations, whistlerLimit);
    REQUIRE(timestepLimit(grid, 1.0, 1.0, 0.0) ==
            Catch::Approx(whistlerLimit).margin(roundoff));
    REQUIRE(exactStabilityLimit(grid, 0.0) ==
            Catch::Approx(whistlerLimit).margin(roundoff));
    // With nu > 0 the boundary radius along the ray of nu k^4 + i omega
    // lies between kIcnMinimumRadius and kIcnOscillationRadius, so the
    // estimate (radius kIcnMinimumRadius) stays below it by at most their
    // ratio.
    for (const double nu : {0.002, 0.005, 0.01}) {
        CAPTURE(nu);
        const double limit = timestepLimit(grid, 1.0, 1.0, nu);
        const double exact = exactStabilityLimit(grid, nu);
        REQUIRE(limit < exact);
        REQUIRE(limit > kIcnMinimumRadius / kIcnOscillationRadius * exact);
    }
}

TEST_CASE("Small hyper-resistivity changes long waves only at the "
          "analytic, slow rate",
          "[hybrid][wave][hyper]") {
    for (const double nu : {1e-4, 1e-3}) {
        for (const WaveBranch branch : {WaveBranch::R, WaveBranch::L}) {
            CAPTURE(nu, branch == WaveBranch::R ? "R" : "L");
            ElectronModel electrons;
            electrons.hyperResistivity = nu;
            HybridScenario scenario =
                waveScenario(branch, kWaveCells, kWaveDt, 2, 0.01, electrons);
            const double k = waveNumber(scenario.wave, scenario.state.grid);
            const int ppc = particlesPerCell(scenario.state);
            const double omega = discreteWaveFrequency(k, kWaveDx, ppc, branch);
            const double mismatch = loadMismatch(k, kWaveDx, ppc, branch);
            const WaveMeasurement measured = measureWave(scenario, 2.0);
            const std::complex<double> dissipative =
                dissipativeFrequency(k, nu, branch);
            // nu also shifts the real frequency, by |Re omega(nu) - omega|
            // (second order in eta' = nu k^2), and the mismatch grows by
            // eta' since the loaded eigenmode is the ideal one.
            const double shift =
                std::fabs(dissipative.real() / waveFrequency(k, branch) - 1.0);
            const double etaPrime = nu * k * k;
            REQUIRE(measured.omega ==
                    Catch::Approx(omega).epsilon(
                        discreteFrequencyTolerance(omega, kWaveDt,
                                                   mismatch + etaPrime, 2.0) +
                        kSafetyFactor * shift));
            REQUIRE(
                measured.amplitudeRatio ==
                Catch::Approx(std::exp(dissipative.imag() * measured.duration))
                    .margin(kSafetyFactor * 2.0 * (mismatch + etaPrime)));
        }
    }
}

TEST_CASE("Small hyper-resistivity removes grid-scale whistlers",
          "[hybrid][wave][hyper]") {
    // Ideal electrons: the grid-scale wave survives with ~83 % of its
    // amplitude, the ICN numerical damping |g|^2 = 1 - y^4/4 + y^6/16 per
    // step at y = omega dt ~ 0.42.
    checkGridScaleAmplification(0.0);
    // Same nu as the long-wave test: the rate nu k^4 ~ 43 leaves ~1 %.
    checkGridScaleAmplification(1e-3);
}

TEST_CASE("A grid-scale perturbation decays below the timestep limit and "
          "grows above it",
          "[hybrid][stability]") {
    REQUIRE(checkNyquistAmplification(0.8, 2) < 1.0);
    REQUIRE(checkNyquistAmplification(1.2, 2) > 1.0);
}

TEST_CASE("With one ICN correction a grid-scale perturbation grows even "
          "below the timestep limit",
          "[hybrid][stability]") {
    for (const double fraction : {0.5, 0.8}) {
        CAPTURE(fraction);
        REQUIRE(checkNyquistAmplification(fraction, 1) > 1.0);
    }
}
