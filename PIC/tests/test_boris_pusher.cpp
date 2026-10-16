/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
Catch2 validation of the pusher against the analytic references.

Validates the Boris pusher against closed-form solutions of the equations
of motion, rather than just checking that the code runs. TESTS.md gives
the physics behind each check and the origin of each tolerance.

The gyration and drift tests use several species (including a negative
charge) in fields of magnitude 2, so that no factor of q, m or |B|
equals 1 and a wrong formula or rotation sense cannot pass by
coincidence.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "boris_pusher.hpp"
#include "scenarios.hpp"
#include "tolerances.hpp"
#include "vector3.hpp"

namespace {

// Floating-point operations of one Boris step (two electric half kicks,
// the magnetic rotation with its two cross products, two half drifts),
// rounded up: the round-off budget of each step.
constexpr double kBorisStepOperations = 100.0;

// Relative error of the Boris gyrofrequency: the rotation turns v by
// 2 atan(omega dt / 2) per step instead of omega dt, too slowly by
// (omega dt)^2 / 12.
double borisFrequencyError(double omegaDt) {
    return omegaDt * omegaDt / 12.0;
}

// Round-off budget of a quantity of magnitude 'scale' carried through
// 'steps' Boris steps.
double borisRoundoff(int steps, double scale) {
    return roundoffTolerance(steps * kBorisStepOperations, scale);
}

// Velocity averaged over the samples of whole gyro-periods (the last
// sample, which repeats the first gyrophase, is left out), which cancels
// the gyration.
Vector3 periodAveragedVelocity(const Trajectory &trajectory) {
    const std::size_t count = trajectory.velocity.size() - 1;
    Vector3 sum{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < count; ++i) {
        sum = sum + trajectory.velocity[i];
    }
    return (1.0 / static_cast<double>(count)) * sum;
}

// Error of the period-averaged drift velocity. The Boris gyration lags
// by (omega dt)^2/12 of a turn per turn, so the run ends that fraction of
// a turn short and the gyration velocity, of size 'gyrationSpeed', does
// not average out completely.
double driftTolerance(double omegaDt, double gyrationSpeed) {
    return kSafetyFactor * borisFrequencyError(omegaDt) * gyrationSpeed;
}

// Unit vector along v.
Vector3 unitVector(const Vector3 &v) {
    return (1.0 / norm(v)) * v;
}

// Component of v perpendicular to the unit vector 'axis'.
Vector3 perpendicularPart(const Vector3 &v, const Vector3 &axis) {
    return v - dot(v, axis) * axis;
}

} // namespace

TEST_CASE(
    "Boris pusher conserves speed exactly in a pure magnetic field, at any dt",
    "[boris][energy-conservation]") {
    // Coarse timestep (8 steps per orbit) on purpose: the magnetic
    // sub-step is an exact rotation of v, so |v| is preserved to
    // floating-point precision independent of how well the orbit itself
    // is resolved.
    const Scenario scenario = cyclotronScenarioCoarseTimestep();
    const Trajectory trajectory = simulateBorisPusher(
        scenario.particle, scenario.fields, scenario.dt, scenario.numSteps);

    const double speed0 = norm(trajectory.velocity.front());
    for (const Vector3 &v : trajectory.velocity) {
        REQUIRE(norm(v) == Catch::Approx(speed0).margin(
                               borisRoundoff(scenario.numSteps, speed0)));
    }
}

TEST_CASE("Boris pusher trajectory matches the analytic cyclotron circle",
          "[boris][cyclotron]") {
    const Scenario scenario = cyclotronScenario();
    const Trajectory trajectory = simulateBorisPusher(
        scenario.particle, scenario.fields, scenario.dt, scenario.numSteps);
    const double omega = cyclotronAngularFrequency(
        scenario.particle.species.charge, scenario.particle.species.mass,
        scenario.fields.B);
    const double omegaDt = omega * scenario.dt;
    // The orbit stays on the exact circle but lags in phase by
    // omega dt * (omega dt)^2/12 per step: after N steps the position is
    // off by r_L times the accumulated lag.
    const double larmorRadius = kV0 / omega;
    const double tolerance = kSafetyFactor * larmorRadius * scenario.numSteps *
                             omegaDt * borisFrequencyError(omegaDt);

    for (std::size_t i = 0; i < trajectory.time.size(); ++i) {
        const Vector3 expected = analyticCyclotronPosition(
            trajectory.time[i], scenario.particle.species.charge,
            scenario.particle.species.mass, kV0, scenario.fields.B);
        const Vector3 &actual = trajectory.position[i];
        REQUIRE(actual.x == Catch::Approx(expected.x).margin(tolerance));
        REQUIRE(actual.y == Catch::Approx(expected.y).margin(tolerance));
        // No force along B = z and vz = 0: z stays exactly 0.
        REQUIRE(actual.z == 0.0);
    }
}

TEST_CASE("Boris pusher reproduces the analytic E x B drift velocity",
          "[boris][exb-drift]") {
    const Scenario scenario = exbDriftScenario();
    const Trajectory trajectory = simulateBorisPusher(
        scenario.particle, scenario.fields, scenario.dt, scenario.numSteps);

    const Vector3 meanVelocity = periodAveragedVelocity(trajectory);
    const Vector3 expectedDrift =
        analyticDriftVelocity(scenario.fields.E, scenario.fields.B);
    const double omega = cyclotronAngularFrequency(
        scenario.particle.species.charge, scenario.particle.species.mass,
        scenario.fields.B);
    const Vector3 gyrationVelocity = scenario.particle.velocity - expectedDrift;
    const double tolerance =
        driftTolerance(omega * scenario.dt, norm(gyrationVelocity));
    REQUIRE(meanVelocity.x == Catch::Approx(expectedDrift.x).margin(tolerance));
    REQUIRE(meanVelocity.y == Catch::Approx(expectedDrift.y).margin(tolerance));
    // No force along B = z and vz = 0: the motion stays in the plane.
    REQUIRE(meanVelocity.z == 0.0);
}

namespace {

// Oblique uniform field, |B| = 2, and an initial velocity with both
// parallel and perpendicular components (helical orbit).
const Vector3 kObliqueB{1.2, 0.0, 1.6};
const Vector3 kHelixVelocity{1.0, 0.5, -0.3};
constexpr int kGyrationStepsPerPeriod = 1000;
constexpr int kGyrationPeriods = 3;
constexpr int kGyrationSteps = kGyrationStepsPerPeriod * kGyrationPeriods;

struct GyrationCase {
    const char *name;
    Species species;
};

const std::vector<GyrationCase> kGyrationCases = {
    {"proton", Species{1.0, 1.0}},
    {"alpha", Species{2.0, 4.0}},
    {"negative light charge", Species{-1.0, 0.05}}};

// Signed rate at which v_perp rotates about B: dv/dt = (q/m) v x B is a
// rotation with angular velocity -(q/m) B, so ions (q > 0) turn
// clockwise seen from the tip of B and negative charges anticlockwise.
double signedGyrofrequency(const Species &species, const Vector3 &B) {
    return -species.charge * norm(B) / species.mass;
}

// Larmor radius m |v_perp| / (|q| |B|).
double larmorRadius(const Species &species, const Vector3 &v,
                    const Vector3 &B) {
    return species.mass * norm(perpendicularPart(v, unitVector(B))) /
           (std::fabs(species.charge) * norm(B));
}

// Guiding center of a particle in a uniform B: r + (m / (q B^2)) v x B.
Vector3 guidingCenter(const Particle3D &particle, const Vector3 &B) {
    const double factor =
        particle.species.mass / (particle.species.charge * dot(B, B));
    return particle.position + factor * cross(particle.velocity, B);
}

// Analytic guiding center at time t of the helical orbit of 'species'
// started at the origin with kHelixVelocity: the initial guiding center
// carried along B at constant v_parallel.
Vector3 expectedGuidingCenter(const Species &species, double t) {
    const Vector3 b = unitVector(kObliqueB);
    const Vector3 initial = guidingCenter(
        Particle3D{Vector3{}, kHelixVelocity, species}, kObliqueB);
    return initial + (dot(kHelixVelocity, b) * t) * b;
}

// Helical orbit of 'species' from the origin in kObliqueB, resolved with
// kGyrationStepsPerPeriod steps per gyro-period.
Trajectory gyrationTrajectory(const Species &species) {
    const double period =
        2.0 * kPi / std::fabs(signedGyrofrequency(species, kObliqueB));
    return simulateBorisPusher(Particle3D{Vector3{}, kHelixVelocity, species},
                               EMFields{Vector3{}, kObliqueB},
                               period / kGyrationStepsPerPeriod,
                               kGyrationSteps);
}

// v_perp as a complex number in the right-handed basis (e1, e2, b) with
// e1 along the initial v_perp, so that its argument is the gyrophase.
std::complex<double> perpendicularPhasor(const Vector3 &v, const Vector3 &e1,
                                         const Vector3 &e2) {
    return {dot(v, e1), dot(v, e2)};
}

// Signed gyrofrequency measured from the gyrophase accumulated over the
// whole trajectory (positive = anticlockwise about B).
double measuredGyrofrequency(const Trajectory &trajectory, const Vector3 &B) {
    const Vector3 b = unitVector(B);
    const Vector3 e1 =
        unitVector(perpendicularPart(trajectory.velocity.front(), b));
    const Vector3 e2 = cross(b, e1);
    double phase = 0.0;
    for (std::size_t i = 1; i < trajectory.velocity.size(); ++i) {
        phase +=
            std::arg(perpendicularPhasor(trajectory.velocity[i], e1, e2) /
                     perpendicularPhasor(trajectory.velocity[i - 1], e1, e2));
    }
    return phase / trajectory.time.back();
}

} // namespace

TEST_CASE("Gyration has the analytic frequency qB/m and rotation sense",
          "[boris][gyration]") {
    const double omegaDt = 2.0 * kPi / kGyrationStepsPerPeriod;
    for (const GyrationCase &gyration : kGyrationCases) {
        CAPTURE(gyration.name);
        const Trajectory trajectory = gyrationTrajectory(gyration.species);
        REQUIRE(measuredGyrofrequency(trajectory, kObliqueB) ==
                Catch::Approx(signedGyrofrequency(gyration.species, kObliqueB))
                    .epsilon(kSafetyFactor * borisFrequencyError(omegaDt)));
    }
}

// The Boris step conserves the guiding center (up to its motion along B)
// and |v_perp| exactly, so the Larmor radius, v_parallel and the
// orbit's angular momentum below are invariants: only round-off moves
// them.
TEST_CASE("Gyration keeps the Larmor radius about a guiding center that "
          "moves along B at constant v_parallel",
          "[boris][gyration]") {
    const Vector3 b = unitVector(kObliqueB);
    const double vParallel = dot(kHelixVelocity, b);
    for (const GyrationCase &gyration : kGyrationCases) {
        CAPTURE(gyration.name);
        const Trajectory trajectory = gyrationTrajectory(gyration.species);
        const double radius =
            larmorRadius(gyration.species, kHelixVelocity, kObliqueB);
        for (std::size_t i = 0; i < trajectory.time.size(); ++i) {
            const Vector3 offset =
                trajectory.position[i] -
                expectedGuidingCenter(gyration.species, trajectory.time[i]);
            REQUIRE(norm(offset) == Catch::Approx(radius).margin(
                                        borisRoundoff(kGyrationSteps, radius)));
            REQUIRE(dot(trajectory.velocity[i], b) ==
                    Catch::Approx(vParallel).margin(
                        borisRoundoff(kGyrationSteps, norm(kHelixVelocity))));
        }
    }
}

TEST_CASE("Gyration is diamagnetic: the orbit's angular momentum opposes "
          "B for ions and follows B for negative charges",
          "[boris][gyration]") {
    const Vector3 b = unitVector(kObliqueB);
    for (const GyrationCase &gyration : kGyrationCases) {
        CAPTURE(gyration.name);
        const Trajectory trajectory = gyrationTrajectory(gyration.species);
        const double radius =
            larmorRadius(gyration.species, kHelixVelocity, kObliqueB);
        const double vPerp = norm(perpendicularPart(kHelixVelocity, b));
        // L_par = ((r - r_gc) x v) . b = -sign(q) r_L v_perp about the
        // actual guiding center: the current loop of the orbit always
        // weakens B.
        const double expected =
            -std::copysign(radius * vPerp, gyration.species.charge);
        for (std::size_t i = 0; i < trajectory.time.size(); ++i) {
            const Vector3 offset =
                trajectory.position[i] -
                expectedGuidingCenter(gyration.species, trajectory.time[i]);
            REQUIRE(dot(cross(offset, trajectory.velocity[i]), b) ==
                    Catch::Approx(expected).margin(
                        borisRoundoff(kGyrationSteps, radius * vPerp)));
        }
    }
}

namespace {

// E x B drift setups with |E| = 0.1 and |B| = 2 (drift speed 0.05). The
// expected drift is written out by hand rather than computed with
// cross(), so that a wrong component of the cross product cannot hide
// in both the pusher and the reference.
struct DriftCase {
    const char *name;
    Vector3 E;
    Vector3 B;
    Vector3 expectedDrift;
};

const std::vector<DriftCase> kDriftCases = {
    {"E_x B_y", {0.1, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 0.05}},
    {"E_y B_x", {0.0, 0.1, 0.0}, {2.0, 0.0, 0.0}, {0.0, 0.0, -0.05}},
    {"E_y B_z", {0.0, 0.1, 0.0}, {0.0, 0.0, 2.0}, {0.05, 0.0, 0.0}},
    {"E_z B_y", {0.0, 0.0, 0.1}, {0.0, 2.0, 0.0}, {-0.05, 0.0, 0.0}},
    {"E_z B_x", {0.0, 0.0, 0.1}, {2.0, 0.0, 0.0}, {0.0, 0.05, 0.0}},
    {"E_x B_z", {0.1, 0.0, 0.0}, {0.0, 0.0, 2.0}, {0.0, -0.05, 0.0}},
    // Every component of E and B mixed; E . B = 0.
    {"oblique", {0.16, 0.1, -0.12}, kObliqueB, {0.04, -0.1, -0.03}}};

// Initial velocity with components along every axis, so that each case
// has both a gyration and a motion along B.
const Vector3 kDriftInitialVelocity{0.2, -0.1, 0.3};
constexpr int kDriftStepsPerPeriod = 1000;
constexpr int kDriftPeriods = 5;

// Orbit of 'species' from the origin in the fields of 'drift', over an
// integer number of gyro-periods.
Trajectory driftTrajectory(const DriftCase &drift, const Species &species) {
    const double period =
        2.0 * kPi / std::fabs(signedGyrofrequency(species, drift.B));
    return simulateBorisPusher(
        Particle3D{Vector3{}, kDriftInitialVelocity, species},
        EMFields{drift.E, drift.B}, period / kDriftStepsPerPeriod,
        kDriftStepsPerPeriod * kDriftPeriods);
}

} // namespace

TEST_CASE("E x B drift has the analytic velocity for every orientation of "
          "E and B, independently of charge and mass",
          "[boris][exb-drift]") {
    const double omegaDt = 2.0 * kPi / kDriftStepsPerPeriod;
    for (const DriftCase &drift : kDriftCases) {
        const Vector3 b = unitVector(drift.B);
        const double vParallel = dot(kDriftInitialVelocity, b);
        // Initial velocity in the frame drifting with the guiding center.
        const Vector3 gyrationVelocity =
            perpendicularPart(kDriftInitialVelocity, b) - drift.expectedDrift;
        const double tolerance =
            driftTolerance(omegaDt, norm(gyrationVelocity));
        for (const GyrationCase &gyration : kGyrationCases) {
            CAPTURE(drift.name, gyration.name);
            const Trajectory trajectory =
                driftTrajectory(drift, gyration.species);
            const Vector3 error = periodAveragedVelocity(trajectory) -
                                  vParallel * b - drift.expectedDrift;
            REQUIRE(norm(error) == Catch::Approx(0.0).margin(tolerance));
            REQUIRE(dot(trajectory.velocity.back(), b) ==
                    Catch::Approx(vParallel).margin(
                        borisRoundoff(kDriftStepsPerPeriod * kDriftPeriods,
                                      norm(kDriftInitialVelocity))));
        }
    }
}
