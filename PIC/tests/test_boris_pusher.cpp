/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
Catch2 validation of the pusher against the analytic references.

Validates the Boris pusher against closed-form solutions of the equations
of motion, rather than just checking that the code runs. See README.md,
section "How we verify correctness", for the reasoning behind each check
and the chosen tolerances.

The gyration tests use several species (including a negative charge) in
an oblique field of magnitude 2, so that no factor of q, m or |B| equals
1 and a wrong formula or rotation sense cannot pass by coincidence.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "boris_pusher.hpp"
#include "scenarios.hpp"
#include "vector3.hpp"

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
        REQUIRE(norm(v) == Catch::Approx(speed0).margin(1e-10));
    }
}

TEST_CASE("Boris pusher trajectory matches the analytic cyclotron circle",
          "[boris][cyclotron]") {
    const Scenario scenario = cyclotronScenario();
    const Trajectory trajectory = simulateBorisPusher(
        scenario.particle, scenario.fields, scenario.dt, scenario.numSteps);

    for (std::size_t i = 0; i < trajectory.time.size(); ++i) {
        const Vector3 expected = analyticCyclotronPosition(
            trajectory.time[i], scenario.particle.species.charge,
            scenario.particle.species.mass, kV0, scenario.fields.B);
        const Vector3 &actual = trajectory.position[i];
        REQUIRE(actual.x == Catch::Approx(expected.x).margin(1e-3));
        REQUIRE(actual.y == Catch::Approx(expected.y).margin(1e-3));
        REQUIRE(actual.z == Catch::Approx(expected.z).margin(1e-12));
    }
}

TEST_CASE("Boris pusher reproduces the analytic E x B drift velocity",
          "[boris][exb-drift]") {
    const Scenario scenario = exbDriftScenario();
    const Trajectory trajectory = simulateBorisPusher(
        scenario.particle, scenario.fields, scenario.dt, scenario.numSteps);

    // Averaging the velocity over an integer number of full gyro-periods
    // cancels the oscillating gyration component, leaving the guiding
    // center drift velocity.
    Vector3 meanVelocity{0.0, 0.0, 0.0};
    for (const Vector3 &v : trajectory.velocity) {
        meanVelocity = meanVelocity + v;
    }
    meanVelocity =
        (1.0 / static_cast<double>(trajectory.velocity.size())) * meanVelocity;

    const Vector3 expectedDrift =
        analyticDriftVelocity(scenario.fields.E, scenario.fields.B);
    REQUIRE(meanVelocity.x == Catch::Approx(expectedDrift.x).margin(1e-3));
    REQUIRE(meanVelocity.y == Catch::Approx(expectedDrift.y).margin(1e-3));
    REQUIRE(meanVelocity.z == Catch::Approx(expectedDrift.z).margin(1e-12));
}

namespace {

// Oblique uniform field, |B| = 2, and an initial velocity with both
// parallel and perpendicular components (helical orbit).
const Vector3 kObliqueB{1.2, 0.0, 1.6};
const Vector3 kHelixVelocity{1.0, 0.5, -0.3};
constexpr int kGyrationStepsPerPeriod = 1000;
constexpr int kGyrationPeriods = 3;

struct GyrationCase {
    const char *name;
    Species species;
};

const std::vector<GyrationCase> kGyrationCases = {
    {"proton", Species{1.0, 1.0}},
    {"alpha", Species{2.0, 4.0}},
    {"negative light charge", Species{-1.0, 0.05}}};

// Unit vector along v.
Vector3 unitVector(const Vector3 &v) {
    return (1.0 / norm(v)) * v;
}

// Component of v perpendicular to the unit vector 'axis'.
Vector3 perpendicularPart(const Vector3 &v, const Vector3 &axis) {
    return v - dot(v, axis) * axis;
}

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
                               kGyrationStepsPerPeriod * kGyrationPeriods);
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
    for (const GyrationCase &gyration : kGyrationCases) {
        CAPTURE(gyration.name);
        const Trajectory trajectory = gyrationTrajectory(gyration.species);
        REQUIRE(measuredGyrofrequency(trajectory, kObliqueB) ==
                Catch::Approx(signedGyrofrequency(gyration.species, kObliqueB))
                    .epsilon(1e-5));
    }
}

TEST_CASE("Gyration keeps the Larmor radius about a guiding center that "
          "moves along B at constant v_parallel",
          "[boris][gyration]") {
    const Vector3 b = unitVector(kObliqueB);
    for (const GyrationCase &gyration : kGyrationCases) {
        CAPTURE(gyration.name);
        const Trajectory trajectory = gyrationTrajectory(gyration.species);
        const double radius =
            larmorRadius(gyration.species, kHelixVelocity, kObliqueB);
        const double vParallel = dot(kHelixVelocity, b);
        for (std::size_t i = 0; i < trajectory.time.size(); ++i) {
            const Vector3 offset =
                trajectory.position[i] -
                expectedGuidingCenter(gyration.species, trajectory.time[i]);
            REQUIRE(norm(offset) == Catch::Approx(radius).epsilon(1e-4));
            REQUIRE(dot(trajectory.velocity[i], b) ==
                    Catch::Approx(vParallel).margin(1e-12));
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
                    Catch::Approx(expected).epsilon(1e-4));
        }
    }
}
