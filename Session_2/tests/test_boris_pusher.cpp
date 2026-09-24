/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
Catch2 validation of the pusher against the analytic references.

Validates the Boris pusher against closed-form solutions of the equations
of motion, rather than just checking that the code runs. See README.md,
section "How we verify correctness", for the reasoning behind each check
and the chosen tolerances.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>

#include "boris_pusher.hpp"
#include "scenarios.hpp"
#include "vector3.hpp"

TEST_CASE("Boris pusher conserves speed exactly in a pure magnetic field, at any dt",
          "[boris][energy-conservation]") {
    // Coarse timestep (8 steps per orbit) on purpose: the magnetic
    // sub-step is an exact rotation of v, so |v| is preserved to
    // floating-point precision independent of how well the orbit itself
    // is resolved.
    const Scenario scenario = cyclotronScenarioCoarseTimestep();
    const Trajectory trajectory = simulateBorisPusher(scenario.particle, scenario.fields, scenario.dt,
                                                        scenario.numSteps);

    const double speed0 = norm(trajectory.velocity.front());
    for (const Vector3 &v : trajectory.velocity) {
        REQUIRE(norm(v) == Catch::Approx(speed0).margin(1e-10));
    }
}

TEST_CASE("Boris pusher trajectory matches the analytic cyclotron circle", "[boris][cyclotron]") {
    const Scenario scenario = cyclotronScenario();
    const Trajectory trajectory = simulateBorisPusher(scenario.particle, scenario.fields, scenario.dt,
                                                        scenario.numSteps);

    for (std::size_t i = 0; i < trajectory.time.size(); ++i) {
        const Vector3 expected =
            analyticCyclotronPosition(trajectory.time[i], scenario.particle.species.charge,
                                      scenario.particle.species.mass, kV0, scenario.fields.B);
        const Vector3 &actual = trajectory.position[i];
        REQUIRE(actual.x == Catch::Approx(expected.x).margin(1e-3));
        REQUIRE(actual.y == Catch::Approx(expected.y).margin(1e-3));
        REQUIRE(actual.z == Catch::Approx(expected.z).margin(1e-12));
    }
}

TEST_CASE("Boris pusher reproduces the analytic E x B drift velocity", "[boris][exb-drift]") {
    const Scenario scenario = exbDriftScenario();
    const Trajectory trajectory = simulateBorisPusher(scenario.particle, scenario.fields, scenario.dt,
                                                        scenario.numSteps);

    // Averaging the velocity over an integer number of full gyro-periods
    // cancels the oscillating gyration component, leaving the guiding
    // center drift velocity.
    Vector3 meanVelocity{0.0, 0.0, 0.0};
    for (const Vector3 &v : trajectory.velocity) {
        meanVelocity = meanVelocity + v;
    }
    meanVelocity = (1.0 / static_cast<double>(trajectory.velocity.size())) * meanVelocity;

    const Vector3 expectedDrift = analyticDriftVelocity(scenario.fields.E, scenario.fields.B);
    REQUIRE(meanVelocity.x == Catch::Approx(expectedDrift.x).margin(1e-3));
    REQUIRE(meanVelocity.y == Catch::Approx(expectedDrift.y).margin(1e-3));
    REQUIRE(meanVelocity.z == Catch::Approx(expectedDrift.z).margin(1e-12));
}
