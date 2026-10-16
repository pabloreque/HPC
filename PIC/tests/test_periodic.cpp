/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Catch2 tests of the periodic boundaries (src/grid.hpp, pusher).
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>

#include "boris_pusher.hpp"
#include "grid.hpp"
#include "scenarios.hpp"

// Distance between two points on a periodic circle of length L.
double periodicDistance(double a, double b, double length) {
    const double d = std::fabs(a - b);
    return std::fmin(d, length - d);
}

TEST_CASE("wrapIndex maps any index onto 0..n-1", "[periodic]") {
    REQUIRE(wrapIndex(0, 8) == 0);
    REQUIRE(wrapIndex(7, 8) == 7);
    REQUIRE(wrapIndex(8, 8) == 0);
    REQUIRE(wrapIndex(17, 8) == 1);
    REQUIRE(wrapIndex(-1, 8) == 7);
    REQUIRE(wrapIndex(-9, 8) == 7);
}

TEST_CASE("wrapPosition maps any position into the box 0 <= x < L",
          "[periodic]") {
    const double L = kPeriodicLength;
    REQUIRE(wrapPosition(1.5, L) == 1.5);
    REQUIRE(wrapPosition(L + 0.25, L) == Catch::Approx(0.25));
    REQUIRE(wrapPosition(-0.25, L) == Catch::Approx(L - 0.25));
    REQUIRE(wrapPosition(L, L) == 0.0);
    REQUIRE(wrapPosition(5.0 * L + 0.1, L) == Catch::Approx(0.1));
    // A tiny negative position must not round to exactly L.
    const double tiny = wrapPosition(-1e-18, L);
    REQUIRE(tiny >= 0.0);
    REQUIRE(tiny < L);
}

TEST_CASE("A periodic trajectory is the unbounded one folded into the box",
          "[periodic][boris]") {
    const Scenario scenario = periodicDriftScenario();
    const Trajectory unbounded = simulateBorisPusher(
        scenario.particle, scenario.fields, scenario.dt, scenario.numSteps);
    const Trajectory periodic =
        simulateBorisPusher(scenario.particle, scenario.fields, scenario.dt,
                            scenario.numSteps, kPeriodicLength);

    // The drift (0.1 per unit time over ~126 time units) must carry the
    // unbounded particle well past the +x boundary.
    REQUIRE(unbounded.position.back().x > 3.0 * kPeriodicLength);
    for (std::size_t i = 0; i < periodic.time.size(); ++i) {
        const Vector3 &p = periodic.position[i];
        const Vector3 &u = unbounded.position[i];
        REQUIRE(p.x >= 0.0);
        REQUIRE(p.x < kPeriodicLength);
        REQUIRE(periodicDistance(p.x, wrapPosition(u.x, kPeriodicLength),
                                 kPeriodicLength) < 1e-9);
        REQUIRE(p.y == u.y);
        REQUIRE(periodic.velocity[i].x == unbounded.velocity[i].x);
    }
}
