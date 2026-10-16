/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Catch2 tests of the linear shape function (src/shape.hpp): moment
deposition and field interpolation on the periodic Yee mesh.
*/

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "grid.hpp"
#include "loader.hpp"
#include "shape.hpp"
#include "tolerances.hpp"

namespace {

const Grid1D kGrid = makeGrid(16, 4.0);

// A proton at rest at x with the given weight.
Particle1D protonAt(double x, double weight) {
    return Particle1D{x, Vector3{}, Species{1.0, 1.0}, weight};
}

// Sum of all node values of a field.
double total(const Field &field) {
    double sum = 0.0;
    for (double value : field.values) {
        sum += value;
    }
    return sum;
}

// Field f(x_i) = a + b*x_i on the nodes of the given centering.
Field linearField(double a, double b, Centering centering) {
    Field f = makeField(centering, kGrid);
    for (int i = 0; i < kGrid.nx; ++i) {
        f.values[i] = a + b * nodePosition(i, centering, kGrid);
    }
    return f;
}

} // namespace

TEST_CASE("A uniform load deposits exactly the requested density",
          "[shape][deposit]") {
    const PopulationSettings population{7, 2.5, 0.0, 1};
    const std::vector<Particle1D> particles =
        loadPlasma(population, WaveSettings{}, kGrid);
    const Moments moments = depositMoments(particles, kGrid);
    // Each node sums the contributions of the particles of its two cells.
    const double tolerance =
        roundoffTolerance(2.0 * population.ppc * kShapeOperations, 2.5);
    for (int i = 0; i < kGrid.nx; ++i) {
        REQUIRE(moments.density.values[i] ==
                Catch::Approx(2.5).margin(tolerance));
        REQUIRE(moments.velocity.y.values[i] == 0.0);
    }
}

TEST_CASE("Deposition conserves the total charge across the periodic seam",
          "[shape][deposit]") {
    const std::vector<Particle1D> particles = {
        protonAt(0.0, 1.0), protonAt(1e-12, 0.5),
        protonAt(kGrid.length - 0.01, 2.0), protonAt(1.37, 0.25)};
    const Moments moments = depositMoments(particles, kGrid);
    const double tolerance = roundoffTolerance(
        (particles.size() * kShapeOperations) + kGrid.nx, 3.75);
    REQUIRE(total(moments.density) == Catch::Approx(3.75).margin(tolerance));
}

TEST_CASE("A particle near x = L shares its charge with nodes nx-1 and 0",
          "[shape][deposit]") {
    const double x = kGrid.length - 0.25 * kGrid.dx;
    const Moments moments = depositMoments({protonAt(x, 1.0)}, kGrid);
    const double tolerance = roundoffTolerance(kShapeOperations, x);
    REQUIRE(moments.density.values[kGrid.nx - 1] ==
            Catch::Approx(0.25).margin(tolerance));
    REQUIRE(moments.density.values[0] == Catch::Approx(0.75).margin(tolerance));
}

TEST_CASE("Interpolation is exact for constant and linear fields",
          "[shape][gather]") {
    for (Centering centering : {Centering::Primal, Centering::Dual}) {
        const Field constant = linearField(1.5, 0.0, centering);
        const Field linear = linearField(0.5, 2.0, centering);
        // Interior points only: a linear field is not periodic, so the
        // seam cell would mix both ends of the ramp.
        for (double x : {0.9, 1.234, 2.0, 3.1}) {
            CAPTURE(x);
            const double expected = 0.5 + 2.0 * x;
            REQUIRE(interpolate(constant, x, kGrid) ==
                    Catch::Approx(1.5).margin(
                        roundoffTolerance(kShapeOperations, 1.5)));
            REQUIRE(interpolate(linear, x, kGrid) ==
                    Catch::Approx(expected).margin(
                        roundoffTolerance(kShapeOperations, expected)));
        }
    }
}
