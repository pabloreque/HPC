// Physical setup shared by the simulation (main.cpp, which exports results
// for plotting) and the test suite (test_boris_pusher.cpp, which checks
// those same results against closed-form solutions). Keeping the scenario
// parameters here, in one place, means the tests can never silently drift
// from what was actually simulated and plotted.
#pragma once

#include <cmath>
#include <string>

#include "boris_pusher.hpp"
#include "vector3.hpp"

constexpr double kPi = 3.14159265358979323846;

// Normalized units: charge = mass = 1, |B| = 1, so the cyclotron angular
// frequency omega_c = q*B/m = 1 and the Larmor radius r_L = v0/omega_c = v0.
constexpr double kCharge = 1.0;
constexpr double kMass = 1.0;
constexpr double kB0 = 1.0;
constexpr double kV0 = 1.0;
constexpr double kE0 = 0.1;

struct Scenario {
    std::string name;
    double charge;
    double mass;
    Vector3 E;
    Vector3 B;
    ParticleState initialState;
    double dt;
    int numSteps;
};

inline double cyclotronAngularFrequency(double charge, double mass, const Vector3 &B) {
    return charge * norm(B) / mass;
}

inline double cyclotronPeriod(double charge, double mass, const Vector3 &B) {
    return 2.0 * kPi / cyclotronAngularFrequency(charge, mass, B);
}

// Fine time resolution (200 steps per orbit), run for 3 orbits: used both to
// export a smooth trajectory for plotting and to check that the numeric
// trajectory tracks the analytic circle closely.
inline Scenario cyclotronScenario() {
    const Vector3 B{0.0, 0.0, kB0};
    const double period = cyclotronPeriod(kCharge, kMass, B);
    // 1000 steps/period keeps the (second-order) discretization error many
    // orders of magnitude below the 1e-3 test margin, even after several
    // full orbits.
    const int stepsPerPeriod = 1000;
    const int numPeriods = 3;
    return Scenario{"cyclotron",
                     kCharge,
                     kMass,
                     Vector3{0.0, 0.0, 0.0},
                     B,
                     ParticleState{Vector3{0.0, 0.0, 0.0}, Vector3{kV0, 0.0, 0.0}},
                     period / stepsPerPeriod,
                     stepsPerPeriod * numPeriods};
}

// Deliberately coarse time resolution (8 steps per orbit) to demonstrate
// that the Boris rotation conserves |v| exactly even when the trajectory
// itself is poorly resolved -- speed conservation does not depend on dt.
inline Scenario cyclotronScenarioCoarseTimestep() {
    const Vector3 B{0.0, 0.0, kB0};
    const double period = cyclotronPeriod(kCharge, kMass, B);
    const int stepsPerPeriod = 8;
    const int numPeriods = 5;
    return Scenario{"cyclotron_coarse",
                     kCharge,
                     kMass,
                     Vector3{0.0, 0.0, 0.0},
                     B,
                     ParticleState{Vector3{0.0, 0.0, 0.0}, Vector3{kV0, 0.0, 0.0}},
                     period / stepsPerPeriod,
                     stepsPerPeriod * numPeriods};
}

// Uniform E perpendicular to B on top of the same cyclotron setup: the
// gyration is unchanged, but the guiding center now drifts at v_d = ExB/B^2.
inline Scenario exbDriftScenario() {
    const Vector3 B{0.0, 0.0, kB0};
    const Vector3 E{kE0, 0.0, 0.0};
    const double period = cyclotronPeriod(kCharge, kMass, B);
    const int stepsPerPeriod = 200;
    const int numPeriods = 10;
    return Scenario{"exb_drift",
                     kCharge,
                     kMass,
                     E,
                     B,
                     ParticleState{Vector3{0.0, 0.0, 0.0}, Vector3{kV0, 0.0, 0.0}},
                     period / stepsPerPeriod,
                     stepsPerPeriod * numPeriods};
}

// Closed-form solution for the cyclotron scenario's initial conditions
// (x0 = 0, v0 = (v0, 0, 0), B = (0, 0, B0)): a circle of radius v0/omega_c
// centered at the guiding center (0, -v0/omega_c, 0).
inline Vector3 analyticCyclotronPosition(double t, double charge, double mass, double v0, const Vector3 &B) {
    const double omega = cyclotronAngularFrequency(charge, mass, B);
    const double r = v0 / omega;
    return Vector3{r * std::sin(omega * t), r * (std::cos(omega * t) - 1.0), 0.0};
}

// Guiding-center E x B drift velocity, independent of charge and mass.
inline Vector3 analyticDriftVelocity(const Vector3 &E, const Vector3 &B) {
    return (1.0 / dot(B, B)) * cross(E, B);
}
