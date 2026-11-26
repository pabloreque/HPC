/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
First-order (linear) shape function: particle <-> grid transfers.

A particle at x_p between nodes x_i and x_{i+1} contributes
(1 - (x_p - x_i)/dx) of its weight to node i and (x_p - x_i)/dx to node
i+1. The same weights interpolate grid fields back to the particle
from whichever centering (primal or dual) each field component lives on.
Indices wrap periodically.

Moments: in normalized units (e = 1), quasi-neutrality n_e = sum_s q_s n_s,
so the deposited ion charge density is the electron density n used by
Ohm's law, and the ion bulk velocity is the charge-weighted mean
u = (sum q w v S) / (sum q w S).
*/
#pragma once

#include <algorithm> // std::max: density floor in safeDensity
#include <cmath>
#include <cstddef>
#include <vector>

#include "boris_pusher.hpp"
#include "grid.hpp"
#include "vector3.hpp"

// Densities below this are clamped before dividing by them, so an empty
// cell does not produce inf/NaN (a sign of too few particles per cell).
constexpr double kDensityFloor = 1e-10;

// Clamps a density to kDensityFloor before it is used as a divisor.
inline double safeDensity(double n) {
    return std::max(n, kDensityFloor);
}

// The two nodes surrounding a position and their linear weights.
struct LinearWeights {
    int leftNode;
    int rightNode;
    double leftWeight;
    double rightWeight;
};

// Linear weights of position x on nodes of the given centering.
inline LinearWeights linearWeights(double x, Centering centering,
                                   const Grid1D &grid) {
    // Position in units of node spacing, counted from the first node of
    // this centering: 2.3 means 30 % of the way from node 2 to node 3.
    const double fractionalIndex = x / grid.dx - centeringOffset(centering);
    const double cell = std::floor(fractionalIndex);
    const double fraction = fractionalIndex - cell;
    const int i = static_cast<int>(cell);
    return LinearWeights{wrapIndex(i, grid.nx), wrapIndex(i + 1, grid.nx),
                         1.0 - fraction, fraction};
}

// Interpolates a grid field at position x (from the field's centering).
inline double interpolate(const Field &field, double x, const Grid1D &grid) {
    const LinearWeights w = linearWeights(x, field.centering, grid);
    return w.leftWeight * field.values[w.leftNode] +
           w.rightWeight * field.values[w.rightNode];
}

// Interpolates every component of a vector field at position x.
inline Vector3 interpolate(const VecField &field, double x,
                           const Grid1D &grid) {
    return Vector3{interpolate(field.x, x, grid), interpolate(field.y, x, grid),
                   interpolate(field.z, x, grid)};
}

// Adds 'amount' to the two nodes of 'w', split by the linear weights.
inline void scatter(Field &field, const LinearWeights &w, double amount) {
    field.values[w.leftNode] += w.leftWeight * amount;
    field.values[w.rightNode] += w.rightWeight * amount;
}

// Ion moments on the grid, the inputs of Ohm's law.
struct Moments {
    // Electron density n = ion charge density (quasi-neutrality), primal.
    Field density;
    // Ion bulk velocity u, all components primal.
    VecField velocity;
};

// Deposits one particle's charge and charge flux onto primal nodes.
inline void depositParticle(Moments &moments, const Particle1D &p,
                            const Grid1D &grid) {
    const LinearWeights w = linearWeights(p.position, Centering::Primal, grid);
    const double charge = p.species.charge * p.weight;
    scatter(moments.density, w, charge);
    scatter(moments.velocity.x, w, charge * p.velocity.x);
    scatter(moments.velocity.y, w, charge * p.velocity.y);
    scatter(moments.velocity.z, w, charge * p.velocity.z);
}

// Turns the deposited flux n*u into the bulk velocity u, node by node.
inline void fluxToVelocity(Moments &moments) {
    for (std::size_t i = 0; i < moments.density.values.size(); ++i) {
        const double n = safeDensity(moments.density.values[i]);
        moments.velocity.x.values[i] /= n;
        moments.velocity.y.values[i] /= n;
        moments.velocity.z.values[i] /= n;
    }
}

// Computes density and bulk velocity on primal nodes from all particles.
inline Moments depositMoments(const std::vector<Particle1D> &particles,
                              const Grid1D &grid) {
    Moments moments{makeField(Centering::Primal, grid),
                    makePrimalVecField(grid)};
    for (const Particle1D &p : particles) {
        depositParticle(moments, p, grid);
    }
    fluxToVelocity(moments);
    return moments;
}

// Mass-weighted sums of one population on the primal nodes, from which
// the ion temperature follows.
struct SecondMomentSums {
    // Sum of w S (number density).
    Field number;
    // Sum of w m S.
    Field mass;
    // Sum of w m v S.
    VecField momentum;
    // Sum of w m |v|^2 S.
    Field energy;
};

// Deposits one particle into the second-moment sums.
inline void depositSecondMoment(SecondMomentSums &sums, const Particle1D &p,
                                const Grid1D &grid) {
    const LinearWeights w = linearWeights(p.position, Centering::Primal, grid);
    const double massWeight = p.weight * p.species.mass;
    scatter(sums.number, w, p.weight);
    scatter(sums.mass, w, massWeight);
    scatter(sums.momentum.x, w, massWeight * p.velocity.x);
    scatter(sums.momentum.y, w, massWeight * p.velocity.y);
    scatter(sums.momentum.z, w, massWeight * p.velocity.z);
    scatter(sums.energy, w, massWeight * dot(p.velocity, p.velocity));
}

// Temperature T = (K - |P|^2/M) / (3 N) at node i from the sums
// K = w m v^2, P = w m v, M = w m, N = w: the kinetic energy left after
// removing the bulk flow, per particle and per degree of freedom.
inline double nodeTemperature(const SecondMomentSums &sums, std::size_t i) {
    const Vector3 momentum{sums.momentum.x.values[i], sums.momentum.y.values[i],
                           sums.momentum.z.values[i]};
    const double thermalEnergy =
        sums.energy.values[i] -
        dot(momentum, momentum) / safeDensity(sums.mass.values[i]);
    return thermalEnergy / (3.0 * safeDensity(sums.number.values[i]));
}

// Isotropic ion temperature T = <m |v - u|^2> / 3 on the primal nodes,
// with u the mass-weighted bulk velocity. A Maxwellian of per-component
// std dev vth gives T = m vth^2. Diagnostic only: Ohm's law does not need
// it, so it is computed at output times, not every step.
inline Field depositTemperature(const std::vector<Particle1D> &particles,
                                const Grid1D &grid) {
    SecondMomentSums sums{
        makeField(Centering::Primal, grid), makeField(Centering::Primal, grid),
        makePrimalVecField(grid), makeField(Centering::Primal, grid)};
    for (const Particle1D &p : particles) {
        depositSecondMoment(sums, p, grid);
    }
    Field temperature = makeField(Centering::Primal, grid);
    for (std::size_t i = 0; i < temperature.values.size(); ++i) {
        temperature.values[i] = nodeTemperature(sums, i);
    }
    return temperature;
}

// Mean of f (1 - f) over the fractional positions f = (q + 1/2)/ppc,
// q = 0 ... ppc-1, of 'particlesPerCell' evenly spaced particles in a
// cell: 1/6 + 1/(12 ppc^2), and 1/6 for a continuum.
inline double meanShapeOverlap(int particlesPerCell) {
    const double ppc = particlesPerCell;
    return 1.0 / 6.0 + 1.0 / (12.0 * ppc * ppc);
}

// Fraction sigma by which a grid mode exp(ikx) shrinks when it is
// interpolated to evenly spaced particles and deposited back with the
// linear shape. A particle at fraction f of a cell sees (1 - f) E_i +
// f E_{i+1} and returns the same weights, so node i gets back
// 1 - 2 f (1 - f)(1 - cos k dx) of the mode: sigma = 4 sin^2(k dx/2)
// <f (1 - f)>. Exact for particles that do not move along x.
inline double shapeResponseDeficit(double k, double dx, int particlesPerCell) {
    const double halfPhase = std::sin(0.5 * k * dx);
    return 4.0 * halfPhase * halfPhase * meanShapeOverlap(particlesPerCell);
}

// Interpolates the Yee-mesh E and B at a particle position x.
inline EMFields gatherFields(const VecField &E, const VecField &B, double x,
                             const Grid1D &grid) {
    return EMFields{interpolate(E, x, grid), interpolate(B, x, grid)};
}
