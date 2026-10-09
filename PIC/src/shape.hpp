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
inline void depositParticle(Moments &moments, const Particle &p,
                            const Grid1D &grid) {
    const LinearWeights w =
        linearWeights(p.position.x, Centering::Primal, grid);
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
inline Moments depositMoments(const std::vector<Particle> &particles,
                              const Grid1D &grid) {
    Moments moments{makeField(Centering::Primal, grid),
                    makePrimalVecField(grid)};
    for (const Particle &p : particles) {
        depositParticle(moments, p, grid);
    }
    fluxToVelocity(moments);
    return moments;
}

// Interpolates the Yee-mesh E and B at a particle position.
inline EMFields gatherFields(const VecField &E, const VecField &B,
                             const Vector3 &position, const Grid1D &grid) {
    return EMFields{interpolate(E, position.x, grid),
                    interpolate(B, position.x, grid)};
}
