/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
Hybrid field equations on the 1D periodic Yee mesh.

Normalized units: B0 = n0 = m_p = e = mu0 = 1, so lengths are in ion
inertial lengths d_i, times in 1/Omega_ci and velocities in v_A.

  Ampere   j = curl B                      (no displacement current)
  Faraday  dB/dt = -curl E
  Ohm      E = -v_e x B - grad(P_e)/n + eta*j - nu*lap(j),
           v_e = u - j/n,  P_e = n*Te

The Hall term (j x B / n, from -v_e x B) makes short waves whistlers,
omega ~ k^2. At grid scale (a few cells per wavelength) they are poorly
resolved, excited by particle noise, and the fastest modes on the mesh
(they set dt ~ dx^2). Real plasmas damp them through electron physics
this model drops (electron inertia, resonances); the hyper-resistive
term -nu*lap(j) stands in for it: curl of it gives dB/dt = -nu d^4B/dx^4,
a damping rate ~nu k^4. That is huge at grid scale but ~nu at k d_i ~ 1,
where it only weakly damps the R mode (and weakly drives the L mode).
nu is numerical: it must be small enough that long-wave results do not
depend on it.

In 1D, d/dy = d/dz = 0, so curl B = (0, -dBz/dx, dBy/dx) and
-curl E = (0, dEz/dx, -dEy/dx). Bx never changes, hence the discrete
div B = dBx/dx stays at its initial value (zero for a uniform Bx).

Yee staggering makes every derivative a centered two-point difference:
a primal quantity differentiates onto dual nodes and vice versa.
*/
#pragma once

#include <cstddef>
#include <utility> // std::move: hand j, v_e, lap(j) to OhmInputs
#include <vector>

#include "grid.hpp"
#include "shape.hpp"

struct ElectronModel {
    // Isothermal electron temperature Te (P_e = n*Te).
    double temperature = 0.0;
    // Resistivity eta (adds eta*j to E); 0 for an ideal electron fluid.
    double resistivity = 0.0;
    // Hyper-resistivity nu (adds -nu*lap(j) to E); 0 to disable.
    double hyperResistivity = 0.0;
};

// d/dx of a field, landing on the opposite centering: each output node
// sits midway between two input nodes, (right - left)/dx.
inline Field derivativeX(const Field &f, const Grid1D &grid) {
    Field out = makeField(opposite(f.centering), grid);
    for (int i = 0; i < grid.nx; ++i) {
        // Primal -> dual: x_{i+1/2} lies between primal i and i+1.
        // Dual -> primal: x_i lies between dual i-1 and i (x_{i-1/2},
        // x_{i+1/2}).
        const int left = f.centering == Centering::Primal ? i : i - 1;
        out.values[i] = (f.at(left + 1) - f.at(left)) / grid.dx;
    }
    return out;
}

// d^2/dx^2 of a field on its own centering: the three-point stencil
// (f[i+1] - 2 f[i] + f[i-1])/dx^2, built from two Yee derivatives.
inline Field laplacianX(const Field &f, const Grid1D &grid) {
    return derivativeX(derivativeX(f, grid), grid);
}

// Component by component d^2/dx^2 of a vector field.
inline VecField laplacianX(const VecField &f, const Grid1D &grid) {
    return VecField{laplacianX(f.x, grid), laplacianX(f.y, grid),
                    laplacianX(f.z, grid)};
}

// Moves a field onto the opposite centering: each output node sits
// midway between two input nodes and takes their mean, (left + right)/2.
inline Field averageToOpposite(const Field &f, const Grid1D &grid) {
    Field out = makeField(opposite(f.centering), grid);
    for (int i = 0; i < grid.nx; ++i) {
        // Same neighbours as in derivativeX.
        const int left = f.centering == Centering::Primal ? i : i - 1;
        out.values[i] = 0.5 * (f.at(left) + f.at(left + 1));
    }
    return out;
}

// Returns a + scale*b, node by node (same centering assumed).
inline Field addScaled(const Field &a, const Field &b, double scale) {
    Field out = a;
    for (std::size_t i = 0; i < out.values.size(); ++i) {
        out.values[i] += scale * b.values[i];
    }
    return out;
}

// Returns -f, node by node.
inline Field negated(const Field &f) {
    Field out = f;
    for (double &value : out.values) {
        value = -value;
    }
    return out;
}

// Returns (a + b)/2 node by node: the time average of one field at two
// time levels (same centering assumed).
inline Field timeAverage(const Field &a, const Field &b) {
    Field out = a;
    for (std::size_t i = 0; i < out.values.size(); ++i) {
        out.values[i] = 0.5 * (a.values[i] + b.values[i]);
    }
    return out;
}

// Returns (a + b)/2 component by component: the time average of two
// field states, used to center the fields at n+1/2.
inline VecField timeAverage(const VecField &a, const VecField &b) {
    return VecField{timeAverage(a.x, b.x), timeAverage(a.y, b.y),
                    timeAverage(a.z, b.z)};
}

// Ampere's law j = curl B: jx = 0 (dual), jy = -dBz/dx, jz = dBy/dx
// (both primal, as the Yee layout of E).
inline VecField ampere(const VecField &B, const Grid1D &grid) {
    VecField j = makeElectricField(grid);
    j.y = negated(derivativeX(B.z, grid));
    j.z = derivativeX(B.y, grid);
    return j;
}

// Faraday's law, one explicit update from Bn with E held fixed:
// By += dt*dEz/dx, Bz -= dt*dEy/dx (dual), Bx unchanged.
inline VecField faraday(const VecField &Bn, const VecField &E, double dt,
                        const Grid1D &grid) {
    return VecField{Bn.x, addScaled(Bn.y, derivativeX(E.z, grid), dt),
                    addScaled(Bn.z, derivativeX(E.y, grid), -dt)};
}

// Electron velocity v_e = u - j/n on primal nodes (jx = 0 in 1D).
inline VecField electronVelocity(const Moments &moments, const VecField &j) {
    VecField ve = moments.velocity;
    for (std::size_t i = 0; i < ve.y.values.size(); ++i) {
        const double n = safeDensity(moments.density.values[i]);
        ve.y.values[i] -= j.y.values[i] / n;
        ve.z.values[i] -= j.z.values[i] / n;
    }
    return ve;
}

// Grid quantities combined by Ohm's law, all at the same time level.
struct OhmInputs {
    const Moments &moments;
    const VecField &B;
    // Current j = curl B and electron velocity v_e = u - j/n.
    VecField j;
    VecField ve;
    // Laplacian of j, for the hyper-resistive term.
    VecField lapJ;
};

// Dissipative part of Ohm's law at node i: eta*j - nu*lap(j).
inline double dissipation(const Field &j, const Field &lapJ, int i,
                          const ElectronModel &electrons) {
    return electrons.resistivity * j.values[i] -
           electrons.hyperResistivity * lapJ.values[i];
}

// Ex on dual nodes: -(ve_y Bz - ve_z By) - Te*(dn/dx)/n + dissipation, with
// v_e and n averaged primal -> dual and By, Bz used where they live.
inline Field ohmX(const OhmInputs &in, const ElectronModel &electrons,
                  const Grid1D &grid) {
    const Field veY = averageToOpposite(in.ve.y, grid);
    const Field veZ = averageToOpposite(in.ve.z, grid);
    const Field n = averageToOpposite(in.moments.density, grid);
    const Field dndx = derivativeX(in.moments.density, grid);
    Field Ex = makeField(electricCentering(Axis::X, Axis::X), grid);
    for (int i = 0; i < grid.nx; ++i) {
        const double ideal = -(veY.values[i] * in.B.z.values[i] -
                               veZ.values[i] * in.B.y.values[i]);
        const double pressure =
            electrons.temperature * dndx.values[i] / safeDensity(n.values[i]);
        Ex.values[i] =
            ideal - pressure + dissipation(in.j.x, in.lapJ.x, i, electrons);
    }
    return Ex;
}

// Ey on primal nodes: -(ve_z Bx - ve_x Bz) + dissipation, Bz averaged to
// primal.
inline Field ohmY(const OhmInputs &in, const ElectronModel &electrons,
                  const Grid1D &grid) {
    const Field Bz = averageToOpposite(in.B.z, grid);
    Field Ey = makeField(electricCentering(Axis::Y, Axis::X), grid);
    for (int i = 0; i < grid.nx; ++i) {
        const double ideal = -(in.ve.z.values[i] * in.B.x.values[i] -
                               in.ve.x.values[i] * Bz.values[i]);
        Ey.values[i] = ideal + dissipation(in.j.y, in.lapJ.y, i, electrons);
    }
    return Ey;
}

// Ez on primal nodes: -(ve_x By - ve_y Bx) + dissipation, By averaged to
// primal.
inline Field ohmZ(const OhmInputs &in, const ElectronModel &electrons,
                  const Grid1D &grid) {
    const Field By = averageToOpposite(in.B.y, grid);
    Field Ez = makeField(electricCentering(Axis::Z, Axis::X), grid);
    for (int i = 0; i < grid.nx; ++i) {
        const double ideal = -(in.ve.x.values[i] * By.values[i] -
                               in.ve.y.values[i] * in.B.x.values[i]);
        Ez.values[i] = ideal + dissipation(in.j.z, in.lapJ.z, i, electrons);
    }
    return Ez;
}

// Generalized Ohm's law: E from the ion moments and B (j from Ampere).
inline VecField ohm(const Moments &moments, const VecField &B,
                    const ElectronModel &electrons, const Grid1D &grid) {
    VecField j = ampere(B, grid);
    VecField ve = electronVelocity(moments, j);
    VecField lapJ = laplacianX(j, grid);
    const OhmInputs in{moments, B, std::move(j), std::move(ve),
                       std::move(lapJ)};
    return VecField{ohmX(in, electrons, grid), ohmY(in, electrons, grid),
                    ohmZ(in, electrons, grid)};
}
