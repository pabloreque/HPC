/*
Pablo Requeijo, October 1 2026.
Session 3, HPC coursework: hybrid PIC.
1D periodic Yee mesh: periodic wrapping, grid geometry, centering and
field storage.

Boundaries are periodic along x: a particle that leaves the box [0, L)
through +x re-enters through -x (and vice versa); grid indices wrap the
same way.

Primal nodes sit at x_i = i*dx and dual nodes at x_{i+1/2} = (i+1/2)*dx,
i = 0 ... nx-1 (arrays only have integer indices, so dual index i stands
for x_{i+1/2}). Both are periodic: node nx is node 0.

Yee layout: each field component has its own centering along each
direction. E (and j) components are dual along their own axis and primal
along the others; B components are the reverse. In 3D this puts E on the
cell edges and B on the face centers:

  Ex (i+1/2, j, k)      Bx (i, j+1/2, k+1/2)
  Ey (i, j+1/2, k)      By (i+1/2, j, k+1/2)
  Ez (i, j, k+1/2)      Bz (i+1/2, j+1/2, k)

This 1D code only uses the x direction (first index): Ex, By, Bz are dual,
Ey, Ez, Bx are primal. The centering functions below take (component,
direction), so a 2D/3D version only has to query y and z as well.
*/
#pragma once

#include <cmath>
#include <cstddef> // std::size_t: node count of the field arrays
#include <stdexcept>
#include <vector>

// Maps any integer index onto [0, count), wrapping negative values too.
inline int wrapIndex(int i, int count) {
    const int remainder = i % count;
    return remainder < 0 ? remainder + count : remainder;
}

// Maps any position onto the periodic interval [0, length).
inline double wrapPosition(double x, double length) {
    const double wrapped = x - length * std::floor(x / length);
    // A tiny negative x (e.g. -1e-17) rounds to exactly 'length' after the
    // shift, outside [0, length) (cell index nx); it is the same point as 0
    // on the periodic circle.
    return wrapped < length ? wrapped : 0.0;
}

enum class Centering { Primal, Dual };

// Cartesian axes, used both as vector components and as directions.
enum class Axis { X, Y, Z };

// Whether E (or j) 'component' sits on primal or dual nodes along
// 'direction': dual if the component points along that direction, primal
// otherwise. In 1D (direction x): Ex dual, Ey and Ez primal.
constexpr Centering electricCentering(Axis component, Axis direction) {
    return component == direction ? Centering::Dual : Centering::Primal;
}

// Whether B 'component' sits on primal or dual nodes along 'direction':
// the reverse of E, primal if the component points along that direction,
// dual otherwise. In 1D (direction x): Bx primal, By and Bz dual.
constexpr Centering magneticCentering(Axis component, Axis direction) {
    return component == direction ? Centering::Primal : Centering::Dual;
}

// Returns the other centering (primal <-> dual).
constexpr Centering opposite(Centering centering) {
    return centering == Centering::Primal ? Centering::Dual : Centering::Primal;
}

struct Grid1D {
    int nx = 0;
    double length = 0.0;
    double dx = 0.0;
};

// Builds a periodic grid of nx cells covering [0, length).
inline Grid1D makeGrid(int nx, double length) {
    if (nx <= 0 || length <= 0.0) {
        throw std::invalid_argument("grid needs nx > 0 and length > 0");
    }
    return Grid1D{nx, length, length / nx};
}

// Offset of a centering in units of dx: 0 for primal, 1/2 for dual.
constexpr double centeringOffset(Centering centering) {
    return centering == Centering::Dual ? 0.5 : 0.0;
}

// Position of node i with the given centering.
inline double nodePosition(int i, Centering centering, const Grid1D &grid) {
    return (i + centeringOffset(centering)) * grid.dx;
}

// One scalar quantity on the grid: nx values plus where they live.
struct Field {
    Centering centering = Centering::Primal;
    std::vector<double> values;

    // Value at a periodic node index (any integer is accepted).
    double at(int i) const {
        return values[wrapIndex(i, static_cast<int>(values.size()))];
    }
};

// Three components of a vector quantity, each with its own centering.
struct VecField {
    Field x;
    Field y;
    Field z;
};

// Allocates a field of zeros with the given centering.
inline Field makeField(Centering centering, const Grid1D &grid) {
    return Field{centering,
                 std::vector<double>(static_cast<std::size_t>(grid.nx))};
}

// Allocates E (or j) with its Yee centering along x, filled with zeros.
inline VecField makeElectricField(const Grid1D &grid) {
    return VecField{makeField(electricCentering(Axis::X, Axis::X), grid),
                    makeField(electricCentering(Axis::Y, Axis::X), grid),
                    makeField(electricCentering(Axis::Z, Axis::X), grid)};
}

// Allocates B with its Yee centering along x, filled with zeros.
inline VecField makeMagneticField(const Grid1D &grid) {
    return VecField{makeField(magneticCentering(Axis::X, Axis::X), grid),
                    makeField(magneticCentering(Axis::Y, Axis::X), grid),
                    makeField(magneticCentering(Axis::Z, Axis::X), grid)};
}

// Allocates a vector field with all components on primal nodes.
inline VecField makePrimalVecField(const Grid1D &grid) {
    return VecField{makeField(Centering::Primal, grid),
                    makeField(Centering::Primal, grid),
                    makeField(Centering::Primal, grid)};
}
