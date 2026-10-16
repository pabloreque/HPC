/*
Pablo Requeijo, October 16 2026.
HPC coursework: PIC.
Shared error model of the Catch2 tests (see TESTS.md).

Every tolerance of the suite is derived from the error it allows, never
fitted to a measured value:

- Round-off: the checked quantity is exact in exact arithmetic (an
  invariant or a closed-form identity), so only rounding can move it.
- Truncation: the leading discretization error of the scheme, times
  kSafetyFactor, which covers the neglected higher-order terms.
- Noise: kNoiseSigmas standard deviations of a statistical estimate.
- Convergence: the error ratio between steps h and h/2 is 2^order, up to
  the next-order terms.
*/
#pragma once

#include <cmath>
#include <limits>

// Relative rounding error of one double-precision operation.
constexpr double kMachineEpsilon = std::numeric_limits<double>::epsilon();

// Margin over a leading-order error estimate: the next-order terms are
// smaller by a further power of the small parameter, so they cannot
// double the estimate.
constexpr double kSafetyFactor = 2.0;

// Allowed deviation of a statistical estimate, in standard deviations: a
// Gaussian exceeds it with probability 6e-7, so a correct code never fails
// by bad luck.
constexpr double kNoiseSigmas = 5.0;

// Floating-point operations of one particle's deposit or one field
// interpolation with the linear shape (cell index, two shape weights, two
// weighted additions), rounded up.
constexpr double kShapeOperations = 10.0;

// Worst-case round-off of a result of magnitude 'scale' computed with
// 'operations' floating-point operations, each rounding by at most one
// machine epsilon (no cancellation, errors added linearly).
inline double roundoffTolerance(double operations, double scale = 1.0) {
    return operations * kMachineEpsilon * scale;
}

// Relative leading truncation error of the second-order spatial operators
// acting on a mode of wavenumber k with node spacing dx. Their Taylor
// coefficients (centered derivative 1/24, two-point average 1/8, linear
// shape smoothing 1/6) are all at most 1/6.
inline double spatialTruncation(double k, double dx) {
    return (k * dx) * (k * dx) / 6.0;
}

// Expected error ratio e(h)/e(h/2) of a method of the given order, when
// the error is measured against a reference computed at h/refinement
// (refinement = infinity for an analytic reference): the reference's own
// error, h^order/refinement^order, is subtracted from both errors.
inline double expectedConvergenceRatio(int order, double refinement) {
    const double p = order;
    return std::pow(2.0, p) * (1.0 - std::pow(refinement, -p)) /
           (1.0 - std::pow(refinement / 2.0, -p));
}

// Margin on that ratio: next-order terms of relative size 'nextOrder' at
// the coarse step change it by about 2^order * nextOrder.
inline double convergenceRatioTolerance(int order, double nextOrder) {
    return kSafetyFactor * std::pow(2.0, order) * nextOrder;
}
