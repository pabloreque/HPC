# Finite Difference Order of Accuracy

Numerically demonstrates that:

- `(f[i+1] - f[i]) / dx` is **first-order accurate** (error ~ O(dx))
- `(f[i+1] - f[i-1]) / (2·dx)` is **second-order accurate** (error ~ O(dx²))

It uses `f(x) = sin(x)` (exact derivative `cos(x)`) evaluated at `x0 = 1.0`,
with a sequence of step sizes `dx` that is halved on every iteration. The
observed order is estimated as `p = log2(error(dx) / error(dx/2))`, since
halving `dx` means an error that behaves like `C·dx^p` produces a ratio of
successive errors equal to `2^p`.

## Project layout

```
Session_1/
├── CMakeLists.txt
├── src/
│   └── main.cpp        # computes forward/central diff errors, writes results.h5
├── analysis.ipynb       # loads results.h5, plots epsilon vs dx, discusses results
└── README.md
```

## Build and run (C++)

```bash
cd Session_1
mkdir build && cd build
cmake ..
cmake --build .
./fd_order
```

This prints a convergence table to the console and writes `build/results.h5`,
an HDF5 file with three datasets at the root:

- `dx` — step sizes used
- `err_forward` — absolute error of the forward difference formula
- `err_central` — absolute error of the central difference formula

## Analysis notebook

`analysis.ipynb` reads `build/results.h5` with `h5py`, plots the absolute
error ε against the step size Δx on a log-log scale (ε on the y-axis, Δx on
the x-axis), overlays reference lines of slope 1 and 2, and computes the
observed order `p = log2(ε(Δx)/ε(Δx/2))` for both formulas.

```bash
cd Session_1
jupyter notebook analysis.ipynb
```

Run the C++ program first so `build/results.h5` exists before opening the
notebook.

## Expected result

`p_forward` converges to ~1.0 and `p_central` to ~2.0 as `dx` decreases, until
floating-point round-off error starts to dominate for very small `dx` values
(expected behavior, not a bug).

## Questions I have

- Should I always prefix standard library names with `std::` (as I do
  throughout `main.cpp`), or is `using namespace std;` (or per-name
  `using std::vector;` declarations) acceptable/preferred in course code?
