# Boris Pusher

Implements and validates the **Boris pusher**, the standard time-centered algorithm for integrating the motion of a charged particle in electric and magnetic fields:

$$\frac{d\vec v}{dt} = \frac{q}{m}\left(\vec E + \vec v \times \vec B\right), \qquad \frac{d\vec x}{dt} = \vec v$$

All quantities use normalized units: $q = m = 1$, $|\vec B| = 1$, so the cyclotron angular frequency $\omega_c = qB/m = 1$ and the Larmor radius $r_L = v_0/\omega_c = v_0$.

## The algorithm

Each step advances $(\vec x, \vec v)$ by `dt` through a symmetric **drift-kick-drift** split:

1. Half position drift: $\vec x_{1/2} = \vec x + \tfrac{1}{2}dt\,\vec v$
2. Half electric kick: $\vec v^- = \vec v + \dfrac{q\,dt}{2m}\vec E$
3. Magnetic rotation (exact, norm-preserving): with $\vec t = \dfrac{q\,dt}{2m}\vec B$, $\vec v' = \vec v^- + \vec v^- \times \vec t$, $\vec s = \dfrac{2\vec t}{1+|\vec t|^2}$, $\vec v^+ = \vec v^- + \vec v' \times \vec s$
4. Half electric kick: $\vec v_{new} = \vec v^+ + \dfrac{q\,dt}{2m}\vec E$
5. Half position drift: $\vec x_{new} = \vec x_{1/2} + \tfrac{1}{2}dt\,\vec v_{new}$

Step 3 is a pure rotation of $\vec v$, so $|\vec v|$ is unchanged by construction, regardless of `dt`. The drift-kick-drift symmetry (steps 1+5 straddling the kick) is what makes the *position* update second-order accurate rather than just the velocity update.

## Simulated scenarios

- **`cyclotron`**: $\vec E = 0$, $\vec B = (0,0,1)$, $\vec v_0 = (1,0,0)$, $\vec x_0 = 0$. Closed-form solution: a circle of radius $r_L=1$ centered at $(0,-1,0)$, $x(t)=\sin(t)$, $y(t)=\cos(t)-1$.
- **`drift`**: same as above plus $\vec E = (0.1,0,0)$. The gyration is unchanged, but the guiding center now drifts at $\vec v_d = (\vec E\times\vec B)/B^2 = (0,-0.1,0)$.

## How we verify correctness

Running the pusher and eyeballing a plausible trajectory is not enough — a subtly wrong implementation can still *look* right. Instead, both the test suite and the notebook check the simulation against independent, closed-form solutions of the same physics:

1. **Speed conservation** (`cyclotron`, deliberately coarse `dt`, 8 steps/orbit): with $\vec E=0$, $|\vec v(t)|$ must equal $|\vec v(0)|$ to floating-point precision, *for any* `dt` — because the rotation step never changes a vector's length. This is checked at coarse resolution specifically to demonstrate it holds independent of how well the orbit itself is resolved. If this drifted, the rotation formula would be broken.
2. **Cyclotron circle** (`cyclotron`, fine `dt`, 1000 steps/orbit): the numeric $(x(t),y(t))$ must track $r_L\sin(\omega_c t)$, $r_L(\cos(\omega_c t)-1)$ within a documented tolerance. This catches sign errors, a wrong radius/frequency, or an inaccurate position update — none of which speed conservation alone would catch (a wrong-radius or wrong-frequency circle still has constant speed).
3. **E×B drift** (`drift`, 10 gyro-periods): averaging the numeric velocity over several full periods cancels the oscillating gyration term and must converge to the analytic $\vec v_d$. This is an independent check that also exercises the electric-kick sub-steps, which the pure-B scenarios never do.

These three checks are implemented as [Catch2](https://github.com/catchorg/Catch2) `TEST_CASE`s in `tests/test_boris_pusher.cpp`, run via `ctest`, and independently visualized in `analysis.ipynb`.

`src/` is a generic pusher: it reads arbitrary initial conditions from a `.txt` file and knows nothing about the test physics. The special-case setups live only in `tests/scenarios.hpp`. The example inputs in `examples/` mirror those scenarios (`cyclotron.txt` ↔ `cyclotronScenario()`, `drift.txt` ↔ `exbDriftScenario()`), so the notebook visualizes exactly what the tests validate — keep the mirrored values in sync if you change either side.

## Project layout

```
Session_2/
├── CMakeLists.txt
├── src/
│   ├── vector3.hpp             # minimal 3D vector type + operators
│   ├── boris_pusher.hpp        # Species, Particle, EMFields, Trajectory, borisStep(), simulateBorisPusher()
│   ├── input.hpp               # .txt parsing: trim, Input, parseInputFile()
│   ├── cli.hpp                 # argv parsing: RunConfig, parseArgs()
│   ├── hdf5_writer.hpp         # HDF5 output: datasets, groups, file handling
│   └── main.cpp                # generic driver: <input.txt> -> pusher -> results.h5 group
├── examples/
│   ├── cyclotron.txt           # input mirroring the cyclotron test scenario
│   └── drift.txt               # input mirroring the E x B drift test scenario
├── tests/
│   ├── scenarios.hpp           # test-only fixtures + analytic reference solutions
│   └── test_boris_pusher.cpp   # Catch2 tests for the three checks above
├── analysis.ipynb              # loads results.h5, plots trajectories against analytic references
└── README.md
```

Input format (`key=value`, `#` = comment): `dt`, `steps`, `E=(Ex Ey Ez)`, `B=(Bx By Bz)`, then one particle after another, each with its own `q=`, `m=`, `x=(x y z)`, `v=(vx vy vz)` lines (a `q=`/`m=` line applies to the next particle, so each particle can be a different species). Assumes uniform static `E`/`B` shared by all particles.

## Build and run (C++)

```bash
cd Session_2
mkdir build && cd build
cmake ..
cmake --build .
```

The first configure fetches [Catch2](https://github.com/catchorg/Catch2) v3.5.4 via CMake `FetchContent` (requires internet access).

Run the simulation (one input file = one run, one group; run twice to regenerate both notebook datasets):

```bash
./session_2_prequeijo ../examples/cyclotron.txt results.h5 /cyclotron
./session_2_prequeijo ../examples/drift.txt results.h5 /drift
```

Full usage: `./session_2_prequeijo <input.txt> [output.h5] [group]` (defaults: `results.h5`, `/trajectory`). If the output file exists the group is replaced, otherwise the file is created. With more than one particle in the input, each trajectory is written to `group/particle_i`.

This prints a per-particle summary (`|v0|`, `|v_final|`, mean velocity — no analytic comparison; that lives in the tests) and writes `build/results.h5`, an HDF5 file with two groups:

- `/cyclotron` — `t`, `x`, `y`, `z`, `vx`, `vy`, `vz` for the pure-B scenario
- `/drift` — same datasets for the E×B scenario

Run the tests:

```bash
ctest --output-on-failure
# or directly:
./session_2_prequeijo_test
```

## Analysis notebook

`analysis.ipynb` reads `build/results.h5` with `h5py` and plots, for each scenario, the numeric trajectory/speed/velocity against the closed-form analytic reference — visually reproducing the same three checks the test suite verifies numerically.

```bash
cd Session_2
jupyter notebook analysis.ipynb
```

Run both `./session_2_prequeijo` commands above first so `build/results.h5` exists before opening the notebook.

## Expected result

- `|v(t)|` for the `cyclotron` scenario is flat at 1.0 to floating-point precision (max deviation ~1e-15).
- The numeric `cyclotron` trajectory overlays the analytic circle almost exactly (max position error ~6e-5 over 3 full orbits at 1000 steps/orbit).
- The running mean of $v_y$ for the `drift` scenario converges to the analytic drift velocity $-0.1$ (final mean ≈ -0.09996) as the gyration term averages out over successive periods.
