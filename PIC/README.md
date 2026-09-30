# Boris Pusher

Integrates the motion of N charged particles, each with its own species (charge, mass), in electric and magnetic fields:

$$\frac{d\vec v}{dt} = \frac{q}{m}\left(\vec E + \vec v \times \vec B\right), \qquad \frac{d\vec x}{dt} = \vec v$$

with a time-centered drift-kick-drift Boris step (see `src/boris_pusher.hpp` for the sequence). The bundled examples use normalized units: $q = m = 1$, $|\vec B| = 1$, so $\omega_c = 1$ and $r_L = v_0$; the pusher itself accepts any species.

Two scenarios (see `examples/` and `tests/scenarios.hpp`):

- **`cyclotron`**: $\vec E = 0$ → circular gyration.
- **`drift`**: $\vec E = (0.1,0,0)$ → gyration plus E×B drift.

Correctness is checked three ways (speed conservation, cyclotron circle, E×B drift) by the Catch2 suite and visualized in the notebook, plus per-particle checks for fixed-seed mixed species.

## Project layout

```
Session_2/
├── CMakeLists.txt
├── src/
│   ├── vector3.hpp             # 3D vector type + operators
│   ├── boris_pusher.hpp        # particle pusher (Species, Particle, EMFields, ...)
│   ├── input.hpp               # .txt parsing
│   ├── cli.hpp                 # argv parsing
│   ├── hdf5_writer.hpp         # HDF5 output
│   └── main.cpp                # generic driver
├── examples/
│   ├── cyclotron.txt           # pure-B input
│   ├── drift.txt               # E×B input
│   ├── multi_cyclotron.txt     # pure-B input, 4 fixed-seed mixed species
│   └── multi_drift.txt         # E×B input, same 4 particles
├── tests/
│   ├── scenarios.hpp           # test fixtures + analytic references
│   └── test_boris_pusher.cpp   # Catch2 tests
├── analysis.ipynb              # loads results.h5, plots vs analytic references
└── README.md
```

Input format (`key=value`, `#` = comment): `dt`, `steps`, `E`, `B`, then per-particle `q=`, `m=`, `x=`, `v=` lines. Uniform static `E`/`B` shared by all particles.

## Build and run (C++)

```bash
cd Session_2
mkdir build && cd build
cmake ..
cmake --build .
```

The first configure fetches [Catch2](https://github.com/catchorg/Catch2) v3.5.4 via CMake `FetchContent` (requires internet access).

```bash
./session_2_prequeijo ../examples/cyclotron.txt results.h5 /cyclotron
./session_2_prequeijo ../examples/drift.txt results.h5 /drift
./session_2_prequeijo ../examples/multi_cyclotron.txt results.h5 /multi_cyclotron
./session_2_prequeijo ../examples/multi_drift.txt results.h5 /multi_drift
```

This prints a per-particle summary and writes `build/results.h5` with four groups (`/cyclotron`, `/drift`, `/multi_cyclotron`, `/multi_drift`), each with `t`, `x`, `y`, `z`, `vx`, `vy`, `vz` (per `particle_N` subgroup for the multi runs).

Full usage: `./session_2_prequeijo <input.txt> [output.h5] [group]` (defaults `results.h5`, `/trajectory`).

Run the tests:

```bash
ctest --output-on-failure
# or directly:
./session_2_prequeijo_test
```

## Analysis notebook

`analysis.ipynb` reads `build/results.h5` with `h5py` and plots each scenario against its analytic reference (including per-particle checks for the mixed-species runs).

```bash
cd Session_2
jupyter notebook analysis.ipynb
```

Run the C++ program first so `build/results.h5` exists before opening the notebook.

## Expected result

- `|v(t)|` flat at 1.0 to floating-point precision (~1e-15).
- Numeric trajectory overlays the analytic circle (max error ~6e-5).
- Mean `vy` converges to the drift velocity −0.1 (≈ −0.09996).
- Mixed species: per-particle `|v|` flat (~1e-15), per-particle trajectories within ~1e-4 of their own analytic curves, common drift for all.

## Next steps

- Yee mesh for the electric and magnetic fields (TODO).
