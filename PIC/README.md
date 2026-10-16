# PIC

Cumulative particle-in-cell project, one HPC session after another:

- **Session 2**: Boris pusher for test particles in uniform static fields.
- **Session 3**: 1D3V hybrid PIC (1D in space, 3D velocities and field components) with self-consistent fields (Ohm's law with resistivity and hyper-resistivity + Faraday on a periodic Yee mesh, iterated Crank–Nicolson), and periodic boundaries along x.

## Session 2: Boris pusher

Integrates the motion of N charged particles, each with its own species (charge, mass), in electric and magnetic fields:

$$\frac{d\vec v}{dt} = \frac{q}{m}\left(\vec E + \vec v \times \vec B\right), \qquad \frac{d\vec x}{dt} = \vec v$$

with a time-centered drift-kick-drift Boris step (see `src/boris_pusher.hpp` for the sequence). The bundled examples use normalized units: $q = m = 1$, $|\vec B| = 1$, so $\omega_c = 1$ and $r_L = v_0$; the pusher itself accepts any species. This mode is fully 3D: positions and velocities have three components, and since the fields are uniform no grid is needed.

Two scenarios (see `examples/` and `tests/scenarios.hpp`):

- **`cyclotron`**: $\vec E = 0$ → circular gyration.
- **`drift`**: $\vec E = (0.1,0,0)$ → gyration plus E×B drift.

Correctness is checked three ways (speed conservation, cyclotron circle, E×B drift) by the Catch2 suite and visualized in the notebook, plus per-particle checks for fixed-seed mixed species.

## Session 3: hybrid PIC and periodic boundaries

**Dimensionality: 1D3V.** The hybrid model is one-dimensional in space but keeps all three components of velocities and fields:

| | Dimensions |
|---|---|
| Grid (`Grid1D`) | 1: $n_x$ nodes along $x$, no $y$/$z$ nodes |
| $\vec E$, $\vec B$, $\vec j$, $n$, $\vec u$ | depend on $x$ only; each vector keeps its 3 components (3 arrays of $n_x$ values) |
| Particle velocity | 3 components |
| Particle position | only $x$ matters |

- The plasma is assumed uniform in $y$ and $z$ ($\partial_y = \partial_z = 0$): each plane $x = $ const is an infinite slab with the same values everywhere. Hence only $x$ derivatives appear, $j_x = 0$ and $B_x$ never changes.
- Velocities must stay 3D: with $\vec B_0 = B_0\hat x$, the $\vec v\times\vec B$ force and the ion gyration live in the $y$–$z$ plane, and so do the transverse waves ($B_y$, $B_z$, $j_y$, $j_z$).
- Particle positions are still `Vector3` (the Session 2 pusher is reused) and the pusher advances $y$ and $z$, but deposit, gather and periodic wrapping only use $x$: the $y$, $z$ coordinates have no effect.
- This is exact for waves propagating along $\vec B_0 = B_0\hat x$ (plane wavefronts), which is what this session studies. It cannot describe oblique propagation, 2D/3D structures (reconnection, turbulence) or instabilities growing across $x$.

**Periodic boundaries.** A particle leaving the box $[0, L)$ through $+x$ re-enters through $-x$ (`src/grid.hpp`), and grid indices wrap the same way. In test-particle mode an optional `length=` line turns this on for the Session 2 pusher (`examples/periodic_drift.txt`).

**Hybrid model** (course slides 130–138; slide numbers refer to `docs/Nicolas_Sessions.pdf`). Ions are macroparticles. Electrons are a massless isothermal fluid. Normalized units: $B_0 = n_0 = m_p = e = \mu_0 = 1$, so lengths are in $d_i$, times in $\Omega_{ci}^{-1}$ and velocities in $v_A$:

$$\vec j = \nabla\times\vec B, \qquad \frac{\partial \vec B}{\partial t} = -\nabla\times\vec E, \qquad \vec E = -\left(\vec u - \frac{\vec j}{n}\right)\times\vec B - \frac{T_e \nabla n}{n} + \eta\,\vec j - \nu\,\nabla^2\vec j$$

$\vec E$ is not evolved in time; it is recomputed from Ohm's law, using the ion moments described below.

**Hyper-resistivity** (`nu`, the $-\nu\nabla^2\nabla\times\vec B$ term of the course's ICN slide, slide 167). In Faraday's law it becomes $-\nu\nabla^4\vec B$, a dissipation rate $\sim\nu k^4$ that removes grid-scale whistlers while barely touching long wavelengths ($\eta$ acts as $\eta k^2$). In 1D, $\nabla^2\vec j = \partial_x^2\vec j$ is the three-point stencil obtained by applying the Yee derivative twice (primal → dual → primal). Two caveats:

- In this Ohm's law, $\eta$ and $\nu$ act on $\vec E$, which the ions also feel, so they do not only dissipate. For cold parallel modes with $\eta' = \eta + \nu k^2$, the dispersion relation (derived below) becomes $\omega^2 + a\omega - a = 0$ with $a = k^2(1 + i\eta')$. At $k d_i = 1$ the R mode is damped ($\gamma_R = -1.17\,\eta'$) but the L mode *grows* slowly ($\gamma_L = +0.17\,\eta'$): the energy budget picks up $\eta'\!\int n\,\vec j\cdot\vec u$, which is positive when the ions move faster than the current, as in the L mode. Both rates are measured to within 1% (R) and 6% (L) in the Catch2 `[hyper]` test.
- The term is explicit, so it does not allow a larger `dt`: it adds its own stability limit (see below). Its purpose is to clean grid-scale noise, e.g. in `thermal_noise`.

**Moments of the distribution function on the grid** (slides 136–137, `src/shape.hpp`). Ohm's law needs the fluid quantities $n$ and $\vec u$ on the grid, but the ions are a set of macroparticles that sample $f(\vec r, \vec v)$. Each particle *deposits* its contribution on the surrounding nodes through a shape function $S$:

$$n_i = \sum_p q_p\, w_p\, S(x_p - x_i), \qquad \vec u_i = \frac{1}{n_i}\sum_p q_p\, w_p\, \vec v_p\, S(x_p - x_i)$$

- $w_p$ is the macroparticle **statistical weight**: how much it contributes to the local density. A uniform load of `ppc` particles per cell uses $w_p = $ `density/ppc`.
- $S$ is the order-1 B-spline, i.e. linear interpolation. A particle between $x_i$ and $x_{i+1}$ gives $\left(1 - \frac{x_p - x_i}{\Delta x}\right) w_p$ to node $i$ and $\frac{x_p - x_i}{\Delta x}\, w_p$ to node $i+1$.
- The charge $q_p$ makes $n$ the ion *charge* density. By quasi-neutrality that is the electron density in Ohm's law. For protons ($q = 1$) it reduces to the plain number density.
- $n$ and $\vec u$ live on primal nodes, and the indices wrap periodically: a particle near $x = L$ also deposits on node 0. An empty cell is clamped to a tiny density floor before dividing.
- The same weights do the reverse operation, interpolating $\vec E$ and $\vec B$ from the grid to each particle at its mid-step position $\vec r^{n+1/2}$ (slide 147). Each field component is interpolated from its own primal or dual centering.

Tested in `tests/test_shape.cpp`: a uniform load gives exactly $n = n_0$ on every node, the total charge is conserved across the periodic seam, and linear fields are interpolated exactly. The moments are saved as `n`, `ux`, `uy`, `uz` in `results.h5`.

**Yee mesh** (slide 157, `src/grid.hpp`): 1D in $x$, $n_x$ cells, periodic. $E_x$, $B_y$, $B_z$ live on dual nodes $x_{i+1/2}$; $E_y$, $E_z$, $B_x$, $n$, $\vec u$ live on primal nodes $x_i$. Every derivative is then a centered two-point difference (second order), and $\partial_t B_x = 0$ keeps the discrete $\nabla\cdot\vec B$ exactly zero. Centering is defined per component *and per direction* (`electricCentering`, `magneticCentering`), so a 2D/3D extension only has to add the $y$/$z$ directions. (Slide 151 writes $\partial_t B_z = +\partial_x E_y - \partial_y E_x$; the curl gives the opposite sign, which is what the code uses.)

**Time scheme** (`src/hybrid_step.hpp`): iterated Crank–Nicolson predictor-corrector, as in PHARE:

```
predict:  B* = Faraday(B^n, E^n),  E* = Ohm(moments^n, B*)
repeat `iterations` times:
    Ebar = (E^n + E*)/2,  Bbar = (B^n + B*)/2           (fields at n+1/2)
    particles* = Boris(particles^n, Ebar, Bbar), x wrapped periodically
    moments*   = deposit(particles*)                     (moments at n+1)
    B* = Faraday(B^n, Ebar),  E* = Ohm(moments*, B*)
commit
```

Measured convergence of $B$ at fixed $t$ against a $dt/16$ reference (Catch2 `[convergence]`):

| corrections | error ratio $e(dt)/e(dt/2)$ | order |
|---|---|---|
| 1 | 2.0 | 1 |
| 2 (default) | 4.00 | 2 |

Stability: the fastest grid mode, $k = 2/dx$, oscillates as a whistler ($\omega \approx k^2 |B_0|/n$) and decays by hyper-resistivity (rate $\nu k^4$). With 2 corrections the ICN amplification factor is $g(z) = 1 - z + z^2/2 - z^3/4$, $z = dt\,(\nu k^4 + i\omega)$, and $|g| \le 1$ holds up to $|z| = 2$ on the imaginary axis but only $|z| \approx 1.77$ in between. The driver therefore warns when $dt > r/|\nu k^4 + i\omega|$, with $r = 2$ for $\nu = 0$ (i.e. $dt < n\,dx^2/(2|B_0|)$) and $r = 1.75$ otherwise (`timestepLimit`, within 20% of the exact boundary). With $dx = 2\pi/64$: for $\nu = 0$ the runs blow up between $dt = 0.004$ and $0.005$, as predicted (limit $0.0048$); with $\nu = 0.005$ the limit drops to $0.0018$, and the `alfven_R` setup ($dt = 0.002$) blows up once `nu=0.005` is added. A naive pure-diffusion estimate $dt < dx^4/(8\nu) = 0.0023$ would miss this.

Why 2 corrections: each correction pushes and deposits every particle again, so with many particles per cell the cost of a step is roughly proportional to `iterations`, and a third correction costs about 50% more. The extra cost buys little, because the corrections converge to the Crank–Nicolson solution, which already has an $O(dt^2)$ error. Two corrections reach that order, and more do not raise it. For a pure oscillation of the field equation alone ($z = i y$, $y = \omega\,dt$, as for a grid-scale whistler), the corrections change $|g|^2$ as follows:

| corrections | $\lvert g\rvert^2$ | |
|---|---|---|
| 1 | $1 + y^4/4$ | unstable: every wave grows |
| 2 (default) | $1 - y^4/4 + y^6/16$ | stable for $y \le 2$, damping $\sim y^4$ |
| 3 | $1 - y^6/16 + y^8/64$ | stable for $y \le 2$, damping $\sim y^6$ |
| $\infty$ (CN) | $1$ | no damping, for $y < 2$ |

The limit $y \le 2$ does not move: it is the convergence condition $|z/2| < 1$ of the fixed-point iteration, and beyond it more corrections only make things worse. The only gain is less numerical damping at grid scale, where hyper-resistivity removes the waves anyway; long waves ($y \ll 1$) barely feel either term. To improve long-wave accuracy, reduce $dt$ (error $\propto dt^2$) rather than add corrections. (The rows for 1 and 3 corrections are computed analytically, not measured.)

**Eigenmodes: R and L waves** (`src/loader.hpp`). The validation waves are not in the course slides; they are the standard parallel electromagnetic waves (R and L waves in F. F. Chen, *Introduction to Plasma Physics and Controlled Fusion*, ch. 4), in the limit of cold ions and massless electrons. Derivation from the hybrid equations above:

1. *Equilibrium*: $n = 1$, $\vec u = 0$, $\vec E = 0$, $\vec B_0 = \hat x$. Nothing moves.
2. *Perturbation*: transverse fields written as complex numbers, $b = B_y + iB_z$, $w = u_y + iu_z$, $e = E_y + iE_z$, all $\propto e^{i(kx - \omega t)}$. The motion is transverse, so the density stays uniform and the pressure term drops out.
3. *Complex identities*: for a transverse vector $\vec v$, $\vec v\times\hat x \leftrightarrow -i\,v$; in 1D, $\nabla\times$ acts as $i\partial_x \to -k$. Hence $j = -k\,b$ (Ampère) and $\partial_t b = k\,e$ (Faraday).
4. *Ohm* (linearized, $\vec v_e = \vec u - \vec j$): $e = -(\vec u - \vec j)\times\hat x \to i(w - j) = i(w + kb)$.
5. *Ion momentum*: $\partial_t w = e - iw = ikb$, so $-i\omega w = ikb$ and
$$w = -\frac{k}{\omega}\,b$$
This is the ion velocity loaded with the wave, so that only one branch is excited.
6. *Faraday*: $-i\omega b = k\,e = i k w + i k^2 b = i k^2 b\left(1 - \tfrac{1}{\omega}\right)$, which gives the dispersion relation
$$\omega^2 + k^2\omega - k^2 = 0 \quad\Rightarrow\quad \omega = \frac{k^2}{2}\left(-1 \pm \sqrt{1 + \frac{4}{k^2}}\right)$$

The $+$ root ($\omega > 0$) is the **L** wave: it rotates with the ions and saturates at the ion cyclotron frequency, $\omega \to \Omega_{ci} = 1$ as $k\to\infty$. The $-$ root ($\omega < 0$) is the **R** wave or whistler: it rotates with the electrons and grows as $\omega \approx -k^2$, because the electron resonance $\Omega_{ce}$ is infinite for massless electrons. For $k d_i \ll 1$ both reduce to the Alfvén wave, $|\omega| = k$. At $k d_i = 1$, $\omega_L = (\sqrt5 - 1)/2 = 0.618$ and $\omega_R = -(\sqrt5 + 1)/2 = -1.618$. Adding a dissipative term $\eta' j$ to Ohm's law (step 4) replaces $k^2$ by $a = k^2(1 + i\eta')$, which gives the hyper-resistive rates above.

**Validation** (Catch2 + `analysis_session_3.ipynb`):

- Periodic folding: the periodic trajectory equals the unbounded one folded into the box.
- Shape function: exact uniform density, charge conservation across the seam, exact interpolation of linear fields.
- Field solver: second-order derivatives, Ampère/Faraday signs, Ohm's ideal, pressure and hyper-resistive terms, and a uniform equilibrium that stays exactly at rest.
- Cold-plasma circularly polarized eigenmodes ($k d_i = 1$): measured $\omega_R = -1.6159$ (analytic $-1.6180$) and $\omega_L = 0.6184$ (analytic $0.6180$). Total energy is conserved to $\sim 10^{-10}$.
- Hyper-resistive eigenmodes ($\nu = 0.002$): measured growth rates $\gamma_R = -0.00234$ (analytic $-0.00234$) and $\gamma_L = 0.00036$ (analytic $0.00034$).
- Hyper-resistivity only acts at grid scale: with $\nu = 10^{-4}$ and $10^{-3}$ the R and L eigenmodes ($k d_i = 1$) keep their analytic frequency (within 0.5%) and amplitude (within 2%), while a grid-scale whistler (4 cells per wavelength, $k = 16$) drops to 1.3% of its amplitude by $t = 0.1$ with $\nu = 10^{-3}$. Without $\nu$ it keeps 83%; that 17% loss is the ICN damping, $|g|^2 = 1 - y^4/4 + y^6/16$ per step at $y = \omega\,dt \approx 0.42$, and it shrinks roughly as $dt^3$ when $dt$ is reduced (down to a ~0.5% residual).
- Timestep limit: below the exact ICN boundary and within 20% of it for $\nu = 0.002$, $0.005$ and $0.01$.
- Thermal noise: the $(k, \omega)$ spectrum of $B_y + iB_z$ follows the analytic R and L branches.

## Project layout

```
PIC/
├── CMakeLists.txt
├── src/
│   ├── vector3.hpp             # 3D vector type + operators
│   ├── boris_pusher.hpp        # particle pusher (Species, Particle, EMFields, ...)
│   ├── grid.hpp                # periodic wrap + 1D Yee mesh: centering, Field, VecField
│   ├── shape.hpp               # linear shape: moment deposit, field gather
│   ├── field_solver.hpp        # Ampere, Faraday, Ohm's law
│   ├── hybrid_step.hpp         # ICN predictor-corrector time step
│   ├── loader.hpp              # uniform population + eigenmode initial wave
│   ├── parse.hpp               # value parsing helpers
│   ├── input.hpp               # .txt parsing (test-particle / hybrid mode)
│   ├── hybrid_input.hpp        # hybrid-mode keys and their validation
│   ├── cli.hpp                 # argv parsing
│   ├── hdf5_writer.hpp         # HDF5 output
│   └── main.cpp                # generic driver
├── examples/
│   ├── cyclotron.txt           # pure-B input
│   ├── drift.txt               # E×B input
│   ├── multi_cyclotron.txt     # pure-B input, 4 fixed-seed mixed species
│   ├── multi_drift.txt         # E×B input, same 4 particles
│   ├── periodic_drift.txt      # E×B drift along x in a periodic box
│   ├── alfven_R.txt            # hybrid: R-mode (whistler) eigenmode
│   ├── alfven_L.txt            # hybrid: L-mode (ion cyclotron) eigenmode
│   └── thermal_noise.txt       # hybrid: warm plasma, dispersion from noise
├── tests/
│   ├── scenarios.hpp           # test fixtures + analytic references
│   ├── test_boris_pusher.cpp   # Catch2 tests of the pusher
│   ├── test_input.cpp          # Catch2 tests of the .txt input parser
│   ├── test_periodic.cpp       # periodic wrap and folded trajectories
│   ├── test_shape.cpp          # deposit / gather
│   ├── test_field_solver.cpp   # Yee derivatives, Ampere, Faraday, Ohm
│   └── test_hybrid_step.cpp    # eigenmodes, hyper-resistivity, ICN order
├── analysis_session_2.ipynb    # Session 2: trajectories vs analytic references
├── analysis_session_3.ipynb    # Session 3: periodic drift, eigenmodes, dispersion
└── README.md
```

Input format (`key=value`, `#` = comment), two modes:

- **Test-particle mode**: `dt`, `steps`, `E`, `B`, then per-particle `q=`, `m=`, `x=`, `v=` lines. Uniform static `E`/`B` shared by all particles. Optional `length=` makes $x$ periodic in $[0, L)$.
- **Hybrid mode**, selected by an `nx=` line. Required: `dt`, `steps`, `B` (initial uniform $\vec B_0$), `length`, `nx`, `ppc` (protons per cell). Optional: `iterations` (ICN corrections, default 2), `output_every` (snapshot period in steps, default 1), `Te`, `eta`, `nu` (electron temperature, resistivity, hyper-resistivity; default 0), `density` (default 1), `vth` (thermal speed, default 0), `seed` (default 1), `wave_mode`, `wave_amplitude`, `wave_branch` (`R`/`L`; needs `B=1 0 0`, `density=1`). `E=` and explicit particles are rejected: $\vec E$ comes from Ohm's law.

## Build and run (C++)

```bash
cd PIC
mkdir -p build && cd build
cmake ..
cmake --build .
```

The first configure fetches [Catch2](https://github.com/catchorg/Catch2) v3.5.4 via CMake `FetchContent` (requires internet access).

```bash
./pic_prequeijo ../examples/cyclotron.txt results.h5 /cyclotron
./pic_prequeijo ../examples/drift.txt results.h5 /drift
./pic_prequeijo ../examples/multi_cyclotron.txt results.h5 /multi_cyclotron
./pic_prequeijo ../examples/multi_drift.txt results.h5 /multi_drift
./pic_prequeijo ../examples/periodic_drift.txt results.h5 /periodic_drift
./pic_prequeijo ../examples/alfven_R.txt results.h5 /alfven_R
./pic_prequeijo ../examples/alfven_L.txt results.h5 /alfven_L
./pic_prequeijo ../examples/thermal_noise.txt results.h5 /thermal_noise
```

Test-particle runs print a per-particle summary and write `t`, `x`, `y`, `z`, `vx`, `vy`, `vz` per group (in `particle_N` subgroups for the multi runs). Hybrid runs print the initial/final energies and $\max|\nabla\cdot\vec B|$, and write `t`, `x_primal`, `x_dual` and one `(snapshots × nx)` dataset per grid quantity: `Ex`, `Ey`, `Ez`, `Bx`, `By`, `Bz`, `n`, `ux`, `uy`, `uz`. Each of these has a `centering` attribute (`primal`/`dual`). `thermal_noise` takes a few seconds and writes ~20 MB.

Full usage: `./pic_prequeijo <input.txt> [output.h5] [group]`, or `--group <group>` (defaults `results.h5`, `/trajectory`).

Run the tests:

```bash
ctest --output-on-failure
# or directly:
./pic_prequeijo_test
```

## Analysis notebooks

`analysis_session_2.ipynb` reads `build/results.h5` with `h5py` and plots each Session 2 scenario against its analytic reference (including per-particle checks for the mixed-species runs). `analysis_session_3.ipynb` plots the periodic drift, the R/L eigenmodes (space-time, measured frequency, energy) and the $(k, \omega)$ spectrum of the thermal-noise run against the analytic dispersion branches.

```bash
cd PIC
jupyter notebook analysis_session_2.ipynb   # or analysis_session_3.ipynb
```

Run the C++ program first so `build/results.h5` exists before opening the notebook.

## Expected result

Session 2:

- `|v(t)|` flat at 1.0 to floating-point precision (~1e-15).
- Numeric trajectory overlays the analytic circle (max error ~6e-5).
- Mean `vy` converges to the drift velocity −0.1 (≈ −0.09996).
- Mixed species: per-particle `|v|` flat (~1e-15), per-particle trajectories within ~1e-4 of their own analytic curves, common drift for all.

Session 3:

- Periodic drift: $x \in [0, 4)$ always, 17 boundary crossings, folded orbit within $r_L = 0.9$ of the guiding center.
- R/L eigenmodes: $\omega$ within ~0.1% of analytic, energy conserved to ~$10^{-10}$, $\nabla\cdot\vec B = 0$ exactly.
- Thermal noise: spectral power on the R and L branches.

## Next steps

- 2D/3D Yee mesh (add the y/z directions to the centering-aware operators).
