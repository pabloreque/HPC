# Tests

The Catch2 suite (`tests/`, 44 test cases) checks the code against physics with a known answer: closed-form orbits, conservation laws, analytic dispersion relations, convergence orders and stability boundaries. Below is the physical motivation of each test and where each of its tolerances comes from.

Run them from `PIC/build`:

```bash
ctest --output-on-failure        # whole suite
./pic_prequeijo_test "[gyration]" # one tag
./pic_prequeijo_test --list-tests # names and tags
```

GitHub Actions runs the suite on every push (`.github/workflows/pic.yml`).

## How tolerances are chosen

No tolerance is fitted to a measured value. Each one is derived from the error it allows, using the helpers in `tests/tolerances.hpp`. There are four kinds:

| Kind | When | Tolerance |
|---|---|---|
| **Round-off** | The checked quantity is exact in exact arithmetic (an invariant or a closed-form identity), so only rounding can move it. | $N_\text{ops}\,\varepsilon\,\lvert x\rvert$: $N_\text{ops}$ floating-point operations of relative error $\varepsilon = 2.2\times10^{-16}$ each, added linearly (worst case). |
| **Truncation** | A discretization error with a known leading term, e.g. $(\omega\,\Delta t)^2/12$ for the Boris rotation or $(k\,\Delta x)^2/6$ for second-order spatial operators. | `kSafetyFactor` × leading term. The next-order terms are smaller by another power of the small parameter, so a factor of 2 covers them. |
| **Noise** | Moments of randomly drawn velocities. | Smoothing error + `kNoiseSigmas` = 5 standard deviations. A Gaussian exceeds $5\sigma$ with probability $6\times10^{-7}$, so a correct code does not fail by bad luck. |
| **Convergence** | The error ratio $e(h)/e(h/2)$ of an order-$p$ method. | Expected ratio $2^p$, corrected for the error of a numerical reference when there is one; margin $2 \cdot 2^p \times$ (relative size of the next-order term). |

Each leading-order estimate was checked against measurements: when the step changes, the measured error must follow the predicted law. For example, halving $\Delta t$ divides the E×B drift error by exactly 4, and the loader's smoothing error divided by its estimate tends to 1.000 as $\Delta x \to 0$. The measured values in the tables below come from the current code.

### Mutation testing

A test is only useful if it fails when the code is wrong. Each Boris test was run against deliberately broken versions of the pusher. Every mutation makes at least one test fail:

| Mutation | Failing test cases (of 8 `[boris]`) |
|---|---|
| Rotation sense reversed ($\vec t \to -\vec t$) | 7 |
| Mass dropped from the rotation ($q\,dt/2$ instead of $q\,dt/2m$) | 4 |
| One sign flipped in the $x$ component of `cross()` | 2 |
| Electric kick drops $E_z$ | 1 |
| Electric kick drops the mass | 1 |

The last three mutations are caught only by the multi-orientation E×B test. The single-geometry tests ($\vec B \parallel \hat z$, $q = m = 1$) cannot see them.

## Boris pusher (`test_boris_pusher.cpp`)

The test particles start from known states in uniform fields, where the exact orbit is known. Most tests use several species, a proton ($q = m = 1$), an alpha particle ($q = 2$, $m = 4$) and a light negative charge ($q = -1$, $m = 0.05$), in fields of magnitude 2. With $q = m = |\vec B| = 1$, a wrong formula such as $m/qB$ instead of $qB/m$ would pass by coincidence.

A property of the Boris step used below: the update $\vec v^{n+1} - \vec v^n = \frac{q\,\Delta t}{m}\,\frac{\vec v^n + \vec v^{n+1}}{2}\times\vec B$ is an exact rotation. It follows that $|\vec v|$, $v_\parallel$, $|\vec v_\perp|$ and the guiding center $\vec r + \frac{m}{qB^2}\vec v\times\vec B$ (apart from its motion along $\vec B$) are conserved exactly, not just to $O(\Delta t^2)$. What the scheme does get wrong is the phase: each step rotates $\vec v$ by $2\arctan(\omega\,\Delta t/2)$ instead of $\omega\,\Delta t$, so the gyration is slow by a relative $(\omega\,\Delta t)^2/12$.

| Test | Physics | Tolerance | Measured |
|---|---|---|---|
| Speed conservation, 8 steps per orbit | The magnetic force does no work. Boris rotates $\vec v$ exactly, so $\lvert\vec v\rvert$ is conserved even on a badly resolved orbit. | Round-off: 100 operations per step × 40 steps → $9\times10^{-13}$ | $6\times10^{-16}$ |
| Cyclotron circle, 1000 steps per orbit, 3 orbits | The orbit is the analytic circle of radius $r_L$. | Truncation: the orbit stays on the circle but lags by $N\,(\omega\Delta t)^3/12$ rad after $N$ steps, so the position error is $r_L$ times that → $1.2\times10^{-4}$. $z$ stays exactly 0 (no force along $\vec B$). | $6.20\times10^{-5}$ (prediction $6.20\times10^{-5}$) |
| E×B drift, scenario `drift` | Over whole gyro-periods the gyration averages out and the guiding center moves at $\vec E\times\vec B/B^2$. | Truncation: the run ends $(\omega\Delta t)^2/12$ of a turn short, so a fraction $(\omega\Delta t)^2/12$ of the gyration speed $\lvert\vec v_\perp - \vec v_E\rvert$ remains in the mean → $1.7\times10^{-4}$ | $8.27\times10^{-5}$ (prediction $8.27\times10^{-5}$) |
| Gyrofrequency and rotation sense, 3 species | $d\vec v/dt = (q/m)\,\vec v\times\vec B$ rotates $\vec v_\perp$ at $-q\vec B/m$: ions clockwise seen from the tip of $\vec B$, negative charges anticlockwise. The measured signed frequency checks the magnitude, the $q/m$ dependence and the sense together. | Truncation: relative $2\times(\omega\Delta t)^2/12 = 6.6\times10^{-6}$ | $3.29\times10^{-6}$ |
| Larmor radius about the guiding center | The distance to the analytic guiding center, which moves along $\vec B$ at constant $v_\parallel$, is $r_L = m v_\perp/(\lvert q\rvert B)$. $v_\parallel$ is constant. | Round-off (both are Boris invariants): 100 operations × 3000 steps | $1.3\times10^{-14}$ (relative) |
| Diamagnetism | The orbit's angular momentum about the guiding center is $-\mathrm{sign}(q)\,r_L v_\perp$ along $\vec B$: the current loop of every gyrating charge weakens $\vec B$. | Round-off (invariant) | $1.6\times10^{-14}$ (relative) |
| E×B drift in 7 geometries × 3 species | The drift $\vec E\times\vec B/B^2$ for the six pairs $E_iB_j$ of perpendicular axes (each component of the cross product with both signs) and one oblique pair mixing all components. The drift does not depend on $q$ or $m$, so the same expected value applies to all species. Expected drifts are written out by hand, not computed with `cross()`, so that a bug in `cross()` cannot cancel out. | Truncation, as for the `drift` scenario, with the gyration speed of each case; $v_\parallel$ to round-off | $\le 1.3\times10^{-6}$ (tolerance $\approx 2\times10^{-6}$) |

## Periodic boundaries (`test_periodic.cpp`)

| Test | Physics | Tolerance |
|---|---|---|
| `wrapIndex`, `wrapPosition` | Any index or position maps into the box, including negative values and values several boxes away, and $-10^{-18}$ must not round to $L$. | Round-off of one subtraction, $\varepsilon\lvert x\rvert$ |
| Folded trajectory | Periodic boundaries must not change the physics: the periodic orbit equals the unbounded one folded into $[0, L)$. The E×B drift scenario crosses the box at least 3 times. | Round-off: the two runs differ only by the wraps, $100\,N_\text{steps}\,\varepsilon\,\lvert x\rvert_{\max}$ |

## Shape function (`test_shape.cpp`)

The deposit and the gather use the linear (order-1 B-spline) shape. These tests check identities that the linear shape satisfies exactly, so all tolerances are round-off.

| Test | Physics |
|---|---|
| Uniform load | Evenly spaced particles of equal weight give every node exactly ppc particles' worth of shape weight, so $n = n_0$. |
| Charge across the seam | The shape weights of each particle add to 1, so the total deposited charge equals the total particle charge, including particles near $x = L$ that deposit on nodes $n_x - 1$ and 0. |
| Particle near $x = L$ | The weights $1 - s$ and $s$ go to the two nearest nodes across the periodic seam. |
| Linear fields | Linear interpolation reproduces constant and linear fields exactly, on both primal and dual nodes. |

## Initial load (`test_loader.cpp`)

The profiles are analytic: $n = 1 + 0.5\sin kx$, $\vec u = (0.1\cos kx,\ 0.3\sin kx,\ -0.2\cos kx)$, $v_{th} = 0.1\,(1 + 0.5\cos kx)$, with one wavelength in the box. They are sampled on the nodes, particles are loaded from them, and the deposited moments are compared with the analytic functions.

**The smoothing error.** Loading interpolates the node values linearly to the particles, and the deposit spreads them back with the same hat function. A profile $f$ therefore comes back as $f + h f''$ with $h = \Delta x^2/6$, and a product of interpolated profiles $fg$ loses $h f'g'$. From this:

- density: $h\,n''$;
- bulk velocity (ratio of the deposits of $n\vec u$ and $n$): $h\,(\vec u'' + n'\vec u'/n)$;
- temperature: $h\,(\lvert\vec u'\rvert^2/3 + v_{th}'^2 + 2v_{th}v_{th}'' + 2n'v_{th}v_{th}'/n)$. The first term is the spread of the bulk velocity across the kernel, which appears as heat.

The measured errors divided by these estimates are 0.987, 0.997, 0.999, 1.000 for $n_x$ = 16, 32, 64, 128 (density), and similarly for $\vec u$.

**The noise.** A node collects ppc particles from each of its two cells with hat weights $S$, so its effective number of independent particles is $(\sum S)^2/\sum S^2 = 1.5\,$ppc. A component of $\vec u$ then has standard deviation $v_{th}/\sqrt{N_\text{eff}}$, and $v_{th}$ (the square root of the sample variance of $3N_\text{eff}$ draws) has $v_{th}/\sqrt{6N_\text{eff}}$. The positions and weights do not depend on the random draws, so the density has no noise.

| Test | Tolerance | Measured |
|---|---|---|
| Cold load, 32 cells | 2 × smoothing: $n$ $6.4\times10^{-3}$, $\vec u$ $4.3\times10^{-3}$ | $3.2\times10^{-3}$, $2.2\times10^{-3}$ |
| Second-order convergence, 16 → 32 cells | Ratio 4, margin $2\cdot4\,(k\Delta x)^2/6 = 0.21$ | 3.94 ($n$), 3.93 ($\vec u$) |
| Weights follow the mesh, charge conserved | Each weight equals $n(x_p)/\text{ppc}$ within the linear-interpolation error $\Delta x^2\lvert n''\rvert/8$ (×2). The total charge equals the sum of the mesh density to round-off. | |
| Warm load, 16 cells, 4000 per cell | 2 × smoothing + $5\sigma$: $n$ $2.6\times10^{-2}$, $\vec u$ $3.3\times10^{-2}$, $v_{th}$ $1.9\times10^{-2}$ | $1.3\times10^{-2}$, $9.3\times10^{-3}$, $7.7\times10^{-3}$ |

## Input parser (`test_input.cpp`)

These tests do not check physics. They check that a malformed input fails loudly instead of running with wrong values. A valid file is read exactly; malformed numbers, unknown keys, missing required keys, inconsistent hybrid settings and bad profile files are rejected with a `path:line` message. A missing `q=` line, for example, must not silently reuse the previous particle's charge. The comparisons are exact, with no tolerances.

## Field solver (`test_field_solver.cpp`)

Each term of the field equations is checked on a sinusoid, whose exact derivative is known. The truncation tolerance is $2\,(k\Delta x)^2/6$ times the amplitude. The Taylor coefficients of the operators involved (centered derivative 1/24, two-point average 1/8, second derivative 1/12) are all at most 1/6.

| Test | Physics |
|---|---|
| Yee derivatives | The centered difference on the staggered mesh is second order, from primal to dual and back. Ratio 4, margin $2\cdot4\,(k\Delta x)^2/6 = 0.05$. |
| Ampère | $\vec j = \nabla\times\vec B$ with the right signs, $j_x = 0$ in 1D, and $\vec j$ lands on the $\vec E$ centering. |
| Faraday | $\partial_t\vec B = -\nabla\times\vec E$ with the right signs; $B_x$ is untouched, so $\nabla\cdot\vec B = 0$ is preserved exactly. |
| Ohm, ideal term | $\vec E = -\vec u\times\vec B$ for a uniform flow, exact to round-off. |
| Ohm, pressure term | $E_x = -T_e\,\partial_x n/n$ on the dual nodes. |
| Ohm, hyper-resistive term | $\vec E = -\nu\,\partial_x^2\vec j$. |
| Uniform equilibrium | A uniform plasma in a uniform $\vec B$ has $\vec E = 0$ and stays at rest for 50 steps, with $\vec B$ unchanged bit for bit. The tolerance on $\vec E$ is round-off. Its largest term is the pressure gradient that the deposit's round-off in $n$ produces over one cell, $T_e\,\delta n/\Delta x$. |

## Hybrid time step (`test_hybrid_step.cpp`)

These tests run the full PIC loop (push, deposit, Ohm, Faraday, ICN corrections) and compare it with the cold-plasma theory derived in the README. In normalized units $B_0 = n_0 = 1$, lengths are in $d_i$.

**Discretization error of a wave.** On the Yee mesh the curl, applied twice, changes $k^2$ by a relative $(k\Delta x)^2/12$. The linear shape smooths the ion response by $(k\Delta x)^2/12$ in each of the gather and the deposit. Since $\lvert d\ln\omega/d\ln k^2\rvert \le 1$, the frequency is off by at most $(k\Delta x)^2/4$. With $k\Delta x = 2\pi/64$ this is $2.4\times10^{-3}$, and the tolerance is twice that.

| Test | Physics | Tolerance | Measured |
|---|---|---|---|
| R mode (whistler), $kd_i = 1$ | The loaded eigenmode oscillates at $\omega_R = -(1+\sqrt5)/2$ and keeps its amplitude (ideal electrons, no dissipation). | Relative $4.8\times10^{-3}$ on $\omega$. Same margin on the amplitude: the loaded continuous eigenvector misses the discrete one by the dispersion error and excites the other branch by that much. | $1.3\times10^{-3}$; amplitude $1\times10^{-4}$ |
| L mode (ion cyclotron), $kd_i = 1$ | $\omega_L = (\sqrt5-1)/2$ | Same | $4.7\times10^{-4}$; amplitude $5.9\times10^{-4}$ |
| ICN convergence order | Two corrections give second order, one gives first order. | The reference ($\Delta t/16$) has its own error, which turns the expected ratio $2^p$ into $2^p(1 - 16^{-p})/(1 - 8^{-p})$: 4.047 and 2.143. Margin $2\cdot2^p\,\omega\Delta t$. | 4.048, 2.145 |
| Hyper-resistive growth rates, $\nu = 0.002$ | With $\eta' = \nu k^2$, the R mode is damped and the L mode grows: $\gamma = \mathrm{Im}\,\omega$ of $\omega^2 + a\omega - a = 0$, $a = k^2(1 + i\eta')$. | The ideal eigenvector misses the dissipative discrete one by a relative $(k\Delta x)^2/4 + \eta'$. This excites the other branch, and the beat (frequency $\omega_L - \omega_R$) biases a least-squares slope over a time $T$ by at most $12A/(WT^2)$. Plus the discretization error of $\gamma$ itself. About $1.2\times10^{-4}$. | R: $8\times10^{-7}$; L: $1.4\times10^{-5}$ |
| Timestep limit | `timestepLimit` uses radius 1.75 for the ICN stability region. The exact boundary along the ray of $\nu k^4 + i\omega$ lies between radii 1.75 and 2, so the estimate must stay below the exact limit and within a factor 1.75/2 of it. The exact limit is found by bisection on $\lvert g(z)\rvert = 1$ to double precision. For $\nu = 0$ it is exactly $\Delta x^2/2$ (round-off). | Exact bounds | |
| Long waves with small $\nu$ ($10^{-4}$, $10^{-3}$) | Hyper-resistivity acts as $\nu k^4$, so at $kd_i = 1$ the waves keep their frequency, and their amplitude changes only at the slow analytic rate $e^{\gamma t}$. | Dispersion error, as for the ideal modes | |
| Grid-scale whistlers (4 cells per wavelength) | Without $\nu$ the wave survives, with only the ICN numerical damping $\lvert g\rvert^N$ (~83 %). $\nu = 10^{-3}$, harmless at $kd_i = 1$, removes it (~1 %). Both are compared with the analytic ICN amplification of the pure whistler, with the Yee wavenumber $(2/\Delta x)\sin(k\Delta x/2)$. | The pure-whistler formula leaves out the ions, which carry a fraction $1/(kd_i)^2$ of the response, so $\ln$ amplitude can differ by $N/(kd_i)^2$. Margin $2N/(kd_i)^2 = 0.48$ in $\ln$. | 0.005 and 0.053 in $\ln$ |
| ICN stability at 2 cells per wavelength | The fastest grid whistler decays at $0.8\times$ the limit and grows at $1.2\times$ (2 corrections), and grows at any $\Delta t$ with 1 correction, each at the analytic amplification $\lvert g(i\omega\Delta t)\rvert^{10}$. | Ion response as above: $2\cdot10\,(\Delta x/2)^2 = 0.048$ in $\ln$. Verified: with 256 instead of 64 cells the deviations drop by 16, as $\Delta x^2$. | $\le 0.025$ in $\ln$ |
