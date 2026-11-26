"""Shared helpers of the PIC analysis notebooks.

Plot style and colors, reading a run back from the HDF5 file of the C++
driver (datasets and the input settings stored as attributes), measuring
the frequency of one Fourier mode, and drawing the analytic R and L
dispersion branches. Only code used by more than one notebook lives here.

From a notebook in ``PIC/``::

    import sys
    sys.path[:0] = ["python", "build"]
    import pic_analysis as pa
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np

import pic_prequeijo as pic

# --- Style and colors ------------------------------------------------------------------

# Categorical colors in a fixed, colorblind-validated order: series take
# them in this order, never cycled past the eighth.
SERIES_COLORS = (
    "#257de0",
    "#f47230",
    "#00b894",
    "#f2a600",
    "#ee7fac",
    "#008a17",
    "#6a4ce2",
    "#ec4647",
)
# Analytic references are drawn in neutral ink, under the numeric series.
ANALYTIC_COLOR = "#666666"
INK_COLOR = "#191919"
GRID_COLOR = "#e7e7e7"
AXIS_COLOR = "#c2c2c2"
L_BRANCH_COLOR = SERIES_COLORS[0]
R_BRANCH_COLOR = SERIES_COLORS[1]
# Signed fields (B_y, u_y): two hues around a neutral midpoint.
DIVERGING_CMAP = "RdBu_r"
# Magnitudes (spectral power): one hue, light to dark.
SEQUENTIAL_CMAP = "Blues"
LINE_WIDTH = 2.0


def apply_style():
    """Sets the matplotlib defaults shared by every notebook figure.

    Thin lines, recessive solid hairline grid and axes, ink-colored text,
    and the categorical colors in their fixed order.
    """
    plt.rcParams.update(
        {
            "axes.prop_cycle": plt.cycler(color=SERIES_COLORS),
            "axes.edgecolor": AXIS_COLOR,
            "axes.labelcolor": INK_COLOR,
            "axes.titlesize": "medium",
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": GRID_COLOR,
            "grid.linewidth": 0.8,
            "grid.linestyle": "-",
            "lines.linewidth": LINE_WIDTH,
            "lines.solid_capstyle": "round",
            "xtick.color": INK_COLOR,
            "ytick.color": INK_COLOR,
            "legend.frameon": False,
            "legend.fontsize": "small",
            "figure.dpi": 110,
        }
    )


# --- Loading results -------------------------------------------------------------------


def load_group(path, group):
    """Reads one run of the C++ driver back from its HDF5 file.

    Parameters
    ----------
    path : str or os.PathLike
        HDF5 file written by ``pic_prequeijo`` (or ``HybridSimulation.save``).
    group : str
        Group of the run, e.g. ``"alfven_R"``.

    Returns
    -------
    data : dict of ndarray
        Every dataset of the group (subgroups are left out).
    settings : dict
        The input settings stored as attributes, named like the .txt keys.
    """
    with h5py.File(path, "r") as file:
        node = file[group]
        data = {
            key: item[()]
            for key, item in node.items()
            if isinstance(item, h5py.Dataset)
        }
        settings = dict(node.attrs)
    return data, settings


def load_particles(path, group):
    """Reads every trajectory of a test-particle run.

    Parameters
    ----------
    path : str or os.PathLike
        HDF5 file written by ``pic_prequeijo``.
    group : str
        Group of the run: flat for one particle, ``particle_N`` subgroups
        otherwise.

    Returns
    -------
    trajectories : list of (dict, dict)
        For each particle, its datasets (t, x, y, z, vx, vy, vz) and its
        attributes (q, m).
    settings : dict
        The run settings (dt, steps, E, B, ...).
    """
    with h5py.File(path, "r") as file:
        names = sorted(
            (key for key in file[group] if key.startswith("particle_")),
            key=lambda name: int(name.split("_")[1]),
        )
    if not names:
        data, settings = load_group(path, group)
        return [(data, settings)], settings
    _, settings = load_group(path, group)
    return [load_group(path, f"{group}/{name}") for name in names], settings


# --- Diagnostics -----------------------------------------------------------------------


def mode_frequency(history, k):
    """Measures the frequency of the Fourier mode k of b = B_y + i B_z.

    The mode b_k(t) = mean over x of b(x, t) exp(-ikx) evolves as
    exp(-i omega t) for a wave exp(i(kx - omega t)); omega is minus the
    least-squares slope of its unwrapped phase.

    Parameters
    ----------
    history : dict of ndarray
        Snapshots of a hybrid run, as read by ``load_group`` or returned by
        ``HybridSimulation.history()``: t, x_dual, By, Bz.
    k : float
        Wavenumber of the mode, in 1/d_i.

    Returns
    -------
    omega : float
        Signed frequency in Omega_ci (> 0 for the L wave, < 0 for R).
    amplitude : ndarray
        |b_k(t)| at every snapshot.
    """
    b = history["By"] + 1j * history["Bz"]
    mode = (b * np.exp(-1j * k * history["x_dual"])).mean(axis=1)
    phase = np.unwrap(np.angle(mode))
    slope = np.polyfit(history["t"], phase, 1)[0]
    return -slope, np.abs(mode)


# --- Plots -----------------------------------------------------------------------------


def plot_wave_branches(ax, k):
    """Draws the analytic cold-plasma L and R branches omega(|k|).

    Both branches are drawn in neutral ink, the L branch solid and the R
    branch dashed, so they stay readable over a colored image.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    k : ndarray
        Wavenumbers, in 1/d_i, sorted and without k = 0; negative values
        are drawn with the frequency of |k|.
    """
    magnitude = np.abs(k)
    ax.plot(
        k,
        pic.wave_frequency(magnitude, "L"),
        color=INK_COLOR,
        linewidth=1.2,
        label="L branch (analytic)",
    )
    ax.plot(
        k,
        pic.wave_frequency(magnitude, "R"),
        color=INK_COLOR,
        linewidth=1.2,
        linestyle="--",
        label="R branch (analytic)",
    )
