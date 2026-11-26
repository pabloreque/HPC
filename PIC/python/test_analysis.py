"""Tests of pic_analysis, the shared helpers of the analysis notebooks."""

import h5py
import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

import pic_analysis as pa  # noqa: E402

# Relative round-off of a least-squares phase slope over a few hundred
# samples: well above eps, far below any physical effect.
FIT_ROUND_OFF = 1e-12


def synthetic_wave(k, omega, nx=32, snapshots=200, length=2 * np.pi):
    """Snapshots of b = exp(i(kx - omega t)) on the dual nodes."""
    dx = length / nx
    x_dual = (np.arange(nx) + 0.5) * dx
    t = np.linspace(0.0, 10.0, snapshots)
    b = np.exp(1j * (k * x_dual[None, :] - omega * t[:, None]))
    return {"t": t, "x_dual": x_dual, "By": b.real, "Bz": b.imag}


@pytest.mark.parametrize("omega", [0.618, -1.618])
def test_mode_frequency_recovers_a_pure_wave(omega):
    k = 1.0
    measured, amplitude = pa.mode_frequency(synthetic_wave(k, omega), k)
    assert measured == pytest.approx(omega, rel=FIT_ROUND_OFF)
    assert np.allclose(amplitude, 1.0)


def test_load_group_reads_datasets_and_settings(tmp_path):
    path = tmp_path / "run.h5"
    with h5py.File(path, "w") as file:
        group = file.create_group("run")
        group["t"] = np.arange(3.0)
        group.attrs["dt"] = 0.5
        group.create_group("particle_0")
    data, settings = pa.load_group(path, "run")
    assert list(data) == ["t"]
    assert settings == {"dt": 0.5}


def test_load_particles_orders_subgroups_numerically(tmp_path):
    path = tmp_path / "run.h5"
    with h5py.File(path, "w") as file:
        group = file.create_group("run")
        group.attrs["dt"] = 0.1
        for index in (0, 2, 10, 1):
            sub = group.create_group(f"particle_{index}")
            sub["t"] = np.full(2, float(index))
            sub.attrs["q"] = float(index)
    trajectories, settings = pa.load_particles(path, "run")
    assert [attributes["q"] for _, attributes in trajectories] == [0, 1, 2, 10]
    assert settings["dt"] == 0.1


def test_plot_wave_branches_draws_both_labelled_branches():
    _, ax = plt.subplots()
    pa.plot_wave_branches(ax, np.linspace(-3, 3, 40))
    labels = [line.get_label() for line in ax.get_lines()]
    assert labels == ["L branch (analytic)", "R branch (analytic)"]
    plt.close("all")
