"""Tests of the pic_prequeijo Python module (pybind11 bindings).

The physics is validated by the Catch2 suite (TESTS.md). These tests only
check the bridge: that Python runs exactly what the C++ driver runs (same
input, same seed, hence bit-identical arrays), that NumPy profiles are
loaded like a profile file, and that bad arguments raise clean errors.
"""

import subprocess

import h5py
import numpy as np
import pytest

import pic_prequeijo as pic
from conftest import PIC_DIR

EXAMPLES_DIR = PIC_DIR / "examples"
# Steps of the shortened example runs: enough for several snapshots.
SHORT_STEPS = 200
# Relative round-off allowed when C++ and NumPy evaluate the same formula.
ROUND_OFF = 4 * np.finfo(float).eps


# --- Helpers --------------------------------------------------------------------------


def shortened_input(tmp_path, example, steps=SHORT_STEPS):
    """Copies an example input into tmp_path with its steps= replaced."""
    lines = (EXAMPLES_DIR / example).read_text().splitlines()
    lines = [f"steps={steps}" if line.startswith("steps=") else line for line in lines]
    path = tmp_path / example
    path.write_text("\n".join(lines) + "\n")
    return path


def run_driver(driver, input_path, output_path, group):
    """Runs the C++ driver on one input into one HDF5 group."""
    subprocess.run(
        [driver, input_path, output_path, group], check=True, capture_output=True
    )


def read_group(path, group):
    """Reads the datasets (not the subgroups) and attributes of one HDF5 group."""
    with h5py.File(path, "r") as file:
        node = file[group]
        datasets = {
            key: item[()]
            for key, item in node.items()
            if isinstance(item, h5py.Dataset)
        }
        attributes = dict(node.attrs)
    return datasets, attributes


def assert_same_arrays(expected, actual):
    """Asserts two dicts hold the same keys and bit-identical values."""
    assert expected.keys() == actual.keys()
    for key in expected:
        assert np.array_equal(expected[key], actual[key]), key


# --- Hybrid runs ----------------------------------------------------------------------


def test_hybrid_history_matches_driver(tmp_path, driver):
    input_path = shortened_input(tmp_path, "alfven_R.txt")
    run_driver(driver, input_path, tmp_path / "driver.h5", "/run")
    driver_datasets, driver_attributes = read_group(tmp_path / "driver.h5", "/run")

    sim = pic.HybridSimulation.from_file(input_path)
    sim.run()
    assert_same_arrays(driver_datasets, sim.history())

    sim.save(tmp_path / "python.h5", "/run")
    python_datasets, python_attributes = read_group(tmp_path / "python.h5", "/run")
    assert_same_arrays(driver_datasets, python_datasets)
    assert_same_arrays(driver_attributes, python_attributes)


def test_run_defaults_to_the_input_steps(tmp_path):
    sim = pic.HybridSimulation.from_file(shortened_input(tmp_path, "alfven_L.txt"))
    sim.run(SHORT_STEPS // 2)
    sim.run()
    assert sim.step == SHORT_STEPS
    sim.run()
    assert sim.step == SHORT_STEPS
    assert sim.time == pytest.approx(SHORT_STEPS * 0.002, rel=ROUND_OFF)


def test_current_state_matches_last_snapshot(tmp_path):
    sim = pic.HybridSimulation.from_file(shortened_input(tmp_path, "alfven_R.txt"))
    sim.run()
    history, fields = sim.history(), sim.fields()
    assert fields["t"] == history["t"][-1]
    for name in ("Ex", "Ey", "Ez", "Bx", "By", "Bz", "n", "ux", "uy", "uz", "Ti"):
        assert np.array_equal(fields[name], history[name][-1]), name
    particles = sim.particles()
    assert particles["v"].shape == (len(particles["x"]), 3)
    energies = sim.energies()
    assert energies["total"] == energies["magnetic"] + energies["kinetic"]


# Small profile setup written both as a .txt + profile file and as NumPy
# arrays: the two loads must give bit-identical runs.
PROFILE_NX = 16
PROFILE_LENGTH = 8.0
PROFILE_SETTINGS = dict(
    nx=PROFILE_NX,
    length=PROFILE_LENGTH,
    dt=0.01,
    steps=20,
    output_every=5,
    B=(1.0, 0.0, 0.0),
    Te=0.5,
    ppc=8,
    seed=7,
)


def profile_columns():
    """Initial n, ux, uy, uz, vth on the primal nodes of the profile setup."""
    phase = 2 * np.pi * np.arange(PROFILE_NX) / PROFILE_NX
    return {
        "n": 1.0 + 0.1 * np.cos(phase),
        "ux": 0.01 * np.sin(phase),
        "uy": 0.02 * np.cos(phase),
        "uz": np.zeros(PROFILE_NX),
        "vth": 0.1 + 0.05 * np.sin(phase),
    }


def write_profile_input(tmp_path):
    """Writes the profile setup as an input file plus its profile file."""
    columns = profile_columns()
    table = np.column_stack(list(columns.values()))
    np.savetxt(
        tmp_path / "profile.txt",
        table,
        header=" ".join(columns),
        comments="",
        fmt="%.17g",
    )
    lines = [f"{key}={value}" for key, value in PROFILE_SETTINGS.items() if key != "B"]
    lines += ["B=1 0 0", "profile_file=profile.txt"]
    path = tmp_path / "input.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_numpy_profiles_match_profile_file(tmp_path):
    from_file = pic.HybridSimulation.from_file(write_profile_input(tmp_path))
    columns = profile_columns()
    from_arrays = pic.HybridSimulation(
        **PROFILE_SETTINGS,
        density=columns["n"],
        vth=columns["vth"],
        bulk_velocity=np.column_stack([columns["ux"], columns["uy"], columns["uz"]]),
    )
    from_file.run()
    from_arrays.run()
    assert_same_arrays(from_file.history(), from_arrays.history())


VALID_HYBRID = dict(nx=8, length=1.0, dt=1e-3, B=(1, 0, 0), ppc=2)


@pytest.mark.parametrize(
    ("changes", "error", "message"),
    [
        ({"density": np.ones(3)}, ValueError, "profile 'n' needs nx=8"),
        ({"vth": -np.ones(8)}, ValueError, "vth must be >= 0"),
        ({"bulk_velocity": np.zeros(8)}, ValueError, r"\(N, 3\) array"),
        ({"ppc": 0}, ValueError, "ppc must be > 0"),
        ({"wave_branch": "X"}, ValueError, "expected R or L"),
        ({"wave_mode": 1, "density": np.ones(8)}, ValueError, "uniform moments"),
        ({"foo": 1}, TypeError, "unexpected keyword argument 'foo'"),
    ],
)
def test_bad_hybrid_arguments_raise(changes, error, message):
    with pytest.raises(error, match=message):
        pic.HybridSimulation(**(VALID_HYBRID | changes))


def test_missing_hybrid_key_raises():
    settings = dict(VALID_HYBRID)
    del settings["ppc"]
    with pytest.raises(ValueError, match="missing required key 'ppc'"):
        pic.HybridSimulation(**settings)


def test_unstable_settings_warn():
    with pytest.warns(RuntimeWarning, match="exceeds the stability limit"):
        sim = pic.HybridSimulation(**(VALID_HYBRID | {"dt": 1.0}))
    assert sim.timestep_limit < 1.0


def test_wrong_mode_file_raises():
    with pytest.raises(ValueError, match="use run_test_particles"):
        pic.HybridSimulation.from_file(EXAMPLES_DIR / "cyclotron.txt")
    with pytest.raises(ValueError, match="use HybridSimulation"):
        pic.run_test_particles(EXAMPLES_DIR / "alfven_R.txt")


# --- Test particles -------------------------------------------------------------------


def test_test_particles_match_driver(tmp_path, driver):
    input_path = shortened_input(tmp_path, "multi_drift.txt")
    run_driver(driver, input_path, tmp_path / "driver.h5", "/run")
    trajectories = pic.run_test_particles(input_path)
    for index, trajectory in enumerate(trajectories):
        datasets, attributes = read_group(
            tmp_path / "driver.h5", f"/run/particle_{index}"
        )
        arrays = {key: trajectory[key] for key in datasets}
        assert_same_arrays(datasets, arrays)
        assert (trajectory["q"], trajectory["m"]) == (attributes["q"], attributes["m"])


def test_test_particle_arrays_match_file(tmp_path, driver):
    input_path = shortened_input(tmp_path, "multi_cyclotron.txt")
    run_driver(driver, input_path, tmp_path / "driver.h5", "/run")
    _, settings = read_group(tmp_path / "driver.h5", "/run")
    from_file = pic.run_test_particles(input_path)
    from_arrays = pic.run_test_particles(
        q=[p["q"] for p in from_file],
        m=[p["m"] for p in from_file],
        x=[[p["x"][0], p["y"][0], p["z"][0]] for p in from_file],
        v=[[p["vx"][0], p["vy"][0], p["vz"][0]] for p in from_file],
        E=settings["E"],
        B=settings["B"],
        dt=settings["dt"],
        steps=settings["steps"],
    )
    for expected, actual in zip(from_file, from_arrays):
        assert_same_arrays(expected, actual)


def test_bad_test_particle_arguments_raise():
    common = dict(E=(0, 0, 0), B=(0, 0, 1), dt=0.1, steps=1)
    with pytest.raises(ValueError, match="m must be > 0"):
        pic.run_test_particles(q=[1], m=[0], x=[[0, 0, 0]], v=[[1, 0, 0]], **common)
    with pytest.raises(ValueError, match="same number of particles"):
        pic.run_test_particles(q=[1, 1], m=[1], x=[[0, 0, 0]], v=[[1, 0, 0]], **common)


# --- Analytic references --------------------------------------------------------------


def test_wave_frequency_at_k_di_1():
    golden_ratio = (1 + np.sqrt(5)) / 2
    assert pic.wave_frequency(1.0, "L") == pytest.approx(
        golden_ratio - 1, rel=ROUND_OFF
    )
    assert pic.wave_frequency(1.0, "R") == pytest.approx(-golden_ratio, rel=ROUND_OFF)
    k = np.linspace(0.1, 3.0, 12).reshape(3, 4)
    assert pic.wave_frequency(k, "R").shape == k.shape
