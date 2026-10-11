"""`Run.run_dynamic()`: fluence-dependent runs match the `lindhard` command,
at any thread count, and `write()` follows the command's publication policy."""

import csv
import json
import subprocess
import warnings

import lindhard as lh
import numpy as np
import pytest
from conftest import EXAMPLES

FILM = EXAMPLES / "dynamic" / "as_1keV_si_film.toml"
IONS = 400
SEED = 3

# A thin film sputtered hard enough that front slabs empty and the slab count
# changes (erosion on).
EROSION = """
[beam]
ion = "As"
energy_ev = 1000.0
[target]
substrate = "Si"
[[target.layers]]
material = "Si"
thickness_nm = 4.0
[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0
[run]
ions = 600
seed = 2
[dynamic]
fluence_cm2 = 4.0e17
ions_per_step = 100
max_change = 0.25
slab_nm = 1.0
erosion = true
[dynamic.energies.As]
e_d_ev = 15.0
"""

STATIC = """
[beam]
ion = "B"
energy_ev = 5000.0
[target]
substrate = "Si"
[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0
[run]
ions = 20
seed = 1
"""


def run_cli(cli, text, out, tmp_path, threads=1):
    src = tmp_path / "in.toml"
    src.write_text(text)
    subprocess.run(
        [str(cli), "run", str(src), "--out", str(out), "--ions", str(IONS),
         "--seed", str(SEED), "--threads", str(threads)],
        check=True, capture_output=True,
    )


def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def strip_run(summary):
    summary = dict(summary)
    summary.pop("run")
    return summary


def fingerprint(res):
    """Everything thread-independent, as bytes (so NaN compares equal)."""
    steps = {k: v.tobytes() for k, v in res.steps.items()}
    comp = [
        {k: (v.tobytes() if hasattr(v, "tobytes") else v) for k, v in c.items()}
        for c in res.composition
    ]
    return steps, comp, strip_run(res.summary())


@pytest.fixture(params=["film", "erosion"])
def text(request):
    return FILM.read_text() if request.param == "film" else EROSION


def test_arrays_equal_the_cli_csv_values(cli, text, tmp_path):
    run_cli(cli, text, tmp_path / "cli", tmp_path)
    res = lh.Run.from_toml(text).run_dynamic(ions=IONS, seed=SEED, threads=1)

    header, rows = read_csv(tmp_path / "cli" / "dynamic_steps.csv")
    cols = dict(zip(header, zip(*rows)))
    s = res.steps
    assert res.n_steps == len(rows) and len(res.composition) == len(rows)
    assert s["step"][0] == 0 and s["ions_done"][0] == 0
    ints = {"step": "step", "first_index": "first_index", "ions": "ions",
            "ions_done": "ions_done", "attempts": "attempts",
            "clamped": "clamped", "removed_slabs": "removed_slabs",
            "n_slabs": "n_slabs", "cum_backscattered": "cum_backscattered",
            "cum_transmitted": "cum_transmitted",
            "cum_stopped_in_target": "cum_stopped_in_target",
            "cum_stopped_in_substrate": "cum_stopped_in_substrate",
            "cum_sputtered": "cum_sputtered",
            "cum_recoils_transmitted": "cum_recoils_transmitted"}
    for key, col in ints.items():
        assert s[key].dtype == np.uint64, key
        assert np.array_equal(s[key], np.array(cols[col], dtype=np.uint64)), key
    for key in ("fluence_cm2", "max_change", "surface_nm", "thickness_nm"):
        assert np.array_equal(s[key], np.array(cols[key], dtype=float)), key
    cli_yield = np.array([np.nan if x == "" else float(x) for x in cols["sputter_yield"]])
    assert np.array_equal(s["sputter_yield"], cli_yield, equal_nan=True)
    for j, sym in enumerate(res.species):
        assert np.array_equal(
            s["cum_sputtered_by_species"][:, j],
            np.array(cols[f"cum_sputtered_{sym}"], dtype=np.uint64),
        )

    header, rows = read_csv(tmp_path / "cli" / "dynamic_composition.csv")
    cols = dict(zip(header, zip(*rows)))
    step = np.array(cols["step"], dtype=np.uint64)
    for i, c in enumerate(res.composition):
        m = step == i
        assert c["step"] == i
        assert c["atoms_per_cm2"].shape == (s["n_slabs"][i], len(res.species))
        assert m.sum() == s["n_slabs"][i]
        for key in ("front_nm", "back_nm", "thickness_nm"):
            assert np.array_equal(c[key], np.array(cols[key], dtype=float)[m]), key
        for j, sym in enumerate(res.species):
            a = np.array(cols[f"atoms_per_cm2_{sym}"], dtype=float)[m]
            f = np.array(cols[f"fraction_{sym}"], dtype=float)[m]
            assert np.array_equal(c["atoms_per_cm2"][:, j], a)
            assert np.array_equal(c["fraction"][:, j], f, equal_nan=True)

    a = strip_run(json.loads((tmp_path / "cli" / "dynamic_summary.json").read_text()))
    assert a == strip_run(res.summary())


def test_finite_film_returns_step_zero_and_all_accepted_steps():
    res = lh.Run.from_toml_file(FILM).run_dynamic(ions=IONS, seed=SEED)
    assert res.species == ["Si", "As"]
    assert list(res.species_z) == sorted(res.species_z)
    steps = res.steps["step"]
    assert list(steps) == list(range(len(steps)))
    assert res.steps["ions_done"][-1] == IONS
    assert res.summary()["totals"]["steps"] == steps[-1]


def test_erosion_fixture_changes_slab_count():
    res = lh.Run.from_toml(EROSION).run_dynamic()
    n = res.steps["n_slabs"]
    assert n[0] == 4 and n.min() < n[0]
    assert res.steps["surface_nm"][-1] > 0
    # Ragged snapshots: each has its own length.
    assert {len(c["front_nm"]) for c in res.composition} == set(n.tolist())


@pytest.mark.parametrize("which", ["film", "erosion"])
def test_results_do_not_depend_on_threads(which):
    run = lh.Run.from_toml(FILM.read_text() if which == "film" else EROSION)
    ref = fingerprint(run.run_dynamic(ions=IONS, seed=SEED, threads=1))
    for n in (2, 4):
        assert fingerprint(run.run_dynamic(ions=IONS, seed=SEED, threads=n)) == ref


def test_write_reproduces_cli_files(cli, text, tmp_path):
    run_cli(cli, text, tmp_path / "cli", tmp_path)
    res = lh.Run.from_toml(text).run_dynamic(ions=IONS, seed=SEED, threads=2)
    res.write(tmp_path / "py")
    names = sorted(p.name for p in (tmp_path / "cli").iterdir())
    assert names == sorted(p.name for p in (tmp_path / "py").iterdir())
    for name in names:
        a, b = (tmp_path / "cli" / name).read_bytes(), (tmp_path / "py" / name).read_bytes()
        if name.endswith(".json"):
            a, b = (strip_run(json.loads(x)) for x in (a, b))
        assert a == b, name
    assert json.loads(res.summary_json()) == res.summary()


def test_write_invalidates_old_summary_when_publication_fails(tmp_path):
    res = lh.Run.from_toml(EROSION).run_dynamic(ions=100, seed=1)
    res.write(tmp_path)
    assert (tmp_path / "dynamic_summary.json").is_file()
    # Obstruct the composition CSV: the replacement fails after the old
    # summary was removed, so no summary describes the half-written set.
    (tmp_path / "dynamic_composition.csv").unlink()
    (tmp_path / "dynamic_composition.csv").mkdir()
    with pytest.raises(OSError, match="dynamic_composition.csv"):
        res.write(tmp_path)
    assert not (tmp_path / "dynamic_summary.json").exists()


def test_missing_dynamic_section_is_an_input_error():
    run = lh.Run.from_toml(STATIC)
    with pytest.raises(lh.InputError, match=r"\[dynamic\]"):
        run.run_dynamic()


def test_threads_zero_is_an_input_error():
    with pytest.raises(lh.InputError, match="threads"):
        lh.Run.from_toml(EROSION).run_dynamic(threads=0)


def test_static_run_still_works_and_rejects_dynamic_input():
    res = lh.Run.from_toml(STATIC).run()
    assert isinstance(res, lh.RunResult)
    with pytest.raises(lh.LindhardError, match="dynamic"):
        lh.Run.from_toml(EROSION).run()


def test_per_call_overrides_do_not_change_the_run():
    run = lh.Run.from_toml(EROSION)
    run.run_dynamic(ions=100, seed=9, threads=1)
    assert (run.ions, run.seed, run.threads) == (600, 2, None)
