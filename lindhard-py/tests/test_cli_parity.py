"""A run matches the `lindhard` command bit for bit for the same seed, at any
thread count."""

import csv
import json
import subprocess

import lindhard as lh
import numpy as np
import pytest
from conftest import EXAMPLES

IONS = 100
SEED = 7


def run_cli(cli, example, out, threads=1):
    subprocess.run(
        [str(cli), "run", str(example), "--out", str(out), "--ions", str(IONS),
         "--seed", str(SEED), "--threads", str(threads)],
        check=True, capture_output=True,
    )


def python_run(example, threads=1):
    return lh.Run.from_toml_file(example).run(ions=IONS, seed=SEED, threads=threads)


def strip_run(summary):
    summary = dict(summary)
    summary.pop("run")
    return summary


def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def floats(col):
    return np.array([np.nan if x == "" else float(x) for x in col])


def test_files_are_byte_identical(cli, example, tmp_path):
    cli_out, py_out = tmp_path / "cli", tmp_path / "py"
    run_cli(cli, example, cli_out)
    python_run(example).write(py_out)
    names = sorted(p.name for p in cli_out.iterdir())
    assert names == sorted(p.name for p in py_out.iterdir())
    for name in names:
        a, b = (cli_out / name).read_bytes(), (py_out / name).read_bytes()
        if name == "summary.json":
            a, b = (strip_run(json.loads(x)) for x in (a, b))
            # `software.git_describe` is compared too. Two builds of the same
            # commit must agree; a wheel built without git (a manylinux
            # container) reports "unknown", so CI passes one
            # LINDHARD_GIT_DESCRIBE to both builds (lindhard-cli/build.rs).
            assert a["software"] == b["software"], (
                "build provenance differs between the CLI and the wheel; build "
                "both with the same LINDHARD_GIT_DESCRIBE"
            )
        assert a == b, name


def test_arrays_equal_the_csv_values(cli, example, tmp_path):
    run_cli(cli, example, tmp_path)
    res = python_run(example)

    _, rows = read_csv(tmp_path / "depth_profile.csv")
    cols = list(zip(*rows))
    h = res.depth
    assert h.counts.dtype == np.uint64
    assert np.array_equal(h.counts, np.array(cols[2], dtype=np.uint64))
    assert np.array_equal(h.edges[:-1], floats(cols[0]))
    assert np.array_equal(h.edges[1:-1], floats(cols[1])[:-1]) and np.isinf(h.edges[-1])
    assert np.array_equal(h.density, floats(cols[3]), equal_nan=True)

    _, rows = read_csv(tmp_path / "lateral_profile.csv")
    for name, hist in (("y", res.lateral_y), ("z", res.lateral_z), ("radial", res.radial)):
        sel = [r for r in rows if r[0] == name]
        body = sel[:-2]
        assert np.array_equal(hist.counts, np.array([r[3] for r in body], dtype=np.uint64))
        assert np.array_equal(hist.edges[:-1], floats([r[1] for r in body]))
        assert np.array_equal(hist.edges[1:], floats([r[2] for r in body]))
        assert np.array_equal(hist.density, floats([r[4] for r in body]))
        assert (hist.underflow, hist.overflow) == (int(sel[-2][3]), int(sel[-1][3]))

    _, rows = read_csv(tmp_path / "damage_profile.csv")
    cols = list(zip(*rows[:-1]))
    for k, hist in ((2, res.vacancies), (3, res.interstitials), (4, res.replacements)):
        assert np.array_equal(hist.counts, np.array(cols[k], dtype=np.uint64))
        assert np.array_equal(hist.density, floats(cols[k + 3]))
    assert np.array_equal(res.vacancies.edges[:-1], floats(cols[0]))

    _, rows = read_csv(tmp_path / "escape_spectra.csv")
    assert res.escapes
    for esc in res.escapes:
        for spec, key in (("energy_ev", "energy"), ("polar_deg", "polar")):
            sel = [r for r in rows if (r[0], r[3], r[4]) == (str(esc["z"]), esc["face"], spec)]
            body = sel[:-2]
            h = esc[key]
            assert np.array_equal(h.counts, np.array([r[7] for r in body], dtype=np.uint64))
            assert np.array_equal(h.edges[:-1], floats([r[5] for r in body]))
            assert np.array_equal(h.density, floats([r[8] for r in body]))


def test_ions_arrays_equal_ions_csv(cli, tmp_path):
    example = EXAMPLES / "ar_1keV_cu.toml"
    run_cli(cli, example, tmp_path)
    ions = python_run(example).ions
    head, rows = read_csv(tmp_path / "ions.csv")
    cols = list(zip(*rows))
    assert np.array_equal(ions["index"], np.array(cols[0], dtype=np.uint64))
    assert [lh.FATE_NAMES[c] for c in ions["fate"]] == list(cols[1])
    assert np.array_equal(ions["position_nm"], np.column_stack([floats(cols[i]) for i in (2, 3, 4)]))
    assert np.array_equal(ions["energy_ev"], floats(cols[5]))
    assert np.array_equal(ions["direction"], np.column_stack([floats(cols[i]) for i in (6, 7, 8)]))
    assert np.array_equal(ions["layer"], np.array(cols[9], dtype=np.uint64))


def test_summary_dict(example):
    res = python_run(example)
    s = res.summary()
    assert s["format"]["name"] == "lindhard-summary"
    assert s["results"]["histories"] == IONS == res.histories
    assert "yields" in s["results"] and "run" in s
    assert json.loads(res.summary_json()) == s


def test_independent_of_thread_count(example):
    one, three = python_run(example, 1), python_run(example, 3)
    assert strip_run(one.summary()) == strip_run(three.summary())
    assert np.array_equal(one.depth.counts, three.depth.counts)
    assert np.array_equal(one.vacancies.counts, three.vacancies.counts)
    assert one.summary()["run"]["threads"] == 1
    assert three.summary()["run"]["threads"] == 3


def test_per_ion_false_has_no_ions():
    run = lh.Run.from_toml_file(EXAMPLES / "ar_1keV_cu.toml")
    t = run.tally
    t.per_ion = False
    run.tally = t
    assert run.run(ions=20).ions is None


def test_same_seed_reproduces_and_seed_matters():
    run = lh.Run.from_toml_file(EXAMPLES / "ar_1keV_cu.toml")
    a, b = run.run(ions=200, seed=1), run.run(ions=200, seed=1)
    c = run.run(ions=200, seed=2)
    assert np.array_equal(a.depth.counts, b.depth.counts)
    assert not np.array_equal(a.ions["energy_ev"], c.ions["energy_ev"])


# Directory reuse and failure behaviour follows the command's output lifecycle.

REUSE_EXAMPLE = EXAMPLES / "ar_1keV_cu.toml"


def small_run(per_ion):
    run = lh.Run.from_toml_file(REUSE_EXAMPLE)
    t = run.tally
    t.per_ion = per_ion
    run.tally = t
    return run.run(ions=20, seed=SEED)


def test_write_removes_stale_ions_csv_and_keeps_unrelated_files(tmp_path):
    small_run(True).write(tmp_path)
    assert (tmp_path / "ions.csv").is_file()
    (tmp_path / "notes.txt").write_text("mine")
    small_run(False).write(tmp_path)
    assert not (tmp_path / "ions.csv").exists()
    assert (tmp_path / "summary.json").is_file()
    assert (tmp_path / "notes.txt").read_text() == "mine"


def test_write_obstructed_csv_on_fresh_dir_names_path_and_leaves_no_summary(tmp_path):
    (tmp_path / "damage_profile.csv").mkdir()
    with pytest.raises(OSError, match="damage_profile.csv"):
        small_run(True).write(tmp_path)
    assert not (tmp_path / "summary.json").exists()


def test_write_failure_on_reused_dir_leaves_no_old_summary(tmp_path):
    res = small_run(True)
    res.write(tmp_path)
    assert (tmp_path / "summary.json").is_file()
    (tmp_path / "escape_spectra.csv").unlink()
    (tmp_path / "escape_spectra.csv").mkdir()
    with pytest.raises(OSError, match="escape_spectra.csv"):
        res.write(tmp_path)
    assert not (tmp_path / "summary.json").exists()


def test_write_failed_summary_invalidation_changes_no_csv(tmp_path):
    names = ["depth_profile.csv", "lateral_profile.csv", "damage_profile.csv",
             "escape_spectra.csv"]
    for n in names:
        (tmp_path / n).write_text("old")
    (tmp_path / "summary.json").mkdir()
    (tmp_path / "summary.json" / "inner").write_text("x")
    with pytest.raises(OSError, match="summary.json"):
        small_run(True).write(tmp_path)
    for n in names:
        assert (tmp_path / n).read_text() == "old", n
