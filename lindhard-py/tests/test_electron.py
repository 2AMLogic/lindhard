"""Electron runs through ElectronRun: CLI parity, thread invariance,
file-relative data paths, optional tallies and error mapping."""

import csv
import json
import shutil
import subprocess

import lindhard as lh
import numpy as np
import pytest
from conftest import EXAMPLES

HISTORIES = 40
SEED = 11
CARTESIAN = """
[electron.tally.cartesian]
x = { lo_nm = 0.0, hi_nm = 1500.0, bins = 5 }
y = { lo_nm = -500.0, hi_nm = 500.0, bins = 4 }
z = { lo_nm = -500.0, hi_nm = 500.0, bins = 3 }
"""


def _small_input() -> str:
    """The CLI example on a coarser table grid, with a Cartesian tally added
    (the example has a cylindrical one). Debug builds build tables slowly."""
    text = (EXAMPLES / "electron" / "e_10keV_si.toml").read_text()
    assert "points_per_decade = 8.0" in text
    text = text.replace("points_per_decade = 8.0", "points_per_decade = 3.0")
    return text.replace("[electron.tally.cylindrical]", CARTESIAN + "\n[electron.tally.cylindrical]")


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    """A directory holding the input and its optical-ELF data file."""
    d = tmp_path_factory.mktemp("electron")
    shutil.copy(EXAMPLES / "electron" / "synthetic_plasmon_elf.toml", d)
    (d / "e.toml").write_text(_small_input())
    return d / "e.toml"


@pytest.fixture(scope="module")
def cli_out(cli, case, tmp_path_factory):
    out = tmp_path_factory.mktemp("cli")
    subprocess.run(
        [str(cli), "run", str(case), "--out", str(out), "--histories", str(HISTORIES),
         "--seed", str(SEED), "--threads", "1"],
        check=True, capture_output=True,
    )
    return out


@pytest.fixture(scope="module")
def res1(case):
    return lh.ElectronRun.from_toml_file(case).run(histories=HISTORIES, seed=SEED, threads=1)


def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def floats(col):
    return np.array([np.nan if x == "" else float(x) for x in col])


def strip_nondeterministic(summary):
    """Drop the measured/thread-dependent parts: the `run` block (threads and
    timings) and the table origin (built or cached)."""
    summary = json.loads(json.dumps(summary))
    summary.pop("run")
    for m in summary["physics"]["materials"]:
        for t in ("elastic_table", "inelastic_table"):
            m[t].pop("source", None)
            m[t].pop("cache", None)
    return summary


def test_summary_matches_cli(cli_out, res1):
    cli = json.loads((cli_out / "electron_summary.json").read_text())
    assert strip_nondeterministic(res1.summary()) == strip_nondeterministic(cli)
    assert res1.summary()["run"]["threads"] == 1


def test_escape_spectra_equal_the_csv(cli_out, res1):
    _, rows = read_csv(cli_out / "electron_escape_spectra.csv")
    sp = res1.escape_spectra
    assert set(sp) == {"front", "back"}
    for face in ("front", "back"):
        for key, spectrum, cls in (
            ("energy", "energy_ev", "all"),
            ("polar_slow", "polar_deg", "slow"),
            ("polar_fast", "polar_deg", "fast"),
        ):
            sel = [r for r in rows if r[:3] == [face, spectrum, cls]]
            h = sp[face][key]
            n = len(h)
            body = list(zip(*sel[:n]))
            assert len(sel) == n + 2
            assert np.array_equal(h.edges[:-1], floats(body[3]))
            assert np.array_equal(h.edges[1:], floats(body[4]))
            assert np.array_equal(h.counts, np.array(body[5], dtype=np.uint64))
            assert np.array_equal(h.density, floats(body[6]))
            assert h.underflow == int(sel[n][5]) and h.overflow == int(sel[n + 1][5])
    assert sp["front"]["energy"].unit == "eV" and sp["front"]["polar_fast"].unit == "deg"


def test_cylindrical_equals_the_csv(cli_out, res1):
    _, rows = read_csv(cli_out / "electron_deposition_cylindrical.csv")
    cols = list(zip(*rows))
    c = res1.cylindrical
    assert c["energy_ev"].shape == (20, 40)
    assert c["r_edges_nm"].shape == (21,) and c["depth_edges_nm"].shape == (41,)
    assert np.array_equal(c["energy_ev"].ravel(), floats(cols[6]))
    assert np.array_equal(c["ev_per_primary_per_nm3"].ravel(), floats(cols[7]))
    assert np.array_equal(c["r_edges_nm"][:-1], floats(cols[2])[::40])
    assert np.array_equal(c["depth_edges_nm"][:-1], floats(cols[4])[:40])
    # The grid plus what fell outside is everything deposited in the layers.
    assert c["outside_ev"] == res1.summary()["results"]["deposition"]["cylindrical"]["outside_ev"]


def test_cartesian_equals_the_csv(cli_out, res1):
    _, rows = read_csv(cli_out / "electron_deposition_cartesian.csv")
    cols = list(zip(*rows))
    c = res1.cartesian
    assert c["energy_ev"].shape == (5, 4, 3)
    assert c["ev_per_primary_per_nm3"].shape == (5, 4, 3)
    # CSV rows run ix, iy, iz with iz fastest: C order.
    assert np.array_equal(c["energy_ev"].ravel(), floats(cols[9]))
    assert np.array_equal(c["ev_per_primary_per_nm3"].ravel(), floats(cols[10]))
    ix = np.array(cols[0], dtype=int).reshape(5, 4, 3)
    assert (ix[:, 0, 0] == np.arange(5)).all()
    assert np.array_equal(c["x_edges_nm"][:-1], floats(cols[3]).reshape(5, 4, 3)[:, 0, 0])
    assert np.array_equal(c["z_edges_nm"][1:], floats(cols[8]).reshape(5, 4, 3)[0, 0, :])


def test_thread_count_does_not_change_results(case, res1):
    four = lh.ElectronRun.from_toml_file(case).run(histories=HISTORIES, seed=SEED, threads=2)
    assert four.summary()["run"]["threads"] == 2
    assert strip_nondeterministic(four.summary()) == strip_nondeterministic(res1.summary())
    assert np.array_equal(four.cylindrical["energy_ev"], res1.cylindrical["energy_ev"])
    assert np.array_equal(four.cartesian["energy_ev"], res1.cartesian["energy_ev"])
    for face in ("front", "back"):
        for key, h in res1.escape_spectra[face].items():
            g = four.escape_spectra[face][key]
            assert np.array_equal(g.counts, h.counts)
            assert np.array_equal(g.density, h.density)


def test_disabled_tallies_are_none(case):
    head = _small_input().split("[electron.tally.cartesian]")[0]
    tail = '[target]\nsubstrate = "Si"\n[run]\nhistories = 5\nseed = 1\n'
    res = lh.ElectronRun.from_toml(head + "\n" + tail, base_dir=case.parent).run(threads=1)
    assert res.cartesian is None and res.cylindrical is None
    files = res.summary()["files"]
    assert files["deposition_cartesian"] is None and files["deposition_cylindrical"] is None
    # The escape spectra are always present.
    assert len(res.escape_spectra["front"]["energy"]) == 100


def test_roundtrip_and_overrides(case):
    run = lh.ElectronRun.from_toml_file(case)
    assert (run.histories, run.seed, run.threads) == (1000, 1, None)
    again = lh.ElectronRun.from_toml(run.to_toml(), base_dir=case.parent)
    assert again.to_toml() == run.to_toml()
    run.histories, run.seed, run.threads = 7, 3, 2
    assert (run.histories, run.seed, run.threads) == (7, 3, 2)
    assert "ElectronRun" in repr(run)


def test_file_relative_data_works_from_another_directory(case, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run = lh.ElectronRun.from_toml_file(case)
    assert run.base_dir == case.parent
    assert isinstance(run.validate(), list)
    # Without the base directory the relative data file is not found.
    with pytest.raises(lh.InputError, match="synthetic_plasmon_elf"):
        lh.ElectronRun.from_toml(case.read_text()).validate()
    # An explicit base directory works for text input.
    assert isinstance(lh.ElectronRun.from_toml(case.read_text(), base_dir=case.parent).validate(), list)


def test_errors(tmp_path, case):
    with pytest.raises(FileNotFoundError, match="missing.toml"):
        lh.ElectronRun.from_toml_file(tmp_path / "missing.toml")
    with pytest.raises(lh.InputError, match="bogus"):
        lh.ElectronRun.from_toml(case.read_text().replace("[run]", "[run]\nbogus = 1"))
    with pytest.raises(lh.InputError):
        lh.ElectronRun.from_toml("not = = toml")
    bad = lh.ElectronRun.from_toml_file(case)
    with pytest.raises(lh.InputError, match="histories"):
        bad.run(histories=0)
    with pytest.raises(lh.InputError, match="threads"):
        bad.run(histories=1, threads=0)
    # An ion input is not an electron input.
    with pytest.raises(lh.InputError):
        lh.ElectronRun.from_toml_file(next(EXAMPLES.glob("*.toml")))


def test_releases_the_gil(case):
    import threading
    import time

    run = lh.ElectronRun.from_toml_file(case)
    ticks = []
    stop = threading.Event()

    def tick():
        while not stop.is_set():
            ticks.append(time.perf_counter())
            time.sleep(0.01)

    t = threading.Thread(target=tick)
    t.start()
    try:
        run.run(histories=2, threads=1)
    finally:
        stop.set()
        t.join()
    assert len(ticks) > 5


def test_missing_data_file_names_it(case, tmp_path):
    (tmp_path / "e.toml").write_text(case.read_text())
    with pytest.raises(lh.InputError, match="synthetic_plasmon_elf"):
        lh.ElectronRun.from_toml_file(tmp_path / "e.toml").run(histories=1)
