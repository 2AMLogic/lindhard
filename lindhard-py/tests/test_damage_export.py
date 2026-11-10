"""The damage-export file: written by a run, read back as NumPy arrays that
equal the run's own arrays exactly."""

import numpy as np
import pytest

import lindhard as lh

BASE = """
[beam]
ion = "B"
energy_ev = 5000.0

[materials.SiO2]
density_g_cm3 = 2.2
elements = [
  { symbol = "Si", atom_fraction = 1.0 },
  { symbol = "O", atom_fraction = 2.0, e_d_ev = 20.0, e_s_ev = 2.0 },
]

[target]
substrate = "Si"
[[target.layers]]
material = "SiO2"
thickness_nm = 20.0

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 100
seed = 9

[tally]
depth_bin_nm = 2.0
depth_bins = 50
per_ion = false
"""

VOXEL = """
[tally.voxel]
origin_nm = [0.0, -10.0, -10.0]
spacing_nm = [4.0, 10.0, 10.0]
dims = [10, 2, 2]
"""

ARRAYS = [
    "implanted", "displacements", "electronic_ev", "lattice_ev", "rest_ev",
    "deposited_ev", "vacancies", "interstitials", "replacements",
]


def same(a, b):
    assert a.keys() == b.keys()
    for k in a:
        if isinstance(a[k], np.ndarray):
            assert a[k].dtype == b[k].dtype and a[k].shape == b[k].shape, k
            assert a[k].tobytes() == b[k].tobytes(), k  # bit for bit
        elif isinstance(a[k], dict):
            same(a[k], b[k])
        else:
            assert a[k] == b[k], k


def test_written_file_reads_back_unchanged(tmp_path):
    result = lh.Run.from_toml(BASE + VOXEL).run(ions=100, seed=9, threads=1)
    live = result.damage_export()
    assert live is not None
    result.write(tmp_path)
    back = lh.read_damage_export(tmp_path / "damage_export.csv")
    same(live, back)

    assert back["dims"] == (10, 2, 2)
    assert back["histories"] == 100
    assert back["species_z"] == [5, 14, 8]
    for name in ARRAYS:
        assert name in back
    assert back["displacements"].shape == (10, 2, 2)
    assert back["vacancies"].shape == (3, 10, 2, 2)
    assert back["displacements"].sum() > 0
    prov = back["provenance"]
    assert prov["engine"] == "ion"
    assert prov["seed"] == 9
    assert prov["lindhard_version"] == lh.__version__
    assert len(prov["input_sha256"]) == 64
    # the channels add up
    ch = back["electronic_ev"] + back["lattice_ev"] + back["rest_ev"]
    assert np.allclose(ch, back["deposited_ev"], rtol=1e-12)


def test_thread_count_does_not_change_the_export():
    run = lh.Run.from_toml(BASE + VOXEL)
    a = run.run(ions=100, seed=9, threads=1).damage_export()
    b = run.run(ions=100, seed=9, threads=2).damage_export()
    same(a, b)


def test_no_voxel_table_no_export(tmp_path):
    result = lh.Run.from_toml(BASE).run(ions=20, seed=1, threads=1)
    assert result.damage_export() is None
    result.write(tmp_path)
    assert not (tmp_path / "damage_export.csv").exists()


def test_reader_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        lh.read_damage_export(tmp_path / "nope.csv")
    bad = tmp_path / "bad.csv"
    bad.write_text("# lindhard voxel tally v1\n")
    with pytest.raises(lh.InputError):
        lh.read_damage_export(bad)
