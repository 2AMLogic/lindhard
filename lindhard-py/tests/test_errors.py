"""Error mapping to Python exceptions."""

import lindhard as lh
import pytest

from test_config import argon_on_copper

GOOD = """
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


def test_hierarchy():
    assert issubclass(lh.InputError, lh.LindhardError)
    assert issubclass(lh.InputError, ValueError)
    assert issubclass(lh.RunError, lh.LindhardError)
    assert issubclass(lh.RunError, RuntimeError)
    assert not issubclass(lh.RunError, ValueError)


def test_unknown_key_names_it():
    with pytest.raises(lh.InputError, match="bogus"):
        lh.Run.from_toml(GOOD.replace("[run]", "[run]\nbogus = 1"))


def test_not_toml():
    with pytest.raises(lh.InputError):
        lh.Run.from_toml("this is = = not toml")


def test_invalid_value_names_the_field():
    run = lh.Run.from_toml(GOOD.replace("5000.0", "-5.0"))
    with pytest.raises(lh.InputError, match="beam"):
        run.run()
    with pytest.raises(ValueError):
        run.validate()


def test_unknown_element_and_material():
    run = lh.Run.from_toml(GOOD.replace('ion = "B"', 'ion = "Xx"'))
    with pytest.raises(lh.InputError):
        run.validate()
    run = lh.Run.from_toml(GOOD.replace('substrate = "Si"', 'substrate = "Nope"'))
    with pytest.raises(lh.InputError, match="target"):
        run.validate()


def test_bad_choice_names_the_field():
    with pytest.raises(lh.InputError, match="physics.potential"):
        lh.Run(
            lh.Beam("Ar", 1000.0),
            lh.Target(substrate="Cu"),
            lh.Physics(2.0, 1.0, potential="not-a-potential"),
            ions=10,
            seed=1,
        )


def test_bad_energy_key_and_material_type():
    with pytest.raises(lh.InputError, match="e_x_ev"):
        lh.Run(
            lh.Beam("Ar", 1000.0),
            lh.Target(substrate="Cu"),
            lh.Physics(2.0, 1.0, energies={"Cu": {"e_x_ev": 1.0}}),
            ions=1,
            seed=1,
        )
    with pytest.raises(lh.InputError, match="material"):
        lh.Layer(3, 1.0)


def test_zero_threads_and_missing_file(tmp_path):
    with pytest.raises(lh.InputError, match="threads"):
        argon_on_copper().run(threads=0)
    with pytest.raises(FileNotFoundError):
        lh.Run.from_toml_file(tmp_path / "missing.toml")


def test_transport_failure_is_a_run_error(tmp_path):
    # A user stopping table that covers the beam energy but not the energies
    # the ion slows through: accepted, but the first history fails.
    (tmp_path / "t.toml").write_text(
        'provenance = "test table"\nion_z = 5\ntarget_z = 14\n'
        "energy_ev = [1.0e3, 2.0e3]\nstopping_ev_1e15_cm2 = [10.0, 12.0]\n"
    )
    run = lh.Run.from_toml(
        GOOD.replace("5000.0", "1500.0") + '[stopping]\ntables = ["t.toml"]\n', base_dir=tmp_path
    )
    with pytest.warns(UserWarning, match="table"):
        with pytest.raises(lh.RunError, match="history"):
            run.run()
