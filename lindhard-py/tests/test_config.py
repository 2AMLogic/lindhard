"""Configuration classes and the TOML round trip."""

import lindhard as lh
import pytest


def argon_on_copper() -> lh.Run:
    return lh.Run(
        beam=lh.Beam("Ar", 1000.0),
        target=lh.Target(substrate="Cu"),
        physics=lh.Physics(
            primary_cutoff_ev=2.0, recoil_cutoff_ev=1.0, potential="kr-c"
        ),
        ions=100,
        seed=1,
        tally=lh.Tally(depth_bin_nm=0.1, depth_bins=100),
    )


def test_classes_round_trip_through_toml():
    run = argon_on_copper()
    text = run.to_toml()
    back = lh.Run.from_toml(text)
    assert back.to_toml() == text
    assert back.beam.ion == "Ar"
    assert back.beam.energy_ev == 1000.0
    assert back.physics.potential == "kr-c"
    assert back.target.substrate == "Cu"
    assert (back.ions, back.seed, back.threads) == (100, 1, None)
    assert back.tally.depth_bins == 100


def test_example_inputs_round_trip(example):
    run = lh.Run.from_toml_file(example)
    text = run.to_toml()
    again = lh.Run.from_toml(text)
    assert again.to_toml() == text
    assert again.validate() == run.validate()
    assert (again.beam.ion, again.ions, again.seed) == (
        run.beam.ion,
        run.ions,
        run.seed,
    )


def test_layers_and_inline_materials():
    sio2 = lh.Material(
        [
            lh.Element("Si", atom_fraction=1.0),
            lh.Element("O", atom_fraction=2.0, e_d_ev=20.0, e_s_ev=2.0),
        ],
        density_g_cm3=2.2,
        name="SiO2",
    )
    run = lh.Run(
        lh.Beam("As", 50e3, tilt_deg=7.0),
        lh.Target([lh.Layer("SiO2", 10.0), lh.Layer(sio2, 5.0)], substrate="Si"),
        lh.Physics(
            5.0, 2.0, energies={"Si": {"e_d_ev": 15.0}}
        ),
        ions=10,
        seed=3,
        materials={"SiO2": sio2},
    )
    back = lh.Run.from_toml(run.to_toml())
    layers = back.target.layers
    assert layers[0].material == "SiO2" and layers[0].thickness_nm == 10.0
    inline = layers[1].material
    assert isinstance(inline, lh.Material)
    assert [e.symbol for e in inline.elements] == ["Si", "O"]
    assert inline.density_g_cm3 == 2.2
    assert back.materials["SiO2"].name == "SiO2"
    assert back.physics.energies == {"Si": {"e_d_ev": 15.0}}
    assert back.validate() == []


def test_properties_are_copies_and_assignable():
    run = argon_on_copper()
    beam = run.beam
    beam.energy_ev = 2000.0
    assert run.beam.energy_ev == 1000.0
    run.beam = beam
    assert run.beam.energy_ev == 2000.0
    run.seed = 9
    run.threads = 2
    assert (run.seed, run.threads) == (9, 2)
    assert "seed = 9" in run.to_toml()
    # threads is not part of the schema's echoed input but is written back.
    assert lh.Run.from_toml(run.to_toml()).threads == 2


def test_validate_returns_warnings_and_accepts_good_input():
    assert argon_on_copper().validate() == []


def test_beam_divergence_round_trips_through_toml():
    run = argon_on_copper()
    run.beam = lh.Beam("Ar", 1000.0, divergence_model="gaussian", divergence_deg=0.5)
    text = run.to_toml()
    assert "[beam.divergence]" in text
    back = lh.Run.from_toml(text)
    assert back.to_toml() == text
    assert (back.beam.divergence_model, back.beam.divergence_deg) == ("gaussian", 0.5)
    assert argon_on_copper().beam.divergence_model is None


def test_invalid_beam_divergence_names_the_field():
    run = argon_on_copper()
    run.beam = lh.Beam(
        "Ar", 1000.0, divergence_model="uniform-cone", divergence_deg=-1.0
    )
    with pytest.raises(Exception, match="beam.divergence.half_angle_deg"):
        run.validate()


def test_unknown_divergence_model_is_rejected():
    run = argon_on_copper()
    bad = lh.Beam("Ar", 1000.0, divergence_model="uniform_cone", divergence_deg=0.5)
    with pytest.raises(lh.InputError, match="beam.divergence_model.*uniform_cone"):
        run.beam = bad
    assert run.beam.divergence_model is None


def test_incomplete_divergence_pair_is_rejected():
    run = argon_on_copper()
    width_only = lh.Beam("Ar", 1000.0, divergence_deg=0.5)
    with pytest.raises(lh.InputError, match="beam.divergence_model"):
        run.beam = width_only
    model_only = lh.Beam("Ar", 1000.0, divergence_model="gaussian")
    with pytest.raises(lh.InputError, match="beam.divergence_deg"):
        run.beam = model_only
    assert run.beam.divergence_model is None
