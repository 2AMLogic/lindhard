//! Tests for the input module.

use std::path::Path;

use super::*;
use crate::elements::element_by_symbol;
use crate::ion::bca::{ElectronicLoss, MeanFreePath};
use crate::ion::potential::{Screening, ScreeningLength};

const B_SI: &str = r#"
[beam]
ion = "B"
energy_ev = 5000.0
tilt_deg = 7.0

[target]
substrate = "Si"

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 10
seed = 1
"#;

fn err(text: &str) -> InputError {
    Input::from_toml_str(text)
        .and_then(|i| i.resolve().map(|_| ()))
        .unwrap_err()
}

fn field(e: &InputError) -> &str {
    match e {
        InputError::Invalid { field, .. } => field,
        InputError::Parse(_) => panic!("expected a value error, got {e}"),
    }
}

#[test]
fn minimal_input_resolves_with_defaults() {
    let i = Input::from_toml_str(B_SI).unwrap();
    let r = i.resolve().unwrap();
    assert_eq!(r.beam.ion.z(), 5);
    assert!((r.beam.polar_rad - 7f64.to_radians()).abs() < 1e-15);
    assert_eq!(r.stack.layers().len(), 1);
    assert_eq!(r.config.mean_free_path, MeanFreePath::Constant);
    assert_eq!(r.config.electronic, ElectronicLoss::NonLocal);
    assert_eq!(r.config.seed, 1);
    assert_eq!(r.screening, Screening::ZblUniversal);
    assert_eq!(r.screening_length, ScreeningLength::Universal);
    assert_eq!(
        r.layers[0]
            .material
            .displacement_energy_ev(14)
            .unwrap()
            .to_bits(),
        15f64.to_bits()
    );
    assert_eq!(
        r.input.physics.screening_length,
        Some(LengthChoice::Universal)
    );
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
}

#[test]
fn echo_round_trips_and_drops_threads() {
    let mut i = Input::from_toml_str(B_SI).unwrap();
    i.run.threads = Some(3);
    let e = i.echo();
    assert_eq!(e.run.threads, None);
    let text = toml::to_string(&e).unwrap();
    assert_eq!(Input::from_toml_str(&text).unwrap(), e);
}

#[test]
fn layers_by_name_and_inline() {
    let text = r#"
[beam]
ion = "As"
energy_ev = 5.0e4

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
thickness_nm = 10.0

[[target.layers]]
material = { elements = [ { symbol = "Si", atom_fraction = 1.0 } ] }
thickness_nm = 5.0

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 1.0
free_path = "energy-dependent"
min_cm_angle_deg = 0.5
stopping = "equipartition-ls-or"
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 10
seed = 2
"#;
    let r = Input::from_toml_str(text).unwrap().resolve().unwrap();
    assert_eq!(r.stack.layers().len(), 3);
    assert_eq!(r.layers[0].source, "SiO2");
    assert_eq!(r.layers[1].source, "inline");
    assert_eq!(r.layers[2].source, "Si");
    // The override reached Si in every layer.
    for l in &r.layers {
        assert_eq!(l.material.displacement_energy_ev(14).unwrap(), 15.0);
    }
    assert_eq!(r.config.electronic, ElectronicLoss::EquipartitionLsOr);
    assert!(matches!(
        r.config.mean_free_path,
        MeanFreePath::EnergyDependent { .. }
    ));
    assert!((r.stack.layers()[1].back_m() - 15e-9).abs() < 1e-20);
    // Echo of a mixed input round-trips through TOML.
    let back = Input::from_toml_str(&toml::to_string(&r.input).unwrap()).unwrap();
    assert_eq!(back, r.input);
}

#[test]
fn weak_collisions_reach_the_engine_and_echo() {
    let r = Input::from_toml_str(&B_SI.replace(
        "primary_cutoff_ev = 5.0",
        "primary_cutoff_ev = 5.0\nweak_collisions = 3",
    ))
    .unwrap()
    .resolve()
    .unwrap();
    assert_eq!(r.config.weak_collisions, 3);
    let echo = toml::to_string(&r.input).unwrap();
    assert!(echo.contains("weak_collisions = 3"), "{echo}");
    // The default is 0 and is not echoed, so older echoes are unchanged.
    let r0 = Input::from_toml_str(B_SI).unwrap().resolve().unwrap();
    assert_eq!(r0.config.weak_collisions, 0);
    assert!(!toml::to_string(&r0.input)
        .unwrap()
        .contains("weak_collisions"));
}

#[test]
fn unknown_key_is_named() {
    let e = err(&B_SI.replace("tilt_deg", "tilt_degrees"));
    assert!(matches!(e, InputError::Parse(_)));
    assert!(e.to_string().contains("tilt_degrees"), "{e}");
}

#[test]
fn missing_displacement_energy_is_named() {
    let e = err(&B_SI.replace("[physics.energies.Si]\ne_d_ev = 15.0\n", ""));
    assert_eq!(field(&e), "target.substrate");
    assert!(e.to_string().contains("e_d_ev"), "{e}");
    assert!(e.to_string().contains("physics.energies.Si"), "{e}");
}

#[test]
fn negative_thickness_is_named() {
    let text = B_SI.replace(
        "[target]\nsubstrate = \"Si\"\n",
        "[[target.layers]]\nmaterial = \"Si\"\nthickness_nm = -3.0\n",
    );
    let e = err(&text);
    assert_eq!(field(&e), "target.layers[0].thickness_nm");
}

#[test]
fn value_errors_name_their_field() {
    let cases = [
        (B_SI.replace("ion = \"B\"", "ion = \"b\""), "beam.ion"),
        (B_SI.replace("5000.0", "-1.0"), "beam.energy_ev"),
        (B_SI.replace("7.0", "90.0"), "beam.tilt_deg"),
        (B_SI.replace("ions = 10", "ions = 0"), "run.ions"),
        (
            B_SI.replace("substrate = \"Si\"", "substrate = \"Unobtainium\""),
            "target.substrate",
        ),
        (
            B_SI.replace("substrate = \"Si\"", "substrate = \"O\""),
            "target.substrate",
        ),
        (
            B_SI.replace("recoil_cutoff_ev = 2.0", "recoil_cutoff_ev = 0.0"),
            "physics.recoil_cutoff_ev",
        ),
        (
            B_SI.replace(
                "primary_cutoff_ev = 5.0",
                "primary_cutoff_ev = 5.0\nmin_cm_angle_deg = 1.0",
            ),
            "physics.min_cm_angle_deg",
        ),
        (
            B_SI.replace(
                "primary_cutoff_ev = 5.0",
                "primary_cutoff_ev = 5.0\nfree_path = \"energy-dependent\"",
            ),
            "physics.min_cm_angle_deg",
        ),
        (
            B_SI.replace(
                "primary_cutoff_ev = 5.0",
                "primary_cutoff_ev = 5.0\nweak_collisions = 4",
            ),
            "physics.weak_collisions",
        ),
        (
            B_SI.replace(
                "primary_cutoff_ev = 5.0",
                "primary_cutoff_ev = 5.0\nweak_collisions = 2\nfree_path = \"energy-dependent\"\nmin_cm_angle_deg = 1.0",
            ),
            "physics.weak_collisions",
        ),
        (
            B_SI.replace("[physics.energies.Si]", "[physics.energies.Ge]"),
            "physics.energies.Ge",
        ),
        (
            B_SI.replace("e_d_ev = 15.0", "e_d_ev = -15.0"),
            "physics.energies.Si.e_d_ev",
        ),
        (
            format!("{B_SI}\n[tally]\ndepth_bin_nm = 0.0\n"),
            "tally.depth_bin_nm",
        ),
        (
            format!("{B_SI}\n[tally]\nlateral_bin_nm = -1.0\n"),
            "tally.lateral_bin_nm",
        ),
        (
            format!("{B_SI}\n[tally]\nlateral_bins = 0\n"),
            "tally.lateral_bins",
        ),
        (
            format!("{B_SI}\n[tally]\nescape_energy_bins = 0\n"),
            "tally.escape_energy_bins",
        ),
        (
            format!("{B_SI}\n[tally]\nescape_polar_bins = 0\n"),
            "tally.escape_polar_bins",
        ),
        (
            format!("{B_SI}\n[tally]\nescape_energy_max_ev = 0.0\n"),
            "tally.escape_energy_max_ev",
        ),
    ];
    for (text, want) in cases {
        let e = err(&text);
        assert_eq!(field(&e), want, "{e}");
    }
}

#[test]
fn unknown_model_name_is_a_parse_error() {
    let e = err(&B_SI.replace(
        "primary_cutoff_ev = 5.0",
        "primary_cutoff_ev = 5.0\npotential = \"zbl2\"",
    ));
    assert!(matches!(e, InputError::Parse(_)));
    assert!(e.to_string().contains("zbl2"), "{e}");
}

#[test]
fn warns_outside_stopping_validity() {
    let text = B_SI.replace(
        "primary_cutoff_ev = 5.0",
        "primary_cutoff_ev = 5.0\nstopping = \"bethe-bloch\"",
    );
    let r = Input::from_toml_str(&text).unwrap().resolve().unwrap();
    assert!(
        r.warnings.iter().any(|w| w.contains("bethe-bloch")),
        "{:?}",
        r.warnings
    );
}

#[test]
fn models_list_every_choice() {
    let r = Input::from_toml_str(B_SI).unwrap().resolve().unwrap();
    let names: Vec<_> = r.models().iter().map(|m| m.name).collect();
    for n in ["zbl-universal", "universal", "lindhard-scharff", "constant"] {
        assert!(names.contains(&n), "{names:?}");
    }
    assert!(r.models().iter().all(|m| !m.citation.is_empty()));
}

// Tuning plumbing. The sets below are synthetic fixtures, not fitted
// coefficients.
const FIX: TuningSet = TuningSet {
    name: "fixture",
    version: 3,
    ions: &["B"],
    energy_range_ev: (1000.0, 10000.0),
    e_s_factors: &[("Si", 1.25), ("Ag", 0.5)],
    provenance: "synthetic test fixture",
};

fn with_tuning(text: &str, name: &str) -> Input {
    let mut i = Input::from_toml_str(text).unwrap();
    i.physics.tuning = name.to_string();
    i
}

fn resolve_fix(i: &Input, sets: &[TuningSet]) -> Result<Resolved, InputError> {
    i.resolve_with_sets(Path::new("."), sets)
}

fn es(r: &Resolved, layer: usize, z: u8) -> f64 {
    r.layers[layer]
        .material
        .surface_binding_energy_ev(z)
        .unwrap()
}

#[test]
fn tuning_none_and_omitted_are_identical() {
    let a = Input::from_toml_str(B_SI).unwrap().resolve().unwrap();
    let b = with_tuning(B_SI, "none").resolve().unwrap();
    assert!(a.tuning.is_none() && b.tuning.is_none());
    assert_eq!(es(&a, 0, 14), es(&b, 0, 14));
    assert_eq!(
        toml::to_string(&a.input).unwrap(),
        toml::to_string(&b.input).unwrap()
    );
    assert!(!toml::to_string(&a.input).unwrap().contains("tuning"));
}

#[test]
fn tuning_multiplies_resolved_default_once() {
    let base = Input::from_toml_str(B_SI).unwrap().resolve().unwrap();
    let r = resolve_fix(&with_tuning(B_SI, "fixture"), &[FIX]).unwrap();
    let t = r.tuning.as_ref().unwrap();
    assert_eq!((t.set.as_str(), t.version), ("fixture", 3));
    assert_eq!(t.components.len(), 1);
    let c = &t.components[0];
    assert_eq!(c.e_s_original_ev, es(&base, 0, 14));
    assert_eq!(c.factor, 1.25);
    assert_eq!(c.e_s_effective_ev, es(&r, 0, 14));
    assert_eq!(c.e_s_effective_ev, c.e_s_original_ev * 1.25);
    assert_eq!(
        toml::to_string(&r.input)
            .unwrap()
            .matches("fixture")
            .count(),
        1
    );
}

#[test]
fn tuning_acts_on_explicit_override_and_each_layer() {
    let text = B_SI.replace("e_d_ev = 15.0", "e_d_ev = 15.0\ne_s_ev = 4.0").replace(
        "[target]\nsubstrate = \"Si\"\n",
        "[target]\nsubstrate = \"Si\"\n[[target.layers]]\nmaterial = \"Si\"\nthickness_nm = 5.0\n",
    );
    let r = resolve_fix(&with_tuning(&text, "fixture"), &[FIX]).unwrap();
    let t = r.tuning.unwrap();
    assert_eq!(t.components.len(), 2);
    for c in &t.components {
        assert_eq!(c.e_s_original_ev, 4.0);
        assert_eq!(c.e_s_effective_ev, 5.0);
    }
}

#[test]
fn tuning_rejections() {
    let fe = |i: &Input, sets: &[TuningSet]| resolve_fix(i, sets).unwrap_err();
    let e = fe(&with_tuning(B_SI, "nope"), &[FIX]);
    assert_eq!(field(&e), "physics.tuning");
    assert!(e.to_string().contains("unknown tuning set"), "{e}");
    // A fixture name is unknown to the shipped registry.
    assert_eq!(
        field(&err(&B_SI.replace("[run]", "[run]\n").replace(
            "recoil_cutoff_ev = 2.0",
            "recoil_cutoff_ev = 2.0\ntuning = \"fixture\""
        ))),
        "physics.tuning"
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let factors: &'static [(&str, f64)] = Box::leak(Box::new([("Si", bad)]));
        let s = TuningSet {
            e_s_factors: factors,
            ..FIX
        };
        assert_eq!(
            field(&fe(&with_tuning(B_SI, "fixture"), &[s])),
            "physics.tuning"
        );
    }
    let unknown_el = TuningSet {
        e_s_factors: &[("Xx", 1.0)],
        ..FIX
    };
    assert!(fe(&with_tuning(B_SI, "fixture"), &[unknown_el])
        .to_string()
        .contains("Xx"));
    let dup = TuningSet {
        e_s_factors: &[("Si", 1.0), ("Si", 2.0)],
        ..FIX
    };
    assert!(fe(&with_tuning(B_SI, "fixture"), &[dup])
        .to_string()
        .contains("twice"));
    // Compound layer.
    let compound = B_SI.replace(
        "[target]\nsubstrate = \"Si\"\n",
        "[materials.SiO2]\ndensity_g_cm3 = 2.2\nelements = [\n  { symbol = \"Si\", atom_fraction = 1.0 },\n  { symbol = \"O\", atom_fraction = 2.0, e_d_ev = 20.0, e_s_ev = 2.0 },\n]\n[target]\nsubstrate = \"SiO2\"\n",
    );
    let e = fe(&with_tuning(&compound, "fixture"), &[FIX]);
    assert!(e.to_string().contains("single-element"), "{e}");
    // Unset energy stays an error, not a default.
    let none_set = TuningSet {
        e_s_factors: &[("Si", 1.0)],
        ..FIX
    };
    let r = resolve_fix(&with_tuning(B_SI, "fixture"), &[none_set]).unwrap();
    assert!(r.warnings.iter().any(|w| w.contains("nothing changed")));
    // A beam the set was not fitted for.
    let e = fe(
        &with_tuning(&B_SI.replace("ion = \"B\"", "ion = \"P\""), "fixture"),
        &[FIX],
    );
    assert_eq!(field(&e), "physics.tuning");
    assert!(e.to_string().contains("beam ion"), "{e}");
    // A target element the set does not list.
    let ge = B_SI
        .replace("substrate = \"Si\"", "substrate = \"Ge\"")
        .replace("[physics.energies.Si]", "[physics.energies.Ge]");
    let e = fe(&with_tuning(&ge, "fixture"), &[FIX]);
    assert!(e.to_string().contains("no factor for Ge"), "{e}");
    // A malformed energy range.
    for range in [(0.0, 1.0), (2.0, 1.0), (1.0, f64::INFINITY)] {
        let s = TuningSet {
            energy_range_ev: range,
            ..FIX
        };
        assert_eq!(
            field(&fe(&with_tuning(B_SI, "fixture"), &[s])),
            "physics.tuning"
        );
    }
}

#[test]
fn tuning_warns_outside_the_fitted_domain() {
    // B_SI: 5 keV (inside 1..10 keV) at 7 deg tilt.
    let r = resolve_fix(&with_tuning(B_SI, "fixture"), &[FIX]).unwrap();
    assert!(r.warnings.iter().any(|w| w.contains("normal incidence")));
    assert!(!r.warnings.iter().any(|w| w.contains("fitted at 1000")));
    let low = B_SI
        .replace("energy_ev = 5000.0", "energy_ev = 500.0")
        .replace("tilt_deg = 7.0", "tilt_deg = 0.0");
    let r = resolve_fix(&with_tuning(&low, "fixture"), &[FIX]).unwrap();
    assert!(r.warnings.iter().any(|w| w.contains("extrapolation")));
    assert!(!r.warnings.iter().any(|w| w.contains("normal incidence")));
    assert_eq!(r.tuning.unwrap().components[0].factor, 1.25);
}

#[test]
fn shipped_sets_are_well_formed() {
    let mut names = std::collections::BTreeSet::new();
    for s in TUNING_SETS {
        assert!(names.insert(s.name), "duplicate set name {}", s.name);
        assert_ne!(s.name, NO_TUNING);
        assert!(!s.ions.is_empty() && !s.e_s_factors.is_empty());
        assert!(!s.provenance.is_empty());
        for &(sym, k) in s.e_s_factors {
            assert!(element_by_symbol(sym).is_some(), "{sym}");
            assert!(k.is_finite() && k > 0.0, "{sym}: {k}");
        }
        // Every listed element resolves for every listed ion, once.
        for ion in s.ions {
            for &(sym, k) in s.e_s_factors {
                let text = format!(
                    "[beam]\nion = \"{ion}\"\nenergy_ev = {}\n[target]\nsubstrate = \"{sym}\"\n\
                     [physics]\nprimary_cutoff_ev = 2.0\nrecoil_cutoff_ev = 1.0\ntuning = \"{}\"\n\
                     [physics.energies.{sym}]\ne_d_ev = 10.0\n\
                     [run]\nions = 1\nseed = 1\n",
                    s.energy_range_ev.0, s.name
                );
                let r = Input::from_toml_str(&text).unwrap().resolve().unwrap();
                let c = &r.tuning.as_ref().unwrap().components[0];
                assert_eq!(c.factor, k);
                assert_eq!(c.e_s_effective_ev, c.e_s_original_ev * k);
                assert_eq!(
                    es(&r, 0, element_by_symbol(sym).unwrap().z),
                    c.e_s_effective_ev
                );
            }
        }
    }
}
