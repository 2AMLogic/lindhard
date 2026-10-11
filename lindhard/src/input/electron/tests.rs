//! Tests for the electron input module.

use std::path::PathBuf;

use super::*;
use crate::electron::inelastic::{FitWeighting, MerminFitOptions, PennAlgorithm, MAX_OSCILLATORS};
use crate::electron::secondary::SecondaryModel;
use crate::electron::transport::BoundaryModel;
use crate::tally::{PsfModel, PsfNormalization};

const ELF: &str = "material = \"synthetic\"\nprovenance = \"synthetic test ELF\"\n\
                   energy_ev = [1.0, 10.0, 100.0]\nelf = [0.1, 1.0, 0.1]\n";

fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "lindhard-electron-input-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

const GOOD: &str = r#"
[electron.beam]
energy_ev = 1000.0

[electron.transport]
cutoff_ev = 20.0

[electron.elastic]
potential = "thomas-fermi-yukawa"

[electron.materials.Si]
optical_elf = "elf.toml"

[target]
substrate = "Si"

[run]
histories = 10
seed = 1
"#;

#[test]
fn detects_electron_documents() {
    assert!(ElectronInput::is_electron_toml(GOOD));
    assert!(!ElectronInput::is_electron_toml("[beam]\nion = \"B\"\n"));
}

#[test]
fn resolves_and_echoes_defaults() {
    let d = dir("good");
    std::fs::write(d.join("elf.toml"), ELF).unwrap();
    let i = ElectronInput::from_toml_str(GOOD).unwrap();
    let r = i.resolve_in(&d).unwrap();
    assert_eq!(r.layers.len(), 1);
    assert_eq!(
        r.materials[0].optical_elf.provenance(),
        "synthetic test ELF"
    );
    assert_eq!(r.materials[0].optical_elf_file.sha256.len(), 64);
    assert_eq!(r.input.electron.tables.max_energy_ev, Some(1000.0));
    assert_eq!(r.input.electron.tally.escape_energy_max_ev, Some(1000.0));
    assert_eq!(r.inelastic, PennAlgorithm::SinglePole);
    assert_eq!(*r.table_energy_ev.last().unwrap(), 1000.0);
    // The echo round-trips through TOML.
    let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
    assert_eq!(back, r.input);
}

#[test]
fn data_without_provenance_is_refused() {
    let d = dir("noprov");
    std::fs::write(d.join("elf.toml"), ELF.replace("synthetic test ELF", " ")).unwrap();
    let e = ElectronInput::from_toml_str(GOOD)
        .unwrap()
        .resolve_in(&d)
        .unwrap_err()
        .to_string();
    assert!(
        e.contains("electron.materials.Si.optical_elf") && e.contains("provenance"),
        "{e}"
    );
    // Band parameters without a provenance are refused too.
    std::fs::write(d.join("elf.toml"), ELF).unwrap();
    let text = GOOD.replace(
        "optical_elf = \"elf.toml\"",
        "optical_elf = \"elf.toml\"\nband = { kind = \"metal\", fermi_ev = 5.0, \
         work_function_ev = 4.0, provenance = \"\" }",
    );
    let e = ElectronInput::from_toml_str(&text)
        .unwrap()
        .resolve_in(&d)
        .unwrap_err()
        .to_string();
    assert!(e.contains("electron.materials.Si.band"), "{e}");
}

#[test]
fn metal_band_takes_an_optional_valence_binding() {
    let d = dir("valence-binding");
    std::fs::write(d.join("elf.toml"), ELF).unwrap();
    let with = |band: &str| {
        GOOD.replace(
            "optical_elf = \"elf.toml\"",
            &format!("optical_elf = \"elf.toml\"\nband = {band}"),
        )
    };
    let resolve = |text: &str| ElectronInput::from_toml_str(text).and_then(|i| i.resolve_in(&d));
    // Absent: no binding, and the echo has no such key.
    let r = resolve(&with(
        "{ kind = \"metal\", fermi_ev = 5.0, work_function_ev = 4.0, provenance = \"p\" }",
    ))
    .unwrap();
    assert_eq!(
        r.materials[0].band.as_ref().unwrap().valence_binding_ev(),
        None
    );
    assert!(!r
        .input
        .to_toml_string()
        .unwrap()
        .contains("valence_binding_ev"));
    for band in [
        "{ kind = \"metal\", fermi_ev = 5.0, work_function_ev = 4.0, valence_binding_ev = 7.5, \
         provenance = \"p\" }",
        "{ kind = \"free-electron-metal\", valence_electrons_per_atom = 1.0, \
         work_function_ev = 4.0, valence_binding_ev = 7.5, provenance = \"p\" }",
    ] {
        let r = resolve(&with(band)).unwrap();
        assert_eq!(
            r.materials[0].band.as_ref().unwrap().valence_binding_ev(),
            Some(7.5)
        );
        let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
        assert_eq!(back, r.input);
    }
    // Refused: a negative binding, and any binding on an insulator.
    let e = resolve(&with(
        "{ kind = \"metal\", fermi_ev = 5.0, work_function_ev = 4.0, valence_binding_ev = -1.0, \
         provenance = \"p\" }",
    ))
    .unwrap_err()
    .to_string();
    assert!(e.contains("electron.materials.Si.band"), "{e}");
    assert!(resolve(&with(
        "{ kind = \"insulator\", valence_band_width_ev = 6.0, band_gap_ev = 3.0, \
         affinity_ev = 1.0, valence_binding_ev = 1.0, provenance = \"p\" }",
    ))
    .is_err());
}

#[test]
fn errors_name_the_field() {
    let d = dir("errors");
    std::fs::write(d.join("elf.toml"), ELF).unwrap();
    let bad = |text: &str| {
        ElectronInput::from_toml_str(text)
            .and_then(|i| i.resolve_in(&d))
            .unwrap_err()
            .to_string()
    };
    let e = bad(&GOOD.replace("cutoff_ev = 20.0", "cutoff_ev = -1.0"));
    assert!(e.contains("electron.transport.cutoff_ev"), "{e}");
    let e = bad(&GOOD.replace("[electron.materials.Si]", "[electron.materials.Ge]"));
    assert!(
        e.contains("target.substrate") && e.contains("electron data"),
        "{e}"
    );
    let e = bad(&GOOD.replace(
        "cutoff_ev = 20.0",
        "cutoff_ev = 20.0\nsecondaries = \"kieft-bosch\"",
    ));
    assert!(e.contains("electron.materials.Si.band"), "{e}");
    let e = bad(&GOOD.replace("potential = \"thomas-fermi-yukawa\"", "potential = \"x\""));
    assert!(e.contains("thomas-fermi-yukawa"), "{e}");
    let e = bad(&GOOD.replace(
        "[electron.elastic]",
        "[electron.inelastic]\nmodel = \"penn\"\n[electron.elastic]",
    ));
    assert!(e.contains("electron.inelastic.model"), "{e}");
    let e = bad(&GOOD.replace("energy_ev = 1000.0", "energy_kev = 1.0"));
    assert!(e.contains("energy_kev"), "{e}");
    let e = bad(&GOOD.replace(
        "potential = \"thomas-fermi-yukawa\"",
        "potential = \"thomas-fermi-yukawa\"\n[electron.elastic.correlation_polarization.polarizability]\n",
    ));
    assert!(e.contains("no polarizability for Si"), "{e}");
}

/// `fermi_energy_ev` sets the inelastic tables of materials without a
/// band; with a band the table takes the band's Fermi energy, so a
/// nonzero value is refused instead of counted twice (#241).
#[test]
fn inelastic_fermi_energy_is_refused_with_a_band() {
    let d = dir("fermi-band");
    std::fs::write(d.join("elf.toml"), ELF).unwrap();
    let fermi = |text: &str, ev: &str| {
        text.replace(
            "[electron.elastic]",
            &format!("[electron.inelastic]\nfermi_energy_ev = {ev}\n[electron.elastic]"),
        )
    };
    let r = ElectronInput::from_toml_str(&fermi(GOOD, "1.5"))
        .unwrap()
        .resolve_in(&d)
        .unwrap();
    assert_eq!(r.inelastic_fermi_ev, 1.5);
    let banded = GOOD.replace(
        "optical_elf = \"elf.toml\"",
        "optical_elf = \"elf.toml\"\nband = { kind = \"metal\", fermi_ev = 5.0, \
         work_function_ev = 4.0, provenance = \"synthetic test band\" }",
    );
    let r = ElectronInput::from_toml_str(&fermi(&banded, "0.0"))
        .unwrap()
        .resolve_in(&d)
        .unwrap();
    assert_eq!(r.inelastic_fermi_ev, 0.0);
    let e = ElectronInput::from_toml_str(&fermi(&banded, "1.5"))
        .unwrap()
        .resolve_in(&d)
        .unwrap_err()
        .to_string();
    assert!(
        e.contains("electron.inelastic.fermi_energy_ev") && e.contains("material Si"),
        "{e}"
    );
}

/// `[electron.materials.<name>.mermin_fit]` (#306): absent gives exactly
/// the default fit options and no echo; given, it resolves, echoes its
/// defaults and round-trips; invalid values, unknown keys and use with
/// another inelastic model are refused naming the key.
#[test]
fn mermin_fit_options_parse_validate_and_resolve() {
    let d = dir("mermin-fit");
    std::fs::write(d.join("elf.toml"), ELF).unwrap();
    let mermin = GOOD.replace(
        "[electron.elastic]",
        "[electron.inelastic]\nmodel = \"mermin-melf\"\n[electron.elastic]",
    );
    let with_fit =
        |base: &str, fit: &str| format!("{base}\n[electron.materials.Si.mermin_fit]\n{fit}\n");
    let resolve = |text: &str| ElectronInput::from_toml_str(text).and_then(|i| i.resolve_in(&d));
    // Absent: the defaults, and nothing echoed.
    for text in [GOOD, mermin.as_str()] {
        let r = resolve(text).unwrap();
        assert_eq!(r.materials[0].mermin_fit, MerminFitOptions::default());
        assert!(r.input.electron.materials["Si"].mermin_fit.is_none());
        assert!(!r.input.to_toml_string().unwrap().contains("mermin_fit"));
    }
    // Explicit defaults resolve to the same options, and the echo fills
    // in the floor.
    for fit in [
        "",
        "oscillators = 3\nweighting = \"relative\"\nrelative_floor = 0.01",
    ] {
        let r = resolve(&with_fit(&mermin, fit)).unwrap();
        assert_eq!(r.materials[0].mermin_fit, MerminFitOptions::default());
        let echo = r.input.electron.materials["Si"].mermin_fit.clone().unwrap();
        assert_eq!(
            echo,
            MerminFitSpec {
                oscillators: 3,
                weighting: MerminWeightingChoice::Relative,
                relative_floor: Some(0.01),
            }
        );
        let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
        assert_eq!(back, r.input);
    }
    let r = resolve(&with_fit(&mermin, "oscillators = 6\nrelative_floor = 0.05")).unwrap();
    assert_eq!(
        r.materials[0].mermin_fit,
        MerminFitOptions {
            n_oscillators: 6,
            weighting: FitWeighting::Relative { floor: 0.05 },
            ..MerminFitOptions::default()
        }
    );
    let r = resolve(&with_fit(
        &mermin,
        "oscillators = 16\nweighting = \"uniform\"",
    ))
    .unwrap();
    assert_eq!(r.materials[0].mermin_fit.n_oscillators, MAX_OSCILLATORS);
    assert_eq!(r.materials[0].mermin_fit.weighting, FitWeighting::Uniform);
    let echo = r.input.electron.materials["Si"].mermin_fit.clone().unwrap();
    assert_eq!(echo.relative_floor, None);
    let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
    assert_eq!(back, r.input);
    // Refused, naming the key.
    let bad = |text: &str| resolve(text).unwrap_err().to_string();
    let key = "electron.materials.Si.mermin_fit";
    for (fit, field) in [
        ("oscillators = 0", "oscillators"),
        ("oscillators = 17", "oscillators"),
        ("relative_floor = 0.0", "relative_floor"),
        ("relative_floor = 1.5", "relative_floor"),
        ("relative_floor = nan", "relative_floor"),
        ("relative_floor = -0.1", "relative_floor"),
        (
            "weighting = \"uniform\"\nrelative_floor = 0.01",
            "relative_floor",
        ),
    ] {
        let e = bad(&with_fit(&mermin, fit));
        assert!(e.contains(&format!("{key}.{field}")), "{fit}: {e}");
    }
    let e = bad(&with_fit(&mermin, "weighting = \"absolute\""));
    assert!(e.contains("weighting") && e.contains("absolute"), "{e}");
    let e = bad(&with_fit(&mermin, "start = 1"));
    assert!(e.contains("start"), "{e}");
    let e = bad(&with_fit(&mermin, "oscillators = -1"));
    assert!(e.contains("oscillators"), "{e}");
    // Any other model: the table would be ignored, so it is refused.
    for model in ["penn-single-pole", "penn-full"] {
        let text = GOOD.replace(
            "[electron.elastic]",
            &format!("[electron.inelastic]\nmodel = \"{model}\"\n[electron.elastic]"),
        );
        let e = bad(&with_fit(&text, "oscillators = 6"));
        assert!(e.contains(key) && e.contains("mermin-melf"), "{model}: {e}");
    }
    let e = bad(&with_fit(GOOD, ""));
    assert!(e.contains(key) && e.contains("penn-single-pole"), "{e}");
}

#[test]
fn psf_section_defaults_echo_and_errors_name_fields() {
    let d = dir("psf");
    std::fs::write(d.join("elf.toml"), ELF).unwrap();
    let psf = "\n[electron.tally.psf]\ndepth_lo_nm = 0.0\ndepth_hi_nm = 5.0\nr_min_nm = 0.5\nr_max_nm = 50.0\nbins = 8\n";
    let text = format!("{GOOD}{psf}");
    let r = ElectronInput::from_toml_str(&text)
        .unwrap()
        .resolve_in(&d)
        .unwrap();
    let c = r.tally.psf.unwrap();
    assert_eq!(c.depth_hi_m, 5e-9);
    assert_eq!(c.radial.bins, 8);
    assert_eq!(
        r.psf_fits,
        vec![PsfModel::DoubleGaussian, PsfModel::TripleGaussian]
    );
    assert_eq!(r.psf_normalization, PsfNormalization::SlabTotal);
    let echo = r.input.electron.tally.psf.clone().unwrap();
    assert_eq!(echo.fits, default_psf_fits());
    assert_eq!(echo.normalization, PsfNormalizationSpec::SlabTotal);
    let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
    assert_eq!(back, r.input);
    // No PSF: nothing resolved.
    let r0 = ElectronInput::from_toml_str(GOOD)
        .unwrap()
        .resolve_in(&d)
        .unwrap();
    assert!(r0.tally.psf.is_none() && r0.psf_fits.is_empty());
    // Empty and explicit fit lists, free normalisation.
    let r = ElectronInput::from_toml_str(&text.replace(
        "bins = 8",
        "bins = 8\nfits = [\"triple\", \"triple\"]\nnormalization = \"free\"",
    ))
    .unwrap()
    .resolve_in(&d)
    .unwrap();
    assert_eq!(r.psf_fits, vec![PsfModel::TripleGaussian]);
    assert_eq!(r.psf_normalization, PsfNormalization::Free);
    let r = ElectronInput::from_toml_str(&text.replace("bins = 8", "bins = 8\nfits = []"))
        .unwrap()
        .resolve_in(&d)
        .unwrap();
    assert!(r.psf_fits.is_empty());
    // Errors.
    let bad = |from: &str, to: &str| {
        ElectronInput::from_toml_str(&text.replace(from, to))
            .and_then(|i| i.resolve_in(&d))
            .unwrap_err()
            .to_string()
    };
    let e = bad("depth_hi_nm = 5.0", "depth_hi_nm = 0.0");
    assert!(e.contains("electron.tally.psf.depth_hi_nm"), "{e}");
    let e = bad("r_max_nm = 50.0", "r_max_nm = 0.5");
    assert!(e.contains("electron.tally.psf.r_max_nm"), "{e}");
    let e = bad("r_min_nm = 0.5", "r_min_nm = 0.0");
    assert!(e.contains("electron.tally.psf.r_min_nm"), "{e}");
    let e = bad("bins = 8", "bins = 0");
    assert!(e.contains("electron.tally.psf.bins"), "{e}");
    let e = bad("bins = 8", "bins = 8\nfits = [\"quadruple\"]");
    assert!(e.contains("quadruple") && e.contains("fits"), "{e}");
    let e = bad("bins = 8", "bins = 8\nnormalization = \"x\"");
    assert!(e.contains("normalization"), "{e}");
    let e = bad("bins = 8", "bins = 8\nbogus = 1");
    assert!(e.contains("bogus"), "{e}");
}

#[test]
fn every_option_resolves_and_is_recorded() {
    let d = dir("rich");
    std::fs::write(d.join("elf.toml"), ELF).unwrap();
    let text = r#"
[electron.beam]
energy_ev = 2000.0
tilt_deg = 30.0
azimuth_deg = 10.0

[electron.transport]
cutoff_ev = 2.0
cutoff_reference = "vacuum-level"
escape_rule = "front-only"
max_events = 1000
secondaries = "kieft-bosch"
instantaneous_momentum = false
boundary = "step-barrier"
refraction = false

[electron.elastic]
potential = "thomas-fermi-yukawa"
exchange = true
[electron.elastic.correlation_polarization]
b_pol_squared = 1.0
polarizability.Si = { bohr3 = 40.0, source = "synthetic test value" }
polarizability.O = { bohr3 = 5.0, source = "synthetic test value" }

[electron.inelastic]
model = "mermin-melf"

[electron.materials.Ox]
optical_elf = "elf.toml"
band = { kind = "insulator", valence_band_width_ev = 6.0, band_gap_ev = 3.0, affinity_ev = 1.0, provenance = "synthetic test band" }
phonon = { preset = "sio2-63mev", temperature_k = 300.0 }
polaron = { c_per_nm = 0.1, gamma_per_ev = 0.5, provenance = "synthetic test polaron" }

[electron.materials.Si]
optical_elf = "elf.toml"
band = { kind = "free-electron-metal", valence_electrons_per_atom = 4.0, work_function_ev = 4.0, provenance = "synthetic test band" }

[electron.tally]
cartesian = { x = { lo_nm = 0.0, hi_nm = 10.0, bins = 2 }, y = { lo_nm = -5.0, hi_nm = 5.0, bins = 2 }, z = { lo_nm = -5.0, hi_nm = 5.0, bins = 2 } }
cylindrical = { r = { lo_nm = 0.0, hi_nm = 5.0, bins = 2 }, depth = { lo_nm = 0.0, hi_nm = 10.0, bins = 2 } }

[materials.Ox]
density_g_cm3 = 2.2
elements = [{ symbol = "Si", atom_fraction = 1.0 }, { symbol = "O", atom_fraction = 2.0 }]

[[target.layers]]
material = "Ox"
thickness_nm = 5.0

[target]
substrate = "Si"

[run]
histories = 10
seed = 3
"#;
    let r = ElectronInput::from_toml_str(text)
        .unwrap()
        .resolve_in(&d)
        .unwrap();
    assert_eq!(r.layer_material, vec![0, 1]);
    assert_eq!(r.materials[0].name, "Ox");
    assert!(r.materials[0].channels.phonon.is_some());
    assert!(r.materials[0].channels.polaron.is_some());
    assert!(r.materials[1].band.as_ref().unwrap().fermi_ev() > 0.0);
    assert_eq!(r.inelastic, PennAlgorithm::Mermin);
    assert_eq!(
        r.config.secondaries,
        SecondaryModel::KieftBosch {
            instantaneous_momentum: false,
            momentum_conservation: true
        }
    );
    assert_eq!(
        r.config.boundary,
        BoundaryModel::StepBarrier {
            quantum_transmission: true,
            refraction: false
        }
    );
    assert!(r.tally.cartesian.is_some() && r.tally.cylindrical.is_some());
    assert_eq!(r.elastic.correlation_polarization.len(), 2);
    let desc = r.elastic.corrections_description();
    assert!(desc.contains("exchange") && desc.contains("synthetic test value"));
    let names: Vec<_> = r.models().iter().map(|m| m.name).collect();
    for n in [
        "furness-mccarthy",
        "salvat-2003",
        "mermin-melf",
        "frohlich",
        "ganachaud-mokrani",
    ] {
        assert!(names.contains(&n), "{n}: {names:?}");
    }
    // The echo has the defaulted options filled in and round-trips.
    let t = &r.input.electron.transport;
    assert_eq!(t.momentum_conservation, Some(true));
    assert_eq!(t.quantum_transmission, Some(true));
    let max = r.input.electron.tables.max_energy_ev.unwrap();
    assert!(
        max > 2000.0,
        "beam energy plus the largest inner potential: {max}"
    );
    let back = ElectronInput::from_toml_str(&r.input.to_toml_string().unwrap()).unwrap();
    assert_eq!(back, r.input);
    // A preset fixes its own source.
    let e = ElectronInput::from_toml_str(&text.replace(
        "preset = \"sio2-63mev\", temperature_k = 300.0",
        "preset = \"sio2-63mev\", temperature_k = 300.0, provenance = \"x\"",
    ))
    .unwrap()
    .resolve_in(&d)
    .unwrap_err()
    .to_string();
    assert!(e.contains("electron.materials.Ox.phonon"), "{e}");
}
