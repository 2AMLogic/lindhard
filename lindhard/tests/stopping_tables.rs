//! The `[stopping]` input table: field-naming errors, composition with the
//! `[physics]` model, and equivalence with the built-in model for a table
//! generated from it. Also `[physics] stopping = "none"` (issue #58): parsing,
//! echo, the models listing, its composition with tables, and zero electronic
//! loss in a run.

use std::path::{Path, PathBuf};

use lindhard::input::{Input, InputError, StoppingChoice};
use lindhard::ion::bca::{Bca, ElectronicLoss, SummaryTally};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::ScatteringTable;
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::none::NoStopping;
use lindhard::ion::stopping::{to_ev_1e15_cm2, ElectronicStopping, Ion};

fn scratch(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("lindhard-stopping-tables")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A table for `ion_z` in `target_z` sampled from our own Lindhard-Scharff
/// model on a log grid: no external data.
fn ls_table(ion_z: u8, target_z: u8, lo: f64, hi: f64, extra: &str) -> String {
    let ls = LindhardScharff::new();
    let ion = Ion::new(ion_z).unwrap();
    let n = 60;
    let (mut e, mut s) = (Vec::new(), Vec::new());
    for i in 0..=n {
        let en = lo * (hi / lo).powf(f64::from(i) / f64::from(n));
        e.push(format!("{en:e}"));
        s.push(format!(
            "{:e}",
            to_ev_1e15_cm2(ls.stopping(&ion, target_z, en).unwrap())
        ));
    }
    format!(
        "provenance = \"Lindhard-Scharff model of this repository sampled on a grid (test)\"\n\
         ion_z = {ion_z}\ntarget_z = {target_z}\n{extra}\nenergy_ev = [{}]\n\
         stopping_ev_1e15_cm2 = [{}]\n",
        e.join(", "),
        s.join(", ")
    )
}

const INPUT: &str = r#"
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
ions = 10
seed = 1
"#;

fn with_tables(tables: &[&str]) -> Input {
    let list = tables
        .iter()
        .map(|t| format!("{t:?}"))
        .collect::<Vec<_>>()
        .join(", ");
    Input::from_toml_str(&format!("{INPUT}\n[stopping]\ntables = [{list}]\n")).unwrap()
}

fn field_of(e: InputError) -> (String, String) {
    match e {
        InputError::Invalid { field, message } => (field, message),
        InputError::Parse(m) => panic!("expected a value error, got {m}"),
    }
}

#[test]
fn absent_table_changes_nothing() {
    let i = Input::from_toml_str(INPUT).unwrap();
    assert!(i.stopping.is_empty());
    // Not echoed, so existing inputs and outputs are unchanged.
    assert!(!toml::to_string(&i.echo()).unwrap().contains("[stopping]"));
    assert!(i.resolve().unwrap().stopping_tables.is_empty());
}

#[test]
fn unknown_key_in_stopping_table_is_a_parse_error_naming_it() {
    let text = format!("{INPUT}\n[stopping]\ntable = \"x.toml\"\n");
    let e = Input::from_toml_str(&text).unwrap_err().to_string();
    assert!(e.contains("table"), "{e}");
}

#[test]
fn errors_name_the_field() {
    let d = scratch("errors");
    let write = |name: &str, text: &str| {
        std::fs::write(d.join(name), text).unwrap();
    };
    let good = ls_table(5, 14, 1.0, 1.0e4, "");

    // Missing file.
    let (f, m) = field_of(with_tables(&["nope.toml"]).resolve_in(&d).unwrap_err());
    assert_eq!(f, "stopping.tables[0]");
    assert!(m.contains("nope.toml"), "{m}");

    // Invalid TOML.
    write("bad.toml", "this is = = not toml");
    let (f, _) = field_of(with_tables(&["bad.toml"]).resolve_in(&d).unwrap_err());
    assert_eq!(f, "stopping.tables[0]");

    // Missing provenance.
    let no_prov = good.lines().skip(1).collect::<Vec<_>>().join("\n");
    write("noprov.toml", &no_prov);
    let (f, m) = field_of(with_tables(&["noprov.toml"]).resolve_in(&d).unwrap_err());
    assert_eq!(f, "stopping.tables[0]");
    assert!(m.contains("provenance"), "{m}");

    // Ion mass differing from the beam ion's.
    write(
        "mass.toml",
        &ls_table(5, 14, 1.0, 1.0e4, "ion_mass_amu = 11.0093"),
    );
    let (f, m) = field_of(with_tables(&["mass.toml"]).resolve_in(&d).unwrap_err());
    assert_eq!(f, "stopping.tables[0]");
    assert!(m.contains("mass"), "{m}");

    // Not covering the beam energy.
    write("range.toml", &ls_table(5, 14, 1.0, 1.0e3, ""));
    let (f, m) = field_of(with_tables(&["range.toml"]).resolve_in(&d).unwrap_err());
    assert_eq!(f, "stopping.tables[0]");
    assert!(m.contains("beam energy"), "{m}");

    // Two tables for one pair.
    write("good.toml", &good);
    write("dup.toml", &good);
    let (f, m) = field_of(
        with_tables(&["good.toml", "dup.toml"])
            .resolve_in(&d)
            .unwrap_err(),
    );
    assert_eq!(f, "stopping.tables[1]");
    assert!(m.contains("stopping.tables[0]"), "{m}");

    // Equipartition ignores the model, so tables would be ignored.
    let mut i = with_tables(&["good.toml"]);
    i.physics.stopping = lindhard::input::StoppingChoice::EquipartitionLsOr;
    let (f, _) = field_of(i.resolve_in(&d).unwrap_err());
    assert_eq!(f, "stopping.tables");
}

#[test]
fn relative_paths_resolve_against_the_base_directory_and_are_recorded() {
    let d = scratch("paths");
    std::fs::create_dir_all(d.join("sub")).unwrap();
    let text = ls_table(5, 14, 1.0, 1.0e4, "");
    std::fs::write(d.join("sub/b.toml"), &text).unwrap();
    let r = with_tables(&["sub/b.toml"]).resolve_in(&d).unwrap();
    let l = &r.stopping_tables[0];
    assert_eq!(l.path, "sub/b.toml");
    assert!(l.resolved_path.is_absolute());
    assert_eq!(l.sha256.len(), 64);
    assert!(l.sha256.bytes().all(|b| b.is_ascii_hexdigit()));
    // Same bytes, same hash; different bytes, different hash.
    std::fs::write(d.join("sub/c.toml"), format!("{text}\n# edited\n")).unwrap();
    let r2 = with_tables(&["sub/c.toml"]).resolve_in(&d).unwrap();
    assert_ne!(l.sha256, r2.stopping_tables[0].sha256);
    // Listed as a model with its provenance.
    let m = r
        .models()
        .into_iter()
        .find(|m| m.name == "user-table")
        .expect("table model listed");
    assert!(m.citation.contains("sub/b.toml") && m.citation.contains("Lindhard-Scharff"));
}

#[test]
fn composite_uses_the_table_for_its_pair_and_the_model_otherwise() {
    let d = scratch("composite");
    let ls = LindhardScharff::new();
    let ion = Ion::new(5).unwrap();
    std::fs::write(d.join("t.toml"), ls_table(5, 14, 1.0, 1.0e4, "")).unwrap();
    let r = with_tables(&["t.toml"]).resolve_in(&d).unwrap();
    let model = r.stopping_model();
    let at_node = model.stopping(&ion, 14, 1.0e3).unwrap();
    let want = ls.stopping(&ion, 14, 1.0e3).unwrap();
    assert!((at_node / want - 1.0).abs() < 1e-3);
    // Other target: the fallback, exactly.
    assert_eq!(
        model.stopping(&ion, 8, 1.0e3).unwrap(),
        ls.stopping(&ion, 8, 1.0e3).unwrap()
    );
    // Out of the table's range is an error, not the fallback.
    assert!(model.stopping(&ion, 14, 5.0e4).is_err());
}

/// With `follow_recoils` on, a table whose ion is a target element also
/// serves that element's recoils, so it is checked against them up front.
#[test]
fn recoil_species_table_starting_above_the_recoil_cutoff_is_an_error() {
    // Judge's reproduction: B into Si with a Si->Si table from 100 eV; the
    // run used to stop at history 0 ("56.3 eV outside table range").
    let d = scratch("recoil-range");
    std::fs::write(d.join("b.toml"), ls_table(5, 14, 1.0, 1.0e4, "")).unwrap();
    std::fs::write(d.join("si.toml"), ls_table(14, 14, 100.0, 1.0e4, "")).unwrap();
    let (f, m) = field_of(
        with_tables(&["b.toml", "si.toml"])
            .resolve_in(&d)
            .unwrap_err(),
    );
    assert_eq!(f, "stopping.tables[1]");
    assert!(
        m.contains("recoil_cutoff_ev") && m.contains("Si recoils"),
        "{m}"
    );

    // Down to the cutoff it is accepted, and serves the recoils.
    std::fs::write(d.join("si.toml"), ls_table(14, 14, 2.0, 1.0e4, "")).unwrap();
    let r = with_tables(&["b.toml", "si.toml"]).resolve_in(&d).unwrap();
    assert!(
        r.warnings.iter().all(|w| !w.contains("stopping.tables")),
        "{:?}",
        r.warnings
    );
    let si = Ion::new(14).unwrap();
    assert!(r.stopping_model().stopping(&si, 14, 2.0).is_ok());
}

#[test]
fn recoil_species_table_ending_below_the_largest_transfer_warns() {
    let d = scratch("recoil-top");
    // B (10.81 u) on Si (28.085 u) at 5 keV can hand a Si atom ~3.98 keV.
    std::fs::write(d.join("si.toml"), ls_table(14, 14, 1.0, 1.0e3, "")).unwrap();
    let r = with_tables(&["si.toml"]).resolve_in(&d).unwrap();
    assert!(
        r.warnings
            .iter()
            .any(|w| w.contains("stopping.tables[0]") && w.contains("largest energy")),
        "{:?}",
        r.warnings
    );
}

#[test]
fn isotopic_self_ion_table_is_rejected_when_recoils_are_followed() {
    // Judge's reproduction: a Si beam of 27.9769 u into Si with a Si->Si
    // table for 27.9769 u. Si recoils have the standard weight, so the run
    // used to stop at history 0 ("asked for 28.085 u").
    let d = scratch("recoil-mass");
    let si_input = |follow: bool| {
        let text = INPUT
            .replace("ion = \"B\"", "ion = \"Si\"\nmass_amu = 27.9769")
            .replace(
                "recoil_cutoff_ev = 2.0",
                &format!("recoil_cutoff_ev = 2.0\nfollow_recoils = {follow}"),
            );
        Input::from_toml_str(&format!("{text}\n[stopping]\ntables = [\"si.toml\"]\n")).unwrap()
    };
    std::fs::write(
        d.join("si.toml"),
        ls_table(14, 14, 1.0, 1.0e4, "ion_mass_amu = 27.9769"),
    )
    .unwrap();
    let (f, m) = field_of(si_input(true).resolve_in(&d).unwrap_err());
    assert_eq!(f, "stopping.tables[0]");
    assert!(
        m.contains("recoils") && m.contains("standard atomic weight"),
        "{m}"
    );
    // A standard-weight table then fails the beam check: no table serves
    // both, and this is reported before the run, not at history 0.
    std::fs::write(d.join("si.toml"), ls_table(14, 14, 1.0, 1.0e4, "")).unwrap();
    let (f, m) = field_of(si_input(true).resolve_in(&d).unwrap_err());
    assert_eq!(f, "stopping.tables[0]");
    assert!(m.contains("beam ion has mass"), "{m}");
    // Without recoils the isotopic table serves the beam alone.
    std::fs::write(
        d.join("si.toml"),
        ls_table(14, 14, 1.0, 1.0e4, "ion_mass_amu = 27.9769"),
    )
    .unwrap();
    assert!(si_input(false).resolve_in(&d).is_ok());
}

#[test]
fn recoil_species_table_without_followed_recoils_is_unused() {
    let d = scratch("recoil-unused");
    // Starts above the recoil cutoff, but it is never queried.
    std::fs::write(d.join("si.toml"), ls_table(14, 14, 100.0, 1.0e4, "")).unwrap();
    let mut i = with_tables(&["si.toml"]);
    i.physics.follow_recoils = false;
    let r = i.resolve_in(&d).unwrap();
    let w = r
        .warnings
        .iter()
        .find(|w| w.contains("stopping.tables[0]"))
        .expect("unused warning");
    assert!(
        w.contains("unused") && w.contains("Si->Si") && w.contains("target elements Si"),
        "{w}"
    );
    assert!(!w.contains('{'), "element symbols, not a set: {w}");
}

fn with_stopping(choice: &str) -> Input {
    Input::from_toml_str(&INPUT.replace(
        "[physics]\n",
        &format!("[physics]\nstopping = \"{choice}\"\n"),
    ))
    .unwrap()
}

#[test]
fn stopping_key_defaults_to_lindhard_scharff() {
    let i = Input::from_toml_str(INPUT).unwrap();
    assert_eq!(i.physics.stopping, StoppingChoice::LindhardScharff);
}

#[test]
fn none_parses_resolves_and_is_echoed_and_listed() {
    let i = with_stopping("none");
    assert_eq!(i.physics.stopping, StoppingChoice::None);
    let r = i.resolve().unwrap();
    // Nonlocal channel with a zero model: no local loss either.
    assert_eq!(r.config.electronic, ElectronicLoss::NonLocal);
    // No energy range, so no advisory out-of-range warning, even far above
    // any model's range.
    assert!(
        r.warnings.iter().all(|w| !w.contains("validity")),
        "{:?}",
        r.warnings
    );
    let far = Input::from_toml_str(
        &INPUT
            .replace("energy_ev = 5000.0", "energy_ev = 5.0e9")
            .replace("[physics]\n", "[physics]\nstopping = \"none\"\n"),
    )
    .unwrap()
    .resolve()
    .unwrap();
    assert!(
        far.warnings.iter().all(|w| !w.contains("validity")),
        "{:?}",
        far.warnings
    );
    // Echoed, and round-trips.
    let echo = toml::to_string(&r.input).unwrap();
    assert!(echo.contains("stopping = \"none\""), "{echo}");
    assert_eq!(Input::from_toml_str(&echo).unwrap(), r.input);
    // Listed as the electronic stopping model; nothing for Bragg to sum.
    let models = r.models();
    let es: Vec<_> = models
        .iter()
        .filter(|m| m.role == "electronic stopping")
        .collect();
    assert_eq!(es.len(), 1);
    assert_eq!(es[0].name, "none");
    assert!(models.iter().all(|m| m.name != "bragg-additivity"));
    // The model is zero.
    let model = r.stopping_model();
    assert_eq!(model.name(), "none");
    assert_eq!(
        model.stopping(&Ion::new(5).unwrap(), 14, 5.0e3).unwrap(),
        0.0
    );
}

#[test]
fn unknown_stopping_choice_is_a_parse_error() {
    let e = Input::from_toml_str(&INPUT.replace("[physics]\n", "[physics]\nstopping = \"off\"\n"))
        .unwrap_err()
        .to_string();
    assert!(e.contains("off") || e.contains("stopping"), "{e}");
}

/// `none` with tables: allowed. The tables serve their pairs, every other
/// pair has zero electronic stopping.
#[test]
fn none_with_tables_uses_the_tables_and_zero_elsewhere() {
    let d = scratch("none-tables");
    let ls = LindhardScharff::new();
    let ion = Ion::new(5).unwrap();
    std::fs::write(d.join("t.toml"), ls_table(5, 14, 1.0, 1.0e4, "")).unwrap();
    let mut i = with_tables(&["t.toml"]);
    i.physics.stopping = StoppingChoice::None;
    let r = i.resolve_in(&d).unwrap();
    let model = r.stopping_model();
    let at_node = model.stopping(&ion, 14, 1.0e3).unwrap();
    assert!((at_node / ls.stopping(&ion, 14, 1.0e3).unwrap() - 1.0).abs() < 1e-3);
    assert_eq!(model.stopping(&ion, 8, 1.0e3).unwrap(), 0.0);
    // Out of the table's range is still an error, not zero.
    assert!(model.stopping(&ion, 14, 5.0e4).is_err());
    let models = r.models();
    assert!(models.iter().any(|m| m.name == "none"));
    assert!(models.iter().any(|m| m.name == "user-table"));
    assert!(models.iter().any(|m| m.name == "bragg-additivity"));
}

/// Si 5 keV -> Si with `stopping = "none"` through the input path: zero
/// electronic loss, and bit-identical to the engine driven directly with
/// [`NoStopping`], the model of the level-1 `range.si5k_si.*` checks.
#[test]
fn none_run_has_zero_electronic_loss_and_matches_the_library_model() {
    let text = INPUT
        .replace("ion = \"B\"", "ion = \"Si\"")
        .replace(
            "recoil_cutoff_ev = 2.0",
            "recoil_cutoff_ev = 2.0\nfollow_recoils = false\nstopping = \"none\"",
        )
        .replace("ions = 10", "ions = 300");
    let r = Input::from_toml_str(&text).unwrap().resolve().unwrap();
    let pot = Potential::new(r.screening, 14.0, 14.0).with_length(r.screening_length);
    assert_eq!(r.screening, Screening::ZblUniversal);
    let table = ScatteringTable::build(&pot, &r.table_spec);
    let run = |stopping: &(dyn ElectronicStopping + Sync)| {
        Bca::new(r.beam, &r.stack, r.config, stopping, &table)
            .unwrap()
            .run(|| SummaryTally::new(1e-9, 100))
            .unwrap()
    };
    let via_input = run(&*r.stopping_model());
    let direct = run(&NoStopping);
    assert_eq!(via_input, direct);
    assert_eq!(via_input.histories, 300);
    assert_eq!(via_input.budget.electronic_nonlocal, 0.0);
    assert_eq!(via_input.budget.electronic_local, 0.0);
    // And Lindhard-Scharff does lose energy to electrons on the same problem.
    let ls = run(&LindhardScharff::new());
    assert!(ls.budget.electronic_nonlocal > 0.0);
}
