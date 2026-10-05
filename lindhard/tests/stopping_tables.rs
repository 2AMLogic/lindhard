//! The `[stopping]` input table: field-naming errors, composition with the
//! `[physics]` model, and equivalence with the built-in model for a table
//! generated from it.

use std::path::{Path, PathBuf};

use lindhard::input::{Input, InputError};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
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
