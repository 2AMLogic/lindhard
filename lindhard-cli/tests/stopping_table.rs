//! `lindhard stopping-table`: files written, the Lindhard-Scharff row against
//! its analytic value, units, batches, compounds, one-point grids, and errors
//! that write nothing.

use std::path::PathBuf;
use std::process::Command;

use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::srim::SrimStopping;
use lindhard::ion::stopping::{ElectronicStopping, Ion};
use lindhard_cli::stopping_table::{CSV_FILE, JSON_FILE};
use serde_json::Value;

fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "lindhard-stopping-table-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(d: &std::path::Path, toml: &str, out: &str) -> std::process::Output {
    std::fs::write(d.join("in.toml"), toml).unwrap();
    Command::new(env!("CARGO_BIN_EXE_lindhard"))
        .arg("stopping-table")
        .arg(d.join("in.toml"))
        .arg("--out")
        .arg(d.join(out))
        .output()
        .unwrap()
}

fn json(d: &std::path::Path, out: &str) -> Value {
    serde_json::from_str(&std::fs::read_to_string(d.join(out).join(JSON_FILE)).unwrap()).unwrap()
}

const BASIC: &str = r#"
targets = ["Si"]
[[ions]]
ion = "B"
[energy]
values_ev = [1000.0, 10000.0]
[units]
stopping = "eV/1e15 atoms/cm^2"
"#;

/// The electronic column is the Lindhard-Scharff value of the library, in the
/// requested unit; the total is the sum; the range is positive and ascending.
#[test]
fn lindhard_scharff_row_matches_the_model() {
    let d = dir("ls");
    let o = run(&d, BASIC, "out");
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let j = json(&d, "out");
    let rows = j["tables"][0]["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 2);
    let ion = Ion::new(5).unwrap();
    for r in rows {
        let e = r["energy_ev"].as_f64().unwrap();
        let want = lindhard::ion::stopping::to_ev_1e15_cm2(
            LindhardScharff::new().stopping(&ion, 14, e).unwrap(),
        );
        let got = r["electronic_stopping"].as_f64().unwrap();
        assert!((got / want - 1.0).abs() < 1e-12, "{got} vs {want}");
        let sum = got + r["nuclear_stopping"].as_f64().unwrap();
        assert!((r["total_stopping"].as_f64().unwrap() / sum - 1.0).abs() < 1e-12);
    }
    let (r0, r1) = (
        rows[0]["csda_range"].as_f64().unwrap(),
        rows[1]["csda_range"].as_f64().unwrap(),
    );
    assert!(0.0 < r0 && r0 < r1);
    assert_eq!(j["units"]["stopping"], "eV/1e15 atoms/cm^2");
    assert!(!j["tables"][0]["models"].as_array().unwrap().is_empty());
    let csv = std::fs::read_to_string(d.join("out").join(CSV_FILE)).unwrap();
    assert_eq!(csv.lines().count(), 3);
    assert!(csv.starts_with("ion,target,energy_ev,electronic_stopping [eV/1e15 atoms/cm^2]"));
    std::fs::remove_dir_all(&d).unwrap();
}

/// Two ions and two targets: four tables in input order, ions outer; a
/// compound (SiO2 by name) and a one-point grid work; output is byte-identical
/// on a rerun.
#[test]
fn batch_compound_and_one_point_grid_are_deterministic() {
    let d = dir("batch");
    let toml = r#"
targets = ["Si", "SiO2"]
[[ions]]
ion = "B"
[[ions]]
ion = "P"
[materials.SiO2]
density_g_cm3 = 2.2
elements = [ { symbol = "Si", atom_fraction = 1.0 }, { symbol = "O", atom_fraction = 2.0 } ]
[energy]
start_ev = 5000.0
points = 1
[range]
enabled = false
"#;
    let o = run(&d, toml, "a");
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(run(&d, toml, "b").status.success());
    let j = json(&d, "a");
    let names: Vec<(String, String)> = j["tables"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| {
            (
                t["ion"]["symbol"].as_str().unwrap().to_string(),
                t["target"]["name"].as_str().unwrap().to_string(),
            )
        })
        .collect();
    let want: Vec<(String, String)> = [("B", "Si"), ("B", "SiO2"), ("P", "Si"), ("P", "SiO2")]
        .iter()
        .map(|(a, b)| (a.to_string(), b.to_string()))
        .collect();
    assert_eq!(names, want);
    let comp = j["tables"][1]["target"]["composition"].as_array().unwrap();
    assert_eq!(comp.len(), 2);
    assert!(j["tables"][0]["rows"][0]["csda_range"].is_null());
    for f in [CSV_FILE, JSON_FILE] {
        assert_eq!(
            std::fs::read(d.join("a").join(f)).unwrap(),
            std::fs::read(d.join("b").join(f)).unwrap(),
            "{f}"
        );
    }
    std::fs::remove_dir_all(&d).unwrap();
}

/// An energy outside the model's hard domain is an error naming the key, and
/// nothing is written.
#[test]
fn out_of_domain_energy_fails_and_writes_nothing() {
    let d = dir("domain");
    let toml = r#"
targets = ["Si"]
[[ions]]
ion = "B"
[physics]
stopping = "bethe-bloch"
[energy]
values_ev = [1.0e8, 1.0e9]
[range]
enabled = false
"#;
    // Bethe-Bloch is meaningful at 100 MeV; lower the grid to make it fail.
    let ok = run(&d, toml, "ok");
    assert!(
        ok.status.success(),
        "{}",
        String::from_utf8_lossy(&ok.stderr)
    );
    let bad = toml.replace("[1.0e8, 1.0e9]", "[1.0, 1.0e9]");
    let o = run(&d, &bad, "bad");
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("energy.values_ev[0]"), "{err}");
    assert!(!d.join("bad").exists());
    std::fs::remove_dir_all(&d).unwrap();
}

/// Unknown keys and a range start above the grid are input errors.
#[test]
fn bad_input_is_rejected_with_the_key() {
    let d = dir("bad");
    let o = run(&d, &BASIC.replace("[units]", "[unitz]"), "x");
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("unitz"));
    let o = run(&d, &format!("{BASIC}\n[range]\nfrom_ev = 5000.0\n"), "x");
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("range.from_ev"));
    let o = run(&d, &BASIC.replace("\"Si\"", "\"Unobtainium\""), "x");
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("targets[0]"));
    std::fs::remove_dir_all(&d).unwrap();
}

/// Synthetic numbers made up for tests; not book values. Mirrors
/// `SRIM_FIXTURE` in `lindhard/tests/stopping_tables.rs` (B in Si).
const SRIM_FIXTURE: &str = r#"
format = "lindhard-srim-coefficients/1"

[[row]]
ion_z = 5
target_z = 14
a1 = 2.0
a2 = 3.0
a3 = 400.0
a4 = 5.0
a5 = 0.01
e_min_kev_per_u = 0.1
e_max_kev_per_u = 1.0e4
citation = "synthetic test data, made up for unit tests"
"#;

/// `physics.stopping = "srim"` resolves and the electronic column is the
/// coefficient file's values (regression: the synthetic run input dropped
/// `stopping.srim_data`, so `srim` could never resolve and the file the user
/// set was reported as missing).
#[test]
fn srim_stopping_resolves_and_uses_the_coefficient_file() {
    let d = dir("srim");
    std::fs::write(d.join("srim.toml"), SRIM_FIXTURE).unwrap();
    let toml = r#"
targets = ["Si"]
[[ions]]
ion = "B"
[physics]
stopping = "srim"
[stopping]
srim_data = "srim.toml"
[energy]
values_ev = [10000.0, 100000.0]
[range]
enabled = false
[units]
stopping = "eV/1e15 atoms/cm^2"
"#;
    let o = run(&d, toml, "out");
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let j = json(&d, "out");
    // The echoed electronic model is the file's, citation and all.
    let srim = j["tables"][0]["models"]
        .as_array()
        .unwrap()
        .iter()
        .find(|m| m["name"] == "srim")
        .expect("the srim model echoed")
        .clone();
    assert!(
        srim["citation"]
            .as_str()
            .unwrap()
            .contains("synthetic test data"),
        "{}",
        srim["citation"]
    );
    // And the row values are that file's coefficients, loaded here directly.
    let model = SrimStopping::from_toml_str(SRIM_FIXTURE).unwrap();
    let ion = Ion::new(5).unwrap();
    for r in j["tables"][0]["rows"].as_array().unwrap() {
        let e = r["energy_ev"].as_f64().unwrap();
        let want = lindhard::ion::stopping::to_ev_1e15_cm2(model.stopping(&ion, 14, e).unwrap());
        let got = r["electronic_stopping"].as_f64().unwrap();
        assert!((got / want - 1.0).abs() < 1e-12, "{got} vs {want}");
    }
    std::fs::remove_dir_all(&d).unwrap();
}

/// A user table for (ion_z in target_z) on [1000, 100000] eV, sampled from
/// the library's own Lindhard-Scharff model: no external data.
fn sampled_table(ion_z: u8, target_z: u8) -> String {
    let ls = LindhardScharff::new();
    let ion = Ion::new(ion_z).unwrap();
    let (mut e, mut s) = (Vec::new(), Vec::new());
    for i in 0..=8 {
        let en = 1000.0 * 100.0f64.powf(f64::from(i) / 8.0);
        e.push(format!("{en:e}"));
        s.push(format!(
            "{:e}",
            lindhard::ion::stopping::to_ev_1e15_cm2(ls.stopping(&ion, target_z, en).unwrap())
        ));
    }
    format!(
        "provenance = \"Lindhard-Scharff model of this repository sampled on a grid (test)\"\n\
         ion_z = {ion_z}\ntarget_z = {target_z}\nenergy_ev = [{}]\nstopping_ev_1e15_cm2 = [{}]\n",
        e.join(", "),
        s.join(", ")
    )
}

/// A user table whose lower bound sits above the synthetic cutoff warns
/// nothing: `physics.primary_cutoff_ev` is not a key of this input, and the
/// resolver's advice is about a slowing down that never happens here
/// (regression: the placeholder fired on most `[stopping] tables` inputs).
#[test]
fn table_lower_bound_does_not_warn_about_the_placeholder_cutoff() {
    let d = dir("cutoff");
    std::fs::write(d.join("b_si.toml"), sampled_table(5, 14)).unwrap();
    let toml = r#"
targets = ["Si"]
[[ions]]
ion = "B"
[stopping]
tables = ["b_si.toml"]
[energy]
values_ev = [1000.0, 10000.0]
[range]
from_ev = 1000.0
[units]
stopping = "eV/1e15 atoms/cm^2"
"#;
    let o = run(&d, toml, "out");
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!err.contains("physics.primary_cutoff_ev"), "{err}");
    std::fs::remove_dir_all(&d).unwrap();
}

/// A user table for a pair the batch never asks for (C in Si, ions = B)
/// warns that the table is unused, in stopping-table terms: it names the
/// ion and the target elements, and neither `physics.follow_recoils`, "run"
/// nor "beam" — those describe a transport run, not this input
/// (regression: the resolver's unused-table warning leaked through the
/// filter with its `run` wording).
#[test]
fn unused_user_table_warns_in_stopping_table_terms() {
    let d = dir("unused");
    std::fs::write(d.join("c_si.toml"), sampled_table(6, 14)).unwrap();
    let toml = r#"
targets = ["Si"]
[[ions]]
ion = "B"
[stopping]
tables = ["c_si.toml"]
[energy]
values_ev = [1000.0, 10000.0]
[units]
stopping = "eV/1e15 atoms/cm^2"
"#;
    let o = run(&d, toml, "out");
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let j = json(&d, "out");
    let w = j["tables"][0]["warnings"].as_array().unwrap();
    let msgs: Vec<&str> = w
        .iter()
        .filter_map(|x| x.as_str())
        .filter(|m| m.starts_with("stopping.tables["))
        .collect();
    assert_eq!(msgs.len(), 1, "{w:?}");
    let msg = msgs[0];
    assert!(
        msg.starts_with("stopping.tables[0]: the table (C->Si, Z1=6, Z2=14)"),
        "{msg}"
    );
    assert!(msg.contains("ion B"), "{msg}");
    assert!(msg.contains("target elements Si"), "{msg}");
    assert!(msg.ends_with("; it is unused"), "{msg}");
    for banned in ["run", "beam", "follow_recoils"] {
        assert!(!msg.contains(banned), "{banned} in: {msg}");
    }
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains(&format!("warning: {msg}")), "{err}");
    std::fs::remove_dir_all(&d).unwrap();
}

/// The same table for the pair it declares (C in Si) is used and emits no
/// warning.
#[test]
fn used_user_table_warns_nothing() {
    let d = dir("used");
    std::fs::write(d.join("c_si.toml"), sampled_table(6, 14)).unwrap();
    let toml = r#"
targets = ["Si"]
[[ions]]
ion = "C"
[stopping]
tables = ["c_si.toml"]
[energy]
values_ev = [1000.0, 10000.0]
[range]
from_ev = 1000.0
[units]
stopping = "eV/1e15 atoms/cm^2"
"#;
    let o = run(&d, toml, "out");
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let j = json(&d, "out");
    let w = j["tables"][0]["warnings"].as_array().unwrap();
    assert!(w.is_empty(), "{w:?}");
    std::fs::remove_dir_all(&d).unwrap();
}
