//! Tests of the steppable [`Simulation`] session (native build).
//!
//! The chunk-invariance tests compare byte-identical JSON between different
//! step sizes on the same build. The native-versus-WASM comparison, which is
//! statistical, is in `parity.rs` and `node/parity.test.mjs`.

use lindhard_wasm::{elements_json, Simulation};

/// `examples/as_50keV_si_sio2.toml` (a compound screen layer on a Si
/// substrate) with the recoil cascade followed, at a reduced ion count.
pub const SCREEN_TOML: &str = include_str!("../../examples/as_50keV_si_sio2.toml");

fn input(ions: u32) -> String {
    let t = SCREEN_TOML.replace("ions = 2000", &format!("ions = {ions}"));
    assert!(t.contains(&format!("ions = {ions}")), "example changed");
    if t.contains("follow_recoils") {
        t
    } else {
        t.replace("[physics]\n", "[physics]\nfollow_recoils = true\n")
    }
}

/// Run `ions` histories in steps of `chunks` (the last repeats until done).
fn run(ions: u32, chunks: &[u32]) -> (String, String) {
    let mut sim = Simulation::new(&input(ions), 20, 400).unwrap();
    let mut trajectories = Vec::new();
    let mut i = 0;
    while !sim.done() {
        let n = chunks[i.min(chunks.len() - 1)];
        i += 1;
        sim.step(n).unwrap();
        trajectories.push(sim.take_trajectories_json());
    }
    // Concatenate the histories of every call into one list.
    let mut all = Vec::new();
    for t in trajectories {
        let v: serde_json::Value = serde_json::from_str(&t).unwrap();
        all.extend(v["histories"].as_array().unwrap().iter().cloned());
    }
    (
        sim.report_json().unwrap(),
        serde_json::to_string(&all).unwrap(),
    )
}

#[test]
fn chunking_does_not_change_the_report_or_the_trajectories() {
    let ions = 300;
    let (one_report, one_traj) = run(ions, &[ions]);
    let (cut_report, cut_traj) = run(ions, &[1, 7, 100, ions]);
    assert_eq!(one_report, cut_report, "report depends on the chunking");
    assert_eq!(one_traj, cut_traj, "trajectories depend on the chunking");
    // The cascade was followed and the trajectories are not empty.
    let traj: serde_json::Value = serde_json::from_str(&one_traj).unwrap();
    assert_eq!(traj.as_array().unwrap().len(), 20);
    let report: serde_json::Value = serde_json::from_str(&one_report).unwrap();
    assert_eq!(report["histories"], ions);
}

#[test]
fn trajectories_are_taken_once() {
    let mut sim = Simulation::new(&input(50), 5, 100).unwrap();
    sim.step(10).unwrap();
    let first: serde_json::Value = serde_json::from_str(&sim.take_trajectories_json()).unwrap();
    assert_eq!(first["histories"].as_array().unwrap().len(), 5);
    let again: serde_json::Value = serde_json::from_str(&sim.take_trajectories_json()).unwrap();
    assert!(again["histories"].as_array().unwrap().is_empty());
}

#[test]
fn step_rounds_up_to_whole_chunks_and_stops_at_the_end() {
    let mut sim = Simulation::new(&input(150), 0, 10).unwrap();
    assert!(!sim.done());
    assert_eq!(sim.step(0).unwrap(), 0);
    assert_eq!(sim.step(1).unwrap(), 64);
    assert_eq!(sim.step(64).unwrap(), 128);
    assert_eq!(sim.step(1000).unwrap(), 150);
    assert!(sim.done());
    assert_eq!(sim.step(5).unwrap(), 150);
}

#[test]
fn input_errors_are_text_not_panics() {
    let err = |toml: &str| Simulation::new(toml, 0, 10).unwrap_err();
    // Bad TOML.
    assert!(!err("this is [not toml").is_empty());
    // Unknown element.
    let unknown = input(10).replace("ion = \"As\"", "ion = \"Xx\"");
    assert!(!err(&unknown).is_empty());
    // Energy out of range.
    let bad_energy = input(10).replace("energy_ev = 50000.0", "energy_ev = -5.0");
    assert!(!err(&bad_energy).is_empty());
    // Unsupported sections and over-large runs.
    let too_many = input(10).replace("ions = 10", "ions = 2000000");
    assert!(err(&too_many).contains("run.ions"));
}

#[test]
fn elements_json_is_generated_from_the_library_table() {
    let v: serde_json::Value = serde_json::from_str(&elements_json()).unwrap();
    let rows = v.as_array().unwrap();
    assert_eq!(rows.len(), lindhard::elements::NUM_ELEMENTS);
    let si = rows.iter().find(|r| r["symbol"] == "Si").unwrap();
    let want = lindhard::elements::element_by_symbol("Si").unwrap();
    assert_eq!(si["z"], 14);
    assert_eq!(si["atomic_weight"], want.atomic_weight);
    assert_eq!(
        si["default_displacement_ev"],
        serde_json::json!(want.default_displacement_ev)
    );
}

// ---- native half of the native-versus-WASM session parity ----------------
//
// `as_50keV_si_sio2.golden.json` holds the native `report_json()` summary of
// `SCREEN_TOML` at the golden file's ion count. `node/session.test.mjs` runs
// the same input through the WASM build and checks the same numbers with the
// same tolerance, which is the one of `parity.rs` and for the same reason
// (wasm32 and native take sin/cos/exp/ln/powf from different libms, so the
// two builds draw the same random numbers into statistically equivalent
// histories): the depth mean and standard deviation within one standard error
// of the golden values, and the stopped count within its binomial scatter (at least one ion).

const GOLDEN: &str = include_str!("as_50keV_si_sio2.golden.json");

fn summary(report: &serde_json::Value) -> [f64; 5] {
    let d = &report["range"]["depth"];
    let f = |v: &serde_json::Value| v.as_f64().unwrap();
    [
        f(&d["mean"]),
        f(&d["std_dev"]),
        f(&d["mean_std_err"]),
        f(&d["std_dev_std_err"]),
        f(&report["range"]["stopped"]),
    ]
}

fn golden_report() -> (serde_json::Value, u32, serde_json::Value) {
    let g: serde_json::Value = serde_json::from_str(GOLDEN).unwrap();
    let ions = u32::try_from(g["ions"].as_u64().unwrap()).unwrap();
    let mut sim = Simulation::new(&input(ions), 0, 10).unwrap();
    while !sim.done() {
        sim.step(1000).unwrap();
    }
    let report = serde_json::from_str(&sim.report_json().unwrap()).unwrap();
    (g, ions, report)
}

#[test]
fn native_session_matches_the_golden_record() {
    let (g, ions, report) = golden_report();
    let [mean, sd, _, _, stopped] = summary(&report);
    let want = |k: &str| g[k].as_f64().unwrap();
    assert!((mean - want("depth_mean_m")).abs() <= want("depth_mean_std_err_m"));
    assert!((sd - want("depth_std_dev_m")).abs() <= want("depth_std_dev_std_err_m"));
    let n = f64::from(ions);
    let sigma = (want("stopped") * (1.0 - want("stopped") / n))
        .sqrt()
        .max(1.0);
    assert!((stopped - want("stopped")).abs() <= sigma);
}

/// Prints the numbers for the golden file: `cargo test -p lindhard-wasm --test
/// session record -- --ignored --nocapture`.
#[test]
#[ignore = "prints the numbers for the golden file"]
fn record() {
    let (_, _, report) = golden_report();
    let [mean, sd, mean_se, sd_se, stopped] = summary(&report);
    println!(
        "depth_mean_m {mean:?}\ndepth_std_dev_m {sd:?}\ndepth_mean_std_err_m {mean_se:?}\n\
         depth_std_dev_std_err_m {sd_se:?}\nstopped {stopped}"
    );
}
