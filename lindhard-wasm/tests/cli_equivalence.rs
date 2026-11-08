//! A setup downloaded from the simulator page runs the same under
//! `lindhard-cli`: for `page_generated.toml` (written by `web/toml_model.js`,
//! and checked against the generator by `node/toml_model.test.mjs`) the
//! session's report is identical to the report of the CLI's driver,
//! `lindhard_cli::sim::simulate`, for the same seed and ion count.

use lindhard::input::Input;
use lindhard_wasm::Simulation;

const PAGE_TOML: &str = include_str!("page_generated.toml");

fn cli_report(toml: &str, threads: usize) -> String {
    let resolved = Input::from_toml_str(toml).unwrap().resolve().unwrap();
    let sim = lindhard_cli::sim::simulate(&resolved, Some(threads)).unwrap();
    serde_json::to_string(&sim.report).unwrap()
}

fn session_report(toml: &str, step: u32) -> String {
    let mut sim = Simulation::new(toml, 0, 10).unwrap();
    while !sim.done() {
        sim.step(step).unwrap();
    }
    sim.report_json().unwrap()
}

#[test]
fn page_generated_input_gives_the_cli_report() {
    let want = cli_report(PAGE_TOML, 1);
    // Whatever the step size and the CLI thread count.
    assert_eq!(session_report(PAGE_TOML, 1000), want);
    assert_eq!(session_report(PAGE_TOML, 7), want);
    assert_eq!(cli_report(PAGE_TOML, 2), want);
    let report: serde_json::Value = serde_json::from_str(&want).unwrap();
    assert_eq!(report["histories"], 200);
    // The cascade is followed, so the damage counts are populated.
    assert!(report["damage"].is_object());
}
