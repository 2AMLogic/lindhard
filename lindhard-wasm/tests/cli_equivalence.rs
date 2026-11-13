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

/// `page_generated.toml` with the single `[beam] energy_ev` replaced by a
/// two-line primary-energy spectrum spanning more than a decade, so a run that
/// ignores the spectrum (transporting every primary at the largest line, which
/// is what the resolver puts in `Beam::energy_ev`) cannot match a run that
/// applies it.
fn spectrum_toml() -> String {
    let t = PAGE_TOML.replace(
        "energy_ev = 50000.0",
        "spectrum = [{ energy_ev = 2000.0, weight = 1.0 }, \
         { energy_ev = 50000.0, weight = 1.0 }]",
    );
    assert!(t.contains("spectrum = ["), "the page input changed");
    t
}

#[test]
fn spectrum_input_gives_the_cli_report() {
    let toml = spectrum_toml();
    // The spectrum is what makes this input different: dropping it (every
    // primary at the largest line) gives a different report, so this test
    // cannot pass on a front end that never reads `Resolved::spectrum`.
    let flat = toml.replace(
        "spectrum = [{ energy_ev = 2000.0, weight = 1.0 }, \
         { energy_ev = 50000.0, weight = 1.0 }]",
        "energy_ev = 50000.0",
    );
    let want = cli_report(&toml, 1);
    assert_ne!(want, cli_report(&flat, 1), "the spectrum changes nothing");
    // The session applies it too, whatever the step size.
    assert_eq!(session_report(&toml, 1000), want);
    assert_eq!(session_report(&toml, 7), want);
}
