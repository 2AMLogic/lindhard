//! A setup downloaded from the simulator page runs the same under
//! `lindhard-cli`: for `page_generated.toml` (written by `web/toml_model.js`,
//! and checked against the generator by `node/toml_model.test.mjs`) the
//! session's report is identical to the report of the CLI's driver,
//! `lindhard_cli::sim::simulate`, for the same seed and ion count.

use lindhard::input::Input;
use lindhard_wasm::Simulation;

const PAGE_TOML: &str = include_str!("page_generated.toml");

/// `examples/b_5keV_si_crystal.toml`: B into a Si substrate transported on the
/// explicit cubic lattice, 7 degrees off the (001) normal.
const CRYSTAL_TOML: &str = include_str!("../../examples/b_5keV_si_crystal.toml");

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

/// The crystal example at `ions`. 150 is in line with the ion counts of the
/// crystal tests in `lindhard/tests/`; it is not a multiple of `chunk_size`
/// (64), so the session's final short chunk is exercised too.
fn crystal_toml(ions: u32) -> String {
    let t = CRYSTAL_TOML.replace("ions = 2000", &format!("ions = {ions}"));
    assert!(t.contains(&format!("ions = {ions}")), "the example changed");
    t
}

/// `toml` with the `[[crystal]]` and `[crystal.thermal]` tables removed: the
/// same beam, target, physics and tally on an amorphous Si substrate.
///
/// Drops every line from a crystal table header to the next table header that
/// is not one. Headers are matched as whole lines, so the `[[crystal]]` the
/// example mentions in its comments is not mistaken for the table itself.
fn without_crystal(toml: &str) -> String {
    let mut out = String::new();
    let mut dropping = false;
    for line in toml.lines() {
        let t = line.trim();
        if t == "[[crystal]]" || t == "[crystal]" || t.starts_with("[crystal.") {
            dropping = true;
        } else if dropping && t.starts_with('[') {
            dropping = false;
        }
        if !dropping {
            out.push_str(line);
            out.push('\n');
        }
    }
    assert!(out.len() < toml.len(), "no crystal table was removed");
    assert!(
        Input::from_toml_str(&out).unwrap().crystal.is_empty(),
        "a crystal table survived the strip"
    );
    out
}

/// A `[[crystal]]` input runs the same under the session as under the CLI's
/// driver. The session needs no crystal-specific path for this: it already
/// builds the engine with `Resolved::bca`, which applies every `[[crystal]]`
/// assignment, and steps it with `Bca::run_range`, which transports the
/// lattice inline in the history.
#[test]
fn crystal_input_gives_the_cli_report() {
    let toml = crystal_toml(150);
    // The lattice is what makes this input different: stripping `[[crystal]]`
    // (the same amorphous Si substrate) gives a different report, so this test
    // cannot pass on a session that silently drops the crystal tables.
    let amorphous = without_crystal(&toml);
    let want = cli_report(&toml, 1);
    assert_ne!(
        want,
        cli_report(&amorphous, 1),
        "the lattice changes nothing"
    );
    // The session applies it too, whatever the step size and the CLI thread
    // count (the determinism-across-thread-counts invariant).
    assert_eq!(session_report(&toml, 1000), want);
    assert_eq!(session_report(&toml, 7), want);
    assert_eq!(cli_report(&toml, 2), want);
    let report: serde_json::Value = serde_json::from_str(&want).unwrap();
    assert_eq!(report["histories"], 150);
}
