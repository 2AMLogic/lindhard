//! The native half of the native-versus-WASM parity test.
//!
//! `b_5keV_si.golden.json` holds the native result for the run in
//! `b_5keV_si.toml` (10^4 ions of 5 keV B into Si, seed 1), reduced to `Rp`,
//! `ΔRp` and the stopped count. This test checks the native build against it
//! by running the fixture through [`Simulation`]; the Node test
//! `tests/node/parity.test.mjs` checks the WASM build against the same file,
//! with the same tolerance, so the two builds agree through the file.
//!
//! # Tolerance
//!
//! Rp and ΔRp must agree within **one standard error** of the golden value
//! (`rp_std_err_nm`, `drp_std_err_nm`), and the stopped count within the
//! binomial scatter of the golden count. Why not bit-for-bit: the transport
//! uses `sin`, `cos`, `exp`, `ln` and `powf`, which `wasm32-unknown-unknown`
//! takes from Rust's bundled libm and a native build takes from the platform
//! libm (or, with fused multiply-add on aarch64, rounds differently). Those
//! differ in the last bit, a collision cascade amplifies a last-bit difference
//! until the histories diverge, and the two builds then draw the same random
//! numbers into different, statistically equivalent trajectories. A
//! tolerance of one standard error accepts that and rejects a systematic
//! difference of the size the statistics can resolve. Observed today:
//! the Linux x86-64 native and the WASM build agree to about 1e-15 relative.

use lindhard_wasm::Simulation;

const INPUT_TOML: &str = include_str!("b_5keV_si.toml");
const GOLDEN: &str = include_str!("b_5keV_si.golden.json");

fn golden() -> serde_json::Value {
    serde_json::from_str(GOLDEN).expect("golden file is JSON")
}

/// Run the whole fixture in engine chunks and return the final report.
fn report() -> serde_json::Value {
    let mut sim = Simulation::new(INPUT_TOML, 0, 10).unwrap();
    while !sim.done() {
        sim.step(1000).unwrap();
    }
    serde_json::from_str(&sim.report_json().unwrap()).expect("report is JSON")
}

#[test]
fn native_matches_the_golden_record() {
    let g = golden();
    let report = report();
    let want = |k: &str| g[k].as_f64().unwrap();
    let depth = &report["range"]["depth"];
    let rp_nm = depth["mean"].as_f64().unwrap() * 1e9;
    let drp_nm = depth["std_dev"].as_f64().unwrap() * 1e9;
    assert!(
        (rp_nm - want("rp_nm")).abs() <= want("rp_std_err_nm"),
        "Rp {rp_nm} vs golden {} (tolerance {})",
        want("rp_nm"),
        want("rp_std_err_nm")
    );
    assert!(
        (drp_nm - want("drp_nm")).abs() <= want("drp_std_err_nm"),
        "dRp {drp_nm} vs golden {} (tolerance {})",
        want("drp_nm"),
        want("drp_std_err_nm")
    );
    let stopped = want("stopped");
    let n = want("ions");
    let sigma = (stopped * (1.0 - stopped / n)).sqrt();
    let got = report["range"]["stopped"].as_f64().unwrap();
    assert!(
        (got - stopped).abs() <= sigma,
        "stopped {got} vs golden {stopped} (tolerance {sigma})"
    );
    assert_eq!(report["histories"].as_u64().unwrap(), 10_000);
    // The depth histogram accounts for every stopped primary.
    let hist = &report["range"]["depth_histogram"];
    let binned: u64 = hist["counts"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_u64().unwrap())
        .sum();
    assert_eq!(
        binned + hist["underflow"].as_u64().unwrap() + hist["overflow"].as_u64().unwrap(),
        got as u64
    );
}

#[test]
fn bad_input_is_an_error_not_a_panic() {
    let err = |toml: &str| Simulation::new(toml, 0, 10).unwrap_err();
    // Unknown beam element.
    let ion = INPUT_TOML.replace("ion = \"B\"", "ion = \"Xx\"");
    assert!(err(&ion).contains("unknown element"), "{}", err(&ion));
    // Unknown target element.
    let target = INPUT_TOML.replace("substrate = \"Si\"", "substrate = \"Qq\"");
    assert!(
        err(&target).contains("unknown material"),
        "{}",
        err(&target)
    );
    // Over-large run (the resolver takes `ions = 0`; the session's cap).
    let too_many = INPUT_TOML.replace("ions = 10000", "ions = 2000000");
    assert!(err(&too_many).contains("run.ions"), "{}", err(&too_many));
}

/// Prints the record for the golden file: `cargo test -p lindhard-wasm --test
/// parity record -- --ignored --nocapture`.
#[test]
#[ignore = "prints the numbers for the golden file"]
fn record() {
    let report = report();
    let depth = &report["range"]["depth"];
    println!(
        "rp_nm {:?}\ndrp_nm {:?}\nrp_std_err_nm {:?}\ndrp_std_err_nm {:?}\nstopped {}",
        depth["mean"].as_f64().unwrap() * 1e9,
        depth["std_dev"].as_f64().unwrap() * 1e9,
        depth["mean_std_err"].as_f64().unwrap() * 1e9,
        depth["std_dev_std_err"].as_f64().unwrap() * 1e9,
        report["range"]["stopped"].as_u64().unwrap(),
    );
}
