//! The native half of the native-versus-WASM parity test.
//!
//! `b_5keV_si.golden.json` holds the native result for 10^4 ions of 5 keV B
//! into Si (seed 1). This test checks the native build against it; the Node
//! test `tests/node/parity.test.mjs` checks the WASM build against the same
//! file, with the same tolerance, so the two builds agree through the file.
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

use lindhard_wasm::{range_profile, Request};

const GOLDEN: &str = include_str!("b_5keV_si.golden.json");

fn golden() -> serde_json::Value {
    serde_json::from_str(GOLDEN).expect("golden file is JSON")
}

fn request(g: &serde_json::Value) -> Request {
    Request {
        ion: g["ion"].as_str().unwrap().to_owned(),
        energy_ev: g["energy_ev"].as_f64().unwrap(),
        target: g["target"].as_str().unwrap().to_owned(),
        ions: u32::try_from(g["ions"].as_u64().unwrap()).unwrap(),
        seed: u32::try_from(g["seed"].as_u64().unwrap()).unwrap(),
    }
}

#[test]
fn native_matches_the_golden_record() {
    let g = golden();
    let p = range_profile(&request(&g)).unwrap();
    let want = |k: &str| g[k].as_f64().unwrap();
    assert!(
        (p.rp_nm - want("rp_nm")).abs() <= want("rp_std_err_nm"),
        "Rp {} vs golden {} (tolerance {})",
        p.rp_nm,
        want("rp_nm"),
        want("rp_std_err_nm")
    );
    assert!(
        (p.drp_nm - want("drp_nm")).abs() <= want("drp_std_err_nm"),
        "dRp {} vs golden {} (tolerance {})",
        p.drp_nm,
        want("drp_nm"),
        want("drp_std_err_nm")
    );
    let stopped = want("stopped");
    let n = want("ions");
    let sigma = (stopped * (1.0 - stopped / n)).sqrt();
    assert!(
        (p.stopped as f64 - stopped).abs() <= sigma,
        "stopped {} vs golden {stopped} (tolerance {sigma})",
        p.stopped
    );
    assert_eq!(p.histories, 10_000);
    assert_eq!(p.counts.iter().sum::<u64>(), p.stopped);
}

#[test]
fn bad_input_is_an_error_not_a_panic() {
    let mut req = request(&golden());
    req.ion = "Xx".into();
    assert!(range_profile(&req).is_err());
    req.ion = "B".into();
    req.target = "Qq".into();
    assert!(range_profile(&req).is_err());
    req.target = "Si".into();
    req.ions = 0;
    assert!(range_profile(&req).is_err());
    req.ions = lindhard_wasm::MAX_IONS + 1;
    assert!(range_profile(&req).is_err());
}

/// Prints the record for the golden file: `cargo test -p lindhard-wasm --test
/// parity record -- --ignored --nocapture`.
#[test]
#[ignore = "prints the numbers for the golden file"]
fn record() {
    let p = range_profile(&request(&golden())).unwrap();
    println!(
        "rp_nm {:?}\ndrp_nm {:?}\nrp_std_err_nm {:?}\ndrp_std_err_nm {:?}\nstopped {}",
        p.rp_nm, p.drp_nm, p.rp_std_err_nm, p.drp_std_err_nm, p.stopped
    );
}
