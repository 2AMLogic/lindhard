//! Level-1 check of `sigma_el` and `sigma_tr1` against published partial-wave
//! values (issue #93; fixture `tests/data/elastic_reference.toml`).
//!
//! **No reference value is committed yet.** The intended source (Jablonski,
//! Salvat & Powell 2004) was read, but tabulates none of these cross sections;
//! its figures cover H, Al, Ni, Ag, Au and Cm (plus some gases), not C, Si or
//! Cu. The DHFS screening
//! coefficients for C, Si and Au are not in the tree either (#130), so the
//! harness below is in place but every case is skipped and reported as such.
//! A skipped case is not a pass: the comparison test reports `checked N of 16`
//! and `INCOMPLETE` while any value is absent, and the `#[ignore]`d test
//! `complete_validation_requires_all_sixteen_values` fails until all 16 are
//! checked (run it with `--ignored`). Filling a case in the fixture activates
//! it; the tolerances are read from the fixture and pinned here so they cannot
//! be loosened without touching this file. The fixture is validated strictly
//! (`parse_fixture`): unknown keys (for example an unsupported
//! `correlation_polarization` setting), blank citations, non-positive or
//! non-finite values, and duplicate or missing cases are rejected.

use lindhard::electron::elastic::corrections::{solve_corrected, Corrections};
use lindhard::electron::elastic::{SalvatDhfs, SolverOptions};
use serde::Deserialize;

const TOL_1KEV: f64 = 0.05;
const TOL_10KEV: f64 = 0.02;

#[derive(Deserialize)]
struct Tolerance {
    rel_at_1000_ev: f64,
    rel_at_10000_ev: f64,
}

const EXPECTED_Z: [u32; 4] = [6, 14, 29, 79];
const EXPECTED_E: [f64; 2] = [1000.0, 10000.0];
const EXPECTED_VALUES: usize = 16;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    z: u32,
    energy_ev: f64,
    sigma_el_a02: Option<f64>,
    sigma_tr1_a02: Option<f64>,
    sigma_el_cite: Option<String>,
    sigma_tr1_cite: Option<String>,
    #[serde(default)]
    exchange: bool,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Reference {
    #[allow(dead_code)]
    source: String,
    #[allow(dead_code)]
    status: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Fixture {
    tolerance: Tolerance,
    #[allow(dead_code)]
    reference: Reference,
    case: Vec<Case>,
}

/// Parses and strictly validates fixture text.
fn parse_fixture(text: &str) -> Result<Fixture, String> {
    let f: Fixture = toml::from_str(text).map_err(|e| format!("parse: {e}"))?;
    if f.tolerance.rel_at_1000_ev != TOL_1KEV || f.tolerance.rel_at_10000_ev != TOL_10KEV {
        return Err("tolerances differ from the pinned 5 % / 2 %".into());
    }
    for z in EXPECTED_Z {
        for e in EXPECTED_E {
            let n = f
                .case
                .iter()
                .filter(|c| c.z == z && c.energy_ev == e)
                .count();
            if n != 1 {
                return Err(format!("Z={z} E={e}: {n} cases, expected exactly 1"));
            }
        }
    }
    if f.case.len() != EXPECTED_Z.len() * EXPECTED_E.len() {
        return Err(format!("{} cases, expected 8", f.case.len()));
    }
    for c in &f.case {
        for (name, v, cite) in [
            ("sigma_el", c.sigma_el_a02, &c.sigma_el_cite),
            ("sigma_tr1", c.sigma_tr1_a02, &c.sigma_tr1_cite),
        ] {
            match (v, cite) {
                (None, None) => {}
                (Some(x), Some(t)) => {
                    if !(x.is_finite() && x > 0.0) {
                        return Err(format!(
                            "Z={} E={} {name}: invalid value {x}",
                            c.z, c.energy_ev
                        ));
                    }
                    if t.trim().is_empty() {
                        return Err(format!(
                            "Z={} E={} {name}: blank citation",
                            c.z, c.energy_ev
                        ));
                    }
                }
                _ => {
                    return Err(format!(
                        "Z={} E={} {name}: value and citation must come together",
                        c.z, c.energy_ev
                    ))
                }
            }
        }
        // Both observables or neither: a half-populated case is not a case.
        if c.sigma_el_a02.is_some() != c.sigma_tr1_a02.is_some() {
            return Err(format!(
                "Z={} E={}: exactly one of sigma_el / sigma_tr1 present",
                c.z, c.energy_ev
            ));
        }
    }
    Ok(f)
}

fn fixture() -> Fixture {
    parse_fixture(include_str!("data/elastic_reference.toml")).expect("fixture is valid")
}

#[test]
fn fixture_covers_the_required_cases_with_the_required_tolerances() {
    fixture();
}

#[cfg(test)]
mod malformed {
    use super::*;

    fn good() -> String {
        include_str!("data/elastic_reference.toml").to_string()
    }

    /// Replaces the first case (Z=6, 1 keV) block with `block`.
    fn with_first_case(block: &str) -> String {
        let t = good();
        let start = t.find("[[case]]").unwrap();
        let end = start + t[start..].find("\n\n").unwrap();
        format!("{}{}{}", &t[..start], block, &t[end..])
    }

    const HEAD: &str = "[[case]]\nz = 6\nenergy_ev = 1000.0\n";

    fn err(text: String) -> String {
        match parse_fixture(&text) {
            Ok(_) => panic!("malformed fixture accepted"),
            Err(e) => e,
        }
    }

    #[test]
    fn rejects_blank_citation() {
        let b = format!("{HEAD}sigma_el_a02 = 1.0\nsigma_el_cite = \"  \"\nsigma_tr1_a02 = 1.0\nsigma_tr1_cite = \"x\"");
        assert!(err(with_first_case(&b)).contains("blank citation"));
    }

    #[test]
    fn rejects_value_without_citation() {
        let b = format!("{HEAD}sigma_el_a02 = 1.0\nsigma_tr1_a02 = 1.0\nsigma_tr1_cite = \"x\"");
        assert!(err(with_first_case(&b)).contains("together"));
    }

    #[test]
    fn rejects_missing_observable() {
        let b = format!("{HEAD}sigma_el_a02 = 1.0\nsigma_el_cite = \"x\"");
        assert!(err(with_first_case(&b)).contains("exactly one"));
    }

    #[test]
    fn rejects_nonpositive_and_nonfinite_values() {
        for v in ["0.0", "-1.0", "nan", "inf"] {
            let b = format!("{HEAD}sigma_el_a02 = {v}\nsigma_el_cite = \"x\"\nsigma_tr1_a02 = 1.0\nsigma_tr1_cite = \"x\"");
            assert!(err(with_first_case(&b)).contains("invalid value"), "{v}");
        }
    }

    #[test]
    fn rejects_duplicate_and_missing_cases() {
        let dup = format!("{}\n[[case]]\nz = 6\nenergy_ev = 1000.0\n", good());
        assert!(err(dup).contains("expected exactly 1"));
        let t = good();
        let start = t.find("[[case]]").unwrap();
        let end = start + t[start..].find("\n\n").unwrap();
        let missing = format!("{}{}", &t[..start], &t[end..]);
        assert!(err(missing).contains("expected exactly 1"));
    }

    #[test]
    fn rejects_unsupported_correction_settings_and_loosened_tolerance() {
        let b = format!("{HEAD}correlation_polarization = true");
        assert!(err(with_first_case(&b)).contains("parse"));
        let loose = good().replace("rel_at_10000_ev = 0.02", "rel_at_10000_ev = 0.05");
        assert!(err(loose).contains("tolerances"));
    }
}

/// Compares every populated case and returns the number of values checked.
/// Needs the DHFS coefficients of `z`; `SalvatDhfs::for_element` fails until
/// #130 lands, in which case a case that has a reference value fails loudly
/// rather than being skipped.
fn run_comparisons(f: &Fixture) -> usize {
    let mut checked = 0;
    for c in &f.case {
        if c.sigma_el_a02.is_none() && c.sigma_tr1_a02.is_none() {
            eprintln!(
                "SKIPPED Z={} E={} eV: no published reference value committed (placeholder)",
                c.z, c.energy_ev
            );
            continue;
        }
        let tol = if c.energy_ev >= 10000.0 {
            TOL_10KEV
        } else {
            TOL_1KEV
        };
        let pot = SalvatDhfs::for_element(c.z)
            .unwrap_or_else(|e| panic!("Z={}: DHFS potential unavailable: {e}", c.z));
        let corr = Corrections {
            exchange: c.exchange,
            correlation_polarization: None,
        };
        let r = solve_corrected(
            &pot,
            &pot,
            c.energy_ev,
            &corr,
            &[],
            SolverOptions::default(),
        )
        .unwrap();
        for (name, ours, theirs, cite) in [
            ("sigma_el", r.sigma_el, c.sigma_el_a02, &c.sigma_el_cite),
            ("sigma_tr1", r.sigma_tr1, c.sigma_tr1_a02, &c.sigma_tr1_cite),
        ] {
            if let Some(t) = theirs {
                let gap = ours / t - 1.0;
                eprintln!(
                    "Z={} E={} {name}: ours {ours:.6e} ref {t:.6e} gap {gap:+.4} ({})",
                    c.z,
                    c.energy_ev,
                    cite.as_deref().unwrap_or("")
                );
                assert!(
                    gap.abs() <= tol,
                    "Z={} E={} {name}: gap {gap:+.4} > {tol}",
                    c.z,
                    c.energy_ev
                );
                checked += 1;
            }
        }
    }
    checked
}

fn populated(f: &Fixture) -> usize {
    f.case
        .iter()
        .map(|c| c.sigma_el_a02.is_some() as usize + c.sigma_tr1_a02.is_some() as usize)
        .sum()
}

/// Every populated value is compared (no silent skip), and the report states
/// whether coverage is complete. Zero comparisons is reported as INCOMPLETE,
/// never as a validation.
#[test]
fn elastic_cross_sections_match_published_partial_waves() {
    let f = fixture();
    let expected = populated(&f);
    let checked = run_comparisons(&f);
    assert_eq!(checked, expected, "a populated value was not compared");
    if checked == EXPECTED_VALUES {
        eprintln!("COMPLETE: {checked} of {EXPECTED_VALUES} reference values checked");
    } else {
        eprintln!(
            "INCOMPLETE: {checked} of {EXPECTED_VALUES} reference values checked; \
             this is NOT a completed validation (issue #93 stays open)"
        );
    }
}

/// Fails until all 16 values (8 cases x 2 observables) are populated and
/// within tolerance. Ignored by default so CI reflects the harness, not the
/// missing data; run with `--ignored` to gate the issue's completion.
#[test]
#[ignore = "fails until all 16 published reference values are committed (#93, #130)"]
fn complete_validation_requires_all_sixteen_values() {
    let f = fixture();
    assert_eq!(
        populated(&f),
        EXPECTED_VALUES,
        "fixture coverage incomplete"
    );
    assert_eq!(run_comparisons(&f), EXPECTED_VALUES);
}
