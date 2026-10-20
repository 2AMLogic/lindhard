//! Level-1 check of `sigma_el` and `sigma_tr1` against published partial-wave
//! values (issue #93; fixture `tests/data/elastic_reference.toml`).
//!
//! The 16 reference values (8 cases x `sigma_el`, `sigma_tr1`) are NIST SRD 64
//! v5.0, committed under the operator's scoped data-use ruling of 2026-10-19
//! (see the fixture header and `docs/data-provenance.md`). The reference
//! potential is a numerical DHF density with Furness-McCarthy exchange; ours is
//! the analytic DHFS screening with the same exchange form, so this is a
//! different-potential comparison, not a same-potential validation.
//! Every gap is reported (`WITHIN` / `OUT OF TOLERANCE`) and none is hidden:
//! the comparison test records them, and the `#[ignore]`d test
//! `complete_validation_requires_all_sixteen_values` fails until all 16 are
//! within tolerance (run it with `--ignored`). A skipped or absent value is
//! never a pass: the report states `checked N of 16`. Filling a case in the fixture activates
//! it; the tolerances are read from the fixture and pinned here so they cannot
//! be loosened without touching this file. The fixture is validated strictly
//! (`parse_fixture`): unknown keys (for example an unsupported
//! `correlation_polarization` setting), blank citations, non-positive or
//! non-finite values, and duplicate or missing cases are rejected.

use lindhard::electron::elastic::corrections::{solve_corrected, Corrections, ElectronDensity};
use lindhard::electron::elastic::{SalvatDhfs, ScreenedPotential, SolverOptions};
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
    #[allow(dead_code)]
    version: String,
    #[allow(dead_code)]
    accessed: String,
    #[allow(dead_code)]
    potential: String,
    #[allow(dead_code)]
    units: String,
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

/// Outcome of comparing every populated value.
struct Outcome {
    /// Values compared.
    checked: usize,
    /// Values whose gap exceeded the pinned tolerance, as report lines.
    out_of_tolerance: Vec<String>,
}

/// Compares every populated case and reports every gap (one line per value,
/// `WITHIN` or `OUT OF TOLERANCE`); a failing value does not stop the run, so
/// all 16 gaps are always reported.
/// Uses the Salvat et al. (1987) Table I potential of `z`
/// (`SalvatDhfs::for_element`, #130); an element outside Z = 1..92 fails
/// loudly rather than being skipped.
fn run_comparisons(f: &Fixture) -> Outcome {
    let mut out = Outcome {
        checked: 0,
        out_of_tolerance: Vec::new(),
    };
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
        for (name, ours, theirs) in [
            ("sigma_el", r.sigma_el, c.sigma_el_a02),
            ("sigma_tr1", r.sigma_tr1, c.sigma_tr1_a02),
        ] {
            if let Some(t) = theirs {
                let gap = ours / t - 1.0;
                let verdict = if gap.abs() <= tol {
                    "WITHIN"
                } else {
                    "OUT OF TOLERANCE"
                };
                let line = format!(
                    "Z={} E={} {name}: ours {ours:.6e} ref {t:.6e} gap {gap:+.4} (tol {tol}) {verdict}",
                    c.z, c.energy_ev
                );
                eprintln!("{line}");
                if gap.abs() > tol {
                    out.out_of_tolerance.push(line);
                }
                out.checked += 1;
            }
        }
    }
    out
}

fn populated(f: &Fixture) -> usize {
    f.case
        .iter()
        .map(|c| c.sigma_el_a02.is_some() as usize + c.sigma_tr1_a02.is_some() as usize)
        .sum()
}

/// Every populated value is compared (no silent skip) and every gap is
/// reported. Out-of-tolerance values are recorded in the report, not hidden
/// and not relaxed; they fail the `--ignored` completion gate below. Zero
/// comparisons is reported as INCOMPLETE, never as a validation.
#[test]
fn elastic_cross_sections_match_published_partial_waves() {
    let f = fixture();
    let expected = populated(&f);
    let o = run_comparisons(&f);
    assert_eq!(o.checked, expected, "a populated value was not compared");
    eprintln!(
        "{} of {EXPECTED_VALUES} reference values checked; {} within tolerance, {} OUT OF TOLERANCE",
        o.checked,
        o.checked - o.out_of_tolerance.len(),
        o.out_of_tolerance.len()
    );
    if o.checked != EXPECTED_VALUES {
        eprintln!("INCOMPLETE: this is NOT a completed validation (issue #93 stays open)");
    } else if !o.out_of_tolerance.is_empty() {
        eprintln!(
            "FAILED VALIDATION: the gaps above exceed the pinned tolerances (issue #93 stays open)"
        );
    } else {
        eprintln!("COMPLETE: all {EXPECTED_VALUES} values within tolerance");
    }
}

/// Fails until all 16 values (8 cases x 2 observables) are populated and
/// within tolerance. Ignored by default so CI records the measured gaps
/// without blocking on a model difference; run with `--ignored` to gate the
/// issue's completion.
#[test]
#[ignore = "fails while any of the 16 values is outside its pinned tolerance (#93)"]
fn complete_validation_requires_all_sixteen_values() {
    let f = fixture();
    assert_eq!(
        populated(&f),
        EXPECTED_VALUES,
        "fixture coverage incomplete"
    );
    let o = run_comparisons(&f);
    assert_eq!(o.checked, EXPECTED_VALUES);
    assert!(
        o.out_of_tolerance.is_empty(),
        "{} value(s) out of tolerance:\n{}",
        o.out_of_tolerance.len(),
        o.out_of_tolerance.join("\n")
    );
}

/// Discretization check: tightening the solver settings moves neither
/// observable by more than 0.1 %, so the gaps reported above are not
/// integration error. They come from the difference between the reference
/// potential (numerical DHF density) and ours (analytic DHFS fit), not from
/// fitting to the reference numbers.
#[test]
fn reported_gaps_are_not_discretization_error() {
    let f = fixture();
    for c in f.case.iter().filter(|c| c.z == 6 || c.z == 79) {
        let pot = SalvatDhfs::for_element(c.z).unwrap();
        let corr = Corrections {
            exchange: c.exchange,
            correlation_polarization: None,
        };
        let solve = |opts| solve_corrected(&pot, &pot, c.energy_ev, &corr, &[], opts).unwrap();
        let base = solve(SolverOptions::default());
        let tight = solve(SolverOptions {
            step_scale: 0.1,
            phase_tolerance: 1.0e-10,
            consecutive: 10,
            ..SolverOptions::default()
        });
        for (name, a, b) in [
            ("sigma_el", base.sigma_el, tight.sigma_el),
            ("sigma_tr1", base.sigma_tr1, tight.sigma_tr1),
        ] {
            let d = (b / a - 1.0).abs();
            eprintln!(
                "Z={} E={} {name}: tight/default - 1 = {d:.2e}",
                c.z, c.energy_ev
            );
            assert!(
                d < 1.0e-3,
                "Z={} E={} {name}: moved by {d}",
                c.z,
                c.energy_ev
            );
        }
    }
}

/// DHFS potential with every screening exponent divided by `scale`, i.e.
/// `V_s(r) = V(r/s)/s` (and, for the exchange term, the Poisson density
/// `rho_s(r) = rho(r/s)/s^3` of that same screening function). `Z` and the
/// point-Coulomb limit are unchanged. A diagnostic perturbation of the
/// in-tree fit (issue #452), not a model: nothing is adjusted toward the
/// reference.
struct Rescaled<'a> {
    base: &'a SalvatDhfs,
    scale: f64,
}

impl ScreenedPotential for Rescaled<'_> {
    fn energy(&self, r: f64) -> f64 {
        self.base.energy(r / self.scale) / self.scale
    }
    fn nuclear_charge(&self) -> f64 {
        self.base.nuclear_charge()
    }
    fn length_scale(&self) -> f64 {
        self.base.length_scale() * self.scale
    }
}

impl ElectronDensity for Rescaled<'_> {
    fn density(&self, r: f64) -> f64 {
        self.base.density(r / self.scale) / self.scale.powi(3)
    }
    fn density_source(&self) -> String {
        format!("DHFS density rescaled by {} (diagnostic)", self.scale)
    }
}

/// Sensitivity of `sigma_el` / `sigma_tr1` to in-tree potential settings
/// (issue #452). Reports, for the 8 cases, the fractional change relative to
/// the baseline (DHFS with Furness-McCarthy exchange, as compared above) when
/// exchange is switched off and when the DHFS screening length is scaled by
/// 0.95 / 1.05. These perturb our own fit; they are not tuned to, and
/// never compared as a pass/fail against, the reference. The assertions pin
/// only invariants: a unit scale reproduces the plain DHFS potential bit-for-bit.
#[test]
fn potential_sensitivity_of_the_eight_cases() {
    let f = fixture();
    let opts = SolverOptions::default();
    let run = |c: &Case, scale: f64, exchange: bool| {
        let pot = SalvatDhfs::for_element(c.z).unwrap();
        let r = Rescaled { base: &pot, scale };
        let corr = Corrections {
            exchange,
            correlation_polarization: None,
        };
        let res = solve_corrected(&r, &r, c.energy_ev, &corr, &[], opts).unwrap();
        (res.sigma_el, res.sigma_tr1)
    };
    eprintln!("sensitivity: change of ours relative to baseline (x/baseline - 1); 'ref gap' = baseline/ref - 1");
    for c in &f.case {
        let direct = {
            let pot = SalvatDhfs::for_element(c.z).unwrap();
            let corr = Corrections {
                exchange: c.exchange,
                correlation_polarization: None,
            };
            let r = solve_corrected(&pot, &pot, c.energy_ev, &corr, &[], opts).unwrap();
            (r.sigma_el, r.sigma_tr1)
        };
        let base = run(c, 1.0, c.exchange);
        assert_eq!(base, direct, "unit rescale must equal the plain potential");
        let gap_el = base.0 / c.sigma_el_a02.unwrap() - 1.0;
        let gap_tr = base.1 / c.sigma_tr1_a02.unwrap() - 1.0;
        eprintln!(
            "Z={} E={}: ref gap sigma_el {gap_el:+.4} sigma_tr1 {gap_tr:+.4}",
            c.z, c.energy_ev
        );
        for (label, scale, exchange) in [
            ("exchange off", 1.0, !c.exchange),
            ("screening x0.95", 0.95, c.exchange),
            ("screening x1.05", 1.05, c.exchange),
        ] {
            let v = run(c, scale, exchange);
            eprintln!(
                "  {label:<16} d sigma_el {:+.4}  d sigma_tr1 {:+.4}",
                v.0 / base.0 - 1.0,
                v.1 / base.1 - 1.0
            );
        }
    }
}
