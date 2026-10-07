//! Level-1 check of `sigma_el` and `sigma_tr1` against published partial-wave
//! values (issue #93; fixture `tests/data/elastic_reference.toml`).
//!
//! **No reference value is committed yet.** The intended source (Jablonski,
//! Salvat & Powell 2004) was read, but tabulates none of these cross sections;
//! its figures cover H, Al, Ni, Ag, Au and Cm (plus some gases), not C, Si or
//! Cu. The DHFS screening
//! coefficients for C, Si and Au are not in the tree either (#130), so the
//! harness below is in place but every case is skipped and reported as such.
//! A skipped case is not a pass. Filling a case in the fixture activates it;
//! the tolerances are read from the fixture and pinned here so they cannot be
//! loosened without touching this file.

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

#[derive(Deserialize)]
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
struct Fixture {
    tolerance: Tolerance,
    case: Vec<Case>,
}

fn fixture() -> Fixture {
    let text = include_str!("data/elastic_reference.toml");
    toml::from_str(text).expect("fixture parses")
}

#[test]
fn fixture_covers_the_required_cases_with_the_required_tolerances() {
    let f = fixture();
    assert_eq!(f.tolerance.rel_at_1000_ev, TOL_1KEV);
    assert_eq!(f.tolerance.rel_at_10000_ev, TOL_10KEV);
    for z in [6, 14, 29, 79] {
        for e in [1000.0, 10000.0] {
            assert_eq!(
                f.case
                    .iter()
                    .filter(|c| c.z == z && c.energy_ev == e)
                    .count(),
                1,
                "Z={z} E={e}"
            );
        }
    }
    assert_eq!(f.case.len(), 8);
    // A value without its table-and-page citation is not allowed.
    for c in &f.case {
        assert_eq!(c.sigma_el_a02.is_some(), c.sigma_el_cite.is_some());
        assert_eq!(c.sigma_tr1_a02.is_some(), c.sigma_tr1_cite.is_some());
    }
}

/// Needs the DHFS coefficients of `z`; `SalvatDhfs::for_element` fails until
/// #130 lands, in which case a case that has a reference value fails loudly
/// rather than being skipped.
#[test]
fn elastic_cross_sections_match_published_partial_waves() {
    let f = fixture();
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
    eprintln!("{checked} reference values checked (0 expected until the source is reachable)");
}
