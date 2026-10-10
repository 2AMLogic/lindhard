//! Inelastic validation (issue #99): our full Penn IMFPs, built from the
//! committed optical ELFs of Al and Cu (`validation/data/optical/`, #98) and
//! of Si (#125),
//! against the published IMFPs of Tanuma, Powell & Penn, Surf. Interface
//! Anal. 43, 689 (2011), doi:10.1002/sia.3522 ("TPP 2011"), Table 4, at the
//! rows nearest 100 eV, 1 keV and 10 keV. Fixture:
//! `tests/data/imfp_tpp2011.toml`, which cites every value to its table, row
//! and manuscript page.
//!
//! What is compared, and what differs by construction:
//!
//! - Ours is [`FullPenn`] at every energy, with the Fermi energy of TPP 2011
//!   Table 1 (kinetic energies above the Fermi level on both sides). TPP 2011
//!   used the full Penn algorithm up to 300 eV and the single-pole
//!   approximation from 330 eV up (manuscript pp. 4-5), and report the two
//!   differ by < 0.2 % at 300 eV for graphite; our own full / single-pole
//!   ratio on a synthetic Drude ELF is 1.002 at 1 keV (`full_penn` module
//!   docs). So the 1 and 10 keV rows compare the same physics to well inside
//!   the tolerance.
//! - The optical inputs differ: ours is Hagemann, Gudat & Kunz (1975) as
//!   committed, served with linear interpolation between knots, which
//!   overshoots the f-sum (Al N_eff +13.8 %, Cu +7.3 %;
//!   `tests/optical_sumrule.rs`). TPP 2011 took Cu from the same Hagemann
//!   measurement for 1-95 eV but from Henke et al. (their Ref. 22) for
//!   101.94 eV to 30 keV (their Table 2, manuscript p. 33); the Al sources
//!   are in Table 2 of their Ref. 10, which was not read. They report f-sum
//!   errors of +0.9 % (Al) and -1.5 % (Cu) (their Table 3, p. 35). A larger
//!   ELF gives a shorter IMFP.
//! - Si: ours is the ELF of Yang et al., Phys. Rev. B 100, 245209 (2019),
//!   digitized from their Fig. 7 (from REELS, not an optical measurement),
//!   which ends at 199 eV: nothing above 199 eV (L-shell tail, K shell)
//!   enters our IMFP, which can only lengthen it at 1 and 10 keV. TPP 2011
//!   took Si from Palik's handbook (their Ref. 23) to 2 keV and Henke et al.
//!   above (their Table 2, manuscript p. 33), and E_F = 12.5 eV (Table 1).
//! - Below 200 eV TPP 2011 expect larger uncertainties (manuscript p. 9), and
//!   our full Penn procedure is our own reading of S2017 (`full_penn` module
//!   docs), not Penn's formulae, so the 100 eV tolerance is wider.
//!
//! The tolerances are fixed by the issue (20 % near 100 eV, 10 % near 1 and
//! 10 keV) and pinned here, so they cannot be loosened by editing the fixture
//! alone. Nothing is tuned: the measured gaps are reported in
//! `docs/validation.md` as they come out. Run with `--nocapture` for the
//! table.
//!
//! The integrals run at a relative tolerance of 1e-3 (as `tests/penn_full.rs`)
//! to keep the unoptimised test build affordable; that is a few times 1e-3 of
//! error, two orders below the tolerances checked here. The cases are
//! evaluated sequentially; each is a pure function of its inputs.

use lindhard::electron::data::OpticalElf;
use lindhard::electron::inelastic::FullPenn;
use serde::Deserialize;
use std::path::PathBuf;

const TOL_100EV: f64 = 0.20;
const TOL_1KEV: f64 = 0.10;
const TOL_10KEV: f64 = 0.10;

/// Relative tolerance of the full Penn integrals in this test.
const INTEGRATION_TOL: f64 = 1e-3;

#[derive(Deserialize)]
struct Tolerance {
    rel_at_100_ev: f64,
    rel_at_1000_ev: f64,
    rel_at_10000_ev: f64,
}

#[derive(Deserialize)]
struct Material {
    symbol: String,
    optical: String,
    fermi_energy_ev: f64,
    fermi_cite: String,
}

#[derive(Deserialize)]
struct Case {
    material: String,
    nominal_energy_ev: f64,
    energy_ev: f64,
    imfp_angstrom: f64,
    cite: String,
}

#[derive(Deserialize)]
struct Fixture {
    tolerance: Tolerance,
    material: Vec<Material>,
    case: Vec<Case>,
}

fn fixture() -> Fixture {
    toml::from_str(include_str!("data/imfp_tpp2011.toml")).expect("fixture parses")
}

fn optical_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../validation/data/optical")
}

fn tolerance_for(nominal_energy_ev: f64) -> f64 {
    if nominal_energy_ev >= 10000.0 {
        TOL_10KEV
    } else if nominal_energy_ev >= 1000.0 {
        TOL_1KEV
    } else {
        TOL_100EV
    }
}

#[test]
fn fixture_covers_al_cu_and_si_at_the_three_energies_with_the_issue_tolerances() {
    let f = fixture();
    assert_eq!(f.tolerance.rel_at_100_ev, TOL_100EV);
    assert_eq!(f.tolerance.rel_at_1000_ev, TOL_1KEV);
    assert_eq!(f.tolerance.rel_at_10000_ev, TOL_10KEV);
    assert_eq!(f.material.len(), 3);
    for m in &f.material {
        assert!(m.fermi_energy_ev > 0.0, "{}", m.symbol);
        assert!(m.fermi_cite.contains("Table 1"), "{}", m.symbol);
    }
    // Table 4 rows nearest the nominal energies (10 % logarithmic grid).
    let rows = [(100.0, 99.5), (1000.0, 992.3), (10000.0, 9897.1)];
    for sym in ["Al", "Cu", "Si"] {
        assert!(f.material.iter().any(|m| m.symbol == sym), "{sym}");
        for (nominal, row) in rows {
            let hits: Vec<&Case> = f
                .case
                .iter()
                .filter(|c| c.material == sym && c.nominal_energy_ev == nominal)
                .collect();
            assert_eq!(hits.len(), 1, "{sym} {nominal}");
            assert_eq!(hits[0].energy_ev, row, "{sym} {nominal}");
        }
    }
    assert_eq!(f.case.len(), 9);
    // Every value carries its table, row and page.
    for c in &f.case {
        assert!(c.imfp_angstrom > 0.0);
        assert!(c.cite.contains("Table 4"), "{}", c.cite);
        assert!(
            c.cite.contains(&format!("column {}", c.material)),
            "{}",
            c.cite
        );
        assert!(
            c.cite.contains(&format!("row {} eV", c.energy_ev)),
            "{}",
            c.cite
        );
        assert!(c.cite.contains(" p. "), "{}", c.cite);
    }
}

#[test]
fn full_penn_imfps_match_tpp_2011_for_al_cu_and_si() {
    if !optical_dir().is_dir() {
        // A packaged crate does not ship validation/data.
        eprintln!("validation/data/optical not found; TPP comparison skipped");
        return;
    }
    let f = fixture();
    let mut failures = Vec::new();
    eprintln!("material  E/eV     ours/A    TPP/A    rel.diff  tol");
    for m in &f.material {
        let elf = OpticalElf::from_toml_file(optical_dir().join(&m.optical)).unwrap();
        let model = FullPenn::new(elf)
            .with_fermi_energy_ev(m.fermi_energy_ev)
            .unwrap()
            .with_relative_tolerance(INTEGRATION_TOL)
            .unwrap();
        for c in f.case.iter().filter(|c| c.material == m.symbol) {
            let ours = model.imfp_m(c.energy_ev).unwrap() * 1e10;
            let gap = ours / c.imfp_angstrom - 1.0;
            let tol = tolerance_for(c.nominal_energy_ev);
            eprintln!(
                "{:<8}  {:<7}  {:>7.3}  {:>7.2}  {:>+8.4}  {:.2}  ({})",
                m.symbol, c.energy_ev, ours, c.imfp_angstrom, gap, tol, c.cite
            );
            if gap.abs() > tol {
                failures.push(format!(
                    "{} at {} eV: ours {ours:.3} A, TPP {} A, gap {gap:+.4} > {tol}",
                    m.symbol, c.energy_ev, c.imfp_angstrom
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
