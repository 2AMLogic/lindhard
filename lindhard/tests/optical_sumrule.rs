//! Sum-rule checks on the committed optical ELF datasets
//! (`validation/data/optical/`, issue #98).
//!
//! For `ELF(E) = Im[-1/eps(E)]` of a material with electron density `n`
//! (`n` counts all `Z` electrons per atom, `Omega_p^2 = n e^2 / (eps0 m)`):
//!
//! - f-sum rule (Bethe sum rule; Shiles et al., Phys. Rev. B 22, 1612
//!   (1980); Hagemann, Gudat and Kunz, DESY report SR-74/7 (1974), Eq. 8):
//!   `N_eff = Z * (2 / (pi * (hbar Omega_p)^2)) * Int E ELF dE`
//!   is the number of electrons per atom taking part up to the table's top.
//! - Perfect-screening sum rule (Shiles et al.): for a metal
//!   `P_eff = (2/pi) * Int ELF / E dE = 1 - 1/eps1(0) = 1`.
//!
//! Reference electron densities come from Z/A and density of NIST X-ray mass
//! attenuation Table 1 (Hubbell and Seltzer, NISTIR 5632),
//! <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html>, read
//! 2026-10-07: Al Z/A 0.48181, 2.699 g/cm3; Cu Z/A 0.45636, 8.960 g/cm3.
//!
//! # Quadrature
//!
//! The integrals are taken over the tabulated knots only (nothing is
//! extrapolated) with **piecewise power-law (log-log) segments**: between
//! knots `f(E) = f_i (E / E_i)^b`, integrated exactly. The trapezoid rule is
//! computed too and printed, but not gated on. The reasons, all measured on
//! these tables:
//!
//! - The published knots are spaced roughly logarithmically over eight
//!   decades (1 meV to 120 keV, several knots per decade), and the integrands
//!   follow power laws between them over most of that range (Drude region,
//!   x-ray tails between edges). A chord across a convex power-law segment
//!   overestimates its integral; on a log-spaced grid the error does not
//!   shrink at high energy.
//! - Measured against the authors' own integration: the report's tables
//!   print a cumulative `N-EFF` from `eps2` (its Eq. 6) on the same rows.
//!   Integrating `eps2 = 2nk` from the transcribed rows, log-log comes closer
//!   to that column than the trapezoid in every stretch of the Al table
//!   (below 16 eV: 2.70 and 2.88 against 2.564; 72.9 eV to 1.55 keV: +9.25 and
//!   +9.34 against +9.056; 1.55 to 120 keV: +1.66 and +1.85 against +1.601),
//!   and in total for Cu (27.79 and 29.78 against 27.505 at 50 keV).
//! - With the trapezoid, Al's `P_eff` comes out 1.100; almost all of the
//!   excess is the two chords on either side of the sharp 15 eV plasmon peak
//!   (0.5 eV knots). With log-log it is 0.999.
//!
//! Over 16 to 72.9 eV, where the knots are dense, both rules reproduce the
//! report's `N-EFF` increment (0.273 and 0.275 against 0.274), which also
//! confirms that the NIST densities above match the authors' normalisation
//! to better than 1 %.
//!
//! The tolerance is 5 %. A material that misses it is listed in
//! `known_failures` with the measured value pinned and the reason; the test
//! then asserts that the failure is still exactly that one, so it can neither
//! be hidden nor drift unnoticed. Run with `--nocapture` to print every
//! number, including the trapezoid values and the authors' Table 10 values.

use lindhard::constants::{AVOGADRO, ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY};
use lindhard::electron::data::OpticalElf;
use std::path::PathBuf;

const TOL: f64 = 0.05;

struct Case {
    file: &'static str,
    z: f64,
    z_over_a: f64,
    density_g_cm3: f64,
    /// `n_eff(Im eps^-1)` of the authors, extrapolated to infinite energy:
    /// DESY report SR-74/7, Table 10 (errata sheet, PDF p. 3). For comparison
    /// only; printed, not asserted.
    authors_n_eff: f64,
    /// Which checks are expected to fail, with the reason: `(check, pinned
    /// value, half-width of the pinned band, explanation)`.
    known_failures: &'static [(&'static str, f64, f64, &'static str)],
}

const AL_N_EFF_EXCESS: &str = "measured: N_eff = 13.98 for Z = 13 (+7.5 %), log-log on the 148 \
     published knots, no extrapolation above 120 keV. The authors' own value from the same \
     function is 13.6 (+4.5 %; DESY report SR-74/7, Table 10), which they attribute to less \
     accurate absorption data from 200 to 500 eV (report p. 15: 'overestimation of n_eff by 0.5 \
     electrons'). The further +2.8 % between 13.98 and 13.6 is quadrature on the printed table, \
     which is a subsample of the authors' working grid: integrating eps2 = 2nk from the same \
     printed rows the same way gives 13.89 against the report's own N-EFF column of 13.495 at \
     120 keV (+2.9 %), the excess lying in the sparsely sampled stretches below 16 eV (+0.14) \
     and from 72.9 eV to 1.55 keV (+0.19). So the data carry the authors' +4.5 % and the knots \
     add about 3 %; neither is a transcription error (29 n, k readings checked against the \
     scanned tables, validation/data/optical/secondread_hagemann1975.json)";

const CASES: [Case; 2] = [
    Case {
        file: "al_elf_hagemann1975.toml",
        z: 13.0,
        z_over_a: 0.48181,
        density_g_cm3: 2.699,
        authors_n_eff: 13.6,
        known_failures: &[("N_eff/Z", 1.075, 0.003, AL_N_EFF_EXCESS)],
    },
    Case {
        file: "cu_elf_hagemann1975.toml",
        z: 29.0,
        z_over_a: 0.45636,
        density_g_cm3: 8.960,
        authors_n_eff: 27.6,
        known_failures: &[],
    },
];

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../validation/data/optical")
}

fn plasma_energy_ev(n_m3: f64) -> f64 {
    HBAR * (n_m3 * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE / (VACUUM_PERMITTIVITY * ELECTRON_MASS))
        .sqrt()
        / ELEMENTARY_CHARGE
}

fn trapezoid(x0: f64, x1: f64, f0: f64, f1: f64) -> f64 {
    0.5 * (x1 - x0) * (f0 + f1)
}

/// Exact integral over `[x0, x1]` of the power law through `(x0, f0)` and
/// `(x1, f1)`: `f0 x0 L (e^u - 1) / u` with `L = ln(x1/x0)`, `u = (b + 1) L`,
/// which stays exact as `b -> -1`. Falls back to the trapezoid if either value
/// is not positive.
fn power_law(x0: f64, x1: f64, f0: f64, f1: f64) -> f64 {
    if f0 <= 0.0 || f1 <= 0.0 {
        return trapezoid(x0, x1, f0, f1);
    }
    let l = (x1 / x0).ln();
    let u = (f1 / f0).ln() + l; // (b + 1) L
    let g = if u.abs() < 1e-12 { 1.0 } else { u.exp_m1() / u };
    f0 * x0 * l * g
}

/// `(Int E ELF dE, Int ELF / E dE)` over the knots with the given segment rule.
fn integrals(e: &[f64], y: &[f64], seg: fn(f64, f64, f64, f64) -> f64) -> (f64, f64) {
    let (mut f_sum, mut p_sum) = (0.0, 0.0);
    for i in 0..e.len() - 1 {
        f_sum += seg(e[i], e[i + 1], e[i] * y[i], e[i + 1] * y[i + 1]);
        p_sum += seg(e[i], e[i + 1], y[i] / e[i], y[i + 1] / e[i + 1]);
    }
    (f_sum, p_sum)
}

#[test]
fn power_law_segment_is_exact_for_power_laws() {
    for b in [-3.5, -1.0, -0.5, 0.0, 1.0, 2.0] {
        let f = |x: f64| 2.0 * x.powf(b);
        let exact = if b == -1.0 {
            2.0 * (5.0f64 / 2.0).ln()
        } else {
            2.0 * (5.0f64.powf(b + 1.0) - 2.0f64.powf(b + 1.0)) / (b + 1.0)
        };
        let got = power_law(2.0, 5.0, f(2.0), f(5.0));
        assert!(
            (got / exact - 1.0).abs() < 1e-12,
            "b = {b}: {got} vs {exact}"
        );
    }
}

#[test]
fn optical_elf_sum_rules() {
    if !data_dir().is_dir() {
        // A packaged crate does not ship validation/data.
        eprintln!("validation/data/optical not found; sum-rule check skipped");
        return;
    }
    for case in &CASES {
        let t = OpticalElf::from_toml_file(data_dir().join(case.file)).unwrap();
        let (e, y) = (t.energy_ev(), t.elf_values());
        let n = case.density_g_cm3 * case.z_over_a * AVOGADRO * 1e6;
        let wp = plasma_energy_ev(n);
        let f_norm = 2.0 / (std::f64::consts::PI * wp * wp);
        let p_norm = 2.0 / std::f64::consts::PI;
        let (f_ll, p_ll) = integrals(e, y, power_law);
        let (f_tr, p_tr) = integrals(e, y, trapezoid);
        let (n_eff_over_z, p_eff) = (f_norm * f_ll, p_norm * p_ll);
        eprintln!(
            "{}: top {:.4e} eV; log-log N_eff = {:.3} of Z = {} (ratio {:.4}), P_eff = {:.4}; \
             trapezoid (not gated) N_eff = {:.3}, P_eff = {:.4}; authors (Table 10) N_eff = {}",
            t.material(),
            e[e.len() - 1],
            n_eff_over_z * case.z,
            case.z,
            n_eff_over_z,
            p_eff,
            f_norm * f_tr * case.z,
            p_norm * p_tr,
            case.authors_n_eff
        );
        for (name, value) in [("N_eff/Z", n_eff_over_z), ("P_eff", p_eff)] {
            let within = (value - 1.0).abs() <= TOL;
            match case.known_failures.iter().find(|k| k.0 == name) {
                None => assert!(
                    within,
                    "{} {name} = {value:.4}, outside 5 % of 1",
                    t.material()
                ),
                Some((_, pinned, band, why)) => {
                    eprintln!("  KNOWN FAILURE {name} = {value:.4}: {why}");
                    assert!(
                        !within,
                        "{} {name} now passes ({value:.4}): update the table",
                        t.material()
                    );
                    assert!(
                        (value - pinned).abs() <= *band,
                        "{} {name} = {value:.4}, pinned {pinned}",
                        t.material()
                    );
                }
            }
        }
    }
}

#[test]
fn optical_elf_datasets_load_and_cite_their_source() {
    if !data_dir().is_dir() {
        return;
    }
    for case in &CASES {
        let t = OpticalElf::from_toml_file(data_dir().join(case.file)).unwrap();
        assert!(t.provenance().contains("10.1364/JOSA.65.000742"));
        assert!(t.provenance().contains("SR-74/7"));
        assert!(
            t.energy_range_ev().1 > 1.0e4,
            "{} covers the x-ray region",
            t.material()
        );
    }
}
