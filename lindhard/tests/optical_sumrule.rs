//! Sum-rule checks on the committed optical ELF datasets
//! (`validation/data/optical/`, issue #98).
//!
//! For `ELF(E) = Im[-1/eps(E)]` of a material with electron density `n`
//! (`n` counts all `Z` electrons per atom, `Omega_p^2 = n e^2 / (eps0 m)`):
//!
//! - f-sum rule (Bethe sum rule; Shiles et al., Phys. Rev. B 22, 1612
//!   (1980)): `N_eff = Z * (2 / (pi * (hbar Omega_p)^2)) * Int E ELF dE`
//!   is the number of electrons per atom taking part up to the table's top.
//! - Perfect-screening sum rule (same reference): for a metal
//!   `P_eff = (2/pi) * Int ELF / E dE = 1 - 1/eps1(0) = 1`.
//!
//! Integrals are trapezoids over the tabulated knots; nothing is extrapolated.
//! Reference electron densities come from Z/A and density of NIST X-ray mass
//! attenuation Table 1 (Hubbell and Seltzer, NISTIR 5632),
//! <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html>, read
//! 2026-10-07: Al Z/A 0.48181, 2.699 g/cm3; Cu Z/A 0.45636, 8.960 g/cm3.
//!
//! The tolerance is 5 %. A material that misses it is listed in `KNOWN_FAILURES`
//! with the measured value pinned and the reason; the test then asserts that
//! the failure is still exactly that one, so it can neither be hidden nor drift
//! unnoticed. Run with `--nocapture` to print every number.

use lindhard::constants::{AVOGADRO, ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY};
use lindhard::electron::data::OpticalElf;
use std::path::PathBuf;

const TOL: f64 = 0.05;

struct Case {
    file: &'static str,
    z: f64,
    z_over_a: f64,
    density_g_cm3: f64,
    /// Which checks are expected to fail, with the reason: `(check, pinned
    /// value, half-width of the pinned band, explanation)`.
    known_failures: &'static [(&'static str, f64, f64, &'static str)],
}

const AL_EXCESS: &str = "measured, cause not established. Cumulative from the table: N_eff = \
     3.003 by 40 eV (the 3 valence electrons), 12.68 by 1.5 keV (at most 11 electrons can act \
     below the K edge), 14.54 at the top (Z = 13); P_eff = 1.060 by 40 eV, 1.100 from 1.5 keV \
     on. The tail above 40 eV lies 3 to 28 % above the Henke et al. f2-derived ELF at 500 eV to \
     20 keV (spot checks, see docs/data-provenance.md). Candidate causes, none tested: the \
     source table's x-ray part, the rounding of its n (4 digits) and k (3 digits), and linear \
     interpolation across the sharp 15 eV plasmon peak (0.5 eV knots)";

const CASES: [Case; 2] = [
    Case {
        file: "al_elf_hagemann1975.toml",
        z: 13.0,
        z_over_a: 0.48181,
        density_g_cm3: 2.699,
        known_failures: &[
            ("N_eff/Z", 1.118, 0.005, AL_EXCESS),
            ("P_eff", 1.100, 0.005, AL_EXCESS),
        ],
    },
    Case {
        file: "cu_elf_hagemann1975.toml",
        z: 29.0,
        z_over_a: 0.45636,
        density_g_cm3: 8.960,
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
        let (mut f_sum, mut p_sum) = (0.0, 0.0);
        for i in 0..e.len() - 1 {
            let h = e[i + 1] - e[i];
            f_sum += 0.5 * h * (e[i] * y[i] + e[i + 1] * y[i + 1]);
            p_sum += 0.5 * h * (y[i] / e[i] + y[i + 1] / e[i + 1]);
        }
        let n = case.density_g_cm3 * case.z_over_a * AVOGADRO * 1e6;
        let wp = plasma_energy_ev(n);
        let n_eff_over_z = 2.0 / (std::f64::consts::PI * wp * wp) * f_sum;
        let p_eff = 2.0 / std::f64::consts::PI * p_sum;
        eprintln!(
            "{}: top {:.4e} eV, N_eff = {:.3} of Z = {} (ratio {:.4}), P_eff = {:.4}",
            t.material(),
            e[e.len() - 1],
            n_eff_over_z * case.z,
            case.z,
            n_eff_over_z,
            p_eff
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
        assert!(
            t.energy_range_ev().1 > 1.0e4,
            "{} covers the x-ray region",
            t.material()
        );
    }
}
