//! Sum-rule checks on the committed optical ELF datasets
//! (`validation/data/optical/`, issues #98, #148 and #125).
//!
//! For `ELF(E) = Im[-1/eps(E)]` of a material with electron density `n`
//! (`n` counts all `Z` electrons per atom, `Omega_p^2 = n e^2 / (eps0 m)`):
//!
//! - f-sum rule (Bethe sum rule; Shiles et al., Phys. Rev. B 22, 1612
//!   (1980); Hagemann, Gudat and Kunz, DESY report SR-74/7 (1974), Eq. 8):
//!   `N_eff = Z * (2 / (pi * (hbar Omega_p)^2)) * Int E ELF dE`
//!   is the number of electrons per atom taking part up to the table's top.
//! - Perfect-screening sum rule (Shiles et al.): for a metal
//!   `P_eff = (2/pi) * Int ELF / E dE = 1 - 1/eps1(0) = 1`. For a
//!   nonconductor `Re[-1/eps(0)] = 1/eps1(0)` is not zero, so the target is
//!   `1 - 1/eps1(0)` (Yang et al., Phys. Rev. B 100, 245209 (2019),
//!   doi:10.1103/PhysRevB.100.245209, Eqs. (28) and (30), p. 6, with
//!   `eps1(0) = n(0)^2`; for Si they take `n(0) = 3.4155` from Palik's
//!   handbook, their Ref. 63, p. 6), which is what Si is checked against.
//!
//! Reference electron densities come from Z/A and density of NIST X-ray mass
//! attenuation Table 1 (Hubbell and Seltzer, NISTIR 5632),
//! <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html>, read
//! 2026-10-07: Al Z/A 0.48181, 2.699 g/cm3; Cu Z/A 0.45636, 8.960 g/cm3;
//! read 2026-10-08: Au Z/A 0.40108, 19.32 g/cm3; C Z/A 0.49954 (the density
//! used for C is the report's 1.5 g/cm3 of glassy carbon, not NIST's graphite);
//! read 2026-10-09: Si Z/A 0.49848, 2.330 g/cm3.
//!
//! # Quadrature
//!
//! The integrals are taken over the tabulated knots only (nothing is
//! extrapolated). The check gated first is the one for the representation a
//! consumer gets: `OpticalElf::elf()` interpolates the ELF **linearly**
//! between knots (`lerp_on_grid`), so `Int E ELF dE` and `Int ELF / E dE` are
//! integrated exactly for a piecewise-linear ELF (`linear_elf`). A
//! **power-law (log-log)** segment rule is computed as a second, separately
//! named gate, a comparison with the source (not what the library serves).
//! The trapezoid rule is printed, not gated. Measured on these tables:
//!
//! - The published knots are spaced roughly logarithmically over eight
//!   decades (1 meV to 120 keV, several knots per decade), and the integrands
//!   follow power laws between them over most of that range (Drude region,
//!   x-ray tails between edges). A chord across a convex power-law segment
//!   overestimates its integral; on a log-spaced grid the error does not
//!   shrink at high energy. So the linear representation overshoots: Al
//!   P_eff 1.094 and N_eff 1.138, Cu N_eff 1.073; the power-law rule gives
//!   0.999 and 1.075 for Al, 0.999 and 0.958 for Cu.
//! - Against the authors' own integration: the report's tables print a
//!   cumulative `N-EFF` from `eps2` (its Eq. 6) on the same rows. Integrating
//!   `eps2 = 2nk` from the transcribed rows, log-log comes closer to that
//!   column than the trapezoid in every stretch of the Al table (below 16 eV:
//!   2.70 and 2.88 against 2.564; 72.9 eV to 1.55 keV: +9.25 and +9.34
//!   against +9.056; 1.55 to 120 keV: +1.66 and +1.85 against +1.601), and in
//!   total for Cu (27.79 and 29.78 against 27.505 at 50 keV). So the
//!   power-law numbers are the better estimate of the underlying data, and
//!   the linear numbers the right description of what `elf()` returns.
//!
//! - Si (Yang et al. 2019, digitized from their Fig. 7) is knotted every
//!   0.25 eV to 25 eV and every 1 eV above (0.25 eV across the L2,3 edge),
//!   so the three rules agree to 0.1 %: P_eff = 0.9135 (linear) against the
//!   nonconductor target 1 - 1/n(0)^2 = 0.9143, where the authors' own
//!   ps-sum, with their ELF extended above 200 eV, is 1.0011 including the
//!   1/n(0)^2 term (their Table II), i.e. 0.9154 for the integral. N_eff
//!   falls short because the table stops at 199 eV (`SI_N_EFF_SHORT`).
//!
//! Over 16 to 72.9 eV, where the knots are dense, the rules agree with the
//! report's `N-EFF` increment (0.273 and 0.275 against 0.274), which also
//! confirms that the NIST densities above match the authors' normalisation
//! to better than 1 %.
//!
//! The tolerance is 5 %. A check that misses it is listed in
//! `known_failures` with the measured value pinned and the reason; the test
//! then asserts that the failure is still exactly that one, so it can neither
//! be hidden nor drift unnoticed. Run with `--nocapture` to print every
//! number, including the trapezoid values and the authors' Table 10 values.

use lindhard::constants::{AVOGADRO, ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY};
use lindhard::electron::data::OpticalElf;
use std::path::PathBuf;

const TOL: f64 = 0.05;

/// Linear `ELF` between knots, exactly integrated: what `OpticalElf::elf()` serves.
const LIBRARY: &str = "library (linear)";
/// Power-law `ELF` between knots: source-side comparison.
const POWER_LAW: &str = "power-law";

struct Case {
    file: &'static str,
    z: f64,
    z_over_a: f64,
    density_g_cm3: f64,
    /// Target of `P_eff`: `1 - 1/eps1(0)`, which is 1 for a metal.
    p_eff_target: f64,
    /// `n_eff(Im eps^-1)` of the authors, extrapolated to infinite energy:
    /// DESY report SR-74/7, Table 10 (errata sheet, PDF p. 3) for the
    /// Hagemann sets; Yang et al. (2019), Table II, for Si (their ELF to
    /// 200 eV extended with other data to 10 MeV, p. 13). For comparison
    /// only; printed, not asserted.
    authors_n_eff: f64,
    /// Strings the dataset's provenance must contain (its source).
    cites: &'static [&'static str],
    /// The table must reach at least this energy (eV).
    min_top_ev: f64,
    /// Which checks are expected to fail, with the reason: `(rule, check,
    /// pinned value, half-width of the pinned band, explanation)`.
    known_failures: &'static [(&'static str, &'static str, f64, f64, &'static str)],
}

const AL_N_EFF_EXCESS: &str =
    "measured with power-law segments: N_eff = 13.98 for Z = 13 (+7.5 %) on the 148 \
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

const AL_LINEAR_EXCESS: &str = "measured: with the linear interpolation `OpticalElf::elf()` \
     serves, integrated exactly over the 148 published knots, Al gives P_eff = 1.094 and N_eff = \
     14.79 for Z = 13 (+13.8 %). The same knots with power-law segments give P_eff = 0.999 and \
     N_eff = 13.98 (AL_N_EFF_EXCESS), so the P_eff excess and about half of the N_eff excess \
     (0.82 of 1.79 electrons) belong to the representation, not the data. By energy, the linear \
     minus power-law N_eff difference is +0.30 electrons from 10 to 30 eV (chords across the \
     sharp 15 eV plasmon peak) and +0.50 from 100 eV to 10 keV (the sparse x-ray knots); \
     nothing elsewhere. The remaining +7.5 % is the data and the printed knots, as explained \
     there. A consumer that integrates this table with `elf()` on a coarse grid inherits the \
     chord overshoot; the fix is a denser table (follow-up), not a different integrator";

const CU_LINEAR_EXCESS: &str = "measured: with the linear interpolation `OpticalElf::elf()` \
     serves, integrated exactly over the 149 published knots, Cu gives N_eff = 31.11 for Z = 29 \
     (+7.3 %), against 27.80 (-4.1 %) with power-law segments and 27.6 in the authors' Table 10. \
     The whole 3.31 electron difference lies above 100 eV (+0.87 from 100 eV to 1 keV, +2.07 \
     from 1 to 10 keV, +0.37 above), where the knots are a few per decade across the L and K \
     edges (about 0.93 and 9 keV) and the integrand falls as a power law between them, so \
     chords overshoot. P_eff is unaffected (1.002), since the screening integral is dominated \
     by the Drude region, which the table samples densely. Not a data error";

const C_N_EFF_EXCESS: &str =
    "measured on glassy carbon (density 1.5 g/cm3, the report's own basis, \
     p. 20): N_eff = 6.76 for Z = 6 (+12.7 %) with power-law segments and 7.01 (+16.9 %) with \
     the linear interpolation `elf()` serves, over the 100 published knots to 30 keV; P_eff \
     passes (0.987 and 0.991). The authors' own N_eff from the same function is 6.2 (+3.3 %; \
     DESY report SR-74/7, Table 10). Integrating eps2 = 2nk from the printed rows with \
     power-law segments and comparing with the report's own N-EFF column (Table 8): 2.384 \
     against 2.317 at 28 eV, 3.790 against 3.653 at 100 eV, 4.196 against 4.017 at 290 eV, \
     6.514 against 5.963 at 1 keV and 6.765 against 6.139 at 30 keV. So about 0.37 of the 0.63 \
     electron excess lies between the K edge (290 eV) and 1 keV, where the printed table has 12 \
     knots, a subsample of the authors' working grid; a further 3 to 4 % appears already \
     below 100 eV, which is not explained here. Not a transcription error: 17 n, k readings \
     checked against the scanned table (secondread_hagemann1975.json), 6 of them from 290 to \
     1000 eV";

const AU_LINEAR_EXCESS: &str = "measured: with the linear interpolation `OpticalElf::elf()` \
     serves, integrated exactly over the 149 published knots, Au gives N_eff = 95.58 for Z = 79 \
     (+21.0 %), against 82.60 (+4.6 %, passes) with power-law segments and 79.0 in the authors' \
     Table 10. Of the 12.98 electron difference, 4.80 lie from 100 eV to 1 keV, 6.66 from 1 to \
     10 keV and 1.51 above, where the knots are a few per decade across the N, M and L edges \
     and chords overshoot. P_eff passes in both (1.002, 1.000). Against the report's own N-EFF \
     column (Table 5), eps2 = 2nk integrated with power-law segments from the printed rows gives \
     7.47 against 7.37 at 22.5 eV, 19.50 against 19.34 at 84 eV, 47.74 against 46.53 at 1 keV \
     and 82.49 against 78.47 at 150 keV, so the density normalisation agrees within about 1 % \
     and the rest is quadrature on the printed knots. Not a data error";

const HAGEMANN: &[&str] = &["10.1364/JOSA.65.000742", "SR-74/7"];

/// `n(0)` of Si used by Yang et al. (2019) for the ps-sum rule (p. 6, from
/// Palik's handbook, their Ref. 63); `eps1(0) = n(0)^2`.
const SI_N0: f64 = 3.4155;

const SI_N_EFF_SHORT: &str = "measured: the table ends at 199 eV (Yang et al. 2019, Fig. 7, \
     digitized), so the Si K shell (binding about 1.84 keV) and the L-shell tail above 199 eV are \
     absent and N_eff counts only the electrons excited below 199 eV. With the linear \
     interpolation `elf()` serves, N_eff = 3.13 to 30 eV, 3.50 to 99 eV (the valence electrons, 4 \
     of the 14, are not all reached below the L2,3 edge) and 7.66 to 199 eV, ratio 0.547 \
     (power-law segments give the same to 0.1 %, the knots being dense). The authors report \
     N_eff = 14.158 (Table II), with their ELF extended above 200 eV by Henke's data to 30 keV \
     and atomic scattering factors to 10 MeV (p. 13). A truncation, not a data error: the \
     perfect-screening sum, which the region above 199 eV hardly feeds, passes";

const CASES: [Case; 5] = [
    Case {
        file: "al_elf_hagemann1975.toml",
        z: 13.0,
        z_over_a: 0.48181,
        density_g_cm3: 2.699,
        p_eff_target: 1.0,
        authors_n_eff: 13.6,
        cites: HAGEMANN,
        min_top_ev: 1.0e4,
        known_failures: &[
            (LIBRARY, "N_eff/Z", 1.138, 0.003, AL_LINEAR_EXCESS),
            (LIBRARY, "P_eff", 1.094, 0.003, AL_LINEAR_EXCESS),
            (POWER_LAW, "N_eff/Z", 1.075, 0.003, AL_N_EFF_EXCESS),
        ],
    },
    Case {
        file: "cu_elf_hagemann1975.toml",
        z: 29.0,
        z_over_a: 0.45636,
        density_g_cm3: 8.960,
        p_eff_target: 1.0,
        authors_n_eff: 27.6,
        cites: HAGEMANN,
        min_top_ev: 1.0e4,
        known_failures: &[(LIBRARY, "N_eff/Z", 1.073, 0.003, CU_LINEAR_EXCESS)],
    },
    Case {
        file: "c_elf_hagemann1975.toml",
        z: 6.0,
        z_over_a: 0.49954,
        // Glassy carbon: the report gives its values "on the basis of
        // rho = 1.5 g cm^-3" (DESY report SR-74/7, p. 20), not the 1.700 of
        // NIST's graphite row; Z/A is NIST's.
        density_g_cm3: 1.5,
        p_eff_target: 1.0,
        authors_n_eff: 6.2,
        cites: HAGEMANN,
        min_top_ev: 1.0e4,
        known_failures: &[
            (LIBRARY, "N_eff/Z", 1.169, 0.003, C_N_EFF_EXCESS),
            (POWER_LAW, "N_eff/Z", 1.127, 0.003, C_N_EFF_EXCESS),
        ],
    },
    Case {
        file: "au_elf_hagemann1975.toml",
        z: 79.0,
        z_over_a: 0.40108,
        density_g_cm3: 19.32,
        p_eff_target: 1.0,
        authors_n_eff: 79.0,
        cites: HAGEMANN,
        min_top_ev: 1.0e4,
        known_failures: &[(LIBRARY, "N_eff/Z", 1.210, 0.003, AU_LINEAR_EXCESS)],
    },
    Case {
        file: "si_elf_yang2019.toml",
        z: 14.0,
        z_over_a: 0.49848,
        density_g_cm3: 2.330,
        p_eff_target: 1.0 - 1.0 / (SI_N0 * SI_N0),
        authors_n_eff: 14.158,
        cites: &["10.1103/PhysRevB.100.245209", "Fig. 7", "DIGITIZED"],
        min_top_ev: 199.0,
        known_failures: &[
            (LIBRARY, "N_eff/Z", 0.547, 0.003, SI_N_EFF_SHORT),
            (POWER_LAW, "N_eff/Z", 0.547, 0.003, SI_N_EFF_SHORT),
        ],
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

/// Exact `(Int E y dE, Int y / E dE)` over `[e0, e1]` for `y` linear in `E`
/// through `(e0, y0)` and `(e1, y1)`: the representation
/// `OpticalElf::elf()` serves (`lerp_on_grid`). With `h = e1 - e0` and
/// `b = (y0 e1 - y1 e0) / h` (the intercept):
/// `Int E y dE = h (y0 (2 e0 + e1) + y1 (e0 + 2 e1)) / 6` and
/// `Int y / E dE = (y1 - y0) + b ln(e1 / e0)`.
fn linear_elf(e0: f64, e1: f64, y0: f64, y1: f64) -> (f64, f64) {
    let h = e1 - e0;
    let b = (y0 * e1 - y1 * e0) / h;
    (
        h * (y0 * (2.0 * e0 + e1) + y1 * (e0 + 2.0 * e1)) / 6.0,
        (y1 - y0) + b * (h / e0).ln_1p(),
    )
}

/// The same pair with a power-law `ELF` between knots (source-side
/// comparison, not what the library serves).
fn power_law_elf(e0: f64, e1: f64, y0: f64, y1: f64) -> (f64, f64) {
    (
        power_law(e0, e1, e0 * y0, e1 * y1),
        power_law(e0, e1, y0 / e0, y1 / e1),
    )
}

/// The same pair with the trapezoid rule on each integrand.
fn trapezoid_elf(e0: f64, e1: f64, y0: f64, y1: f64) -> (f64, f64) {
    (
        trapezoid(e0, e1, e0 * y0, e1 * y1),
        trapezoid(e0, e1, y0 / e0, y1 / e1),
    )
}

type Rule = fn(f64, f64, f64, f64) -> (f64, f64);

/// `(Int E ELF dE, Int ELF / E dE)` over the knots with the given segment rule.
fn integrals(e: &[f64], y: &[f64], seg: Rule) -> (f64, f64) {
    let (mut f_sum, mut p_sum) = (0.0, 0.0);
    for i in 0..e.len() - 1 {
        let (f, p) = seg(e[i], e[i + 1], y[i], y[i + 1]);
        f_sum += f;
        p_sum += p;
    }
    (f_sum, p_sum)
}

#[test]
fn linear_segment_is_exact_for_linear_elf() {
    for (a, b) in [(0.0, 3.0), (2.0, 0.0), (-0.5, 4.0), (1.5, -2.0)] {
        let (e0, e1) = (2.0f64, 5.0f64);
        let (f, p) = linear_elf(e0, e1, a * e0 + b, a * e1 + b);
        let f_exact = a * (e1.powi(3) - e0.powi(3)) / 3.0 + b * (e1 * e1 - e0 * e0) / 2.0;
        let p_exact = a * (e1 - e0) + b * (e1 / e0).ln();
        assert!(
            (f / f_exact - 1.0).abs() < 1e-12,
            "{a} {b}: {f} vs {f_exact}"
        );
        assert!(
            (p / p_exact - 1.0).abs() < 1e-12,
            "{a} {b}: {p} vs {p_exact}"
        );
    }
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
        let rules: [(&str, Rule, bool); 3] = [
            (LIBRARY, linear_elf, true),
            (POWER_LAW, power_law_elf, true),
            ("trapezoid", trapezoid_elf, false),
        ];
        for (rule, seg, gated) in rules {
            let (f, p) = integrals(e, y, seg);
            let (n_eff_over_z, p_eff) = (f_norm * f, p_norm * p);
            eprintln!(
                "{} [{rule}{}]: top {:.4e} eV; N_eff = {:.3} of Z = {} (ratio {:.4}), \
                 P_eff = {:.4} (target {:.4}); authors' N_eff = {}",
                t.material(),
                if gated { "" } else { ", not gated" },
                e[e.len() - 1],
                n_eff_over_z * case.z,
                case.z,
                n_eff_over_z,
                p_eff,
                case.p_eff_target,
                case.authors_n_eff
            );
            if !gated {
                continue;
            }
            for (name, value, target) in [
                ("N_eff/Z", n_eff_over_z, 1.0),
                ("P_eff", p_eff, case.p_eff_target),
            ] {
                let within = (value / target - 1.0).abs() <= TOL;
                let known = case
                    .known_failures
                    .iter()
                    .find(|k| k.0 == rule && k.1 == name);
                match known {
                    None => assert!(
                        within,
                        "{} [{rule}] {name} = {value:.4}, outside 5 % of {target:.4}",
                        t.material()
                    ),
                    Some((_, _, pinned, band, why)) => {
                        eprintln!("  KNOWN FAILURE [{rule}] {name} = {value:.4}: {why}");
                        assert!(
                            !within,
                            "{} [{rule}] {name} now passes ({value:.4}): update the table",
                            t.material()
                        );
                        assert!(
                            (value - pinned).abs() <= *band,
                            "{} [{rule}] {name} = {value:.4}, pinned {pinned}",
                            t.material()
                        );
                    }
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
        for cite in case.cites {
            assert!(t.provenance().contains(cite), "{}: {cite}", case.file);
        }
        assert!(
            t.energy_range_ev().1 >= case.min_top_ev,
            "{} reaches {} eV",
            t.material(),
            case.min_top_ev
        );
    }
}
