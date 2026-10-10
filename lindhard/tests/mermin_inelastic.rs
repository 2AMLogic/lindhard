//! Acceptance tests for the Mermin-ELF model (`lindhard::electron::inelastic::
//! {MerminPenn, MerminGas, fit_mermin_oscillators}`).
//!
//! Every ELF but one is **synthetic** (Drude-Lorentz oscillators with
//! made-up parameters; not a material), kinematics are nonrelativistic. The
//! reference values are limits and sum rules of the models themselves and the
//! other two models of the selector. The exception is the Cu section at the
//! end (#300), which pins what the default fit does to the committed Cu ELF.

use lindhard::electron::data::OpticalElf;
use lindhard::electron::inelastic::{
    fit_mermin_oscillators, DrudeLorentz, DrudeLorentzOscillator, MerminFitOptions, MerminPenn,
    PennAlgorithm, PennInelastic, SinglePolePenn,
};
use std::path::PathBuf;

const TOL: f64 = 1e-4;

fn fixture() -> OpticalElf {
    DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])
        .unwrap()
        .to_optical_elf("synthetic Drude plasmon", 0.05, 2e4, 240)
        .unwrap()
}

fn two_oscillator_fixture() -> OpticalElf {
    DrudeLorentz::new(vec![
        DrudeLorentzOscillator {
            strength_ev2: 200.0,
            energy_ev: 17.0,
            width_ev: 3.0,
        },
        DrudeLorentzOscillator {
            strength_ev2: 150.0,
            energy_ev: 120.0,
            width_ev: 80.0,
        },
    ])
    .unwrap()
    .to_optical_elf("synthetic two oscillators", 0.05, 2e4, 300)
    .unwrap()
}

fn model(elf: OpticalElf, n: usize) -> MerminPenn {
    MerminPenn::fit(
        elf,
        &MerminFitOptions {
            n_oscillators: n,
            ..Default::default()
        },
    )
    .unwrap()
    .with_relative_tolerance(TOL)
    .unwrap()
}

#[test]
fn the_fit_reproduces_the_optical_elf_and_reports_the_sum_rules() {
    let elf = two_oscillator_fixture();
    let fit = fit_mermin_oscillators(
        &elf,
        &MerminFitOptions {
            n_oscillators: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(fit.rms_residual < 1e-5, "{fit}");
    assert!(fit.oscillators.iter().all(|o| o.strength_ev2 >= 0.0));
    let exact = DrudeLorentz::new(vec![
        DrudeLorentzOscillator {
            strength_ev2: 200.0,
            energy_ev: 17.0,
            width_ev: 3.0,
        },
        DrudeLorentzOscillator {
            strength_ev2: 150.0,
            energy_ev: 120.0,
            width_ev: 80.0,
        },
    ])
    .unwrap();
    assert!((fit.f_sum_ev2 / exact.f_sum_ev2() - 1.0).abs() < 1e-3);
    assert!((fit.p_eff / exact.p_eff() - 1.0).abs() < 1e-3);
    // the model's q = 0 loss function is the fitted ELF
    let m = model(elf, 2);
    for w in [3.0, 17.0, 60.0, 400.0] {
        let want = m.fit_report().model_elf(w);
        assert!((m.loss_function(0.0, w) - want).abs() <= 1e-15 + 1e-12 * want);
    }
}

#[test]
fn loss_function_decays_and_broadens_with_momentum() {
    let m = model(fixture(), 1);
    // plasmon at 20 eV: the peak height falls with q (Landau damping and
    // dispersion take over)
    let peak = |q: f64| {
        (100..=400)
            .map(|i| m.loss_function(q, i as f64 * 0.1))
            .fold(0.0_f64, f64::max)
    };
    assert!(peak(1e8) > peak(1e10));
    assert!(peak(1e10) > peak(5e10));
}

#[test]
fn high_energy_limit_agrees_with_the_single_pole_model() {
    // Both reach the Bethe limit, fixed by the f-sum and the mean excitation
    // energy of the ELF; compare the stopping power at 10 keV.
    let elf = fixture();
    let mermin = model(elf.clone(), 1);
    let spa = SinglePolePenn::new(elf);
    let e = 1e4;
    let sm = mermin.imfp_and_stopping(e).unwrap();
    let ss = spa.imfp_and_stopping(e).unwrap();
    let r = sm.stopping_ev_per_m / ss.stopping_ev_per_m;
    assert!((r - 1.0).abs() < 0.03, "stopping ratio {r}");
    let ri = sm.inverse_imfp_per_m / ss.inverse_imfp_per_m;
    assert!((ri - 1.0).abs() < 0.10, "inverse IMFP ratio {ri}");
}

#[test]
fn diimfp_integrates_to_the_inverse_imfp_and_the_stopping_power() {
    let m = model(fixture(), 1);
    let e = 500.0;
    let p = m.imfp_and_stopping(e).unwrap();
    // independent trapezoid in ln w of the DIIMFP
    let (lo, hi, n) = ((1e-6f64 * e).ln(), e.ln(), 400usize);
    let (mut inv, mut stop) = (0.0, 0.0);
    for i in 0..=n {
        let w = (lo + (hi - lo) * i as f64 / n as f64).exp();
        let wt = if i == 0 || i == n { 0.5 } else { 1.0 };
        let d = m.diimfp_per_m_ev(e, w).unwrap() * w;
        inv += wt * d;
        stop += wt * d * w;
    }
    let h = (hi - lo) / n as f64;
    assert!((inv * h / p.inverse_imfp_per_m - 1.0).abs() < 2e-2);
    assert!((stop * h / p.stopping_ev_per_m - 1.0).abs() < 2e-2);
    assert_eq!(m.diimfp_per_m_ev(e, e * 1.01).unwrap(), 0.0);
}

#[test]
fn the_model_is_selectable_beside_single_pole_and_full_and_recorded() {
    let elf = fixture();
    let mm = PennInelastic::try_new(PennAlgorithm::Mermin, elf.clone()).unwrap();
    let sp = PennInelastic::new(PennAlgorithm::SinglePole, elf.clone());
    let fp = PennInelastic::new(PennAlgorithm::Full, elf.clone());
    assert_eq!(mm.algorithm(), PennAlgorithm::Mermin);
    let (idm, ids, idf) = (
        mm.model_identity(),
        sp.model_identity(),
        fp.model_identity(),
    );
    assert!(idm.starts_with("mermin-melf ("), "{idm}");
    assert!(idm.contains("oscillators") && idm.contains("nonrelativistic"));
    assert!(ids.starts_with("penn-single-pole;") && idf.starts_with("penn-full;"));
    assert_ne!(idm, ids);
    assert_ne!(idm, idf);
    assert_eq!(
        "mermin-melf".parse::<PennAlgorithm>().unwrap(),
        PennAlgorithm::Mermin
    );
    let fermi = mm.with_fermi_energy_ev(10.0).unwrap();
    assert!(fermi.model_identity().contains("E_F = 10 eV"));
    assert_eq!(fermi.optical_elf().material(), "synthetic Drude plasmon");
    // chosen fit options flow into the identity
    let two = PennInelastic::mermin(
        elf,
        &MerminFitOptions {
            n_oscillators: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(two.model_identity().contains("2 oscillators"));
}

#[test]
fn rejects_invalid_inputs() {
    let m = model(fixture(), 1);
    assert!(m.imfp_and_stopping(0.0).is_err());
    assert!(m.imfp_and_stopping(f64::NAN).is_err());
    assert!(m.diimfp_per_m_ev(-1.0, 1.0).is_err());
    assert!(m.clone().with_relative_tolerance(0.5).is_err());
    assert!(m.clone().with_fermi_energy_ev(-1.0).is_err());
    assert_eq!(m.loss_function(1e9, -1.0), 0.0);
}

// ---- Cu: the default fit of the committed ELF (#300) ----
//
// The `mermin-melf` IMFP of Cu is 15 to 41 % longer than TPP 2011 and 18 to
// 34 % longer than the single-pole model on the same table (`docs/validation.md`, "Cu: the
// Mermin IMFP and the default oscillator fit"). The tests below pin the
// finding: the default fit (3 oscillators) carries 57 % of the table's f-sum
// and 67 % of its `P_eff`, most of the missing `P_eff` below 20 eV; the
// single-pole model on that *fitted* ELF is longer still than the Mermin
// model on it; so the excess is the fit, not the Mermin dielectric function
// or the momentum integral. de Vera et al., Int. J. Mol. Sci. 23, 6121
// (2022), doi:10.3390/ijms23116121, section 2.1.1 (opened 2026-10-09):
// "The consistency of the fitting procedure is checked by fulfilling the
// Kramers-Kronig and f-sum rules"; this fit does not pass that check.
// Model Fermi energy 0, as in the delta(E) runs. Values measured
// 2026-10-09 (release build); the allowances cover the debug build and
// nothing else. Nothing is tuned.

/// TPP 2011 Table 4 (continued), column Cu, rows 54.6, 99.5, 200.3, 492.7
/// and 992.3 eV (manuscript p. 41 of 60; S. Tanuma, C. J. Powell and D. R.
/// Penn, Surf. Interface Anal. 43, 689 (2011), doi:10.1002/sia.3522, read in
/// the authors' manuscript, NIMS MDR doi:10.48505/nims.3238), Å.
const CU_TPP2011: [(f64, f64); 5] = [
    (54.6, 4.94),
    (99.5, 5.00),
    (200.3, 6.29),
    (492.7, 10.3),
    (992.3, 16.6),
];

/// Measured outputs of this code (2026-10-09), Å, at the energies of
/// [`CU_TPP2011`]: `mermin-melf` with the default fit, single-pole on the
/// table, single-pole on the fitted ELF sampled at 4000 log-spaced points.
const CU_MERMIN: [f64; 5] = [6.967, 6.282, 7.430, 11.807, 19.213];
const CU_SINGLE_POLE: [f64; 5] = [5.186, 5.065, 6.285, 10.021, 16.164];
const CU_SINGLE_POLE_ON_FIT: [f64; 5] = [9.286, 7.533, 8.504, 12.869, 20.606];

/// Relative allowance on the pinned values above.
const PIN: f64 = 5e-3;

fn cu_elf() -> Option<OpticalElf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../validation/data/optical/cu_elf_hagemann1975.toml");
    if !p.is_file() {
        // A packaged crate does not ship validation/data.
        eprintln!("{} not found; Cu Mermin checks skipped", p.display());
        return None;
    }
    Some(OpticalElf::from_toml_file(p).unwrap())
}

fn close(got: f64, want: f64, what: &str) {
    assert!(
        (got / want - 1.0).abs() < PIN,
        "{what}: {got} vs pinned {want}"
    );
}

#[test]
fn cu_default_fit_misses_the_sum_rules_of_its_table() {
    let Some(elf) = cu_elf() else { return };
    let fit = fit_mermin_oscillators(&elf, &MerminFitOptions::default()).unwrap();
    assert_eq!(fit.oscillators.len(), 3, "{fit}");
    assert!(fit.converged, "{fit}");
    // the table: P_eff 1.002 (a metal: perfect screening), f-sum 5721 eV^2
    close(fit.data_p_eff, 1.0017, "data P_eff");
    close(fit.data_f_sum_ev2, 5720.9, "data f-sum");
    // the fit: 57 % of the f-sum, 67 % of P_eff
    close(fit.f_sum_ev2 / fit.data_f_sum_ev2, 0.5668, "f-sum ratio");
    close(fit.p_eff, 0.6703, "fit P_eff");
    close(fit.weighted_rms, 0.4029, "weighted rms");
    // most of the missing P_eff is below 20 eV: (2/pi) ∫ ELF dW/W there
    let peff_below = |f: &dyn Fn(f64) -> f64| {
        let (a, b, n) = (elf.energy_range_ev().0.ln(), 20f64.ln(), 20000);
        let du = (b - a) / n as f64;
        let s: f64 = (0..n).map(|i| f((a + (i as f64 + 0.5) * du).exp())).sum();
        2.0 / std::f64::consts::PI * s * du
    };
    let data = peff_below(&|w| elf.elf(w).unwrap());
    let model = peff_below(&|w| fit.model_elf(w));
    close(data, 0.5353, "data P_eff below 20 eV");
    close(model, 0.2433, "fit P_eff below 20 eV");
    let missing = fit.data_p_eff - fit.p_eff;
    assert!(data - model > 0.85 * missing, "{data} {model} {missing}");
}

#[test]
fn cu_mermin_imfp_excess_is_the_fit_not_the_momentum_extension() {
    let Some(elf) = cu_elf() else { return };
    let fit = fit_mermin_oscillators(&elf, &MerminFitOptions::default()).unwrap();
    let (lo, hi) = elf.energy_range_ev();
    let fitted_elf = DrudeLorentz::new(fit.oscillators.clone())
        .unwrap()
        .to_optical_elf("Cu, default Mermin fit", lo, hi, 4000)
        .unwrap();
    let mermin = MerminPenn::new(elf.clone(), fit).unwrap();
    let sp = SinglePolePenn::new(elf);
    let sp_fit = SinglePolePenn::new(fitted_elf);
    eprintln!("E/eV    TPP    Mermin   single-pole  single-pole(fit ELF)");
    for (i, &(e, tpp)) in CU_TPP2011.iter().enumerate() {
        let m = mermin.imfp_m(e).unwrap() * 1e10;
        let s = sp.imfp_m(e).unwrap() * 1e10;
        let sf = sp_fit.imfp_m(e).unwrap() * 1e10;
        eprintln!("{e:<6}  {tpp:<5}  {m:.3}  {s:.3}  {sf:.3}");
        close(m, CU_MERMIN[i], "Mermin IMFP");
        close(s, CU_SINGLE_POLE[i], "single-pole IMFP");
        close(
            sf,
            CU_SINGLE_POLE_ON_FIT[i],
            "single-pole IMFP on the fitted ELF",
        );
        // the gap of the issue: Mermin 15 to 41 % longer than TPP 2011,
        // single-pole within 6 % (+5.0 % at 54.6 eV)
        assert!(m / tpp - 1.0 > 0.14, "{e}: Mermin {m} vs TPP {tpp}");
        assert!(
            (s / tpp - 1.0).abs() < 0.06,
            "{e}: single-pole {s} vs TPP {tpp}"
        );
        // the fitted ELF alone lengthens the single-pole IMFP by more than
        // the whole Mermin excess, and the Mermin extension of that same
        // ELF shortens it: the excess is the fit
        assert!(sf > m && m > s, "{e}: {sf} > {m} > {s}");
    }
}

// ---- Cu: fits with more oscillators and the width floor (#307) ----
//
// Without a lower bound on the widths, the 10- and 16-oscillator fits of the
// committed Cu ELF placed an oscillator of width 3e-4 eV (1e-4 eV) between
// the knots at 999.9 and 1100.1 eV (1100.1 and 1200.2 eV), which the
// residuals at the knots do not see and the sum rules do: f-sum 6771 (16467)
// times the table's, `P_eff` 22.7 (46.4), Mermin IMFP 1.6 to 45 times shorter
// than with 6 oscillators. The fit now bounds every width below by the local
// knot spacing at the oscillator's energy (our choice, documented on
// `fit_mermin_oscillators`). Values below measured 2026-10-09 (release and
// debug builds agree); the bounds are stated next to each assertion.

/// Mermin IMFP of Cu with the 6-oscillator fit, Å, at the energies of
/// [`CU_TPP2011`] (measured 2026-10-09; the width floor does not bind in that
/// fit, which is bitwise the same as before #307).
const CU_MERMIN_6: [f64; 5] = [6.002, 5.731, 6.974, 11.274, 18.466];

/// The local knot spacing at `e` as `fit_mermin_oscillators` documents it:
/// the larger of the two gaps next to the knot nearest to `e` (the single
/// gap of an end knot; ties to the lower knot). Written again here, so the
/// test does not read the floor from the code under test.
fn knot_spacing(w: &[f64], e: f64) -> f64 {
    let n = w.len();
    let k = match w.iter().position(|&x| x >= e) {
        None => n - 1,
        Some(0) => 0,
        Some(i) => {
            if e - w[i - 1] <= w[i] - e {
                i - 1
            } else {
                i
            }
        }
    };
    let left = if k > 0 { w[k] - w[k - 1] } else { 0.0 };
    let right = if k + 1 < n { w[k + 1] - w[k] } else { 0.0 };
    left.max(right)
}

#[test]
fn cu_fits_keep_widths_above_the_knot_spacing_and_the_sum_rules_bounded() {
    let Some(elf) = cu_elf() else { return };
    let w = elf.energy_ev().to_vec();
    for n in [3, 6, 10, 16] {
        let fit = fit_mermin_oscillators(
            &elf,
            &MerminFitOptions {
                n_oscillators: n,
                ..Default::default()
            },
        )
        .unwrap();
        for o in &fit.oscillators {
            let floor = knot_spacing(&w, o.energy_ev);
            assert!(
                o.width_ev >= floor,
                "n = {n}: width {} eV below the knot spacing {floor} eV at {} eV\n{fit}",
                o.width_ev,
                o.energy_ev
            );
        }
        if n < 10 {
            continue;
        }
        // f-sum: measured 0.652 (n = 10) and 0.650 (n = 16) of the table's,
        // against 0.807 with 6 oscillators and 6771 / 16467 without the
        // floor. Bound: between half the table's f-sum and the table's.
        let ratio = fit.f_sum_ev2 / fit.data_f_sum_ev2;
        assert!(
            ratio > 0.5 && ratio <= 1.0,
            "n = {n}: f-sum ratio {ratio}\n{fit}"
        );
        // P_eff: measured 0.782 and 0.795 against the table's 1.0017 (22.7 and
        // 46.4 without the floor). Bound: within a factor 2 of the table's.
        let r = fit.p_eff / fit.data_p_eff;
        assert!(r > 0.5 && r < 2.0, "n = {n}: P_eff {}\n{fit}", fit.p_eff);
        // Mermin IMFP: measured within 4 % of the 6-oscillator values
        // (0.62 to 0.022 of them without the floor). Bound: within 10 %.
        let m = MerminPenn::new(elf.clone(), fit).unwrap();
        for (&(e, _), &six) in CU_TPP2011.iter().zip(&CU_MERMIN_6) {
            let l = m.imfp_m(e).unwrap() * 1e10;
            assert!(
                (l / six - 1.0).abs() < 0.10,
                "n = {n}, {e} eV: Mermin IMFP {l} vs {six} with 6 oscillators"
            );
        }
    }
}
