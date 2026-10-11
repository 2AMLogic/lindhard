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

// ---- Cu and C: the width floor (#307) ----
//
// Without a lower bound on the widths, the 10- and 16-oscillator fits of the
// committed Cu ELF placed an oscillator of width 3e-4 eV (1e-4 eV) between
// the knots at 999.9 and 1100.1 eV (1100.1 and 1200.2 eV), which the
// residuals at the knots do not see and the sum rules do: f-sum 6771 (16467)
// times the table's. The fit now bounds every width below by the local knot
// spacing at the oscillator's energy (our choice, documented on
// `fit_mermin_oscillators`), and these tests check that bound. The same fits
// still place a very wide oscillator (4e5 to 8e5 eV) whose f-sum lies above
// the table, 188 and 106 times the table's (measured 2026-10-09, release
// build, unconverged): a different defect, under the default weighting, that
// the floor does not address (#311). No bound on the f-sum, `P_eff` or IMFP
// of those fits is asserted here.

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

fn c_elf() -> Option<OpticalElf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../validation/data/optical/c_elf_hagemann1975.toml");
    if !p.is_file() {
        eprintln!("{} not found; C Mermin checks skipped", p.display());
        return None;
    }
    Some(OpticalElf::from_toml_file(p).unwrap())
}

fn assert_widths_above_the_knot_spacing(elf: &OpticalElf, n: usize) {
    let w = elf.energy_ev();
    let fit = fit_mermin_oscillators(
        elf,
        &MerminFitOptions {
            n_oscillators: n,
            ..Default::default()
        },
    )
    .unwrap();
    for o in &fit.oscillators {
        let floor = knot_spacing(w, o.energy_ev);
        assert!(
            o.width_ev >= floor,
            "n = {n}: width {} eV below the knot spacing {floor} eV at {} eV\n{fit}",
            o.width_ev,
            o.energy_ev
        );
    }
}

#[test]
fn cu_fits_keep_widths_above_the_knot_spacing() {
    // fails on the code before #307 at n = 10 (width 3e-4 eV at 1052.6 eV,
    // where the knot spacing is 100.3 eV)
    let Some(elf) = cu_elf() else { return };
    for n in [3, 6, 10, 16] {
        assert_widths_above_the_knot_spacing(&elf, n);
    }
}

#[test]
fn c_default_fit_starts_on_the_floor_and_leaves_it() {
    // The one committed ELF where the floor binds at the default options:
    // the default start puts an oscillator on the 0.2 eV peak with width
    // 0.05 eV, under the 0.1 eV knot spacing there. The projection raises
    // it to the floor and the fit moves it off again, to the same fit as
    // without the floor: (2.000, 5.020, 28.457) eV, weighted rms 0.29194
    // (measured 2026-10-09 with and without the floor, which agree to about
    // 1e-7). A floor that froze the width gave (0.161, 2.81, 28.31) eV and
    // rms 0.3030.
    let Some(elf) = c_elf() else { return };
    assert_widths_above_the_knot_spacing(&elf, 3);
    let fit = fit_mermin_oscillators(&elf, &MerminFitOptions::default()).unwrap();
    assert!(fit.converged, "{fit}");
    for (o, want) in fit.oscillators.iter().zip([2.000, 5.020, 28.457]) {
        close(o.energy_ev, want, "C oscillator energy");
    }
    close(fit.weighted_rms, 0.29194, "C weighted rms");
}

// ---- Cu: the fit options of the electron input (#306) ----

/// A Cu electron input with `mermin-melf` and the given `mermin_fit` table
/// (none if `None`), resolved against the committed Cu ELF.
fn cu_input(fit: Option<&str>) -> Option<lindhard::input::electron::ResolvedElectron> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../validation/data/optical");
    if !dir.join("cu_elf_hagemann1975.toml").is_file() {
        eprintln!("{} not found; Cu input checks skipped", dir.display());
        return None;
    }
    let fit = fit.map_or(String::new(), |f| {
        format!("[electron.materials.Cu.mermin_fit]\n{f}\n")
    });
    let text = format!(
        "[electron.beam]\nenergy_ev = 1000.0\n[electron.transport]\ncutoff_ev = 20.0\n\
         [electron.elastic]\npotential = \"thomas-fermi-yukawa\"\n\
         [electron.inelastic]\nmodel = \"mermin-melf\"\n\
         [electron.materials.Cu]\noptical_elf = \"cu_elf_hagemann1975.toml\"\n{fit}\
         [target]\nsubstrate = \"Cu\"\n[run]\nhistories = 1\nseed = 1\n"
    );
    Some(
        lindhard::input::electron::ElectronInput::from_toml_str(&text)
            .unwrap()
            .resolve_in(&dir)
            .unwrap(),
    )
}

#[test]
fn cu_input_with_six_oscillators_resolves_and_fits_six() {
    let Some(r) = cu_input(Some("oscillators = 6")) else {
        return;
    };
    let m = &r.materials[0];
    assert_eq!(m.mermin_fit.n_oscillators, 6);
    let model = PennInelastic::mermin(m.optical_elf.clone(), &m.mermin_fit).unwrap();
    assert_eq!(model.algorithm(), PennAlgorithm::Mermin);
    let id = model.model_identity();
    assert!(id.contains("(6 oscillators "), "{id}");
    // Without the table the input resolves to the default options, and the
    // model is the one `try_new` builds.
    let r0 = cu_input(None).unwrap();
    assert_eq!(r0.materials[0].mermin_fit, MerminFitOptions::default());
    let default = PennInelastic::mermin(
        r0.materials[0].optical_elf.clone(),
        &r0.materials[0].mermin_fit,
    )
    .unwrap();
    let try_new =
        PennInelastic::try_new(PennAlgorithm::Mermin, r0.materials[0].optical_elf.clone()).unwrap();
    assert_eq!(default.model_identity(), try_new.model_identity());
    assert!(default.model_identity().contains("(3 oscillators "));
}

// ---------------------------------------------------------------------------
// The loss window of a Mermin table row (#342)
// ---------------------------------------------------------------------------

fn si_elf() -> Option<OpticalElf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../validation/data/optical/si_elf_yang2019.toml");
    if !p.is_file() {
        // A packaged crate does not ship validation/data.
        eprintln!("{} not found; Si Mermin table checks skipped", p.display());
        return None;
    }
    Some(OpticalElf::from_toml_file(p).unwrap())
}

/// The fraction of the model's inverse IMFP at `t` from losses in
/// `(w_lo, t]`: an independent trapezoid in `ln ω` of the public DIIMFP.
fn rate_fraction_above(m: &PennInelastic, t: f64, w_lo: f64) -> f64 {
    let (lo, hi, n) = (w_lo.ln(), t.ln(), 400usize);
    let h = (hi - lo) / n as f64;
    let mut sum = 0.0;
    for i in 0..=n {
        let w = (lo + h * i as f64).exp();
        let wt = if i == 0 || i == n { 0.5 } else { 1.0 };
        sum += wt * m.diimfp_per_m_ev(t, w).unwrap() * w;
    }
    sum * h / m.imfp_and_stopping(t).unwrap().inverse_imfp_per_m
}

/// Issue #342. The Mermin rate integrates the DIIMFP of its fitted ELF from
/// `ω = 1e-8 T`, below the lowest tabulated ELF energy `elf_min` (0.5 eV for
/// the committed Si ELF). The table's sampled loss distribution must cover
/// the same window: the fraction of draws below `elf_min` is the fraction of
/// the rate there, the stopping power of the stored rows is the model's, and
/// rows with `T <= elf_min` carry the model's rate. Before #342 the window
/// opened at `elf_min`: no draw fell below it, the row at `T = 0.6 eV` gave
/// 1.42 times the model's stopping power, and the row at 0.3 eV was empty.
#[test]
fn si_mermin_rows_sample_the_window_of_their_rate() {
    use lindhard::electron::inelastic::table::{
        build_inelastic_table_for_model, stopping_power_ev_per_m, InelasticTableOptions,
    };
    use lindhard::material::Material;
    use lindhard::rng::stream;
    use rand_core::Rng;

    let Some(elf) = si_elf() else { return };
    let elf_min = elf.energy_ev()[0];
    assert_eq!(elf_min, 0.5);
    // the band of validation/experiments/se_yield.py (E_F = W_v + E_g)
    let m = PennInelastic::try_new(PennAlgorithm::Mermin, elf)
        .unwrap()
        .with_fermi_energy_ev(13.46)
        .unwrap();
    let ts = vec![0.3, 0.6, 1.0, 2.0, 5.0];
    let mat = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    let table =
        build_inelastic_table_for_model(&m, &mat, &InelasticTableOptions::new(ts.clone())).unwrap();
    let n = 100_000u64;
    for (i, &t) in ts.iter().enumerate() {
        let model = m.imfp_and_stopping(t).unwrap();
        assert_eq!(
            table.inverse_mfp_per_m()[i],
            model.inverse_imfp_per_m,
            "T = {t}"
        );
        let q = table.quantiles(i).expect("non-empty row");
        // the window is (1e-8 T, T]
        assert!((q[0] / (1e-8 * t) - 1.0).abs() < 1e-12, "T = {t}: {}", q[0]);
        assert!(q.iter().all(|&w| w > 0.0 && w <= t), "T = {t}");
        // the stopping power of the row is the model's (density tolerance
        // 1e-5 and probability refinement 1e-4 of the mean loss)
        let s = stopping_power_ev_per_m(&table, i).unwrap();
        assert!(
            (s / model.stopping_ev_per_m - 1.0).abs() < 1e-3,
            "T = {t}: stopping {s} vs model {}",
            model.stopping_ev_per_m
        );
        // the fraction of draws below elf_min is the fraction of the rate
        // there (5 binomial standard deviations plus 1e-3 for the trapezoid)
        let below_rate = if t > elf_min {
            1.0 - rate_fraction_above(&m, t, elf_min)
        } else {
            1.0
        };
        let below = (0..n)
            .filter(|&k| {
                let u = (stream(342, k).next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64);
                table.inverse_cdf(i, u).unwrap() < elf_min
            })
            .count() as f64
            / n as f64;
        let sigma = (below_rate * (1.0 - below_rate) / n as f64).sqrt();
        eprintln!("T = {t}: draws below elf_min {below:.4}, rate below {below_rate:.4}");
        assert!(
            (below - below_rate).abs() <= 5.0 * sigma + 1e-3,
            "T = {t}: draws below elf_min {below} vs rate fraction {below_rate}"
        );
    }
}
