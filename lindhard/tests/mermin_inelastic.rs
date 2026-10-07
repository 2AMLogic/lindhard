//! Acceptance tests for the Mermin-ELF model (`lindhard::electron::inelastic::
//! {MerminPenn, MerminGas, fit_mermin_oscillators}`).
//!
//! Every ELF is **synthetic** (Drude-Lorentz oscillators with made-up
//! parameters; not a material), kinematics are nonrelativistic. The reference
//! values are limits and sum rules of the models themselves and the other two
//! models of the selector.

use lindhard::electron::data::OpticalElf;
use lindhard::electron::inelastic::{
    fit_mermin_oscillators, DrudeLorentz, DrudeLorentzOscillator, MerminFitOptions, MerminPenn,
    PennAlgorithm, PennInelastic, SinglePolePenn,
};

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
