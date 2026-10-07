//! Acceptance tests for the full Penn algorithm and its comparison with the
//! single-pole model (`lindhard::electron::inelastic::{FullPenn,
//! SinglePolePenn}`).
//!
//! As in `penn_inelastic.rs`, every ELF is the **synthetic** Drude-Lorentz
//! plasmon (E_p = 20 eV, width 5 eV; not a material), and kinematics are
//! nonrelativistic. The reference values are limits and sum rules of the
//! models themselves, and the single-pole model, whose own tests are in
//! `penn_inelastic.rs`. The numerical tolerance of the full model is set to
//! 1e-3 here to keep the (nested, three-dimensional) integrals affordable in
//! an unoptimised test build; its default is 1e-4.

use lindhard::constants::BOHR_RADIUS;
use lindhard::electron::data::OpticalElf;
use lindhard::electron::inelastic::{
    DrudeLorentz, DrudeLorentzOscillator, FullPenn, PennAlgorithm, PennInelastic, SinglePolePenn,
    SumRuleReport,
};
use rayon::prelude::*;

const TOL: f64 = 1e-3;

fn fixture() -> OpticalElf {
    DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])
        .unwrap()
        .to_optical_elf("synthetic Drude plasmon", 0.05, 2e4, 240)
        .unwrap()
}

fn full() -> FullPenn {
    FullPenn::new(fixture())
        .with_relative_tolerance(TOL)
        .unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

/// Composite Simpson rule in `ln x` on `[lo, hi]` with `n` (even) intervals.
fn simpson_log(f: &impl Fn(f64) -> f64, lo: f64, hi: f64, n: usize) -> f64 {
    let (a, b) = (lo.ln(), hi.ln());
    let h = (b - a) / n as f64;
    let mut s = 0.0;
    for i in 0..=n {
        let x = (a + h * i as f64).exp().clamp(lo, hi);
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        s += w * f(x) * x;
    }
    s * h / 3.0
}

// ---------------------------------------------------------------------------
// The extended ELF
// ---------------------------------------------------------------------------

#[test]
fn loss_function_reduces_to_the_optical_elf_at_q_zero() {
    let elf = fixture();
    let fpa = full();
    for w in [0.5, 10.0, 20.0, 37.0, 900.0] {
        assert!(
            rel(fpa.loss_function(0.0, w), elf.elf(w).unwrap()) < 1e-12,
            "{w}"
        );
    }
    // and continuously so: small q (the plasmon barely disperses)
    let q = 1e-2 / BOHR_RADIUS;
    for w in [15.0, 20.0, 25.0] {
        let want = elf.elf(w).unwrap();
        assert!(rel(fpa.loss_function(q, w), want) < 2e-2, "{w}");
    }
}

#[test]
fn expanded_elf_keeps_the_optical_f_sum_at_every_momentum() {
    // Each Lindhard gas carries (pi/2) w_p^2 at every q, so the expanded ELF
    // keeps the optical f-sum (module docs of `full_penn`).
    let f_opt = SumRuleReport::new(&fixture()).f_sum_ev2;
    let fpa = full();
    for q_au in [0.05, 0.3, 0.5, 1.0, 2.0, 6.0] {
        let f = fpa.f_sum_ev2_at(q_au / BOHR_RADIUS).unwrap();
        assert!(rel(f, f_opt) < 1e-3, "q = {q_au} a.u.: {f} vs {f_opt}");
    }
}

#[test]
fn integrating_the_loss_function_over_omega_gives_the_f_sum() {
    // The same sum rule through `loss_function` itself (plasmon term and
    // single-electron term of S2017 eq. (6) together), integrated over the
    // energy loss.
    let f_opt = SumRuleReport::new(&fixture()).f_sum_ev2;
    let fpa = FullPenn::new(fixture())
        .with_relative_tolerance(1e-2)
        .unwrap();
    for q_au in [0.4, 1.5] {
        let q = q_au / BOHR_RADIUS;
        let f = simpson_log(&|w| w * fpa.loss_function(q, w), 0.05, 2e4, 1200);
        assert!(rel(f, f_opt) < 1e-3, "q = {q_au} a.u.: {f} vs {f_opt}");
    }
}

#[test]
fn full_loss_function_is_wider_than_the_single_pole_at_large_q() {
    // The single pole sits on one curve w = w_q; the Lindhard continuum
    // spreads the strength (S2017 section on the ELF at q > 0).
    let (fpa, spa) = (full(), SinglePolePenn::new(fixture()));
    let q = 2.0 / BOHR_RADIUS;
    assert_eq!(spa.loss_function(q, 25.0), 0.0);
    assert!(fpa.loss_function(q, 25.0) > 0.0);
}

// ---------------------------------------------------------------------------
// DIIMFP, IMFP and stopping
// ---------------------------------------------------------------------------

#[test]
fn imfp_and_stopping_are_the_integrals_of_the_diimfp() {
    let fpa = full();
    let e = 300.0;
    let pt = fpa.imfp_and_stopping(e).unwrap();
    // Simpson in s = sqrt(e - w)? The DIIMFP has a plasmon peak and a
    // sqrt-type end; use a fine log grid plus a final sqrt-mapped piece.
    let split = 0.5 * e;
    let p = |w: f64| fpa.diimfp_per_m_ev(e, w).unwrap();
    let lo = 0.05;
    let part_lo = simpson_log(&p, lo, split, 120);
    let part_lo_s = simpson_log(&|w| w * p(w), lo, split, 120);
    let (a, b) = (0.0, (e - split).sqrt());
    let m = 24;
    let h = (b - a) / m as f64;
    let (mut hi, mut hi_s) = (0.0, 0.0);
    for i in 0..=m {
        let s = a + h * i as f64;
        let w = e - s * s;
        let c = if i == 0 || i == m {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let v = p(w) * 2.0 * s;
        hi += c * v;
        hi_s += c * v * w;
    }
    let (hi, hi_s) = (hi * h / 3.0, hi_s * h / 3.0);
    let (inv, s) = (part_lo + hi, part_lo_s + hi_s);
    assert!(
        rel(inv, pt.inverse_imfp_per_m) < 1e-2,
        "{inv:e} vs {:e}",
        pt.inverse_imfp_per_m
    );
    assert!(
        rel(s, pt.stopping_ev_per_m) < 1e-2,
        "{s:e} vs {:e}",
        pt.stopping_ev_per_m
    );
    assert_eq!(fpa.diimfp_per_m_ev(100.0, 150.0).unwrap(), 0.0);
    assert_eq!(fpa.diimfp_per_m_ev(100.0, 0.0).unwrap(), 0.0);
    assert!(fpa.diimfp_per_m_ev(f64::NAN, 1.0).is_err());
}

#[test]
fn full_and_single_pole_agree_within_1_percent_at_10_to_50_kev() {
    // Both reach the Bethe limit (same f-sum and mean excitation energy);
    // the constants of their inverse IMFPs differ only by the single-electron
    // continuum, which is below 1 % here.
    let (fpa, spa) = (full(), SinglePolePenn::new(fixture()));
    for e in [1e4, 2e4, 5e4] {
        let (a, b) = (
            fpa.imfp_and_stopping(e).unwrap(),
            spa.imfp_and_stopping(e).unwrap(),
        );
        assert!(
            rel(a.stopping_ev_per_m, b.stopping_ev_per_m) < 1e-2,
            "{e} eV: S {:e} vs {:e}",
            a.stopping_ev_per_m,
            b.stopping_ev_per_m
        );
        assert!(
            rel(a.inverse_imfp_per_m, b.inverse_imfp_per_m) < 1e-2,
            "{e} eV: 1/lambda {:e} vs {:e}",
            a.inverse_imfp_per_m,
            b.inverse_imfp_per_m
        );
        // and the Bethe stopping power of the same ELF
        let sb = spa.bethe_stopping_ev_per_m(e).unwrap();
        assert!(
            rel(a.stopping_ev_per_m, sb) < 2e-2,
            "{e} eV: {:e} vs {sb:e}",
            a.stopping_ev_per_m
        );
    }
}

#[test]
fn below_1_kev_the_models_differ_by_the_documented_amounts() {
    // Side by side (printed with --nocapture), and the differences bounded
    // as documented in the module docs of `full_penn`: at 100 eV to 1 keV
    // the stopping powers differ by at most 5 % and the inverse IMFPs by at
    // most 3 %; at 30 eV the full model has the larger inverse IMFP (the
    // single-electron excitations the single pole omits), by 1.4 to 2.5 times.
    let (fpa, spa) = (full(), SinglePolePenn::new(fixture()));
    println!("   E/eV   S_full/S_spa  invIMFP_full/invIMFP_spa");
    for e in [30.0, 50.0, 100.0, 200.0, 500.0, 1000.0] {
        let (a, b) = (
            fpa.imfp_and_stopping(e).unwrap(),
            spa.imfp_and_stopping(e).unwrap(),
        );
        let (rs, ri) = (
            a.stopping_ev_per_m / b.stopping_ev_per_m,
            a.inverse_imfp_per_m / b.inverse_imfp_per_m,
        );
        println!("{e:7.1}  {rs:12.4}  {ri:12.4}");
        if e >= 100.0 {
            assert!((rs - 1.0).abs() < 0.05, "{e} eV: S ratio {rs}");
            assert!((ri - 1.0).abs() < 0.03, "{e} eV: 1/lambda ratio {ri}");
        }
        if e == 30.0 {
            assert!(ri > 1.4 && ri < 2.5, "{e} eV: 1/lambda ratio {ri}");
        }
    }
}

#[test]
fn tables_are_bit_identical_across_thread_counts() {
    let fpa = full();
    let grid = [40.0, 300.0, 3e4];
    let serial = fpa.tabulate(&grid).unwrap();
    for threads in [1, 3] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let parallel: Vec<_> = pool.install(|| {
            grid.par_iter()
                .map(|&e| fpa.imfp_and_stopping(e).unwrap())
                .collect()
        });
        for (a, b) in serial.iter().zip(&parallel) {
            assert_eq!(
                a.inverse_imfp_per_m.to_bits(),
                b.inverse_imfp_per_m.to_bits()
            );
            assert_eq!(a.stopping_ev_per_m.to_bits(), b.stopping_ev_per_m.to_bits());
        }
    }
}

// ---------------------------------------------------------------------------
// Selection per material, recorded in the metadata
// ---------------------------------------------------------------------------

#[test]
fn the_algorithm_is_selectable_per_material_and_recorded() {
    let elf = fixture();
    let sp = PennInelastic::new(PennAlgorithm::SinglePole, elf.clone());
    let fp = PennInelastic::new(PennAlgorithm::Full, elf.clone());
    assert_eq!(sp.algorithm(), PennAlgorithm::SinglePole);
    assert_eq!(fp.algorithm(), PennAlgorithm::Full);
    let (ids, idf) = (sp.model_identity(), fp.model_identity());
    assert!(ids.starts_with("penn-single-pole;"), "{ids}");
    assert!(idf.starts_with("penn-full;"), "{idf}");
    assert_ne!(ids, idf);
    assert!(idf.contains("nonrelativistic") && idf.contains("E_F = 0 eV"));
    // the selector gives the underlying models' numbers
    let e = 5e3;
    assert_eq!(
        sp.imfp_and_stopping(e).unwrap(),
        SinglePolePenn::new(elf).imfp_and_stopping(e).unwrap()
    );
    let fermi = sp.with_fermi_energy_ev(10.0).unwrap();
    assert!(fermi.model_identity().contains("E_F = 10 eV"));
    assert_eq!(fermi.optical_elf().material(), "synthetic Drude plasmon");
}

#[test]
fn rejects_invalid_inputs() {
    let fpa = full();
    assert!(fpa.imfp_and_stopping(0.0).is_err());
    assert!(fpa.imfp_and_stopping(f64::NAN).is_err());
    assert!(fpa.diimfp_per_m_ev(-1.0, 1.0).is_err());
    assert!(fpa.f_sum_ev2_at(0.0).is_err());
    assert!(FullPenn::new(fixture())
        .with_relative_tolerance(0.5)
        .is_err());
    assert!(FullPenn::new(fixture()).with_fermi_energy_ev(-1.0).is_err());
    assert_eq!(fpa.loss_function(1e9, -1.0), 0.0);
}
