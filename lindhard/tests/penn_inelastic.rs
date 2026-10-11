//! Acceptance tests for the single-pole Penn model
//! (`lindhard::electron::inelastic`).
//!
//! Every ELF here is **synthetic**: a Drude-Lorentz model with chosen
//! parameters (`DrudeLorentz`), whose sum rules are known in closed form. No
//! optical data of a real material is used. The reference values are either
//! closed forms of the fixture (module docs of `inelastic::drude`) or limits
//! of the model itself (module docs of `inelastic::penn` and
//! `inelastic::bethe`).
//!
//! Kinematics are nonrelativistic throughout, as in the model.

use lindhard::electron::data::OpticalElf;
use lindhard::electron::inelastic::{
    DrudeLorentz, DrudeLorentzOscillator, SinglePolePenn, SumRuleReport,
};
use rayon::prelude::*;

/// Plasmon energy and width of the main fixture, eV (synthetic values).
const PLASMON_EV: f64 = 20.0;
const WIDTH_EV: f64 = 5.0;

fn plasmon_model() -> DrudeLorentz {
    DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(PLASMON_EV, WIDTH_EV)]).unwrap()
}

/// The fixture as a dense table: 4000 log-spaced samples from 1 meV to 100
/// keV. The truncation of the f-sum above 100 keV is `2γ/(π W_max)` relative
/// (3e-5) and the interpolation error near the peak about 1e-5.
fn dense_fixture() -> OpticalElf {
    plasmon_model()
        .to_optical_elf("synthetic Drude plasmon", 1e-3, 1e5, 4000)
        .unwrap()
}

/// A coarser table of the same model, for the costlier direct DIIMFP
/// integrations.
fn coarse_fixture() -> OpticalElf {
    plasmon_model()
        .to_optical_elf("synthetic Drude plasmon", 0.05, 2e4, 240)
        .unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

/// `n` log-spaced points from `lo` to `hi`.
fn log_grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| lo * (hi / lo).powf(i as f64 / (n - 1) as f64))
        .collect()
}

/// Composite Simpson rule in `ln x` on `[lo, hi]` with `n` (even) intervals.
fn simpson_log(f: &impl Fn(f64) -> f64, lo: f64, hi: f64, n: usize) -> f64 {
    let (a, b) = (lo.ln(), hi.ln());
    let h = (b - a) / n as f64;
    let mut s = 0.0;
    for i in 0..=n {
        let u = a + h * i as f64;
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let x = u.exp().clamp(lo, hi);
        s += w * f(x) * x;
    }
    s * h / 3.0
}

/// Composite Simpson rule on `[a, b]` with `m` (even) intervals.
fn simpson(f: &impl Fn(f64) -> f64, a: f64, b: f64, m: usize) -> f64 {
    let h = (b - a) / m as f64;
    let mut s = f(a) + f(b);
    for i in 1..m {
        s += if i % 2 == 1 { 4.0 } else { 2.0 } * f(a + h * i as f64);
    }
    s * h / 3.0
}

/// `∫ f(ω) dω` from the first table knot to `t`, for a DIIMFP-like `f`:
/// Simpson's rule on each interval between table knots (the DIIMFP has kinks
/// only at the knots), and on the last interval in `s = sqrt(t - ω)`
/// (`dω = 2s ds`), which smooths the `sqrt(t - ω)` behaviour at `ω -> t`.
fn integrate_losses(f: impl Fn(f64) -> f64, knots: &[f64], t: f64) -> f64 {
    let mut b: Vec<f64> = knots.iter().copied().filter(|&w| w < t).collect();
    b.push(t);
    let n = b.len();
    let mut total = 0.0;
    for k in 0..n - 2 {
        total += simpson(&f, b[k], b[k + 1], 8);
    }
    let a = b[n - 2];
    total
        + simpson(
            &|s: f64| f((t - s * s).max(a)) * 2.0 * s,
            0.0,
            (t - a).sqrt(),
            32,
        )
}

// ---------------------------------------------------------------------------
// Sum rules against the closed forms
// ---------------------------------------------------------------------------

#[test]
fn sum_rule_report_reproduces_the_drude_closed_forms() {
    let model = plasmon_model();
    let r = SumRuleReport::new(&dense_fixture());
    let f = rel(r.f_sum_ev2, model.f_sum_ev2());
    let ep = rel(r.plasma_energy_ev, model.plasma_energy_ev());
    let p = rel(r.p_eff, model.p_eff());
    assert!(f < 1e-3, "f-sum off by {f:e}");
    assert!(ep < 1e-3, "plasma energy off by {ep:e}");
    assert!(
        p < 1e-3,
        "P_eff = {} vs {} (rel {p:e})",
        r.p_eff,
        model.p_eff()
    );
    // A free-electron plasmon: P_eff = 1 and the plasma energy is the input.
    assert!((model.p_eff() - 1.0).abs() < 1e-15);
    assert!((model.plasma_energy_ev() - PLASMON_EV).abs() < 1e-12);
}

#[test]
fn sum_rules_of_a_two_oscillator_insulator_like_fixture() {
    // A bound oscillator (A != E^2): P_eff = sum A/E^2 < 1 here.
    let model = DrudeLorentz::new(vec![
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
    let elf = model
        .to_optical_elf("synthetic two-oscillator", 1e-3, 1e6, 6000)
        .unwrap();
    let r = SumRuleReport::new(&elf);
    assert!(
        rel(r.f_sum_ev2, model.f_sum_ev2()) < 1e-3,
        "{} vs {}",
        r.f_sum_ev2,
        model.f_sum_ev2()
    );
    assert!(
        rel(r.p_eff, model.p_eff()) < 1e-3,
        "{} vs {}",
        r.p_eff,
        model.p_eff()
    );
    assert!(model.p_eff() < 0.75);
}

#[test]
fn effective_electron_count_with_a_target_density() {
    let model = plasmon_model();
    let elf = dense_fixture();
    let n_e = SumRuleReport::new(&elf).electron_density_per_m3;
    // Closed form: the electron density whose plasma energy is E_p,
    // n = E_p^2 eps0 m_e / hbar^2 (E_p in J).
    use lindhard::constants::{ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY};
    let ep_j = model.plasma_energy_ev() * ELEMENTARY_CHARGE;
    let n_closed = ep_j * ep_j * VACUUM_PERMITTIVITY * ELECTRON_MASS
        / (HBAR * HBAR * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE);
    assert!(rel(n_e, n_closed) < 1e-3, "{n_e:e} vs {n_closed:e}");

    // Four electrons per target unit, by construction.
    let r = SumRuleReport::with_target_density(&elf, n_closed / 4.0).unwrap();
    let neff = r.effective_electrons().unwrap();
    assert!((neff / 4.0 - 1.0).abs() < 1e-3, "N_eff = {neff}");
    // N_eff(W) is non-decreasing, exact at the knots, and saturates.
    let mut last = 0.0;
    for w in log_grid(1e-3, 1e5, 300) {
        let x = r.effective_electrons_up_to(w).unwrap();
        assert!(x >= last - 1e-12, "N_eff({w}) = {x} < {last}");
        last = x;
    }
    for (k, w) in r.energy_ev.iter().enumerate().step_by(97) {
        let a = r.electron_density_up_to(*w);
        let b = r.cumulative_electron_density_per_m3[k];
        assert!((a - b).abs() <= 1e-12 * b.max(1.0), "{w}: {a} vs {b}");
    }
    // Most of the oscillator strength lies below 10 E_p: the fraction above
    // W is about 2γ/(πW) for W >> E_p, i.e. 1.6 % at W = 10 E_p.
    let frac = r.effective_electrons_up_to(10.0 * PLASMON_EV).unwrap() / neff;
    assert!(frac > 0.97 && frac < 0.995, "{frac}");
    assert!(SumRuleReport::with_target_density(&elf, 0.0).is_err());
    assert!(r.effective_electrons_up_to(0.0).unwrap() == 0.0);
    let text = r.to_string();
    assert!(text.contains("P_eff") && text.contains("N_eff"));
}

#[test]
fn mean_excitation_energy_matches_direct_quadrature() {
    let model = plasmon_model();
    let r = SumRuleReport::new(&dense_fixture());
    let num = simpson_log(&|w| w * model.elf(w) * w.ln(), 1e-6, 1e8, 40000);
    let den = simpson_log(&|w| w * model.elf(w), 1e-6, 1e8, 40000);
    let want = (num / den).exp();
    assert!(
        rel(r.mean_excitation_energy_ev, want) < 1e-3,
        "I = {} vs {want}",
        r.mean_excitation_energy_ev
    );
}

// ---------------------------------------------------------------------------
// The momentum-dependent loss function
// ---------------------------------------------------------------------------

#[test]
fn loss_function_reduces_to_the_optical_elf_at_q_zero() {
    let elf = dense_fixture();
    let penn = SinglePolePenn::new(elf.clone());
    for w in [0.5, 10.0, 20.0, 37.0, 900.0] {
        let want = elf.elf(w).unwrap();
        assert!(rel(penn.loss_function(0.0, w), want) < 1e-12, "{w}");
        // and continuously so: q = 1e-3 / a0
        let q = 1e-3 / lindhard::constants::BOHR_RADIUS;
        assert!(rel(penn.loss_function(q, w), want) < 1e-4, "{w}");
    }
}

#[test]
fn loss_function_keeps_the_f_sum_at_every_momentum() {
    // Each plasmon pole carries the f-sum of its Lindhard function, so
    // int W Im[-1/eps(q, W)] dW is the optical f-sum for every q (S2017
    // eqs. (4), (5), (7)). This checks the Jacobian of eq. (13).
    let elf = dense_fixture();
    let f_opt = SumRuleReport::new(&elf).f_sum_ev2;
    let penn = SinglePolePenn::new(elf);
    let a0 = lindhard::constants::BOHR_RADIUS;
    for q_au in [0.05, 0.5, 2.0, 6.0] {
        let q = q_au / a0;
        let f = simpson_log(&|w| w * penn.loss_function(q, w), 1e-3, 3e5, 60000);
        // the f-sum above the table maximum is missing on both sides
        assert!(rel(f, f_opt) < 2e-3, "q = {q_au} a.u.: {f} vs {f_opt}");
    }
}

// ---------------------------------------------------------------------------
// DIIMFP, IMFP and stopping
// ---------------------------------------------------------------------------

#[test]
fn imfp_and_stopping_are_the_integrals_of_the_diimfp() {
    let penn = SinglePolePenn::new(coarse_fixture());
    for e in [60.0, 2000.0] {
        let pt = penn.imfp_and_stopping(e).unwrap();
        let p = |w: f64| penn.diimfp_per_m_ev(e, w).unwrap();
        let knots = penn.optical_elf().energy_ev();
        let inv = integrate_losses(p, knots, e);
        let s = integrate_losses(|w| w * p(w), knots, e);
        let (di, ds) = (
            rel(inv, pt.inverse_imfp_per_m),
            rel(s, pt.stopping_ev_per_m),
        );
        assert!(di < 1e-5, "{e} eV: {inv:e} vs {:e}", pt.inverse_imfp_per_m);
        assert!(ds < 1e-5, "{e} eV: {s:e} vs {:e}", pt.stopping_ev_per_m);
    }
    assert_eq!(penn.diimfp_per_m_ev(100.0, 150.0).unwrap(), 0.0);
    assert_eq!(penn.diimfp_per_m_ev(100.0, 0.0).unwrap(), 0.0);
    assert_eq!(penn.diimfp_per_m_ev(100.0, -1.0).unwrap(), 0.0);
    assert!(penn.diimfp_per_m_ev(f64::NAN, 1.0).is_err());
}

#[test]
fn integration_error_is_below_the_stated_tolerance() {
    let elf = dense_fixture();
    let coarse = SinglePolePenn::new(elf.clone());
    let fine = SinglePolePenn::new(elf)
        .with_relative_tolerance(1e-11)
        .unwrap();
    for e in [15.0, 80.0, 1e3, 3e4] {
        let a = coarse.imfp_and_stopping(e).unwrap();
        let b = fine.imfp_and_stopping(e).unwrap();
        let tol = 2.0 * coarse.relative_tolerance();
        assert!(rel(a.inverse_imfp_per_m, b.inverse_imfp_per_m) < tol, "{e}");
        assert!(rel(a.stopping_ev_per_m, b.stopping_ev_per_m) < tol, "{e}");
    }
}

#[test]
fn stopping_approaches_nonrelativistic_bethe_within_2_percent() {
    let penn = SinglePolePenn::new(dense_fixture());
    for e in [1e4, 2e4, 5e4] {
        let s = penn.stopping_power_ev_per_m(e).unwrap();
        let sb = penn.bethe_stopping_ev_per_m(e).unwrap();
        let d = rel(s, sb);
        assert!(d < 0.02, "{e} eV: S = {s:e}, Bethe {sb:e} (rel {d:e})");
    }
}

#[test]
fn imfp_approaches_the_bethe_asymptotic_form_within_3_percent() {
    let penn = SinglePolePenn::new(dense_fixture());
    for e in [1e4, 2e4, 5e4] {
        let inv = penn.imfp_and_stopping(e).unwrap().inverse_imfp_per_m;
        let inv_b = penn.bethe_inverse_imfp_per_m(e).unwrap();
        let d = rel(inv, inv_b);
        assert!(
            d < 0.03,
            "{e} eV: 1/lambda = {inv:e}, Bethe {inv_b:e} (rel {d:e})"
        );
    }
    // The Bethe slope is fixed by the optical ELF alone (S2017 eqs. (28)-(29)):
    // E/lambda -> A ln E + const with A = (1/2pi) int ELF dW (per a0, in Ha).
    let (e1, e2): (f64, f64) = (2e4, 5e4);
    let a = (e2 * penn.bethe_inverse_imfp_per_m(e2).unwrap()
        - e1 * penn.bethe_inverse_imfp_per_m(e1).unwrap())
        / (e2 / e1).ln();
    let elf = penn.optical_elf();
    let hartree = lindhard::constants::HARTREE_ENERGY / lindhard::constants::ELEMENTARY_CHARGE;
    let int_elf = simpson_log(&|w| elf.elf(w).unwrap_or(0.0), 1e-3, 1e5, 40000) / hartree;
    let a_want =
        int_elf / (2.0 * std::f64::consts::PI) / lindhard::constants::BOHR_RADIUS * hartree;
    assert!(rel(a, a_want) < 1e-4, "{a:e} vs {a_want:e}");
}

#[test]
fn low_energy_imfp_is_finite_and_positive_with_a_minimum_between_30_and_300_ev() {
    let penn = SinglePolePenn::new(dense_fixture());
    let grid = log_grid(1e-3, 5e4, 160);
    let table = penn.tabulate(&grid).unwrap();
    let first = table
        .iter()
        .position(|p| p.inverse_imfp_per_m > 0.0)
        .expect("loss is allowed somewhere");
    assert!(grid[first] < 0.01, "first allowed energy {}", grid[first]);
    for p in &table[first..] {
        let l = p.imfp_m();
        assert!(l.is_finite() && l > 0.0, "{} eV: lambda = {l}", p.energy_ev);
        assert!(p.stopping_ev_per_m.is_finite() && p.stopping_ev_per_m > 0.0);
    }
    for p in &table[..first] {
        assert!(p.imfp_m().is_infinite());
    }
    let min = table[first..]
        .iter()
        .min_by(|a, b| a.imfp_m().total_cmp(&b.imfp_m()))
        .unwrap();
    assert!(
        min.energy_ev > 30.0 && min.energy_ev < 300.0,
        "IMFP minimum at {} eV",
        min.energy_ev
    );
}

#[test]
fn tables_are_bit_identical_across_thread_counts() {
    let penn = SinglePolePenn::new(coarse_fixture());
    let grid = log_grid(10.0, 5e4, 40);
    let serial = penn.tabulate(&grid).unwrap();
    for threads in [1, 2, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let parallel: Vec<_> = pool.install(|| {
            grid.par_iter()
                .map(|&e| penn.imfp_and_stopping(e).unwrap())
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

#[test]
fn fermi_energy_limits_losses_to_the_energy_above_the_fermi_level() {
    let penn = SinglePolePenn::new(coarse_fixture())
        .with_fermi_energy_ev(10.0)
        .unwrap();
    assert_eq!(penn.fermi_energy_ev(), 10.0);
    assert_eq!(penn.diimfp_per_m_ev(30.0, 30.5).unwrap(), 0.0);
    assert!(penn.diimfp_per_m_ev(30.0, 25.0).unwrap() > 0.0);
    // The stopping power is the integral of the DIIMFP up to T, not T'.
    let pt = penn.imfp_and_stopping(30.0).unwrap();
    let s = integrate_losses(
        |w| w * penn.diimfp_per_m_ev(30.0, w).unwrap(),
        penn.optical_elf().energy_ev(),
        30.0,
    );
    assert!(
        rel(s, pt.stopping_ev_per_m) < 1e-5,
        "{s:e} vs {:e}",
        pt.stopping_ev_per_m
    );
    assert!(SinglePolePenn::new(coarse_fixture())
        .with_fermi_energy_ev(-1.0)
        .is_err());
}

#[test]
fn rejects_invalid_inputs() {
    let penn = SinglePolePenn::new(coarse_fixture());
    assert!(penn.imfp_and_stopping(0.0).is_err());
    assert!(penn.imfp_and_stopping(f64::INFINITY).is_err());
    assert!(penn.bethe_inverse_imfp_per_m(-5.0).is_err());
    assert!(SinglePolePenn::new(coarse_fixture())
        .with_relative_tolerance(0.5)
        .is_err());
    assert_eq!(penn.loss_function(f64::NAN, 10.0), 0.0);
}
