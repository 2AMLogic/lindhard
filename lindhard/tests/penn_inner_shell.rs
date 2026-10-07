//! Tests for the shell-resolved inner-shell channels and the optional
//! Born-Ochkur exchange correction of the single-pole Penn model
//! (`lindhard::electron::inelastic::{ShellResolvedChannels,
//! ExchangeCorrection, born_ochkur_factor}`).
//!
//! Every input is **synthetic**: a Drude-Lorentz valence ELF, made-up shell
//! ELFs, and a made-up `SubshellBindingTable` with arbitrary binding energies
//! (provenance strings say so). None is physical data; the real EADL table is
//! not needed.

use lindhard::constants::{ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR};
use lindhard::electron::data::{
    AtomBindings, ElectronDataError, OpticalElf, ShellBinding, Subshell, SubshellBindingTable,
};
use lindhard::electron::inelastic::{
    born_ochkur_factor, Channel, DrudeLorentz, DrudeLorentzOscillator, ExchangeCorrection,
    InnerShell, ShellResolvedChannels, SinglePolePenn,
};
use rayon::prelude::*;

const SYNTHETIC: &str = "synthetic test fixture, not physical data";

fn fixture_elf() -> OpticalElf {
    DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])
        .unwrap()
        .to_optical_elf("synthetic Drude plasmon", 0.05, 2e4, 240)
        .unwrap()
}

fn shell(label: &str, ebi: f64, occ: f64) -> ShellBinding {
    ShellBinding::new(Subshell::from_label(label).unwrap(), ebi, occ).unwrap()
}

/// A made-up 14-electron atom (the binding energies are arbitrary).
fn synthetic_table() -> SubshellBindingTable {
    let atom = AtomBindings::new(
        14,
        vec![
            shell("K", 1800.0, 2.0),
            shell("L1", 150.0, 2.0),
            shell("L2", 100.0, 2.0),
            shell("L3", 99.0, 4.0),
            shell("M1", 15.0, 2.0),
            shell("M2", 8.0, 1.0),
            shell("M3", 8.0, 1.0),
        ],
    )
    .unwrap();
    SubshellBindingTable::new(SYNTHETIC, vec![atom]).unwrap()
}

fn inner(label: &str) -> InnerShell {
    InnerShell::from_table(&synthetic_table(), 14, Subshell::from_label(label).unwrap()).unwrap()
}

/// A made-up shell ELF on the knots of the valence grid at and above `B`:
/// `c (1 - E0/E) (E0/E)^2.5`, `E0` the first such knot, so it is zero at its
/// first knot (and a sum with the valence table on the full grid is the
/// same piecewise-linear function).
fn shell_elf(b: f64, c: f64) -> OpticalElf {
    let grid: Vec<f64> = fixture_elf()
        .energy_ev()
        .iter()
        .copied()
        .filter(|&e| e >= b)
        .collect();
    let e0 = grid[0];
    let v = grid
        .iter()
        .map(|&e| c * (1.0 - e0 / e) * (e0 / e).powf(2.5))
        .collect();
    OpticalElf::new("synthetic shell", SYNTHETIC, grid, v).unwrap()
}

/// L3 (99 eV), L1 (150 eV) and K (1800 eV) with made-up ELFs.
fn shells() -> Vec<(InnerShell, OpticalElf)> {
    [("L3", 0.08), ("L1", 0.05), ("K", 0.02)]
        .iter()
        .map(|&(l, c)| {
            let s = inner(l);
            (s, shell_elf(s.binding_energy_ev, c))
        })
        .collect()
}

fn channels(valence: SinglePolePenn) -> ShellResolvedChannels {
    ShellResolvedChannels::new(valence, shells()).unwrap()
}

fn with_exchange() -> SinglePolePenn {
    SinglePolePenn::new(fixture_elf()).with_exchange(ExchangeCorrection::new(2e4).unwrap())
}

/// The valence ELF plus every shell ELF, on the valence grid.
fn summed_elf() -> OpticalElf {
    let v = fixture_elf();
    let parts: Vec<OpticalElf> = shells().into_iter().map(|(_, e)| e).collect();
    let sum: Vec<f64> = v
        .energy_ev()
        .iter()
        .zip(v.elf_values())
        .map(|(&e, &x)| {
            x + parts
                .iter()
                .map(|p| match p.energy_ev().iter().position(|&g| g == e) {
                    Some(k) => p.elf_values()[k],
                    None => 0.0,
                })
                .sum::<f64>()
        })
        .collect();
    OpticalElf::new("synthetic sum", SYNTHETIC, v.energy_ev().to_vec(), sum).unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

#[test]
fn inner_shells_take_their_binding_energy_from_the_table() {
    assert_eq!(inner("K").binding_energy_ev, 1800.0);
    assert_eq!(inner("L3").binding_energy_ev, 99.0);
    let t = synthetic_table();
    assert!(InnerShell::from_table(&t, 8, Subshell::from_label("K").unwrap()).is_err());
    assert!(InnerShell::from_table(&t, 14, Subshell::from_label("N1").unwrap()).is_err());
    let c = channels(SinglePolePenn::new(fixture_elf()));
    let eb: Vec<f64> = c
        .inner_shells()
        .iter()
        .map(|s| s.binding_energy_ev)
        .collect();
    assert_eq!(eb, vec![99.0, 150.0, 1800.0]);
}

/// Without exchange the SPA is linear in the ELF: the channels add up to the
/// model of the summed ELF.
#[test]
fn channels_add_up_to_the_summed_elf_model_and_respect_the_edges() {
    let c = channels(SinglePolePenn::new(fixture_elf()));
    let total = SinglePolePenn::new(summed_elf());
    for e in [90.0, 250.0, 4000.0, 8000.0] {
        let got = c.inverse_imfps(e).unwrap().total_per_m();
        let want = total.imfp_and_stopping(e).unwrap().inverse_imfp_per_m;
        assert!(rel(got, want) < 1e-6, "{e} eV: {got:e} vs {want:e}");
        let (l, d) = (37.0, c.diimfps_per_m_ev(e, 37.0).unwrap());
        let want = total.diimfp_per_m_ev(e, l).unwrap();
        assert!(rel(d.total_per_m_ev(), want) < 1e-6, "{e} eV");
    }
    // A loss cannot exceed T: shells whose edge is above T do not open.
    let r = c.inverse_imfps(145.0).unwrap();
    assert_eq!(r.shells_per_m[1], 0.0); // L1, 150 eV
    assert_eq!(r.shells_per_m[2], 0.0); // K
    let r = c.inverse_imfps(400.0).unwrap();
    assert!(r.shells_per_m[0] > 0.0 && r.shells_per_m[1] > 0.0); // L3, L1
    assert_eq!(r.shells_per_m[2], 0.0); // K, 1800 eV
    let r = c.inverse_imfps(60.0).unwrap();
    assert!(r.valence_per_m > 0.0 && r.shells_per_m.iter().all(|&x| x == 0.0));
    // The SPA pole needs T well above the edge (its dispersion): the K
    // shell (1800 eV) is still closed at 4 keV and open at 8 keV.
    assert_eq!(c.inverse_imfps(4000.0).unwrap().shells_per_m[2], 0.0);
    let r = c.inverse_imfps(8000.0).unwrap();
    assert!(r.shells_per_m.iter().all(|&x| x > 0.0));
    // Below its edge a shell's DIIMFP is zero.
    let d = c.diimfps_per_m_ev(500.0, 140.0).unwrap();
    assert!(d.shells_per_m_ev[0] > 0.0 && d.shells_per_m_ev[1] == 0.0);
    assert!(c.inverse_imfps(-1.0).is_err());
}

#[test]
fn secondary_energy_is_the_loss_minus_the_binding_energy() {
    let c = channels(SinglePolePenn::new(fixture_elf()));
    let k = Channel::InnerShell(2); // K, 1800 eV
    assert_eq!(c.secondary_energy_ev(k, 1900.0), Some(100.0));
    assert_eq!(c.secondary_energy_ev(k, 1800.0), Some(0.0));
    assert_eq!(c.secondary_energy_ev(k, 1799.0), None);
    assert_eq!(c.secondary_energy_ev(Channel::InnerShell(7), 1900.0), None);
    assert_eq!(c.secondary_energy_ev(Channel::Valence, 1900.0), None);
}

#[test]
fn sampling_follows_the_channel_diimfps() {
    let c = channels(SinglePolePenn::new(fixture_elf()));
    let (e, w) = (1000.0, 200.0);
    let d = c.diimfps_per_m_ev(e, w).unwrap();
    let total = d.total_per_m_ev();
    let n = 10_000;
    let mut counts = [0usize; 4];
    for k in 0..n {
        let u = (k as f64 + 0.5) / n as f64;
        match c.sample_channel(e, w, u).unwrap().unwrap() {
            Channel::Valence => counts[3] += 1,
            Channel::InnerShell(i) => counts[i] += 1,
        }
    }
    assert_eq!(counts[2], 0); // K closed at 200 eV
    for (i, p) in d.shells_per_m_ev.iter().enumerate().take(2) {
        assert!(
            (counts[i] as f64 / n as f64 - p / total).abs() < 1e-3,
            "{i}"
        );
    }
    assert!((counts[3] as f64 / n as f64 - d.valence_per_m_ev / total).abs() < 1e-3);
    // No loss is possible above T.
    assert_eq!(c.sample_channel(100.0, 150.0, 0.5).unwrap(), None);
}

/// The Born-Ochkur factor of de Vera et al. (2022), text below eq. (32):
/// `F = -x + x²`, `x = (k²/2m)/(T - W)`, `W = E - B` the emitted energy.
#[test]
fn born_ochkur_factor_uses_the_emitted_energy() {
    // T = 400 eV, B = 150 eV, E = 200 eV: W = 50 eV, T - W = 350 eV.
    let x: f64 = 175.0 / 350.0;
    assert_eq!(born_ochkur_factor(175.0, 400.0, 200.0, 150.0), -x + x * x);
    assert_eq!(born_ochkur_factor(175.0, 400.0, 200.0, 150.0), -0.25);
    // With B = 0 the denominator is T - E = 200 eV.
    let x0: f64 = 175.0 / 200.0;
    assert_eq!(born_ochkur_factor(175.0, 400.0, 200.0, 0.0), -x0 + x0 * x0);
    assert!(born_ochkur_factor(1.0, 100.0, 100.0, 0.0).is_nan());
}

/// `ħ²/2m` in eV m², from CODATA SI constants (independent of the atomic
/// units used inside the model).
fn recoil_coefficient_ev_m2() -> f64 {
    HBAR * HBAR / (2.0 * ELECTRON_MASS) / ELEMENTARY_CHARGE
}

/// `∫ dq/q L(q, ω) (1 + F(q))` over `q± = (sqrt(2mT) ± sqrt(2m(T - ω)))/ħ`
/// by composite Simpson in `ln q`, with `F` built on the denominator `d_ev`
/// (`None`: no exchange). Independent of the model's own quadrature.
fn q_integral(model: &SinglePolePenn, t: f64, w: f64, d_ev: Option<f64>) -> f64 {
    let a = recoil_coefficient_ev_m2();
    let (qm, qp) = (
        (t.sqrt() - (t - w).sqrt()) / a.sqrt(),
        (t.sqrt() + (t - w).sqrt()) / a.sqrt(),
    );
    let n = 200_000;
    let (u0, u1) = (qm.ln(), qp.ln());
    let h = (u1 - u0) / n as f64;
    let f = |u: f64| {
        let q = u.exp();
        let l = model.loss_function(q, w);
        match d_ev {
            Some(d) => {
                let x = a * q * q / d;
                l * (1.0 - x + x * x)
            }
            None => l,
        }
    };
    let mut s = f(u0) + f(u1);
    for k in 1..n {
        s += f(u0 + k as f64 * h) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    s * h / 3.0
}

/// Regression for the inner-shell exchange: at T = 400 eV, B = 150 eV and a
/// loss of 200 eV the exchange-corrected DIIMFP must use `T - W = 350 eV`
/// (eq. (32)), not `T - E = 200 eV`. Checked against an independent `q`
/// integral of the closed form.
#[test]
fn inner_shell_exchange_diimfp_matches_the_closed_form_with_nonzero_binding() {
    let (t, w, b) = (400.0, 200.0, 150.0);
    let plain = SinglePolePenn::new(fixture_elf());
    let ex = with_exchange();
    let ratio =
        ex.diimfp_with_binding_per_m_ev(t, w, b).unwrap() / plain.diimfp_per_m_ev(t, w).unwrap();
    let direct = q_integral(&plain, t, w, None);
    let want = q_integral(&plain, t, w, Some(t - (w - b))) / direct;
    let wrong = q_integral(&plain, t, w, Some(t - w)) / direct;
    assert!(rel(ratio, want) < 1e-6, "{ratio} vs {want}");
    // The test discriminates: the loss-based denominator is far off.
    assert!(rel(wrong, want) > 1e-2, "{wrong} vs {want}");
    // B = 0 is the plain model's exchange.
    let r0 = ex.diimfp_per_m_ev(t, w).unwrap() / plain.diimfp_per_m_ev(t, w).unwrap();
    assert!(rel(r0, wrong) < 1e-6, "{r0} vs {wrong}");
    assert_eq!(
        ex.diimfp_with_binding_per_m_ev(t, w, 0.0).unwrap(),
        ex.diimfp_per_m_ev(t, w).unwrap()
    );
    // Without exchange the binding energy does not enter.
    assert_eq!(
        plain.diimfp_with_binding_per_m_ev(t, w, b).unwrap(),
        plain.diimfp_per_m_ev(t, w).unwrap()
    );
    assert!(ex.diimfp_with_binding_per_m_ev(t, w, -1.0).is_err());
    assert!(ex.diimfp_with_binding_per_m_ev(t, w, f64::NAN).is_err());
}

/// The loss limit is `ω <= (T' + B)/2` (eq. (35) text), and `-1/4 <= F <= 0`
/// below it for every binding energy.
#[test]
fn inner_shell_exchange_limit_and_bound() {
    let plain = SinglePolePenn::new(fixture_elf());
    let ex = with_exchange();
    let (t, b) = (400.0, 150.0);
    assert_eq!(ex.diimfp_with_binding_per_m_ev(t, 276.0, b).unwrap(), 0.0);
    assert!(ex.diimfp_with_binding_per_m_ev(t, 274.0, b).unwrap() > 0.0);
    assert_eq!(ex.diimfp_with_binding_per_m_ev(t, 201.0, 0.0).unwrap(), 0.0);
    for w in [160.0, 200.0, 250.0, 275.0] {
        let a = plain.diimfp_per_m_ev(t, w).unwrap();
        let x = ex.diimfp_with_binding_per_m_ev(t, w, b).unwrap();
        assert!(x <= a && x >= 0.75 * a, "{w} eV: {x:e} vs {a:e}");
    }
}

#[test]
fn exchange_is_optional_and_off_by_default() {
    let plain = SinglePolePenn::new(fixture_elf());
    assert!(plain.exchange().is_none());
    assert!(!plain.exchange_applies_at(50.0));
    let ex =
        SinglePolePenn::new(fixture_elf()).with_exchange(ExchangeCorrection::new(300.0).unwrap());
    assert_eq!(ex.exchange().unwrap().applies_below_ev(), 300.0);
    assert!(ex.exchange_applies_at(299.0) && !ex.exchange_applies_at(300.0));
    // At and above the stated energy the model is the plain one, bit for bit.
    for e in [300.0, 1000.0] {
        assert_eq!(
            ex.imfp_and_stopping(e).unwrap(),
            plain.imfp_and_stopping(e).unwrap()
        );
    }
    assert_ne!(
        ex.imfp_m(100.0).unwrap().to_bits(),
        plain.imfp_m(100.0).unwrap().to_bits()
    );
    assert!(ExchangeCorrection::new(0.0).is_err());
    assert!(ExchangeCorrection::new(f64::NAN).is_err());
}

/// The effect on lambda(E): large below 200 eV, falling with energy, and
/// under 0.5 % at 10 keV (it is 0.17 % there for this fixture).
#[test]
fn exchange_lengthens_the_imfp_at_low_energy_and_vanishes_at_high_energy() {
    let plain = SinglePolePenn::new(fixture_elf());
    let ex = with_exchange();
    let ratio = |e: f64| ex.imfp_m(e).unwrap() / plain.imfp_m(e).unwrap();
    let low: Vec<f64> = [30.0, 50.0, 100.0, 150.0]
        .iter()
        .map(|&e| ratio(e))
        .collect();
    // Exchange removes inverse path: lambda grows, by more than 10 % below
    // 150 eV, and the effect falls with energy.
    assert!(low.iter().all(|&r| r > 1.10), "{low:?}");
    assert!(low.windows(2).all(|w| w[0] > w[1]), "{low:?}");
    let r200 = ratio(200.0);
    assert!(r200 > 1.05 && r200 < low[3], "{r200}");
    // Vanishes above 10 keV to 0.5 %.
    for e in [1e4, 1.5e4] {
        let r = ratio(e);
        assert!((r - 1.0).abs() < 5e-3, "{e} eV: {r}");
    }
    // The stopping power is reduced too, never increased.
    for e in [50.0, 200.0, 5000.0] {
        let a = plain.stopping_power_ev_per_m(e).unwrap();
        let b = ex.stopping_power_ev_per_m(e).unwrap();
        assert!(b < a && b > 0.0, "{e}: {b} vs {a}");
    }
}

#[test]
fn exchange_diimfp_is_cut_at_half_the_energy_and_never_larger() {
    let plain = SinglePolePenn::new(fixture_elf());
    let ex = with_exchange();
    let t = 100.0;
    assert_eq!(ex.diimfp_per_m_ev(t, 51.0).unwrap(), 0.0);
    for w in [5.0, 20.0, 40.0, 50.0] {
        let (a, b) = (
            plain.diimfp_per_m_ev(t, w).unwrap(),
            ex.diimfp_per_m_ev(t, w).unwrap(),
        );
        // -1/4 <= F <= 0 for losses up to T/2.
        assert!(b <= a && b >= 0.75 * a, "{w} eV: {b:e} vs {a:e}");
    }
}

#[test]
fn channels_with_exchange_are_smaller_and_use_their_binding_energy() {
    let a = channels(SinglePolePenn::new(fixture_elf()));
    let b = channels(with_exchange());
    let e = 400.0;
    let (ra, rb) = (a.inverse_imfps(e).unwrap(), b.inverse_imfps(e).unwrap());
    assert!(rb.valence_per_m < ra.valence_per_m);
    // L3 and L1 are open at 400 eV and lose strength; K is closed.
    assert!(rb.shells_per_m[..2]
        .iter()
        .zip(&ra.shells_per_m)
        .all(|(x, y)| x < y));
    assert_eq!(rb.shells_per_m[2], 0.0);
    // Each shell channel is its model's exchange integral with its own B,
    // which differs from the B = 0 one (larger denominator, higher limit).
    for i in 0..2 {
        let m = b.shell_model(i);
        let eb = b.inner_shells()[i].binding_energy_ev;
        let with_b = m
            .imfp_and_stopping_with_binding(e, eb)
            .unwrap()
            .inverse_imfp_per_m;
        assert_eq!(rb.shells_per_m[i], with_b);
        let without = m.imfp_and_stopping(e).unwrap().inverse_imfp_per_m;
        assert!(with_b > without, "{i}: {with_b:e} vs {without:e}");
    }
}

#[test]
fn results_are_bit_identical_across_thread_counts() {
    let c = channels(with_exchange());
    let energies = [40.0, 90.0, 130.0, 300.0];
    let seq: Vec<_> = energies
        .iter()
        .map(|&e| c.inverse_imfps(e).unwrap())
        .collect();
    for threads in [1usize, 2] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let par: Vec<_> = pool.install(|| {
            energies
                .par_iter()
                .map(|&e| c.inverse_imfps(e).unwrap())
                .collect()
        });
        assert_eq!(seq, par, "{threads} threads");
    }
}

#[test]
fn rejects_invalid_inputs() {
    let v = || SinglePolePenn::new(fixture_elf());
    let l3 = inner("L3");
    // A shell ELF that starts below its edge.
    assert!(matches!(
        ShellResolvedChannels::new(v(), vec![(l3, shell_elf(50.0, 0.1))]),
        Err(ElectronDataError::Invalid { .. })
    ));
    // The same shell twice.
    assert!(ShellResolvedChannels::new(
        v(),
        vec![(l3, shell_elf(99.0, 0.1)), (l3, shell_elf(99.0, 0.1))]
    )
    .is_err());
    // A non-positive or non-finite binding energy.
    for b in [0.0, -5.0, f64::NAN] {
        let s = InnerShell {
            binding_energy_ev: b,
            ..l3
        };
        assert!(ShellResolvedChannels::new(v(), vec![(s, shell_elf(99.0, 0.1))]).is_err());
    }
    // No inner shell at all is allowed (valence only).
    let c = ShellResolvedChannels::new(v(), vec![]).unwrap();
    let r = c.inverse_imfps(300.0).unwrap();
    assert_eq!(
        r.total_per_m(),
        v().imfp_and_stopping(300.0).unwrap().inverse_imfp_per_m
    );
}
