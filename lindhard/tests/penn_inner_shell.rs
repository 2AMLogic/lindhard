//! Tests for the inner-shell channel attribution and the optional Born-Ochkur
//! exchange correction of the single-pole Penn model
//! (`lindhard::electron::inelastic::{ChannelPartition, ExchangeCorrection}`).
//!
//! Every input is **synthetic**: a Drude-Lorentz ELF and a made-up
//! `SubshellBindingTable` with arbitrary binding energies (provenance string
//! says so). Neither is physical data; the real EADL table is not needed.

use lindhard::electron::data::{
    AtomBindings, ElectronDataError, OpticalElf, ShellBinding, Subshell, SubshellBindingTable,
};
use lindhard::electron::inelastic::{
    Channel, ChannelPartition, DrudeLorentz, DrudeLorentzOscillator, ExchangeCorrection,
    SinglePolePenn,
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

/// A made-up 14-electron atom: two valence subshells below 50 eV and four
/// deeper ones (the binding energies are arbitrary).
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

fn partition() -> ChannelPartition {
    ChannelPartition::unsourced_occupancy_weighted(&synthetic_table(), &[(14, 1.0)], 50.0).unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

#[test]
fn inner_shells_are_those_at_or_above_the_cutoff_in_energy_order() {
    let p = partition();
    let eb: Vec<f64> = p
        .inner_shells()
        .iter()
        .map(|s| s.binding_energy_ev)
        .collect();
    assert_eq!(eb, vec![99.0, 100.0, 150.0, 1800.0]);
    let labels: Vec<&str> = p
        .inner_shells()
        .iter()
        .map(|s| s.subshell.label())
        .collect();
    assert_eq!(labels, vec!["L3", "L2", "L1", "K"]);
}

#[test]
fn fractions_follow_the_occupancy_rule_and_sum_to_one() {
    let p = partition();
    // Below every edge: all valence.
    assert_eq!(p.valence_fraction(5.0), 1.0);
    for i in 0..4 {
        assert_eq!(p.shell_fraction(5.0, i), 0.0);
    }
    // At 20 eV only M1, M2, M3 (4 electrons) are open: still all valence.
    assert!((p.valence_fraction(20.0) - 1.0).abs() < 1e-15);
    // At 99.5 eV: open = 4 valence electrons + L3 (4): L3 takes 1/2.
    assert!((p.shell_fraction(99.5, 0) - 0.5).abs() < 1e-15);
    assert!((p.valence_fraction(99.5) - 0.5).abs() < 1e-15);
    assert_eq!(p.shell_fraction(99.5, 1), 0.0);
    // At 120 eV: 4 + 4 + 2 = 10 electrons.
    assert!((p.shell_fraction(120.0, 0) - 0.4).abs() < 1e-15);
    assert!((p.shell_fraction(120.0, 1) - 0.2).abs() < 1e-15);
    assert!((p.valence_fraction(120.0) - 0.4).abs() < 1e-15);
    // Above the K edge: all 14 electrons.
    assert!((p.shell_fraction(2000.0, 3) - 2.0 / 14.0).abs() < 1e-15);
    for w in [
        1.0, 20.0, 99.0, 99.5, 100.0, 149.0, 150.0, 1799.0, 1800.0, 5e3,
    ] {
        let s: f64 = p.valence_fraction(w) + (0..4).map(|i| p.shell_fraction(w, i)).sum::<f64>();
        assert!((s - 1.0).abs() < 1e-14, "{w}: {s}");
    }
}

#[test]
fn compound_weights_use_atoms_per_formula_unit() {
    let a = AtomBindings::new(
        1,
        vec![shell("K", 60.0, 1.0)], // arbitrary, with the cutoff below
    )
    .unwrap();
    let b = AtomBindings::new(2, vec![shell("K", 90.0, 2.0)]).unwrap();
    let t = SubshellBindingTable::new(SYNTHETIC, vec![a, b]).unwrap();
    // H2He: weights 2*1 (Z=1, 60 eV) and 1*2 (Z=2, 90 eV); cutoff 50: both inner.
    let p =
        ChannelPartition::unsourced_occupancy_weighted(&t, &[(1, 2.0), (2, 1.0)], 50.0).unwrap();
    assert_eq!(p.inner_shells().len(), 2);
    assert!((p.shell_fraction(70.0, 0) - 1.0).abs() < 1e-15);
    assert!((p.shell_fraction(100.0, 0) - 0.5).abs() < 1e-15);
    assert!((p.shell_fraction(100.0, 1) - 0.5).abs() < 1e-15);
}

#[test]
fn secondary_energy_is_the_loss_minus_the_binding_energy() {
    let p = partition();
    let k = Channel::InnerShell(3); // K, 1800 eV
    assert_eq!(p.secondary_energy_ev(k, 1900.0), Some(100.0));
    assert_eq!(p.secondary_energy_ev(k, 1800.0), Some(0.0));
    assert_eq!(p.secondary_energy_ev(k, 1799.0), None);
    assert_eq!(p.secondary_energy_ev(Channel::Valence, 1900.0), None);
}

#[test]
fn sampling_reproduces_the_fractions() {
    let p = partition();
    let n = 10_000;
    let w = 120.0;
    let mut counts = [0usize; 5];
    for k in 0..n {
        let u = (k as f64 + 0.5) / n as f64;
        match p.sample_channel(w, u) {
            Channel::Valence => counts[4] += 1,
            Channel::InnerShell(i) => counts[i] += 1,
        }
    }
    // Shells L3, L2 open at 120 eV; L1 (150) and K are closed.
    assert_eq!(counts[2], 0);
    assert_eq!(counts[3], 0);
    assert!((counts[0] as f64 / n as f64 - 0.4).abs() < 1e-3);
    assert!((counts[1] as f64 / n as f64 - 0.2).abs() < 1e-3);
    assert!((counts[4] as f64 / n as f64 - 0.4).abs() < 1e-3);
    // Below every inner edge the draw is always valence.
    assert_eq!(p.sample_channel(60.0, 0.999), Channel::Valence);
}

#[test]
fn channels_add_up_to_the_total_and_respect_the_edges() {
    let penn = SinglePolePenn::new(fixture_elf());
    let p = partition();
    for e in [90.0, 250.0] {
        let c = p.inverse_imfps(&penn, e).unwrap();
        let total = penn.imfp_and_stopping(e).unwrap().inverse_imfp_per_m;
        assert!(
            rel(c.total_per_m(), total) < 1e-4,
            "{e} eV: {:e} vs {total:e}",
            c.total_per_m()
        );
    }
    // A loss cannot exceed T: shells whose edge is above T do not open.
    let c = p.inverse_imfps(&penn, 120.0).unwrap();
    assert!(c.shells_per_m[0] > 0.0 && c.shells_per_m[1] > 0.0); // L3, L2
    assert_eq!(c.shells_per_m[2], 0.0); // L1, 150 eV
    assert_eq!(c.shells_per_m[3], 0.0); // K
                                        // Below the first inner edge everything is valence.
    let c = p.inverse_imfps(&penn, 60.0).unwrap();
    assert!(c.valence_per_m > 0.0);
    assert!(c.shells_per_m.iter().all(|&x| x == 0.0));
    // Above all edges all channels contribute, and deeper edges less.
    let c = p.inverse_imfps(&penn, 4000.0).unwrap();
    assert!(c.shells_per_m.iter().all(|&x| x > 0.0));
    assert!(c.shells_per_m[3] < c.shells_per_m[2]);
    assert!(p.inverse_imfps(&penn, -1.0).is_err());
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
    let ex =
        SinglePolePenn::new(fixture_elf()).with_exchange(ExchangeCorrection::new(2e4).unwrap());
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
    let ex =
        SinglePolePenn::new(fixture_elf()).with_exchange(ExchangeCorrection::new(2e4).unwrap());
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
fn channels_with_exchange_are_smaller() {
    let plain = SinglePolePenn::new(fixture_elf());
    let ex =
        SinglePolePenn::new(fixture_elf()).with_exchange(ExchangeCorrection::new(2e4).unwrap());
    let p = partition();
    let e = 400.0;
    let a = p.inverse_imfps(&plain, e).unwrap();
    let b = p.inverse_imfps(&ex, e).unwrap();
    assert!(b.valence_per_m < a.valence_per_m);
    for (x, y) in b.shells_per_m.iter().zip(&a.shells_per_m) {
        assert!(x <= y, "{x:e} vs {y:e}");
    }
    // L3, L2 and L1 are open at 400 eV and lose strength; K is closed.
    assert!(b.shells_per_m[..3]
        .iter()
        .zip(&a.shells_per_m)
        .all(|(x, y)| x < y));
    assert!(b.total_per_m() < a.total_per_m());
    // Shell i is limited to omega <= (T + E_B)/2: the K shell (1800 eV) is
    // closed at 400 eV either way.
    assert_eq!(b.shells_per_m[3], 0.0);
}

#[test]
fn results_are_bit_identical_across_thread_counts() {
    let ex =
        SinglePolePenn::new(fixture_elf()).with_exchange(ExchangeCorrection::new(2e4).unwrap());
    let p = partition();
    let energies = [40.0, 90.0, 130.0, 300.0];
    let seq: Vec<_> = energies
        .iter()
        .map(|&e| p.inverse_imfps(&ex, e).unwrap())
        .collect();
    for threads in [1usize, 2] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let par: Vec<_> = pool.install(|| {
            energies
                .par_iter()
                .map(|&e| p.inverse_imfps(&ex, e).unwrap())
                .collect()
        });
        assert_eq!(seq, par, "{threads} threads");
    }
}

#[test]
fn rejects_invalid_inputs() {
    let t = synthetic_table();
    assert!(ChannelPartition::unsourced_occupancy_weighted(&t, &[], 50.0).is_err());
    assert!(ChannelPartition::unsourced_occupancy_weighted(&t, &[(14, 1.0)], f64::NAN).is_err());
    assert!(ChannelPartition::unsourced_occupancy_weighted(&t, &[(14, 1.0)], -1.0).is_err());
    assert!(ChannelPartition::unsourced_occupancy_weighted(&t, &[(14, 0.0)], 50.0).is_err());
    assert!(
        ChannelPartition::unsourced_occupancy_weighted(&t, &[(14, 1.0), (14, 1.0)], 50.0).is_err()
    );
    assert!(matches!(
        ChannelPartition::unsourced_occupancy_weighted(&t, &[(8, 1.0)], 50.0),
        Err(ElectronDataError::Invalid { .. })
    ));
}
