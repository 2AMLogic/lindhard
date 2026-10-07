//! Insulator channels of the electron transport loop: Fröhlich LO-phonon
//! emission and absorption, and polaron trapping
//! (`lindhard::electron::phonon`).
//!
//! Every phonon and polaron parameter here is **synthetic** (round numbers
//! chosen to give useful rates), not material data. The analytic references
//! are the formulas cited in the module docs of `electron::phonon`, written
//! out again here independently of the library code.

use lindhard::constants::{BOHR_RADIUS, BOLTZMANN, ELEMENTARY_CHARGE};
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::phonon::{
    cos_theta_cdf, FrohlichPhonon, InsulatorChannels, PolaronTrapping,
};
use lindhard::electron::transport::{
    ElectronState, ElectronTally, Face, Fate, LayerTables, PhononEvent, Primary, SummaryTally,
    Transport, TransportConfig,
};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::tally::{Binning, ElectronTallyConfig, FullElectronTally};

const GRID: [f64; 5] = [0.01, 1.0, 10.0, 100.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];

/// Synthetic mode: `ħω = 0.1 eV`, `ε₀ = 4`, `ε∞ = 2`, `T = 600 K`.
const HW: f64 = 0.1;
const EPS0: f64 = 4.0;
const EPSINF: f64 = 2.0;
const TEMP: f64 = 600.0;

fn material() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

fn table(axis: SamplingAxis, rate: f64, row: impl Fn(f64, f64) -> f64) -> CrossSectionTable {
    let quantiles = GRID
        .iter()
        .map(|&e| {
            if rate == 0.0 {
                vec![]
            } else {
                PROB.iter().map(|&p| row(e, p)).collect()
            }
        })
        .collect();
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "computed in tests/electron_phonon.rs".into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: vec![rate; GRID.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

/// No elastic and no inelastic scattering: only the insulator channels act.
fn no_tables() -> LayerTables {
    LayerTables {
        elastic: table(SamplingAxis::ElasticPolarAngle, 0.0, |_, _| 0.0),
        inelastic: table(SamplingAxis::InelasticEnergyLoss, 0.0, |_, _| 0.0),
    }
}

fn mode() -> FrohlichPhonon {
    FrohlichPhonon::new(
        HW,
        EPS0,
        EPSINF,
        TEMP,
        "synthetic (tests/electron_phonon.rs)",
    )
    .unwrap()
}

fn polaron(c: f64, gamma: f64) -> PolaronTrapping {
    PolaronTrapping::new(c, gamma, "synthetic (tests/electron_phonon.rs)").unwrap()
}

fn phonon_only() -> InsulatorChannels {
    InsulatorChannels {
        phonon: Some(mode()),
        polaron: None,
    }
}

/// `n(T) = 1/(exp(ħω/k_B T) - 1)`.
fn bose() -> f64 {
    1.0 / ((HW * ELEMENTARY_CHARGE / (BOLTZMANN * TEMP)).exp() - 1.0)
}

/// The emission inverse mean free path, written out as in the cited
/// equation (no rearrangement of the logarithm).
fn emission_ref(e: f64) -> f64 {
    let x = HW / e;
    let r = (1.0 - x).sqrt();
    (1.0 / BOHR_RADIUS) * (EPS0 - EPSINF) / (EPS0 * EPSINF) * x * (bose() + 1.0) / 2.0
        * ((1.0 + r) / (1.0 - r)).ln()
}

/// The absorption inverse mean free path, as written in the cited equation.
fn absorption_ref(e: f64) -> f64 {
    let x = HW / e;
    let q = (1.0 + x).sqrt();
    (1.0 / BOHR_RADIUS) * (1.0 / EPSINF - 1.0 / EPS0) * x * bose() / 2.0
        * ((1.0 + q) / (q - 1.0)).ln()
}

// ---------------------------------------------------------------------------
// Analytic inverse mean free paths
// ---------------------------------------------------------------------------

#[test]
fn phonon_and_polaron_imfp_match_the_analytic_formulas_at_three_energies() {
    let p = mode();
    assert!((p.occupation() / bose() - 1.0).abs() < 1e-12);
    let pol = polaron(2.0e8, 0.7);
    for e in [0.25, 3.0, 40.0] {
        let (em, ab) = (
            p.emission_inverse_mfp_per_m(e),
            p.absorption_inverse_mfp_per_m(e),
        );
        assert!(
            (em / emission_ref(e) - 1.0).abs() < 1e-6,
            "emission at {e} eV"
        );
        assert!(
            (ab / absorption_ref(e) - 1.0).abs() < 1e-6,
            "absorption at {e} eV"
        );
        let tr = 2.0e8 * (-0.7 * e).exp();
        assert!(
            (pol.inverse_mfp_per_m(e) / tr - 1.0).abs() < 1e-6,
            "polaron at {e} eV"
        );
    }
}

// ---------------------------------------------------------------------------
// Detailed balance by Monte Carlo
// ---------------------------------------------------------------------------

/// The first flight of each history: its length and which phonon event ended
/// it (the run stops after one collision).
#[derive(Default)]
struct FirstEvent {
    path_m: f64,
    flights: u64,
    emissions: u64,
    absorptions: u64,
    first: bool,
}

impl ElectronTally for FirstEvent {
    fn begin_history(&mut self, _i: u64, _s: &ElectronState) {
        self.first = true;
    }
    fn step(&mut self, _from: [f64; 3], _end: &ElectronState, length_m: f64) {
        if self.first {
            self.path_m += length_m;
            self.flights += 1;
        }
    }
    fn phonon(&mut self, _a: &ElectronState, kind: PhononEvent, _hw: f64, _t: f64) {
        if self.first {
            match kind {
                PhononEvent::Emission => self.emissions += 1,
                PhononEvent::Absorption => self.absorptions += 1,
            }
            self.first = false;
        }
    }
    fn merge(&mut self, o: Self) {
        self.path_m += o.path_m;
        self.flights += o.flights;
        self.emissions += o.emissions;
        self.absorptions += o.absorptions;
    }
}

fn first_events(energy_ev: f64, seed: u64, n: u64) -> FirstEvent {
    let mut cfg = TransportConfig::new(0.001);
    cfg.max_events = 1;
    let t = Transport::new(Stack::semi_infinite(material()), vec![no_tables()], cfg)
        .unwrap()
        .with_insulator_channels(0, phonon_only())
        .unwrap();
    t.run(
        seed,
        n,
        4096,
        &Primary::normal(energy_ev),
        FirstEvent::default,
    )
    .unwrap()
    .tally
}

#[test]
fn emission_over_absorption_follows_detailed_balance() {
    // Emission E -> E' and absorption E' -> E with E' = E - ħω: the measured
    // rates obey E λ₊⁻¹(E) / (E' λ₋⁻¹(E')) = (n + 1)/n.
    let (e, ep) = (1.0, 1.0 - HW);
    let n = 300_000;
    let hi = first_events(e, 0xD1, n);
    let lo = first_events(ep, 0xD2, n);
    assert_eq!(hi.flights, n);
    assert_eq!(lo.flights, n);
    let em_rate = hi.emissions as f64 / hi.path_m;
    let ab_rate = lo.absorptions as f64 / lo.path_m;
    // Each measured rate also matches its formula.
    let sig_em = (1.0 / hi.emissions as f64 + 1.0 / n as f64).sqrt();
    let sig_ab = (1.0 / lo.absorptions as f64 + 1.0 / n as f64).sqrt();
    assert!(
        (em_rate / emission_ref(e) - 1.0).abs() < 4.0 * sig_em,
        "emission rate {em_rate} vs {}",
        emission_ref(e)
    );
    assert!(
        (ab_rate / absorption_ref(ep) - 1.0).abs() < 4.0 * sig_ab,
        "absorption rate {ab_rate} vs {}",
        absorption_ref(ep)
    );
    let ratio = e * em_rate / (ep * ab_rate);
    let expect = (bose() + 1.0) / bose();
    let sigma = (sig_em * sig_em + sig_ab * sig_ab).sqrt();
    assert!(
        (ratio / expect - 1.0).abs() < 4.0 * sigma,
        "detailed balance: measured {ratio}, (n+1)/n = {expect}, 1 sigma {sigma}"
    );
}

#[test]
fn no_absorption_at_zero_temperature() {
    let cold = FrohlichPhonon::new(HW, EPS0, EPSINF, 0.0, "synthetic").unwrap();
    let mut cfg = TransportConfig::new(0.05);
    cfg.max_events = 1000;
    let t = Transport::new(Stack::semi_infinite(material()), vec![no_tables()], cfg)
        .unwrap()
        .with_insulator_channels(
            0,
            InsulatorChannels {
                phonon: Some(cold),
                polaron: None,
            },
        )
        .unwrap();
    let r = t
        .run(5, 2000, 64, &Primary::normal(1.0), SummaryTally::default)
        .unwrap()
        .tally;
    assert_eq!(r.phonon_absorptions, 0);
    assert!(r.phonon_emissions > 0);
}

// ---------------------------------------------------------------------------
// Angular distribution
// ---------------------------------------------------------------------------

#[derive(Default)]
struct FirstCos {
    emission_cos: Vec<f64>,
    first: bool,
}

impl ElectronTally for FirstCos {
    fn begin_history(&mut self, _i: u64, _s: &ElectronState) {
        self.first = true;
    }
    fn phonon(&mut self, a: &ElectronState, kind: PhononEvent, _hw: f64, theta: f64) {
        if self.first && kind == PhononEvent::Emission {
            // Normal incidence: the new x direction is cos θ.
            assert!((a.dir[0] - theta.cos()).abs() < 1e-12);
            self.emission_cos.push(a.dir[0]);
        }
        self.first = false;
    }
    fn merge(&mut self, o: Self) {
        self.emission_cos.extend(o.emission_cos);
    }
}

#[test]
fn phonon_deflection_follows_the_frohlich_angular_distribution() {
    let e = 0.5;
    let mut cfg = TransportConfig::new(0.001);
    cfg.max_events = 1;
    let t = Transport::new(Stack::semi_infinite(material()), vec![no_tables()], cfg)
        .unwrap()
        .with_insulator_channels(0, phonon_only())
        .unwrap();
    let r = t
        .run(0xA9, 100_000, 4096, &Primary::normal(e), FirstCos::default)
        .unwrap()
        .tally;
    let n = r.emission_cos.len() as f64;
    assert!(n > 50_000.0);
    for mu in [-0.5, 0.0, 0.5, 0.9] {
        let p = cos_theta_cdf(e, e - HW, mu);
        let got = r.emission_cos.iter().filter(|&&c| c <= mu).count() as f64 / n;
        let sigma = (p * (1.0 - p) / n).sqrt();
        assert!((got - p).abs() < 4.0 * sigma, "F({mu}): {got} vs {p}");
    }
}

// ---------------------------------------------------------------------------
// Polaron trapping
// ---------------------------------------------------------------------------

#[derive(Default)]
struct TrapLog {
    trapped: u64,
    escaped_back: u64,
    other: u64,
    deposited_ev: f64,
    max_depth_m: f64,
}

impl ElectronTally for TrapLog {
    fn polaron_trapped(&mut self, at: &ElectronState) {
        self.deposited_ev += at.energy_ev;
        self.max_depth_m = self.max_depth_m.max(at.pos[0]);
    }
    fn end_history(&mut self, _i: u64, fate: Fate) {
        match fate {
            Fate::PolaronTrapped => self.trapped += 1,
            Fate::Escaped(Face::Back) => self.escaped_back += 1,
            _ => self.other += 1,
        }
    }
    fn merge(&mut self, o: Self) {
        self.trapped += o.trapped;
        self.escaped_back += o.escaped_back;
        self.other += o.other;
        self.deposited_ev += o.deposited_ev;
        self.max_depth_m = self.max_depth_m.max(o.max_depth_m);
    }
}

#[test]
fn trapping_probability_over_a_fixed_path_is_one_minus_exp() {
    let (c, gamma, e, s) = (1.0e8, 0.5, 2.0, 2.0e-8);
    let t = Transport::new(
        Stack::new(vec![(material(), s)], None).unwrap(),
        vec![no_tables()],
        TransportConfig::new(0.01),
    )
    .unwrap()
    .with_insulator_channels(
        0,
        InsulatorChannels {
            phonon: None,
            polaron: Some(polaron(c, gamma)),
        },
    )
    .unwrap();
    let n = 100_000u64;
    let r = t
        .run(0x7A, n, 1000, &Primary::normal(e), TrapLog::default)
        .unwrap()
        .tally;
    assert_eq!(r.other, 0);
    assert_eq!(r.trapped + r.escaped_back, n);
    let p = 1.0 - (-s * c * (-gamma * e).exp()).exp();
    let got = r.trapped as f64 / n as f64;
    let sigma = (p * (1.0 - p) / n as f64).sqrt();
    assert!((got - p).abs() < 4.0 * sigma, "trapped {got} vs {p}");
    // Each trapped electron deposits its whole energy, inside the slab.
    assert!((r.deposited_ev - e * r.trapped as f64).abs() < 1e-9 * r.deposited_ev);
    assert!(r.max_depth_m <= s);
}

#[test]
fn a_trapped_electron_ends_its_history_and_is_tallied() {
    // A trap so likely that every electron is trapped on its first collision.
    let t = Transport::new(
        Stack::semi_infinite(material()),
        vec![no_tables()],
        TransportConfig::new(0.01),
    )
    .unwrap()
    .with_insulator_channels(
        0,
        InsulatorChannels {
            phonon: None,
            polaron: Some(polaron(1e9, 0.0)),
        },
    )
    .unwrap();
    let r = t
        .run(3, 500, 64, &Primary::normal(7.0), SummaryTally::default)
        .unwrap()
        .tally;
    assert_eq!(r.histories, 500);
    assert_eq!(r.polaron_trapped, 500);
    assert_eq!(r.stopped + r.trapped + r.event_capped, 0);
    assert!((r.polaron_deposited_ev - 3500.0).abs() < 1e-9);
    assert_eq!(r.rest_energy_ev, 0.0);
}

// ---------------------------------------------------------------------------
// Opt-in, metadata, energy and determinism
// ---------------------------------------------------------------------------

fn mixed(channels: bool) -> Transport {
    let layer = |el: f64, inel: f64, frac: f64| LayerTables {
        elastic: table(SamplingAxis::ElasticPolarAngle, 1.0 / el, |_, p| {
            (1.0 - 2.0 * p).acos()
        }),
        inelastic: table(
            SamplingAxis::InelasticEnergyLoss,
            1.0 / inel,
            move |e, p| frac * e * p,
        ),
    };
    let t = Transport::new(
        Stack::new(vec![(material(), 5e-9)], Some(material())).unwrap(),
        vec![layer(2e-9, 4e-9, 0.2), layer(1e-9, 3e-9, 0.3)],
        TransportConfig::new(0.3),
    )
    .unwrap();
    if !channels {
        return t;
    }
    // The first layer is the "insulator"; the substrate stays without
    // insulator channels.
    t.with_insulator_channels(
        0,
        InsulatorChannels {
            phonon: Some(mode()),
            polaron: Some(polaron(3e8, 0.4)),
        },
    )
    .unwrap()
}

fn mixed_run(threads: usize, channels: bool) -> SummaryTally {
    let t = mixed(channels);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let primary = Primary {
        energy_ev: 60.0,
        direction: [1.0, 0.2, 0.0],
    };
    pool.install(|| {
        t.run(0xC0DE, 800, 9, &primary, SummaryTally::default)
            .unwrap()
            .tally
    })
}

#[test]
fn insulator_channels_are_bit_identical_across_thread_counts() {
    let a = mixed_run(1, true);
    assert!(a.phonon_emissions > 0 && a.phonon_absorptions > 0 && a.polaron_trapped > 0);
    for n in [2, 8] {
        let b = mixed_run(n, true);
        assert_eq!(a.path_m.to_bits(), b.path_m.to_bits());
        assert_eq!(a.phonon_emitted_ev.to_bits(), b.phonon_emitted_ev.to_bits());
        assert_eq!(
            a.polaron_deposited_ev.to_bits(),
            b.polaron_deposited_ev.to_bits()
        );
        assert_eq!(a, b, "tally differs on {n} threads");
    }
}

#[test]
fn energy_is_conserved_with_insulator_channels() {
    let s = mixed_run(2, true);
    assert_eq!(s.event_capped, 0);
    let incident = 800.0 * 60.0;
    let out = s.inelastic_loss_ev
        + s.phonon_emitted_ev
        + s.escaped_energy_ev
        + s.rest_energy_ev
        + s.polaron_deposited_ev;
    let inp = incident + s.phonon_absorbed_ev;
    assert!((out - inp).abs() < 1e-9 * inp, "{out} vs {inp}");
}

#[test]
fn full_tally_balances_energy_with_phonon_and_polaron_channels() {
    let t = mixed(true);
    let config = ElectronTallyConfig::new(
        Binning::new(0.0, 60.0, 30).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    let proto = FullElectronTally::new(&t, config).unwrap();
    let primary = Primary {
        energy_ev: 60.0,
        direction: [1.0, 0.2, 0.0],
    };
    let r = t
        .run(0xC0DE, 800, 9, &primary, || proto.clone())
        .unwrap()
        .tally
        .report();
    let b = r.budget;
    assert!(r.fates.polaron_trapped > 0 && b.polaron_ev > 0.0);
    assert!(b.phonon_emitted_ev > 0.0 && b.phonon_absorbed_ev > 0.0);
    assert_eq!(r.fates.event_capped, 0);
    assert!(b.relative_imbalance < 1e-9, "{}", b.relative_imbalance);
    // The deposit includes the phonon and polaron parts.
    assert!(b.deposited_ev >= b.phonon_emitted_ev + b.polaron_ev);
    let layers: f64 = r.deposition.per_layer_ev.iter().sum();
    assert!((layers - b.deposited_ev).abs() < 1e-9 * b.deposited_ev);
    // The same run with the summary tally agrees on the polaron deposit.
    let s = mixed_run(2, true);
    assert_eq!(r.fates.polaron_trapped, s.polaron_trapped);
    assert!((b.polaron_ev - s.polaron_deposited_ev).abs() < 1e-9 * b.polaron_ev);
}

#[test]
fn channels_are_off_by_default_and_the_choice_is_in_the_metadata() {
    let off = mixed(false);
    assert!(off.insulator_channels().iter().all(|c| !c.any()));
    // Explicitly switching nothing on changes nothing, bit for bit.
    let none = mixed(false)
        .with_insulator_channels(0, InsulatorChannels::none())
        .unwrap();
    let primary = Primary::normal(60.0);
    let a = off
        .run(1, 300, 16, &primary, SummaryTally::default)
        .unwrap();
    let b = none
        .run(1, 300, 16, &primary, SummaryTally::default)
        .unwrap();
    assert_eq!(a.tally, b.tally);
    assert_eq!(a.tally.phonon_emissions + a.tally.polaron_trapped, 0);
    assert!(a
        .metadata
        .layers
        .iter()
        .all(|l| l.phonon.is_none() && l.polaron.is_none()));

    let on = mixed(true)
        .run(1, 10, 16, &primary, SummaryTally::default)
        .unwrap();
    let l0 = &on.metadata.layers[0];
    assert_eq!(l0.phonon.as_ref().unwrap().hbar_omega_ev(), HW);
    assert_eq!(l0.phonon.as_ref().unwrap().temperature_k(), TEMP);
    assert_eq!(l0.polaron.as_ref().unwrap().c_per_m(), 3e8);
    assert!(on.metadata.layers[1].phonon.is_none() && on.metadata.layers[1].polaron.is_none());
    let json = serde_json::to_string(&on.metadata).unwrap();
    assert!(json.contains("\"hbar_omega_ev\":0.1"), "{json}");
    assert!(json.contains("\"gamma_per_ev\":0.4"), "{json}");
    assert!(json.contains("\"polaron\":null"), "{json}");
}

#[test]
fn opting_in_a_layer_out_of_range_is_rejected() {
    assert!(mixed(false)
        .with_insulator_channels(2, phonon_only())
        .is_err());
}
