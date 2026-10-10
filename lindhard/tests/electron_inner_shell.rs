//! Inner-shell ionisation channels in the electron loop (#273), on synthetic
//! cross-section tables and synthetic band parameters (made-up numbers chosen
//! to exercise the code, not material data): valence-only compatibility,
//! channel frequencies and conditional losses, shell edges and rows that
//! straddle them, per-event and per-run energy closure with and without the
//! secondary model, determinism across thread counts, and the rejection of
//! incompatible channel tables.

use lindhard::electron::boundary::{BandModel, BandStructure};
use lindhard::electron::data::{
    CrossSectionTable, CrossSectionTableParts, SamplingAxis, ShellChannelTable, Subshell,
};
use lindhard::electron::inelastic::InnerShell;
use lindhard::electron::secondary::{SecondaryEvent, SecondaryModel};
use lindhard::electron::transport::{
    BoundaryModel, CutoffReference, ElectronState, ElectronTally, LayerTables, Primary,
    SummaryTally, Transport, TransportConfig,
};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::tally::{Binning, ElectronTallyConfig, FullElectronTally};

const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];
const PROVENANCE: &str = "synthetic, tests/electron_inner_shell.rs";
const MATERIAL: &str = "synthetic";

fn material() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

/// A table on `GRID` with the given rate per row and `row(e, p)` as the
/// inverse CDF of every row with a positive rate.
fn table(
    axis: SamplingAxis,
    material: &str,
    rates: [f64; 5],
    row: impl Fn(f64, f64) -> f64,
) -> CrossSectionTable {
    let quantiles = GRID
        .iter()
        .zip(rates)
        .map(|(&e, r)| {
            if r == 0.0 {
                vec![]
            } else {
                PROB.iter().map(|&p| row(e, p)).collect()
            }
        })
        .collect();
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: material.into(),
        provenance: PROVENANCE.into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: rates.to_vec(),
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

fn elastic_none() -> CrossSectionTable {
    table(
        SamplingAxis::ElasticPolarAngle,
        MATERIAL,
        [0.0; 5],
        |_, _| 0.0,
    )
}

fn elastic_isotropic(lambda_m: f64) -> CrossSectionTable {
    table(
        SamplingAxis::ElasticPolarAngle,
        MATERIAL,
        [1.0 / lambda_m; 5],
        |_, p| (1.0 - 2.0 * p).acos(),
    )
}

/// Valence losses uniform on `[0, min(10 eV, E)]`, at a constant rate.
fn valence(rate: f64) -> CrossSectionTable {
    table(
        SamplingAxis::InelasticEnergyLoss,
        MATERIAL,
        [rate; 5],
        |e, p| (10.0 * p).min(e),
    )
}

/// Valence losses uniform on `[0, frac E]`.
fn valence_frac(lambda_m: f64, frac: f64) -> CrossSectionTable {
    table(
        SamplingAxis::InelasticEnergyLoss,
        MATERIAL,
        [1.0 / lambda_m; 5],
        move |e, p| frac * e * p,
    )
}

/// A shell channel of binding energy `b` whose rows with a positive rate
/// draw `ω` uniformly on `[lo, hi]`, cut to the row energy (rows must start
/// at or above `b`).
fn shell(label: &str, b: f64, rates: [f64; 5], lo: f64, hi: f64) -> ShellChannelTable {
    shell_for(MATERIAL, label, b, rates, lo, hi)
}

fn shell_for(
    material: &str,
    label: &str,
    b: f64,
    rates: [f64; 5],
    lo: f64,
    hi: f64,
) -> ShellChannelTable {
    let t = table(
        SamplingAxis::InelasticEnergyLoss,
        material,
        rates,
        |e, p| (lo + (hi - lo) * p).min(e),
    );
    ShellChannelTable::new(
        14,
        Subshell::from_label(label).unwrap(),
        b,
        "synthetic binding energy, not physical data",
        t,
    )
    .unwrap()
}

/// Synthetic metal: `E_F = 5 eV`, `Φ = 4 eV`, so `U = 9 eV`.
fn metal() -> BandStructure {
    BandStructure::new(
        BandModel::Metal {
            fermi_ev: 5.0,
            work_function_ev: 4.0,
        },
        PROVENANCE,
    )
    .unwrap()
}

/// Synthetic insulator: `E_F = 7.5 eV`, `E_g = 3 eV`, `U = 10 eV`.
fn insulator() -> BandStructure {
    BandStructure::new(
        BandModel::Insulator {
            valence_band_width_ev: 6.0,
            band_gap_ev: 3.0,
            affinity_ev: 1.0,
        },
        PROVENANCE,
    )
    .unwrap()
}

fn semi_infinite() -> Stack {
    Stack::new(vec![], Some(material())).unwrap()
}

fn off(cutoff: f64) -> TransportConfig {
    TransportConfig::new(cutoff)
}

/// Secondaries and the step barrier, cutoff 0.5 eV above the vacuum level
/// (thresholds 9.5 eV in the metal, 10.5 eV in the insulator).
fn full_physics() -> TransportConfig {
    let mut cfg = TransportConfig::new(0.5);
    cfg.secondaries = SecondaryModel::KIEFT_BOSCH;
    cfg.boundary = BoundaryModel::STEP_BARRIER;
    cfg.cutoff_reference = CutoffReference::VacuumLevel;
    cfg
}

/// A 6 nm metal on an insulator substrate, valence-only.
fn metal_on_insulator() -> Transport {
    Transport::with_band_structures(
        Stack::new(vec![(material(), 6e-9)], Some(material())).unwrap(),
        vec![
            LayerTables {
                elastic: elastic_isotropic(3e-9),
                inelastic: valence_frac(4e-9, 0.3),
            },
            LayerTables {
                elastic: elastic_isotropic(2e-9),
                inelastic: valence_frac(3e-9, 0.25),
            },
        ],
        vec![metal(), insulator()],
        full_physics(),
    )
    .unwrap()
}

/// The shells of [`metal_on_insulator`]: in the metal an L3-like shell
/// (100 eV, losses `[100, 110]`, but exactly 100 at the 100 eV row, where it
/// is closed) and an L1-like shell whose every loss equals its binding
/// energy (60 eV); in the insulator a shell of 150 eV with rows only from
/// 1 keV, so the rate between the 100 eV and 1 keV rows interpolates from an
/// inactive row.
fn with_shells(t: Transport) -> Transport {
    let r = 1.0 / 5e-9;
    t.with_inner_shells(
        0,
        vec![
            shell("L3", 100.0, [0.0, 0.0, r, r, r], 100.0, 110.0),
            shell("L1", 60.0, [0.0, 0.0, r, r, r], 60.0, 60.0),
        ],
    )
    .unwrap()
    .with_inner_shells(
        1,
        vec![shell("K", 150.0, [0.0, 0.0, 0.0, r, r], 150.0, 200.0)],
    )
    .unwrap()
}

// ---------------------------------------------------------------------------
// Valence-only compatibility
// ---------------------------------------------------------------------------

#[test]
fn no_shells_and_inactive_shells_leave_runs_bit_identical() {
    let run = |t: &Transport| {
        t.run(5, 200, 16, &Primary::normal(400.0), SummaryTally::default)
            .unwrap()
            .tally
    };
    let base = metal_on_insulator();
    let s0 = run(&base);
    assert!(s0.inelastic_events > 1_000 && s0.secondaries > 100);
    assert_eq!(s0.inner_shell_events, 0);
    // An empty list, and a channel whose rate is zero everywhere.
    assert_eq!(run(&base.clone().with_inner_shells(0, vec![]).unwrap()), s0);
    let inactive = base
        .clone()
        .with_inner_shells(0, vec![shell("K", 50.0, [0.0; 5], 50.0, 50.0)])
        .unwrap()
        .with_inner_shells(1, vec![shell("K", 50.0, [0.0; 5], 50.0, 50.0)])
        .unwrap();
    assert_eq!(run(&inactive), s0);
    // The metadata of a valence-only run serializes as before.
    let m = base.metadata(5, 200, 16, &Primary::normal(400.0));
    let json = serde_json::to_string(&m).unwrap();
    assert!(!json.contains("inner_shells"), "{json}");
    let m = inactive.metadata(5, 200, 16, &Primary::normal(400.0));
    assert_eq!(m.layers[0].inner_shells.len(), 1);
    assert_eq!(m.layers[0].inner_shells[0].subshell, "K");
    assert_eq!(
        m.layers[0].inner_shells[0].binding_provenance,
        "synthetic binding energy, not physical data"
    );
    assert!(serde_json::to_string(&m).unwrap().contains("inner_shells"));
}

// ---------------------------------------------------------------------------
// Channel frequencies and conditional losses
// ---------------------------------------------------------------------------

/// Every inelastic event: the energy before, the loss, and the binding
/// energy of its shell (`None` for valence).
#[derive(Default)]
struct Events {
    before: f64,
    list: Vec<(f64, f64, Option<f64>)>,
}

impl ElectronTally for Events {
    fn step(&mut self, _from: [f64; 3], end: &ElectronState, _l: f64) {
        self.before = end.energy_ev;
    }
    fn inelastic(&mut self, after: &ElectronState, w: f64) {
        assert_eq!(after.energy_ev, self.before - w);
        self.list.push((self.before, w, None));
    }
    fn inner_shell(&mut self, _after: &ElectronState, shell: &InnerShell, w: f64) {
        let last = self.list.last_mut().expect("inelastic comes first");
        assert_eq!(last.1, w);
        assert!(last.2.is_none());
        last.2 = Some(shell.binding_energy_ev);
    }
    fn merge(&mut self, o: Self) {
        self.list.extend(o.list);
    }
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

#[test]
fn channels_are_chosen_by_rate_and_losses_come_from_the_chosen_channel() {
    // Constant rates over the energies visited (the 1 and 10 keV rows are
    // equal), unequal and with disjoint loss ranges; one channel inactive.
    let (rv, ra, rb) = (1e9, 5e8, 2.5e8);
    let t = Transport::new(
        semi_infinite(),
        vec![LayerTables {
            elastic: elastic_none(),
            inelastic: valence(rv),
        }],
        off(2_000.0),
    )
    .unwrap()
    .with_inner_shells(
        0,
        vec![
            shell("L3", 100.0, [0.0, 0.0, 0.0, ra, ra], 100.0, 110.0),
            shell("L1", 50.0, [0.0; 5], 50.0, 50.0),
            shell("K", 300.0, [0.0, 0.0, 0.0, rb, rb], 300.0, 320.0),
        ],
    )
    .unwrap();
    let ev = t
        .run(17, 300, 16, &Primary::normal(9_000.0), Events::default)
        .unwrap()
        .tally
        .list;
    let n = ev.len() as f64;
    assert!(n > 20_000.0, "{n}");
    let by = |b: Option<f64>| -> Vec<f64> { ev.iter().filter(|e| e.2 == b).map(|e| e.1).collect() };
    let (v, a, k) = (by(None), by(Some(100.0)), by(Some(300.0)));
    assert!(by(Some(50.0)).is_empty(), "an inactive channel was chosen");
    assert_eq!(v.len() + a.len() + k.len(), ev.len());
    let total = rv + ra + rb;
    for (got, rate) in [(v.len(), rv), (a.len(), ra), (k.len(), rb)] {
        let p = rate / total;
        let f = got as f64 / n;
        let sigma = (p * (1.0 - p) / n).sqrt();
        assert!((f - p).abs() < 5.0 * sigma, "fraction {f} vs {p}");
    }
    // Conditional losses: each channel's own uniform range.
    assert!(v.iter().all(|&w| (0.0..=10.0).contains(&w)));
    assert!(a.iter().all(|&w| (100.0..=110.0).contains(&w)));
    assert!(k.iter().all(|&w| (300.0..=320.0).contains(&w)));
    // Uniform on a range of width d: the mean's standard error is
    // d / sqrt(12 n).
    for (w, centre, d) in [(&v, 5.0, 10.0), (&a, 105.0, 10.0), (&k, 310.0, 20.0)] {
        let se = d / (12.0 * w.len() as f64).sqrt();
        assert!(
            (mean(w) - centre).abs() < 5.0 * se,
            "{} vs {centre}",
            mean(w)
        );
    }
}

#[test]
fn a_shell_opens_only_above_its_edge_across_an_inactive_row() {
    // The shell's 1 keV row is inactive and its 10 keV row active, so between
    // them the interpolated rate is positive at every energy, also below the
    // 4 keV edge: the transport must close the channel there itself, and
    // keep every loss within [B, E].
    let b = 4_000.0;
    let t = Transport::new(
        semi_infinite(),
        vec![LayerTables {
            elastic: elastic_none(),
            inelastic: valence(1e9),
        }],
        off(50.0),
    )
    .unwrap()
    .with_inner_shells(
        0,
        vec![shell("K", b, [0.0, 0.0, 0.0, 0.0, 1e9], b, 4_100.0)],
    )
    .unwrap();
    let ev = t
        .run(23, 200, 16, &Primary::normal(9_000.0), Events::default)
        .unwrap()
        .tally
        .list;
    let shells: Vec<_> = ev.iter().filter(|e| e.2.is_some()).collect();
    assert!(shells.len() > 50, "{}", shells.len());
    for &&(before, w, _) in &shells {
        assert!(before > b, "shell event at {before} eV, below its edge");
        assert!(
            w >= b && w <= before.min(4_100.0),
            "loss {w} at {before} eV"
        );
    }
    // Some electrons do pass the edge with valence losses only.
    assert!(ev.iter().any(|e| e.2.is_none() && e.0 <= b));
}

// ---------------------------------------------------------------------------
// Energy closure
// ---------------------------------------------------------------------------

/// Checks every inelastic event of the Kieft-Bosch model, valence or shell.
struct EventCheck {
    bands: Vec<BandStructure>,
    before: f64,
    w: f64,
    shell: Option<f64>,
    shell_events: u64,
    tracked: u64,
    untracked: u64,
    at_edge: u64,
}

impl EventCheck {
    fn new(bands: Vec<BandStructure>) -> Self {
        Self {
            bands,
            before: 0.0,
            w: 0.0,
            shell: None,
            shell_events: 0,
            tracked: 0,
            untracked: 0,
            at_edge: 0,
        }
    }
}

impl ElectronTally for EventCheck {
    fn step(&mut self, _from: [f64; 3], end: &ElectronState, _l: f64) {
        self.before = end.energy_ev;
    }
    fn inelastic(&mut self, after: &ElectronState, w: f64) {
        assert_eq!(after.energy_ev, self.before - w);
        self.w = w;
        self.shell = None;
    }
    fn inner_shell(&mut self, after: &ElectronState, shell: &InnerShell, w: f64) {
        let ef = self.bands[after.layer].fermi_ev();
        let b = shell.binding_energy_ev;
        assert_eq!(w, self.w);
        // The channel's kinematic limits on the band-bottom axis.
        assert!(w >= b && w <= self.before - ef, "{w} at {}", self.before);
        assert!(after.energy_ev >= ef - 1e-12 * self.before);
        self.shell = Some(b);
    }
    fn secondary(&mut self, p: &ElectronState, e: &SecondaryEvent, c: Option<&ElectronState>) {
        let ef = self.bands[p.layer].fermi_ev();
        assert_eq!(e.loss_ev, self.w);
        let scale = ef + e.loss_ev;
        let residual = (e.loss_ev - (e.secondary_ev + e.binding_ev + e.deposited_ev)).abs();
        assert!(residual <= 4.0 * f64::EPSILON * scale, "{e:?}");
        let Some(b) = self.shell.take() else {
            return;
        };
        self.shell_events += 1;
        assert!(e.liberated);
        assert_eq!(e.binding_ev, b - ef);
        // Verduin Eq. 3.86 with the shell's B: E_SE = E_F + ω - B >= E_F.
        let e_se = ef + e.loss_ev - b;
        assert!(e_se >= ef - 1e-12 * scale);
        if e.loss_ev == b {
            self.at_edge += 1;
        }
        match c {
            Some(s) => {
                self.tracked += 1;
                assert_eq!(
                    (e.secondary_ev, s.energy_ev, e.deposited_ev),
                    (e_se, e_se, 0.0)
                );
            }
            None => {
                self.untracked += 1;
                assert_eq!((e.secondary_ev, e.deposited_ev), (0.0, e_se));
            }
        }
    }
    fn merge(&mut self, o: Self) {
        self.shell_events += o.shell_events;
        self.tracked += o.tracked;
        self.untracked += o.untracked;
        self.at_edge += o.at_edge;
    }
}

#[test]
fn shell_events_conserve_energy_per_event_with_the_secondary_model() {
    let t = with_shells(metal_on_insulator());
    let r = t
        .run(11, 400, 16, &Primary::normal(1_500.0), || {
            EventCheck::new(vec![metal(), insulator()])
        })
        .unwrap()
        .tally;
    assert!(r.shell_events > 1_000, "{}", r.shell_events);
    // Some shell secondaries are followed, some are below the threshold
    // (every one of the 60 eV shell, whose loss equals its binding energy).
    assert!(
        r.tracked > 100 && r.untracked > 100,
        "{} / {}",
        r.tracked,
        r.untracked
    );
    assert!(r.at_edge > 100, "{}", r.at_edge);
}

fn summary_balance(s: &SummaryTally, n: u64, e0: f64) {
    let source = n as f64 * e0 + s.secondary_energy_ev;
    let sink = s.inelastic_loss_ev + s.escaped_energy_ev + s.rest_energy_ev + s.barrier_ev;
    assert!(
        (source - sink).abs() < 1e-9 * source,
        "source {source} vs sink {sink}"
    );
    assert!(s.inner_shell_events > 0 && s.inner_shell_loss_ev > 0.0);
    assert!(s.inner_shell_loss_ev < s.inelastic_loss_ev);
}

#[test]
fn runs_with_shells_close_the_energy_balance_with_and_without_secondaries() {
    let (n, e0) = (300u64, 1_500.0);
    // Kieft-Bosch, step barrier, tracked and untracked shell secondaries.
    let t = with_shells(metal_on_insulator());
    let s = t
        .run(3, n, 32, &Primary::normal(e0), SummaryTally::default)
        .unwrap()
        .tally;
    summary_balance(&s, n, e0);
    let events = s.secondary_energy_ev + s.binding_ev + s.deposited_ev;
    assert!((events - s.inelastic_loss_ev).abs() < 1e-9 * s.inelastic_loss_ev);
    // The full tally's budget, with the binding energies kept in the solid.
    let mut cfg = ElectronTallyConfig::new(
        Binning::new(0.0, e0, 30).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    cfg.cylindrical = None;
    let proto = FullElectronTally::new(&t, cfg).unwrap();
    let report = t
        .run(3, n, 32, &Primary::normal(e0), || proto.clone())
        .unwrap()
        .tally
        .report();
    assert!(
        report.budget.relative_imbalance < 1e-9,
        "{:?}",
        report.budget
    );
    // Secondaries off, no band parameters: the loss is all deposited.
    let off_t = Transport::new(
        semi_infinite(),
        vec![LayerTables {
            elastic: elastic_isotropic(3e-9),
            inelastic: valence_frac(4e-9, 0.3),
        }],
        off(5.0),
    )
    .unwrap()
    .with_inner_shells(
        0,
        vec![shell("L3", 100.0, [0.0, 0.0, 1e8, 1e8, 1e8], 100.0, 110.0)],
    )
    .unwrap();
    let s = off_t
        .run(4, n, 32, &Primary::normal(e0), SummaryTally::default)
        .unwrap()
        .tally;
    assert_eq!((s.secondaries, s.secondary_energy_ev), (0, 0.0));
    summary_balance(&s, n, e0);
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

#[test]
fn shell_runs_are_bit_identical_across_thread_counts() {
    let t = with_shells(metal_on_insulator());
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                t.run(9, 120, 8, &Primary::normal(1_200.0), SummaryTally::default)
                    .unwrap()
                    .tally
            })
    };
    let one = run(1);
    assert!(one.inner_shell_events > 0);
    assert_eq!(one, run(2));
    assert_eq!(one, run(4));
}

// ---------------------------------------------------------------------------
// Rejection of incompatible channels
// ---------------------------------------------------------------------------

#[test]
fn incompatible_shell_channels_are_rejected() {
    let t = metal_on_insulator();
    let r = [0.0, 0.0, 1e8, 1e8, 1e8];
    // Out-of-range layer.
    assert!(t
        .clone()
        .with_inner_shells(2, vec![shell("K", 100.0, r, 100.0, 110.0)])
        .is_err());
    // A table built for another target.
    assert!(t
        .clone()
        .with_inner_shells(0, vec![shell_for("other", "K", 100.0, r, 100.0, 110.0)])
        .is_err());
    // The same shell twice.
    assert!(t
        .clone()
        .with_inner_shells(
            0,
            vec![
                shell("K", 100.0, r, 100.0, 110.0),
                shell("K", 100.0, r, 100.0, 110.0)
            ]
        )
        .is_err());
    // A binding energy at or below the Fermi level (5 eV in the metal).
    let rows = [0.0, 1e8, 1e8, 1e8, 1e8];
    assert!(t
        .clone()
        .with_inner_shells(0, vec![shell("M1", 5.0, rows, 5.0, 6.0)])
        .is_err());
    assert!(t
        .with_inner_shells(0, vec![shell("M1", 5.5, rows, 5.5, 6.0)])
        .is_ok());
}
