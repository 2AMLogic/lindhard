//! Integration tests for the amorphous BCA engine (`lindhard::ion::bca`).

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaError, BcaTally, Beam, ElectronicChannel, ElectronicLoss, Face,
    LatticeDeposit, MeanFreePath, Particle, SummaryTally,
};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::{ElectronicStopping, Ion, StoppingError, ValidityRange};
use lindhard::material::{EnergyKind, Material};

const NM: f64 = 1e-9;

/// One ZBL table shared by every test in this binary (the angle depends on
/// the screening function only, so one table serves every pair).
fn table() -> &'static ScatteringTable {
    static T: OnceLock<ScatteringTable> = OnceLock::new();
    T.get_or_init(|| {
        ScatteringTable::build(
            &Potential::new(Screening::ZblUniversal, 14.0, 14.0),
            &TableSpec {
                eps_min: 1e-6,
                eps_max: 1e4,
                beta_min: 1e-5,
                beta_max: 1e2,
                per_decade: 16,
            },
        )
    })
}

/// Set `E_d` = 15 eV for every element and `E_s` = 2 eV where the element
/// table has no default. These are test parameters, not data: the engine only
/// needs some value, and the tests do not compare against measurements.
fn ready(mut m: Material) -> Material {
    for (z, _) in m.atom_fractions() {
        m.set_displacement_energy_ev(z, 15.0).unwrap();
        if m.surface_binding_energy_ev(z).is_err() {
            m.set_surface_binding_energy_ev(z, 2.0).unwrap();
        }
    }
    m
}

fn elemental(z: u8) -> Material {
    ready(Material::from_atom_fractions(&[(z, 1.0)], None).unwrap())
}

fn si() -> Material {
    elemental(14)
}

fn sio2() -> Material {
    ready(Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap())
}

fn run_summary(bca: &Bca, bin_m: f64, bins: usize) -> SummaryTally {
    bca.run(|| SummaryTally::new(bin_m, bins)).unwrap()
}

fn with_threads<R: Send>(n: usize, f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .unwrap()
        .install(f)
}

#[test]
fn energy_is_conserved_per_history_with_cascades() {
    let ls = LindhardScharff::new();
    // Every code path: two finite layers of different materials, front and
    // back faces, cascades, both free-path conventions, both electronic modes,
    // with and without weak collisions.
    let stack = Stack::new(vec![(sio2(), 3.0 * NM), (si(), 4.0 * NM)], None).unwrap();
    for (mfp, el, weak) in [
        (MeanFreePath::Constant, ElectronicLoss::NonLocal, 0),
        (
            MeanFreePath::EnergyDependent {
                min_cm_angle_rad: 0.01,
            },
            ElectronicLoss::EquipartitionLsOr,
            0,
        ),
        (MeanFreePath::Constant, ElectronicLoss::NonLocal, 2),
        (MeanFreePath::Constant, ElectronicLoss::EquipartitionLsOr, 3),
    ] {
        let mut cfg = BcaConfig::new(5.0, 1.0);
        cfg.mean_free_path = mfp;
        cfg.electronic = el;
        cfg.weak_collisions = weak;
        cfg.seed = 7;
        let beam = Beam {
            ion: Ion::new(18).unwrap(),
            energy_ev: 3.0e3,
            polar_rad: 0.6,
            azimuth_rad: 0.3,
            count: 300,
        };
        let bca = Bca::new(beam, &stack, cfg, &ls, table()).unwrap();
        let t = run_summary(&bca, NM, 20);
        assert_eq!(t.histories, 300);
        assert!(
            t.max_relative_residual < 1e-9,
            "{mfp:?} {el:?} weak {weak}: residual {}",
            t.max_relative_residual
        );
        // Every channel was exercised.
        assert!(t.recoils > 0 && t.sputtered > 0, "{t:?}");
        assert!(t.budget.electronic_nonlocal > 0.0 && t.budget.lattice > 0.0);
        assert!(t.budget.surface_barrier > 0.0 && t.budget.rest > 0.0);
        if el == ElectronicLoss::EquipartitionLsOr {
            assert!(t.budget.electronic_local > 0.0);
        }
        assert_eq!(
            t.primaries_stopped + t.backscattered + t.transmitted,
            t.histories
        );
    }
}

/// Two-sample Kolmogorov-Smirnov statistic between two histograms.
fn ks_statistic(a: &[u64], b: &[u64]) -> f64 {
    let (na, nb) = (a.iter().sum::<u64>() as f64, b.iter().sum::<u64>() as f64);
    let (mut ca, mut cb, mut d) = (0.0, 0.0, 0.0f64);
    for (x, y) in a.iter().zip(b) {
        ca += *x as f64 / na;
        cb += *y as f64 / nb;
        d = d.max((ca - cb).abs());
    }
    d
}

#[test]
fn split_layer_is_equivalent() {
    // A semi-infinite Si target and the same target split at 4 nm (inside
    // the B range), 1e4 ions, fixed seed: the profiles must agree.
    let ls = LindhardScharff::new();
    let one = Stack::semi_infinite(si());
    let two = Stack::new(vec![(si(), 4.0 * NM)], Some(si())).unwrap();
    for (mfp, weak) in [
        (MeanFreePath::Constant, 0),
        (
            MeanFreePath::EnergyDependent {
                min_cm_angle_rad: 0.02,
            },
            0,
        ),
        (MeanFreePath::Constant, 3),
    ] {
        let mut cfg = BcaConfig::new(5.0, 2.0);
        cfg.follow_recoils = false;
        cfg.mean_free_path = mfp;
        cfg.weak_collisions = weak;
        cfg.seed = 11;
        let beam = Beam {
            ion: Ion::new(5).unwrap(),
            energy_ev: 2.0e3,
            polar_rad: 0.3,
            azimuth_rad: 0.0,
            count: 10_000,
        };
        let a = run_summary(
            &Bca::new(beam, &one, cfg, &ls, table()).unwrap(),
            0.5 * NM,
            60,
        );
        let b = run_summary(
            &Bca::new(beam, &two, cfg, &ls, table()).unwrap(),
            0.5 * NM,
            60,
        );
        let n = a.primaries_stopped.min(b.primaries_stopped) as f64;
        assert!(n > 9000.0, "{mfp:?}: too few stopped ({n})");
        // KS critical value at alpha = 0.001 for equal sample sizes.
        let d = ks_statistic(&a.depth_hist, &b.depth_hist);
        let d_crit = 1.95 * (2.0 / n).sqrt();
        assert!(d < d_crit, "{mfp:?}: KS {d} >= {d_crit}");
        // Means within 3 standard errors of their difference.
        let se = (a.depth_std().powi(2) / a.primaries_stopped as f64
            + b.depth_std().powi(2) / b.primaries_stopped as f64)
            .sqrt();
        assert!(
            (a.mean_depth() - b.mean_depth()).abs() < 3.0 * se,
            "{mfp:?}: {} vs {} (se {se})",
            a.mean_depth(),
            b.mean_depth()
        );
        assert_eq!(a.backscattered + a.primaries_stopped, 10_000);

        // Negative control: the same criteria must reject a real difference,
        // a substrate 10 % denser than the film.
        let mut dense = si();
        dense.set_mass_density(1.1 * dense.mass_density()).unwrap();
        let three = Stack::new(vec![(si(), 4.0 * NM)], Some(dense)).unwrap();
        let c = run_summary(
            &Bca::new(beam, &three, cfg, &ls, table()).unwrap(),
            0.5 * NM,
            60,
        );
        let d = ks_statistic(&a.depth_hist, &c.depth_hist);
        assert!(d > d_crit, "{mfp:?}: control not detected, KS {d}");
        assert!((a.mean_depth() - c.mean_depth()).abs() > 3.0 * se);
    }
}

/// Checks that every nonlocal electronic segment lies inside the layer the
/// particle is in, and that local losses and rest positions do too.
struct LayerCheck {
    bounds: Vec<(f64, f64)>,
    segments: u64,
    violations: u64,
}

impl LayerCheck {
    fn inside(&self, layer: usize, x: f64) -> bool {
        let (lo, hi) = self.bounds[layer];
        x >= lo && x <= hi
    }
}

impl BcaTally for LayerCheck {
    fn electronic(&mut self, p: &Particle, from: [f64; 3], _c: ElectronicChannel, _e: f64) {
        self.segments += 1;
        if !(self.inside(p.layer, from[0]) && self.inside(p.layer, p.pos[0])) {
            self.violations += 1;
        }
    }
    fn stopped(&mut self, p: &Particle) {
        if !self.inside(p.layer, p.pos[0]) {
            self.violations += 1;
        }
    }
    fn merge(&mut self, o: Self) {
        self.segments += o.segments;
        self.violations += o.violations;
    }
}

#[test]
fn flights_never_leave_their_layer() {
    let ls = LindhardScharff::new();
    // Thin alternating layers (comparable to a free path) over a substrate.
    let stack = Stack::new(
        vec![
            (si(), 0.3 * NM),
            (sio2(), 0.25 * NM),
            (si(), 0.4 * NM),
            (sio2(), 0.3 * NM),
        ],
        Some(si()),
    )
    .unwrap();
    let bounds: Vec<_> = stack
        .layers()
        .iter()
        .map(|l| (l.front_m(), l.back_m()))
        .collect();
    for mfp in [
        MeanFreePath::Constant,
        MeanFreePath::EnergyDependent {
            min_cm_angle_rad: 0.01,
        },
    ] {
        let mut cfg = BcaConfig::new(5.0, 1.0);
        cfg.mean_free_path = mfp;
        let beam = Beam {
            ion: Ion::new(18).unwrap(),
            energy_ev: 2.0e3,
            polar_rad: 1.0,
            azimuth_rad: 0.0,
            count: 200,
        };
        let bca = Bca::new(beam, &stack, cfg, &ls, table()).unwrap();
        let t = bca
            .run(|| LayerCheck {
                bounds: bounds.clone(),
                segments: 0,
                violations: 0,
            })
            .unwrap();
        assert!(t.segments > 10_000, "{}", t.segments);
        assert_eq!(t.violations, 0, "{mfp:?}");
    }
}

#[test]
fn unset_energy_is_rejected_before_running() {
    let ls = LindhardScharff::new();
    // Oxygen has no default surface binding energy; everything else is set.
    let mut raw = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
    raw.set_displacement_energy_ev(14, 15.0).unwrap();
    raw.set_displacement_energy_ev(8, 15.0).unwrap();
    let stack = Stack::new(vec![(si(), NM)], Some(raw)).unwrap();
    let r = Bca::new(
        Beam::normal(Ion::new(5).unwrap(), 1e3, 10),
        &stack,
        BcaConfig::new(5.0, 1.0),
        &ls,
        table(),
    );
    match r {
        Err(BcaError::EnergyNotSet { layer, z, kind }) => {
            assert_eq!((layer, z, kind), (1, 8, EnergyKind::SurfaceBinding));
        }
        Err(e) => panic!("wrong error {e}"),
        Ok(_) => panic!("accepted a material with an unset energy"),
    }
}

/// A stopping model that fails above a threshold energy.
struct FailsAbove(f64);

impl ElectronicStopping for FailsAbove {
    fn name(&self) -> &'static str {
        "fails-above"
    }
    fn stopping(&self, _ion: &Ion, _z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        if energy_ev > self.0 {
            Err(StoppingError::InvalidEnergy(energy_ev))
        } else {
            Ok(0.0)
        }
    }
    fn validity(&self, _ion: &Ion) -> ValidityRange {
        ValidityRange {
            min_energy_ev: 0.0,
            max_energy_ev: self.0,
        }
    }
}

#[test]
fn stopping_errors_report_the_lowest_failing_history() {
    let model = FailsAbove(500.0);
    let stack = Stack::semi_infinite(si());
    let mut cfg = BcaConfig::new(5.0, 1.0);
    cfg.chunk_size = 3;
    let bca = Bca::new(
        Beam::normal(Ion::new(5).unwrap(), 1e3, 20),
        &stack,
        cfg,
        &model,
        table(),
    )
    .unwrap();
    for threads in [1, 4] {
        match with_threads(threads, || bca.run(|| SummaryTally::new(NM, 4))) {
            Err(BcaError::Stopping { index, .. }) => assert_eq!(index, 0),
            other => panic!("expected a stopping error, got {other:?}"),
        }
    }
}

fn summary_bits(t: &SummaryTally) -> Vec<u64> {
    let b = &t.budget;
    let mut v = vec![
        t.histories,
        t.primaries_stopped,
        t.backscattered,
        t.transmitted,
        t.sputtered,
        t.recoils_transmitted,
        t.recoils,
        t.depth_sum.to_bits(),
        t.depth_sq_sum.to_bits(),
        t.max_relative_residual.to_bits(),
    ];
    for x in [
        b.incident,
        b.electronic_nonlocal,
        b.electronic_local,
        b.lattice,
        b.surface_barrier,
        b.backscattered,
        b.sputtered,
        b.transmitted,
        b.rest,
    ] {
        v.push(x.to_bits());
    }
    v.extend(&t.depth_hist);
    v
}

#[test]
fn bit_identical_across_thread_counts() {
    let ls = LindhardScharff::new();
    let stack = Stack::new(vec![(sio2(), 2.0 * NM)], Some(si())).unwrap();
    let mut cfg = BcaConfig::new(5.0, 1.0);
    cfg.seed = 0xC0FFEE;
    cfg.chunk_size = 16;
    let beam = Beam {
        ion: Ion::new(18).unwrap(),
        energy_ev: 2.0e3,
        polar_rad: 0.4,
        azimuth_rad: 0.0,
        count: 400,
    };
    for weak in [0, 3] {
        cfg.weak_collisions = weak;
        let bca = Bca::new(beam, &stack, cfg, &ls, table()).unwrap();
        let r1 = with_threads(1, || run_summary(&bca, NM, 16));
        assert!(r1.recoils > 0);
        let b1 = summary_bits(&r1);
        for n in [2, 8] {
            let r = with_threads(n, || run_summary(&bca, NM, 16));
            assert_eq!(b1, summary_bits(&r), "weak {weak}: differs on {n} threads");
        }
    }
}

/// Weak collisions (Moller and Eckstein, IPP 9/64 (1988), p. 14) add the
/// nuclear loss of impact parameters beyond `p_max`: in a self-ion cascade
/// the lattice gets more and the electrons less, the weak transfers are
/// reported as `LatticeDeposit::Weak` and never make recoils.
#[test]
fn weak_collisions_move_cascade_energy_from_electrons_to_the_lattice() {
    struct Kinds {
        weak: u64,
        weak_ev: f64,
        recoils: u64,
    }
    impl BcaTally for Kinds {
        fn lattice(&mut self, _at: [f64; 3], _l: usize, kind: LatticeDeposit, e: f64) {
            if kind == LatticeDeposit::Weak {
                self.weak += 1;
                self.weak_ev += e;
            }
        }
        fn recoil(&mut self, _r: &Particle) {
            self.recoils += 1;
        }
        fn merge(&mut self, o: Self) {
            self.weak += o.weak;
            self.weak_ev += o.weak_ev;
            self.recoils += o.recoils;
        }
    }
    let ls = LindhardScharff::new();
    let mut cu = elemental(29);
    cu.set_displacement_energy_ev(29, 5.0).unwrap();
    let stack = Stack::semi_infinite(cu);
    let beam = Beam::normal(Ion::new(29).unwrap(), 1.0e3, 300);
    let mut shares = Vec::new();
    for weak in [0, 1, 3] {
        let mut cfg = BcaConfig::new(2.0, 2.0);
        cfg.seed = 64;
        cfg.weak_collisions = weak;
        let bca = Bca::new(beam, &stack, cfg, &ls, table()).unwrap();
        let t = run_summary(&bca, NM, 8);
        assert!(t.max_relative_residual < 1e-9);
        shares.push(t.budget.electronic_nonlocal / t.budget.incident);
        let k = bca
            .run(|| Kinds {
                weak: 0,
                weak_ev: 0.0,
                recoils: 0,
            })
            .unwrap();
        assert_eq!(k.recoils, t.recoils);
        if weak == 0 {
            assert_eq!(k.weak, 0);
        } else {
            assert!(k.weak > 0 && k.weak_ev > 0.0);
            assert!(k.weak_ev < t.budget.lattice);
        }
    }
    assert!(
        shares[0] > shares[1] && shares[1] > shares[2],
        "electronic share should fall with K: {shares:?}"
    );
}

#[test]
fn escape_channels_behave_physically() {
    let ls = LindhardScharff::new();
    let cfg = BcaConfig::new(5.0, 1.0);
    // A film much thinner than the range transmits most ions; a
    // semi-infinite target transmits none.
    let thin = Stack::new(vec![(si(), 1.0 * NM)], None).unwrap();
    let beam = Beam::normal(Ion::new(5).unwrap(), 5.0e3, 200);
    let t = run_summary(&Bca::new(beam, &thin, cfg, &ls, table()).unwrap(), NM, 4);
    assert!(t.transmitted > 150, "{}", t.transmitted);
    let thick = Stack::semi_infinite(si());
    let t = run_summary(&Bca::new(beam, &thick, cfg, &ls, table()).unwrap(), NM, 4);
    assert_eq!(t.transmitted, 0);
    // Light ions backscatter from a heavy target far more than heavy ions
    // from a light one (kinematics: a heavy projectile cannot be turned
    // back in one collision with a light atom).
    let au = Stack::semi_infinite(elemental(79));
    let c = Stack::semi_infinite(elemental(6));
    let mut no_cascade = cfg;
    no_cascade.follow_recoils = false;
    let he = Beam::normal(Ion::new(2).unwrap(), 2.0e3, 400);
    let xe = Beam::normal(Ion::new(54).unwrap(), 2.0e3, 400);
    let he_au = run_summary(&Bca::new(he, &au, no_cascade, &ls, table()).unwrap(), NM, 4);
    let xe_c = run_summary(&Bca::new(xe, &c, no_cascade, &ls, table()).unwrap(), NM, 4);
    assert!(
        he_au.backscattered > 10 * xe_c.backscattered.max(1),
        "He/Au {} vs Xe/C {}",
        he_au.backscattered,
        xe_c.backscattered
    );
    // Sputtering needs the cascade: Ar on Si sputters with cascades on and
    // cannot without them.
    let ar = Beam::normal(Ion::new(18).unwrap(), 1.0e3, 300);
    let on = run_summary(&Bca::new(ar, &thick, cfg, &ls, table()).unwrap(), NM, 4);
    let off = run_summary(
        &Bca::new(ar, &thick, no_cascade, &ls, table()).unwrap(),
        NM,
        4,
    );
    assert!(on.sputtered > 0);
    assert_eq!(off.sputtered, 0);
}

#[test]
fn free_path_conventions_agree_on_the_range() {
    // With a small minimum angle the energy-dependent free path keeps nearly
    // all of the nuclear stopping, so the mean range should match the
    // constant convention closely.
    let ls = LindhardScharff::new();
    let stack = Stack::semi_infinite(si());
    let beam = Beam::normal(Ion::new(15).unwrap(), 5.0e3, 2000);
    let mut cfg = BcaConfig::new(5.0, 1.0);
    cfg.follow_recoils = false;
    let a = run_summary(&Bca::new(beam, &stack, cfg, &ls, table()).unwrap(), NM, 40);
    cfg.mean_free_path = MeanFreePath::EnergyDependent {
        min_cm_angle_rad: 0.005,
    };
    let b = run_summary(&Bca::new(beam, &stack, cfg, &ls, table()).unwrap(), NM, 40);
    let rel = (a.mean_depth() / b.mean_depth() - 1.0).abs();
    assert!(
        rel < 0.05,
        "constant {} m vs energy-dependent {} m",
        a.mean_depth(),
        b.mean_depth()
    );
}

#[test]
fn escaped_particles_are_outside_the_target() {
    struct Check {
        back: f64,
        bad: u64,
        n: u64,
    }
    impl BcaTally for Check {
        fn escaped(&mut self, p: &Particle, face: Face) {
            self.n += 1;
            let ok = match face {
                Face::Front => p.pos[0] == 0.0 && p.dir[0] < 0.0,
                Face::Back => p.pos[0] == self.back && p.dir[0] > 0.0,
            };
            if !ok {
                self.bad += 1;
            }
        }
        fn merge(&mut self, o: Self) {
            self.bad += o.bad;
            self.n += o.n;
        }
    }
    let ls = LindhardScharff::new();
    let stack = Stack::new(vec![(si(), 2.0 * NM)], None).unwrap();
    let back = stack.back_face_m().unwrap();
    let beam = Beam::normal(Ion::new(18).unwrap(), 2.0e3, 100);
    let bca = Bca::new(beam, &stack, BcaConfig::new(5.0, 1.0), &ls, table()).unwrap();
    let t = bca.run(|| Check { back, bad: 0, n: 0 }).unwrap();
    assert!(t.n > 0);
    assert_eq!(t.bad, 0);
}
