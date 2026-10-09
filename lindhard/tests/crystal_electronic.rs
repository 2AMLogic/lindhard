//! The local (Oen-Robinson) electronic loss in the crystal flight model
//! (issue #226): validation against the amorphous engine, the channeling
//! reduction, and the run-metadata flag for the unverified constants.
//!
//! The local loss itself is implemented in `Bca::crystal_collide`
//! (`lindhard::ion::bca::crystal`): under `ElectronicLoss::EquipartitionLsOr`
//! every lattice site within `p_max` of the path that becomes a partner takes
//! the Oen-Robinson loss at its closest approach. These tests check that the
//! split between the local and the nonlocal half is right.
//!
//! # Metric
//!
//! `R = E_local / E_nonlocal`, summed over histories, for the **primary ion
//! only** (recoils are not followed), over the first stretch of its path:
//! the loss events are counted until the primary's energy first falls below
//! [`WINDOW`] times the beam energy (the event that crosses it included). A
//! ratio, not a range, so the nuclear scattering differences between the
//! crystal and the amorphous model drop out: per unit path length the
//! nonlocal half is the same in both, and the local half counts the partners
//! within `p_max`.
//!
//! The error of `R` is the ratio-estimator (delta method) standard error over
//! histories, `sigma_R = sqrt(Var(L_i - R N_i) / n) / mean(N_i)`.
//!
//! # Comparator
//!
//! Both models cut the local half at a finite radius, so the comparison uses
//! the same radius on both sides: the amorphous constant-path radius
//! `p_max = (pi N^(2/3))^(-1/2)` (1.53 Angstrom for Si; it follows the
//! density) is set explicitly as the crystal's `p_max_m`. Same species,
//! material, `BcaConfig` (no weak collisions; they are amorphous-only) and
//! seed on both sides.
//!
//! The statistical checks need thousands of ions and a release build, so they
//! are `#[ignore]`d; run them with
//! `cargo test --release -p lindhard --test crystal_electronic -- --ignored --nocapture`.
//! The measured numbers are in the module docs of `lindhard::ion::bca::crystal`.

use std::f64::consts::PI;
use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, CrystalTarget, ElectronicChannel, ElectronicLoss, EnergyBudget,
    Particle, Thermal,
};
use lindhard::ion::crystal::debye::THETA_D_SI;
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::oen_robinson::OR_CONSTANTS_UNVERIFIED;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;

/// The loss events of a primary are counted while its energy is at least this
/// fraction of the beam energy.
const WINDOW: f64 = 0.8;

/// Histories per side in the random-direction check (the issue asks for at
/// least 2000).
const N_HISTORIES: u64 = 4000;

/// Histories per run in the axial check: a channeled ion makes few local
/// losses, so `R` is noisier there (<110> needs more than 4000 for a 1 %
/// error).
const N_AXIAL: u64 = 8000;

/// Largest accepted relative standard error of `R` on each side.
const MAX_REL_ERROR: f64 = 0.01;

/// Tolerance floor of the crystal/amorphous comparison (not to be widened).
const FLOOR: f64 = 0.05;

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

/// Silicon at the tabulated (crystal) density with `E_d` = 15 eV, as in
/// `tests/crystal_bca.rs`.
fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    if m.surface_binding_energy_ev(14).is_err() {
        m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    }
    m
}

/// The amorphous constant-path radius of `m`, m.
fn amorphous_p_max(m: &Material) -> f64 {
    1.0 / (PI * m.atom_number_density().powf(2.0 / 3.0)).sqrt()
}

/// Nearest-neighbour distance of silicon, m.
fn si_nn() -> f64 {
    Lattice::silicon().lattice_constant() * 3f64.sqrt() / 4.0
}

fn beam(z: u8, e: f64, tilt_deg: f64, twist_deg: f64, count: u64) -> Beam {
    Beam {
        ion: Ion::new(z).unwrap(),
        energy_ev: e,
        polar_rad: tilt_deg.to_radians(),
        azimuth_rad: twist_deg.to_radians(),
        count,
    }
}

/// Si cut on (100), reference [010], matching `beam(.., tilt, twist, ..)`,
/// with search radius `p_max`.
fn crystal_si(tilt_deg: f64, twist_deg: f64, p_max: f64, thermal: bool) -> CrystalTarget {
    let lat = Lattice::silicon();
    let o = Orientation::new(
        &lat,
        [1, 0, 0],
        [0, 1, 0],
        tilt_deg.to_radians(),
        twist_deg.to_radians(),
        0.0,
    )
    .unwrap();
    let mut c = CrystalTarget::new(lat, o);
    c.p_max_m = Some(p_max);
    if thermal {
        c = c.with_thermal(Thermal::new(300.0, THETA_D_SI));
    }
    c
}

/// Per-history local and nonlocal loss of the primary inside the window,
/// accumulated as the sums the ratio estimator needs. `merge` is a plain sum,
/// so the result is bit-identical at any thread count.
#[derive(Debug, Clone, Default, PartialEq)]
struct RatioTally {
    /// Energy below which the window closes, eV.
    threshold_ev: f64,
    open: bool,
    cur_local: f64,
    cur_nonlocal: f64,
    local_events: u64,
    n: u64,
    sum_l: f64,
    sum_n: f64,
    sum_ll: f64,
    sum_nn: f64,
    sum_ln: f64,
}

impl RatioTally {
    fn new(beam_ev: f64) -> Self {
        Self {
            threshold_ev: WINDOW * beam_ev,
            ..Self::default()
        }
    }

    /// `R` and its standard error.
    fn ratio(&self) -> (f64, f64) {
        let n = self.n as f64;
        let r = self.sum_l / self.sum_n;
        let mean_n = self.sum_n / n;
        // Sample variance of L_i - R N_i (its mean is zero by construction).
        let var = (self.sum_ll - 2.0 * r * self.sum_ln + r * r * self.sum_nn) / (n - 1.0);
        (r, (var.max(0.0) / n).sqrt() / mean_n)
    }
}

impl BcaTally for RatioTally {
    fn begin_history(&mut self, _index: u64) {
        self.open = true;
        self.cur_local = 0.0;
        self.cur_nonlocal = 0.0;
    }

    fn electronic(
        &mut self,
        p: &Particle,
        _from: [f64; 3],
        channel: ElectronicChannel,
        energy_ev: f64,
    ) {
        if !(p.is_primary() && self.open) {
            return;
        }
        match channel {
            ElectronicChannel::Local => {
                self.cur_local += energy_ev;
                self.local_events += 1;
            }
            ElectronicChannel::NonLocal => self.cur_nonlocal += energy_ev,
        }
        if p.energy_ev < self.threshold_ev {
            self.open = false;
        }
    }

    fn end_history(&mut self, _index: u64, _budget: &EnergyBudget) {
        let (l, n) = (self.cur_local, self.cur_nonlocal);
        self.n += 1;
        self.sum_l += l;
        self.sum_n += n;
        self.sum_ll += l * l;
        self.sum_nn += n * n;
        self.sum_ln += l * n;
    }

    fn merge(&mut self, o: Self) {
        self.local_events += o.local_events;
        self.n += o.n;
        self.sum_l += o.sum_l;
        self.sum_n += o.sum_n;
        self.sum_ll += o.sum_ll;
        self.sum_nn += o.sum_nn;
        self.sum_ln += o.sum_ln;
    }
}

/// One run of `R` (local + nonlocal, primary only), amorphous when `crystal`
/// is `None`.
fn run_ratio(b: Beam, crystal: Option<CrystalTarget>, threads: Option<usize>) -> RatioTally {
    let st = Stack::semi_infinite(si());
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 1;
    cfg.follow_recoils = false;
    cfg.electronic = ElectronicLoss::EquipartitionLsOr;
    let e0 = b.energy_ev;
    let mut bca = Bca::new(b, &st, cfg, &ls, table()).unwrap();
    if let Some(c) = crystal {
        bca = bca.with_crystal(c, &[0]).unwrap();
    }
    let go = || bca.run(|| RatioTally::new(e0)).unwrap();
    match threads {
        None => go(),
        Some(n) => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(go),
    }
}

/// `R` with its error, checked to the accepted relative error.
fn measured(label: &str, t: &RatioTally) -> (f64, f64) {
    let (r, s) = t.ratio();
    println!(
        "{label}: R = {r:.5} +- {s:.5} ({:.2} %), n = {}, local events {}",
        100.0 * s / r,
        t.n,
        t.local_events
    );
    assert!(t.n >= 2000, "{label}: n = {}", t.n);
    assert!(
        s < MAX_REL_ERROR * r,
        "{label}: standard error {s} is not below {MAX_REL_ERROR} of R = {r}; raise n"
    );
    (r, s)
}

/// `R_crystal / R_amorphous` and its standard error.
fn quotient(c: (f64, f64), a: (f64, f64)) -> (f64, f64) {
    let q = c.0 / a.0;
    (q, q * ((c.1 / c.0).powi(2) + (a.1 / a.0).powi(2)).sqrt())
}

/// Ar 20 keV into Si in a random direction: the crystal should split the
/// electronic loss between the local and the nonlocal half as the amorphous
/// model does, at the same `p_max`, within `max(3 sigma, 5 %)`. "Random" is
/// the direction 30 degrees tilt, 17 degrees twist, far from every low-index
/// axis and plane (the off-axis case of `tests/crystal_bca.rs`), asserted for
/// the static lattice and at 300 K. Two cases are measured and printed but not
/// asserted: the 7/22 orientation, which looks random but is 2.6 degrees from a
/// {100} and a {110} plane, and the crystal at its default `p_max` (nn), to
/// show how much the result depends on the cutoff.
///
/// **This check fails today: known gap, #250.** The crystal ratio is 5-8 %
/// above the amorphous one (30/17 static: 1.076 +- 0.004). The 5 % floor is
/// the acceptance criterion of #226 and is deliberately not widened; #250
/// records the measurement and the diagnosis so far.
#[test]
#[ignore = "statistical, and fails today (known gap #250); run with --release -- --ignored"]
fn random_direction_local_to_nonlocal_ratio_matches_amorphous() {
    let p_max = amorphous_p_max(&si());
    println!(
        "p_max: amorphous {:.4} A, crystal {:.4} A (set equal); nn {:.4} A",
        p_max * 1e10,
        p_max * 1e10,
        si_nn() * 1e10
    );
    let (z, e) = (18u8, 2.0e4);
    let run = |tilt, twist, crystal| run_ratio(beam(z, e, tilt, twist, N_HISTORIES), crystal, None);
    // Measure everything first, so that a failing run still prints it all.
    let am = measured("amorphous 30/17", &run(30.0, 17.0, None));
    let mut checks = Vec::new();
    for thermal in [false, true] {
        let label = if thermal {
            "30/17, 300 K"
        } else {
            "30/17, static"
        };
        let cr = measured(
            &format!("crystal {label}"),
            &run(30.0, 17.0, Some(crystal_si(30.0, 17.0, p_max, thermal))),
        );
        checks.push((label, quotient(cr, am)));
    }
    let am722 = measured("amorphous 7/22", &run(7.0, 22.0, None));
    let cr722 = measured(
        "crystal 7/22, static",
        &run(7.0, 22.0, Some(crystal_si(7.0, 22.0, p_max, false))),
    );
    let (q, sq) = quotient(cr722, am722);
    println!("crystal/amorphous 7/22, static (not asserted): {q:.4} +- {sq:.4}");
    let nn = measured(
        "crystal 30/17, static, p_max = nn",
        &run(30.0, 17.0, Some(crystal_si(30.0, 17.0, si_nn(), false))),
    );
    let (q, sq) = quotient(nn, am);
    println!(
        "crystal (p_max = nn) / amorphous (p_max = {:.3} A), 30/17 (not asserted): \
         {q:.4} +- {sq:.4}",
        p_max * 1e10
    );

    for (label, (q, sq)) in checks {
        let tol = (3.0 * sq).max(FLOOR);
        println!("crystal/amorphous {label}: {q:.4} +- {sq:.4} (tolerance {tol:.4})");
        assert!(
            (q - 1.0).abs() <= tol,
            "{label}: R_crystal / R_amorphous = {q} +- {sq}, tolerance {tol} (#250)"
        );
    }
}

/// Along <110> and <100> of Si the channeled ion passes fewer sites within
/// `p_max`, so the local half is a smaller share of the loss: `R_axial` is
/// below `R_random` (the crystal at 30/17) by more than three combined
/// standard errors. Same metric, ion and `p_max` as the random-direction
/// check; static lattice and 300 K.
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn axial_local_to_nonlocal_ratio_is_below_random() {
    let p_max = amorphous_p_max(&si());
    let (z, e) = (18u8, 2.0e4);
    for thermal in [false, true] {
        let state = if thermal { "300 K" } else { "static" };
        let random = measured(
            &format!("crystal 30/17, {state}"),
            &run_ratio(
                beam(z, e, 30.0, 17.0, N_AXIAL),
                Some(crystal_si(30.0, 17.0, p_max, thermal)),
                None,
            ),
        );
        // <110>: [100] normal tilted 45 degrees towards the in-plane [010];
        // <100>: along the [100] surface normal.
        for (axis, tilt, twist) in [("<110>", 45.0, 0.0), ("<100>", 0.0, 0.0)] {
            let label = format!("{axis}, {state}");
            let ax = measured(
                &format!("crystal {label}"),
                &run_ratio(
                    beam(z, e, tilt, twist, N_AXIAL),
                    Some(crystal_si(tilt, twist, p_max, thermal)),
                    None,
                ),
            );
            let sigma = (ax.1 * ax.1 + random.1 * random.1).sqrt();
            println!(
                "{label}: R_axial / R_random = {:.4}, below by {:.1} sigma",
                ax.0 / random.0,
                (random.0 - ax.0) / sigma
            );
            assert!(
                random.0 - ax.0 > 3.0 * sigma,
                "{label}: R_axial = {} +- {}, R_random = {} +- {}",
                ax.0,
                ax.1,
                random.0,
                random.1
            );
        }
    }
}

/// The ratio tally of a crystal run with the local loss (static and 300 K) is
/// bit-identical on 1, 2 and 8 threads, and both channels are sampled.
#[test]
fn local_loss_ratio_is_identical_on_1_2_8_threads() {
    let p_max = amorphous_p_max(&si());
    for thermal in [false, true] {
        let run = |threads| {
            run_ratio(
                beam(18, 2.0e4, 30.0, 17.0, 100),
                Some(crystal_si(30.0, 17.0, p_max, thermal)),
                Some(threads),
            )
        };
        let r1 = run(1);
        assert!(r1.sum_l > 0.0 && r1.sum_n > 0.0, "{r1:?}");
        for n in [2, 8] {
            assert_eq!(r1, run(n), "thermal {thermal}: {n} threads differ");
        }
    }
}

/// The run metadata flags the unverified Oen-Robinson constants for a run
/// that uses them, and only then; the flag follows the one constant in
/// `oen_robinson.rs`.
#[test]
fn metadata_flags_unverified_oen_robinson_constants() {
    let st = Stack::semi_infinite(si());
    let ls = LindhardScharff::new();
    let meta = |electronic| {
        let mut cfg = BcaConfig::new(5.0, 2.0);
        cfg.electronic = electronic;
        Bca::new(beam(18, 2.0e4, 7.0, 22.0, 1), &st, cfg, &ls, table())
            .unwrap()
            .with_crystal(crystal_si(7.0, 22.0, si_nn(), false), &[0])
            .unwrap()
            .crystal_metadata()
    };
    let local = meta(ElectronicLoss::EquipartitionLsOr);
    assert_eq!(local.len(), 1);
    assert_eq!(
        local[0].electronic_constants_unverified,
        OR_CONSTANTS_UNVERIFIED
    );
    let json = serde_json::to_string(&local).unwrap();
    assert!(
        json.contains("\"electronic_constants_unverified\":true"),
        "{json}"
    );

    let nonlocal = meta(ElectronicLoss::NonLocal);
    assert!(!nonlocal[0].electronic_constants_unverified);
    let json = serde_json::to_string(&nonlocal).unwrap();
    assert!(!json.contains("electronic_constants_unverified"), "{json}");
}
