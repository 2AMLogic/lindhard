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
//! [`WINDOW`] times the beam energy (the event that crosses it included;
//! [`Window::Energy`]). A
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
//! # Fixed direction and direction average (#250)
//!
//! Along one off-axis direction the crystal's `R` is above the amorphous one
//! (by 7.6 % at 30/17). Only the average over directions equals a random
//! medium: the rule of angular averages of J. Lindhard, Mat. Fys. Medd. Dan.
//! Vid. Selsk. 34, no. 14 (1965), section 5. So there are three checks:
//!
//! * `direction_averaged_ratio_matches_amorphous`: the 5 % criterion of #226,
//!   on the average over directions, counted over a fixed path
//!   ([`Window::Path`]) instead of the energy window;
//! * `fixed_direction_local_to_nonlocal_ratio_is_above_amorphous`: the
//!   recorded quotient along 30/17 in the energy window;
//! * `fixed_direction_partners_lean_to_small_impact_parameters`: the
//!   impact parameters of the partners along 30/17 ([`PartnerTally`]).
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

/// Tolerance floor of the crystal/amorphous comparison (not to be widened):
/// the 5 % of #226.
const FLOOR: f64 = 0.05;

/// Recorded `R_crystal / R_amorphous` of Ar 20 keV into Si along the fixed
/// direction 30/17 at the same `p_max`, static lattice and 300 K
/// ([`N_HISTORIES`] per side, seed 1, energy window; aarch64 macOS).
const FIXED_DIRECTION_STATIC: f64 = 1.076;
const FIXED_DIRECTION_300_K: f64 = 1.072;

/// Accepted distance from the recorded fixed-direction quotients. The
/// quotient has a standard error of 0.004, and the crystal path is not bit
/// identical across platforms (`docs/architecture.md`), so another platform
/// draws an independent sample: three standard errors of the difference of
/// two samples are 0.017, rounded up.
const FIXED_DIRECTION_TOLERANCE: f64 = 0.02;

/// Path window of the direction-averaged and impact-parameter checks, m:
/// about the path over which Ar 20 keV loses a fifth of its energy in Si
/// (the energy window of #226).
const PATH_WINDOW_M: f64 = 100.0e-10;

/// Histories per direction in the direction average.
const N_PER_DIRECTION: u64 = 16;

/// Points of the spiral that [`isotropic_directions`] filters.
const N_SPIRAL: usize = 1500;

/// Bins of `(p / p_max)^2` in [`PartnerTally`]: a uniform disc fills them
/// equally.
const P2_BINS: usize = 10;

/// Largest accepted deviation of a bin's share from the uniform disc, and of
/// the partners per path from the amorphous value, in the direction average.
const UNIFORM_TOLERANCE: f64 = 0.05;

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

/// When the loss events of a primary stop being counted. The event that
/// crosses the limit is included.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Window {
    /// Until the primary's energy first falls below this, eV (the metric of
    /// #226).
    Energy(f64),
    /// Until the primary's path first exceeds this, m. Every history then
    /// contributes the same path whatever it loses, which the average over
    /// directions needs (see `direction_averaged_ratio_matches_amorphous`).
    Path(f64),
}

/// Per-history local and nonlocal loss of the primary inside the window,
/// accumulated as the sums the ratio estimator needs. `merge` is a plain sum,
/// so the result is bit-identical at any thread count.
#[derive(Debug, Clone, PartialEq)]
struct RatioTally {
    window: Window,
    open: bool,
    cur_local: f64,
    cur_nonlocal: f64,
    cur_path: f64,
    local_events: u64,
    /// Path inside the window, m, all histories.
    path_m: f64,
    n: u64,
    sum_l: f64,
    sum_n: f64,
    sum_ll: f64,
    sum_nn: f64,
    sum_ln: f64,
}

impl RatioTally {
    /// The window of #226: down to [`WINDOW`] times the beam energy.
    fn new(beam_ev: f64) -> Self {
        Self::with_window(Window::Energy(WINDOW * beam_ev))
    }

    fn with_window(window: Window) -> Self {
        Self {
            window,
            open: false,
            cur_local: 0.0,
            cur_nonlocal: 0.0,
            cur_path: 0.0,
            local_events: 0,
            path_m: 0.0,
            n: 0,
            sum_l: 0.0,
            sum_n: 0.0,
            sum_ll: 0.0,
            sum_nn: 0.0,
            sum_ln: 0.0,
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
        self.cur_path = 0.0;
    }

    fn electronic(
        &mut self,
        p: &Particle,
        from: [f64; 3],
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
            ElectronicChannel::NonLocal => {
                self.cur_nonlocal += energy_ev;
                let d2: f64 = (0..3).map(|k| (p.pos[k] - from[k]).powi(2)).sum();
                self.cur_path += d2.sqrt();
            }
        }
        let closed = match self.window {
            Window::Energy(threshold_ev) => p.energy_ev < threshold_ev,
            Window::Path(limit_m) => self.cur_path > limit_m,
        };
        if closed {
            self.open = false;
        }
    }

    fn end_history(&mut self, _index: u64, _budget: &EnergyBudget) {
        let (l, n) = (self.cur_local, self.cur_nonlocal);
        self.path_m += self.cur_path;
        self.n += 1;
        self.sum_l += l;
        self.sum_n += n;
        self.sum_ll += l * l;
        self.sum_nn += n * n;
        self.sum_ln += l * n;
    }

    fn merge(&mut self, o: Self) {
        self.local_events += o.local_events;
        self.path_m += o.path_m;
        self.n += o.n;
        self.sum_l += o.sum_l;
        self.sum_n += o.sum_n;
        self.sum_ll += o.sum_ll;
        self.sum_nn += o.sum_nn;
        self.sum_ln += o.sum_ln;
    }
}

/// One run (local + nonlocal loss, primary only, recoils not followed) into
/// the tally made by `new_tally`; amorphous when `crystal` is `None`.
fn run_tally<T: BcaTally>(
    b: Beam,
    crystal: Option<CrystalTarget>,
    threads: Option<usize>,
    seed: u64,
    new_tally: impl Fn() -> T + Sync,
) -> T {
    let st = Stack::semi_infinite(si());
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = seed;
    cfg.follow_recoils = false;
    cfg.electronic = ElectronicLoss::EquipartitionLsOr;
    let mut bca = Bca::new(b, &st, cfg, &ls, table()).unwrap();
    if let Some(c) = crystal {
        bca = bca.with_crystal(c, &[0]).unwrap();
    }
    let go = || bca.run(&new_tally).unwrap();
    match threads {
        None => go(),
        Some(n) => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(go),
    }
}

/// One run of `R` in the energy window of #226, seed 1.
fn run_ratio(b: Beam, crystal: Option<CrystalTarget>, threads: Option<usize>) -> RatioTally {
    let e0 = b.energy_ev;
    run_tally(b, crystal, threads, 1, || RatioTally::new(e0))
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

/// Ar 20 keV into Si along one fixed off-axis direction, 30 degrees tilt and
/// 17 degrees twist (far from every low-index axis and plane; the off-axis
/// case of `tests/crystal_bca.rs`): the crystal's local share is **above**
/// the amorphous one at the same `p_max`, by the recorded
/// [`FIXED_DIRECTION_STATIC`] (static lattice) and [`FIXED_DIRECTION_300_K`].
///
/// #226 asked for agreement within `max(3 sigma, 5 %)` here, and #250 found
/// why a single direction does not give it. The crystal equals the random
/// medium only in the average over directions (Lindhard's rule of angular
/// averages; `direction_averaged_ratio_matches_amorphous` asserts that at the
/// 5 % of #226). The channeling directions lie below that average, so a
/// direction outside every channel lies above it. This test is therefore a
/// recorded-value regression check, not agreement with a reference.
///
/// Two more cases are measured and printed but not asserted: the 7/22
/// orientation, which is 2.6 degrees from a {100} and a {110} plane, and the
/// crystal at its default `p_max` (nn), to show how much the result depends
/// on the cutoff.
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn fixed_direction_local_to_nonlocal_ratio_is_above_amorphous() {
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
    for (thermal, recorded) in [
        (false, FIXED_DIRECTION_STATIC),
        (true, FIXED_DIRECTION_300_K),
    ] {
        let label = if thermal {
            "30/17, 300 K"
        } else {
            "30/17, static"
        };
        let cr = measured(
            &format!("crystal {label}"),
            &run(30.0, 17.0, Some(crystal_si(30.0, 17.0, p_max, thermal))),
        );
        checks.push((label, recorded, quotient(cr, am)));
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

    for (label, recorded, (q, sq)) in checks {
        println!(
            "crystal/amorphous {label}: {q:.4} +- {sq:.4} \
             (recorded {recorded}, tolerance {FIXED_DIRECTION_TOLERANCE})"
        );
        assert!(
            (q - recorded).abs() <= FIXED_DIRECTION_TOLERANCE,
            "{label}: R_crystal / R_amorphous = {q} +- {sq}, recorded {recorded} \
             +- {FIXED_DIRECTION_TOLERANCE}"
        );
        assert!(
            q - 1.0 > 3.0 * sq,
            "{label}: R_crystal / R_amorphous = {q} +- {sq} is not above 1"
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

/// The impact parameters of the collision partners (#250), on top of a
/// [`RatioTally`] with the same window: for the primary inside the window,
/// every partner reported by [`BcaTally::partner`], binned in
/// `(p / p_max)^2`, with the local loss it took. `merge` is a plain sum.
#[derive(Debug, Clone, PartialEq)]
struct PartnerTally {
    ratio: RatioTally,
    p_max_m: f64,
    /// Bins of the partners of the current step whose local loss has not been
    /// reported yet (the engine reports the partners of a step, then their
    /// local losses in the same order).
    queue: Vec<usize>,
    /// Partners, and their summed local loss in eV, per bin.
    partners: [u64; P2_BINS],
    local_ev: [f64; P2_BINS],
    /// Partners beyond `p_max` (none are expected).
    outside: u64,
    /// Collision steps with 1, 2, 3 and 4 or more partners.
    steps: [u64; 4],
}

impl PartnerTally {
    fn new(window: Window, p_max_m: f64) -> Self {
        Self {
            ratio: RatioTally::with_window(window),
            p_max_m,
            queue: Vec::new(),
            partners: [0; P2_BINS],
            local_ev: [0.0; P2_BINS],
            outside: 0,
            steps: [0; 4],
        }
    }

    fn n_partners(&self) -> u64 {
        self.partners.iter().sum()
    }

    /// Partners in bins `range`, as a share of all partners relative to the
    /// uniform disc (1 for a uniform disc), with its binomial standard error.
    fn share(&self, range: std::ops::Range<usize>) -> (f64, f64) {
        let n = self.n_partners() as f64;
        let uniform = range.len() as f64 / P2_BINS as f64;
        let f = self.partners[range].iter().sum::<u64>() as f64 / n;
        (f / uniform, (f * (1.0 - f) / n).sqrt() / uniform)
    }

    /// Partners per unit path, in units of the amorphous `N pi p_max^2`.
    fn partners_per_path(&self) -> f64 {
        let n_pi_p2 = si().atom_number_density() * PI * self.p_max_m * self.p_max_m;
        self.n_partners() as f64 / self.ratio.path_m / n_pi_p2
    }

    /// Mean local loss of a partner in bin `k`, eV.
    fn local_per_partner(&self, k: usize) -> f64 {
        self.local_ev[k] / self.partners[k] as f64
    }
}

impl BcaTally for PartnerTally {
    fn begin_history(&mut self, index: u64) {
        self.ratio.begin_history(index);
        self.queue.clear();
    }

    fn partner(&mut self, p: &Particle, impact_parameter_m: f64, partners: usize) {
        if !(p.is_primary() && self.ratio.open) {
            return;
        }
        if self.queue.len() >= partners {
            // A new step: a partner of the previous one took no local loss.
            self.queue.clear();
        }
        if self.queue.is_empty() {
            self.steps[partners.clamp(1, 4) - 1] += 1;
        }
        let u = (impact_parameter_m / self.p_max_m).powi(2);
        if u > 1.0 + 1e-9 {
            self.outside += 1;
        }
        let bin = ((u * P2_BINS as f64) as usize).min(P2_BINS - 1);
        self.partners[bin] += 1;
        self.queue.push(bin);
    }

    fn electronic(
        &mut self,
        p: &Particle,
        from: [f64; 3],
        channel: ElectronicChannel,
        energy_ev: f64,
    ) {
        if p.is_primary() && self.ratio.open {
            match channel {
                ElectronicChannel::Local if !self.queue.is_empty() => {
                    let bin = self.queue.remove(0);
                    self.local_ev[bin] += energy_ev;
                }
                ElectronicChannel::Local => {}
                // A flight separates two collision steps.
                ElectronicChannel::NonLocal => self.queue.clear(),
            }
        }
        self.ratio.electronic(p, from, channel, energy_ev);
    }

    fn end_history(&mut self, index: u64, budget: &EnergyBudget) {
        self.ratio.end_history(index, budget);
    }

    fn merge(&mut self, o: Self) {
        self.ratio.merge(o.ratio);
        for k in 0..P2_BINS {
            self.partners[k] += o.partners[k];
            self.local_ev[k] += o.local_ev[k];
        }
        for k in 0..4 {
            self.steps[k] += o.steps[k];
        }
        self.outside += o.outside;
    }
}

/// Print the impact-parameter diagnostic of one run, per bin against
/// `reference` (the amorphous run) when there is one.
fn print_partners(label: &str, t: &PartnerTally, reference: Option<&PartnerTally>) {
    let (r, s) = t.ratio.ratio();
    let steps = t.steps.iter().sum::<u64>() as f64;
    println!(
        "{label}: R = {r:.4} +- {s:.4}, {} partners ({} beyond p_max), partners per path / \
         (N pi p_max^2) = {:.4}; steps with 1 / 2 / 3 / 4+ partners: {:.4} / {:.4} / {:.4} / {:.4}",
        t.n_partners(),
        t.outside,
        t.partners_per_path(),
        t.steps[0] as f64 / steps,
        t.steps[1] as f64 / steps,
        t.steps[2] as f64 / steps,
        t.steps[3] as f64 / steps
    );
    for k in 0..P2_BINS {
        let (f, sf) = t.share(k..k + 1);
        let vs = reference.map_or(String::new(), |a| {
            format!(
                " ({:.4} x amorphous)",
                t.local_per_partner(k) / a.local_per_partner(k)
            )
        });
        println!(
            "  (p/p_max)^2 in [{:.1}, {:.1}): share / uniform {f:.4} +- {sf:.4}, local loss per \
             partner {:.3} eV{vs}",
            k as f64 / P2_BINS as f64,
            (k + 1) as f64 / P2_BINS as f64,
            t.local_per_partner(k)
        );
    }
}

/// Beam directions (tilt, twist in degrees) spread evenly over the solid
/// angle of the directions that are closer to [100] than to any other cube
/// axis. That region is one sixth of the sphere and, by the cubic symmetry, a
/// fair sample of all directions. The points are a golden-ratio spiral in
/// `(cos tilt, twist)`, uniform in solid angle over the cone out to <111>
/// (54.7 degrees), of which those inside the region are kept: a deterministic
/// set, with no random draw.
fn isotropic_directions(n_spiral: usize) -> Vec<(f64, f64)> {
    let cos_max = 1.0 / 3f64.sqrt();
    let golden = (5f64.sqrt() - 1.0) / 2.0;
    (0..n_spiral)
        .filter_map(|k| {
            let cos_t = 1.0 - (k as f64 + 0.5) / n_spiral as f64 * (1.0 - cos_max);
            let twist = 2.0 * PI * (k as f64 * golden).fract();
            let tan_t = (1.0 - cos_t * cos_t).sqrt() / cos_t;
            // Closer to the surface normal [100] than to either in-plane axis.
            (tan_t * twist.cos().abs().max(twist.sin().abs()) <= 1.0)
                .then(|| (cos_t.acos().to_degrees(), twist.to_degrees()))
        })
        .collect()
}

/// **The random-direction criterion of #226, in the form that holds.** Ar
/// 20 keV into Si, averaged over all beam directions
/// ([`isotropic_directions`], [`N_PER_DIRECTION`] histories each, a different
/// seed per direction, the same seeds on the amorphous side) and counted over
/// a fixed path ([`PATH_WINDOW_M`]): the crystal splits the electronic loss
/// between the local and the nonlocal half as the amorphous model does at the
/// same `p_max`, within `max(3 sigma, 5 %)`; its partners fill the disc of
/// radius `p_max` uniformly, and there are as many of them per unit path.
/// Static lattice and 300 K.
///
/// This is the rule of angular averages of J. Lindhard, Mat. Fys. Medd. Dan.
/// Vid. Selsk. 34, no. 14 (1965), section 5, eqs. (5.1)-(5.3): a quantity
/// linear in the probability of being at a point of the crystal "has the same
/// angular average as in a random system", a system "with the same density,
/// but without directional effects". The paper limits the rule to small
/// energy loss, hence the fixed, short path here: the energy window of #226
/// lets a channeled ion travel farther, so it weights the channels more.
///
/// The standard error pools all histories about the common `R`, so it also
/// contains the spread between directions and is an upper bound.
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn direction_averaged_ratio_matches_amorphous() {
    let p_max = amorphous_p_max(&si());
    let (z, e) = (18u8, 2.0e4);
    let window = Window::Path(PATH_WINDOW_M);
    let dirs = isotropic_directions(N_SPIRAL);
    assert!(dirs.len() > 1000, "{} directions", dirs.len());
    let run = |k: usize, crystal: Option<CrystalTarget>| {
        let (tilt, twist) = dirs[k];
        run_tally(
            beam(z, e, tilt, twist, N_PER_DIRECTION),
            crystal,
            None,
            k as u64 + 1,
            || PartnerTally::new(window, p_max),
        )
    };
    let mut am = PartnerTally::new(window, p_max);
    for k in 0..dirs.len() {
        am.merge(run(k, None));
    }
    println!(
        "{} directions, {N_PER_DIRECTION} histories each",
        dirs.len()
    );
    print_partners("amorphous, all directions", &am, None);
    let r_am = measured("amorphous, all directions", &am.ratio);

    for thermal in [false, true] {
        let state = if thermal { "300 K" } else { "static" };
        let mut cr = PartnerTally::new(window, p_max);
        // R of each direction, to show the spread (16 histories: noisy).
        let mut each = Vec::new();
        for (k, &(tilt, twist)) in dirs.iter().enumerate() {
            let t = run(k, Some(crystal_si(tilt, twist, p_max, thermal)));
            each.push((t.ratio.ratio().0, tilt, twist));
            cr.merge(t);
        }
        let label = format!("crystal {state}, all directions");
        print_partners(&label, &cr, Some(&am));
        each.sort_by(|a, b| a.0.total_cmp(&b.0));
        let below = each.iter().filter(|x| x.0 < r_am.0).count();
        println!(
            "  R per direction: below the amorphous R in {below} of {}; quantiles",
            each.len()
        );
        for q in [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0] {
            let (r, tilt, twist) = each[((each.len() - 1) as f64 * q) as usize];
            println!("    {q:.2}: {r:.3} (at {tilt:.1}/{twist:.1})");
        }

        let (q, sq) = quotient(measured(&label, &cr.ratio), r_am);
        let tol = (3.0 * sq).max(FLOOR);
        println!(
            "crystal/amorphous {state}, all directions: {q:.4} +- {sq:.4} (tolerance {tol:.4})"
        );
        assert!(
            (q - 1.0).abs() <= tol,
            "{state}: direction-averaged R_crystal / R_amorphous = {q} +- {sq}, tolerance {tol}"
        );
        assert_eq!(cr.outside, 0, "{state}: partners beyond p_max");
        let per_path = cr.partners_per_path() / am.partners_per_path();
        assert!(
            (per_path - 1.0).abs() <= UNIFORM_TOLERANCE,
            "{state}: partners per path, crystal / amorphous = {per_path}"
        );
        for k in 0..P2_BINS {
            let (f, sf) = cr.share(k..k + 1);
            assert!(
                (f - 1.0).abs() <= (3.0 * sf).max(UNIFORM_TOLERANCE),
                "{state}: bin {k} of (p/p_max)^2 holds {f} +- {sf} of its uniform share"
            );
        }
    }
}

/// The diagnostic of #250: which impact parameters occur along one fixed
/// direction (30/17, static lattice, path window). The crystal makes as many
/// collisions per unit path as the amorphous model, and a partner at a given
/// impact parameter takes the same local loss in both (the two engines share
/// the loss function). But the partners are not spread uniformly over the
/// disc of radius `p_max`: the inner three tenths of its area hold more than
/// their share and the outer three tenths less. That shift is the whole
/// excess of `R` along this direction. It vanishes in the average over
/// directions (`direction_averaged_ratio_matches_amorphous`), so it is where
/// the ions go in the lattice along this direction, not a counting error.
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn fixed_direction_partners_lean_to_small_impact_parameters() {
    let p_max = amorphous_p_max(&si());
    let window = Window::Path(PATH_WINDOW_M);
    let run = |crystal| {
        run_tally(
            beam(18, 2.0e4, 30.0, 17.0, N_HISTORIES),
            crystal,
            None,
            1,
            || PartnerTally::new(window, p_max),
        )
    };
    let am = run(None);
    print_partners("amorphous 30/17", &am, None);
    let cr = run(Some(crystal_si(30.0, 17.0, p_max, false)));
    print_partners("crystal 30/17, static", &cr, Some(&am));

    // The amorphous disc is uniform, bin by bin.
    for k in 0..P2_BINS {
        let (f, sf) = am.share(k..k + 1);
        assert!(
            (f - 1.0).abs() <= 4.0 * sf,
            "amorphous bin {k}: {f} +- {sf}"
        );
    }
    // Same collision rate, same loss at a given impact parameter.
    assert_eq!(cr.outside, 0);
    let per_path = cr.partners_per_path() / am.partners_per_path();
    println!("partners per path, crystal / amorphous: {per_path:.4}");
    assert!(
        (per_path - 1.0).abs() <= 0.03,
        "partners per path {per_path}"
    );
    for k in 0..P2_BINS {
        let q = cr.local_per_partner(k) / am.local_per_partner(k);
        assert!(
            (q - 1.0).abs() <= 0.05,
            "bin {k}: local loss per partner {q}"
        );
    }
    // The crystal's partners lean inwards.
    let (inner, s_inner) = cr.share(0..3);
    let (outer, s_outer) = cr.share(P2_BINS - 3..P2_BINS);
    println!(
        "crystal share / uniform: inner three bins {inner:.4} +- {s_inner:.4}, outer three \
         {outer:.4} +- {s_outer:.4}"
    );
    assert!(inner - 1.0 > 3.0 * s_inner, "inner {inner} +- {s_inner}");
    assert!(1.0 - outer > 3.0 * s_outer, "outer {outer} +- {s_outer}");
    let (q, sq) = quotient(cr.ratio.ratio(), am.ratio.ratio());
    println!("crystal/amorphous 30/17, static, path window: {q:.4} +- {sq:.4}");
    assert!(q - 1.0 > 3.0 * sq, "R_crystal / R_amorphous = {q} +- {sq}");
}
