//! The Lindhard partition of a self-ion cascade, solved deterministically
//! from the integral equation of Lindhard, Nielsen, Scharff and Thomsen with
//! the engine's own inputs, against the engine's electronic share (issue #64).
//!
//! # The equation
//!
//! J. Lindhard, V. Nielsen, M. Scharff and P. V. Thomsen, "Integral equations
//! governing radiation effects", Mat. Fys. Medd. Dan. Vid. Selsk. 33 (10)
//! (1963), read at
//! <http://gymarkiv.sdu.dk/MFM/kdvs/mfm%2030-39/mfm-33-10.pdf>. For an atom
//! of the medium (`Z1 = Z2`) of energy `E`, the energy `nu(E)` finally given
//! to atomic motion obeys the homogeneous part of their eq. (2.7) (p. 14;
//! "`nu(E)` ... is normally a solution of the homogeneous part of equation
//! (2.1)", p. 20, with the boundary condition `nu(E)/E -> 1` as `E -> 0`):
//!
//! ```text
//! S_e(E) nu'(E) = Int dsigma_n(E, T) { nu(E - T) - nu(E) + nu(T) }
//! ```
//!
//! (approximations (A) to (D) of their § 2: electrons make no recoils,
//! no binding loss, small electronic transfers, nuclear and electronic
//! collisions separate). The electronic share is `1 - nu(E)/E`.
//!
//! The equation is solved here with the **same inputs as the engine**: the
//! engine's ZBL scattering table for `theta(eps, beta)`, Lindhard-Scharff
//! `S_e`, the density, and the engine's thresholds: an atom below the cutoff
//! `E_c` stops (`nu(E) = E` for `E < E_c`), and a transfer `T <= E_d` makes no
//! recoil (it stays in the lattice, `nu -> T`). Two cross sections are used:
//!
//! * [`Kernel::Full`]: every impact parameter, `p` up to `beta = 100`
//!   screening lengths (the engine table's edge; the part of the nuclear
//!   stopping beyond it is below `1e-12` at 1 eV for Cu and Si). This is the
//!   physical partition for these inputs.
//! * [`Kernel::Weak`]`(K)`: what the engine samples with `K` weak collisions
//!   (`ion::bca` module docs): recoiling collisions for `p <= p_max`, and for
//!   `p_max < p <= p_max sqrt(K + 1)` collisions whose transfer stays in the
//!   lattice (`nu -> T`, no recoil). `K = 0` is the engine without weak
//!   collisions.
//!
//! # Discretisation
//!
//! As in `lss.rs`: an implicit march upward on a geometric energy grid from
//! `E_c`, `nu` interpolated linearly in `E` between known values (linear in
//! the unknown `nu(E_i)` in the top interval, which both `E - T` and `T` can
//! fall in), `nu'` a backward difference, Simpson's rule in `ln p`. Two grids
//! (100 and 200 points per decade) are Richardson-extrapolated and their
//! difference is the reported discretisation uncertainty. The solver is
//! verified against the closed result of LNST for power-law scattering
//! ([`powerlaw_coefficient`]).

use std::f64::consts::PI;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, ElectronicChannel, EnergyBudget, Face, LatticeDeposit, Particle,
};
use lindhard::ion::damage::damage_energy_ev;
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{theta_quadrature, ScatteringTable};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::{ElectronicStopping, Ion};
use lindhard::material::Material;
use lindhard::units::J_PER_EV;

use crate::report::{num, pct, Check};

/// Which nuclear cross section the partition equation uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kernel {
    /// All impact parameters (to the table edge), every transfer above `E_d`
    /// a recoil.
    Full,
    /// The engine's sampling with this many weak collisions.
    Weak(u8),
}

/// A self-ion problem in a monatomic target, in the engine's conventions.
pub struct Problem<'a> {
    pub z: u8,
    /// Atomic mass, u (projectile and target).
    pub m: f64,
    /// Atom density, m^-3.
    pub n: f64,
    /// Cutoff energy `E_c`, eV; also the displacement energy `E_d`.
    pub cutoff_ev: f64,
    /// Electronic stopping cross section per atom, eV m^2, at an energy in eV.
    pub se: &'a dyn Fn(f64) -> f64,
    pub table: &'a ScatteringTable,
}

impl Problem<'_> {
    fn p_max(&self) -> f64 {
        1.0 / (PI * self.n.powf(2.0 / 3.0)).sqrt()
    }
}

/// Simpson nodes `(p, weight)` for `Int 2 pi p dp` over `[lo, hi]` in `ln p`,
/// at least `per_decade` intervals per decade.
fn simpson_p(lo: f64, hi: f64, per_decade: usize) -> Vec<(f64, f64)> {
    let mut m = (((hi / lo).log10() * per_decade as f64).ceil() as usize).max(2);
    m += m % 2;
    let du = (hi / lo).ln() / m as f64;
    (0..=m)
        .map(|k| {
            let p = lo * (k as f64 * du).exp();
            let s = if k == 0 || k == m {
                1.0
            } else if k % 2 == 1 {
                4.0
            } else {
                2.0
            };
            (p, s * du / 3.0 * 2.0 * PI * p * p)
        })
        .collect()
}

/// March `S_e nu' = Int dsigma {nu(E - T) - nu(E) + h(T)}` from `cutoff` to
/// `e0` on `per_decade` points per decade and return `nu(e0)`.
/// `collisions(E)` lists `(T, weight, recoil)`: `h(T) = nu(T)` if `recoil`
/// and `T > e_d`, else `T`. Weights are cross sections (any unit, the same as
/// `se` divided by an energy).
fn march<C, S>(e0: f64, cutoff: f64, e_d: f64, per_decade: usize, se: S, collisions: C) -> f64
where
    C: Fn(f64) -> Vec<(f64, f64, bool)>,
    S: Fn(f64) -> f64,
{
    let k = ((e0 / cutoff).log10() * per_decade as f64).ceil() as usize;
    let ln_r = (e0 / cutoff).ln() / k as f64;
    let grid: Vec<f64> = (0..=k).map(|i| cutoff * (i as f64 * ln_r).exp()).collect();
    let mut nu = vec![0.0f64; k + 1];
    nu[0] = cutoff;
    for i in 1..=k {
        let e = grid[i];
        // nu(x) as (coefficient of the unknown nu(E_i), known part).
        let value = |x: f64, nu: &[f64]| -> (f64, f64) {
            if x < cutoff {
                return (0.0, x);
            }
            let j = (((x / cutoff).ln() / ln_r).floor() as usize).min(i - 1);
            let s = ((x - grid[j]) / (grid[j + 1] - grid[j])).clamp(0.0, 1.0);
            if j + 1 == i {
                (s, (1.0 - s) * nu[j])
            } else {
                (0.0, (1.0 - s) * nu[j] + s * nu[j + 1])
            }
        };
        let (mut diag, mut known) = (0.0, 0.0);
        for (t, w, recoil) in collisions(e) {
            let (c1, k1) = value(e - t, &nu);
            let (c2, k2) = if recoil && t > e_d {
                value(t, &nu)
            } else {
                (0.0, t)
            };
            diag += w * (1.0 - c1 - c2);
            known += w * (k1 + k2);
        }
        let de = grid[i] - grid[i - 1];
        let s = se(e);
        nu[i] = (known + s * nu[i - 1] / de) / (diag + s / de);
    }
    nu[k]
}

/// `nu(e0)` on one grid for the ZBL problem.
fn solve_grid(pr: &Problem, e0: f64, kernel: Kernel, per_decade: usize) -> f64 {
    let (z, m) = (f64::from(pr.z), pr.m);
    let pot = Potential::new(Screening::ZblUniversal, z, z);
    let a = pot.screening_length();
    let eps_per_ev = pot.reduced_energy(Potential::cm_energy(J_PER_EV, m, m));
    let p_max = pr.p_max();
    let mut nodes: Vec<(f64, f64, bool)> = Vec::new();
    match kernel {
        Kernel::Full => {
            let hi = pr.table.spec().beta_max * a;
            nodes.extend(
                simpson_p(1e-6 * p_max, hi, 40)
                    .into_iter()
                    .map(|(p, w)| (p, w, true)),
            );
        }
        Kernel::Weak(k) => {
            nodes.extend(
                simpson_p(1e-6 * p_max, p_max, 40)
                    .into_iter()
                    .map(|(p, w)| (p, w, true)),
            );
            if k > 0 {
                let hi = p_max * f64::from(k + 1).sqrt();
                nodes.extend(
                    simpson_p(p_max, hi, 400)
                        .into_iter()
                        .map(|(p, w)| (p, w, false)),
                );
            }
        }
    }
    let collisions = |e: f64| -> Vec<(f64, f64, bool)> {
        let eps = e * eps_per_ev;
        nodes
            .iter()
            .map(|&(p, w, recoil)| {
                let beta = p / a;
                let theta = pr
                    .table
                    .theta(eps, beta)
                    .unwrap_or_else(|| theta_quadrature(Screening::ZblUniversal, eps, beta));
                let s = (0.5 * theta).sin();
                // Equal masses: T = E sin^2(theta/2).
                ((e * s * s).min(e), w, recoil)
            })
            .collect()
    };
    march(
        e0,
        pr.cutoff_ev,
        pr.cutoff_ev,
        per_decade,
        pr.se,
        collisions,
    )
}

/// Electronic share `1 - nu(e0)/e0`, Richardson-extrapolated, and its
/// discretisation uncertainty.
pub fn electronic_share(pr: &Problem, e0: f64, kernel: Kernel) -> (f64, f64) {
    let coarse = solve_grid(pr, e0, kernel, 100);
    let fine = solve_grid(pr, e0, kernel, 200);
    let nu = 2.0 * fine - coarse;
    (1.0 - nu / e0, ((fine - coarse) / e0).abs())
}

/// Solver verification on the one case LNST solve in closed form for their
/// eq. (2.7): power-law scattering with `s = 2` (`dsigma = C E^(-1/2)
/// T^(-3/2) dT`, so `S_n` is constant) and `S_e = S_n xi`,
/// `xi = (E/E_xi)^(1/2)`. Their p. 23: the solution of (2.7) is the series
/// `eta(E) = E - nu(E) = a1 E^(3/2) E_xi^(-1/2) - ...` with
/// `a1 = 4/(3 pi - 6) = 1.17`, while approximations (E), (E') and eq. (4.4)
/// give `a1 = 1`, `8/7 = 1.14` and `16/13 = 1.23`. Returns `(eta/E)/xi`
/// from the solver at `E = 1e-6 E_xi` (`xi = 1e-3`, cutoff `1e-12 E_xi`), and
/// `a1`.
pub fn powerlaw_coefficient() -> (f64, f64) {
    let a1 = 4.0 / (3.0 * PI - 6.0);
    let (e0, cutoff) = (1e-6, 1e-12);
    // T nodes as fractions x = T/E, Simpson in ln x over [1e-10, 1].
    let m = 400usize;
    let lo: f64 = 1e-10;
    let du = (1.0 / lo).ln() / m as f64;
    let fr: Vec<(f64, f64)> = (0..=m)
        .map(|k| {
            let s = if k == 0 || k == m {
                1.0
            } else if k % 2 == 1 {
                4.0
            } else {
                2.0
            };
            (lo * (k as f64 * du).exp(), s * du / 3.0)
        })
        .collect();
    // S_n = Int T dsigma = 1 with dsigma = (1/2) E^(-1/2) T^(-3/2) dT,
    // i.e. weight (1/2) E^(-1/2) T^(-1/2) per unit ln T.
    let collisions = |e: f64| -> Vec<(f64, f64, bool)> {
        fr.iter()
            .map(|&(x, sw)| {
                let t = x * e;
                (t, sw * 0.5 / (e * t).sqrt(), true)
            })
            .collect()
    };
    let se = |e: f64| e.sqrt();
    let solve = |n: usize| march(e0, cutoff, 0.0, n, se, collisions);
    let nu = 2.0 * solve(200) - solve(100);
    ((1.0 - nu / e0) / e0.sqrt(), a1)
}

/// Per-run sums for the engine side: the electronic share and the measured
/// size of the BCA's step conventions (see [`checks`]).
#[derive(Debug, Default, Clone)]
struct ShareTally {
    cutoff: f64,
    histories: u64,
    electronic: f64,
    kept: f64,
    /// Per-history share, for the standard error.
    sum_f: f64,
    sum_f2: f64,
    hist_el: f64,
    /// Nonlocal loss re-evaluated at the energy at the end of its step.
    end_of_step: f64,
    /// Electronic loss charged below the cutoff by flights that cross it.
    overshoot: f64,
    /// (energy at the start of the pending flight, its loss).
    pending: Option<(f64, f64)>,
    /// The last recoil reported (energy, position), to tell a recoil that is
    /// not followed (stopped at once) from the end of the moving particle.
    last_recoil: Option<(f64, [f64; 3])>,
    /// Depth, m, within which a weak partner can lie beyond the front
    /// surface: the outer radius of the weak rings, `p_max sqrt(K + 1)`.
    near_depth_m: f64,
    /// Weak-collision transfers made at depths below `near_depth_m`.
    near_weak: f64,
}

impl ShareTally {
    fn new(cutoff: f64, near_depth_m: f64) -> Self {
        Self {
            cutoff,
            near_depth_m,
            ..Self::default()
        }
    }

    /// Close the pending flight: its step ends with energy `e_end`.
    /// Lindhard-Scharff `S_e` is proportional to `E^(1/2)`.
    fn close(&mut self, e_end: f64) {
        if let Some((e_start, de)) = self.pending.take() {
            self.end_of_step += de * (e_end.max(0.0) / e_start).sqrt().min(1.0);
        }
    }
}

impl BcaTally for ShareTally {
    fn begin_history(&mut self, _i: u64) {
        self.pending = None;
        self.last_recoil = None;
        self.hist_el = 0.0;
    }
    fn electronic(&mut self, p: &Particle, _from: [f64; 3], c: ElectronicChannel, de: f64) {
        self.hist_el += de;
        if c != ElectronicChannel::NonLocal {
            self.end_of_step += de;
            return;
        }
        let e_start = p.energy_ev + de;
        // A truncated flight continues; its next segment starts where this
        // one ended, which closes this one at the right energy.
        self.close(e_start);
        self.pending = Some((e_start, de));
        if e_start >= self.cutoff && p.energy_ev < self.cutoff {
            self.overshoot += self.cutoff - p.energy_ev;
        }
    }
    fn recoil(&mut self, r: &Particle) {
        self.last_recoil = Some((r.energy_ev, r.pos));
    }
    fn lattice(&mut self, at: [f64; 3], _layer: usize, kind: LatticeDeposit, energy_ev: f64) {
        if kind == LatticeDeposit::Weak && at[0] < self.near_depth_m {
            self.near_weak += energy_ev;
        }
    }
    fn stopped(&mut self, p: &Particle) {
        if let Some((e, pos)) = self.last_recoil.take() {
            if p.generation > 0 && e == p.energy_ev && pos == p.pos {
                return;
            }
        }
        self.close(p.energy_ev);
    }
    fn escaped(&mut self, p: &Particle, _face: Face) {
        self.close(p.energy_ev);
    }
    fn end_history(&mut self, _i: u64, b: &EnergyBudget) {
        self.close(0.0);
        let kept = b.incident - b.backscattered - b.sputtered - b.transmitted - b.surface_barrier;
        let el = b.electronic_nonlocal + b.electronic_local;
        self.histories += 1;
        self.electronic += el;
        self.kept += kept;
        let f = el / kept;
        self.sum_f += f;
        self.sum_f2 += f * f;
    }
    fn merge(&mut self, o: Self) {
        self.histories += o.histories;
        self.electronic += o.electronic;
        self.kept += o.kept;
        self.sum_f += o.sum_f;
        self.sum_f2 += o.sum_f2;
        self.end_of_step += o.end_of_step;
        self.overshoot += o.overshoot;
        self.near_weak += o.near_weak;
    }
}

impl ShareTally {
    fn share(&self) -> f64 {
        self.electronic / self.kept
    }
    fn se(&self) -> f64 {
        let n = self.histories as f64;
        let m = self.sum_f / n;
        ((self.sum_f2 / n - m * m) * n / (n - 1.0) / n).sqrt()
    }
}

/// Weak collisions used by the asserted rows: TRIDYN's maximum (IPP 9/64,
/// p. 14).
const WEAK: u8 = 3;
/// Cutoffs and displacement energy, eV: low, so atoms are set in motion and
/// followed down to where the BCA stops being meaningful.
const CUTOFF_EV: f64 = 1.0;

struct Case {
    id: &'static str,
    sym: &'static str,
    z: u8,
    e0: f64,
}

const CASES: [Case; 4] = [
    Case {
        id: "damage.cascade_electronic_share.cu_1k",
        sym: "Cu",
        z: 29,
        e0: 1e3,
    },
    Case {
        id: "damage.cascade_electronic_share.cu_10k",
        sym: "Cu",
        z: 29,
        e0: 1e4,
    },
    Case {
        id: "damage.cascade_electronic_share.si_1k",
        sym: "Si",
        z: 14,
        e0: 1e3,
    },
    Case {
        id: "damage.cascade_electronic_share.si_10k",
        sym: "Si",
        z: 14,
        e0: 1e4,
    },
];

pub fn checks(table: &ScatteringTable, quick: bool) -> Vec<Check> {
    let mut out = Vec::new();
    let ls = LindhardScharff::new();

    let (got, a1) = powerlaw_coefficient();
    out.push(Check::at_most(
        "damage.partition_solver",
        format!(
            "Partition solver, power-law scattering s = 2, S_e/S_n = xi = 1e-3: (eta/E)/xi = {} vs LNST a1 = 4/(3 pi - 6) = {} (|rel. dev.|)",
            num(got),
            num(a1)
        ),
        (got / a1 - 1.0).abs(),
        0.01,
        pct,
        "LNST (1963) p. 23, solution of eq. (2.7); 1 % is under half the gap to the nearest other approximation, (E') with 8/7 (2.2 %), so the check shows the solver solves (2.7)",
    ));

    let mut robinson = Vec::new();
    let mut control_ok = true;
    let mut control_worst = f64::INFINITY;
    for c in &CASES {
        let mut mat = Material::from_atom_fractions(&[(c.z, 1.0)], None).unwrap();
        let n = mat.atom_number_density();
        let m = lindhard::elements::element(c.z).unwrap().atomic_weight;
        let ion = Ion::new(c.z).unwrap();
        let se = |e: f64| ls.stopping(&ion, c.z, e).unwrap() / J_PER_EV;
        let pr = Problem {
            z: c.z,
            m,
            n,
            cutoff_ev: CUTOFF_EV,
            se: &se,
            table,
        };
        let (full, d_full) = electronic_share(&pr, c.e0, Kernel::Full);
        let (weak, d_weak) = electronic_share(&pr, c.e0, Kernel::Weak(WEAK));

        // Engine: E_d = cutoffs = 1 eV, E_b = 0, and surface barriers so high
        // that nothing leaves: a particle reaching a face is reflected
        // specularly, which for an amorphous medium is the mirror image of
        // its path in an infinite medium, the setting of the equation. One
        // exception: a weak partner beyond the face is skipped (TRIDYN's
        // surface test, `ion::bca` module docs), where the mirror image would
        // have one. That can only happen within the outer ring radius
        // p_max sqrt(K + 1) of the face, and is bounded by the `surface` term
        // of the tolerance below.
        mat.set_displacement_energy_ev(c.z, CUTOFF_EV).unwrap();
        mat.set_lattice_binding_energy_ev(c.z, 0.0).unwrap();
        mat.set_surface_binding_energy_ev(c.z, 1e9).unwrap();
        let stack = Stack::semi_infinite(mat);
        let count = match (quick, c.e0 > 5e3) {
            (true, false) => 200,
            (true, true) => 40,
            (false, false) => 2000,
            (false, true) => 400,
        };
        let near_depth_m = pr.p_max() * f64::from(WEAK + 1).sqrt();
        let run = |weak_collisions: u8| {
            let mut cfg = BcaConfig::new(CUTOFF_EV, CUTOFF_EV);
            cfg.weak_collisions = weak_collisions;
            cfg.primary_surface_binding_ev = 1e9;
            cfg.seed = 0x64 + u64::from(c.z);
            Bca::new(Beam::normal(ion, c.e0, count), &stack, cfg, &ls, table)
                .unwrap()
                .run(|| ShareTally::new(CUTOFF_EV, near_depth_m))
                .unwrap()
        };
        let t = run(WEAK);
        let share = t.share();
        let kept = t.kept / (c.e0 * t.histories as f64);
        // Model allowance, every term computed or measured on this run:
        // * truncation: the weak rings stop at p_max sqrt(K + 1) and make no
        //   recoils, |LNST(K) - LNST(full)|;
        // * step: the engine charges a flight's electronic loss at the energy
        //   the flight starts with, before the step's collisions, where the
        //   equation's continuous slowing down would charge it along the way;
        //   the same losses at the energy the step ends with bracket it from
        //   below. Plus the loss charged below the cutoff by flights that
        //   cross it, and the primary's first flight R l (on average l/2 of
        //   loss missing, N S_e(E0) l / 2);
        // * surface: the weak transfers the skip drops. A partner is at most
        //   p_max sqrt(K + 1) from the path, and its azimuth is uniform, so
        //   for a collision at depth x >= 0 the chance that it lies beyond the
        //   face is at most 1/2: the dropped transfers are at most those made
        //   within that depth, which the tally measures. Dropped nuclear loss
        //   stays with the particle, and at most all of it later goes to
        //   electrons, so to first order it moves the share by at most that
        //   sum over the kept energy. (Full-statistics development runs:
        //   switching the skip off moves the share by 0.0002 to 0.0020, 10
        //   to 25 times less than this bound.)
        let trunc = (weak - full).abs() + d_weak;
        let first = n * se(c.e0) * n.powf(-1.0 / 3.0) * 0.5 / c.e0;
        let step = (t.electronic - t.end_of_step + t.overshoot) / t.kept + first;
        let surface = t.near_weak / t.kept;
        let tol = 3.0 * t.se() + d_full + trunc + step + surface;
        out.push(Check::at_most(
            c.id,
            format!(
                "{s} {e} keV -> {s}, K = {WEAK} weak collisions: engine electronic share {} vs LNST partition {} (|abs. diff.|)",
                num(share),
                num(full),
                s = c.sym,
                e = c.e0 / 1e3
            ),
            (share - full).abs(),
            tol,
            num,
            format!(
                "share = (el_nl + el_loc)/(incident - backscattered - sputtered - transmitted - surface_barrier) = {} kept, {count} ions; LNST eq. (2.7) with the engine's ZBL, LS, E_c = E_d = 1 eV; tol = 3 SE {} + solver {} + truncation {} (LNST with K = {WEAK} rings: {}) + step {} + surface {}",
                num(kept),
                num(3.0 * t.se()),
                num(d_full),
                num(trunc),
                num(weak),
                num(step),
                num(surface)
            ),
        ));

        // Negative control: the engine without weak collisions (bit for bit
        // the engine before #64) against the same reference and tolerance.
        let t0 = run(0);
        let dev0 = (t0.share() - full).abs();
        control_ok &= dev0 > tol;
        control_worst = control_worst.min(dev0 / tol);

        robinson.push(format!(
            "{} {} keV: {} (LNST {})",
            c.sym,
            c.e0 / 1e3,
            num(1.0 - damage_energy_ev(c.e0, f64::from(c.z), m, f64::from(c.z), m) / c.e0),
            num(full)
        ));
    }
    out.push(Check::holds(
        "damage.cascade_electronic_share.control",
        "Negative control: without weak collisions (the engine before #64) every case above misses its tolerance",
        control_ok,
        format!("smallest |deviation| / tolerance {}", num(control_worst)),
        "shows the partition check detects the dropped nuclear loss beyond p_max",
    ));
    out.push(Check::info(
        "damage.cascade_electronic_share.robinson",
        "Context: Robinson's fit of the Lindhard partition (Thomas-Fermi cross section, no cutoffs), 1 - T_dam/T",
        robinson.join("; "),
        "not comparable: with ZBL S_n and LS S_e the ratio S_e/S_n rises again below about 100 eV, while LNST's Thomas-Fermi cross section behaves as power-law scattering with s = 3 at low energy (p. 26), so S_e/S_n keeps falling (p. 21: it tends to zero for s < 4); the cascade tail spends much of its electronic loss there",
    ));
    out
}
