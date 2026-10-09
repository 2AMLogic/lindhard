//! Mean ranges from the LSS first-moment transport equation, solved
//! deterministically, against the Monte Carlo engine.
//!
//! # The equation
//!
//! For an ion of energy `E` in an infinite amorphous medium of atom density
//! `N`, the mean projected range `R1(E)` (the mean penetration along the
//! initial direction) obeys the first-moment backward transport equation of
//! J. Lindhard, M. Scharff and H. E. Schiott, "Range concepts and heavy ion
//! ranges", Mat. Fys. Medd. Dan. Vid. Selsk. 33 (14) (1963):
//!
//! ```text
//! 1 = N Int dsigma(p) [R1(E) - cos(psi) R1(E - T)] + N S_e(E) dR1/dE
//! ```
//!
//! with `T` the energy transfer and `psi` the laboratory deflection of the
//! ion in a collision at impact parameter `p`. Replacing `cos(psi)` by 1 gives
//! the mean total path length `L(E)`. The derivation is one line: over a
//! short step `dx` the ion either travels on (losing `N S_e dx`) or collides
//! (probability `N dsigma dx`), after which its remaining mean penetration
//! along the old axis is `cos(psi) R1(E - T)` by rotational symmetry.
//!
//! The equation is solved here with the **same physics inputs as the
//! engine**: the engine's ZBL scattering table for `theta(eps, beta)`, `p`
//! over the constant-free-path disc `p <= p_max`, Lindhard-Scharff `S_e`, the
//! density, and `R1 = 0` below the primary cutoff energy. The two-body
//! kinematics that turn `theta` into `cos(psi)` and `T` are *not* taken from
//! the engine; they are written out independently in [`cos_psi_lab`] and
//! [`energy_transfer`], so a bug in the engine's kinematics does not cancel.
//! So it tests the transport (deflection kinematics, flight, rotation,
//! bookkeeping), not the physics inputs: a disagreement means a transport
//! bug, an agreement says nothing about whether the inputs match experiment.
//!
//! # Discretisation
//!
//! Implicit march upward on a geometric energy grid from the cutoff: at each
//! grid energy `E_i`, `R1(E - T)` is interpolated linearly in `E` between
//! already known values, except in the top interval where it is linear in the
//! unknown `R1(E_i)`, and `dR1/dE` is a backward difference. That makes each
//! step one linear equation in `R1(E_i)`. The scheme is first order in the
//! grid spacing; two grids (`n` and `2n` points per decade) are solved and
//! Richardson-extrapolated, and their difference is reported as the
//! discretisation uncertainty. The impact-parameter integral is Simpson's rule
//! in `ln p` from `1e-6 p_max` to `p_max` (the omitted head-on disc carries a
//! fraction `1e-12` of the cross section).
//!
//! # What the engine does differently
//!
//! * Flights have the constant length `l = N^(-1/3)` (first flight `R l`)
//!   rather than exponential lengths with the same mean; for the first moment
//!   this matters only at the very end of the path (order `l`).
//! * The target is semi-infinite, so backscattered ions leave and are not in
//!   the engine's mean; the equation keeps them at negative depth. Cases are
//!   chosen with backscatter below 1 % where the bias is negligible, and the
//!   backscatter fraction is reported.

use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, BcaTally, Beam, Face, Particle};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{theta_quadrature, ScatteringTable};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::{ElectronicStopping, Ion, StoppingError, ValidityRange};
use lindhard::material::Material;
use lindhard::units::J_PER_EV;
use std::f64::consts::PI;

use crate::report::{num, pct, spct, Check};

const NM: f64 = 1e-9;

/// Which first moment to solve for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Moment {
    /// Mean projected range (along the incident direction).
    Projected,
    /// Mean total path length.
    Path,
    /// Negative control: the projected range with the centre-of-mass angle
    /// used as the laboratory deflection (a classic BCA bug). The engine must
    /// disagree with this one.
    ProjectedCmAngleBug,
}

/// One ion-target problem, in the engine's conventions.
pub(crate) struct Problem<'a> {
    pub ion: Ion,
    /// Target atomic number and mass (u); monatomic targets only.
    pub z2: u8,
    pub m2: f64,
    /// Atom density, m^-3.
    pub n: f64,
    /// Primary cutoff, eV.
    pub cutoff_ev: f64,
    /// Electronic stopping cross section per atom, eV m^2, at an energy in eV.
    pub se: &'a dyn Fn(f64) -> f64,
    pub table: &'a ScatteringTable,
    /// Weak collisions per collision step in the engine
    /// (`BcaConfig::weak_collisions`).
    pub weak: u8,
}

impl Problem<'_> {
    /// Constant-free-path `p_max = (pi N^(2/3))^(-1/2)` (as in the engine).
    pub(crate) fn p_max(&self) -> f64 {
        1.0 / (PI * self.n.powf(2.0 / 3.0)).sqrt()
    }

    /// Largest impact parameter the engine samples: `p_max`, or
    /// `p_max sqrt(K + 1)` with `K` weak collisions, whose partners fill the
    /// annuli out to that radius uniformly, one atom per annulus per flight
    /// (Moller and Eckstein, IPP 9/64 (1988), p. 14, eq. (26)). For the
    /// primary's range only the deflection and the transfer matter, so the
    /// weak collisions enter the equation as ordinary collisions.
    pub(crate) fn p_cut(&self) -> f64 {
        self.p_max() * f64::from(self.weak + 1).sqrt()
    }

    fn theta(&self, eps: f64, beta: f64) -> f64 {
        self.table
            .theta(eps, beta)
            .unwrap_or_else(|| theta_quadrature(Screening::ZblUniversal, eps, beta))
    }
}

/// Simpson nodes `(p, weight)` for `Int 2 pi p dp` over `[1e-6 p_max, p_cut]`
/// in `ln p`.
fn p_nodes(p_max: f64, p_cut: f64, per_decade: usize) -> Vec<(f64, f64)> {
    let lo = p_max * 1e-6;
    let decades = (p_cut / lo).log10();
    let mut m = (decades * per_decade as f64).ceil() as usize;
    m += m % 2;
    let du = (p_cut / lo).ln() / m as f64;
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

/// Cosine of the laboratory deflection `psi` of a projectile of mass `m1`
/// scattered elastically through centre-of-mass angle `theta` by an
/// initially stationary atom of mass `m2`.
///
/// Classical two-body elastic kinematics (H. Goldstein, *Classical
/// Mechanics*, 2nd ed. (1980), Sec. 3.11; W. Eckstein, *Computer
/// Simulation of Ion-Solid Interactions*, Springer (1991), Ch. 2). With
/// `A = m2 / m1`,
///
/// ```text
/// cos(psi) = (1 + A cos(theta)) / sqrt(1 + 2 A cos(theta) + A^2)
/// ```
///
/// Written out here, deliberately not imported from the engine, so that the
/// LSS solution is independent of the engine's deflection kinematics.
pub(crate) fn cos_psi_lab(theta: f64, m1: f64, m2: f64) -> f64 {
    let a = m2 / m1;
    let c = theta.cos();
    (1.0 + a * c) / (1.0 + 2.0 * a * c + a * a).sqrt()
}

/// Energy `T` transferred to an initially stationary atom of mass `m2` by a
/// projectile of mass `m1` and energy `e` scattered through centre-of-mass
/// angle `theta`: `T = 4 m1 m2 / (m1 + m2)^2 * E sin^2(theta / 2)`
/// (same references as [`cos_psi_lab`]). Clamped to `e` against rounding.
pub(crate) fn energy_transfer(e: f64, theta: f64, m1: f64, m2: f64) -> f64 {
    let s = (0.5 * theta).sin();
    (4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2)) * e * s * s).min(e)
}

/// Solve on one grid of `per_decade` energies per decade; returns `f(E0)`, m.
pub(crate) fn solve_grid(pr: &Problem, e0: f64, moment: Moment, per_decade: usize) -> f64 {
    let m1 = pr.ion.mass_amu();
    let pot = Potential::new(
        Screening::ZblUniversal,
        f64::from(pr.ion.z()),
        f64::from(pr.z2),
    );
    let a = pot.screening_length();
    let eps_per_ev = pot.reduced_energy(Potential::cm_energy(J_PER_EV, m1, pr.m2));
    let nodes = p_nodes(pr.p_max(), pr.p_cut(), 40);

    let k = ((e0 / pr.cutoff_ev).log10() * per_decade as f64).ceil() as usize;
    let ln_r = (e0 / pr.cutoff_ev).ln() / k as f64;
    let grid: Vec<f64> = (0..=k)
        .map(|i| pr.cutoff_ev * (i as f64 * ln_r).exp())
        .collect();
    let mut f = vec![0.0f64; k + 1];
    for i in 1..=k {
        let e = grid[i];
        let eps = e * eps_per_ev;
        let (mut diag, mut known) = (0.0, 0.0);
        for &(p, w) in &nodes {
            let theta = pr.theta(eps, p / a);
            let t = energy_transfer(e, theta, m1, pr.m2);
            let c = match moment {
                Moment::Projected => cos_psi_lab(theta, m1, pr.m2),
                Moment::Path => 1.0,
                Moment::ProjectedCmAngleBug => theta.cos(),
            };
            diag += w;
            let ep = e - t;
            if ep < pr.cutoff_ev {
                continue;
            }
            let j = (((ep / pr.cutoff_ev).ln() / ln_r).floor() as usize).min(i - 1);
            let s = ((ep - grid[j]) / (grid[j + 1] - grid[j])).clamp(0.0, 1.0);
            if j + 1 == i {
                diag -= w * c * s;
                known += w * c * (1.0 - s) * f[j];
            } else {
                known += w * c * ((1.0 - s) * f[j] + s * f[j + 1]);
            }
        }
        let de = grid[i] - grid[i - 1];
        let nse = pr.n * (pr.se)(e);
        f[i] = (1.0 + pr.n * known + nse * f[i - 1] / de) / (pr.n * diag + nse / de);
    }
    f[k]
}

/// Richardson-extrapolated solution and its discretisation uncertainty, m.
pub(crate) fn solve(pr: &Problem, e0: f64, moment: Moment) -> (f64, f64) {
    let coarse = solve_grid(pr, e0, moment, 200);
    let fine = solve_grid(pr, e0, moment, 400);
    (2.0 * fine - coarse, (fine - coarse).abs())
}

/// Continuous-slowing-down path length `Int dE / (N (S_n + S_e))` from the
/// cutoff, with `S_n` integrated over the same `p <= p_cut` disc. Not a mean
/// path: it ignores the fluctuation of the nuclear energy loss.
pub(crate) fn csda_path(pr: &Problem, e0: f64) -> f64 {
    let m1 = pr.ion.mass_amu();
    let pot = Potential::new(
        Screening::ZblUniversal,
        f64::from(pr.ion.z()),
        f64::from(pr.z2),
    );
    let a = pot.screening_length();
    let eps_per_ev = pot.reduced_energy(Potential::cm_energy(J_PER_EV, m1, pr.m2));
    let nodes = p_nodes(pr.p_max(), pr.p_cut(), 40);
    let sn = |e: f64| -> f64 {
        nodes
            .iter()
            .map(|&(p, w)| w * energy_transfer(e, pr.theta(e * eps_per_ev, p / a), m1, pr.m2))
            .sum()
    };
    // Simpson in ln E: Int E dE/(E N S) d(ln E).
    let m = 2 * ((e0 / pr.cutoff_ev).log10() * 100.0).ceil() as usize;
    let du = (e0 / pr.cutoff_ev).ln() / m as f64;
    let mut acc = 0.0;
    for k in 0..=m {
        let e = pr.cutoff_ev * (k as f64 * du).exp();
        let g = e / (pr.n * (sn(e) + (pr.se)(e)));
        let s = if k == 0 || k == m {
            1.0
        } else if k % 2 == 1 {
            4.0
        } else {
            2.0
        };
        acc += s * g;
    }
    acc * du / 3.0
}

/// Primary-only range tally: projected depth and path length of stopped
/// primaries. The path is the polyline through the primary's collision sites
/// (every collision with `T > 0` reports its site through `lattice` or
/// `recoil`; a `T = 0` collision does not deflect, so skipping it does not
/// change the polyline). Only valid with `follow_recoils = false`, where the
/// primary is the only moving particle.
#[derive(Debug, Default, Clone)]
pub(crate) struct RangeTally {
    last: [f64; 3],
    path: f64,
    pub stopped: u64,
    pub backscattered: u64,
    pub sum_x: f64,
    pub sum_x2: f64,
    pub sum_l: f64,
    pub sum_l2: f64,
    pub sum_xl: f64,
}

impl RangeTally {
    fn vertex(&mut self, at: [f64; 3]) {
        let d: f64 = (0..3).map(|k| (at[k] - self.last[k]).powi(2)).sum();
        self.path += d.sqrt();
        self.last = at;
    }

    pub(crate) fn mean_x(&self) -> f64 {
        self.sum_x / self.stopped as f64
    }
    pub(crate) fn mean_l(&self) -> f64 {
        self.sum_l / self.stopped as f64
    }
    fn var(n: f64, s: f64, s2: f64) -> f64 {
        (s2 / n - (s / n).powi(2)) * n / (n - 1.0)
    }
    /// Standard error of the mean projected depth.
    pub(crate) fn se_x(&self) -> f64 {
        let n = self.stopped as f64;
        (Self::var(n, self.sum_x, self.sum_x2) / n).sqrt()
    }
    /// Standard error of the mean path.
    pub(crate) fn se_l(&self) -> f64 {
        let n = self.stopped as f64;
        (Self::var(n, self.sum_l, self.sum_l2) / n).sqrt()
    }
    /// Standard error of the ratio of means `mean_x / mean_l` (delta method,
    /// with the covariance).
    pub(crate) fn se_ratio(&self) -> f64 {
        let n = self.stopped as f64;
        let (mx, ml) = (self.mean_x(), self.mean_l());
        let vx = Self::var(n, self.sum_x, self.sum_x2) / n;
        let vl = Self::var(n, self.sum_l, self.sum_l2) / n;
        let cxl = (self.sum_xl / n - mx * ml) * n / (n - 1.0) / n;
        let r = mx / ml;
        (r * r * (vx / (mx * mx) + vl / (ml * ml) - 2.0 * cxl / (mx * ml))).sqrt()
    }
}

impl BcaTally for RangeTally {
    fn begin_history(&mut self, _index: u64) {
        self.last = [0.0; 3];
        self.path = 0.0;
    }
    fn lattice(
        &mut self,
        at: [f64; 3],
        _layer: usize,
        _kind: lindhard::ion::bca::LatticeDeposit,
        _e: f64,
    ) {
        self.vertex(at);
    }
    fn recoil(&mut self, r: &Particle) {
        self.vertex(r.pos);
    }
    fn stopped(&mut self, p: &Particle) {
        if p.is_primary() {
            self.vertex(p.pos);
            let (x, l) = (p.pos[0], self.path);
            self.stopped += 1;
            self.sum_x += x;
            self.sum_x2 += x * x;
            self.sum_l += l;
            self.sum_l2 += l * l;
            self.sum_xl += x * l;
        }
    }
    fn escaped(&mut self, p: &Particle, face: Face) {
        if p.is_primary() && face == Face::Front {
            self.backscattered += 1;
        }
    }
    fn merge(&mut self, o: Self) {
        self.stopped += o.stopped;
        self.backscattered += o.backscattered;
        self.sum_x += o.sum_x;
        self.sum_x2 += o.sum_x2;
        self.sum_l += o.sum_l;
        self.sum_l2 += o.sum_l2;
        self.sum_xl += o.sum_xl;
    }
}

/// No electronic stopping (nuclear-only runs).
struct NoElectronic;

impl ElectronicStopping for NoElectronic {
    fn name(&self) -> &'static str {
        "none"
    }
    fn stopping(&self, _ion: &Ion, _z: u8, _e: f64) -> Result<f64, StoppingError> {
        Ok(0.0)
    }
    fn validity(&self, _ion: &Ion) -> ValidityRange {
        ValidityRange {
            min_energy_ev: 0.0,
            max_energy_ev: f64::INFINITY,
        }
    }
}

/// Amorphous Si at its tabulated density. `E_d` = 15 eV is a test
/// parameter (the engine requires one; primary-only runs do not depend on it).
pub(crate) fn silicon() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m
}

const CUTOFF_EV: f64 = 5.0;

fn run_engine(
    table: &ScatteringTable,
    stopping: &(dyn ElectronicStopping + Sync),
    z1: u8,
    e0: f64,
    count: u64,
    seed: u64,
    weak: u8,
) -> RangeTally {
    let stack = Stack::semi_infinite(silicon());
    let mut cfg = BcaConfig::new(CUTOFF_EV, 1.0);
    cfg.follow_recoils = false;
    cfg.seed = seed;
    cfg.weak_collisions = weak;
    let beam = Beam::normal(Ion::new(z1).unwrap(), e0, count);
    Bca::new(beam, &stack, cfg, stopping, table)
        .unwrap()
        .run(RangeTally::default)
        .unwrap()
}

fn weak_label(weak: u8) -> String {
    if weak == 0 {
        String::new()
    } else {
        format!(" (K = {weak} weak collisions)")
    }
}

/// Systematic allowance for the engine-vs-equation comparison, relative: the
/// constant (not exponential) flight length and the random first flight move
/// the end of the path by a fraction of `l = 0.27 nm` in Si, and the solver's
/// discretisation uncertainty is added separately. 0.5 % of a 30-70 nm range
/// is 0.15-0.35 nm.
const SYSTEMATIC: f64 = 0.005;

/// Range checks: engine vs the first-moment equation.
pub(crate) fn checks(table: &ScatteringTable, quick: bool) -> Vec<Check> {
    let si = silicon();
    let n = si.atom_number_density();
    let m2 = lindhard::elements::element(14).unwrap().atomic_weight;
    let ls = LindhardScharff::new();
    let count = if quick { 10_000 } else { 100_000 };
    let mut out = Vec::new();

    // Electronic + nuclear: P and As at 50 keV in Si (backscatter < 1 %).
    // The last case runs the engine with K = 3 weak collisions and the
    // equation with p out to p_max sqrt(4) (issue #64).
    for (id_rp, id_ratio, sym, z1, e0, weak) in [
        (
            "range.p50k_si.rp",
            "range.p50k_si.rp_over_l",
            "P",
            15u8,
            50e3,
            0u8,
        ),
        (
            "range.as50k_si.rp",
            "range.as50k_si.rp_over_l",
            "As",
            33,
            50e3,
            0,
        ),
        (
            "range.p50k_si_weak3.rp",
            "range.p50k_si_weak3.rp_over_l",
            "P",
            15,
            50e3,
            3,
        ),
    ] {
        let ion = Ion::new(z1).unwrap();
        let se = |e: f64| ls.stopping(&ion, 14, e).unwrap() / J_PER_EV;
        let pr = Problem {
            ion,
            z2: 14,
            m2,
            n,
            cutoff_ev: CUTOFF_EV,
            se: &se,
            table,
            weak,
        };
        let (rp_lss, drp) = solve(&pr, e0, Moment::Projected);
        let (l_lss, dl) = solve(&pr, e0, Moment::Path);
        let t = run_engine(table, &ls, z1, e0, count, 0x5EED + u64::from(z1), weak);
        let bs = t.backscattered as f64 / count as f64;
        let dev = t.mean_x() / rp_lss - 1.0;
        let tol = 3.0 * t.se_x() / rp_lss + drp / rp_lss + SYSTEMATIC;
        out.push(Check::at_most(
            id_rp,
            format!(
                "{sym} {} keV -> Si{}: engine mean projected range {:.2} nm vs LSS first-moment equation {:.2} nm (|rel. dev.|)",
                e0 / 1e3,
                weak_label(weak),
                t.mean_x() / NM,
                rp_lss / NM
            ),
            dev.abs(),
            tol,
            pct,
            format!(
                "ZBL + LS, primary only, {count} ions; tol = 3 SE + solver + 0.5 % systematic; backscatter {:.2} % excluded from engine mean",
                100.0 * bs
            ),
        ));
        if z1 == 15 && weak == 0 {
            // Power of the check: the same comparison against the equation
            // solved with a wrong deflection must fail by a wide margin.
            let (rp_bug, _) = solve(&pr, e0, Moment::ProjectedCmAngleBug);
            let dev_bug = (t.mean_x() / rp_bug - 1.0).abs();
            out.push(Check::holds(
                "range.p50k_si.control",
                format!(
                    "Negative control: equation with the CM angle used as the lab angle gives {:.2} nm; engine deviates from it by more than the tolerance above",
                    rp_bug / NM
                ),
                dev_bug > tol,
                spct(t.mean_x() / rp_bug - 1.0),
                "shows the Rp check can detect a deflection bug, which an energy-bookkeeping (CSDA) check cannot",
            ));
        }
        let r_eng = t.mean_x() / t.mean_l();
        let r_lss = rp_lss / l_lss;
        let tol_r = 3.0 * t.se_ratio() + r_lss * (drp / rp_lss + dl / l_lss) + SYSTEMATIC * r_lss;
        out.push(Check::at_most(
            id_ratio,
            format!(
                "{sym} {} keV -> Si{}: Rp / path ratio, engine {} vs equation {} (|abs. diff.|)",
                e0 / 1e3,
                weak_label(weak),
                num(r_eng),
                num(r_lss)
            ),
            (r_eng - r_lss).abs(),
            tol_r,
            num,
            format!(
                "mean path {:.2} nm (engine) vs {:.2} nm (equation)",
                t.mean_l() / NM,
                l_lss / NM
            ),
        ));
    }

    // Nuclear stopping only: the reduced range rho(eps) regime. Si -> Si
    // (M1 = M2) at 5 keV, eps ~ 0.1; again with K = 3 weak collisions.
    for (id, weak) in [
        ("range.si5k_si.nuclear_rho", 0u8),
        ("range.si5k_si_weak3.nuclear_rho", 3),
    ] {
        let ion = Ion::new(14).unwrap();
        let zero = |_e: f64| 0.0;
        let pr = Problem {
            ion,
            z2: 14,
            m2,
            n,
            cutoff_ev: CUTOFF_EV,
            se: &zero,
            table,
            weak,
        };
        let e0 = 5e3;
        let (l_lss, dl) = solve(&pr, e0, Moment::Path);
        let (rp_lss, drp) = solve(&pr, e0, Moment::Projected);
        let t = run_engine(table, &NoElectronic, 14, e0, count, 0xA11CE, weak);
        // Reduced path rho = N L pi a^2 4 M1 M2 / (M1 + M2)^2 (LSS 1963).
        let pot = Potential::new(Screening::ZblUniversal, 14.0, 14.0);
        let a = pot.screening_length();
        let m1 = ion.mass_amu();
        let to_rho = n * PI * a * a * 4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2));
        let eps0 = pot.reduced_energy(Potential::cm_energy(e0 * J_PER_EV, m1, m2));
        let dev = t.mean_l() / l_lss - 1.0;
        out.push(Check::at_most(
            id,
            format!(
                "Si 5 keV -> Si{}, nuclear stopping only (eps = {eps0:.3}): engine reduced path rho {:.3} vs LSS first-moment equation {:.3} (|rel. dev.|)",
                weak_label(weak),
                t.mean_l() * to_rho,
                l_lss * to_rho
            ),
            dev.abs(),
            3.0 * t.se_l() / l_lss + dl / l_lss + SYSTEMATIC,
            pct,
            format!(
                "{count} ions; backscatter {:.1} % excluded from engine mean (path is less sensitive to it than Rp)",
                100.0 * t.backscattered as f64 / count as f64
            ),
        ));
        if weak > 0 {
            continue;
        }
        let csda = csda_path(&pr, e0);
        out.push(Check::info(
            "range.si5k_si.csda",
            "Si 5 keV -> Si, nuclear only: LSS mean path / continuous-slowing-down path",
            spct(l_lss / csda - 1.0),
            "context: the CSDA integral ignores the fluctuation of the nuclear energy loss, so it is not the mean path; not a pass/fail quantity",
        ));
        out.push(Check::info(
            "range.si5k_si.rp",
            format!(
                "Si 5 keV -> Si, nuclear only: engine Rp {:.2} nm vs equation {:.2} nm",
                t.mean_x() / NM,
                rp_lss / NM
            ),
            spct(t.mean_x() / rp_lss - 1.0),
            format!(
                "not asserted: the {:.1} % of ions that backscatter (equal masses) leave the engine's semi-infinite target but stay in the equation; solver uncertainty {:.3} %",
                100.0 * t.backscattered as f64 / count as f64,
                100.0 * drp / rp_lss
            ),
        ));
    }

    // B 10 keV: the case the PR #46 review solved independently (44.4 nm in
    // an infinite medium). Reported, because B backscatters several percent
    // from Si and the engine's semi-infinite mean drops those ions.
    {
        let ion = Ion::new(5).unwrap();
        let se = |e: f64| ls.stopping(&ion, 14, e).unwrap() / J_PER_EV;
        let pr = Problem {
            ion,
            z2: 14,
            m2,
            n,
            cutoff_ev: CUTOFF_EV,
            se: &se,
            table,
            weak: 0,
        };
        let (rp_lss, _) = solve(&pr, 10e3, Moment::Projected);
        let t = run_engine(table, &ls, 5, 10e3, count, 0xB0B, 0);
        out.push(Check::info(
            "range.b10k_si.rp",
            format!(
                "B 10 keV -> Si: engine Rp {:.2} nm vs equation {:.2} nm",
                t.mean_x() / NM,
                rp_lss / NM
            ),
            spct(t.mean_x() / rp_lss - 1.0),
            format!(
                "not asserted: {:.1} % backscatter excluded from the engine mean (semi-infinite target) but kept by the infinite-medium equation",
                100.0 * t.backscattered as f64 / count as f64
            ),
        ));
    }
    out
}
