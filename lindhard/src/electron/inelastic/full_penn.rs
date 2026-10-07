//! The full Penn algorithm (FPA): the optical energy-loss function expanded
//! over the Lindhard loss functions of free-electron gases, and the DIIMFP,
//! inelastic mean free path and stopping power built from it.
//!
//! # Sources
//!
//! The algorithm is Penn's, D. R. Penn, Phys. Rev. B 35, 482 (1987),
//! doi:10.1103/PhysRevB.35.482, which is closed access and **was not read**.
//! As for the single-pole model ([`super::penn`]), the equations are those of
//! the account by Penn and co-authors in H. Shinotsuka, B. Da, S. Tanuma, H.
//! Yoshikawa, C. J. Powell, D. R. Penn, Surf. Interface Anal. 49, 238 (2017),
//! doi:10.1002/sia.6123 (NIST author manuscript, PMC5524379, "S2017"), section
//! "Full Penn algorithm (FPA)", eqs. (4)-(6):
//!
//! ```text
//! Im[-1/ε(q, ω)] = ∫_0^∞ g(ω_p) Im[-1/ε_L(q, ω; ω_p)] dω_p,   (S2017 eq. (4))
//! g(ω_p) = (2/(π ω_p)) Im[-1/ε(ω_p)],                         (S2017 eq. (5))
//! ```
//!
//! where `ε_L` is the Lindhard dielectric function of a free-electron gas of
//! plasma frequency `ω_p` ([`super::lindhard_gas`]; `ω_p² = 4π n`), and the
//! optical ELF `Im[-1/ε(ω)]` is the user's table. S2017 eq. (6) splits the
//! result into a plasmon-pole part and a single-electron part,
//!
//! ```text
//! Im[-1/ε(q, ω)] = Im[-1/ε(q, ω)]_pl + Im[-1/ε(q, ω)]_se.
//! ```
//!
//! S2017 defers "the details of the calculation procedures for each term" to
//! its reference 22 (Shinotsuka et al., Surf. Interface Anal. 47, 871
//! (2015)), which was not read either. The procedure used here is our own, a
//! direct evaluation of eq. (4): for every `ω_p` the Lindhard loss function
//! is its electron-hole continuum `Im ε/|ε|²` plus the undamped plasmon
//! delta function of weight `π/(∂Re ε/∂ω)` (exactly the pole of the SPA,
//! S2017 eq. (7), at `q -> 0`). Integrating the delta function over `ω_p`
//! gives a smooth plasmon term `g(ω₀) W/|dω_pl/dω_p|` at the `ω_p = ω₀`
//! whose plasmon sits at the loss `ω`. Below the plasmon cutoff the plasmon
//! term is the analogue of the SPA's S2017 eq. (13); above it, and at all
//! `q` outside the plasmon, the continuum term is the "single-electron
//! excitations" that the SPA neglects. The two are the `pl` and `se` of
//! eq. (6) as far as we can tell from the account in S2017; nothing here was
//! checked against Penn's own formulae.
//!
//! The kinematics, units and the definitions of the DIIMFP, IMFP and
//! stopping power are those of the single-pole module (S2017 eqs. (2)-(3),
//! nonrelativistic; Hartree units internally, eV and metres in the API).
//!
//! # Sum rules
//!
//! Each Lindhard gas satisfies `∫ ω Im[-1/ε_L] dω = (π/2) ω_p²` at every `q`
//! (plasmon included; tested in [`super::lindhard_gas`]), so with eq. (5) the
//! expanded ELF has `∫ ω Im[-1/ε(q, ω)] dω = ∫ ω_p Im[-1/ε(ω_p)] dω_p`, the
//! optical f-sum, at every `q` ([`FullPenn::f_sum_ev2_at`] and the
//! integration tests check this).
//!
//! # How the integrals are done
//!
//! **IMFP and stopping power** ([`FullPenn::imfp_and_stopping`]) are
//! evaluated with the order of integration `ω_p`, `q`, `ω`:
//!
//! ```text
//! λ⁻¹ = (1/(π T')) ∫ dω_p g(ω_p) ∫ (dq/q) ∫_0^{ω_c(q)} Im[-1/ε_L] dω,
//! S   = (1/(π T')) ∫ dω_p g(ω_p) ∫ (dq/q) ∫_0^{ω_c(q)} ω Im[-1/ε_L] dω,
//! ω_c(q) = min(T, T' - (q - k)²/2),   k = sqrt(2T'),   0 < q < 2k,
//! ```
//!
//! (`q- <= q <= q+` is `ω <= T' - (q - k)²/2`, S2017 eq. (2)). The innermost
//! integral is split at the kinks of the Lindhard continuum and includes the
//! plasmon delta function if its energy is below `ω_c(q)`; the `q` integral
//! is in `ln q`, split at the kinematic landmarks and at the momenta where
//! the plasmon enters (`ω_pl = ω_c`) and leaves (cutoff) the allowed region;
//! the outer integral is split at every knot of the
//! table (the ELF is piecewise linear, and a feature narrower than the
//! spacing of a subset of knots could be missed entirely by the adaptive
//! quadrature, so no knot is dropped). All three use the adaptive Gauss-Legendre quadrature of
//! `super::quadrature`, each with the tolerance of the model
//! ([`DEFAULT_FULL_TOLERANCE`] = 1e-4): the result has a relative error of a
//! few times that. The DIIMFP is the same integrals in the order `q`, `ω_p`
//! at fixed `ω`.
//!
//! No Fermi-energy, exchange or relativistic correction is included, as for
//! the single-pole model.
//!
//! # Comparison with the single-pole model (synthetic Drude fixture)
//!
//! For the Drude plasmon of the tests and examples (`E_p = 20 eV`, `γ = 5 eV`,
//! synthetic, tabulated at 1000 energies from 10 meV to 100 keV; default
//! tolerances; `cargo run --release -p lindhard --example
//! penn_full_vs_single_pole`), ratios full/single-pole of the stopping power
//! `S` and of the IMFP `λ`:
//!
//! ```text
//!    E / eV     S ratio    λ ratio
//!      20.0      3.7971     0.2830
//!      30.0      1.6561     0.5489
//!      50.0      0.7988     1.1365
//!     100.0      0.9605     1.0197
//!     200.0      0.9860     1.0067
//!     500.0      0.9953     1.0028
//!    1000.0      0.9977     1.0018
//!   10000.0      0.9998     1.0010
//!   50000.0      1.0000     1.0008
//! ```
//!
//! **Above 10 keV the two agree to 1e-3 in both quantities** (inside the 1 %
//! requirement of the issue); both reach the Bethe limit, which fixes `S`
//! through the f-sum and the mean excitation energy of the optical ELF, and
//! their IMFP constants differ only by the single-electron continuum, a
//! 1e-3 effect at these energies. At **100 eV to 1 keV** the differences are
//! at most 4 % (`S`) and 2 % (`λ`). **Below about 100 eV they are large and
//! of both signs**: at 20 to 30 eV the full model's stopping power and
//! inverse IMFP are several times those of the single pole (the full model
//! has single-electron excitations that the single pole does not, which is
//! presumably the cause; it was not isolated), while at 50 eV the single pole is the larger. That the
//! full-algorithm IMFP is smaller than the single-pole one at low energies
//! is the direction reported for water in S2017 (section on Fig. 2); the
//! sign change near 50 eV here depends on the shape of the synthetic ELF,
//! and no statement about a real material follows. Nothing was tuned to
//! make the two agree: these are the numbers of the implementation as it
//! stands, and the regression bounds of `tests/penn_full.rs` follow them.

use super::lindhard_gas::{kf_coefficient, LindhardGas};
use super::penn::InelasticPoint;
use super::quadrature::{integrate_segments, GaussLegendre};
use crate::constants::{BOHR_RADIUS, ELEMENTARY_CHARGE, HARTREE_ENERGY};
use crate::electron::data::{ElectronDataError, OpticalElf};
use std::f64::consts::PI;

type Result<T> = std::result::Result<T, ElectronDataError>;

/// Default relative tolerance of the integrals of the full Penn model.
pub const DEFAULT_FULL_TOLERANCE: f64 = 1.0e-4;

const GL_ORDER: usize = 5;

fn hartree_ev() -> f64 {
    HARTREE_ENERGY / ELEMENTARY_CHARGE
}

fn check_energy(what: &'static str, v: f64) -> Result<()> {
    if v.is_finite() && v > 0.0 {
        Ok(())
    } else {
        Err(ElectronDataError::Invalid {
            what,
            reason: format!("must be finite and positive, got {v} eV"),
        })
    }
}

/// Geometric bisection for the point where `pred` switches, on `[lo, hi]`,
/// with `pred(lo) != pred(hi)`. Returns the end where `pred` equals
/// `pred(hi)`.
fn bisect_switch(mut lo: f64, mut hi: f64, pred: &impl Fn(f64) -> bool) -> f64 {
    let at_hi = pred(hi);
    for _ in 0..200 {
        let mid = if lo > 0.0 {
            (lo * hi).sqrt()
        } else {
            0.5 * (lo + hi)
        };
        if !(mid > lo && mid < hi) {
            break;
        }
        if pred(mid) == at_hi {
            hi = mid;
        } else {
            lo = mid;
        }
        if hi - lo <= 1e-14 * hi {
            break;
        }
    }
    hi
}

/// The full Penn model built from one optical ELF. See the module docs.
#[derive(Debug, Clone)]
pub struct FullPenn {
    elf: OpticalElf,
    /// Tabulated energies, Hartree.
    w: Vec<f64>,
    /// Tabulated ELF values.
    e: Vec<f64>,
    /// The break points of the outer integral (every knot of `w`).
    outer: Vec<f64>,
    fermi: f64,
    rel_tol: f64,
    gl: GaussLegendre,
}

impl FullPenn {
    /// The model for `elf`, with Fermi energy zero and the default tolerance.
    pub fn new(elf: OpticalElf) -> Self {
        let h = hartree_ev();
        let w: Vec<f64> = elf.energy_ev().iter().map(|x| x / h).collect();
        let e = elf.elf_values().to_vec();
        let outer: Vec<f64> = w.clone();
        Self {
            elf,
            w,
            e,
            outer,
            fermi: 0.0,
            rel_tol: DEFAULT_FULL_TOLERANCE,
            gl: GaussLegendre::new(GL_ORDER),
        }
    }

    /// Set the Fermi energy (eV, finite and non-negative), as for the
    /// single-pole model: energies are kinetic energies above the Fermi
    /// level and the kinematics use `T' = T + E_F`.
    pub fn with_fermi_energy_ev(mut self, fermi_ev: f64) -> Result<Self> {
        if !(fermi_ev.is_finite() && fermi_ev >= 0.0) {
            return Err(ElectronDataError::Invalid {
                what: "Fermi energy",
                reason: format!("must be finite and non-negative, got {fermi_ev} eV"),
            });
        }
        self.fermi = fermi_ev / hartree_ev();
        Ok(self)
    }

    /// Set the relative integration tolerance, in `[1e-9, 1e-2]`.
    pub fn with_relative_tolerance(mut self, rel_tol: f64) -> Result<Self> {
        if !(1e-9..=1e-2).contains(&rel_tol) {
            return Err(ElectronDataError::Invalid {
                what: "integration tolerance",
                reason: format!("must be in [1e-9, 1e-2], got {rel_tol}"),
            });
        }
        self.rel_tol = rel_tol;
        Ok(self)
    }

    /// The optical ELF the model is built from.
    pub fn optical_elf(&self) -> &OpticalElf {
        &self.elf
    }

    /// The Fermi energy, eV.
    pub fn fermi_energy_ev(&self) -> f64 {
        self.fermi * hartree_ev()
    }

    /// The relative integration tolerance.
    pub fn relative_tolerance(&self) -> f64 {
        self.rel_tol
    }

    /// The optical ELF at `x` Hartree (linear interpolation, zero outside
    /// the table).
    fn elf_au(&self, x: f64) -> f64 {
        let n = self.w.len();
        if !(x >= self.w[0] && x <= self.w[n - 1]) {
            return 0.0;
        }
        let i = self.w.partition_point(|&g| g <= x).clamp(1, n - 1);
        let (x0, x1) = (self.w[i - 1], self.w[i]);
        let t = (x - x0) / (x1 - x0);
        (self.e[i - 1] * (1.0 - t) + self.e[i] * t).max(0.0)
    }

    /// `g(ω_p)` of S2017 eq. (5).
    fn g_au(&self, x: f64) -> f64 {
        2.0 * self.elf_au(x) / (PI * x)
    }

    fn table_range(&self) -> (f64, f64) {
        (self.w[0], self.w[self.w.len() - 1])
    }

    /// `∫ dω_p g(ω_p) [∫ ω^j Im[-1/ε_L(q, ω)] dω over (0, ∞)]` is the
    /// optical f-sum at `j = 1`: the `ω` integral of `ω Im[-1/ε(q, ω)]` at
    /// momentum `q` (m⁻¹), eV², evaluated Lindhard gas by Lindhard gas.
    /// The integral over the Lindhard continuum is cut off at
    /// `ω_max = 40 (q k_F + q²/2 + ω_p)` for each gas. The optical f-sum
    /// `∫ ω Im[-1/ε(ω)] dω` of the same table (the value this must equal,
    /// [`super::sum_rules::SumRuleReport::f_sum_ev2`]) is independent of `q`.
    pub fn f_sum_ev2_at(&self, q_per_m: f64) -> Result<f64> {
        if !(q_per_m.is_finite() && q_per_m > 0.0) {
            return Err(ElectronDataError::Invalid {
                what: "momentum transfer",
                reason: format!("must be finite and positive, got {q_per_m} m^-1"),
            });
        }
        let q = q_per_m * BOHR_RADIUS;
        let breaks = &self.outer;
        let r = integrate_segments(
            &self.gl,
            &mut |x: f64| {
                let weight = self.g_au(x);
                if weight == 0.0 {
                    return [0.0];
                }
                let gas = LindhardGas::from_plasma_frequency(x);
                let wcut = 40.0 * (gas.continuum_upper_edge(q) + x);
                [weight * gas.loss_moments(&self.gl, q, wcut, self.rel_tol)[1]]
            },
            breaks,
            self.rel_tol,
        )[0];
        let h = hartree_ev();
        Ok(r * h * h)
    }

    // ---- the expanded loss function at fixed (q, ω) ----

    /// The smooth plasmon term of `Im[-1/ε(q, ω)]`:
    /// `g(ω₀) W(q; ω₀)/(dω_pl/dω_p)` at the `ω₀` whose plasmon is at `ω`.
    fn plasmon_term(&self, q: f64, w: f64) -> f64 {
        let (wmin, wmax) = self.table_range();
        let gas_at = |x: f64| LindhardGas::from_plasma_frequency(x);
        let energy = |x: f64| gas_at(x).plasmon(q).map(|p| p.energy);
        // the plasmon exists for ω_p above a threshold x_c(q)
        let exists = |x: f64| energy(x).is_some();
        if !exists(wmax) {
            return 0.0;
        }
        let x_c = if exists(wmin) {
            wmin
        } else {
            bisect_switch(wmin, wmax, &exists)
        };
        let e_c = energy(x_c).unwrap_or(f64::INFINITY);
        let e_top = energy(wmax).unwrap_or(0.0);
        if w < e_c || w > e_top {
            return 0.0;
        }
        // ω_pl(ω_p) increases with ω_p
        let reached = |x: f64| energy(x).is_some_and(|e| e >= w);
        let x0 = if reached(x_c) {
            x_c
        } else {
            bisect_switch(x_c, wmax, &reached)
        };
        let weight = self.g_au(x0);
        if weight == 0.0 {
            return 0.0;
        }
        let Some(p0) = gas_at(x0).plasmon(q) else {
            return 0.0;
        };
        let h = 1e-6 * x0;
        let (lo, hi) = (x0 - h, x0 + h);
        let slope = match (energy(lo), energy(hi)) {
            (Some(a), Some(b)) => (b - a) / (2.0 * h),
            (_, Some(b)) => (b - p0.energy) / h,
            (Some(a), None) => (p0.energy - a) / h,
            _ => return 0.0,
        };
        if slope <= 0.0 {
            return 0.0;
        }
        weight * p0.weight / slope
    }

    /// The expanded loss function `Im[-1/ε(q, ω)]` in atomic units.
    fn loss_au(&self, q: f64, w: f64) -> f64 {
        let (wmin, wmax) = self.table_range();
        let a = w / q - 0.5 * q;
        let x_a = (a.abs() / kf_coefficient()).powf(1.5).max(wmin);
        let mut cont = 0.0;
        if x_a < wmax {
            let mut breaks = vec![x_a];
            if a + q > 0.0 {
                let x = ((a + q) / kf_coefficient()).powf(1.5);
                if x > x_a && x < wmax {
                    breaks.push(x);
                }
            }
            for &x in &self.outer {
                if x > x_a && x < wmax {
                    breaks.push(x);
                }
            }
            breaks.push(wmax);
            breaks.sort_by(f64::total_cmp);
            breaks.dedup();
            cont = integrate_segments(
                &self.gl,
                &mut |x: f64| {
                    let weight = self.g_au(x);
                    if weight == 0.0 {
                        [0.0]
                    } else {
                        [weight * LindhardGas::from_plasma_frequency(x).loss_continuum(q, w)]
                    }
                },
                &breaks,
                self.rel_tol,
            )[0];
        }
        cont + self.plasmon_term(q, w)
    }

    /// `Im[-1/ε(q, ω)]` of the full Penn algorithm (S2017 eq. (4)), at
    /// momentum transfer `q` (m⁻¹) and energy loss `ω` (eV). At `q = 0` it is
    /// the optical ELF.
    pub fn loss_function(&self, q_per_m: f64, loss_ev: f64) -> f64 {
        if !(q_per_m.is_finite() && q_per_m >= 0.0 && loss_ev.is_finite() && loss_ev > 0.0) {
            return 0.0;
        }
        let w = loss_ev / hartree_ev();
        if q_per_m == 0.0 {
            return self.elf_au(w);
        }
        self.loss_au(q_per_m * BOHR_RADIUS, w)
    }

    /// The DIIMFP `p(T, ω)` (S2017 eq. (2), nonrelativistic) at kinetic
    /// energy `energy_ev` above the Fermi level and energy loss `loss_ev`, in
    /// m⁻¹ eV⁻¹. Zero for `ω <= 0` or `ω > T`.
    pub fn diimfp_per_m_ev(&self, energy_ev: f64, loss_ev: f64) -> Result<f64> {
        check_energy("electron energy", energy_ev)?;
        if !loss_ev.is_finite() {
            return Err(ElectronDataError::Invalid {
                what: "energy loss",
                reason: format!("must be finite, got {loss_ev} eV"),
            });
        }
        let h = hartree_ev();
        let p_au = self.diimfp_au(energy_ev / h, loss_ev / h);
        Ok(p_au / (BOHR_RADIUS * h))
    }

    fn diimfp_au(&self, t: f64, w: f64) -> f64 {
        if !(w > 0.0 && w <= t) {
            return 0.0;
        }
        let tp = t + self.fermi;
        let k = (2.0 * tp).sqrt();
        let s = (2.0 * (tp - w)).max(0.0).sqrt();
        let (q_minus, q_plus) = (2.0 * w / (k + s), k + s);
        if q_plus.is_nan() || q_plus <= q_minus {
            return 0.0;
        }
        let (lo, hi) = (q_minus.ln(), q_plus.ln());
        let breaks: Vec<f64> = (0..=8).map(|i| lo + (hi - lo) * i as f64 / 8.0).collect();
        let integral = integrate_segments(
            &self.gl,
            &mut |u: f64| [self.loss_au(u.exp(), w)],
            &breaks,
            self.rel_tol,
        )[0];
        integral / (PI * tp)
    }

    // ---- IMFP and stopping power ----

    /// `[∫ dq/q ∫_0^{ω_c} Im[-1/ε_L] dω, same with an extra ω]` for the gas
    /// of plasma frequency `x`, at kinetic energy `t` (Hartree).
    fn gas_moments(&self, x: f64, t: f64) -> [f64; 2] {
        let tp = t + self.fermi;
        let k = (2.0 * tp).sqrt();
        let gas = LindhardGas::from_plasma_frequency(x);
        let kf = gas.fermi_wavenumber();
        let wcut = |q: f64| t.min(tp - 0.5 * (q - k) * (q - k));
        let q_lo = 1e-3 * (x / k).min(k);
        let q_hi = 2.0 * k;
        if q_lo.is_nan() || q_lo >= q_hi {
            return [0.0; 2];
        }
        let mut marks = vec![2.0 * kf, k - kf, k + kf, k];
        // plasmon landmarks: cutoff, entry into and exit from the allowed region
        let qc = gas.plasmon_cutoff();
        marks.push(qc);
        let inc = |q: f64| gas.plasmon(q).is_some_and(|p| p.energy < wcut(q));
        let q_top = qc.min(k) * (1.0 - 1e-9);
        if q_lo < q_top && inc(q_top) {
            marks.push(if inc(q_lo) {
                q_lo
            } else {
                bisect_switch(q_lo, q_top, &inc)
            });
        }
        let q_out = qc * (1.0 - 1e-9);
        if qc > k && q_out < q_hi && !inc(q_out) && inc(k) {
            marks.push(bisect_switch(k, q_out, &inc));
        }
        let mut breaks = vec![q_lo.ln(), q_hi.ln()];
        for m in marks {
            if m > q_lo && m < q_hi {
                breaks.push(m.ln());
            }
        }
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        integrate_segments(
            &self.gl,
            &mut |u: f64| {
                let q = u.exp();
                gas.loss_moments(&self.gl, q, wcut(q), self.rel_tol)
            },
            &breaks,
            self.rel_tol,
        )
    }

    /// `(λ⁻¹, S)` at kinetic energy `t` (Hartree) above the Fermi level, in
    /// atomic units.
    fn imfp_and_stopping_au(&self, t: f64) -> (f64, f64) {
        if t.is_nan() || t <= 0.0 {
            return (0.0, 0.0);
        }
        let tp = t + self.fermi;
        let r = integrate_segments(
            &self.gl,
            &mut |x: f64| {
                let weight = self.g_au(x);
                if weight == 0.0 {
                    return [0.0, 0.0];
                }
                let m = self.gas_moments(x, t);
                [weight * m[0], weight * m[1]]
            },
            &self.outer,
            self.rel_tol,
        );
        (r[0] / (PI * tp), r[1] / (PI * tp))
    }

    /// The inverse IMFP (m⁻¹) and the stopping power (eV/m) at kinetic
    /// energy `energy_ev` above the Fermi level.
    pub fn imfp_and_stopping(&self, energy_ev: f64) -> Result<InelasticPoint> {
        check_energy("electron energy", energy_ev)?;
        let h = hartree_ev();
        let (inv, s) = self.imfp_and_stopping_au(energy_ev / h);
        Ok(InelasticPoint {
            energy_ev,
            inverse_imfp_per_m: inv / BOHR_RADIUS,
            stopping_ev_per_m: s * h / BOHR_RADIUS,
        })
    }

    /// The inelastic mean free path `λ(E)`, m.
    pub fn imfp_m(&self, energy_ev: f64) -> Result<f64> {
        Ok(self.imfp_and_stopping(energy_ev)?.imfp_m())
    }

    /// The stopping power `S(E)`, eV/m.
    pub fn stopping_power_ev_per_m(&self, energy_ev: f64) -> Result<f64> {
        Ok(self.imfp_and_stopping(energy_ev)?.stopping_ev_per_m)
    }

    /// IMFP and stopping power on a grid of energies, in order. Sequential
    /// and deterministic.
    pub fn tabulate(&self, energies_ev: &[f64]) -> Result<Vec<InelasticPoint>> {
        energies_ev
            .iter()
            .map(|&e| self.imfp_and_stopping(e))
            .collect()
    }
}
