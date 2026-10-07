//! The single-pole Penn algorithm: the momentum-dependent energy-loss
//! function, the differential inverse inelastic mean free path (DIIMFP), the
//! inelastic mean free path (IMFP) and the stopping power, built from an
//! optical ELF.
//!
//! # Sources
//!
//! The algorithm is Penn's, D. R. Penn, Phys. Rev. B 35, 482 (1987),
//! doi:10.1103/PhysRevB.35.482 (the "simple Penn algorithm" or single-pole
//! approximation, SPA). That paper is closed access and **was not read** for
//! this implementation. The equations below are taken from the account of the
//! SPA by Penn and co-authors in H. Shinotsuka, B. Da, S. Tanuma, H.
//! Yoshikawa, C. J. Powell, D. R. Penn, Surf. Interface Anal. 49, 238 (2017),
//! doi:10.1002/sia.6123 (NIST author manuscript, open copy PMC5524379), and
//! every equation number in this module refers to **that** paper ("S2017").
//! The plasmon-pole idea for extending an optical ELF into momentum transfer
//! goes back to Lindhard, Mat. Fys. Medd. Dan. Vid. Selsk. 28, no. 8 (1954)
//! (the dielectric function each pole stands for) and Ritchie and Howie,
//! Philos. Mag. 36, 463 (1977), doi:10.1080/14786437708244948 (neither read
//! here; cited for attribution only).
//!
//! Hartree atomic units (`m_e = e = ħ = 1`) are used internally, as in S2017;
//! the public API is in eV and metres.
//!
//! # The model
//!
//! The optical ELF `Im[-1/ε(ω)]` is written as a superposition of free-electron
//! (Lindhard) loss functions of plasma frequency `ω_p` with weight
//! `g(ω_p) = (2/(π ω_p)) Im[-1/ε(ω_p)]` (S2017 eqs. (4)-(5)), and in the SPA
//! each Lindhard loss function is replaced by a single plasmon pole
//! (S2017 eq. (7)),
//!
//! ```text
//! Im[-1/ε_L(q, ω; ω_p)] ≈ (π/2) (ω_p² / ω_q(ω_p)) δ(ω - ω_q(ω_p)),
//! ```
//!
//! with Penn's dispersion (S2017 eq. (8))
//!
//! ```text
//! ω_q²(ω_p) = ω_p² + (1/3) (k_F(ω_p) q)² + q⁴/4,    k_F(ω_p) = (3π/4)^(1/3) ω_p^(2/3).
//! ```
//!
//! The `ω_p` integration then gives (S2017 eqs. (9)-(13))
//!
//! ```text
//! Im[-1/ε(q, ω)] = Im[-1/ε(ω₀)] / (1 + π q² / (6 k_F(ω₀))),
//! ```
//!
//! where `ω₀(q, ω)` is the plasma frequency whose pole sits at `ω`, the real
//! positive root `x = ω₀^(2/3)` of the cubic (S2017 eq. (11))
//! `x³ + a(q) x² + b(q, ω) = 0`, `a = (π²/48)^(1/3) q²`, `b = q⁴/4 - ω²`. There
//! is a root only for `b < 0`, i.e. `q < sqrt(2ω)`; otherwise the loss
//! function is zero. ([`SinglePolePenn::loss_function`].)
//!
//! **Kinematics are nonrelativistic** throughout: the relativistic factors of
//! S2017 eq. (2) are set to 1 (`c -> ∞`). Relativistic kinematics is a later
//! option. The DIIMFP is (S2017 eq. (2), nonrelativistic)
//!
//! ```text
//! p(T, ω) = (1 / (π T')) ∫_{q-}^{q+} (dq/q) Im[-1/ε(q, ω)],
//! q± = sqrt(2T') ± sqrt(2(T' - ω)),
//! ```
//!
//! with `T` the kinetic energy measured from the Fermi level, `T' = T + E_F`,
//! and losses `0 < ω <= T` (S2017 eq. (3), `ω_max = T' - E_F`). The Fermi
//! energy `E_F` defaults to zero, so `T = T'` is the kinetic energy. The IMFP
//! is `λ = [∫ p dω]⁻¹` (S2017 eq. (3)) and the stopping power
//! `S = ∫ ω p dω`.
//!
//! # Optional exchange correction (Born-Ochkur)
//!
//! [`ExchangeCorrection`] makes the primary and the struck electron
//! indistinguishable. The source is the open account in P. de Vera, S.
//! Taioli, P. E. Trevisanutto, M. Dapor, I. Abril, S. Simonucci and R.
//! Garcia-Molina, "Energy Deposition around Swift Carbon-Ion Tracks in Liquid
//! Water", Int. J. Mol. Sci. 23, 6121 (2022), doi:10.3390/ijms23116121 (open
//! access, PMC9181504), Section 2.4.3, eqs. (32)-(35) ("dV2022" below). It
//! cites Ochkur and Rudge for the factor (their refs. 122-125) and its own
//! earlier work (their ref. 42) for the setup; none of these was read. From
//! dV2022, for a loss `E = W + B_α` that ionises a shell of binding energy
//! `B_α` and emits an electron of kinetic energy `W`:
//!
//! * the exchange term of the DIIMFP is `(1/T) ∫ dk/k F Im[-1/ε(k, E)]`
//!   (eq. (32), last two terms), with the Born-Ochkur factor
//!   `F = -x + x²`, `x = (k²/2m) / (T - W)` (text below eq. (32); `ħ = m = 1`
//!   here). The denominator is `T - W = T - E + B_α`: it contains the
//!   **emitted** energy `W`, not the loss `E`;
//! * its `k` limits are `sqrt(2mT) ± sqrt(2m(T - E))` (eq. (34)), the limits of
//!   the plain first Born DIIMFP at loss `E`;
//! * the emitted energy is limited by indistinguishability to
//!   `W <= (T - B_α)/2` (text below eq. (35)), i.e. the loss to
//!   `E <= (T + B_α)/2`.
//!
//! Here that term is `(1/(π T')) ∫ dq/q F Im[-1/ε(q, ω)]` with `ω = E` and
//! `x = (q²/2) / (T' - ω + B)` ([`SinglePolePenn::diimfp_with_binding_per_m_ev`]),
//! and the loss is limited to `ω <= (T' + B)/2`. `B` is the binding energy of
//! the channel the loss is attributed to (`super::inner_shell`). The plain
//! model ([`SinglePolePenn::diimfp_per_m_ev`], [`SinglePolePenn::imfp_and_stopping`])
//! and the valence channel use `B = 0`, so `x = (q²/2)/(T' - ω)` and
//! `ω <= T'/2`, which is also the `ω_max` of the SSPA* variant of S2017. The SPA
//! loss function vanishes for `q² >= 2ω` (S2017 eq. (11)), and
//! `ω <= (T' + B)/2` gives `T' - ω + B >= ω`, so `0 <= x <= 1` wherever the
//! integrand is nonzero and `-1/4 <= F <= 0` for every channel.
//!
//! **Where this departs from dV2022 eq. (32), deliberately:**
//!
//! * *Direct term.* dV2022 applies a Coulomb-field correction (their ref. 126,
//!   not read) to the direct term only: `T -> T + 2B_α` in its prefactor
//!   `1/(T + 2B_α)` and in its `k` limits (eq. (33)). We do **not** apply it:
//!   the direct term here is the plain first Born (S2017 eq. (2)) DIIMFP with
//!   prefactor `1/T'` and the eq. (34) limits for every channel, which is what
//!   dV2022 eq. (32) reduces to for `B_α = 0`. So for an inner shell only the
//!   exchange term, its denominator and the loss limit follow eq. (32); the
//!   direct term does not. Adding the correction would also break the property
//!   that the channels add up to the plain-model total.
//! * *Valence binding energy.* dV2022 uses an energy-dependent mean binding
//!   energy `B(T)` of the outer shells, fitted for liquid water in their ref.
//!   42 (not read). No such function is available for a general target, so
//!   the valence channel here has `B = 0`: its exchange denominator is
//!   `T' - ω` and its limit `ω <= T'/2`.
//! * *Fermi energy.* dV2022 has none; here `T' = T + E_F` stands for its `T`
//!   (zero by default).
//!
//! The correction is optional and applied only for kinetic energies below the
//! energy given to [`ExchangeCorrection::new`]; at and above it the model is the
//! plain SPA (no `ω_max` cut either), so `λ(E)` has a step of the size of the
//! correction at that energy. The total is then evaluated by direct numerical
//! integration over `ω` of the DIIMFP (no closed form), with the outer
//! tolerance `100 × rel_tol` since each DIIMFP value is itself an integral.
//!
//! # How the integrals are done
//!
//! **DIIMFP** ([`SinglePolePenn::diimfp_per_m_ev`]): the `q` integral of
//! S2017 eq. (2) is done in `ln q` over `[q-, q+]`, restricted to the `q`
//! where `ω₀(q, ω)` lies inside the tabulated range of the optical ELF (the
//! ELF is zero outside it: nothing is extrapolated), and split at the `q`
//! where `ω₀` crosses a table knot, so that each piece is smooth. `ω₀` is the
//! root of the S2017 eq. (11) cubic, found by Newton's method from
//! `x = ω^(2/3)` (the cubic is convex and increasing for `x > 0`, so the
//! iteration converges monotonically).
//!
//! **IMFP and stopping power** ([`SinglePolePenn::imfp_and_stopping`]): the
//! double integral over `(ω, q)` is done in the other order. Because the SPA
//! puts the whole weight of each `ω_p` on one curve `ω = ω_q(ω_p)`, inserting
//! S2017 eqs. (4), (5), (7) into `∫∫ dω (dq/q)` gives exactly
//!
//! ```text
//! λ⁻¹ = (1/(π T')) ∫ dω_p Im[-1/ε(ω_p)] ω_p ∫_{qa}^{qb} dq / (q ω_q(ω_p)),
//! S   = (1/(π T')) ∫ dω_p Im[-1/ε(ω_p)] ω_p ln(qb / qa),
//! ```
//!
//! where `[qa, qb]` is the set of `q` whose pole `ω = ω_q(ω_p)` is allowed:
//! `q- <= q <= q+` is equivalent to `ω <= q k - q²/2` (`k = sqrt(2T')`), and
//! `ω <= T` to `q <= q_T` with `ω_q(q_T) = T`. The function
//! `φ(q) = q k - q²/2 - ω_q(q)` is concave (the numerator of `d²ω_q/dq²` is
//! `c ω_p² + (3/2) q² ω_p² + (3/4) c q⁴ + q⁶/8 > 0`, `c = k_F²/3`), so the
//! allowed set is one interval; its ends are found by Newton's method from
//! `q = 0` and `q = k`, each converging monotonically. The inner `q`
//! integral of `λ⁻¹` is elementary (`y = q²`,
//! `∫ dy/(y sqrt(y²/4 + c y + ω_p²)) = -(1/ω_p) ln((2ω_p² + c y + 2 ω_p ω_q)/y)`).
//! This reduction (the order of integration and the closed-form inner
//! integral) is our own algebra on the published model; the tests check it
//! against the direct `∫ p(T, ω) dω` to 1e-5.
//!
//! The outer integrals use adaptive Gauss-Legendre quadrature split at the
//! table knots (`super::quadrature`): the estimated relative error of each
//! integral is at most twice the tolerance (default
//! [`DEFAULT_RELATIVE_TOLERANCE`] = 1e-7), plus the error of the Newton
//! roots (relative 1e-15).
//!
//! # High-energy (Bethe) limits
//!
//! As `T -> ∞`, `qa -> ω_p/k` and `qb -> k`. The stopping power tends to
//! `S -> (Ω_p²/(2T)) ln(2T/I)`, with `Ω_p² = (2/π) ∫ ω Im[-1/ε] dω` and
//! `ln I = ∫ ω Im[-1/ε] ln ω dω / ∫ ω Im[-1/ε] dω`: the nonrelativistic Bethe
//! formula with `W_max = T` (see [`super::bethe`]). The inverse IMFP tends to
//!
//! ```text
//! λ⁻¹ -> (1/(2π T')) ∫ dω_p Im[-1/ε(ω_p)] ln(8 T' / (ω_p + k_F(ω_p)²/3)),
//! ```
//!
//! which is the Bethe form `T/λ = A ln(B T)` of S2017 eq. (27) (with
//! `α(T) = 1`): its slope `A = (1/(2π)) ∫ Im[-1/ε] dω = Ω_p² M_tot² / 2`
//! (`Ω_p` and `M_tot²` as in S2017 eq. (29)) is fixed by the optical ELF
//! alone (S2017 eqs. (28)-(29)), and its constant
//! `ln B = ln 8 - ⟨ln(ω_p + k_F²/3)⟩` (average weighted by the ELF) by the
//! SPA dispersion. That limit is our own expansion of the model
//! ([`SinglePolePenn::bethe_inverse_imfp_per_m`]); the tests check that the
//! model approaches it.

use super::bethe;
use super::quadrature::{integrate_segments, GaussLegendre};
use super::sum_rules::SumRuleReport;
use crate::constants::{BOHR_RADIUS, ELEMENTARY_CHARGE, HARTREE_ENERGY};
use crate::electron::data::{ElectronDataError, OpticalElf};
use std::f64::consts::PI;

type Result<T> = std::result::Result<T, ElectronDataError>;

/// Default relative tolerance of the IMFP, stopping and DIIMFP integrals.
pub const DEFAULT_RELATIVE_TOLERANCE: f64 = 1.0e-7;

/// The outer (loss) integral of the exchange-corrected total is this many
/// times looser than the tolerance of the inner integrals.
const EXCHANGE_OUTER_TOLERANCE_FACTOR: f64 = 100.0;

/// Gauss-Legendre order of the adaptive integrals.
const GL_ORDER: usize = 5;

/// The Hartree energy in eV.
pub(crate) fn hartree_ev() -> f64 {
    HARTREE_ENERGY / ELEMENTARY_CHARGE
}

/// `(3π/4)^(1/3)`: `k_F(ω_p) = KF · ω_p^(2/3)` (S2017, below eq. (8)).
fn kf_coefficient() -> f64 {
    (0.75 * PI).cbrt()
}

/// `k_F(ω_p)` in atomic units.
fn fermi_wavenumber(wp: f64) -> f64 {
    kf_coefficient() * wp.cbrt() * wp.cbrt()
}

/// `c = k_F(ω_p)² / 3`, the coefficient of `q²` in S2017 eq. (8).
fn dispersion_c(wp: f64) -> f64 {
    let kf = fermi_wavenumber(wp);
    kf * kf / 3.0
}

/// `ω_q(ω_p)` of S2017 eq. (8), given `c = k_F²/3`.
fn pole(wp: f64, c: f64, q: f64) -> f64 {
    let q2 = q * q;
    (wp * wp + c * q2 + 0.25 * q2 * q2).sqrt()
}

/// `q²` at which the pole of `ω_p` reaches the energy `w > ω_p`:
/// `q⁴/4 + c q² + ω_p² - w² = 0`, written without cancellation.
fn q2_at_pole_energy(wp: f64, c: f64, w: f64) -> f64 {
    let d = (w - wp) * (w + wp);
    2.0 * d / (c + (c * c + d).sqrt())
}

/// `ω₀(q, ω)`: the plasma frequency whose pole is at `ω` for momentum `q`,
/// the positive root `x = ω₀^(2/3)` of S2017 eq. (11). `None` if `q² >= 2ω`.
fn pole_plasma_frequency(q: f64, w: f64) -> Option<f64> {
    let q2 = q * q;
    let a = (PI * PI / 48.0).cbrt() * q2;
    let b = 0.25 * q2 * q2 - w * w;
    if b >= 0.0 {
        return None;
    }
    // f(x) = x³ + a x² + b is convex and increasing for x > 0, and
    // f(ω^(2/3)) = a ω^(4/3) + q⁴/4 >= 0: Newton from there converges
    // monotonically from the right.
    let mut x = w.cbrt() * w.cbrt();
    for _ in 0..200 {
        let f = x * x * x + a * x * x + b;
        let d = 3.0 * x * x + 2.0 * a * x;
        if d <= 0.0 {
            break;
        }
        let dx = f / d;
        x -= dx;
        if dx.abs() <= 1e-15 * x.abs() {
            break;
        }
    }
    (x > 0.0).then(|| x * x.sqrt())
}

/// The interval `[qa, qb]` of momenta for which the pole of `ω_p` is an
/// allowed loss at kinetic energy `t` above the Fermi level, with
/// `k = sqrt(2T')`. See the module docs.
fn allowed_q(wp: f64, c: f64, t: f64, k: f64) -> Option<(f64, f64)> {
    if wp.is_nan() || wp >= t {
        return None;
    }
    let phi = |q: f64| {
        let wq = pole(wp, c, q);
        (
            q * k - 0.5 * q * q - wq,
            k - q - (c * q + 0.5 * q * q * q) / wq,
        )
    };
    // Left end: Newton from q = 0, monotone from the left on a concave φ.
    let mut qa = 0.0;
    let mut converged = false;
    for _ in 0..200 {
        let (f, d) = phi(qa);
        if f >= 0.0 {
            converged = true;
            break;
        }
        if d <= 0.0 {
            return None; // past the maximum of φ with φ < 0: no allowed q
        }
        let dq = -f / d;
        qa += dq;
        if dq <= 1e-15 * qa {
            converged = true;
            break;
        }
    }
    if !converged {
        return None;
    }
    // Right end: Newton from q = k, monotone from the right. φ(k) < 0 and
    // φ'(k) < 0 because ω_q(k) > k²/2 for ω_p > 0, so k lies right of the
    // allowed interval.
    let mut qb = k;
    for _ in 0..200 {
        let (f, d) = phi(qb);
        if f >= 0.0 || d >= 0.0 {
            break;
        }
        let dq = -f / d;
        qb += dq;
        if dq.abs() <= 1e-15 * qb {
            break;
        }
    }
    let qt = q2_at_pole_energy(wp, c, t).sqrt();
    let qb = qb.min(qt);
    (qb > qa).then_some((qa, qb))
}

/// `∫_{qa}^{qb} dq / (q ω_q)` times `ω_p` (see the module docs).
fn imfp_kernel(wp: f64, c: f64, qa: f64, qb: f64) -> f64 {
    let (ya, yb) = (qa * qa, qb * qb);
    let g = |y: f64, q: f64| 2.0 * wp * wp + c * y + 2.0 * wp * pole(wp, c, q);
    0.5 * ((g(ya, qa) / ya) / (g(yb, qb) / yb)).ln()
}

/// The Born-Ochkur exchange factor `F = x² - x`, `x = recoil / (T - W)`
/// (dV2022, text below eq. (32); module docs), for a recoil energy
/// `q²/2` and `T - W` in the same unit. `T - W` must be positive.
fn ochkur_factor(recoil: f64, t_minus_emitted: f64) -> f64 {
    let x = recoil / t_minus_emitted;
    x * x - x
}

/// The Born-Ochkur exchange factor `F = -x + x²` of dV2022 (text below eq.
/// (32); module docs), with `x = (ħ²q²/2m) / (T - W)` and `W = E - B` the
/// kinetic energy of the emitted electron, for a recoil energy
/// `ħ²q²/2m = recoil_ev`, a kinetic energy `T = energy_ev` (use `T' = T + E_F`
/// with a Fermi energy), a loss `E = loss_ev` and a binding energy
/// `B = binding_ev`, all in eV. This is the factor the DIIMFP integrand of
/// [`SinglePolePenn::diimfp_with_binding_per_m_ev`] is multiplied by
/// (`1 + F`). `NaN` unless `T - W > 0`.
pub fn born_ochkur_factor(recoil_ev: f64, energy_ev: f64, loss_ev: f64, binding_ev: f64) -> f64 {
    let t_minus_emitted = energy_ev - (loss_ev - binding_ev);
    if t_minus_emitted > 0.0 {
        ochkur_factor(recoil_ev, t_minus_emitted)
    } else {
        f64::NAN
    }
}

/// Optional exchange correction (module docs, "Optional exchange
/// correction"), applied for kinetic energies strictly below
/// `applies_below_ev`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExchangeCorrection {
    applies_below_ev: f64,
}

impl ExchangeCorrection {
    /// A correction applied for `T < applies_below_ev` (finite, positive).
    pub fn new(applies_below_ev: f64) -> Result<Self> {
        check_energy("exchange cutoff energy", applies_below_ev)?;
        Ok(Self { applies_below_ev })
    }

    /// The energy below which the correction is applied, eV.
    pub fn applies_below_ev(&self) -> f64 {
        self.applies_below_ev
    }
}

/// One point of an IMFP and stopping-power table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InelasticPoint {
    /// Kinetic energy above the Fermi level, eV.
    pub energy_ev: f64,
    /// Inverse inelastic mean free path `λ⁻¹`, m⁻¹ (zero where no loss is
    /// allowed).
    pub inverse_imfp_per_m: f64,
    /// Stopping power `S = -dE/dx`, eV/m.
    pub stopping_ev_per_m: f64,
}

impl InelasticPoint {
    /// The inelastic mean free path `λ`, m (infinite where no loss is
    /// allowed).
    pub fn imfp_m(&self) -> f64 {
        1.0 / self.inverse_imfp_per_m
    }
}

/// The single-pole Penn model built from one optical ELF. See the module docs
/// for the equations, the units and the integration tolerance.
#[derive(Debug, Clone)]
pub struct SinglePolePenn {
    elf: OpticalElf,
    /// Tabulated energies, Hartree.
    w: Vec<f64>,
    /// Tabulated ELF values.
    e: Vec<f64>,
    /// Fermi energy, Hartree.
    fermi: f64,
    rel_tol: f64,
    gl: GaussLegendre,
    exchange: Option<ExchangeCorrection>,
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

fn check_binding(v: f64) -> Result<()> {
    if v.is_finite() && v >= 0.0 {
        Ok(())
    } else {
        Err(ElectronDataError::Invalid {
            what: "binding energy",
            reason: format!("must be finite and non-negative, got {v} eV"),
        })
    }
}

impl SinglePolePenn {
    /// The model for `elf`, with Fermi energy zero and the default tolerance.
    pub fn new(elf: OpticalElf) -> Self {
        let h = hartree_ev();
        let w = elf.energy_ev().iter().map(|x| x / h).collect();
        let e = elf.elf_values().to_vec();
        Self {
            elf,
            w,
            e,
            fermi: 0.0,
            rel_tol: DEFAULT_RELATIVE_TOLERANCE,
            gl: GaussLegendre::new(GL_ORDER),
            exchange: None,
        }
    }

    /// A model for another ELF with this model's Fermi energy, tolerance and
    /// exchange setting.
    pub(crate) fn same_settings_for(&self, elf: OpticalElf) -> Self {
        let mut m = Self::new(elf);
        m.fermi = self.fermi;
        m.rel_tol = self.rel_tol;
        m.exchange = self.exchange;
        m
    }

    /// Switch the exchange correction on (it is off by default).
    pub fn with_exchange(mut self, exchange: ExchangeCorrection) -> Self {
        self.exchange = Some(exchange);
        self
    }

    /// The exchange correction, if enabled.
    pub fn exchange(&self) -> Option<ExchangeCorrection> {
        self.exchange
    }

    /// Whether the exchange correction applies at kinetic energy `energy_ev`.
    pub fn exchange_applies_at(&self, energy_ev: f64) -> bool {
        self.exchange
            .is_some_and(|e| energy_ev < e.applies_below_ev)
    }

    /// Set the Fermi energy (eV, finite and non-negative). Energies passed to
    /// the model are then kinetic energies above the Fermi level, and the
    /// kinematics use `T' = T + E_F` (S2017 eq. (2)).
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

    /// Set the relative integration tolerance, in `[1e-12, 1e-2]`.
    pub fn with_relative_tolerance(mut self, rel_tol: f64) -> Result<Self> {
        if !(1e-12..=1e-2).contains(&rel_tol) {
            return Err(ElectronDataError::Invalid {
                what: "integration tolerance",
                reason: format!("must be in [1e-12, 1e-2], got {rel_tol}"),
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

    /// The optical ELF at `x` Hartree: linear in segment `i`
    /// (`w[i] <= x <= w[i+1]`), clamped to be non-negative.
    fn elf_in_segment(&self, i: usize, x: f64) -> f64 {
        let (x0, x1) = (self.w[i], self.w[i + 1]);
        let t = (x - x0) / (x1 - x0);
        (self.e[i] * (1.0 - t) + self.e[i + 1] * t).max(0.0)
    }

    /// The optical ELF at `x` Hartree; zero outside the table.
    fn elf_au(&self, x: f64) -> f64 {
        let n = self.w.len();
        if !(x >= self.w[0] && x <= self.w[n - 1]) {
            return 0.0;
        }
        let i = self.w.partition_point(|&g| g <= x).clamp(1, n - 1);
        self.elf_in_segment(i - 1, x)
    }

    /// `Im[-1/ε(q, ω)]` in the SPA (S2017 eq. (13)), at momentum transfer
    /// `q` (as a wavenumber, m⁻¹, so the momentum is `ħq`) and energy loss
    /// `ω` (eV). Zero where S2017 eq. (11) has no root (`ħ²q²/2m >= ω`) and
    /// where `ω₀` falls outside the tabulated range. At `q = 0` it is the
    /// optical ELF.
    pub fn loss_function(&self, q_per_m: f64, loss_ev: f64) -> f64 {
        if !(q_per_m.is_finite() && q_per_m >= 0.0 && loss_ev.is_finite() && loss_ev > 0.0) {
            return 0.0;
        }
        let q = q_per_m * BOHR_RADIUS;
        let w = loss_ev / hartree_ev();
        if q == 0.0 {
            return self.elf_au(w);
        }
        match pole_plasma_frequency(q, w) {
            Some(w0) => self.elf_au(w0) / (1.0 + PI * q * q / (6.0 * fermi_wavenumber(w0))),
            None => 0.0,
        }
    }

    /// The DIIMFP `p(T, ω)` (S2017 eq. (2), nonrelativistic) at kinetic
    /// energy `energy_ev` above the Fermi level and energy loss `loss_ev`,
    /// in m⁻¹ eV⁻¹. Zero for `ω <= 0` or `ω > T`, and, with the exchange
    /// correction applying at this energy, for `ω > T'/2` (module docs; no
    /// binding energy, `B = 0`).
    pub fn diimfp_per_m_ev(&self, energy_ev: f64, loss_ev: f64) -> Result<f64> {
        self.diimfp_with_binding_per_m_ev(energy_ev, loss_ev, 0.0)
    }

    /// The DIIMFP at kinetic energy `energy_ev` above the Fermi level and
    /// energy loss `loss_ev`, m⁻¹ eV⁻¹, for a loss that ionises a shell of
    /// binding energy `binding_ev` (finite, non-negative) and so emits an
    /// electron of kinetic energy `W = ω - B`. `B` enters only through the
    /// exchange correction, where it applies: the Born-Ochkur denominator is
    /// `T' - W = T' - ω + B` and the loss is limited to `ω <= (T' + B)/2`
    /// (dV2022 eqs. (32), (34), (35); module docs). Without exchange the
    /// result does not depend on `B`; with `B = 0` it is
    /// [`Self::diimfp_per_m_ev`].
    pub fn diimfp_with_binding_per_m_ev(
        &self,
        energy_ev: f64,
        loss_ev: f64,
        binding_ev: f64,
    ) -> Result<f64> {
        check_energy("electron energy", energy_ev)?;
        if !loss_ev.is_finite() {
            return Err(ElectronDataError::Invalid {
                what: "energy loss",
                reason: format!("must be finite, got {loss_ev} eV"),
            });
        }
        check_binding(binding_ev)?;
        let h = hartree_ev();
        let p_au = self.diimfp_au(energy_ev / h, loss_ev / h, binding_ev / h);
        Ok(p_au / (BOHR_RADIUS * h))
    }

    /// `p(T, ω)` in atomic units (per bohr per hartree) for a channel of
    /// binding energy `b`, with the exchange correction (and its
    /// `ω <= (T' + B)/2` limit) where it applies.
    fn diimfp_au(&self, t: f64, w: f64, b: f64) -> f64 {
        let exchange = self.exchange_applies_at(t * hartree_ev());
        if exchange && w > 0.5 * (t + self.fermi + b) {
            return 0.0;
        }
        self.diimfp_core_au(t, w, exchange, b)
    }

    /// `p(T, ω)` in atomic units, with the exchange factor `1 + F` in the
    /// integrand if `exchange`, `F` built with the emitted energy
    /// `W = ω - b` (`b` the channel's binding energy, Hartree); no limit on
    /// `ω` other than `ω <= T`.
    pub(crate) fn diimfp_core_au(&self, t: f64, w: f64, exchange: bool, b: f64) -> f64 {
        if !(w > 0.0 && w <= t) {
            return 0.0;
        }
        // T' - W, the Born-Ochkur denominator (module docs).
        let t_minus_emitted = t + self.fermi - w + b;
        if exchange && t_minus_emitted.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return 0.0;
        }
        let n = self.w.len();
        let (wmin, wmax) = (self.w[0], self.w[n - 1]);
        if w <= wmin {
            return 0.0; // ω₀ < ω for q > 0, and ω₀ must be in the table
        }
        let tp = t + self.fermi;
        let k = (2.0 * tp).sqrt();
        let s = (2.0 * (tp - w)).max(0.0).sqrt();
        let (q_minus, q_plus) = (2.0 * w / (k + s), k + s);
        // q where ω₀ = table knot w_j (decreasing in w_j). Knots below ω
        // only; ω₀ <= ω always.
        let q_of = |wj: f64| q2_at_pole_energy(wj, dispersion_c(wj), w).sqrt();
        let q_lo_table = if wmax < w { q_of(wmax) } else { 0.0 };
        let q_hi_table = q_of(wmin);
        let (lo, hi) = (q_minus.max(q_lo_table), q_plus.min(q_hi_table));
        if hi.is_nan() || hi <= lo {
            return 0.0;
        }
        // Breakpoints in ln q at the knots strictly inside (lo, hi).
        let mut breaks = vec![lo.ln()];
        let first_knot_below_w = self.w.partition_point(|&x| x < w);
        for j in (1..first_knot_below_w).rev() {
            let qj = q_of(self.w[j]);
            if qj > lo && qj < hi {
                breaks.push(qj.ln());
            }
        }
        breaks.push(hi.ln());
        let integral = integrate_segments(
            &self.gl,
            &mut |u: f64| {
                let q = u.exp();
                match pole_plasma_frequency(q, w) {
                    Some(w0) => {
                        let base =
                            self.elf_au(w0) / (1.0 + PI * q * q / (6.0 * fermi_wavenumber(w0)));
                        if exchange {
                            [base * (1.0 + ochkur_factor(0.5 * q * q, t_minus_emitted))]
                        } else {
                            [base]
                        }
                    }
                    None => [0.0],
                }
            },
            &breaks,
            self.rel_tol,
        )[0];
        integral / (PI * tp)
    }

    /// `(λ⁻¹, S)` at kinetic energy `t` (Hartree) above the Fermi level, in
    /// atomic units (per bohr, hartree per bohr).
    /// `b` is the channel's binding energy (Hartree), used only by the
    /// exchange correction.
    fn imfp_and_stopping_au(&self, t: f64, b: f64) -> (f64, f64) {
        let n = self.w.len();
        if t.is_nan() || t <= self.w[0] {
            return (0.0, 0.0);
        }
        if self.exchange_applies_at(t * hartree_ev()) {
            return self.imfp_and_stopping_exchange_au(t, b);
        }
        let tp = t + self.fermi;
        let k = (2.0 * tp).sqrt();
        let top = self.w[n - 1].min(t);
        let mut breaks: Vec<f64> = self.w.iter().copied().filter(|&x| x < top).collect();
        breaks.push(top);
        let r = integrate_segments(
            &self.gl,
            &mut |x: f64| {
                let c = dispersion_c(x);
                match allowed_q(x, c, t, k) {
                    Some((qa, qb)) => {
                        let weight = self.elf_au(x);
                        [
                            weight * imfp_kernel(x, c, qa, qb),
                            weight * x * (qb / qa).ln(),
                        ]
                    }
                    None => [0.0, 0.0],
                }
            },
            &breaks,
            self.rel_tol,
        );
        (r[0] / (PI * tp), r[1] / (PI * tp))
    }

    /// `(λ⁻¹, S)` with the exchange correction: the direct integral of the
    /// DIIMFP over `0 < ω <= min((T' + B)/2, T)` (module docs).
    fn imfp_and_stopping_exchange_au(&self, t: f64, b: f64) -> (f64, f64) {
        let top = (0.5 * (t + self.fermi + b)).min(t);
        if top <= self.w[0] {
            return (0.0, 0.0);
        }
        let mut breaks: Vec<f64> = self.w.iter().copied().filter(|&x| x < top).collect();
        breaks.push(top);
        let r = integrate_segments(
            &self.gl,
            &mut |x: f64| {
                let p = self.diimfp_core_au(t, x, true, b);
                [p, x * p]
            },
            &breaks,
            (self.rel_tol * EXCHANGE_OUTER_TOLERANCE_FACTOR).min(1e-2),
        );
        (r[0], r[1])
    }

    /// The inverse IMFP (m⁻¹) and the stopping power (eV/m) at kinetic
    /// energy `energy_ev` above the Fermi level. Both are zero where no loss
    /// is allowed (below the first energy at which a pole of the tabulated
    /// range fits inside the kinematic limits).
    pub fn imfp_and_stopping(&self, energy_ev: f64) -> Result<InelasticPoint> {
        self.imfp_and_stopping_with_binding(energy_ev, 0.0)
    }

    /// As [`Self::imfp_and_stopping`], for losses that ionise a shell of
    /// binding energy `binding_ev` (finite, non-negative): with the exchange
    /// correction applying, the Born-Ochkur denominator is `T' - ω + B` and
    /// the losses run up to `min((T' + B)/2, T)`
    /// ([`Self::diimfp_with_binding_per_m_ev`]). Without exchange it does not
    /// depend on `B`. The loss range is not cut at `B` from below: an ELF that
    /// should only act above the edge must start there (the model is zero
    /// below the first tabulated energy).
    pub fn imfp_and_stopping_with_binding(
        &self,
        energy_ev: f64,
        binding_ev: f64,
    ) -> Result<InelasticPoint> {
        check_energy("electron energy", energy_ev)?;
        check_binding(binding_ev)?;
        let h = hartree_ev();
        let (inv, s) = self.imfp_and_stopping_au(energy_ev / h, binding_ev / h);
        Ok(InelasticPoint {
            energy_ev,
            inverse_imfp_per_m: inv / BOHR_RADIUS,
            stopping_ev_per_m: s * h / BOHR_RADIUS,
        })
    }

    /// The inelastic mean free path `λ(E)`, m; infinite where no loss is
    /// allowed.
    pub fn imfp_m(&self, energy_ev: f64) -> Result<f64> {
        Ok(self.imfp_and_stopping(energy_ev)?.imfp_m())
    }

    /// The stopping power `S(E)`, eV/m.
    pub fn stopping_power_ev_per_m(&self, energy_ev: f64) -> Result<f64> {
        Ok(self.imfp_and_stopping(energy_ev)?.stopping_ev_per_m)
    }

    /// IMFP and stopping power on a grid of energies (eV above the Fermi
    /// level), in order. Sequential and deterministic.
    pub fn tabulate(&self, energies_ev: &[f64]) -> Result<Vec<InelasticPoint>> {
        energies_ev
            .iter()
            .map(|&e| self.imfp_and_stopping(e))
            .collect()
    }

    /// The model's own high-energy (Bethe) limit of the inverse IMFP, m⁻¹
    /// (see the module docs): `(1/(2π T')) ∫ dω_p ELF(ω_p) ln(8T'/(ω_p +
    /// k_F²/3))`, over the whole table.
    pub fn bethe_inverse_imfp_per_m(&self, energy_ev: f64) -> Result<f64> {
        check_energy("electron energy", energy_ev)?;
        let h = hartree_ev();
        let tp = energy_ev / h + self.fermi;
        let r = integrate_segments(
            &self.gl,
            &mut |x: f64| [self.elf_au(x) * (8.0 * tp / (x + dispersion_c(x))).ln()],
            &self.w,
            self.rel_tol,
        )[0];
        Ok(r / (2.0 * PI * tp) / BOHR_RADIUS)
    }

    /// The nonrelativistic Bethe stopping power with `W_max = T`, the
    /// electron density from the f-sum of the same ELF and the mean
    /// excitation energy `I` of the same ELF ([`SumRuleReport`],
    /// [`bethe::stopping_power_ev_per_m`]), eV/m. With a Fermi energy, the
    /// kinematic energy is `T' = T + E_F` and `W_max = T`.
    pub fn bethe_stopping_ev_per_m(&self, energy_ev: f64) -> Result<f64> {
        check_energy("electron energy", energy_ev)?;
        let r = SumRuleReport::new(&self.elf);
        Ok(bethe::stopping_power_ev_per_m(
            r.electron_density_per_m3,
            r.mean_excitation_energy_ev,
            energy_ev + self.fermi_energy_ev(),
            energy_ev,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closed_form_imfp_kernel_matches_quadrature() {
        let gl = GaussLegendre::new(GL_ORDER);
        for (wp, qa, qb) in [(0.5, 1e-3, 3.0), (2.0, 0.2, 40.0), (0.05, 1e-5, 0.4)] {
            let c = dispersion_c(wp);
            let want = wp
                * integrate_segments(
                    &gl,
                    &mut |u: f64| [1.0 / pole(wp, c, u.exp())],
                    &[f64::ln(qa), f64::ln(qb)],
                    1e-13,
                )[0];
            let got = imfp_kernel(wp, c, qa, qb);
            assert!((got / want - 1.0).abs() < 1e-11, "{wp}: {got} vs {want}");
        }
    }

    #[test]
    fn cubic_root_satisfies_the_dispersion() {
        for (q, w) in [(0.1, 1.0), (1.0, 0.9), (1e-4, 3.0), (2.0, 2.0001)] {
            let w0 = pole_plasma_frequency(q, w).unwrap();
            let back = pole(w0, dispersion_c(w0), q);
            assert!((back / w - 1.0).abs() < 1e-13, "{q} {w}: {back}");
        }
        assert!(pole_plasma_frequency(2.0, 1.9).is_none());
    }

    #[test]
    fn allowed_interval_ends_satisfy_the_kinematics() {
        let (t, wp): (f64, f64) = (50.0, 0.7);
        let k = (2.0 * t).sqrt();
        let c = dispersion_c(wp);
        let (qa, qb) = allowed_q(wp, c, t, k).unwrap();
        let phi = |q: f64| q * k - 0.5 * q * q - pole(wp, c, q);
        assert!(phi(qa).abs() < 1e-12, "{}", phi(qa));
        // qb is the kinematic root or q_T, whichever is smaller.
        let qt = q2_at_pole_energy(wp, c, t).sqrt();
        assert!(phi(qb).abs() < 1e-10 || (qb - qt).abs() < 1e-12 * qt);
        assert!(allowed_q(wp, c, 0.5, 1.0).is_none());
    }

    /// Eq. (32) of de Vera et al.: with `T = 400`, `B = 150` and `ω = 200` eV
    /// (`E_F = 0`), `W = ω - B = 50` eV and `T - W = 350` eV, not the `200` eV
    /// of `T' - ω`.
    #[test]
    fn exchange_factor_uses_primary_energy_less_emitted_energy() {
        let h = hartree_ev();
        let (t, b, w) = (400.0 / h, 150.0 / h, 200.0 / h);
        let t_minus_w = t - (w - b);
        assert!((t_minus_w * h - 350.0).abs() < 1e-9);
        for q in [0.5, 2.0, 5.0] {
            let x = 0.5 * q * q / (350.0 / h);
            let want = -x + x * x;
            let got = ochkur_factor(0.5 * q * q, t_minus_w);
            assert!((got - want).abs() <= 1e-14 * want.abs(), "{q}: {got}");
            // The wrong denominator would give a different factor.
            let wrong = ochkur_factor(0.5 * q * q, t - w);
            assert!((wrong - want).abs() > 1e-3 * want.abs(), "{q}");
        }
    }

    #[test]
    fn diimfp_exchange_depends_on_binding_energy() {
        use crate::electron::inelastic::{DrudeLorentz, DrudeLorentzOscillator};
        let elf = DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])
            .unwrap()
            .to_optical_elf("synthetic Drude plasmon", 0.05, 2e4, 240)
            .unwrap();
        let m = SinglePolePenn::new(elf);
        let h = hartree_ev();
        let (t, w, b) = (400.0 / h, 200.0 / h, 150.0 / h);
        let direct = m.diimfp_core_au(t, w, false, 0.0);
        let b0 = m.diimfp_core_au(t, w, true, 0.0);
        let bj = m.diimfp_core_au(t, w, true, b);
        assert!(direct > 0.0 && b0 > 0.0 && bj > 0.0);
        assert!((bj / b0 - 1.0).abs() > 1e-3, "{bj} vs {b0}");
        // Without exchange the binding energy does not enter.
        assert_eq!(m.diimfp_core_au(t, w, false, b), direct);
    }
}
