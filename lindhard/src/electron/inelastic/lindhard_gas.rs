//! The Lindhard dielectric function of a free-electron gas at zero
//! temperature (random-phase approximation), in Hartree atomic units
//! (`m_e = e = ħ = 1`), and the loss-function integrals built from it.
//!
//! This is the building block of the full Penn algorithm
//! ([`super::full_penn`]): Penn expands the optical energy-loss function over
//! the loss functions of free-electron gases of all plasma frequencies
//! (Shinotsuka et al., Surf. Interface Anal. 49, 238 (2017),
//! doi:10.1002/sia.6123, eq. (4), with `ε_L` "the dielectric function from
//! the Lindhard model of a free electron gas", citing J. Lindhard, Mat. Fys.
//! Medd. Dan. Vid. Selsk. 28, no. 8 (1954)).
//!
//! # Sources and what was and was not read
//!
//! Lindhard (1954) was **not opened** (no open copy was found). The function
//! is implemented in its standard closed form, written here in the variables
//! `z = q/(2 k_F)` and `u = ω/(q k_F)`:
//!
//! ```text
//! Re ε(q, ω) = 1 + (k_TF²/q²) [ 1/2 + (A(z - u) + A(z + u)) / (8 z) ],
//! A(s) = (1 - s²) ln|(s + 1)/(s - 1)|,      k_TF² = 4 k_F/π,
//! Im ε(q, ω) = A_c(q, ω) / q³,
//! A_c = (k_F² - a²)₊ - (k_F² - (a + q)²)₊,    a = ω/q - q/2,    ω > 0.
//! ```
//!
//! These were compared with an opened source: A. V. Latyshev and A. A.
//! Yushkanov, arXiv:1212.6260v1 (open-access preprint, PDF pages read
//! 2026-10-07), section 5 "Degenerate plasma": the integral `B(k, z)` on
//! p. 14 (`∝ ∫_{-1}^{1} (1 - P_x²) dP_x / ((P_x - z/k)² - (k/2)²)`, with
//! `k = q/k_F`, `z = (ω + iν)/(k_F v_F)`), its closed form `b(k, z)` on
//! pp. 14-15, and the dielectric function eq. (5.5″) on p. 15,
//! `ε = 1 - (3 x_p²/4k²)(x + iy) b(k,z) b(k,0) / (x b(k,0) + i y b(k,z))`
//! with `x_p = ω_p/(k_F v_F)`. Its abstract states that for collision
//! frequency `ν -> 0` (`y -> 0`) the result is Lindhard's formula, which is
//! then `ε = 1 - (3 x_p²/4k²) b(k, x)`. Compared numerically (atomic units,
//! `v_F = k_F`, `ω_p² = 4 k_F³/(3π)`) at six `(k_F, q, ω)` points:
//!
//! * the real part above equals the closed form `b(k, z)` of p. 15 to
//!   1e-10 (the `-2` and the two `ln` terms; their `ln` arguments are the
//!   ones of `A(z ∓ u)` combined);
//! * the imaginary part equals the residue (`+i0`) contribution of the
//!   p. 14 integral at its two poles `P_x = z/k ± k/2` when they lie in
//!   `[-1, 1]`, i.e. our `A_c/q³`, to 1e-15 (the printed closed form `b`
//!   gives only the real part under principal-branch logarithms, so the
//!   imaginary part is compared through the integral, not the closed form).
//!
//! This is a comparison with a secondary source that derives the formula, not
//! with Lindhard's paper; the original equation numbers are not known here.
//! The further checks are:
//!
//! * the imaginary part is derived here directly from the Fermi golden rule
//!   (`Im ε = (4π/q²) π · 2 ∫ d³k/(2π)³ f_k (1 - f_{k+q}) δ(ω - ΔE)`, an
//!   elementary Fermi-sphere area), and the unit tests check the real part
//!   against it by a **Kramers-Kronig integral**;
//! * the **static limit** `ε(q, 0) -> 1 + k_TF²/q²` for `q -> 0`
//!   (Thomas-Fermi screening) and the exact static form
//!   `1 + (k_TF²/q²)(1/2 + (1 - z²)/(4z) ln|(1+z)/(1-z)|)`;
//! * the **long-wavelength limit** `ε(q -> 0, ω) -> 1 - ω_p²/ω²`;
//! * the **f-sum rule** `∫ ω Im[-1/ε] dω = (π/2) ω_p²` at every `q`
//!   (plasmon delta function included), the property Penn's algorithm relies
//!   on (it is what makes the expanded ELF satisfy the optical f-sum rule).
//!
//! The plasma frequency and Fermi wavenumber are related by
//! `ω_p² = 4π n`, `n = k_F³/(3π²)`, i.e. `k_F = (3π/4)^(1/3) ω_p^(2/3)`
//! (S2017 below eq. (8)).
//!
//! # The loss function
//!
//! For `ω > 0`, `Im[-1/ε] = Im ε / |ε|²` inside the electron-hole continuum
//! (`|a| < k_F`) and, above the continuum, a delta function at the plasmon
//! energy `ω_pl(q)` (the root of `Re ε = 0` above the continuum edge
//! `ω_+ = q k_F + q²/2`) of weight `π / (∂Re ε/∂ω)`. That weight is the one
//! of S2017 eq. (7) at `q -> 0`: with `Re ε = 1 - ω_p²/ω²`, `∂Re ε/∂ω =
//! 2/ω_p` at `ω = ω_p` and the weight is `π ω_p/2`.

use super::quadrature::{integrate_segments, GaussLegendre};
use std::f64::consts::PI;

/// `(3π/4)^(1/3)`: `k_F = KF_COEFFICIENT · ω_p^(2/3)`.
pub(crate) fn kf_coefficient() -> f64 {
    (0.75 * PI).cbrt()
}

/// `k_F(ω_p)` in atomic units.
pub fn fermi_wavenumber_au(plasma_frequency: f64) -> f64 {
    kf_coefficient() * plasma_frequency.cbrt() * plasma_frequency.cbrt()
}

/// `(1 - s²) ln|(s + 1)/(s - 1)|`, zero at `s = ±1` (the limit).
fn a_fn(s: f64) -> f64 {
    let d = 1.0 - s * s;
    if d == 0.0 {
        0.0
    } else {
        d * ((s + 1.0) / (s - 1.0)).abs().ln()
    }
}

/// The bracket `1/2 + (A(z - u) + A(z + u))/(8z)` of `Re ε`.
///
/// For `u - z > 4` (far above the continuum) the direct form loses all its
/// digits to cancellation (the bracket is `O(1/u²)` but its terms are
/// `O(u/z)`), so the expansion of `A(s) = -2s + Σ_k 4 s^-(2k+1)/((2k+1)(2k+3))`
/// (from `ln|(s+1)/(s-1)| = 2 Σ_k s^-(2k+1)/(2k+1)`, `|s| > 1`) is used; with
/// `A` odd, the bracket is `-(1/4) Σ_k c_k S_m/(u² - z²)^m`, `m = 2k+1`,
/// `c_k = 4/((2k+1)(2k+3))` and `S_m = Σ_j (u+z)^(m-1-j) (u-z)^j`, every term
/// positive. The unit tests compare both forms where they overlap.
fn shape(z: f64, u: f64) -> f64 {
    if u - z > 4.0 {
        shape_series(z, u)
    } else {
        shape_direct(z, u)
    }
}

fn shape_direct(z: f64, u: f64) -> f64 {
    0.5 + (a_fn(z - u) + a_fn(z + u)) / (8.0 * z)
}

fn shape_series(z: f64, u: f64) -> f64 {
    let (hi, lo) = (u + z, u - z);
    let d = hi * lo;
    let mut sum = 0.0;
    let mut dpow = d; // (u² - z²)^m
    for k in 0..200 {
        let m = 2 * k + 1;
        let mut sm = 0.0;
        for j in 0..m {
            sm += hi.powi(m - 1 - j) * lo.powi(j);
        }
        let ck = 4.0 / ((2 * k + 1) as f64 * (2 * k + 3) as f64);
        let term = ck * sm / dpow;
        sum += term;
        if term <= 1e-17 * sum {
            break;
        }
        dpow *= d * d;
    }
    -0.25 * sum
}

/// `∂(shape)/∂u = (A'(z + u) - A'(z - u))/(8z)`; for `u - z > 4` the
/// series `(1/4) Σ_k c_k m S_{m+1}/(u² - z²)^(m+1)` (all terms positive, from
/// `A'(s) = -2 - Σ_k c_k m s^-(m+1)`), otherwise the direct form.
fn shape_du(z: f64, u: f64) -> f64 {
    if u - z <= 4.0 {
        return (a_prime(z + u) - a_prime(z - u)) / (8.0 * z);
    }
    let (hi, lo) = (u + z, u - z);
    let d = hi * lo;
    let mut sum = 0.0;
    let mut dpow = d * d; // (u² - z²)^(m+1)
    for k in 0..200 {
        let m = 2 * k + 1;
        let n = m + 1;
        let mut sn = 0.0;
        for j in 0..n {
            sn += hi.powi(n - 1 - j) * lo.powi(j);
        }
        let ck = 4.0 / ((2 * k + 1) as f64 * (2 * k + 3) as f64);
        let term = ck * m as f64 * sn / dpow;
        sum += term;
        if term <= 1e-17 * sum {
            break;
        }
        dpow *= d * d;
    }
    0.25 * sum
}

/// `dA/ds = 2 - 2 s ln|(s+1)/(s-1)|`.
fn a_prime(s: f64) -> f64 {
    2.0 - 2.0 * s * ((s + 1.0) / (s - 1.0)).abs().ln()
}

/// The undamped plasmon of a Lindhard gas at one momentum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LindhardPlasmon {
    /// Plasmon energy `ω_pl(q)`, Hartree.
    pub energy: f64,
    /// Weight of the delta function in `Im[-1/ε]`: `π/(∂Re ε/∂ω)`, Hartree.
    pub weight: f64,
}

/// A free-electron gas of Fermi wavenumber `k_F` (bohr⁻¹), zero temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LindhardGas {
    kf: f64,
}

impl LindhardGas {
    /// The gas with Fermi wavenumber `kf > 0` (atomic units).
    pub fn from_fermi_wavenumber(kf: f64) -> Self {
        debug_assert!(kf > 0.0 && kf.is_finite());
        Self { kf }
    }

    /// The gas with plasma frequency `wp > 0` (Hartree; `ω_p² = 4π n`).
    pub fn from_plasma_frequency(wp: f64) -> Self {
        Self::from_fermi_wavenumber(fermi_wavenumber_au(wp))
    }

    /// `k_F`, bohr⁻¹.
    pub fn fermi_wavenumber(&self) -> f64 {
        self.kf
    }

    /// `ω_p = sqrt(4 k_F³/(3π))`, Hartree.
    pub fn plasma_frequency(&self) -> f64 {
        (4.0 * self.kf * self.kf * self.kf / (3.0 * PI)).sqrt()
    }

    /// `k_TF² = 4 k_F/π`, bohr⁻².
    pub fn thomas_fermi_wavenumber_sq(&self) -> f64 {
        4.0 * self.kf / PI
    }

    /// The upper edge of the electron-hole continuum, `q k_F + q²/2`.
    pub fn continuum_upper_edge(&self, q: f64) -> f64 {
        q * self.kf + 0.5 * q * q
    }

    /// `Re ε(q, ω)` for `q > 0`, `ω >= 0`.
    pub fn re_epsilon(&self, q: f64, w: f64) -> f64 {
        let z = q / (2.0 * self.kf);
        let u = w / (q * self.kf);
        1.0 + self.thomas_fermi_wavenumber_sq() / (q * q) * shape(z, u)
    }

    /// `∂Re ε/∂ω` (valid off the singular points `|z ± u| = 1`).
    fn d_re_epsilon_dw(&self, q: f64, w: f64) -> f64 {
        let z = q / (2.0 * self.kf);
        let u = w / (q * self.kf);
        let df = shape_du(z, u);
        self.thomas_fermi_wavenumber_sq() / (q * q) * df / (q * self.kf)
    }

    /// `Im ε(q, ω)` for `q > 0`, `ω > 0`: zero outside the continuum.
    pub fn im_epsilon(&self, q: f64, w: f64) -> f64 {
        let kf2 = self.kf * self.kf;
        let a = w / q - 0.5 * q;
        let b = w / q + 0.5 * q;
        let area = (kf2 - a * a).max(0.0) - (kf2 - b * b).max(0.0);
        area / (q * q * q)
    }

    /// `ε(q, ω)` as `(Re, Im)`.
    pub fn epsilon(&self, q: f64, w: f64) -> (f64, f64) {
        (self.re_epsilon(q, w), self.im_epsilon(q, w))
    }

    /// The continuum part of the loss function `Im[-1/ε(q, ω)]` (zero
    /// outside the electron-hole continuum; the plasmon delta function is
    /// separate, see [`Self::plasmon`]).
    pub fn loss_continuum(&self, q: f64, w: f64) -> f64 {
        let im = self.im_epsilon(q, w);
        if im <= 0.0 {
            return 0.0;
        }
        let re = self.re_epsilon(q, w);
        im / (re * re + im * im)
    }

    /// Whether [`Self::plasmon`] finds a plasmon at `q`: the sign test at the
    /// head of that method, without the root search (a single evaluation of
    /// `Re ε`).
    pub(crate) fn plasmon_exists(&self, q: f64) -> bool {
        if !(q > 0.0 && q.is_finite()) {
            return false;
        }
        let lo0 = self.continuum_upper_edge(q) * (1.0 + 1e-12);
        self.re_epsilon(q, lo0) < 0.0
    }

    /// Whether the plasmon at `q` exists and has energy at least `w`, without
    /// solving for the plasmon energy. Above the continuum edge `Re ε` rises
    /// monotonically from a negative value to 1, so the plasmon energy is
    /// `>= w` exactly when `Re ε(q, w) <= 0` (and always, for `w` at or below
    /// the edge).
    pub(crate) fn plasmon_reaches(&self, q: f64, w: f64) -> bool {
        if !self.plasmon_exists(q) {
            return false;
        }
        w <= self.continuum_upper_edge(q) * (1.0 + 1e-12) || self.re_epsilon(q, w) <= 0.0
    }

    /// Whether the plasmon at `q` exists and has energy below `w`; the
    /// complement of [`Self::plasmon_reaches`] among the `q` with a plasmon.
    pub(crate) fn plasmon_below(&self, q: f64, w: f64) -> bool {
        self.plasmon_exists(q) && !self.plasmon_reaches(q, w)
    }

    /// The plasmon at momentum `q`, if it exists (above the continuum).
    pub fn plasmon(&self, q: f64) -> Option<LindhardPlasmon> {
        if !(q > 0.0 && q.is_finite()) {
            return None;
        }
        let edge = self.continuum_upper_edge(q);
        let lo0 = edge * (1.0 + 1e-12);
        if self.re_epsilon(q, lo0) >= 0.0 {
            return None;
        }
        // Re ε increases from negative at the edge to 1 at infinity.
        let mut hi = 2.0 * (edge + self.plasma_frequency());
        for _ in 0..80 {
            if self.re_epsilon(q, hi) > 0.0 {
                break;
            }
            hi *= 2.0;
        }
        let (mut lo, mut x) = (lo0, 0.5 * (lo0 + hi));
        for _ in 0..200 {
            let f = self.re_epsilon(q, x);
            if f < 0.0 {
                lo = x;
            } else {
                hi = x;
            }
            let d = self.d_re_epsilon_dw(q, x);
            let newton = x - f / d;
            let next = if newton > lo && newton < hi {
                newton
            } else {
                0.5 * (lo + hi)
            };
            let done = (next - x).abs() <= 1e-15 * next.abs();
            x = next;
            if done || hi - lo <= 1e-15 * hi {
                break;
            }
        }
        let d = self.d_re_epsilon_dw(q, x);
        (d > 0.0).then(|| LindhardPlasmon {
            energy: x,
            weight: PI / d,
        })
    }

    /// The plasmon cutoff momentum `q_c`: the plasmon exists for `q < q_c`
    /// (it enters the continuum there and is Landau damped).
    pub fn plasmon_cutoff(&self) -> f64 {
        let exists =
            |q: f64| self.re_epsilon(q, self.continuum_upper_edge(q) * (1.0 + 1e-12)) < 0.0;
        let mut hi = 1e-3 * self.kf;
        while exists(hi) && hi < 1e6 * self.kf {
            hi *= 2.0;
        }
        let mut lo = 0.5 * hi;
        if !exists(lo) {
            return lo; // pathological: no plasmon at all
        }
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if exists(mid) {
                lo = mid;
            } else {
                hi = mid;
            }
            if hi - lo <= 1e-15 * hi {
                break;
            }
        }
        0.5 * (lo + hi)
    }

    /// Break points in `ω` (ascending, inside `(0, wcut)`) where the
    /// continuum loss function has a kink.
    fn continuum_breaks(&self, q: f64, wcut: f64) -> Vec<f64> {
        let h = 0.5 * q * q;
        let qk = q * self.kf;
        let mut v = vec![0.0];
        for p in [h - qk, qk - h, h + qk] {
            if p > 0.0 && p < wcut {
                v.push(p);
            }
        }
        v.push(wcut);
        v.sort_by(f64::total_cmp);
        v.dedup();
        v
    }

    /// `[∫_0^{wcut} Im[-1/ε] dω, ∫_0^{wcut} ω Im[-1/ε] dω]` at momentum `q`,
    /// the plasmon delta function included, with the adaptive integrator's
    /// relative tolerance `rel_tol`.
    pub(crate) fn loss_moments(
        &self,
        gl: &GaussLegendre,
        q: f64,
        wcut: f64,
        rel_tol: f64,
    ) -> [f64; 2] {
        if wcut.is_nan() || wcut <= 0.0 {
            return [0.0; 2];
        }
        let breaks = self.continuum_breaks(q, wcut);
        let mut r = integrate_segments(
            gl,
            &mut |w: f64| {
                let l = self.loss_continuum(q, w);
                [l, w * l]
            },
            &breaks,
            rel_tol,
        );
        // the root search for the plasmon only if it is inside the range
        if self.plasmon_below(q, wcut) {
            if let Some(p) = self.plasmon(q) {
                if p.energy < wcut {
                    r[0] += p.weight;
                    r[1] += p.weight * p.energy;
                }
            }
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gas() -> LindhardGas {
        LindhardGas::from_plasma_frequency(0.7)
    }

    #[test]
    fn kf_and_plasma_frequency_are_consistent() {
        let wp = 0.7;
        let g = LindhardGas::from_plasma_frequency(wp);
        assert!((g.plasma_frequency() / wp - 1.0).abs() < 1e-14);
        // ω_p² = 4π n with n = k_F³/(3π²)
        let n = g.kf.powi(3) / (3.0 * PI * PI);
        assert!((wp * wp / (4.0 * PI * n) - 1.0).abs() < 1e-14);
        // k_TF² = 3 ω_p²/k_F² (= 4 k_F/π)
        let tf = 3.0 * wp * wp / (g.kf * g.kf);
        assert!((g.thomas_fermi_wavenumber_sq() / tf - 1.0).abs() < 1e-14);
    }

    /// Static limit: Thomas-Fermi screening, and the exact static form.
    #[test]
    fn static_limit_is_thomas_fermi_screening() {
        let g = gas();
        let q = 1e-3;
        let tf = 1.0 + g.thomas_fermi_wavenumber_sq() / (q * q);
        let re = g.re_epsilon(q, 0.0);
        assert!((re / tf - 1.0).abs() < 1e-5, "{re} vs {tf}");
        // exact static form, several q (including q > 2 k_F)
        for q in [0.3, 1.0, 1.7, 3.0, 5.0] {
            let z = q / (2.0 * g.kf);
            let want = 1.0
                + g.thomas_fermi_wavenumber_sq() / (q * q)
                    * (0.5 + (1.0 - z * z) / (4.0 * z) * ((1.0 + z) / (1.0 - z)).abs().ln());
            let got = g.re_epsilon(q, 0.0);
            assert!((got / want - 1.0).abs() < 1e-12, "q={q}: {got} vs {want}");
        }
        assert!(g.im_epsilon(0.5, 1e-12) < 1e-9);
    }

    /// Long-wavelength limit: ε -> 1 - ω_p²/ω² (Drude).
    #[test]
    fn long_wavelength_limit_is_drude() {
        let g = gas();
        let wp = g.plasma_frequency();
        for w in [1.5 * wp, 3.0 * wp, 0.8 * wp] {
            let q = 1e-3;
            let want = 1.0 - wp * wp / (w * w);
            let got = g.re_epsilon(q, w);
            assert!((got - want).abs() < 1e-5, "w={w}: {got} vs {want}");
            assert_eq!(g.im_epsilon(q, w), 0.0);
        }
        // the plasmon sits at ω_p with weight π ω_p/2 as q -> 0
        let p = g.plasmon(1e-3).unwrap();
        assert!((p.energy / wp - 1.0).abs() < 1e-4, "{}", p.energy);
        assert!(
            (p.weight / (0.5 * PI * wp) - 1.0).abs() < 1e-3,
            "{}",
            p.weight
        );
    }

    /// The first correction to the long-wavelength limit (Bohm-Pines):
    /// ω_pl² = ω_p² + (3/5) k_F² q² + O(q⁴).
    #[test]
    fn plasmon_dispersion_starts_as_bohm_pines() {
        let g = gas();
        let (wp, kf) = (g.plasma_frequency(), g.kf);
        let q = 0.02;
        let p = g.plasmon(q).unwrap();
        let want = (wp * wp + 0.6 * kf * kf * q * q).sqrt();
        assert!(
            (p.energy / want - 1.0).abs() < 2e-4,
            "{} vs {want}",
            p.energy
        );
    }

    #[test]
    fn large_u_series_matches_the_direct_form_where_both_are_accurate() {
        for (z, u) in [(0.3, 4.5), (0.5, 6.0), (1.0, 5.5), (0.05, 5.0)] {
            let (a, b) = (shape_series(z, u), shape_direct(z, u));
            assert!(
                (a - b).abs() < 1e-10 * a.abs().max(1e-3),
                "{z} {u}: {a} vs {b}"
            );
        }
    }

    #[test]
    fn shape_derivative_matches_a_finite_difference() {
        for (z, u) in [(0.3, 4.5), (0.5, 6.0), (0.05, 2.0), (0.8, 1.2), (0.2, 9.0)] {
            let h = 1e-5 * u;
            let fd = (shape(z, u + h) - shape(z, u - h)) / (2.0 * h);
            let an = shape_du(z, u);
            assert!((fd / an - 1.0).abs() < 1e-6, "{z} {u}: {fd} vs {an}");
        }
    }

    #[test]
    fn continuum_has_the_known_low_frequency_form() {
        // ω -> 0, q < 2 k_F: Im ε = 2ω/q³ (the Fermi-sphere area is q(2a+q) = 2ω).
        let g = gas();
        let (q, w) = (0.4, 1e-4);
        assert!((g.im_epsilon(q, w) / (2.0 * w / q.powi(3)) - 1.0).abs() < 1e-12);
        // no absorption above the upper edge
        assert_eq!(g.im_epsilon(q, g.continuum_upper_edge(q) * 1.0001), 0.0);
    }

    /// Kramers-Kronig: Re ε(q, ω) - 1 = (2/π) P∫ ω' Im ε(q, ω')/(ω'² - ω²) dω'.
    /// This checks the closed-form real part against the Fermi-golden-rule
    /// imaginary part, independently.
    #[test]
    fn real_part_obeys_kramers_kronig_against_the_imaginary_part() {
        let g = gas();
        let gl = GaussLegendre::new(8);
        for (q, w) in [(0.3, 0.2), (0.8, 1.9), (1.5, 0.4), (2.5, 3.0), (0.5, 0.05)] {
            // Subtract the singular point: P∫ f(ω')/(ω'²-ω²) with
            // f = ω' Im ε(ω') -> ∫ [f(ω') - f(ω)]/(ω'²-ω²) dω' + f(ω) P∫ 1/(ω'²-ω²),
            // and P∫_0^Λ dω'/(ω'²-ω²) = (1/(2ω)) ln|(Λ-ω)/(Λ+ω)|.
            let lam = g.continuum_upper_edge(q) * 1.0 + 1e-300;
            let f = |x: f64| x * g.im_epsilon(q, x);
            let fw = f(w);
            let mut breaks = g.continuum_breaks(q, lam);
            if w < lam && !breaks.contains(&w) {
                breaks.push(w);
                breaks.sort_by(f64::total_cmp);
            }
            let sub = integrate_segments(
                &gl,
                &mut |x: f64| {
                    let d = x * x - w * w;
                    if d.abs() < 1e-14 {
                        [0.0]
                    } else {
                        [(f(x) - fw) / d]
                    }
                },
                &breaks,
                1e-10,
            )[0];
            let pv = if w < lam {
                fw * (1.0 / (2.0 * w)) * ((lam - w) / (lam + w)).ln()
            } else {
                0.0
            };
            let kk = 1.0 + (2.0 / PI) * (sub + pv);
            let re = g.re_epsilon(q, w);
            assert!(
                (kk - re).abs() < 2e-5 * (1.0 + re.abs()),
                "q={q} w={w}: KK {kk} vs {re}"
            );
        }
    }

    /// f-sum rule of one Lindhard gas, plasmon included, at several q.
    #[test]
    fn single_gas_f_sum_rule_holds_at_every_q() {
        let g = gas();
        let gl = GaussLegendre::new(8);
        let want = 0.5 * PI * g.plasma_frequency().powi(2);
        // q well below, near and above the plasmon cutoff, and above 2 k_F
        let qc = g.plasmon_cutoff();
        assert!(qc > 0.0 && qc < 2.0 * g.kf, "qc = {qc}");
        for q in [
            0.01,
            0.3 * qc,
            0.9 * qc,
            1.2 * qc,
            2.0 * qc,
            2.5 * g.kf,
            5.0,
        ] {
            let wcut = 20.0 * (g.continuum_upper_edge(q) + g.plasma_frequency());
            let [_, s1] = g.loss_moments(&gl, q, wcut, 1e-10);
            assert!((s1 / want - 1.0).abs() < 1e-6, "q = {q}: {s1} vs {want}");
        }
    }

    /// The sign tests used by the full Penn searches agree with the root
    /// search of `plasmon` (away from the root itself).
    #[test]
    fn plasmon_predicates_agree_with_the_root_search() {
        let g = gas();
        let qc = g.plasmon_cutoff();
        for q in [1e-3, 0.2 * qc, 0.7 * qc, 0.99 * qc, 1.01 * qc, 3.0 * qc] {
            let p = g.plasmon(q);
            assert_eq!(g.plasmon_exists(q), p.is_some(), "q = {q}");
            for f in [0.5, 0.9, 0.999, 1.001, 1.1, 2.0] {
                let w = f * p.map_or(g.continuum_upper_edge(q), |p| p.energy);
                let reaches = p.is_some_and(|p| p.energy >= w);
                assert_eq!(g.plasmon_reaches(q, w), reaches, "q = {q}, w = {w}");
                assert_eq!(
                    g.plasmon_below(q, w),
                    p.is_some_and(|p| p.energy < w),
                    "q = {q}, w = {w}"
                );
            }
        }
        assert!(!g.plasmon_exists(0.0) && !g.plasmon_exists(f64::NAN));
    }

    #[test]
    fn plasmon_cutoff_is_where_the_plasmon_meets_the_continuum() {
        let g = gas();
        let qc = g.plasmon_cutoff();
        assert!(g.plasmon(0.999 * qc).is_some());
        assert!(g.plasmon(1.001 * qc).is_none());
        let p = g.plasmon(0.999999 * qc).unwrap();
        assert!((p.energy / g.continuum_upper_edge(qc) - 1.0).abs() < 1e-3);
    }
}
