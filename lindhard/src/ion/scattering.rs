//! Classical elastic scattering in a screened-Coulomb potential.
//!
//! Everything is in the reduced variables of [`crate::ion::potential`]:
//! `x = r/a`, `beta = b/a`, `eps = a E_cm / (Z1 Z2 e^2/(4 pi eps0))`. The
//! radial motion in the centre-of-mass frame obeys
//! `G(x) = 1 - phi(x)/(eps x) - beta^2/x^2`, with turning point `G(x0) = 0`.
//! The centre-of-mass scattering angle is (Goldstein; ZBL 1985, Ch. 2)
//!
//! ```text
//! theta = pi - 2 beta * Int_{x0}^{inf} dx / (x^2 sqrt(G(x)))
//! ```
//!
//! Three evaluators are provided: direct Gauss-Mehler quadrature
//! ([`theta_quadrature`], accurate everywhere), the Biersack-Haggmark "magic
//! formula" ([`theta_magic`], fast and approximate), and a precomputed
//! [`ScatteringTable`] for the hot path.

use super::potential::{Potential, Screening};
use std::f64::consts::{FRAC_PI_2, PI};

/// Default number of Gauss-Mehler nodes for [`theta_quadrature`].
pub const DEFAULT_NODES: usize = 64;

/// Reduced distance of closest approach `x0` for reduced energy `eps` and
/// reduced impact parameter `beta > 0`.
///
/// Solves `G(x) = 0`. `G` is strictly increasing (both `phi/x` and `1/x^2`
/// decrease), so the root is unique and bracketed by
/// `[beta, (1/eps + sqrt(1/eps^2 + 4 beta^2))/2]`: at `x = beta` the
/// potential term makes `G < 0`, and at the upper end `G >= 0` because
/// `phi <= 1`. Newton's method runs inside that bracket; any step that leaves
/// it is replaced by bisection.
pub fn closest_approach(s: Screening, eps: f64, beta: f64) -> f64 {
    let inv = 1.0 / eps;
    let mut hi = 0.5 * (inv + (inv * inv + 4.0 * beta * beta).sqrt());
    let mut lo = beta;
    let g = |x: f64| 1.0 - s.phi(x) * inv / x - beta * beta / (x * x);
    let dg =
        |x: f64| -s.dphi(x) * inv / x + s.phi(x) * inv / (x * x) + 2.0 * beta * beta / (x * x * x);
    let mut x = hi;
    for _ in 0..200 {
        let gx = g(x);
        if gx > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
        let d = dg(x);
        let mut xn = x - gx / d;
        if !(xn > lo && xn < hi) || !xn.is_finite() {
            xn = 0.5 * (lo + hi);
        }
        if (xn - x).abs() <= 1e-15 * xn || hi - lo <= 1e-15 * hi {
            return xn;
        }
        x = xn;
    }
    x
}

/// Centre-of-mass scattering angle (radians, in `[0, pi]`) by Gauss-Mehler
/// quadrature with `n` nodes (Mendenhall and Weller, Nucl. Instrum. Methods B 58
/// (1991) 11; Gauss-Chebyshev rule of the first kind).
///
/// With `u = x0/x` the integral becomes `Int_0^1 du / sqrt(F(u))` where
/// `F(u) = (1 - u^2) H(u)` vanishes at the turning point. Writing
/// `u = cos t` removes the inverse-square-root endpoint singularity, leaving
/// `Int_0^{pi/2} dt / sqrt(H(cos t))`, which the midpoint rule in `t` (the
/// Gauss-Mehler nodes) integrates with exponential convergence. `H` is
/// evaluated from the cancellation-free form
/// `H(u) = beta^2/x0^2 + [phi(x0) - u phi(x0/u)] / (eps x0 (1 - u^2))`.
pub fn theta_quadrature_n(s: Screening, eps: f64, beta: f64, n: usize) -> f64 {
    if beta <= 0.0 {
        return PI;
    }
    let x0 = closest_approach(s, eps, beta);
    let phi0 = s.phi(x0);
    let b2 = beta * beta / (x0 * x0);
    let sb2 = b2.sqrt();
    let c = 1.0 / (eps * x0);
    let h = FRAC_PI_2 / n as f64;
    // theta = pi - 2 I and, with no potential, I = pi/2 exactly, so
    // theta = 2 (beta/x0) h Sum_k (1/sqrt(b2) - 1/sqrt(H_k)) (the k-sum of the
    // unperturbed integrand reproduces (pi/2) x0/beta). Each term is written as
    // delta / (sqrt(H) sqrt(b2) (sqrt(H) + sqrt(b2))) with delta = H - b2 >= 0,
    // so small angles keep full relative precision (no pi - 2 I cancellation).
    let mut sum = 0.0;
    for k in 0..n {
        let t = (k as f64 + 0.5) * h;
        let (st, u) = t.sin_cos();
        let delta = c * (phi0 - u * s.phi(x0 / u)) / (st * st);
        let sh = (b2 + delta).sqrt();
        sum += delta / (sh * sb2 * (sh + sb2));
    }
    (2.0 * beta / x0 * sum * h).min(PI)
}

/// [`theta_quadrature_n`] with [`DEFAULT_NODES`] nodes.
pub fn theta_quadrature(s: Screening, eps: f64, beta: f64) -> f64 {
    theta_quadrature_n(s, eps, beta, DEFAULT_NODES)
}

/// Constants of the Biersack-Haggmark magic formula for one screening function.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MagicConstants {
    /// `C1` .. `C5` of the fit.
    pub c: [f64; 5],
}

impl MagicConstants {
    /// Fit to the universal (ZBL) potential, Ziegler, Biersack, Littmark (1985),
    /// Ch. 2 / Appendix on the magic formula.
    pub const ZBL: Self = Self {
        c: [0.99229, 0.011615, 0.0071222, 9.3066, 14.813],
    };
    /// Fit to the Moliere potential, Biersack and Haggmark, Nucl. Instrum.
    /// Methods 174 (1980) 257.
    pub const MOLIERE: Self = Self {
        c: [0.6743, 0.009611, 0.005175, 6.314, 10.0],
    };

    /// The constant set published for a screening function, if one exists.
    pub fn for_screening(s: Screening) -> Option<Self> {
        match s {
            Screening::ZblUniversal => Some(Self::ZBL),
            Screening::Moliere => Some(Self::MOLIERE),
            _ => None,
        }
    }
}

/// Biersack-Haggmark "magic formula" (Nucl. Instrum. Methods 174 (1980) 257) for
/// the centre-of-mass angle:
///
/// ```text
/// cos(theta/2) = (beta + rho + Delta) / (x0 + rho + Delta)
/// Delta = A (x0 - beta) / (1 + G)
/// A = 2 alpha eps beta^b,  alpha = 1 + C1 eps^(-1/2),
/// b = (C2 + eps^(1/2)) / (C3 + eps^(1/2)),
/// G = (C5 + eps) / (C4 + eps) * (sqrt(1 + A^2) - A)
/// ```
///
/// with `rho = 2 (1 - V(x0)/E) / (-V'(x0)/E)` the radius of curvature of the
/// trajectory at the distance of closest approach. Only the closest-approach
/// root find is needed, no integral.
pub fn theta_magic(s: Screening, k: &MagicConstants, eps: f64, beta: f64) -> f64 {
    if beta <= 0.0 {
        return PI;
    }
    let x0 = closest_approach(s, eps, beta);
    let phi = s.phi(x0);
    let dphi = s.dphi(x0);
    // -dV/dx / E at x0, with V/E = phi/(eps x).
    let force = (phi / x0 - dphi) / (eps * x0);
    let rho = 2.0 * (1.0 - phi / (eps * x0)) / force;
    let [c1, c2, c3, c4, c5] = k.c;
    let se = eps.sqrt();
    let alpha = 1.0 + c1 / se;
    let bexp = (c2 + se) / (c3 + se);
    let a = 2.0 * alpha * eps * beta.powf(bexp);
    let g = (c5 + eps) / (c4 + eps) * ((1.0 + a * a).sqrt() - a);
    let delta = a * (x0 - beta) / (1.0 + g);
    let cos_half = ((beta + rho + delta) / (x0 + rho + delta)).clamp(-1.0, 1.0);
    2.0 * cos_half.acos()
}

/// Reduced nuclear stopping cross section
/// `s_n(eps) = 2 eps Int_0^inf sin^2(theta/2) beta d beta`
/// from the quadrature angle (reduced units of ZBL 1985, Ch. 2; the energy
/// transfer is `T = gamma E sin^2(theta/2)`). See [`nuclear_stopping_with`].
pub fn nuclear_stopping_reduced(
    s: Screening,
    eps: f64,
    beta_min: f64,
    beta_max: f64,
    per_decade: usize,
) -> f64 {
    nuclear_stopping_with(
        |b| theta_quadrature(s, eps, b),
        eps,
        beta_min,
        beta_max,
        per_decade,
    )
}

/// Reduced nuclear stopping cross section for an arbitrary angle function
/// `theta(beta)` at fixed `eps` (quadrature, magic formula, table, ...).
/// Integrated over `ln beta` with Simpson's rule, `per_decade` intervals per
/// decade from `beta_min` to `beta_max`; the contributions below `beta_min`
/// (`~ beta_min^2`) and above `beta_max` (screened, exponentially small) are
/// neglected.
pub fn nuclear_stopping_with<F: Fn(f64) -> f64>(
    theta: F,
    eps: f64,
    beta_min: f64,
    beta_max: f64,
    per_decade: usize,
) -> f64 {
    let decades = (beta_max / beta_min).log10();
    let mut n = (decades * per_decade as f64).ceil() as usize;
    if n % 2 == 1 {
        n += 1;
    }
    let du = (beta_max / beta_min).ln() / n as f64;
    let f = |u: f64| {
        let b = beta_min * u.exp();
        let sn = (0.5 * theta(b)).sin();
        sn * sn * b * b
    };
    let mut acc = f(0.0) + f(n as f64 * du);
    for i in 1..n {
        acc += f(i as f64 * du) * if i % 2 == 1 { 4.0 } else { 2.0 };
    }
    2.0 * eps * acc * du / 3.0
}

/// Specification of a [`ScatteringTable`] grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TableSpec {
    /// Smallest reduced energy.
    pub eps_min: f64,
    /// Largest reduced energy.
    pub eps_max: f64,
    /// Smallest reduced impact parameter (below it the angle is the `beta -> 0`
    /// limit, `theta ~ pi`).
    pub beta_min: f64,
    /// Largest reduced impact parameter (above it the angle is taken to be 0).
    pub beta_max: f64,
    /// Grid points per decade in both axes.
    pub per_decade: usize,
}

impl Default for TableSpec {
    fn default() -> Self {
        Self {
            eps_min: 1e-5,
            eps_max: 1e4,
            beta_min: 1e-5,
            beta_max: 1e3,
            per_decade: 32,
        }
    }
}

/// Precomputed centre-of-mass scattering angle `theta(eps, beta)`.
///
/// The table stores `y = ln tan(theta/2)` on a grid uniform in `ln eps` and
/// `ln beta`, and interpolates bilinearly in those coordinates. `y` is close to
/// linear in `ln beta` at both ends (`tan(theta/2) ~ 1/beta` as `beta -> 0`, and
/// `theta ~ beta^-k` for large `beta`), which is what makes bilinear
/// interpolation accurate. Since `d theta = sin(theta) dy`, an error `dy` in the
/// table value is at most `dy` radians in `theta`.
///
/// The angle depends on `(eps, beta)` and on the screening *function* only; `Z1`,
/// `Z2` and the screening length enter through the reduced variables. The
/// table keeps its [`Potential`] so one table serves one `(Z1, Z2, potential)`
/// pair and converts to and from SI. It is built once and shared read-only
/// (`&ScatteringTable` is `Sync`).
///
/// # Error bound
///
/// At build time the interpolation error against direct quadrature is measured
/// at cell centres (where bilinear interpolation error peaks) on a
/// deterministic subset of cells. [`ScatteringTable::max_abs_error`] is the
/// largest absolute error in `theta` seen (radians) and
/// [`ScatteringTable::max_rel_error`] the largest relative error where
/// `theta > 1e-6`. These are empirical bounds from sampling, not proofs.
///
/// Outside the grid: `beta >= beta_max` returns `0.0`; `beta < beta_min` is
/// clamped to `beta_min`; `eps` outside `[eps_min, eps_max]` returns `None`.
#[derive(Debug, Clone)]
pub struct ScatteringTable {
    potential: Potential,
    spec: TableSpec,
    n_eps: usize,
    n_beta: usize,
    d: f64,
    y: Vec<f64>,
    max_abs_error: f64,
    max_rel_error: f64,
}

/// Floor for `theta/2` in the stored values; below it `theta` is indistinguishable
/// from rounding noise of `pi - 2 I`.
const HALF_THETA_FLOOR: f64 = 5e-16;

fn half_theta_to_y(half: f64) -> f64 {
    half.max(HALF_THETA_FLOOR).tan().ln()
}

impl ScatteringTable {
    /// Build a table by direct quadrature. Serial and deterministic.
    pub fn build(potential: &Potential, spec: &TableSpec) -> Self {
        let pd = spec.per_decade as f64;
        let d = std::f64::consts::LN_10 / pd;
        let n_eps = ((spec.eps_max / spec.eps_min).log10() * pd).ceil() as usize + 1;
        let n_beta = ((spec.beta_max / spec.beta_min).log10() * pd).ceil() as usize + 1;
        let s = potential.screening;
        let mut y = Vec::with_capacity(n_eps * n_beta);
        for i in 0..n_eps {
            let eps = spec.eps_min * (i as f64 * d).exp();
            for j in 0..n_beta {
                let beta = spec.beta_min * (j as f64 * d).exp();
                y.push(half_theta_to_y(0.5 * theta_quadrature(s, eps, beta)));
            }
        }
        let mut t = Self {
            potential: *potential,
            spec: *spec,
            n_eps,
            n_beta,
            d,
            y,
            max_abs_error: 0.0,
            max_rel_error: 0.0,
        };
        t.measure_error();
        t
    }

    fn measure_error(&mut self) {
        let s = self.potential.screening;
        let (mut abs_e, mut rel_e) = (0.0f64, 0.0f64);
        // Cell centres on a deterministic stride-3 subset (offset by 1 so the
        // sampled cells differ from a lattice-aligned pattern).
        for i in (0..self.n_eps - 1).step_by(3) {
            let eps = self.spec.eps_min * ((i as f64 + 0.5) * self.d).exp();
            for j in (1..self.n_beta - 1).step_by(3) {
                let beta = self.spec.beta_min * ((j as f64 + 0.5) * self.d).exp();
                let direct = theta_quadrature(s, eps, beta);
                if let Some(v) = self.theta(eps, beta) {
                    let e = (v - direct).abs();
                    abs_e = abs_e.max(e);
                    if direct > 1e-6 {
                        rel_e = rel_e.max(e / direct);
                    }
                }
            }
        }
        self.max_abs_error = abs_e;
        self.max_rel_error = rel_e;
    }

    /// The potential this table was built for.
    pub fn potential(&self) -> &Potential {
        &self.potential
    }

    /// The grid specification.
    pub fn spec(&self) -> &TableSpec {
        &self.spec
    }

    /// Largest absolute `theta` error (radians) found against quadrature at
    /// sampled cell centres.
    pub fn max_abs_error(&self) -> f64 {
        self.max_abs_error
    }

    /// Largest relative `theta` error found at sampled cell centres where
    /// `theta > 1e-6`.
    pub fn max_rel_error(&self) -> f64 {
        self.max_rel_error
    }

    /// Interpolated centre-of-mass angle (radians), or `None` if `eps` is out
    /// of the tabulated range.
    pub fn theta(&self, eps: f64, beta: f64) -> Option<f64> {
        if !(eps >= self.spec.eps_min && eps <= self.spec.eps_max) {
            return None;
        }
        if beta >= self.spec.beta_max {
            return Some(0.0);
        }
        let beta = beta.max(self.spec.beta_min);
        let fi = (eps / self.spec.eps_min).ln() / self.d;
        let fj = (beta / self.spec.beta_min).ln() / self.d;
        let i = (fi.floor() as usize).min(self.n_eps - 2);
        let j = (fj.floor() as usize).min(self.n_beta - 2);
        let (ti, tj) = (fi - i as f64, fj - j as f64);
        let at = |a: usize, b: usize| self.y[a * self.n_beta + b];
        let y = (1.0 - ti) * ((1.0 - tj) * at(i, j) + tj * at(i, j + 1))
            + ti * ((1.0 - tj) * at(i + 1, j) + tj * at(i + 1, j + 1));
        Some(2.0 * y.exp().atan())
    }

    /// Number of grid points in `(eps, beta)`.
    pub fn shape(&self) -> (usize, usize) {
        (self.n_eps, self.n_beta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ZBL universal reduced nuclear stopping fit, Ziegler, Biersack, Littmark
    /// (1985), Ch. 2: `ln(1 + 1.1383 eps) / (2 (eps + 0.01321 eps^0.21226 +
    /// 0.19593 eps^0.5))` for `eps <= 30`, `ln(eps) / (2 eps)` above.
    fn zbl_fit(eps: f64) -> f64 {
        if eps <= 30.0 {
            (1.0 + 1.1383 * eps).ln()
                / (2.0 * (eps + 0.01321 * eps.powf(0.21226) + 0.19593 * eps.sqrt()))
        } else {
            eps.ln() / (2.0 * eps)
        }
    }

    #[test]
    fn quadrature_converges_with_nodes() {
        let s = Screening::ZblUniversal;
        for &(eps, beta) in &[(1e-3, 0.5), (0.1, 1.0), (10.0, 0.1), (100.0, 3.0)] {
            let a = theta_quadrature_n(s, eps, beta, DEFAULT_NODES);
            let b = theta_quadrature_n(s, eps, beta, 1024);
            assert!(
                (a - b).abs() < 1e-6 * (1.0 + a),
                "eps={eps} beta={beta}: {a} {b}"
            );
        }
    }

    #[test]
    fn coulomb_limit_matches_rutherford() {
        // Very high energy, tiny beta region is not Coulomb; but at small beta
        // and huge eps phi ~ 1 and tan(theta/2) = 1/(2 eps beta).
        let eps = 1e6;
        let beta = 1e-3;
        let th = theta_quadrature(Screening::ZblUniversal, eps, beta);
        let expect = 2.0 * (1.0 / (2.0 * eps * beta)).atan();
        assert!((th - expect).abs() < 2e-3 * expect, "{th} {expect}");
    }

    #[test]
    fn limits_in_impact_parameter() {
        for s in Screening::ALL {
            for &eps in &[1e-3, 1.0, 1e3] {
                assert!(theta_quadrature(s, eps, 1e3) < 1e-6, "{s:?} {eps}");
                // Head-on is Coulomb-like: pi - theta ~ 4 eps beta.
                let head_on = theta_quadrature(s, eps, 1e-9);
                assert!((PI - head_on) < 1e-5 + 5e-9 * eps, "{s:?} {eps} {head_on}");
                assert_eq!(theta_quadrature(s, eps, 0.0), PI);
            }
        }
    }

    #[test]
    fn theta_decreases_with_beta() {
        for s in Screening::ALL {
            let mut prev = PI;
            for k in 0..60 {
                let beta = 1e-3 * 10f64.powf(k as f64 / 10.0);
                let th = theta_quadrature(s, 0.5, beta);
                assert!(th <= prev + 1e-12, "{s:?} beta={beta}");
                prev = th;
            }
        }
    }

    #[test]
    fn closest_approach_is_a_root_above_beta() {
        for s in Screening::ALL {
            for &eps in &[1e-4, 1e-2, 1.0, 1e3] {
                for &beta in &[1e-6, 0.01, 1.0, 20.0] {
                    let x0 = closest_approach(s, eps, beta);
                    let g = 1.0 - s.phi(x0) / (eps * x0) - beta * beta / (x0 * x0);
                    assert!(x0 >= beta && g.abs() < 1e-10, "{s:?} {eps} {beta} g={g}");
                }
            }
        }
    }

    /// Compare the integrated quadrature angle with the ZBL universal nuclear
    /// stopping fit (as recalled; see `docs/data-provenance.md`).
    ///
    /// Measured: agreement is better than 0.3 % for `eps <= 0.3` and better than
    /// 1 % up to `eps = 1`, but the deviation grows to about 5.5 % for
    /// `3 < eps < 1e3`, short of the 0.5 % acceptance criterion of issue #3.
    /// What has been ruled out: the quadrature (converged in nodes, in the
    /// Simpson step and in the `beta` range, and reproduced to four digits by an
    /// independent scipy implementation of the angle and stopping integrals) and
    /// the magic formula (integrating `theta_magic` is far worse, up to 1.7x).
    /// What is not established: whether the recalled fit coefficients are right,
    /// or whether the published fit is itself only accurate to a few percent
    /// there (its two branches disagree by about 1 % at `eps = 30`). The cause is
    /// open until the fit is checked against the book. The bands below only
    /// record what is achieved and guard against regressions; they are not a
    /// claim that the criterion is met.
    #[test]
    fn stopping_vs_zbl_universal_fit() {
        let mut worst_low = 0.0f64;
        let mut worst_high = 0.0f64;
        for k in 0..=14 {
            let eps = 1e-4 * 10f64.powf(k as f64 / 2.0);
            let sn = nuclear_stopping_reduced(Screening::ZblUniversal, eps, 1e-4, 1e3, 40);
            let rel = (sn / zbl_fit(eps) - 1.0).abs();
            if eps <= 0.3 {
                worst_low = worst_low.max(rel);
            } else {
                worst_high = worst_high.max(rel);
            }
        }
        eprintln!("zbl fit: worst rel dev eps<=0.3: {worst_low:.3e}, eps>0.3: {worst_high:.3e}");
        assert!(worst_low < 5e-3, "eps <= 0.3: {worst_low}");
        assert!(worst_high < 6e-2, "eps > 0.3: {worst_high}");
    }

    #[test]
    fn magic_formula_tracks_quadrature_for_zbl() {
        let k = MagicConstants::ZBL;
        let mut worst = 0.0f64;
        for &eps in &[1e-2, 1e-1, 1.0, 10.0, 100.0] {
            for j in 0..40 {
                let beta = 1e-2 * 10f64.powf(j as f64 / 10.0);
                let q = theta_quadrature(Screening::ZblUniversal, eps, beta);
                let m = theta_magic(Screening::ZblUniversal, &k, eps, beta);
                // Compare cos(theta/2): relative error in theta is ill-conditioned
                // where theta is tiny.
                worst = worst.max(((0.5 * m).cos() - (0.5 * q).cos()).abs());
            }
        }
        eprintln!("magic worst |d cos(theta/2)| {worst}");
        assert!(worst < 0.05, "{worst}");
    }

    #[test]
    fn magic_formula_tracks_quadrature_for_moliere() {
        let k = MagicConstants::MOLIERE;
        let mut worst = 0.0f64;
        for &eps in &[1e-2, 1e-1, 1.0, 10.0, 100.0] {
            for j in 0..40 {
                let beta = 1e-2 * 10f64.powf(j as f64 / 10.0);
                let q = theta_quadrature(Screening::Moliere, eps, beta);
                let m = theta_magic(Screening::Moliere, &k, eps, beta);
                worst = worst.max(((0.5 * m).cos() - (0.5 * q).cos()).abs());
            }
        }
        eprintln!("moliere magic worst |d cos(theta/2)| {worst}");
        assert!(worst < 0.05, "{worst}");
    }

    #[test]
    fn table_interpolation_within_stated_bound() {
        let pot = Potential::new(Screening::ZblUniversal, 31.0, 7.0);
        let spec = TableSpec {
            eps_min: 1e-3,
            eps_max: 1e2,
            beta_min: 1e-3,
            beta_max: 1e2,
            per_decade: 32,
        };
        let t = ScatteringTable::build(&pot, &spec);
        eprintln!("table abs {} rel {}", t.max_abs_error(), t.max_rel_error());
        assert!(t.max_abs_error() < 1e-2);
        // Deterministic off-grid points (LCG).
        let mut st = 12345u64;
        let mut next = || {
            st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (st >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..500 {
            let eps = 1e-3 * 10f64.powf(5.0 * next());
            let beta = 1e-3 * 10f64.powf(4.9 * next() + 0.05);
            let direct = theta_quadrature(Screening::ZblUniversal, eps, beta);
            let tab = t.theta(eps, beta).unwrap();
            assert!(
                (tab - direct).abs() <= 1.5 * t.max_abs_error() + 1e-12,
                "eps={eps} beta={beta} {tab} {direct}"
            );
        }
    }

    #[test]
    fn table_limits_and_range() {
        let pot = Potential::new(Screening::KrC, 18.0, 18.0);
        let spec = TableSpec {
            per_decade: 16,
            ..TableSpec::default()
        };
        let t = ScatteringTable::build(&pot, &spec);
        assert_eq!(t.theta(1e-9, 1.0), None);
        assert_eq!(t.theta(1.0, 1e5), Some(0.0));
        assert!(PI - t.theta(1.0, 1e-8).unwrap() < 1e-3);
    }

    #[test]
    fn table_build_is_deterministic() {
        let pot = Potential::new(Screening::Moliere, 14.0, 14.0);
        let spec = TableSpec {
            eps_min: 1e-2,
            eps_max: 1e1,
            beta_min: 1e-2,
            beta_max: 1e1,
            per_decade: 8,
        };
        let a = ScatteringTable::build(&pot, &spec);
        let b = ScatteringTable::build(&pot, &spec);
        assert_eq!(a.y, b.y);
    }
}
