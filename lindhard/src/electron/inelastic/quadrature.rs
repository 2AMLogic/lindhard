//! Gauss-Legendre quadrature with adaptive bisection, for the inelastic
//! integrals.
//!
//! The nodes and weights are computed at construction, not tabulated: the
//! nodes are the roots of the Legendre polynomial `P_n`, found by Newton's
//! method on the three-term (Bonnet) recurrence
//! `k P_k(x) = (2k - 1) x P_{k-1}(x) - (k - 1) P_{k-2}(x)`, and the weights are
//! `w_i = 2 / ((1 - x_i²) P_n'(x_i)²)` with
//! `P_n'(x) = n (x P_n(x) - P_{n-1}(x)) / (x² - 1)`. These are identities of
//! the Legendre polynomials, so no constant is entered by hand; the unit tests
//! check that an `n`-point rule integrates polynomials of degree `2n - 1`
//! exactly.
//!
//! # Error control
//!
//! [`integrate_segments`] integrates over a list of breakpoints (the places
//! where the integrand has a kink, e.g. the knots of a piecewise-linear
//! table). Each segment is first integrated with one `n`-point rule; the sum of
//! the absolute values of these first estimates, `Σ|I_i|`, sets the scale. A
//! segment is then bisected recursively until the two halves agree with the
//! whole to within `rel_tol · max(|I_i|, Σ|I| / N)` (the budget is halved at
//! each level). The difference between successive estimates is the error
//! estimate, so the estimated error of the total is at most
//! `2 · rel_tol · Σ|I|`. Bisection stops at depth [`MAX_DEPTH`].
//!
//! Every operation is a fixed sequence of floating-point operations, so the
//! result does not depend on the thread that computes it.

/// Maximum bisection depth of [`integrate_segments`] (a segment is never cut
/// into pieces smaller than `2^-MAX_DEPTH` of its width).
pub(crate) const MAX_DEPTH: u32 = 40;

/// An `n`-point Gauss-Legendre rule on `[-1, 1]`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GaussLegendre {
    nodes: Vec<f64>,
    weights: Vec<f64>,
}

/// `(P_n(x), P_{n-1}(x))` by the Bonnet recurrence.
fn legendre_pair(n: usize, x: f64) -> (f64, f64) {
    let (mut p_prev, mut p) = (1.0, x);
    for k in 2..=n {
        let kf = k as f64;
        let next = ((2.0 * kf - 1.0) * x * p - (kf - 1.0) * p_prev) / kf;
        p_prev = p;
        p = next;
    }
    (p, p_prev)
}

impl GaussLegendre {
    /// The `n`-point rule, `n >= 2`.
    pub(crate) fn new(n: usize) -> Self {
        assert!(n >= 2, "a Gauss-Legendre rule needs at least 2 points");
        let nf = n as f64;
        let mut nodes = Vec::with_capacity(n);
        let mut weights = Vec::with_capacity(n);
        for i in 0..n {
            // Starting guess: the i-th root lies near cos(pi (i + 3/4) / (n + 1/2)).
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
            for _ in 0..100 {
                let (p, p1) = legendre_pair(n, x);
                let dx = p / (nf * (x * p - p1) / (x * x - 1.0));
                x -= dx;
                if dx.abs() <= 1e-16 {
                    break;
                }
            }
            let (p, p1) = legendre_pair(n, x);
            let dp = nf * (x * p - p1) / (x * x - 1.0);
            nodes.push(x);
            weights.push(2.0 / ((1.0 - x * x) * dp * dp));
        }
        Self { nodes, weights }
    }

    /// `∫_a^b f`, one application of the rule, for `N` integrands at once.
    pub(crate) fn integrate<const N: usize>(
        &self,
        f: &mut impl FnMut(f64) -> [f64; N],
        a: f64,
        b: f64,
    ) -> [f64; N] {
        let half = 0.5 * (b - a);
        let mid = 0.5 * (a + b);
        let mut acc = [0.0; N];
        for (x, w) in self.nodes.iter().zip(&self.weights) {
            let v = f(mid + half * x);
            for (a, v) in acc.iter_mut().zip(v) {
                *a += w * v;
            }
        }
        acc.map(|a| a * half)
    }
}

fn add<const N: usize>(a: [f64; N], b: [f64; N]) -> [f64; N] {
    let mut out = a;
    for (o, b) in out.iter_mut().zip(b) {
        *o += b;
    }
    out
}

/// Recursive bisection of `[a, b]` whose one-rule estimate is `whole`, until
/// every component of the halves agrees with the whole to within `tol`.
fn refine<const N: usize>(
    gl: &GaussLegendre,
    f: &mut impl FnMut(f64) -> [f64; N],
    (a, b): (f64, f64),
    whole: [f64; N],
    tol: [f64; N],
    depth: u32,
) -> [f64; N] {
    let m = 0.5 * (a + b);
    let left = gl.integrate(f, a, m);
    let right = gl.integrate(f, m, b);
    let both = add(left, right);
    let converged = both
        .iter()
        .zip(&whole)
        .zip(&tol)
        .all(|((x, w), t)| (x - w).abs() <= *t);
    if converged || depth == 0 || !(m > a && m < b) {
        return both;
    }
    let half_tol = tol.map(|t| 0.5 * t);
    add(
        refine(gl, f, (a, m), left, half_tol, depth - 1),
        refine(gl, f, (m, b), right, half_tol, depth - 1),
    )
}

/// `∫ f` over consecutive segments `[breaks[i], breaks[i + 1]]`, with the
/// error control described in the module docs. `breaks` must be
/// non-decreasing; empty segments contribute zero.
pub(crate) fn integrate_segments<const N: usize>(
    gl: &GaussLegendre,
    f: &mut impl FnMut(f64) -> [f64; N],
    breaks: &[f64],
    rel_tol: f64,
) -> [f64; N] {
    if breaks.len() < 2 {
        return [0.0; N];
    }
    let wholes: Vec<[f64; N]> = breaks
        .windows(2)
        .map(|w| {
            if w[1] > w[0] {
                gl.integrate(f, w[0], w[1])
            } else {
                [0.0; N]
            }
        })
        .collect();
    let mut scale = [0.0; N];
    for w in &wholes {
        for (s, x) in scale.iter_mut().zip(w) {
            *s += x.abs();
        }
    }
    let per_segment = scale.map(|s| s / wholes.len() as f64);
    let mut total = [0.0; N];
    for (w, whole) in breaks.windows(2).zip(&wholes) {
        if w[1] <= w[0] {
            continue;
        }
        let mut tol = [0.0; N];
        for k in 0..N {
            tol[k] = rel_tol * whole[k].abs().max(per_segment[k]);
        }
        total = add(total, refine(gl, f, (w[0], w[1]), *whole, tol, MAX_DEPTH));
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rule_is_exact_for_polynomials_of_degree_2n_minus_1() {
        for n in [2usize, 5, 8, 12] {
            let gl = GaussLegendre::new(n);
            let w: f64 = gl.weights.iter().sum();
            assert!((w - 2.0).abs() < 1e-14, "n = {n}: weights sum to {w}");
            for deg in 0..(2 * n) {
                // ∫_0^2 x^deg dx = 2^(deg+1) / (deg + 1)
                let got = gl.integrate(&mut |x: f64| [x.powi(deg as i32)], 0.0, 2.0)[0];
                let want = 2f64.powi(deg as i32 + 1) / (deg as f64 + 1.0);
                assert!(
                    ((got - want) / want).abs() < 1e-13,
                    "n = {n}, degree {deg}: {got} vs {want}"
                );
            }
        }
    }

    #[test]
    fn adaptive_integration_meets_its_tolerance_across_a_kink() {
        let gl = GaussLegendre::new(8);
        // |x - 0.3| on [0, 1] with no breakpoint at the kink: 0.045 + 0.245.
        let got = integrate_segments(&gl, &mut |x: f64| [(x - 0.3).abs()], &[0.0, 1.0], 1e-10)[0];
        assert!((got - 0.29).abs() < 1e-9, "{got}");
        // A peaked integrand: ∫_0^10 1/(1 + 100 (x-3)^2) dx = (atan(70)+atan(30))/10.
        let want = (70f64.atan() + 30f64.atan()) / 10.0;
        let got = integrate_segments(
            &gl,
            &mut |x: f64| [1.0 / (1.0 + 100.0 * (x - 3.0) * (x - 3.0))],
            &[0.0, 10.0],
            1e-9,
        )[0];
        assert!(((got - want) / want).abs() < 1e-8, "{got} vs {want}");
    }
}
