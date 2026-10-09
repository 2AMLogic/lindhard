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

    /// `∫_a^b f` (one application of the rule) and the values of `f` at the
    /// nodes, in node order.
    fn integrate_with_values(
        &self,
        f: &mut impl FnMut(f64) -> f64,
        a: f64,
        b: f64,
    ) -> (f64, Vec<f64>) {
        let half = 0.5 * (b - a);
        let mid = 0.5 * (a + b);
        let mut acc = 0.0;
        let mut values = Vec::with_capacity(self.nodes.len());
        for (x, w) in self.nodes.iter().zip(&self.weights) {
            let v = f(mid + half * x);
            acc += w * v;
            values.push(v);
        }
        (acc * half, values)
    }

    /// `∫_{-1}^{τ} P(t) dt`, `-1 <= τ <= 1`, of the polynomial `P` of degree
    /// `n - 1` that interpolates `values` at the nodes (Lagrange form). The
    /// rule integrates `P` exactly on `[-1, τ]` (degree `n - 1 <= 2n - 1`).
    fn interpolant_integral(&self, values: &[f64], tau: f64) -> f64 {
        let half = 0.5 * (tau + 1.0);
        let mut acc = 0.0;
        for (t, w) in self.nodes.iter().zip(&self.weights) {
            let x = -1.0 + half * (1.0 + t);
            let mut p = 0.0;
            for (k, (tk, vk)) in self.nodes.iter().zip(values).enumerate() {
                let mut l = 1.0;
                for (j, tj) in self.nodes.iter().enumerate() {
                    if j != k {
                        l *= (x - tj) / (tk - tj);
                    }
                }
                p += vk * l;
            }
            acc += w * p;
        }
        acc * half
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

/// The accepted panels of one adaptive integration ([`integrate_panels`]),
/// kept so that the integral over any sub-interval can be taken later
/// without new evaluations of the integrand.
///
/// On every panel the integrand is represented by the polynomial of degree
/// `n - 1` through its values at the `n` Gauss-Legendre nodes, the
/// interpolant whose integral over the whole panel is the panel's rule. The
/// integral up to a point inside a panel is that polynomial's integral
/// (exact, see [`GaussLegendre::interpolant_integral`]); whole panels count
/// with their rule values, so the full range returns the adaptive integral
/// itself.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Panels {
    /// Panel edges, ascending: panel `i` is `[edges[i], edges[i + 1]]`.
    edges: Vec<f64>,
    /// `n` node values per panel.
    values: Vec<f64>,
    /// `cum[i]`: the integral over panels `0..i`.
    cum: Vec<f64>,
}

impl Panels {
    /// The integral from the lower end to `x`, `x` clamped to the range.
    pub(crate) fn integral_to(&self, gl: &GaussLegendre, x: f64) -> f64 {
        let np = self.edges.len() - 1;
        if np == 0 || x.is_nan() || x <= self.edges[0] {
            return 0.0;
        }
        if x >= self.edges[np] {
            return self.cum[np];
        }
        // panel i with edges[i] <= x < edges[i + 1]
        let i = (self.edges.partition_point(|&e| e <= x) - 1).min(np - 1);
        let (a, b) = (self.edges[i], self.edges[i + 1]);
        let n = gl.nodes.len();
        let tau = (2.0 * (x - a) / (b - a) - 1.0).clamp(-1.0, 1.0);
        let part = 0.5 * (b - a) * gl.interpolant_integral(&self.values[i * n..(i + 1) * n], tau);
        self.cum[i] + part
    }

    /// The integral over `[lo, hi]` (clamped to the range; zero if
    /// `hi <= lo`).
    pub(crate) fn integral(&self, gl: &GaussLegendre, lo: f64, hi: f64) -> f64 {
        if hi.is_nan() || lo.is_nan() || hi <= lo {
            return 0.0;
        }
        self.integral_to(gl, hi) - self.integral_to(gl, lo)
    }

    /// Number of panels.
    pub(crate) fn len(&self) -> usize {
        self.edges.len() - 1
    }
}

/// [`refine`] for one integrand, recording the accepted panels (the two
/// halves of every leaf) in order.
fn refine_recorded(
    gl: &GaussLegendre,
    f: &mut impl FnMut(f64) -> f64,
    (a, b): (f64, f64),
    whole: f64,
    tol: f64,
    depth: u32,
    out: &mut Vec<(f64, f64, f64, Vec<f64>)>,
) -> f64 {
    let m = 0.5 * (a + b);
    let (left, lv) = gl.integrate_with_values(f, a, m);
    let (right, rv) = gl.integrate_with_values(f, m, b);
    let both = left + right;
    if (both - whole).abs() <= tol || depth == 0 || !(m > a && m < b) {
        out.push((a, m, left, lv));
        out.push((m, b, right, rv));
        return both;
    }
    let half_tol = 0.5 * tol;
    let l = refine_recorded(gl, f, (a, m), left, half_tol, depth - 1, out);
    let r = refine_recorded(gl, f, (m, b), right, half_tol, depth - 1, out);
    l + r
}

/// [`integrate_segments`] for one integrand, returning the accepted panels
/// instead of only their sum. `breaks` must be strictly increasing. The
/// tolerance of segment `i` is `rel_tol · max(|I_i|, floor · Σ|I|/N)`: with
/// `floor = 1` this is exactly the error control (and so the integral, up to
/// the order of the final sum) of [`integrate_segments`]; a smaller `floor`
/// holds every segment to `rel_tol` of its own integral down to that
/// fraction of the mean, so that a sub-range made of a few segments is as
/// accurate, relative to itself, as the whole.
pub(crate) fn integrate_panels(
    gl: &GaussLegendre,
    f: &mut impl FnMut(f64) -> f64,
    breaks: &[f64],
    rel_tol: f64,
    floor: f64,
) -> Panels {
    let mut leaves: Vec<(f64, f64, f64, Vec<f64>)> = Vec::new();
    if breaks.len() >= 2 {
        let wholes: Vec<f64> = breaks
            .windows(2)
            .map(|w| gl.integrate(&mut |x: f64| [f(x)], w[0], w[1])[0])
            .collect();
        let scale: f64 = wholes.iter().map(|x| x.abs()).sum();
        let per_segment = floor * scale / wholes.len() as f64;
        for (w, whole) in breaks.windows(2).zip(&wholes) {
            let tol = rel_tol * whole.abs().max(per_segment);
            refine_recorded(gl, f, (w[0], w[1]), *whole, tol, MAX_DEPTH, &mut leaves);
        }
    }
    let mut edges = Vec::with_capacity(leaves.len() + 1);
    let mut values = Vec::with_capacity(leaves.len() * gl.nodes.len());
    let mut cum = Vec::with_capacity(leaves.len() + 1);
    let mut acc = 0.0;
    for (a, _, integral, v) in &leaves {
        edges.push(*a);
        cum.push(acc);
        acc += integral;
        values.extend_from_slice(v);
    }
    edges.push(
        leaves
            .last()
            .map_or(breaks.first().copied().unwrap_or(0.0), |l| l.1),
    );
    cum.push(acc);
    Panels { edges, values, cum }
}

/// Number of consecutive segments a smoothness check of
/// [`integrate_weighted`] spans at first.
const GROUP: usize = 16;
/// Fraction of the tolerance the smoothness check of a group may use
/// (see [`integrate_weighted`]).
const GROUP_FRACTION: f64 = 0.1;

/// `∫ w(x) k(x) dx` over consecutive segments `[breaks[i], breaks[i + 1]]`,
/// where the weight `w` is non-negative and *linear* on every segment (the
/// breaks are the knots of a piecewise-linear table, times a smooth factor
/// absorbed in `k`) and the kernel `k` is smooth except at the breaks.
///
/// This is [`integrate_segments`] with a cheaper error control when there are
/// many segments and the kernel varies slowly over most of them, which is the
/// case of the full Penn integrals over the knots of a measured ELF. Every
/// segment gets one `n`-point rule (`n` evaluations of `k`). Then the
/// *kernel alone* is checked on groups of up to [`GROUP`] consecutive
/// segments: one `n`-point rule on the whole group is compared with the sum
/// of the rules on its segments, which are already known. If they differ by
/// no more than [`GROUP_FRACTION`] of the tolerance `rel_tol · max(Σ|I_s|,
/// (group size) · Σ|I|/N)` (times the largest weight in the group, since
/// `w` multiplies the error), the group is accepted with the per-segment
/// rules. The group rule has the error of the rule on a span `G` times
/// wider, so it is orders of magnitude less accurate than the per-segment
/// sum, and agreement within a fraction of the tolerance means that the
/// per-segment errors are far below it. The comparison is of signed
/// integrals, so in principle errors of opposite sign on the segments of a
/// group can cancel and let a group pass that a check of each segment would
/// fail. This is the trade-off for evaluating the kernel only once more per
/// group. On the full Penn integrals it did not show: with this routine in
/// the `ω_p` integral, the Al inverse IMFP and stopping power at the 57
/// energies of the #169 grid stayed within 3e-8 and 2e-7 relative of the
/// earlier per-segment bisection (#256, `docs/validation.md`). A group that
/// fails is halved, down
/// to single segments, and a single failing segment is refined by the
/// bisection of [`integrate_segments`] with the same tolerance, starting from
/// its first estimate. Segments on which `w` vanishes at both ends are
/// skipped, as `k` is not evaluated there.
///
/// Every operation is a fixed sequence of floating-point operations, so the
/// result does not depend on the thread that computes it.
pub(crate) fn integrate_weighted<const N: usize>(
    gl: &GaussLegendre,
    weight: &impl Fn(f64) -> f64,
    kernel: &mut impl FnMut(f64) -> [f64; N],
    breaks: &[f64],
    rel_tol: f64,
) -> [f64; N] {
    if breaks.len() < 2 {
        return [0.0; N];
    }
    let nseg = breaks.len() - 1;
    // per-segment first estimates of the integral and of the kernel alone
    let mut whole = vec![[0.0; N]; nseg];
    let mut kernel_whole = vec![[0.0; N]; nseg];
    let mut active = vec![false; nseg];
    let at_breaks: Vec<f64> = breaks.iter().map(|&b| weight(b)).collect();
    for (i, w) in breaks.windows(2).enumerate() {
        if w[1] <= w[0] || (at_breaks[i] <= 0.0 && at_breaks[i + 1] <= 0.0) {
            continue;
        }
        active[i] = true;
        let half = 0.5 * (w[1] - w[0]);
        let mid = 0.5 * (w[0] + w[1]);
        for (x, wt) in gl.nodes.iter().zip(&gl.weights) {
            let xx = mid + half * x;
            let k = kernel(xx);
            let g = weight(xx);
            for c in 0..N {
                kernel_whole[i][c] += wt * k[c];
                whole[i][c] += wt * g * k[c];
            }
        }
        for c in 0..N {
            kernel_whole[i][c] *= half;
            whole[i][c] *= half;
        }
    }
    let mut scale = [0.0; N];
    for w in &whole {
        for (s, x) in scale.iter_mut().zip(w) {
            *s += x.abs();
        }
    }
    let per_segment = scale.map(|s| s / nseg as f64);
    let mut total = [0.0; N];
    let mut start = 0;
    while start < nseg {
        if !active[start] {
            start += 1;
            continue;
        }
        let mut end = start;
        while end < nseg && end < start + GROUP && active[end] {
            end += 1;
        }
        total = add(
            total,
            weighted_group(
                gl,
                weight,
                kernel,
                breaks,
                (&whole, &kernel_whole, &at_breaks),
                &per_segment,
                (start, end),
                rel_tol,
            ),
        );
        start = end;
    }
    total
}

/// Segments `lo..hi` of [`integrate_weighted`] (all active): accept, split or
/// refine.
#[allow(clippy::too_many_arguments)]
fn weighted_group<const N: usize>(
    gl: &GaussLegendre,
    weight: &impl Fn(f64) -> f64,
    kernel: &mut impl FnMut(f64) -> [f64; N],
    breaks: &[f64],
    (whole, kernel_whole, at_breaks): (&[[f64; N]], &[[f64; N]], &[f64]),
    per_segment: &[f64; N],
    (lo, hi): (usize, usize),
    rel_tol: f64,
) -> [f64; N] {
    let n = hi - lo;
    if n == 1 {
        let mut f = |x: f64| {
            let g = weight(x);
            if g == 0.0 {
                [0.0; N]
            } else {
                kernel(x).map(|k| g * k)
            }
        };
        let mut tol = [0.0; N];
        for c in 0..N {
            tol[c] = rel_tol * whole[lo][c].abs().max(per_segment[c]);
        }
        return refine(
            gl,
            &mut f,
            (breaks[lo], breaks[lo + 1]),
            whole[lo],
            tol,
            MAX_DEPTH,
        );
    }
    let group = gl.integrate(kernel, breaks[lo], breaks[hi]);
    let mut gmax: f64 = 0.0;
    for &g in &at_breaks[lo..=hi] {
        gmax = gmax.max(g);
    }
    let mut smooth = true;
    let mut sum = [0.0; N];
    for c in 0..N {
        let composite: f64 = kernel_whole[lo..hi].iter().map(|k| k[c]).sum();
        let mass: f64 = whole[lo..hi].iter().map(|w| w[c].abs()).sum();
        sum[c] = whole[lo..hi].iter().map(|w| w[c]).sum();
        let tol = GROUP_FRACTION * rel_tol * mass.max(n as f64 * per_segment[c]);
        if gmax * (group[c] - composite).abs() > tol {
            smooth = false;
        }
    }
    if smooth {
        return sum;
    }
    let mid = lo + n / 2;
    let args = (whole, kernel_whole, at_breaks);
    add(
        weighted_group(
            gl,
            weight,
            kernel,
            breaks,
            args,
            per_segment,
            (lo, mid),
            rel_tol,
        ),
        weighted_group(
            gl,
            weight,
            kernel,
            breaks,
            args,
            per_segment,
            (mid, hi),
            rel_tol,
        ),
    )
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
    fn interpolant_integral_is_exact_for_polynomials_of_degree_n_minus_1() {
        let gl = GaussLegendre::new(5);
        // P(t) = 1 + 2t - t^3 + 0.5 t^4: ∫_{-1}^{τ} P = F(τ) - F(-1)
        let p = |t: f64| 1.0 + 2.0 * t - t.powi(3) + 0.5 * t.powi(4);
        let big_f = |t: f64| t + t * t - 0.25 * t.powi(4) + 0.1 * t.powi(5);
        let values: Vec<f64> = gl.nodes.iter().map(|&t| p(t)).collect();
        for tau in [-1.0, -0.6, 0.0, 0.3, 0.99, 1.0] {
            let got = gl.interpolant_integral(&values, tau);
            let want = big_f(tau) - big_f(-1.0);
            assert!((got - want).abs() < 1e-13, "tau {tau}: {got} vs {want}");
        }
    }

    #[test]
    fn panels_give_the_adaptive_integral_and_its_partial_integrals() {
        let gl = GaussLegendre::new(5);
        // A peaked integrand with a closed-form antiderivative.
        let f = |x: f64| 1.0 / (1.0 + 100.0 * (x - 3.0) * (x - 3.0));
        let big_f = |x: f64| (10.0 * (x - 3.0)).atan() / 10.0;
        let breaks = [0.0, 2.5, 5.0, 7.5, 10.0];
        let p = integrate_panels(&gl, &mut |x: f64| f(x), &breaks, 1e-8, 1.0);
        let whole = integrate_segments(&gl, &mut |x: f64| [f(x)], &breaks, 1e-8)[0];
        // the same adaptive integral, summed in another order
        assert!((p.integral(&gl, 0.0, 10.0) / whole - 1.0).abs() < 1e-14);
        assert!(p.len() > breaks.len() - 1);
        for (lo, hi) in [
            (0.0, 10.0),
            (2.9, 3.05),
            (0.1, 9.7),
            (3.0, 3.0001),
            (-5.0, 1.0),
        ] {
            let got = p.integral(&gl, lo, hi);
            let want = big_f(hi.min(10.0)) - big_f(lo.max(0.0));
            assert!(
                (got - want).abs() < 1e-8 * (big_f(10.0) - big_f(0.0)),
                "[{lo}, {hi}]: {got} vs {want}"
            );
        }
        assert_eq!(p.integral(&gl, 4.0, 3.0), 0.0);
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
