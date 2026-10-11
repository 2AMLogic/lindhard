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
    /// Cells of the bisection accepted without meeting their tolerance
    /// (see [`integrate_panels`]).
    unconverged: u64,
    /// The sum of the error estimates of those cells.
    unconverged_error: f64,
    /// Evaluations of the integrand whose own integral did not converge,
    /// as recorded by the caller (see [`Self::with_inner_unconverged`]).
    inner_unconverged: u64,
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

    /// Records that `n` evaluations of the integrand were themselves
    /// integrals that did not converge (#298).
    pub(crate) fn with_inner_unconverged(mut self, n: u64) -> Self {
        self.inner_unconverged = n;
        self
    }

    /// Evaluations of the integrand that did not converge, as recorded by
    /// [`Self::with_inner_unconverged`] (zero if never recorded).
    pub(crate) fn inner_unconverged(&self) -> u64 {
        self.inner_unconverged
    }

    /// Number of panels.
    pub(crate) fn len(&self) -> usize {
        self.edges.len() - 1
    }

    /// Number of cells of the bisection that were accepted without meeting
    /// their tolerance (zero when the integration converged).
    pub(crate) fn unconverged_panels(&self) -> u64 {
        self.unconverged
    }

    /// The sum of the error estimates of the cells counted in
    /// [`Self::unconverged_panels`], relative to the magnitude of the whole
    /// integral (zero when it converged, or when the integral is zero).
    pub(crate) fn unconverged_relative_error(&self) -> f64 {
        let total = self.cum[self.cum.len() - 1].abs();
        if self.unconverged_error > 0.0 && total > 0.0 {
            self.unconverged_error / total
        } else {
            0.0
        }
    }
}

/// The bisection of one segment of [`integrate_panels`], recording the
/// accepted panels in order. `[a, b]` has the one-rule estimate `whole.0`
/// with node values `whole.1`. A cell is accepted, and its two halves
/// recorded, when they agree with the whole to within the larger of `tol`
/// (halved at each level) and `rel_tol` times their own integral: the rule
/// and the bound of [`integrate_weighted`] (#298), so that noise in the
/// integrand near a jump or a narrow peak cannot drive the bisection to
/// [`MAX_DEPTH`] over a whole neighbourhood. `remaining` is the segment's
/// budget of rule applications, spent depth first; a cell that cannot be
/// bisected (budget spent, depth 0, midpoint not representable) is recorded
/// as it stands and counted in `unconverged`, as in [`refine_bounded`].
#[allow(clippy::too_many_arguments)]
fn refine_recorded(
    gl: &GaussLegendre,
    f: &mut impl FnMut(f64) -> f64,
    (a, b): (f64, f64),
    whole: (f64, Vec<f64>),
    (tol, rel_tol): (f64, f64),
    depth: u32,
    (remaining, parent_err): (&mut u64, f64),
    out: &mut Vec<(f64, f64, f64, Vec<f64>)>,
    unconverged: &mut Weighted<1>,
) -> f64 {
    if *remaining < 2 {
        unconverged.note_unconverged([parent_err]);
        out.push((a, b, whole.0, whole.1));
        return whole.0;
    }
    *remaining -= 2;
    unconverged.rules += 2;
    let m = 0.5 * (a + b);
    let (left, lv) = gl.integrate_with_values(f, a, m);
    let (right, rv) = gl.integrate_with_values(f, m, b);
    let both = left + right;
    let err = (both - whole.0).abs();
    let converged = err <= tol.max(rel_tol * both.abs());
    if converged || depth == 0 || !(m > a && m < b) || *remaining < 2 {
        if !converged {
            unconverged.note_unconverged([err]);
        }
        out.push((a, m, left, lv));
        out.push((m, b, right, rv));
        return both;
    }
    let half_tol = 0.5 * tol;
    let l = refine_recorded(
        gl,
        f,
        (a, m),
        (left, lv),
        (half_tol, rel_tol),
        depth - 1,
        (&mut *remaining, err),
        out,
        unconverged,
    );
    let r = refine_recorded(
        gl,
        f,
        (m, b),
        (right, rv),
        (half_tol, rel_tol),
        depth - 1,
        (&mut *remaining, err),
        out,
        unconverged,
    );
    l + r
}

/// [`integrate_segments`] for one integrand, returning the accepted panels
/// instead of only their sum. `breaks` must be strictly increasing. The
/// tolerance of segment `i` is `rel_tol · max(|I_i|, floor · Σ|I|/N)`
/// (halved at each level of the bisection); a smaller `floor` holds every
/// segment to `rel_tol` of its own integral down to that fraction of the
/// mean, so that a sub-range made of a few segments is as accurate,
/// relative to itself, as the whole. A cell is also accepted when it meets
/// `rel_tol` relative to its own integral, and each segment's bisection has
/// a budget of [`SEGMENT_BUDGET`] rule applications
/// ([`refine_recorded`], #298); the panels left without meeting their
/// tolerance are counted in [`Panels::unconverged_panels`]. (Before #298, with
/// `floor = 1` this was exactly the error control of [`integrate_segments`];
/// the relative acceptance makes it accept some cells earlier, and widens the
/// bound on the estimated error by the sum of the accepted cells' own
/// integrals times `rel_tol`: see the contract of [`integrate_weighted`].)
pub(crate) fn integrate_panels(
    gl: &GaussLegendre,
    f: &mut impl FnMut(f64) -> f64,
    breaks: &[f64],
    rel_tol: f64,
    floor: f64,
) -> Panels {
    let mut leaves: Vec<(f64, f64, f64, Vec<f64>)> = Vec::new();
    let mut stats = Weighted {
        value: [0.0],
        rules: 0,
        unconverged_cells: 0,
        unconverged_error: [0.0],
        largest_segment_rules: 0,
    };
    if breaks.len() >= 2 {
        let wholes: Vec<(f64, Vec<f64>)> = breaks
            .windows(2)
            .map(|w| gl.integrate_with_values(f, w[0], w[1]))
            .collect();
        stats.rules += wholes.len() as u64;
        let scale: f64 = wholes.iter().map(|x| x.0.abs()).sum();
        let per_segment = floor * scale / wholes.len() as f64;
        for (w, whole) in breaks.windows(2).zip(wholes) {
            let tol = rel_tol * whole.0.abs().max(per_segment);
            let mut budget = SEGMENT_BUDGET;
            refine_recorded(
                gl,
                f,
                (w[0], w[1]),
                whole,
                (tol, rel_tol),
                MAX_DEPTH,
                (&mut budget, 0.0),
                &mut leaves,
                &mut stats,
            );
            stats.largest_segment_rules = stats.largest_segment_rules.max(SEGMENT_BUDGET - budget);
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
    Panels {
        edges,
        values,
        cum,
        unconverged: stats.unconverged_cells,
        unconverged_error: stats.unconverged_error[0],
        inner_unconverged: 0,
    }
}

/// Budget of rule applications for the bisection of one segment of
/// [`integrate_weighted`] and of [`integrate_panels`] (each application is
/// `n` evaluations of the integrand). A segment that has not met its
/// tolerance when the budget is spent is not reported as converged: see
/// [`Weighted`].
///
/// Chosen from measurements (#298) on the full Penn loss function, model
/// tolerance 1e-4: the most any one knot interval used was 78 at the two
/// Si edge-resonance points of `tests/penn_full.rs` (1.2 and 5.8 million
/// before the relative acceptance), and at most 62 on scans of 2000 momenta
/// from 0.005 to 40 a.u. at 11 losses each on the committed Si and Al ELFs.
/// The budget is some 50 times that, so it acts only as a backstop against
/// an integrand that cannot meet the tolerance: on the Si resonance at
/// tolerances of 1e-7 and below, the floating-point noise of the integrand
/// exhausts it, and the evaluation says so (`docs/validation.md`).
pub(crate) const SEGMENT_BUDGET: u64 = 4096;

/// The result of [`integrate_weighted`] with what it cost and whether every
/// cell met its tolerance (also the bookkeeping of [`integrate_panels`],
/// whose `value` is unused).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Weighted<const N: usize> {
    /// The integral.
    pub(crate) value: [f64; N],
    /// Applications of the `n`-point rule, all stages (first estimates,
    /// group checks and bisection). Deterministic for a given integrand.
    pub(crate) rules: u64,
    /// Cells accepted *without* meeting their tolerance: at [`MAX_DEPTH`],
    /// where the midpoint is no longer representable, or when the segment's
    /// [`SEGMENT_BUDGET`] ran out. Zero when the integration converged.
    pub(crate) unconverged_cells: u64,
    /// The sum of `|halves - whole|` over those cells: the error estimate of
    /// the part of `value` that did not converge (zero when it did).
    pub(crate) unconverged_error: [f64; N],
    /// The most rule applications the bisection of any one segment used (at
    /// most [`SEGMENT_BUDGET`]).
    pub(crate) largest_segment_rules: u64,
}

impl<const N: usize> Weighted<N> {
    /// Counts one cell accepted without meeting its tolerance, with its
    /// error estimate.
    fn note_unconverged(&mut self, err: [f64; N]) {
        self.unconverged_cells += 1;
        self.unconverged_error = add(self.unconverged_error, err);
    }
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
/// to single segments, and a single failing segment is refined by bisection
/// ([`refine_bounded`]) starting from its first estimate. Segments on which
/// `w` vanishes at both ends are skipped, as `k` is not evaluated there.
///
/// # The bisection of a segment, and its bound (#298)
///
/// The bisection is that of [`integrate_segments`] (halves compared with the
/// whole, the absolute tolerance `rel_tol · max(|I_s|, Σ|I|/N)` halved at
/// each level) with one more way to accept a cell: its halves agree with
/// the whole to within `rel_tol` of the halves' own integral. A cell
/// accepted either way has an estimated error at most `rel_tol` times the
/// larger of its share of the halved tolerance and its own integral, so the
/// estimated error of the total is at most
/// `rel_tol · (2 Σ|I| + Σ_cells|I_c|)`: the absolute tolerance of a segment,
/// `rel_tol · max(|I_s|, Σ|I|/N)`, sums over the segments to at most
/// `2 · rel_tol · Σ|I|` (the contract of [`integrate_segments`] on its own),
/// and the cells accepted by the relative test add at most
/// `rel_tol · Σ_cells|I_c|`. That is at most `3 · rel_tol` times the integral
/// of `|w k|` (up to the first estimates of `Σ|I|`).
///
/// Without the relative acceptance, a narrow peak that carries much of a
/// segment's integral cannot converge: a cell at depth `d` must meet
/// `2^-d` of the segment's tolerance, so the cells on the peak are asked for
/// an accuracy relative to their own integral of about `rel_tol 2^-d` divided
/// by the peak's share, which falls below the floating-point noise of the
/// integrand long before [`MAX_DEPTH`]. The bisection then fills the peak
/// with cells down to [`MAX_DEPTH`]: on the full Penn loss function of the
/// Si ELF at two `(q, ω)` points this was 1.2 and 5.8 million rule
/// applications for one segment, with an error of the exhausted cells of
/// 1e-11 of the result (`docs/validation.md`, "Electron oracles").
///
/// Each segment's bisection also has a budget of
/// [`SEGMENT_BUDGET`] rule applications. Cells left without meeting
/// their tolerance (budget spent, [`MAX_DEPTH`] or a midpoint that is not
/// representable) keep their two-half estimate and are counted, with their
/// error estimates, in [`Weighted`], so that an exhausted refinement is never
/// reported as converged.
///
/// Every operation is a fixed sequence of floating-point operations, so the
/// result does not depend on the thread that computes it.
pub(crate) fn integrate_weighted<const N: usize>(
    gl: &GaussLegendre,
    weight: &impl Fn(f64) -> f64,
    kernel: &mut impl FnMut(f64) -> [f64; N],
    breaks: &[f64],
    rel_tol: f64,
) -> Weighted<N> {
    let mut out = Weighted {
        value: [0.0; N],
        rules: 0,
        unconverged_cells: 0,
        unconverged_error: [0.0; N],
        largest_segment_rules: 0,
    };
    if breaks.len() < 2 {
        return out;
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
        out.rules += 1;
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
                &mut out,
            ),
        );
        start = end;
    }
    out.value = total;
    out
}

/// Segments `lo..hi` of [`integrate_weighted`] (all active): accept, split or
/// refine. The work and any unconverged cells are added to `stats` (whose
/// `value` is not touched).
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
    stats: &mut Weighted<N>,
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
        let mut budget = SEGMENT_BUDGET;
        let r = refine_bounded(
            gl,
            &mut f,
            (breaks[lo], breaks[lo + 1]),
            whole[lo],
            (tol, rel_tol),
            MAX_DEPTH,
            (&mut budget, [0.0; N]),
            stats,
        );
        stats.largest_segment_rules = stats.largest_segment_rules.max(SEGMENT_BUDGET - budget);
        return r;
    }
    stats.rules += 1;
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
    let left = weighted_group(
        gl,
        weight,
        kernel,
        breaks,
        args,
        per_segment,
        (lo, mid),
        rel_tol,
        stats,
    );
    let right = weighted_group(
        gl,
        weight,
        kernel,
        breaks,
        args,
        per_segment,
        (mid, hi),
        rel_tol,
        stats,
    );
    add(left, right)
}

/// The bisection of one segment of [`integrate_weighted`] (see there for the
/// acceptance rule and its error bound). `[a, b]` has the one-rule estimate
/// `whole`; a cell is accepted when every component of its halves agrees
/// with its whole to within the larger of `tol` (halved at each level) and
/// `rel_tol` times the halves' own integral.
///
/// `remaining` is the segment's budget of rule applications, spent depth
/// first (left half before right half). A cell that cannot be bisected
/// (budget spent, depth 0, or a midpoint that is not representable) keeps
/// its best estimate and is counted in `stats` with its error estimate: the
/// disagreement of its halves with its whole, or, for a cell reached after
/// the budget ran out (only its one-rule estimate `whole` is known), the
/// disagreement `parent_err` of the cell it was cut from (an overestimate,
/// on the safe side).
#[allow(clippy::too_many_arguments)]
fn refine_bounded<const N: usize>(
    gl: &GaussLegendre,
    f: &mut impl FnMut(f64) -> [f64; N],
    (a, b): (f64, f64),
    whole: [f64; N],
    (tol, rel_tol): ([f64; N], f64),
    depth: u32,
    (remaining, parent_err): (&mut u64, [f64; N]),
    stats: &mut Weighted<N>,
) -> [f64; N] {
    if *remaining < 2 {
        stats.note_unconverged(parent_err);
        return whole;
    }
    *remaining -= 2;
    stats.rules += 2;
    let m = 0.5 * (a + b);
    let left = gl.integrate(f, a, m);
    let right = gl.integrate(f, m, b);
    let both = add(left, right);
    let mut err = [0.0; N];
    let mut converged = true;
    for c in 0..N {
        err[c] = (both[c] - whole[c]).abs();
        if err[c] > tol[c].max(rel_tol * both[c].abs()) {
            converged = false;
        }
    }
    if converged {
        return both;
    }
    if depth == 0 || !(m > a && m < b) || *remaining < 2 {
        stats.note_unconverged(err);
        return both;
    }
    let half_tol = tol.map(|t| 0.5 * t);
    let l = refine_bounded(
        gl,
        f,
        (a, m),
        left,
        (half_tol, rel_tol),
        depth - 1,
        (&mut *remaining, err),
        stats,
    );
    let r = refine_bounded(
        gl,
        f,
        (m, b),
        right,
        (half_tol, rel_tol),
        depth - 1,
        (&mut *remaining, err),
        stats,
    );
    add(l, r)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A deterministic pseudo-random number in `[-1, 1]` from the bits of
    /// `x` (a 64-bit finaliser), standing in for floating-point noise.
    fn noise(x: f64) -> f64 {
        let mut n = x.to_bits();
        n ^= n >> 33;
        n = n.wrapping_mul(0xff51_afd7_ed55_8ccd);
        n ^= n >> 33;
        n = n.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        n ^= n >> 33;
        (n >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    }

    #[test]
    fn weighted_bisection_of_a_noisy_narrow_peak_is_bounded() {
        // The case of #298 in miniature: a Lorentzian of width 1e-7 that
        // carries most of the integral of one segment, evaluated with a
        // relative noise of 1e-9 (the cancellation in Re ε near its zero).
        // A cell at depth d is held to 2^-d of the segment tolerance, which
        // the noise does not allow near the peak; the acceptance relative to
        // the cell's own integral lets the bisection stop after a few
        // hundred rules, converged, and within the tolerance.
        let gl = GaussLegendre::new(5);
        let (x0, g) = (0.3 + 1.234e-7, 1e-7);
        let lorentz = |x: f64| g / ((x - x0) * (x - x0) + g * g);
        let exact = |lo: f64, hi: f64| ((hi - x0) / g).atan() - ((lo - x0) / g).atan();
        let breaks = [0.3, 0.31, 0.32, 0.33, 0.34];
        let r = integrate_weighted(
            &gl,
            &|_| 1.0,
            &mut |x: f64| [lorentz(x) * (1.0 + 1e-9 * noise(x))],
            &breaks,
            1e-4,
        );
        let want = exact(0.3, 0.34);
        assert!(
            (r.value[0] / want - 1.0).abs() < 1e-4,
            "{:?} vs {want}",
            r.value
        );
        assert_eq!(r.unconverged_cells, 0, "{r:?}");
        assert!(r.largest_segment_rules < 400, "{r:?}");
        assert!(r.rules < 500, "{r:?}");
    }

    #[test]
    fn weighted_budget_exhaustion_is_reported_not_converged() {
        // A kernel whose noise (1e-3 relative) is far above the tolerance
        // can never converge: every segment spends its budget, the cells
        // left are counted with their error estimates, and the value is
        // still the integral to within the noise.
        let gl = GaussLegendre::new(5);
        let breaks = [0.0, 1.0, 2.0];
        let r = integrate_weighted(
            &gl,
            &|_| 1.0,
            &mut |x: f64| [1.0 + 1e-3 * noise(x)],
            &breaks,
            1e-8,
        );
        assert!(r.unconverged_cells > 0);
        assert!(r.unconverged_error[0] > 0.0);
        assert_eq!(r.largest_segment_rules, SEGMENT_BUDGET);
        // first estimates, one group check, and two budgets
        assert!(r.rules <= 2 + 1 + 2 * SEGMENT_BUDGET, "{r:?}");
        assert!((r.value[0] - 2.0).abs() < 1e-3, "{:?}", r.value);
    }

    #[test]
    fn weighted_depth_exhaustion_at_an_interior_jump_is_reported() {
        // A jump that is not at a break point: the cell holding it never
        // meets its tolerance, so the bisection follows it to MAX_DEPTH
        // (a few rules per level) and counts that last cell, whose error
        // estimate is the jump times its width (2^-40). At 1/3 the jump sits
        // at a third or two thirds of every cell that holds it, between
        // Gauss nodes, where every rule sees it.
        let gl = GaussLegendre::new(5);
        let jump = 1.0 / 3.0;
        let r = integrate_weighted(
            &gl,
            &|_| 1.0,
            &mut |x: f64| [if x < jump { 1.0 } else { 2.0 }],
            &[0.0, 1.0],
            1e-6,
        );
        let want = jump + 2.0 * (1.0 - jump);
        assert!((r.value[0] - want).abs() < 1e-9, "{:?} vs {want}", r.value);
        assert!(r.unconverged_cells >= 1, "{r:?}");
        assert!(r.unconverged_error[0] < 1e-10, "{r:?}");
        assert!(r.largest_segment_rules <= 6 * u64::from(MAX_DEPTH), "{r:?}");
    }

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
        // the same integral to the tolerance (the panels also accept a cell
        // that meets the tolerance relative to itself, #298)
        assert!((p.integral(&gl, 0.0, 10.0) / whole - 1.0).abs() < 1e-8);
        assert_eq!(p.unconverged_panels(), 0);
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
