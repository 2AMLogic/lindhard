//! The rate refinement of the inelastic table's energy grid (parent module
//! docs, "Rate refinement"): [`RateRefinement`], the rate of a row and the
//! extra rows that hold its linear interpolation to a tolerance.

use super::build::InelasticTableOptions;
use super::density::LossModel;
use super::{invalid, EnergyAxis, InelasticTableError};
use rayon::prelude::*;

/// The tolerance of the rate refinement (#339) and of the `rate` gate of
/// `lindhard-cli/examples/inelastic_low_energy.rs` (#291), fixed there before
/// any value was measured: the largest relative error of the interpolated
/// `1/λ` over [`RATE_WINDOW_EV`] above the Fermi level must not exceed 1 %.
/// It is the statistical error of the most precise committed `default` δ
/// entries, rounded up: in `validation/experiments/se_yield_results.json`
/// the Poisson floor of Al `default` is 0.0585 on δ = 6.85 at 600 eV
/// (0.85 %) and 0.0571 on δ = 6.52 at 400 eV (0.88 %). Other configurations
/// are more precise: Al `barrier-off` reaches 0.60 % (600 eV), and at 400 eV
/// the Al `barrier-off`, `cutoff-band-bottom`, `phi-low` and `phi-high`
/// entries are 0.63 to 0.90 %; the Au, Cu and Si `default` entries at 400 eV
/// are 1.4, 1.7 and 2.8 %. A grid error below the tolerance cannot be told
/// from the noise of those tables if δ responds at most in proportion to the
/// rate in that window; that response is not measured here.
pub const RATE_TOLERANCE: f64 = 0.01;

/// The window of the rate refinement and of the `rate` gate (#291, #339):
/// `E - E_F` from 2 to 20 eV. The upper end covers the steep rise of the
/// rate above the Fermi level, which for the Al single pole peaks in error at
/// 15 to 17 eV. The lower end lies below the transport cutoff of every
/// `default` δ run of `validation/experiments/se_yield.py` (1 eV above the
/// vacuum level, at least 5 eV above `E_F`); the `cutoff-band-bottom` runs
/// read the rate down to 1 eV above `E_F`. Below 2 eV the rate starts from
/// zero at `E_F`, and a relative tolerance there would bisect towards that
/// node without end.
pub const RATE_WINDOW_EV: (f64, f64) = (2.0, 20.0);

/// The default bound on the rows the rate refinement may add (a safety
/// bound against a rate that fails the test in many cells at once;
/// not a tuned value): 1000 rows, so cells of 18 meV on average across the
/// 18 eV of [`RATE_WINDOW_EV`], far below the 0.1 eV step of the `rate`
/// gate. A step in the rate splits one cell per level and is stopped
/// earlier, when its cell is exhausted to rounding (module docs, "Rate
/// refinement").
pub const DEFAULT_RATE_REFINEMENT_MAX_ROWS: usize = 1000;

/// Extra rows that resolve the inelastic rate `λ⁻¹(E)` near the Fermi level
/// on the band-bottom axis (module docs, "Rate refinement"; #339).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RateRefinement {
    /// The window of `E - E_F` (eV) whose cells are refined: a cell is
    /// refined when it overlaps it.
    pub window_ev: (f64, f64),
    /// The relative tolerance of the linear interpolation of the rate at a
    /// cell's test points (its quartiles and midpoint).
    pub tolerance: f64,
    /// The most rows the refinement may add; more is an error.
    pub max_rows: usize,
}

impl Default for RateRefinement {
    /// [`RATE_WINDOW_EV`], [`RATE_TOLERANCE`] and
    /// [`DEFAULT_RATE_REFINEMENT_MAX_ROWS`].
    fn default() -> Self {
        Self {
            window_ev: RATE_WINDOW_EV,
            tolerance: RATE_TOLERANCE,
            max_rows: DEFAULT_RATE_REFINEMENT_MAX_ROWS,
        }
    }
}

/// The inverse mean free path of the row at model energy `t` (`T`, module
/// docs, "Energy axis"): zero where the row is stored empty (`T <= 0` on the
/// band-bottom axis, or an empty loss window), else the model's, floored at
/// zero. The rate every row of a table stores.
pub(super) fn row_rate<M: LossModel>(
    penn: &M,
    axis: EnergyAxis,
    t: f64,
) -> Result<f64, InelasticTableError> {
    if axis == EnergyAxis::BandBottom && t <= 0.0 {
        // A band-bottom energy at or below the Fermi level: no loss is
        // allowed. (On the model's own axis the model rejects such an
        // energy.)
        return Ok(0.0);
    }
    if t <= penn.loss_min_ev(t) {
        // The sampled loss window (loss_min, T] is empty (#340): no loss can
        // be drawn, so the row is empty, whatever rate the model reports
        // below its lowest loss (module docs, "Rows below the ELF table").
        return Ok(0.0);
    }
    let inv = penn.inverse_imfp(t)?;
    Ok(if inv > 0.0 { inv } else { 0.0 })
}

/// The points at which the rate refinement tests the cell `a..b` (module
/// docs, "Rate refinement"): its first quartile, midpoint and third
/// quartile, each formed as the midpoint of a half, so that a quartile is
/// bitwise the midpoint of the half a bisection would make.
///
/// The interpolation error of a function with constant curvature peaks at
/// the midpoint; the band-bottom rate's curvature changes across a cell,
/// and the relative error divides by a rate that changes across it, so the
/// peak moves off the midpoint (#339 measured misses of up to 1.43 % in
/// cells whose midpoint missed by less than 1 %). The quartiles are the
/// points of the next level of the bisection itself, a fixed dyadic set
/// rather than a tuned count; they bracket the peak, which
/// `rate_error_peak` then locates.
fn rate_test_points(a: f64, b: f64) -> [f64; 3] {
    let m = 0.5 * (a + b);
    [0.5 * (a + m), m, 0.5 * (m + b)]
}

/// The estimated position of the largest relative interpolation error in
/// the cell `a..b`, from its signed relative errors `e` at the quartiles and
/// midpoint (`rate_test_points`; the error is zero at the cell's ends): the
/// vertex of the parabola through the largest of the five samples and its
/// two neighbours (parabolic interpolation of an extremum, the step of
/// Brent's method; R. P. Brent, *Algorithms for Minimization without
/// Derivatives*, Prentice-Hall 1973, ch. 5). The largest sample's
/// neighbours do not exceed it in magnitude, so the vertex lies between
/// them. `None` where every sample is zero or the three are collinear.
fn rate_error_peak(a: f64, b: f64, e: [f64; 3]) -> Option<f64> {
    let s = [0.0, e[0], e[1], e[2], 0.0];
    // The largest sample in magnitude (the midpoint on a tie).
    let k = [2, 1, 3]
        .into_iter()
        .fold(2, |k, j| if s[j].abs() > s[k].abs() { j } else { k });
    let (l, c, r) = (s[k - 1], s[k], s[k + 1]);
    let curvature = l - 2.0 * c + r;
    if c == 0.0 || curvature == 0.0 || !curvature.is_finite() {
        return None;
    }
    let h = 0.25 * (b - a);
    let xk = a + h * k as f64;
    let v = xk + 0.5 * h * (l - r) / curvature;
    (v > a && v < b && v.is_finite()).then_some(v)
}

/// The miss of the interpolated rate at a test point of the rate refinement.
#[derive(Debug, Clone, Copy)]
struct RateMiss {
    /// `(interpolated - exact) / exact` (infinite where the exact rate is
    /// zero and the interpolated one is not).
    rel: f64,
    /// The energy of the test point, eV.
    at: f64,
    /// The interpolated rate, 1/m.
    interpolated: f64,
    /// The exact rate, 1/m.
    exact: f64,
}

/// The energy grid of a table with the rows of the rate refinement of
/// `options` merged in (module docs, "Rate refinement"), or `None` without
/// one. `options` is checked.
fn refined_energy_grid<M: LossModel>(
    penn: &M,
    options: &InelasticTableOptions,
) -> Result<Option<Vec<f64>>, InelasticTableError> {
    let Some(r) = options.rate_refinement else {
        return Ok(None);
    };
    let ef = penn.fermi_energy_ev();
    let (lo, hi) = (ef + r.window_ev.0, ef + r.window_ev.1);
    let grid = &options.energy_ev;
    let overlaps = |a: f64, b: f64| b > lo && a < hi;
    // The rate's onset, where it is at or above the window start (#357):
    // a cell below it, or straddling it, is not tested (module docs). An
    // onset below the window start changes nothing.
    let onset = penn.rate_onset_ev().map(|w| ef + w).filter(|&o| o >= lo);
    let tested = |a: f64, b: f64| overlaps(a, b) && onset.is_none_or(|o| a >= o);
    // The cells that overlap the window are contiguous: [first, last].
    let Some(first) = grid.windows(2).position(|w| overlaps(w[0], w[1])) else {
        return Ok(Some(grid.clone()));
    };
    let last = grid
        .windows(2)
        .rposition(|w| overlaps(w[0], w[1]))
        .expect("one cell");
    let rate = |e: f64| row_rate(penn, options.axis, e - ef);
    // The nodes of the refined span and their exact rates, in parallel and
    // collected in grid order; the lowest failing energy is the error.
    let mut x: Vec<f64> = grid[first..=last + 1].to_vec();
    let mut y: Vec<f64> = x.par_iter().map(|&e| rate(e)).collect::<Result<_, _>>()?;
    // Cells still to test: those that overlap the window, at or above the
    // onset.
    let mut todo: Vec<bool> = x.windows(2).map(|w| tested(w[0], w[1])).collect();
    let mut added = 0usize;
    loop {
        let cells: Vec<usize> = (0..todo.len()).filter(|&i| todo[i]).collect();
        // A cell still to test is the child of a cell that failed. If it
        // can no longer be bisected (its midpoint rounds to an end), the
        // rate does not meet the tolerance there at any width, as at a
        // step: an error naming the lowest such cell, never a grid that
        // leaves it unresolved.
        if let Some(&i) = cells.iter().find(|&&i| {
            let m = 0.5 * (x[i] + x[i + 1]);
            !(m > x[i] && m < x[i + 1])
        }) {
            return Err(invalid(
                "rate refinement",
                format!(
                    "cannot hold the rate to {} between E - E_F = {} and {} eV: cell {}..{} eV \
                     is exhausted to rounding and its parent still failed at a test point (end \
                     rates {:e} and {:e} 1/m)",
                    r.tolerance,
                    r.window_ev.0,
                    r.window_ev.1,
                    x[i],
                    x[i + 1],
                    y[i],
                    y[i + 1]
                ),
            ));
        }
        if cells.is_empty() {
            break;
        }
        // The test points of each cell, in cell order: its quartiles and
        // midpoint (`rate_test_points`), then the estimated peak of its
        // error (`rate_error_peak`) where those pass.
        let points: Vec<[f64; 3]> = cells
            .iter()
            .map(|&i| rate_test_points(x[i], x[i + 1]))
            .collect();
        let exact: Vec<f64> = points
            .par_iter()
            .flat_map_iter(|p| p.iter().copied())
            .map(rate)
            .collect::<Result<_, _>>()?;
        // The miss of the transport's linear interpolation in cell `i` at
        // `e`, where the exact rate is `m`: signed, relative to `m` (a zero
        // exact rate with a nonzero interpolant is an infinite miss), with
        // the interpolated rate.
        let miss = |i: usize, e: f64, m: f64| {
            let lin = y[i] + (e - x[i]) / (x[i + 1] - x[i]) * (y[i + 1] - y[i]);
            let d = lin - m;
            let rel = if d == 0.0 {
                0.0
            } else if m > 0.0 {
                d / m
            } else {
                f64::INFINITY
            };
            RateMiss {
                rel,
                at: e,
                interpolated: lin,
                exact: m,
            }
        };
        let mut worst: Vec<RateMiss> = Vec::with_capacity(cells.len());
        let mut peaks: Vec<(usize, f64)> = Vec::new();
        for (k, (&i, (p, ex))) in cells
            .iter()
            .zip(points.iter().zip(exact.chunks_exact(3)))
            .enumerate()
        {
            let m: Vec<RateMiss> = p.iter().zip(ex).map(|(&e, &v)| miss(i, e, v)).collect();
            let w = m
                .iter()
                .fold(m[1], |a, &b| if b.rel.abs() > a.rel.abs() { b } else { a });
            if w.rel.abs() <= r.tolerance {
                if let Some(e) = rate_error_peak(x[i], x[i + 1], [m[0].rel, m[1].rel, m[2].rel]) {
                    peaks.push((k, e));
                }
            }
            worst.push(w);
        }
        // The estimated peaks of the cells that passed at their quartiles and
        // midpoint, in parallel and in cell order.
        let peak_exact: Vec<f64> = peaks
            .par_iter()
            .map(|&(_, e)| rate(e))
            .collect::<Result<_, _>>()?;
        for (&(k, e), &v) in peaks.iter().zip(&peak_exact) {
            let p = miss(cells[k], e, v);
            if p.rel.abs() > worst[k].rel.abs() {
                worst[k] = p;
            }
        }
        // Split at the midpoint where the miss at any test point exceeds the
        // tolerance.
        let split: Vec<bool> = worst.iter().map(|w| w.rel.abs() > r.tolerance).collect();
        let n_split = split.iter().filter(|&&s| s).count();
        if n_split == 0 {
            break;
        }
        if added + n_split > r.max_rows {
            let k = split.iter().position(|&s| s).expect("one split");
            let i = cells[k];
            let RateMiss {
                at,
                interpolated: lin,
                exact: m,
                ..
            } = worst[k];
            return Err(invalid(
                "rate refinement",
                format!(
                    "more than {} rows needed to hold the rate to {} between E - E_F = {} and \
                     {} eV (cell {}..{} eV still fails at {} eV: interpolated {:e}, exact {:e} \
                     1/m)",
                    r.max_rows,
                    r.tolerance,
                    r.window_ev.0,
                    r.window_ev.1,
                    x[i],
                    x[i + 1],
                    at,
                    lin,
                    m
                ),
            ));
        }
        added += n_split;
        let mut nx = Vec::with_capacity(x.len() + n_split);
        let mut ny = Vec::with_capacity(x.len() + n_split);
        let mut ntodo = Vec::with_capacity(todo.len() + n_split);
        let mut k = 0;
        for i in 0..todo.len() {
            nx.push(x[i]);
            ny.push(y[i]);
            if k < cells.len() && cells[k] == i {
                if split[k] {
                    // The midpoint, a test point: its exact rate is known.
                    let m = points[k][1];
                    nx.push(m);
                    ny.push(exact[3 * k + 1]);
                    ntodo.push(tested(x[i], m));
                    ntodo.push(tested(m, x[i + 1]));
                } else {
                    ntodo.push(false);
                }
                k += 1;
            } else {
                ntodo.push(false);
            }
        }
        nx.push(x[x.len() - 1]);
        ny.push(y[y.len() - 1]);
        x = nx;
        y = ny;
        todo = ntodo;
    }
    let mut out = Vec::with_capacity(grid.len() + added);
    out.extend_from_slice(&grid[..first]);
    out.extend_from_slice(&x);
    out.extend_from_slice(&grid[last + 2..]);
    Ok(Some(out))
}

/// `options` with the rows of its rate refinement merged into the energy
/// grid (and the refinement then spent), or `options` itself without one.
pub(super) fn with_refined_grid<'a, M: LossModel>(
    penn: &M,
    options: &'a InelasticTableOptions,
) -> Result<(std::borrow::Cow<'a, InelasticTableOptions>, usize), InelasticTableError> {
    Ok(match refined_energy_grid(penn, options)? {
        None => (std::borrow::Cow::Borrowed(options), 0),
        Some(energy_ev) => {
            let added = energy_ev.len() - options.energy_ev.len();
            let mut o = options.clone();
            o.energy_ev = energy_ev;
            (std::borrow::Cow::Owned(o), added)
        }
    })
}

/// The provenance text of the rate refinement: empty without one (so tables
/// built without it keep their bytes).
pub(super) fn refinement_provenance(options: &InelasticTableOptions, added: usize) -> String {
    match &options.rate_refinement {
        None => String::new(),
        Some(r) => format!(
            "; rate refinement (#339): {added} rows added where the linear interpolation of \
             the rate missed it by more than {} at a cell's quartiles, midpoint or estimated \
             error peak, E - E_F = {}..{} eV",
            r.tolerance, r.window_ev.0, r.window_ev.1
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electron::data::ElectronDataError;
    use crate::electron::inelastic::mermin::LOWEST_LOSS_FRACTION;
    use crate::electron::inelastic::model::PennInelastic;

    /// A loss model with Fermi energy 10 eV whose rate is `rate(T)`, and a
    /// flat DIIMFP: synthetic.
    struct SteepRate(fn(f64) -> f64);

    impl LossModel for SteepRate {
        fn elf_min_ev(&self) -> f64 {
            0.05
        }
        fn fermi_energy_ev(&self) -> f64 {
            10.0
        }
        fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
            Ok(if loss_ev > 0.0 && loss_ev <= energy_ev {
                1.0
            } else {
                0.0
            })
        }
        fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
            Ok((self.0)(energy_ev))
        }
    }

    fn band_options(grid: Vec<f64>) -> InelasticTableOptions {
        InelasticTableOptions::new(grid)
            .with_axis(EnergyAxis::BandBottom)
            .with_rate_refinement(RateRefinement::default())
    }

    /// A rate that climbs steeply above `E_F`, like the band-bottom rates of
    /// #291 (synthetic).
    fn steep(t: f64) -> f64 {
        1e6 * t.powi(3) / (1.0 + (t / 15.0).powi(4))
    }

    #[test]
    fn rate_refinement_adds_rows_only_in_the_window_and_meets_the_test() {
        let model = SteepRate(steep);
        let grid = vec![5.0, 9.0, 14.0, 22.0, 35.0, 60.0, 100.0];
        let o = band_options(grid.clone());
        let refined = refined_energy_grid(&model, &o).unwrap().unwrap();
        assert!(refined.len() > grid.len());
        assert!(refined.windows(2).all(|w| w[1] > w[0]));
        // The given grid is a subsequence of the refined one.
        let mut it = refined.iter();
        assert!(grid.iter().all(|g| it.any(|r| r == g)));
        let (lo, hi) = (10.0 + RATE_WINDOW_EV.0, 10.0 + RATE_WINDOW_EV.1);
        // Rows are added only inside cells of the grid that overlap the
        // window (here 9..35 eV), never in a cell outside it.
        for e in refined.iter().filter(|e| !grid.contains(e)) {
            let k = grid.partition_point(|g| g < e);
            let (a, b) = (grid[k - 1], grid[k]);
            assert!(b > lo && a < hi, "{e} in the cell {a}..{b}");
        }
        // Every final cell overlapping the window passes the test at its
        // quartiles and midpoint.
        let rate = |e: f64| row_rate(&model, EnergyAxis::BandBottom, e - 10.0).unwrap();
        for w in refined.windows(2) {
            if w[1] > lo && w[0] < hi {
                for e in rate_test_points(w[0], w[1]) {
                    let m = rate(e);
                    let lin = rate(w[0]) + (e - w[0]) / (w[1] - w[0]) * (rate(w[1]) - rate(w[0]));
                    assert!((lin - m).abs() <= RATE_TOLERANCE * m, "{w:?} at {e}");
                }
            }
        }
        // Without the option the grid is the given one.
        let mut off = o.clone();
        off.rate_refinement = None;
        assert_eq!(refined_energy_grid(&model, &off).unwrap(), None);
    }

    #[test]
    fn rate_refinement_is_the_same_on_any_thread_count() {
        let model = SteepRate(steep);
        let o = band_options(vec![5.0, 9.0, 14.0, 22.0, 35.0, 60.0, 100.0]);
        let on = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| refined_energy_grid(&model, &o).unwrap())
        };
        let one = on(1);
        assert_eq!(on(3), one);
        assert_eq!(on(8), one);
    }

    #[test]
    fn rate_refinement_of_a_step_stops_at_the_row_bound() {
        // A step in the rate inside the window cannot be resolved by
        // bisection: an error naming the cell, not a runaway.
        let model = SteepRate(|t| if t < 7.3 { 1e6 } else { 2e6 });
        let mut o = band_options(vec![5.0, 9.0, 14.0, 22.0, 35.0, 60.0, 100.0]);
        o.rate_refinement.as_mut().unwrap().max_rows = 30;
        let err = refined_energy_grid(&model, &o).unwrap_err().to_string();
        assert!(err.contains("more than 30 rows"), "{err}");
        assert!(err.contains("cell 17."), "{err}");
    }

    #[test]
    fn rate_refinement_of_a_step_is_an_error_at_the_default_bound() {
        // A step splits one cell per level, so at the default bound the
        // cell around it reaches rounding (about 51 levels) long before the
        // row bound: still an error naming the cell, never a grid with the
        // step unresolved.
        let model = SteepRate(|t| if t < 7.3 { 1e6 } else { 2e6 });
        let o = band_options(vec![5.0, 9.0, 14.0, 22.0, 35.0, 60.0, 100.0]);
        assert_eq!(
            o.rate_refinement.unwrap().max_rows,
            DEFAULT_RATE_REFINEMENT_MAX_ROWS
        );
        let err = refined_energy_grid(&model, &o).unwrap_err().to_string();
        assert!(err.contains("exhausted to rounding"), "{err}");
        assert!(err.contains("cell 17."), "{err}");
    }

    #[test]
    fn rate_test_points_are_the_quartiles_and_midpoint() {
        let p = rate_test_points(2.0, 6.0);
        assert_eq!(p, [3.0, 4.0, 5.0]);
        // A quartile is bitwise the midpoint of the half a bisection makes.
        let (a, b) = (0.1, 0.7);
        let p = rate_test_points(a, b);
        assert_eq!(p[0], rate_test_points(a, p[1])[1]);
        assert_eq!(p[2], rate_test_points(p[1], b)[1]);
    }

    #[test]
    fn rate_error_peak_is_the_vertex_of_the_parabola_at_the_largest_sample() {
        // Samples of 1.69 - (x - 1.3)^2 at the quartiles and midpoint of
        // 0..4 (zero at x = 0): the largest is at x = 1, and the parabola
        // through it and its neighbours x = 0 and 2 has its vertex at 1.3.
        let g = |x: f64| 1.69 - (x - 1.3) * (x - 1.3);
        let v = rate_error_peak(0.0, 4.0, [g(1.0), g(2.0), g(3.0)]).unwrap();
        assert!((v - 1.3).abs() < 1e-12, "{v}");
        // The same with the sign reversed (an interpolant below the rate).
        let v = rate_error_peak(0.0, 4.0, [-g(1.0), -g(2.0), -g(3.0)]).unwrap();
        assert!((v - 1.3).abs() < 1e-12, "{v}");
        // A symmetric profile peaks at the midpoint; no samples, no peak.
        let v = rate_error_peak(0.0, 4.0, [0.75, 1.0, 0.75]).unwrap();
        assert!((v - 2.0).abs() < 1e-12, "{v}");
        assert_eq!(rate_error_peak(0.0, 4.0, [0.0; 3]), None);
    }

    #[test]
    fn a_window_above_the_grid_adds_no_rows() {
        let model = SteepRate(steep);
        let grid = vec![1.0, 5.0, 11.0];
        let o = band_options(grid.clone());
        assert_eq!(refined_energy_grid(&model, &o).unwrap(), Some(grid));
    }

    /// A loss model with Fermi energy 10 eV, lowest tabulated ELF energy
    /// `elf_min`, rate `rate(T)` and a flat DIIMFP, whose rate counts losses
    /// down to `loss_min_fraction · T` where given (Mermin-style, a positive
    /// rate below `elf_min`) and otherwise from `elf_min` (Penn-style), and
    /// whose reported rate onset is `onset`: synthetic.
    struct OnsetRate {
        elf_min: f64,
        loss_min_fraction: Option<f64>,
        onset: Option<f64>,
        rate: fn(f64) -> f64,
    }

    impl OnsetRate {
        /// Penn-style: zero rate at and below `elf_min`, reported as the
        /// onset.
        fn penn(elf_min: f64, rate: fn(f64) -> f64) -> Self {
            Self {
                elf_min,
                loss_min_fraction: None,
                onset: Some(elf_min),
                rate,
            }
        }
    }

    impl LossModel for OnsetRate {
        fn elf_min_ev(&self) -> f64 {
            self.elf_min
        }
        fn loss_min_ev(&self, energy_ev: f64) -> f64 {
            self.loss_min_fraction
                .map_or(self.elf_min, |f| f * energy_ev)
        }
        fn rate_onset_ev(&self) -> Option<f64> {
            self.onset
        }
        fn fermi_energy_ev(&self) -> f64 {
            10.0
        }
        fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
            Ok(if loss_ev > 0.0 && loss_ev <= energy_ev {
                1.0
            } else {
                0.0
            })
        }
        fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
            Ok((self.rate)(energy_ev))
        }
    }

    /// A rate that is zero up to `T = 6` eV and climbs steeply above it,
    /// linearly at the onset (synthetic).
    fn onset_steep(t: f64) -> f64 {
        let s = (t - 6.0).max(0.0);
        1e6 * s * (1.0 + s * s) / (1.0 + (s / 15.0).powi(4))
    }

    /// Every final cell of `refined` that overlaps the window and starts at
    /// or above `from` passes the test at its quartiles and midpoint.
    fn assert_cells_pass<M: LossModel>(model: &M, refined: &[f64], from: f64) {
        let ef = model.fermi_energy_ev();
        let (lo, hi) = (ef + RATE_WINDOW_EV.0, ef + RATE_WINDOW_EV.1);
        let rate = |e: f64| row_rate(model, EnergyAxis::BandBottom, e - ef).unwrap();
        for w in refined.windows(2) {
            if w[1] > lo && w[0] < hi && w[0] >= from {
                for e in rate_test_points(w[0], w[1]) {
                    let m = rate(e);
                    let lin = rate(w[0]) + (e - w[0]) / (w[1] - w[0]) * (rate(w[1]) - rate(w[0]));
                    assert!((lin - m).abs() <= RATE_TOLERANCE * m, "{w:?} at {e}");
                }
            }
        }
    }

    /// The given grid is a subsequence of the refined one.
    fn assert_subsequence(grid: &[f64], refined: &[f64]) {
        assert!(refined.windows(2).all(|w| w[1] > w[0]));
        let mut it = refined.iter();
        assert!(grid.iter().all(|g| it.any(|r| r == g)));
    }

    #[test]
    fn rate_refinement_starts_at_an_onset_above_the_window() {
        // Issue #357: an ELF tabulated from 6 eV, so the Penn-style rate is
        // zero up to E = E_F + 6 = 16 eV, inside the window 12..30 eV. The
        // cell 14..17 straddles the onset and is left as given; the cells
        // from 17 eV are refined to the tolerance.
        let model = OnsetRate::penn(6.0, onset_steep);
        let onset = 16.0;
        for grid in [
            vec![5.0, 9.0, 14.0, 17.0, 22.0, 35.0, 60.0, 100.0],
            // a node on the onset: the cell from it is tested
            vec![5.0, 9.0, 14.0, 16.0, 22.0, 35.0, 60.0, 100.0],
        ] {
            let o = band_options(grid.clone());
            let refined = refined_energy_grid(&model, &o).unwrap().unwrap();
            assert!(refined.len() > grid.len(), "{grid:?}");
            assert_subsequence(&grid, &refined);
            // Rows are added only in cells of the grid at or above the onset.
            for e in refined.iter().filter(|e| !grid.contains(e)) {
                let k = grid.partition_point(|g| g < e);
                assert!(
                    grid[k - 1] >= onset,
                    "{e} in the cell {}..{}",
                    grid[k - 1],
                    grid[k]
                );
            }
            assert_cells_pass(&model, &refined, onset);
        }
        // Without the onset the straddling cell is bisected towards the step
        // until it is exhausted to rounding (the error before #357).
        let unhooked = OnsetRate {
            onset: None,
            ..OnsetRate::penn(6.0, onset_steep)
        };
        let o = band_options(vec![5.0, 9.0, 14.0, 17.0, 22.0, 35.0, 60.0, 100.0]);
        let err = refined_energy_grid(&unhooked, &o).unwrap_err().to_string();
        assert!(err.contains("exhausted to rounding"), "{err}");
        assert!(err.contains("cell 15.") || err.contains("cell 16"), "{err}");
    }

    #[test]
    fn an_onset_below_the_window_leaves_the_refinement_unchanged() {
        // elf_min = 0.5 eV (the Si ELF; every committed ELF starts at or
        // below it), so the onset E_F + 0.5 lies below the window start
        // E_F + 2: the grid is the one refined without an onset, including
        // rows in the cell 10.3..14 that straddles both.
        let grid = vec![5.0, 9.0, 10.3, 14.0, 22.0, 35.0, 60.0, 100.0];
        let o = band_options(grid.clone());
        let with = OnsetRate::penn(0.5, steep);
        let without = OnsetRate {
            onset: None,
            ..OnsetRate::penn(0.5, steep)
        };
        let refined = refined_energy_grid(&with, &o).unwrap().unwrap();
        assert_eq!(refined, refined_energy_grid(&without, &o).unwrap().unwrap());
        assert!(refined.iter().any(|&e| e > 10.3 && e < 14.0), "{refined:?}");
        // An onset exactly at the window start is the lowest that applies.
        let at_start = OnsetRate::penn(RATE_WINDOW_EV.0, onset_steep);
        let refined = refined_energy_grid(&at_start, &o).unwrap().unwrap();
        assert!(
            !refined.iter().any(|&e| e > 10.3 && e < 14.0),
            "{refined:?}"
        );
        assert_subsequence(&grid, &refined);
    }

    #[test]
    fn a_rate_without_an_onset_keeps_the_full_window() {
        // Mermin-style: elf_min = 6 eV, but the rate counts losses down to
        // 1e-8 T and is positive below elf_min, so it has no onset and the
        // cells below E_F + 6 are refined as before.
        let model = OnsetRate {
            elf_min: 6.0,
            loss_min_fraction: Some(LOWEST_LOSS_FRACTION),
            onset: None,
            rate: steep,
        };
        let grid = vec![5.0, 9.0, 14.0, 17.0, 22.0, 35.0, 60.0, 100.0];
        let refined = refined_energy_grid(&model, &band_options(grid.clone()))
            .unwrap()
            .unwrap();
        assert_subsequence(&grid, &refined);
        assert!(refined.iter().any(|&e| e > 12.0 && e < 16.0), "{refined:?}");
        assert_cells_pass(&model, &refined, f64::NEG_INFINITY);
    }

    #[test]
    fn rate_refinement_from_an_onset_is_the_same_on_any_thread_count() {
        let model = OnsetRate::penn(6.0, onset_steep);
        let o = band_options(vec![5.0, 9.0, 14.0, 17.0, 22.0, 35.0, 60.0, 100.0]);
        let on = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| refined_energy_grid(&model, &o).unwrap())
        };
        let one = on(1);
        assert_eq!(on(3), one);
        assert_eq!(on(8), one);
    }

    #[test]
    fn the_single_pole_model_has_a_rate_onset_and_the_full_and_mermin_models_none() {
        // A synthetic ELF (one Drude-Lorentz plasmon) tabulated from 6 eV.
        let elf = crate::electron::inelastic::drude::DrudeLorentz::new(vec![
            crate::electron::inelastic::drude::DrudeLorentzOscillator::plasmon(15.0, 4.0),
        ])
        .unwrap()
        .to_optical_elf("synthetic", 6.0, 2000.0, 200)
        .unwrap();
        use crate::electron::inelastic::model::PennAlgorithm;
        for (algorithm, onset) in [
            (PennAlgorithm::SinglePole, Some(6.0)),
            // #361: the full Penn DIIMFP is positive below elf_min, so is
            // its rate
            (PennAlgorithm::Full, None),
            (PennAlgorithm::Mermin, None),
        ] {
            let model = PennInelastic::try_new(algorithm, elf.clone()).unwrap();
            assert_eq!(model.rate_onset_ev(), onset, "{algorithm:?}");
            // single pole: the rate is zero at the onset, as the hook says;
            // full Penn and Mermin: positive below elf_min.
            let below = row_rate(&model, EnergyAxis::ModelFermiLevel, 5.0).unwrap();
            assert_eq!(below > 0.0, onset.is_none(), "{algorithm:?}: {below}");
        }
        // The single-pole Penn model on the band-bottom axis builds a refined
        // grid with the default refinement, rows only from the onset.
        let model = PennInelastic::new(PennAlgorithm::SinglePole, elf)
            .with_fermi_energy_ev(10.0)
            .unwrap();
        let grid = vec![5.0, 9.0, 14.0, 17.0, 22.0, 35.0, 60.0, 100.0];
        let refined = refined_energy_grid(&model, &band_options(grid.clone()))
            .unwrap()
            .unwrap();
        assert!(refined.len() > grid.len(), "{refined:?}");
        assert_subsequence(&grid, &refined);
        for e in refined.iter().filter(|e| !grid.contains(e)) {
            assert!(*e > 17.0, "{e}");
        }
        assert_cells_pass(&model, &refined, 16.0);
    }
}
