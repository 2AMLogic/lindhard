//! The loss-density machinery of the inelastic table (parent module docs):
//! the piecewise-linear density and its adaptive refinement, the
//! [`LossModel`] interface of the ELF-extension models, the rows' loss
//! densities and quantiles, and the shared probability-grid refinement.

use super::super::full_penn::DiimfpGrid;
#[cfg(doc)]
use super::super::full_penn::FullPenn;
use super::super::full_penn::FULL_LOWEST_LOSS_FRACTION;
use super::super::mermin::LOWEST_LOSS_FRACTION;
use super::super::model::PennInelastic;
use super::super::penn::SinglePolePenn;
use super::{InelasticTableError, INITIAL_LOSS_PANELS, MAX_NODES};
use crate::electron::data::{CrossSectionTable, ElectronDataError};
use crate::electron::elastic::table::{logistic, logit};
use rayon::prelude::*;

// ---------------------------------------------------------------------------
// Piecewise-linear density in one variable

/// A density `g(x)` that is linear between nodes, with its cumulative
/// integral; `x` is `ln W` or `ln q`.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct LinearDensity {
    pub(super) x: Vec<f64>,
    g: Vec<f64>,
    /// Cumulative integral at each node.
    cum: Vec<f64>,
}

impl LinearDensity {
    fn from_nodes(x: Vec<f64>, g: Vec<f64>) -> Self {
        let mut cum = Vec::with_capacity(x.len());
        let mut acc = 0.0;
        cum.push(0.0);
        for i in 1..x.len() {
            acc += 0.5 * (x[i] - x[i - 1]) * (g[i] + g[i - 1]);
            cum.push(acc);
        }
        Self { x, g, cum }
    }

    pub(super) fn total(&self) -> f64 {
        *self.cum.last().expect("non-empty")
    }

    /// `∫ e^x g dx` of the piecewise-linear model (closed form per panel).
    fn first_moment(&self) -> f64 {
        (1..self.x.len())
            .map(|i| panel_exp_moment(self.x[i - 1], self.x[i], self.g[i - 1], self.g[i]))
            .sum()
    }

    /// The inverse CDF: the `x` with cumulative probability `u`.
    pub(super) fn quantile(&self, u: f64) -> f64 {
        let n = self.x.len();
        if u <= 0.0 {
            return self.x[0];
        }
        if u >= 1.0 {
            return self.x[n - 1];
        }
        let target = u * self.total();
        // First node with cum >= target.
        let j = self.cum.partition_point(|&c| c < target).clamp(1, n - 1);
        let (x0, x1) = (self.x[j - 1], self.x[j]);
        let (g0, g1) = (self.g[j - 1], self.g[j]);
        let h = x1 - x0;
        let rem = (target - self.cum[j - 1]).max(0.0);
        // g0 t + (g1 - g0) t² / (2h) = rem, t in [0, h].
        let disc = (g0 * g0 + 2.0 * (g1 - g0) * rem / h).max(0.0);
        let den = g0 + disc.sqrt();
        let t = if den > 0.0 { 2.0 * rem / den } else { 0.0 };
        x0 + t.clamp(0.0, h)
    }
}

/// `∫_{x0}^{x1} e^x g(x) dx` for `g` linear from `g0` to `g1`.
fn panel_exp_moment(x0: f64, x1: f64, g0: f64, g1: f64) -> f64 {
    // ∫ e^x (g0 + m (x - x0)) dx = e^{x0} [ g0 (E - 1) + m (h E - (E - 1)) ]
    // with h = x1 - x0 and E = e^h.
    let h = x1 - x0;
    if h <= 0.0 {
        return 0.0;
    }
    let m = (g1 - g0) / h;
    let em1 = h.exp_m1();
    x0.exp() * (g0 * em1 + m * (h * (em1 + 1.0) - em1))
}

/// Build the density model of `f` by adaptive panel splitting (module docs):
/// `breaks` are fixed nodes, each segment is first cut into `sub` equal
/// panels. `rough` lists ranges of `x` (ascending) where `f` is known not to
/// be smooth below the width `rough_width` (the unresolved cells of a
/// [`DiimfpGrid`], module docs): a panel within `rough_width` of one of them
/// is not split once it is no wider than `rough_width`. With no ranges the
/// refinement is the plain one. A panel whose midpoint is not a floating-point
/// number strictly inside it is never split (module docs, "Floating-point
/// resolution").
///
/// Returns the density and whether the refinement stopped at [`MAX_NODES`]
/// with a panel still failing the split test.
pub(super) fn adaptive_density(
    f: &(impl Fn(f64) -> f64 + ?Sized),
    breaks: &[f64],
    sub: usize,
    tol: f64,
    first_moment: bool,
    rough: &[(f64, f64)],
    rough_width: f64,
) -> (LinearDensity, bool) {
    let at_resolution = |a: f64, b: f64| {
        b - a <= rough_width
            && rough
                .iter()
                .any(|&(lo, hi)| a <= hi + rough_width && b >= lo - rough_width)
    };
    let mut x: Vec<f64> = Vec::new();
    for w in breaks.windows(2) {
        for i in 0..sub {
            x.push(w[0] + (w[1] - w[0]) * i as f64 / sub as f64);
        }
    }
    x.push(*breaks.last().expect("two breaks"));
    let mut g: Vec<f64> = x.iter().map(|&xi| f(xi)).collect();
    // Cached midpoint values of the panels, once evaluated.
    let mut mids: Vec<Option<f64>> = vec![None; x.len() - 1];
    let span = x[x.len() - 1] - x[0];
    loop {
        for (i, m) in mids.iter_mut().enumerate() {
            if m.is_none() {
                *m = Some(f(0.5 * (x[i] + x[i + 1])));
            }
        }
        let d = LinearDensity::from_nodes(x.clone(), g.clone());
        let i0 = d.total();
        let i1 = if first_moment { d.first_moment() } else { 0.0 };
        let mut split = vec![false; mids.len()];
        let mut any = false;
        let mut capped = false;
        for i in 0..mids.len() {
            let h = x[i + 1] - x[i];
            // The midpoint of a panel one ulp wide rounds onto an end: there
            // is nothing left to split (module docs, "Floating-point
            // resolution").
            let xm = 0.5 * (x[i] + x[i + 1]);
            if !(xm > x[i] && xm < x[i + 1]) {
                continue;
            }
            let gm = mids[i].expect("filled");
            let allowed = h / span * tol;
            let e0 = (2.0 * h / 3.0 * (gm - 0.5 * (g[i] + g[i + 1]))).abs();
            let mut bad = e0 > allowed * i0;
            if first_moment {
                let w = |xx: f64, gg: f64| xx.exp() * gg;
                let lin = 0.5 * (w(x[i], g[i]) + w(x[i + 1], g[i + 1]));
                let e1 = (2.0 * h / 3.0 * (w(xm, gm) - lin)).abs();
                bad |= e1 > allowed * i1;
            }
            if bad && !at_resolution(x[i], x[i + 1]) {
                if x.len() + 1 < MAX_NODES {
                    split[i] = true;
                    any = true;
                } else {
                    capped = true;
                }
            }
        }
        if !any {
            return (d, capped);
        }
        let mut nx = Vec::with_capacity(x.len() + 16);
        let mut ng = Vec::with_capacity(x.len() + 16);
        let mut nm = Vec::with_capacity(mids.len() + 16);
        for i in 0..mids.len() {
            nx.push(x[i]);
            ng.push(g[i]);
            if split[i] {
                nx.push(0.5 * (x[i] + x[i + 1]));
                ng.push(mids[i].expect("filled"));
                nm.push(None);
                nm.push(None);
            } else {
                nm.push(mids[i]);
            }
        }
        nx.push(x[x.len() - 1]);
        ng.push(g[g.len() - 1]);
        x = nx;
        g = ng;
        mids = nm;
    }
}

/// What the table builder needs from an ELF-extension model: the lowest
/// tabulated ELF energy, the lowest loss of its rate integral, the DIIMFP
/// and the inverse IMFP. Implemented for [`SinglePolePenn`] and for every
/// model of [`PennInelastic`]; the table is built the same way for each
/// (module docs).
pub(super) trait LossModel: Sync {
    fn elf_min_ev(&self) -> f64;
    /// The lowest loss (eV) the model's inverse IMFP at `energy_ev` counts,
    /// which opens the sampled loss window of the row (module docs, "Rows
    /// below the ELF table"). By default the lowest tabulated ELF energy,
    /// below which the DIIMFP is zero.
    fn loss_min_ev(&self, _energy_ev: f64) -> f64 {
        self.elf_min_ev()
    }
    /// The kinetic energy `T` (eV, above the model's Fermi level) at and
    /// below which the model's rate is identically zero and above which it
    /// is positive, if it has one: the step of the rate refinement (module
    /// docs, "Rate refinement"; #357). By default the lowest tabulated ELF
    /// energy, the default [`LossModel::loss_min_ev`], below which the DIIMFP
    /// is zero. A model whose rate counts losses below it (Mermin, full
    /// Penn) has no onset.
    fn rate_onset_ev(&self) -> Option<f64> {
        Some(self.elf_min_ev())
    }
    /// A loss strictly inside the window `(loss_min, energy_ev)` at which
    /// the row's initial panels are split into two ranges, so that the
    /// part above it keeps its resolution however far below it the window
    /// opens: the lowest ELF energy for the full Penn model (#361). None by
    /// default.
    fn window_break_ev(&self, _energy_ev: f64) -> Option<f64> {
        None
    }
    fn fermi_energy_ev(&self) -> f64;
    fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError>;
    fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError>;
    /// A tabulation of the DIIMFP for the rows at `energy_ev`, if the model
    /// has one ([`FullPenn::diimfp_grid`]); the rows then read the DIIMFP
    /// from it wherever it covers the loss.
    fn diimfp_grid(&self, _energy_ev: &[f64]) -> Result<Option<DiimfpGrid>, ElectronDataError> {
        Ok(None)
    }
}

impl LossModel for SinglePolePenn {
    fn elf_min_ev(&self) -> f64 {
        self.optical_elf().energy_ev()[0]
    }
    fn fermi_energy_ev(&self) -> f64 {
        SinglePolePenn::fermi_energy_ev(self)
    }
    fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
        self.diimfp_per_m_ev(energy_ev, loss_ev)
    }
    fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
        Ok(self.imfp_and_stopping(energy_ev)?.inverse_imfp_per_m)
    }
}

impl LossModel for PennInelastic {
    fn elf_min_ev(&self) -> f64 {
        self.optical_elf().energy_ev()[0]
    }
    fn loss_min_ev(&self, energy_ev: f64) -> f64 {
        match self {
            // the lower limit of the Mermin rate integral (#342); the fitted
            // ELF has no cutoff at the lowest tabulated energy
            PennInelastic::Mermin(_) => LOWEST_LOSS_FRACTION * energy_ev,
            // the full Penn DIIMFP is not zero below the lowest ELF energy
            // (the electron-hole continuum of each Lindhard gas reaches
            // omega -> 0): the floor of its model-integrated rate (#361)
            PennInelastic::Full(_) => FULL_LOWEST_LOSS_FRACTION * energy_ev,
            // the single-pole DIIMFP is zero below it
            PennInelastic::SinglePole(_) => self.elf_min_ev(),
        }
    }
    fn window_break_ev(&self, energy_ev: f64) -> Option<f64> {
        let b = self.elf_min_ev();
        match self {
            PennInelastic::Full(_) if b > self.loss_min_ev(energy_ev) && b < energy_ev => Some(b),
            _ => None,
        }
    }
    fn rate_onset_ev(&self) -> Option<f64> {
        match self {
            // a positive rate at every T > 0, from losses in (1e-8 T, T]
            // (#342): no onset
            PennInelastic::Mermin(_) => None,
            // a positive rate at every T > 0, from losses in (1e-6 T, T]
            // (#361): no onset
            PennInelastic::Full(_) => None,
            // the rate counts losses above elf_min only (`loss_min_ev`), so
            // it is zero for T <= elf_min
            PennInelastic::SinglePole(_) => Some(self.elf_min_ev()),
        }
    }
    fn fermi_energy_ev(&self) -> f64 {
        PennInelastic::fermi_energy_ev(self)
    }
    fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
        self.diimfp_per_m_ev(energy_ev, loss_ev)
    }
    fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
        Ok(self.imfp_and_stopping(energy_ev)?.inverse_imfp_per_m)
    }
    fn diimfp_grid(&self, energy_ev: &[f64]) -> Result<Option<DiimfpGrid>, ElectronDataError> {
        match self {
            PennInelastic::Full(m) => Ok(Some(m.diimfp_grid(energy_ev)?)),
            _ => Ok(None),
        }
    }
}

/// The DIIMFP of `penn` at `(energy_ev, loss_ev)`, from `grid` where it
/// covers the point; zero where the model fails (only for non-finite input).
fn row_diimfp<M: LossModel>(
    penn: &M,
    grid: Option<&DiimfpGrid>,
    energy_ev: f64,
    loss_ev: f64,
) -> f64 {
    grid.and_then(|g| g.diimfp_per_m_ev(energy_ev, loss_ev))
        .unwrap_or_else(|| penn.diimfp(energy_ev, loss_ev).unwrap_or(0.0))
}

/// The density model of the row at `energy_ev`, on `ln W` from the model's
/// lowest loss ([`LossModel::loss_min_ev`]) to `energy_ev` (module docs),
/// and whether its refinement stopped at [`MAX_NODES`].
fn row_density<M: LossModel>(
    penn: &M,
    grid: Option<&DiimfpGrid>,
    energy_ev: f64,
    tol: f64,
) -> (LinearDensity, bool) {
    let w_lo = penn.loss_min_ev(energy_ev);
    let f = |s: f64| {
        let w = s.exp().min(energy_ev);
        row_diimfp(penn, grid, energy_ev, w) * w
    };
    // The unresolved cells of the grid, on the row's `ln W` axis (module
    // docs, "Unresolved grid cells"); none without a grid.
    let (rough, rough_width) = match grid {
        Some(g) => (g.unresolved_ln_loss_ev(), g.min_cell_width()),
        None => (Vec::new(), 0.0),
    };
    let mut breaks = vec![w_lo.ln(), energy_ev.ln()];
    if let Some(b) = penn.window_break_ev(energy_ev) {
        breaks.insert(1, b.ln());
    }
    adaptive_density(
        &f,
        &breaks,
        INITIAL_LOSS_PANELS,
        tol,
        true,
        &rough,
        rough_width,
    )
}

/// The number of nodes of the loss-density model of the row at `energy_ev`
/// (eV, on the model's own axis) of a table of `model` built with density
/// tolerance `density_tolerance`, reading the DIIMFP from `grid` where it
/// covers the loss as a table build does (pass the grid of
/// [`FullPenn::diimfp_grid`] for the table's energies, or `None`). The
/// count is bounded by [`MAX_NODES`] (a few more: the last refinement pass
/// splits every failing panel); a row that reaches it is refused by the
/// table build ([`InelasticTableError::DensityNodeLimit`]). Each node not
/// covered by the grid costs one direct DIIMFP evaluation.
pub fn loss_density_node_count(
    model: &PennInelastic,
    grid: Option<&DiimfpGrid>,
    energy_ev: f64,
    density_tolerance: f64,
) -> usize {
    row_density(model, grid, energy_ev, density_tolerance)
        .0
        .x
        .len()
}

/// The loss density of one energy.
pub(super) fn loss_density<M: LossModel>(
    penn: &M,
    grid: Option<&DiimfpGrid>,
    energy_ev: f64,
    inverse_mfp: f64,
    tol: f64,
) -> Result<LinearDensity, InelasticTableError> {
    let (d, capped) = row_density(penn, grid, energy_ev, tol);
    if !(d.total() > 0.0 && d.total().is_finite()) {
        return Err(InelasticTableError::EmptyDistribution {
            energy_ev,
            inverse_mfp_per_m: inverse_mfp,
            integral: d.total(),
        });
    }
    if capped {
        return Err(InelasticTableError::DensityNodeLimit {
            energy_ev,
            max_nodes: MAX_NODES,
        });
    }
    Ok(d)
}

/// The loss `W` at cumulative probability `u` of a row, clamped to
/// `[0, energy_ev]`.
pub(super) fn loss_quantile(d: &LinearDensity, energy_ev: f64, u: f64) -> f64 {
    d.quantile(u).exp().min(energy_ev)
}

/// `∫_0^1 W(u) du` by the trapezoid rule of the table's linear
/// interpolation.
pub fn mean_loss_ev(probability: &[f64], quantiles: &[f64]) -> f64 {
    probability
        .windows(2)
        .zip(quantiles.windows(2))
        .map(|(p, q)| 0.5 * (p[1] - p[0]) * (q[0] + q[1]))
        .sum()
}

/// The stopping power `λ⁻¹ ⟨W⟩` (eV/m) of row `i` of an energy-loss table,
/// with `⟨W⟩` from [`mean_loss_ev`]; zero for a zero-rate row. `None` if `i`
/// is out of range.
pub fn stopping_power_ev_per_m(table: &CrossSectionTable, i: usize) -> Option<f64> {
    let inv = *table.inverse_mfp_per_m().get(i)?;
    Some(match table.quantiles(i) {
        Some(q) => inv * mean_loss_ev(table.probability(), q),
        None => 0.0,
    })
}

/// Adaptive refinement of the shared probability grid (module docs).
pub(super) fn refine_probability(
    rows: &[Option<LinearDensity>],
    energy_ev: &[f64],
    mut grid: Vec<f64>,
    mut quantiles: Vec<Vec<f64>>,
    tol: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let means: Vec<f64> = quantiles
        .iter()
        .map(|q| {
            if q.is_empty() {
                0.0
            } else {
                mean_loss_ev(&grid, q)
            }
        })
        .collect();
    // The two segments touching 0 and 1 are left alone.
    let mut todo: Vec<bool> = (0..grid.len() - 1)
        .map(|s| s > 0 && s + 2 < grid.len())
        .collect();
    loop {
        let segs: Vec<usize> = (0..todo.len()).filter(|&s| todo[s]).collect();
        if segs.is_empty() || grid.len() + segs.len() > MAX_NODES {
            return (grid, quantiles);
        }
        let umid: Vec<f64> = segs
            .iter()
            .map(|&s| logistic(0.5 * (logit(grid[s]) + logit(grid[s + 1]))))
            .collect();
        // Exact quantiles at the midpoints, one row per task.
        let exact: Vec<Vec<f64>> = rows
            .par_iter()
            .zip(energy_ev.par_iter())
            .map(|(r, &e)| match r {
                Some(d) => umid.iter().map(|&u| loss_quantile(d, e, u)).collect(),
                None => Vec::new(),
            })
            .collect();
        let mut split = vec![false; segs.len()];
        for (k, &s) in segs.iter().enumerate() {
            let du = grid[s + 1] - grid[s];
            if !(umid[k] > grid[s] && umid[k] < grid[s + 1]) {
                continue; // grid exhausted to rounding
            }
            // The table interpolates linearly in u, at the logistic midpoint
            // as at any other point.
            let t = (umid[k] - grid[s]) / du;
            for (i, q) in quantiles.iter().enumerate() {
                if q.is_empty() {
                    continue;
                }
                let lin = q[s] + t * (q[s + 1] - q[s]);
                if (lin - exact[i][k]).abs() * du > tol * means[i] * du.max(1.0e-4) {
                    split[k] = true;
                    break;
                }
            }
        }
        if !split.iter().any(|&b| b) {
            return (grid, quantiles);
        }
        let mut ngrid = Vec::with_capacity(grid.len() + segs.len());
        let mut nq: Vec<Vec<f64>> = quantiles
            .iter()
            .map(|q| {
                Vec::with_capacity(if q.is_empty() {
                    0
                } else {
                    q.len() + segs.len()
                })
            })
            .collect();
        let mut ntodo = Vec::with_capacity(todo.len() + segs.len());
        let mut k = 0;
        for s in 0..todo.len() {
            ngrid.push(grid[s]);
            for (i, q) in quantiles.iter().enumerate() {
                if !q.is_empty() {
                    nq[i].push(q[s]);
                }
            }
            let mut halves = false;
            if todo[s] {
                if split[k] {
                    ngrid.push(umid[k]);
                    for (i, q) in nq.iter_mut().enumerate() {
                        if !quantiles[i].is_empty() {
                            q.push(exact[i][k]);
                        }
                    }
                    halves = true;
                }
                k += 1;
            }
            ntodo.push(halves);
            if halves {
                ntodo.push(true);
            }
        }
        ngrid.push(grid[grid.len() - 1]);
        for (i, q) in quantiles.iter().enumerate() {
            if !q.is_empty() {
                nq[i].push(q[q.len() - 1]);
            }
        }
        grid = ngrid;
        quantiles = nq;
        todo = ntodo;
    }
}

#[cfg(test)]
mod tests {
    use super::super::build::{loss_rows, InelasticTableOptions};
    use super::super::refine::row_rate;
    use super::super::EnergyAxis;
    use super::*;

    /// A row-like integrand in `x = ln W` (a smooth loss peak) that, inside
    /// the cell `[c, c + 1.5e-4]`, is replaced by a second model differing
    /// by 1e-3 relative and carrying noise at 1e-4: the shape of a full-Penn
    /// row at an unresolved [`DiimfpGrid`] cell, where the grid gives way to
    /// the direct model (module docs, "Unresolved grid cells"). Synthetic,
    /// not a material.
    fn stepped_row(c: f64) -> impl Fn(f64) -> f64 {
        move |x: f64| {
            let smooth = (-(x - 3.0).powi(2)).exp();
            if (c..=c + 1.5e-4).contains(&x) {
                smooth * (1.0 + 1e-3 + 1e-4 * (1e7 * x).sin())
            } else {
                smooth
            }
        }
    }

    #[test]
    fn a_step_at_an_unresolved_cell_does_not_exhaust_the_nodes() {
        // A starting node (x = 3, the 25th of 64 panels on [0, 8]) falls in
        // the cell, at the peak, so the step is sampled and far above the
        // split threshold (about 3 tol of the mean).
        let c = 3.0 - 0.5e-4;
        let f = stepped_row(c);
        let breaks = [0.0, 8.0];
        let cells = [(c, c + 1.5e-4)];
        let (d, capped) =
            adaptive_density(&f, &breaks, INITIAL_LOSS_PANELS, 1e-5, true, &cells, 1e-4);
        assert!(!capped);
        let (smooth, _) = adaptive_density(
            &|x: f64| (-(x - 3.0).powi(2)).exp(),
            &breaks,
            INITIAL_LOSS_PANELS,
            1e-5,
            true,
            &[],
            0.0,
        );
        let (n, n0) = (d.x.len(), smooth.x.len());
        eprintln!("nodes: {n} with the stepped cell declared, {n0} for the smooth row");
        // The cell costs a few bisection levels around its two edges, not
        // the ~MAX_NODES of a bisection to floating-point resolution.
        assert!(n < n0 + 100, "{n} nodes vs {n0} (MAX_NODES = {MAX_NODES})");
        // The density still matches the smooth row's integral to far below
        // the model tolerance: the step adds ~1.5e-7 of it.
        let rel = (d.total() / smooth.total() - 1.0).abs();
        assert!(rel < 1e-5, "{rel:e}");
        // No panel was bisected towards floating-point resolution.
        let narrowest = d.x.windows(2).map(|w| w[1] - w[0]).fold(f64::MAX, f64::min);
        assert!(narrowest > 1e-5, "{narrowest:e}");
    }

    #[test]
    fn without_unresolved_cells_the_refinement_is_unchanged() {
        // A smooth row passes no rough ranges and is refined exactly as
        // before; declaring a cell far from the row's structure changes no
        // node either, since no panel near it ever gets that narrow.
        let f = |x: f64| (-(x - 3.0).powi(2)).exp();
        let breaks = [0.0, 8.0];
        let plain = adaptive_density(&f, &breaks, INITIAL_LOSS_PANELS, 1e-5, true, &[], 0.0);
        let far = adaptive_density(
            &f,
            &breaks,
            INITIAL_LOSS_PANELS,
            1e-5,
            true,
            &[(7.5, 7.5 + 1.5e-4)],
            1e-4,
        );
        assert_eq!(plain, far);
    }

    /// The smooth row of [`stepped_row`] with jumps of a fixed relative size
    /// `rel` at `at` (up, then down at the next, and so on): the shape of the
    /// quadrature noise of the Mermin DIIMFP at high `T` (module docs,
    /// "Floating-point resolution"). Synthetic, not a material.
    fn jumped_row(at: &'static [f64], rel: f64) -> impl Fn(f64) -> f64 {
        move |x: f64| {
            let up = at.iter().filter(|&&a| x >= a).count() % 2;
            (-(x - 3.0).powi(2)).exp() * (1.0 + rel * up as f64)
        }
    }

    #[test]
    fn a_jump_is_bisected_to_floating_point_resolution_and_no_further() {
        // Two pairs of jumps of 2e-5 of the density near its peak: each fails
        // the split test at every width (about 3 tol of the mean). Before
        // #345 the panel around each was bisected until it was one ulp wide
        // and then padded with zero-width panels up to MAX_NODES.
        const AT: &[f64] = &[2.8765432, 2.8812345, 3.1234567, 3.1301234];
        let breaks = [0.0, 8.0];
        let (d, capped) = adaptive_density(
            &jumped_row(AT, 2e-5),
            &breaks,
            INITIAL_LOSS_PANELS,
            1e-5,
            true,
            &[],
            0.0,
        );
        let (smooth, _) = adaptive_density(
            &jumped_row(&[], 0.0),
            &breaks,
            INITIAL_LOSS_PANELS,
            1e-5,
            true,
            &[],
            0.0,
        );
        let (n, n0) = (d.x.len(), smooth.x.len());
        eprintln!("nodes: {n} with four jumps, {n0} for the smooth row");
        assert!(!capped);
        // About 50 bisection levels per jump, not MAX_NODES.
        assert!(
            n < n0 + 4 * 60,
            "{n} nodes vs {n0} (MAX_NODES = {MAX_NODES})"
        );
        // No zero-width panel; the narrowest are one ulp wide, at the jumps.
        let narrowest = d.x.windows(2).map(|w| w[1] - w[0]).fold(f64::MAX, f64::min);
        assert!(narrowest > 0.0, "{narrowest:e}");
        assert!(narrowest <= 3.2 * f64::EPSILON, "{narrowest:e}");
        // The jumps add about 1e-7 of the integral.
        let rel = (d.total() / smooth.total() - 1.0).abs();
        assert!(rel < 1e-5, "{rel:e}");
    }

    /// A loss model whose DIIMFP carries relative noise of 1e-2 at every
    /// loss, which no panel width resolves: synthetic.
    struct NoisyRow;

    impl LossModel for NoisyRow {
        fn elf_min_ev(&self) -> f64 {
            0.5
        }
        fn fermi_energy_ev(&self) -> f64 {
            0.0
        }
        fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
            if !(loss_ev > 0.0 && loss_ev <= energy_ev) {
                return Ok(0.0);
            }
            // a deterministic hash of the loss, in [0, 1)
            let h = loss_ev.to_bits().wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 11;
            Ok(1.0 + 1e-2 * (h as f64 / (1u64 << 53) as f64))
        }
        fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
            Ok(energy_ev - 0.5)
        }
    }

    #[test]
    fn a_row_that_reaches_the_node_bound_is_refused() {
        let opts = InelasticTableOptions::new(vec![10.0]);
        match loss_rows(&NoisyRow, &opts) {
            Err(InelasticTableError::DensityNodeLimit {
                energy_ev,
                max_nodes,
            }) => {
                assert_eq!(energy_ev, 10.0);
                assert_eq!(max_nodes, MAX_NODES);
            }
            Err(e) => panic!("unexpected error: {e}"),
            Ok(_) => panic!("a row at the node bound was stored"),
        }
        let msg = loss_rows(&NoisyRow, &opts).unwrap_err().to_string();
        assert!(msg.contains("within 20000 nodes"), "{msg}");
    }

    /// The committed Si ELF (Yang et al. 2019, `validation/data/optical`,
    /// tabulated from 0.5 eV), or `None` in a packaged crate without it.
    fn si_elf() -> Option<crate::electron::data::OpticalElf> {
        let p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../validation/data/optical/si_elf_yang2019.toml");
        if !p.is_file() {
            eprintln!("{} not found; Si check skipped", p.display());
            return None;
        }
        Some(crate::electron::data::OpticalElf::from_toml_file(p).unwrap())
    }

    #[test]
    fn mermin_rows_below_the_elf_table_open_at_the_rate_lower_limit() {
        // Issues #340 and #342: the Mermin model of Si on the band-bottom
        // axis (E_F = W_v + E_g = 13.46 eV, the band of
        // validation/experiments/se_yield.py) with rows at T = E - E_F inside
        // (0, elf_min]. The Mermin rate there is positive (its fitted ELF has
        // no lower cutoff). Before #340 these rows were refused (a reversed
        // ln W window); from #340 to #342 they were stored empty; now the
        // window opens at the lower limit of the rate integral, 1e-8 T, and
        // they carry the model's rate.
        let Some(elf) = si_elf() else { return };
        let ef = 13.46;
        let model = PennInelastic::try_new(
            crate::electron::inelastic::model::PennAlgorithm::Mermin,
            elf,
        )
        .unwrap()
        .with_fermi_energy_ev(ef)
        .unwrap();
        let elf_min = model.elf_min_ev();
        assert_eq!(elf_min, 0.5);
        let t_low = [0.2085, 0.2530, 0.4216];
        let mut grid: Vec<f64> = t_low.iter().map(|t| ef + t).collect();
        // the T = elf_min row, exactly on the band-bottom axis, and one above
        grid.push(ef + elf_min);
        grid.push(ef + 1.0);
        let opts = InelasticTableOptions::new(grid.clone()).with_axis(EnergyAxis::BandBottom);
        let (inv, _, q) = loss_rows(&model, &opts).expect("rows inside (0, elf_min] build");
        for (k, &e) in grid.iter().enumerate() {
            let t = e - ef;
            let w_min = model.loss_min_ev(t);
            assert_eq!(w_min, LOWEST_LOSS_FRACTION * t);
            assert_eq!(inv[k], model.inverse_imfp(t).unwrap(), "T = {t}");
            assert!(inv[k] > 0.0, "T = {t}");
            let row = &q[k];
            // the window is (1e-8 T, T]: the lowest quantile is its lower
            // edge, far below elf_min, and the highest is T
            assert!((row[0] / w_min - 1.0).abs() < 1e-12, "T = {t}: {}", row[0]);
            assert!(row.iter().all(|&w| w > 0.0 && w <= t), "T = {t}");
            assert!((row[row.len() - 1] - t).abs() <= 1e-12 * t, "T = {t}");
            assert!(row.windows(2).all(|w| w[1] >= w[0]), "T = {t}");
        }
    }

    #[test]
    fn the_si_mermin_row_at_1000_ev_stops_at_floating_point_resolution() {
        // Issue #345. The Mermin row at T = 1000 eV of the committed Si ELF
        // (model Fermi energy 13.46 eV, density tolerance 1e-5) reached
        // MAX_NODES (20 099 nodes, 14 307 of them zero-width panels) on 18
        // jumps of its DIIMFP's quadrature noise; the rows at T <= 100 eV
        // need 2 000 to 3 300 nodes.
        let Some(elf) = si_elf() else { return };
        let model = PennInelastic::try_new(
            crate::electron::inelastic::model::PennAlgorithm::Mermin,
            elf,
        )
        .unwrap()
        .with_fermi_energy_ev(13.46)
        .unwrap();
        let (low, _) = row_density(&model, None, 100.0, 1e-5);
        let (d, capped) = row_density(&model, None, 1000.0, 1e-5);
        let (n, n_low) = (d.x.len(), low.x.len());
        eprintln!("nodes: {n} at T = 1000 eV, {n_low} at T = 100 eV");
        assert!(!capped);
        assert!(
            n < 3 * n_low && n < MAX_NODES / 2,
            "{n} nodes vs {n_low} at 100 eV (MAX_NODES = {MAX_NODES})"
        );
        assert_eq!(
            n,
            loss_density_node_count(&model, None, 1000.0, 1e-5),
            "the public count is the row's"
        );
        assert!(d.x.windows(2).all(|w| w[1] > w[0]), "a zero-width panel");
        // The row's mean loss is the model's, S / λ⁻¹, to the tolerances.
        let p = model.imfp_and_stopping(1000.0).unwrap();
        let mean = d.first_moment() / d.total();
        let want = p.stopping_ev_per_m / p.inverse_imfp_per_m;
        let rel = (mean / want - 1.0).abs();
        eprintln!("mean loss {mean} eV vs model {want} eV: {rel:e}");
        assert!(rel < 1e-5, "{rel:e}");
    }

    /// A loss model whose rate counts losses down to `loss_min < elf_min`,
    /// with a flat DIIMFP on `(loss_min, T]`: synthetic.
    struct WideWindow;

    impl LossModel for WideWindow {
        fn elf_min_ev(&self) -> f64 {
            0.5
        }
        fn loss_min_ev(&self, energy_ev: f64) -> f64 {
            0.01 * energy_ev
        }
        fn rate_onset_ev(&self) -> Option<f64> {
            // a positive rate below elf_min
            None
        }
        fn fermi_energy_ev(&self) -> f64 {
            0.0
        }
        fn diimfp(&self, energy_ev: f64, loss_ev: f64) -> Result<f64, ElectronDataError> {
            Ok(if loss_ev > 0.0 && loss_ev <= energy_ev {
                1.0
            } else {
                0.0
            })
        }
        fn inverse_imfp(&self, energy_ev: f64) -> Result<f64, ElectronDataError> {
            Ok(energy_ev * 0.99)
        }
    }

    #[test]
    fn a_row_opens_at_the_model_lowest_loss() {
        // T = 0.3 lies in (loss_min, elf_min]: the row is built, not empty,
        // and both rows open at loss_min = 0.01 T, not at elf_min.
        let opts = InelasticTableOptions::new(vec![0.3, 2.0]);
        let (inv, prob, q) = loss_rows(&WideWindow, &opts).unwrap();
        assert_eq!(inv, vec![0.3 * 0.99, 2.0 * 0.99]);
        for (row, t) in q.iter().zip([0.3, 2.0]) {
            assert!((row[0] / (0.01 * t) - 1.0).abs() < 1e-12, "{}", row[0]);
            assert!((row[row.len() - 1] - t).abs() <= 1e-12 * t);
            // a flat DIIMFP: W(u) = loss_min + u (T - loss_min)
            for (&u, &w) in prob.iter().zip(row) {
                let want = 0.01 * t + u * 0.99 * t;
                assert!((w - want).abs() <= 1e-4 * t, "u {u}: {w} vs {want}");
            }
        }
    }

    /// A loss model with a positive rate and no DIIMFP anywhere: synthetic.
    struct RateWithoutDensity;

    impl LossModel for RateWithoutDensity {
        fn elf_min_ev(&self) -> f64 {
            0.5
        }
        fn fermi_energy_ev(&self) -> f64 {
            0.0
        }
        fn diimfp(&self, _: f64, _: f64) -> Result<f64, ElectronDataError> {
            Ok(0.0)
        }
        fn inverse_imfp(&self, _: f64) -> Result<f64, ElectronDataError> {
            Ok(1.0e6)
        }
    }

    #[test]
    fn a_positive_rate_without_density_above_the_elf_table_is_still_refused() {
        // Below or at elf_min (the default lowest loss, as for the Penn
        // models) the row is empty; above it, a positive rate with a density
        // that integrates to zero is a defect, refused.
        let opts = InelasticTableOptions::new(vec![0.25, 0.5]);
        let (inv, _, q) = loss_rows(&RateWithoutDensity, &opts).unwrap();
        assert_eq!(inv, vec![0.0, 0.0]);
        assert!(q.iter().all(|r| r.is_empty()));
        let opts = InelasticTableOptions::new(vec![0.25, 2.0]);
        match loss_rows(&RateWithoutDensity, &opts) {
            Err(InelasticTableError::EmptyDistribution { energy_ev, .. }) => {
                assert_eq!(energy_ev, 2.0)
            }
            other => panic!("expected EmptyDistribution, got {other:?}"),
        }
    }

    #[test]
    fn full_penn_loss_window_opens_at_a_fraction_of_the_energy_and_has_no_onset() {
        // #361: the window of a full Penn row is (1e-6 T, T] at every T > 0,
        // below elf_min too; the single pole window opens at elf_min. The
        // initial panels are split at elf_min (the part above it keeps its
        // resolution) when elf_min lies inside the window.
        let elf = crate::electron::inelastic::drude::DrudeLorentz::new(vec![
            crate::electron::inelastic::drude::DrudeLorentzOscillator::plasmon(15.0, 4.0),
        ])
        .unwrap()
        .to_optical_elf("synthetic", 6.0, 2000.0, 100)
        .unwrap();
        use crate::electron::inelastic::model::PennAlgorithm;
        let full = PennInelastic::try_new(PennAlgorithm::Full, elf.clone()).unwrap();
        let single = PennInelastic::try_new(PennAlgorithm::SinglePole, elf).unwrap();
        for t in [0.5, 6.0, 40.0] {
            assert_eq!(full.loss_min_ev(t), FULL_LOWEST_LOSS_FRACTION * t);
            assert_eq!(single.loss_min_ev(t), 6.0);
        }
        assert_eq!(full.window_break_ev(5.0), None);
        assert_eq!(full.window_break_ev(6.0), None);
        assert_eq!(full.window_break_ev(40.0), Some(6.0));
        assert_eq!(single.window_break_ev(40.0), None);
        // a full Penn row at T <= elf_min is not empty, and a single-pole
        // one is
        assert!(row_rate(&full, EnergyAxis::ModelFermiLevel, 5.0).unwrap() > 0.0);
        assert_eq!(
            row_rate(&single, EnergyAxis::ModelFermiLevel, 5.0).unwrap(),
            0.0
        );
        assert_eq!(full.rate_onset_ev(), None);
        // With no onset the refinement tests every cell that overlaps its
        // window, the cell straddling the old onset E_F + elf_min included
        // (pinned for a rate of this kind by
        // `a_rate_without_an_onset_keeps_the_full_window`; a refinement of
        // the real model is not repeated here, its rate costs seconds per
        // point in an unoptimised build).
    }
}
