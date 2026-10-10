//! Backscatter coefficient η(E) per target: lindhard against measured data.
//!
//! Input is `validation/experiments/backscatter_results.json` (format
//! `lindhard-backscatter-results/1`), written by
//! `validation/experiments/backscatter.py`. Per target it holds the measured
//! groups (per comparison energy: median, min and max over the stored sets,
//! which are transcribed from the compilation and measurements cited in
//! `docs/data-provenance.md`) and the lindhard runs (η and its binomial σ).
//! The figure shows, per target, the measured min-max band, the measured
//! median, the pass band `median ± tolerance` over the energies where the
//! pass rule applies (`E ≥ pass_min_kev`), and lindhard's η with ±σ bars.
//! The pass rule itself is pinned in `backscatter.py` and documented in
//! `docs/validation.md`; this module only draws what the results file says.

use std::collections::BTreeMap;

use anyhow::{bail, ensure, Context, Result};
use rizzma::artist::Line2D;
use rizzma::core::Rgba;
use rizzma::{Figure, LegendLocation};
use serde::Deserialize;

/// The results-file format this module reads.
pub const FORMAT: &str = "lindhard-backscatter-results/1";

/// The parts of the results file the figure uses; other fields are ignored.
#[derive(Debug, Deserialize)]
struct RawResults {
    format: String,
    lindhard: String,
    histories: u64,
    tolerance: f64,
    pass_min_kev: f64,
    targets: Vec<RawTarget>,
}

#[derive(Debug, Deserialize)]
struct RawTarget {
    target: String,
    groups: Vec<RawGroup>,
    runs: BTreeMap<String, RawRun>,
}

#[derive(Debug, Deserialize)]
struct RawGroup {
    energy_kev: f64,
    median: f64,
    min: f64,
    max: f64,
}

#[derive(Debug, Deserialize)]
struct RawRun {
    eta: f64,
    eta_se: f64,
}

/// Everything the figure draws, in drawing order.
#[derive(Debug, Clone, PartialEq)]
pub struct Results {
    /// The engine version string recorded with the runs.
    pub lindhard: String,
    /// Primaries per run.
    pub histories: u64,
    /// Pass tolerance on |η − measured median|, absolute.
    pub tolerance: f64,
    /// Lowest energy (keV) at which the pass rule applies.
    pub pass_min_kev: f64,
    /// One entry per target, in file order.
    pub targets: Vec<TargetSeries>,
}

/// The series for one target, each sorted by ascending energy.
#[derive(Debug, Clone, PartialEq)]
pub struct TargetSeries {
    /// Target label (element symbol).
    pub target: String,
    /// Measured groups: energy (keV), median, min, max.
    pub measured: Vec<MeasuredPoint>,
    /// lindhard runs: energy (keV), η, σ.
    pub lindhard: Vec<RunPoint>,
}

/// Measured η at one comparison energy, over the stored sets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeasuredPoint {
    /// Comparison energy, keV.
    pub energy_kev: f64,
    /// Median over the sets.
    pub median: f64,
    /// Smallest value over the sets.
    pub min: f64,
    /// Largest value over the sets.
    pub max: f64,
}

/// lindhard's η at one energy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RunPoint {
    /// Primary energy, keV.
    pub energy_kev: f64,
    /// Backscatter coefficient.
    pub eta: f64,
    /// Its standard error (binomial).
    pub eta_se: f64,
}

/// Parse a results file into the series the figure draws.
///
/// Runs are keyed by energy strings in the file (`"1"`, `"10"`, ...); they are
/// sorted numerically here, not lexically. Measured groups are sorted the same
/// way.
///
/// # Errors
///
/// Fails on malformed JSON, an unexpected `format`, a run key that is not a
/// number, a non-finite value, or a target with no groups or no runs.
pub fn parse(json: &str) -> Result<Results> {
    let raw: RawResults = serde_json::from_str(json).context("parsing backscatter results")?;
    if raw.format != FORMAT {
        bail!("expected format {FORMAT:?}, found {:?}", raw.format);
    }
    ensure!(
        raw.tolerance.is_finite() && raw.pass_min_kev.is_finite(),
        "non-finite tolerance or pass_min_kev"
    );
    let mut targets = Vec::with_capacity(raw.targets.len());
    for t in raw.targets {
        ensure!(
            !t.groups.is_empty(),
            "target {}: no measured groups",
            t.target
        );
        ensure!(!t.runs.is_empty(), "target {}: no runs", t.target);
        let mut measured: Vec<MeasuredPoint> = t
            .groups
            .iter()
            .map(|g| MeasuredPoint {
                energy_kev: g.energy_kev,
                median: g.median,
                min: g.min,
                max: g.max,
            })
            .collect();
        let mut lindhard = Vec::with_capacity(t.runs.len());
        for (key, run) in &t.runs {
            let energy_kev: f64 = key.parse().with_context(|| {
                format!("target {}: run key {key:?} is not an energy", t.target)
            })?;
            lindhard.push(RunPoint {
                energy_kev,
                eta: run.eta,
                eta_se: run.eta_se,
            });
        }
        let finite = measured.iter().all(|m| {
            [m.energy_kev, m.median, m.min, m.max]
                .iter()
                .all(|v| v.is_finite())
        }) && lindhard.iter().all(|r| {
            [r.energy_kev, r.eta, r.eta_se]
                .iter()
                .all(|v| v.is_finite())
        });
        ensure!(finite, "target {}: non-finite value", t.target);
        measured.sort_by(|a, b| a.energy_kev.total_cmp(&b.energy_kev));
        lindhard.sort_by(|a, b| a.energy_kev.total_cmp(&b.energy_kev));
        targets.push(TargetSeries {
            target: t.target,
            measured,
            lindhard,
        });
    }
    ensure!(!targets.is_empty(), "no targets");
    Ok(Results {
        lindhard: raw.lindhard,
        histories: raw.histories,
        tolerance: raw.tolerance,
        pass_min_kev: raw.pass_min_kev,
        targets,
    })
}

// Fixed styling. Colors are from the Okabe-Ito colour-blind-safe palette.
const BAND: Rgba = Rgba::new(0.6, 0.6, 0.6, 0.35);
const MEDIAN: Rgba = Rgba::new(0.0, 0.0, 0.0, 1.0);
const PASS: Rgba = Rgba::new(0.0, 0.447, 0.698, 1.0); // Okabe-Ito blue
const LINDHARD: Rgba = Rgba::new(0.835, 0.369, 0.0, 1.0); // Okabe-Ito vermillion
const FIG_WIDTH_IN: f64 = 12.0;
const FIG_HEIGHT_IN: f64 = 7.5;
const COLUMNS: usize = 3;
/// Upper y limit is the largest drawn value rounded up to this step.
const Y_STEP: f64 = 0.05;

/// Render the figure as an SVG document.
///
/// One panel per target in a grid of three columns, in file order; a spare
/// panel, when the grid has one, holds the legend. The output depends only on
/// `results`: fixed canvas size, colours and font, no timestamps.
#[must_use]
pub fn render_svg(results: &Results) -> String {
    // rizzma writes the document on one line; one element per line keeps a
    // regenerated figure's `git diff` readable. Glyph paths contain no `><`.
    let mut svg = figure(results).to_svg().replace("><", ">\n<");
    svg.push('\n');
    svg
}

fn figure(results: &Results) -> Figure {
    let n = results.targets.len();
    // Leave room for a legend panel when the targets fill the grid exactly.
    let rows = n / COLUMNS + 1;
    let mut fig = Figure::new(FIG_WIDTH_IN, FIG_HEIGHT_IN);
    fig.suptitle(format!(
        "Backscatter coefficient \u{3b7} vs primary energy E: {}, {} primaries per run",
        results.lindhard, results.histories
    ));
    let x_max = results
        .targets
        .iter()
        .flat_map(|t| {
            t.measured
                .iter()
                .map(|m| m.energy_kev)
                .chain(t.lindhard.iter().map(|r| r.energy_kev))
        })
        .fold(0.0_f64, f64::max);
    let x_lim = (0.0, x_max * 1.05);
    for (i, t) in results.targets.iter().enumerate() {
        let ax = fig.add_subplot(rows, COLUMNS, i + 1);
        draw_target(ax, t, results.tolerance, results.pass_min_kev, x_lim);
    }
    let legend = fig.add_subplot(rows, COLUMNS, rows * COLUMNS);
    legend.set_axis_off();
    legend.legend_at(
        vec![
            (BAND, "measured: min-max over the sets".to_owned()),
            (MEDIAN, "measured: median over the sets".to_owned()),
            (
                PASS,
                format!(
                    "pass band: median \u{b1} {} (E \u{2265} {} keV)",
                    results.tolerance, results.pass_min_kev
                ),
            ),
            (LINDHARD, "lindhard \u{3b7} \u{b1} \u{3c3}".to_owned()),
        ],
        LegendLocation::UpperLeft,
    );
    fig
}

fn draw_target(
    ax: &mut rizzma::Axes,
    t: &TargetSeries,
    tolerance: f64,
    pass_min_kev: f64,
    x_lim: (f64, f64),
) {
    let e: Vec<f64> = t.measured.iter().map(|m| m.energy_kev).collect();
    let lo: Vec<f64> = t.measured.iter().map(|m| m.min).collect();
    let hi: Vec<f64> = t.measured.iter().map(|m| m.max).collect();
    let med: Vec<f64> = t.measured.iter().map(|m| m.median).collect();
    ax.fill_between_with_color(&e, &lo, &hi, BAND);
    ax.add_line(
        Line2D::new(e.clone(), med)
            .with_color(MEDIAN)
            .with_linewidth(1.2),
    );

    // Pass band: median ± tolerance where the rule applies.
    let pass: Vec<&MeasuredPoint> = t
        .measured
        .iter()
        .filter(|m| m.energy_kev >= pass_min_kev)
        .collect();
    let pe: Vec<f64> = pass.iter().map(|m| m.energy_kev).collect();
    for sign in [-1.0, 1.0] {
        let y: Vec<f64> = pass.iter().map(|m| m.median + sign * tolerance).collect();
        ax.add_line(
            Line2D::new(pe.clone(), y)
                .with_color(PASS)
                .with_linewidth(1.0)
                .with_dashes(Some((0.0, vec![4.0, 3.0]))),
        );
    }

    // lindhard η with ±σ bars (σ is small at 1e5 primaries; the bars may be
    // shorter than the marker).
    let cap = (x_lim.1 - x_lim.0) * 0.008;
    for r in &t.lindhard {
        let (y0, y1) = (r.eta - r.eta_se, r.eta + r.eta_se);
        ax.add_line(Line2D::new(vec![r.energy_kev; 2], vec![y0, y1]).with_color(LINDHARD));
        for y in [y0, y1] {
            ax.add_line(
                Line2D::new(vec![r.energy_kev - cap, r.energy_kev + cap], vec![y, y])
                    .with_color(LINDHARD),
            );
        }
    }
    let x: Vec<f64> = t.lindhard.iter().map(|r| r.energy_kev).collect();
    let y: Vec<f64> = t.lindhard.iter().map(|r| r.eta).collect();
    // `scatter` takes the next cycle colour; restyle the collection it returns.
    let markers = ax.scatter(&x, &y);
    *markers = markers.clone().with_facecolors(vec![LINDHARD]);

    let y_top = t
        .measured
        .iter()
        .map(|m| {
            let band = if m.energy_kev >= pass_min_kev {
                m.median + tolerance
            } else {
                m.median
            };
            m.max.max(band)
        })
        .chain(t.lindhard.iter().map(|r| r.eta + r.eta_se))
        .fold(0.0_f64, f64::max);
    let y_max = (y_top / Y_STEP).ceil() * Y_STEP;
    ax.set_xlim(x_lim.0, x_lim.1);
    ax.set_ylim(0.0, y_max);
    ax.set_title(t.target.clone());
    ax.set_xlabel("E (keV)");
    ax.set_ylabel("\u{3b7}");
    ax.grid(true);
}
