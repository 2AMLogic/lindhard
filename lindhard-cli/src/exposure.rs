//! The `lindhard exposure` command: the dose of a shot list from a
//! Gaussian-sum PSF, with optional iso-dose contours. Model and conventions:
//! [`lindhard::exposure`] and [`lindhard::contour`]; files: `docs/cli.md`.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::Args;
use lindhard::contour::{contours, Contour};
use lindhard::exposure::{
    compare_with_convolution, comparison_csv, direct_dose, parse_shots_csv, uniform_area_dose,
    DirectConfig, Exposure, Grid, DEFAULT_SUBSAMPLES,
};
use lindhard::input::electron::ElectronInput;
use lindhard::tally::GaussianPsf;
use serde_json::json;

use crate::publish::OutputSet;

/// File names written in the output directory.
pub const DOSE_FILE: &str = "exposure_dose.csv";
/// Contour vertices (with `--threshold` or `--threshold-fraction`).
pub const CONTOURS_FILE: &str = "exposure_contours.csv";
/// One row per contour: area and centroid.
pub const FEATURES_FILE: &str = "exposure_features.csv";
/// With `--direct`: the transport dose and its standard error per point.
pub const DIRECT_FILE: &str = "exposure_direct.csv";
/// With `--direct`: convolution and direct dose side by side, with the
/// deviation in combined standard errors and the flag.
pub const COMPARISON_FILE: &str = "exposure_comparison.csv";
/// The inputs echoed, the area dose and the contour count.
pub const SUMMARY_FILE: &str = "exposure_summary.json";

/// Arguments of `lindhard exposure`. Lengths are in metres, doses in eV/m².
#[derive(Args, Debug)]
pub struct ExposureArgs {
    /// Shot list CSV (columns x_m, y_m, weight; optional blur_m, shape,
    /// size_m, size2_m).
    #[arg(long)]
    pub shots: PathBuf,
    /// Read the PSF from an `electron_psf_parameters.csv` (see --model).
    #[arg(long, value_name = "CSV", conflicts_with_all = ["total_ev", "alpha_m", "beta_m", "eta", "gamma_m", "nu"])]
    pub psf_parameters: Option<PathBuf>,
    /// Which fit of the parameters file to use (`double_gaussian` or
    /// `triple_gaussian`). Default: the triple fit if present, else double.
    #[arg(long, requires = "psf_parameters")]
    pub model: Option<String>,
    /// PSF energy E per primary, eV (typed-in PSF).
    #[arg(long)]
    pub total_ev: Option<f64>,
    /// Forward width alpha, m.
    #[arg(long)]
    pub alpha_m: Option<f64>,
    /// Backscatter width beta, m.
    #[arg(long)]
    pub beta_m: Option<f64>,
    /// Backscatter ratio eta.
    #[arg(long)]
    pub eta: Option<f64>,
    /// Third width gamma, m (triple Gaussian).
    #[arg(long, requires = "nu")]
    pub gamma_m: Option<f64>,
    /// Third ratio nu (triple Gaussian).
    #[arg(long, requires = "gamma_m")]
    pub nu: Option<f64>,
    /// Beam blur b in exp(-r^2/b^2), m.
    #[arg(long, default_value_t = 0.0)]
    pub blur_m: f64,
    /// Sub-samples per axis of a shaped spot.
    #[arg(long, default_value_t = DEFAULT_SUBSAMPLES)]
    pub subsamples: usize,
    /// Grid: x of the first point, m.
    #[arg(long, allow_negative_numbers = true)]
    pub x0_m: f64,
    /// Grid: y of the first point, m.
    #[arg(long, allow_negative_numbers = true)]
    pub y0_m: f64,
    /// Grid spacing in x, m.
    #[arg(long)]
    pub dx_m: f64,
    /// Grid spacing in y, m.
    #[arg(long)]
    pub dy_m: f64,
    /// Grid points in x.
    #[arg(long)]
    pub nx: usize,
    /// Grid points in y.
    #[arg(long)]
    pub ny: usize,
    /// Iso-dose threshold, eV/m^2.
    #[arg(long, conflicts_with = "threshold_fraction")]
    pub threshold: Option<f64>,
    /// Iso-dose threshold as a fraction of the largest dose on the grid.
    #[arg(long)]
    pub threshold_fraction: Option<f64>,
    /// Also run the exposure by electron transport (direct Monte Carlo) on
    /// the target of this electron input (`[electron]` TOML) and compare it
    /// with the convolution. Only the input's layers, tables and beam energy
    /// are used; the shot list sets where the primaries land. The target is
    /// the laterally uniform layer stack of the input.
    #[arg(long, value_name = "TOML")]
    pub direct: Option<PathBuf>,
    /// Direct mode: primaries per shot, summed over the replicas.
    #[arg(long, default_value_t = 80_000, requires = "direct")]
    pub direct_histories: u64,
    /// Direct mode: independent replicas for the statistical error (at least 2).
    #[arg(long, default_value_t = 8, requires = "direct")]
    pub direct_replicas: usize,
    /// Direct mode: base random seed.
    #[arg(long, default_value_t = 1, requires = "direct")]
    pub direct_seed: u64,
    /// Direct mode: lower depth of the film slab, m.
    #[arg(long, default_value_t = 0.0, requires = "direct")]
    pub slab_lo_m: f64,
    /// Direct mode: upper depth of the film slab, m (required).
    #[arg(long, requires = "direct")]
    pub slab_hi_m: Option<f64>,
    /// Direct mode: flag a point when the two doses differ by more than this
    /// many combined standard errors.
    #[arg(long, default_value_t = 5.0, requires = "direct")]
    pub tolerance_sigmas: f64,
    /// Direct mode: relative uncertainty assigned to the convolution dose.
    #[arg(long, default_value_t = 0.0, requires = "direct")]
    pub conv_rel_err: f64,
    /// Direct mode: points below this fraction of the largest dose are not
    /// judged (the run does not resolve them).
    #[arg(long, default_value_t = 0.01, requires = "direct")]
    pub floor_fraction: f64,
    /// Output directory (created if missing).
    #[arg(long, short)]
    pub out: PathBuf,
    /// Worker threads. Never changes the results.
    #[arg(long)]
    pub threads: Option<usize>,
}

/// The PSF fit named `model` (or the best available) in a parameters CSV.
pub fn psf_from_parameters_csv(text: &str, model: Option<&str>) -> Result<GaussianPsf> {
    let rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|l| l.split(',').map(str::trim).collect())
        .filter(|r: &Vec<&str>| r.len() >= 3)
        .collect();
    let get = |m: &str, name: &str| -> Result<Option<f64>> {
        rows.iter()
            .find(|r| r[0] == m && r[1] == name)
            .map(|r| {
                r[2].parse::<f64>()
                    .with_context(|| format!("{m} {name}: `{}` is not a number", r[2]))
            })
            .transpose()
    };
    let model = match model {
        Some(m) => m,
        None if get("triple_gaussian", "total_ev")?.is_some() => "triple_gaussian",
        None => "double_gaussian",
    };
    let need = |name: &str| {
        get(model, name)?.with_context(|| format!("no `{model}` parameter `{name}` in the file"))
    };
    let psf = match model {
        "double_gaussian" => GaussianPsf::double(
            need("total_ev")?,
            need("alpha_m")?,
            need("beta_m")?,
            need("eta")?,
        ),
        "triple_gaussian" => GaussianPsf::triple(
            need("total_ev")?,
            need("alpha_m")?,
            need("beta_m")?,
            need("eta")?,
            need("gamma_m")?,
            need("nu")?,
        ),
        other => bail!("unknown model `{other}` (double_gaussian or triple_gaussian)"),
    };
    Ok(psf?)
}

fn read_psf(a: &ExposureArgs) -> Result<GaussianPsf> {
    if let Some(p) = &a.psf_parameters {
        let text =
            std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
        return psf_from_parameters_csv(&text, a.model.as_deref())
            .with_context(|| format!("{}: PSF parameters", p.display()));
    }
    let (Some(e), Some(al), Some(be), Some(et)) = (a.total_ev, a.alpha_m, a.beta_m, a.eta) else {
        bail!("give --psf-parameters, or --total-ev, --alpha-m, --beta-m and --eta");
    };
    Ok(match (a.gamma_m, a.nu) {
        (Some(g), Some(n)) => GaussianPsf::triple(e, al, be, et, g, n)?,
        _ => GaussianPsf::double(e, al, be, et)?,
    })
}

/// Runs the command.
pub fn run(a: &ExposureArgs) -> Result<()> {
    if a.threads == Some(0) {
        bail!("--threads must be at least 1");
    }
    let psf = read_psf(a)?;
    let text = std::fs::read_to_string(&a.shots)
        .with_context(|| format!("reading {}", a.shots.display()))?;
    let shots = parse_shots_csv(&text).with_context(|| format!("{}", a.shots.display()))?;
    let n_shots = shots.len();
    let direct_shots = a.direct.is_some().then(|| shots.clone());
    let exposure = Exposure::with_subsamples(psf, shots, a.blur_m, a.subsamples)?;
    let grid = Grid {
        x0_m: a.x0_m,
        y0_m: a.y0_m,
        dx_m: a.dx_m,
        dy_m: a.dy_m,
        nx: a.nx,
        ny: a.ny,
    };
    let dose = match a.threads {
        Some(n) => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()?
            .install(|| exposure.dose_on_grid(&grid))?,
        None => exposure.dose_on_grid(&grid)?,
    };
    std::fs::create_dir_all(&a.out).with_context(|| format!("creating {}", a.out.display()))?;

    let max = dose
        .dose_ev_per_m2
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let threshold = match (a.threshold, a.threshold_fraction) {
        (Some(t), _) => Some(t),
        (None, Some(f)) => Some(f * max),
        _ => None,
    };
    // Every output is computed before the first write, so the summary can be
    // serialized up front and the whole set published through `OutputSet`.
    let mut features = serde_json::Value::Null;
    let mut contour_rows = None;
    let mut feature_rows = None;
    if let Some(t) = threshold {
        let c = contours(&dose, t);
        contour_rows = Some(Contour::csv(&c));
        feature_rows = Some(features_csv(&c));
        features = json!({
            "threshold_ev_per_m2": t,
            "contours": c.len(),
            "closed": c.iter().filter(|c| c.closed).count(),
        });
    }
    let (direct_rows, comparison_rows, direct) = match (&a.direct, direct_shots) {
        (Some(toml), Some(shots)) => {
            let d = run_direct(a, toml, &shots, &grid, &exposure)?;
            (Some(d.dose_csv), Some(d.comparison_csv), d.summary)
        }
        _ => (None, None, serde_json::Value::Null),
    };
    // The area dose per unit weight density: the interior value of a large
    // uniform exposure is `weight_per_m2 * E` (the forward part is its alpha share).
    let unit = uniform_area_dose(&psf, 1.0);
    let summary = json!({
        "shots": n_shots,
        "psf": psf,
        "blur_m": a.blur_m,
        "subsamples": a.subsamples,
        "grid": { "x0_m": a.x0_m, "y0_m": a.y0_m, "dx_m": a.dx_m, "dy_m": a.dy_m,
                  "nx": a.nx, "ny": a.ny },
        "max_dose_ev_per_m2": max,
        "area_dose_per_unit_weight_density": {
            "total_ev_per_m2": unit.total_ev_per_m2,
            "forward_ev_per_m2": unit.forward_ev_per_m2,
            "long_range_ev_per_m2": unit.long_range_ev_per_m2,
        },
        "contours": features,
        "direct": direct,
    });
    let set = OutputSet::begin(
        &a.out,
        SUMMARY_FILE,
        serde_json::to_string_pretty(&summary)? + "\n",
    )?;
    set.write(DOSE_FILE, dose.csv())?;
    set.reconcile_optional(CONTOURS_FILE, contour_rows)?;
    set.reconcile_optional(FEATURES_FILE, feature_rows)?;
    set.reconcile_optional(DIRECT_FILE, direct_rows)?;
    set.reconcile_optional(COMPARISON_FILE, comparison_rows)?;
    // Last, so the summary describes the completed output set.
    set.publish()?;
    eprintln!(
        "{n_shots} shots on a {}x{} grid; wrote {}",
        a.nx,
        a.ny,
        a.out.display()
    );
    Ok(())
}

/// Direct mode's outputs: the bodies of [`DIRECT_FILE`] and
/// [`COMPARISON_FILE`], and the `direct` block of the summary.
struct Direct {
    dose_csv: String,
    comparison_csv: String,
    summary: serde_json::Value,
}

/// Direct mode: transport the shot list on the target of `toml` and return the
/// direct dose, the comparison and the summary block. Nothing is written here;
/// the caller publishes the whole output set.
fn run_direct(
    a: &ExposureArgs,
    toml: &Path,
    shots: &[lindhard::exposure::Shot],
    grid: &Grid,
    exposure: &Exposure,
) -> Result<Direct> {
    let Some(slab_hi_m) = a.slab_hi_m else {
        bail!("--direct needs --slab-hi-m (the film's lower and upper depth bound the tally)");
    };
    let text =
        std::fs::read_to_string(toml).with_context(|| format!("reading {}", toml.display()))?;
    let input = ElectronInput::from_toml_str(&text)
        .with_context(|| format!("{}: invalid input", toml.display()))?;
    let base = toml
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let resolved = input
        .resolve_in(base)
        .with_context(|| format!("{}: invalid input", toml.display()))?;
    for w in &resolved.warnings {
        eprintln!("warning: {w}");
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(a.threads.unwrap_or(0))
        .build()?;
    let (transport, _tables, _secs) = crate::electron::build_transport(&resolved, &pool, None)?;
    let cfg = DirectConfig {
        histories_per_shot: a.direct_histories,
        replicas: a.direct_replicas,
        chunk_size: crate::electron::CHUNK_SIZE,
        seed: a.direct_seed,
        depth_lo_m: a.slab_lo_m,
        depth_hi_m: slab_hi_m,
    };
    let d =
        pool.install(|| direct_dose(&transport, &resolved.primary, shots, a.blur_m, grid, &cfg))?;
    let cmp = compare_with_convolution(
        &d,
        exposure,
        a.conv_rel_err,
        a.tolerance_sigmas,
        a.floor_fraction,
    )?;
    let flagged = cmp.iter().filter(|c| c.flagged).count();
    eprintln!(
        "direct: {flagged} of {} points differ from the convolution by more than {} combined standard errors",
        cmp.len(),
        a.tolerance_sigmas
    );
    Ok(Direct {
        dose_csv: d.csv(),
        comparison_csv: comparison_csv(&cmp),
        summary: json!({
            "target": toml.display().to_string(),
            "histories_per_shot": a.direct_histories,
            "replicas": d.replicas,
            "histories_per_shot_per_replica": d.histories_per_shot_per_replica,
            "seed": a.direct_seed,
            "slab_m": [a.slab_lo_m, slab_hi_m],
            "tolerance_sigmas": a.tolerance_sigmas,
            "conv_rel_err": a.conv_rel_err,
            "floor_fraction": a.floor_fraction,
            "points": cmp.len(),
            "flagged": flagged,
        }),
    })
}

fn features_csv(c: &[Contour]) -> String {
    let f = |o: Option<f64>| o.map(|v| v.to_string()).unwrap_or_default();
    let mut s = String::from("contour,closed,vertices,signed_area_m2,centroid_x_m,centroid_y_m\n");
    for (k, c) in c.iter().enumerate() {
        s.push_str(&format!(
            "{k},{},{},{},{},{}\n",
            u8::from(c.closed),
            c.points.len(),
            f(c.signed_area_m2),
            f(c.centroid_m.map(|m| m[0])),
            f(c.centroid_m.map(|m| m[1])),
        ));
    }
    s
}
