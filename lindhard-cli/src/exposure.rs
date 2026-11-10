//! The `lindhard exposure` command: the dose of a shot list from a
//! Gaussian-sum PSF, with optional iso-dose contours. Model and conventions:
//! [`lindhard::exposure`] and [`lindhard::contour`]; files: `docs/cli.md`.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::Args;
use lindhard::contour::{contours, Contour};
use lindhard::exposure::{parse_shots_csv, uniform_area_dose, Exposure, Grid, DEFAULT_SUBSAMPLES};
use lindhard::tally::GaussianPsf;
use serde_json::json;

/// File names written in the output directory.
pub const DOSE_FILE: &str = "exposure_dose.csv";
/// Contour vertices (with `--threshold` or `--threshold-fraction`).
pub const CONTOURS_FILE: &str = "exposure_contours.csv";
/// One row per contour: area and centroid.
pub const FEATURES_FILE: &str = "exposure_features.csv";
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
    write(&a.out, DOSE_FILE, &dose.csv())?;

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
    let mut features = serde_json::Value::Null;
    if let Some(t) = threshold {
        let c = contours(&dose, t);
        write(&a.out, CONTOURS_FILE, &Contour::csv(&c))?;
        write(&a.out, FEATURES_FILE, &features_csv(&c))?;
        features = json!({
            "threshold_ev_per_m2": t,
            "contours": c.len(),
            "closed": c.iter().filter(|c| c.closed).count(),
        });
    }
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
    });
    write(
        &a.out,
        SUMMARY_FILE,
        &(serde_json::to_string_pretty(&summary)? + "\n"),
    )?;
    eprintln!(
        "{n_shots} shots on a {}x{} grid; wrote {}",
        a.nx,
        a.ny,
        a.out.display()
    );
    Ok(())
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

fn write(dir: &Path, name: &str, text: &str) -> Result<()> {
    let p = dir.join(name);
    std::fs::write(&p, text).with_context(|| format!("writing {}", p.display()))
}
