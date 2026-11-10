//! Multi-shot exposure: the energy density deposited by a list of shots, from
//! a fitted Gaussian-sum point-spread function ([`GaussianPsf`]).
//!
//! # The model
//!
//! A shot is a beam landing at `(x, y)` with a weight (the number of
//! primaries, or any relative dwell) and an optional spot shape and beam
//! blur. The deposited energy density at a point is the sum over shots of the
//! shot's weight times the PSF, `E f(r)` (eV/m²), centred on the shot, where
//! `E` is [`GaussianPsf::total_ev`], the energy one primary deposits in the
//! film. The PSF is the sum of Gaussians documented in
//! [`crate::tally::psf`]:
//!
//! ```text
//! f(r) = 1/(π Σc) Σ_t (c_t / w_t²) exp(-r²/w_t²)
//! ```
//!
//! with `(w_t, c_t) = (α, 1), (β, η)` and, for the triple form, `(γ, ν)`.
//!
//! **Beam blur.** A Gaussian beam profile of width `b` in the same
//! convention, `exp(-r²/b²)/(π b²)`, is convolved with the PSF. The
//! convolution of two normalised 2D Gaussians is the normalised Gaussian
//! whose `w²` is the sum (the variances of independent Gaussians add; here
//! `w² = 2σ²`, so the `w²` add too):
//!
//! ```text
//! ∫ exp(-|r-s|²/w²)/(π w²) · exp(-s²/b²)/(π b²) d²s = exp(-r²/(w²+b²))/(π (w²+b²))
//! ```
//!
//! This is textbook (the convolution theorem for Gaussians; it follows from
//! completing the square in the exponent), so a blurred point shot needs no
//! numerical convolution: every term's `w_t²` becomes `w_t² + b²`, and the
//! weights `c_t` and the normalisation are unchanged. The blur of a shot is
//! `b² = b_global² + b_shot²`.
//!
//! **Shaped spots.** A uniform disc or rectangle is a fixed sub-sampling
//! into point shots: the bounding box of the spot is cut into
//! `n × n` equal cells (`n` is [`Exposure::subsamples`], default
//! [`DEFAULT_SUBSAMPLES`]); each cell centre inside the spot becomes a point
//! shot with an equal share of the weight (the shares sum to the shot's
//! weight). The rectangle is the midpoint rule on the exact shape, so its
//! error is `O(h²)` in the cell pitch `h` relative to the narrowest Gaussian
//! width `sqrt(w_t² + b²)`; it is accurate when `h` is below that width. A
//! disc has the same quadrature error and, in addition, the boundary of the
//! sampled cells approximates the circle to `O(h)` in area; raise
//! `subsamples` for a spot that is large compared with the PSF's `α`.
//!
//! # Determinism
//!
//! The dose at one point sums the point sources in list order (the shots in
//! the order given, each shot's sub-samples in row-major order). A grid
//! parallelises over grid points only, so each point's sum is the same fixed
//! sequence of floating-point operations: the result is bit-identical at any
//! thread count.
//!
//! # Accuracy
//!
//! The dose is only as good as the PSF it uses. A fit that does not describe
//! the simulated profile (see the reduced `χ²` of [`PsfFit`]) gives a dose
//! that does not describe the simulated exposure either; the fit's validation
//! is tracked separately (issue #151). Superposition also assumes the film is
//! laterally uniform.
//!
//! [`PsfFit`]: crate::tally::PsfFit

use std::f64::consts::PI;

use rayon::prelude::*;
use thiserror::Error;

use crate::tally::GaussianPsf;

/// Default sub-sampling per axis of a shaped spot.
pub const DEFAULT_SUBSAMPLES: usize = 16;

/// Errors of building or reading an exposure.
#[derive(Debug, Error)]
pub enum ExposureError {
    /// A parameter is out of range.
    #[error("invalid exposure: {0}")]
    Invalid(String),
    /// The shot-list CSV could not be read.
    #[error("shot list line {line}: {message}")]
    Csv {
        /// 1-based line number.
        line: usize,
        /// What is wrong.
        message: String,
    },
}

/// The shape of a shot's spot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpotShape {
    /// A point (a Gaussian beam of the blur width only).
    Point,
    /// A uniform disc of radius `radius_m`.
    Disc {
        /// Radius, m.
        radius_m: f64,
    },
    /// A uniform axis-aligned rectangle centred on the shot.
    Rect {
        /// Extent in x, m.
        width_m: f64,
        /// Extent in y, m.
        height_m: f64,
    },
}

/// One shot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Shot {
    /// Position, m.
    pub x_m: f64,
    /// Position, m.
    pub y_m: f64,
    /// Number of primaries (or relative dwell): the PSF is scaled by it.
    pub weight: f64,
    /// Additional Gaussian beam width `b` of this shot, m, in the
    /// `exp(-r²/b²)` convention; adds in quadrature to the exposure's blur.
    pub blur_m: f64,
    /// Spot shape.
    pub shape: SpotShape,
}

impl Shot {
    /// A point shot with no extra blur.
    pub fn point(x_m: f64, y_m: f64, weight: f64) -> Self {
        Self {
            x_m,
            y_m,
            weight,
            blur_m: 0.0,
            shape: SpotShape::Point,
        }
    }

    fn validate(&self) -> Result<(), ExposureError> {
        let bad = |what: &str| Err(ExposureError::Invalid(format!("shot {what}: {self:?}")));
        if !(self.x_m.is_finite() && self.y_m.is_finite() && self.weight.is_finite()) {
            return bad("position and weight must be finite");
        }
        if !(self.blur_m.is_finite() && self.blur_m >= 0.0) {
            return bad("blur must be finite and non-negative");
        }
        let positive = |v: f64| v.is_finite() && v > 0.0;
        let ok = match self.shape {
            SpotShape::Point => true,
            SpotShape::Disc { radius_m } => positive(radius_m),
            SpotShape::Rect { width_m, height_m } => positive(width_m) && positive(height_m),
        };
        if ok {
            Ok(())
        } else {
            bad("spot size must be positive and finite")
        }
    }
}

/// One Gaussian point source after expanding shaped spots.
#[derive(Debug, Clone, Copy)]
struct Source {
    x: f64,
    y: f64,
    weight: f64,
    /// Blur variance `b²`, m².
    b2: f64,
}

/// A regular grid of points: point `(i, j)` is at
/// `(x0 + i dx, y0 + j dy)`, stored row-major (`j * nx + i`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grid {
    /// x of point `(0, 0)`, m.
    pub x0_m: f64,
    /// y of point `(0, 0)`, m.
    pub y0_m: f64,
    /// Spacing in x, m.
    pub dx_m: f64,
    /// Spacing in y, m.
    pub dy_m: f64,
    /// Points in x (at least 1).
    pub nx: usize,
    /// Points in y (at least 1).
    pub ny: usize,
}

impl Grid {
    /// Checks the spacings are positive and finite and the counts non-zero.
    pub fn validate(&self) -> Result<(), ExposureError> {
        let ok = self.x0_m.is_finite()
            && self.y0_m.is_finite()
            && self.dx_m.is_finite()
            && self.dx_m > 0.0
            && self.dy_m.is_finite()
            && self.dy_m > 0.0
            && self.nx >= 1
            && self.ny >= 1;
        if ok {
            Ok(())
        } else {
            Err(ExposureError::Invalid(format!("grid: {self:?}")))
        }
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.nx * self.ny
    }

    /// True if the grid has no points.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Position of point `(i, j)`, m.
    pub fn point(&self, i: usize, j: usize) -> [f64; 2] {
        [
            self.x0_m + i as f64 * self.dx_m,
            self.y0_m + j as f64 * self.dy_m,
        ]
    }
}

/// The dose on a [`Grid`], eV/m², row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct DoseGrid {
    /// The grid.
    pub grid: Grid,
    /// Deposited energy density per point, eV/m².
    pub dose_ev_per_m2: Vec<f64>,
}

impl DoseGrid {
    /// The dose at point `(i, j)`.
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.dose_ev_per_m2[j * self.grid.nx + i]
    }

    /// CSV with columns `x_m,y_m,dose_ev_per_m2`, one row per point, rows in
    /// grid order.
    pub fn csv(&self) -> String {
        let mut s = String::from("x_m,y_m,dose_ev_per_m2\n");
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let [x, y] = self.grid.point(i, j);
                s.push_str(&format!("{x},{y},{}\n", self.at(i, j)));
            }
        }
        s
    }
}

/// A list of shots on a PSF, ready to evaluate.
#[derive(Debug, Clone)]
pub struct Exposure {
    psf: GaussianPsf,
    blur_m: f64,
    subsamples: usize,
    shots: Vec<Shot>,
    sources: Vec<Source>,
    /// `(w_t², c_t / Σc)` of the PSF.
    terms: Vec<(f64, f64)>,
}

impl Exposure {
    /// An exposure of `shots` on `psf` with a global beam blur `blur_m`
    /// (`b` in `exp(-r²/b²)`, 0 for none) and [`DEFAULT_SUBSAMPLES`].
    pub fn new(psf: GaussianPsf, shots: Vec<Shot>, blur_m: f64) -> Result<Self, ExposureError> {
        Self::with_subsamples(psf, shots, blur_m, DEFAULT_SUBSAMPLES)
    }

    /// As [`Exposure::new`], with `subsamples` cells per axis for shaped
    /// spots (at least 1).
    pub fn with_subsamples(
        psf: GaussianPsf,
        shots: Vec<Shot>,
        blur_m: f64,
        subsamples: usize,
    ) -> Result<Self, ExposureError> {
        psf.validate()
            .map_err(|e| ExposureError::Invalid(e.to_string()))?;
        if !(blur_m.is_finite() && blur_m >= 0.0) {
            return Err(ExposureError::Invalid(format!(
                "blur must be finite and non-negative, got {blur_m}"
            )));
        }
        if subsamples == 0 {
            return Err(ExposureError::Invalid(
                "subsamples must be at least 1".into(),
            ));
        }
        for s in &shots {
            s.validate()?;
        }
        let terms = psf_terms(&psf);
        let sources = expand(&shots, blur_m, subsamples);
        Ok(Self {
            psf,
            blur_m,
            subsamples,
            shots,
            sources,
            terms,
        })
    }

    /// The PSF.
    pub fn psf(&self) -> &GaussianPsf {
        &self.psf
    }

    /// The global beam blur `b`, m.
    pub fn blur_m(&self) -> f64 {
        self.blur_m
    }

    /// Sub-samples per axis of a shaped spot.
    pub fn subsamples(&self) -> usize {
        self.subsamples
    }

    /// The shots, in the order they are summed.
    pub fn shots(&self) -> &[Shot] {
        &self.shots
    }

    /// The deposited energy density at `(x, y)`, eV/m².
    pub fn dose_at(&self, x_m: f64, y_m: f64) -> f64 {
        let mut sum = 0.0;
        for s in &self.sources {
            let d2 = (x_m - s.x) * (x_m - s.x) + (y_m - s.y) * (y_m - s.y);
            let mut f = 0.0;
            for &(w2, c) in &self.terms {
                let v = w2 + s.b2;
                f += c / v * (-d2 / v).exp();
            }
            sum += s.weight * f;
        }
        self.psf.total_ev * sum / PI
    }

    /// The dose at each of `points`, in order. Parallel over points; the
    /// result does not depend on the thread count.
    pub fn dose_at_points(&self, points: &[[f64; 2]]) -> Vec<f64> {
        points
            .par_iter()
            .map(|p| self.dose_at(p[0], p[1]))
            .collect()
    }

    /// The dose on a grid. Parallel over grid points; the result does not
    /// depend on the thread count.
    pub fn dose_on_grid(&self, grid: &Grid) -> Result<DoseGrid, ExposureError> {
        grid.validate()?;
        let nx = grid.nx;
        let dose: Vec<f64> = (0..grid.len())
            .into_par_iter()
            .map(|k| {
                let [x, y] = grid.point(k % nx, k / nx);
                self.dose_at(x, y)
            })
            .collect();
        Ok(DoseGrid {
            grid: *grid,
            dose_ev_per_m2: dose,
        })
    }
}

/// `(w_t², c_t / Σc)` of a PSF, the `α` term first.
fn psf_terms(psf: &GaussianPsf) -> Vec<(f64, f64)> {
    let mut t = vec![(psf.alpha_m, 1.0), (psf.beta_m, psf.eta)];
    if let (Some(g), Some(n)) = (psf.gamma_m, psf.nu) {
        t.push((g, n));
    }
    let sum: f64 = t.iter().map(|x| x.1).sum();
    t.into_iter().map(|(w, c)| (w * w, c / sum)).collect()
}

/// Shots to point sources, in list order (see the module docs).
fn expand(shots: &[Shot], global_blur_m: f64, n: usize) -> Vec<Source> {
    let mut out = Vec::with_capacity(shots.len());
    for s in shots {
        let b2 = global_blur_m * global_blur_m + s.blur_m * s.blur_m;
        let (hx, hy) = match s.shape {
            SpotShape::Point => {
                out.push(Source {
                    x: s.x_m,
                    y: s.y_m,
                    weight: s.weight,
                    b2,
                });
                continue;
            }
            SpotShape::Disc { radius_m } => (radius_m, radius_m),
            SpotShape::Rect { width_m, height_m } => (width_m / 2.0, height_m / 2.0),
        };
        let mut cells = Vec::with_capacity(n * n);
        for j in 0..n {
            for i in 0..n {
                let ox = -hx + (i as f64 + 0.5) * (2.0 * hx / n as f64);
                let oy = -hy + (j as f64 + 0.5) * (2.0 * hy / n as f64);
                let inside = match s.shape {
                    SpotShape::Disc { radius_m } => ox * ox + oy * oy <= radius_m * radius_m,
                    _ => true,
                };
                if inside {
                    cells.push((ox, oy));
                }
            }
        }
        let share = s.weight / cells.len() as f64;
        out.extend(cells.into_iter().map(|(ox, oy)| Source {
            x: s.x_m + ox,
            y: s.y_m + oy,
            weight: share,
            b2,
        }));
    }
    out
}

/// The dose far inside a large uniform exposure, split by PSF term.
#[derive(Debug, Clone, PartialEq)]
pub struct AreaDose {
    /// The total, eV/m²: `weight_per_m2 · E`. Independent of every width, as
    /// the PSF integrates to `E` (beam blur does not change it either).
    pub total_ev_per_m2: f64,
    /// The forward-scattering (`α`) term's share, eV/m².
    pub forward_ev_per_m2: f64,
    /// The long-range terms' share (`β`, and `γ` for the triple form),
    /// eV/m²: the backscatter contribution.
    pub long_range_ev_per_m2: f64,
}

/// The area dose of a uniform exposure of `weight_per_m2` primaries per m²
/// over an area much larger than the PSF's longest width, with the
/// forward-only part reported next to the total.
///
/// Each normalised Gaussian term integrates to 1 over the plane, so a term
/// with share `c_t/Σc` contributes `weight_per_m2 · E · c_t/Σc` and the
/// total is `weight_per_m2 · E`.
pub fn uniform_area_dose(psf: &GaussianPsf, weight_per_m2: f64) -> AreaDose {
    let t = psf_terms(psf);
    let total = weight_per_m2 * psf.total_ev;
    let forward = total * t[0].1;
    AreaDose {
        total_ev_per_m2: total,
        forward_ev_per_m2: forward,
        long_range_ev_per_m2: total - forward,
    }
}

/// Reads a shot list.
///
/// A CSV with a header row. Required columns `x_m`, `y_m`, `weight`;
/// optional `blur_m` (default 0), `shape` (`point`, the default, `disc` or
/// `rect`), `size_m` (disc radius, or rectangle width) and `size2_m`
/// (rectangle height, default `size_m`). Blank lines and lines starting with
/// `#` are skipped. Lengths are in metres.
pub fn parse_shots_csv(text: &str) -> Result<Vec<Shot>, ExposureError> {
    let err = |line: usize, message: String| ExposureError::Csv { line, message };
    let mut lines = text
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.trim()))
        .filter(|(_, l)| !l.is_empty() && !l.starts_with('#'));
    let (hline, header) = lines.next().ok_or_else(|| err(1, "no header row".into()))?;
    let names: Vec<&str> = header.split(',').map(str::trim).collect();
    let col = |n: &str| names.iter().position(|c| *c == n);
    let need = |n: &str| col(n).ok_or_else(|| err(hline, format!("missing column `{n}`")));
    let (cx, cy, cw) = (need("x_m")?, need("y_m")?, need("weight")?);
    let (cb, cs, c1, c2) = (col("blur_m"), col("shape"), col("size_m"), col("size2_m"));
    let mut shots = Vec::new();
    for (n, line) in lines {
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        if f.len() != names.len() {
            return Err(err(
                n,
                format!("expected {} fields, found {}", names.len(), f.len()),
            ));
        }
        let num = |i: usize| {
            f[i].parse::<f64>()
                .map_err(|_| err(n, format!("`{}` is not a number", f[i])))
        };
        let opt = |c: Option<usize>| -> Result<Option<f64>, ExposureError> {
            match c {
                Some(i) if !f[i].is_empty() => num(i).map(Some),
                _ => Ok(None),
            }
        };
        let size = opt(c1)?;
        let size2 = opt(c2)?;
        let need_size = |what: &str| size.ok_or_else(|| err(n, format!("{what} needs size_m")));
        let shape = match cs.map(|i| f[i]).unwrap_or("point") {
            "" | "point" => SpotShape::Point,
            "disc" => SpotShape::Disc {
                radius_m: need_size("disc")?,
            },
            "rect" => {
                let w = need_size("rect")?;
                SpotShape::Rect {
                    width_m: w,
                    height_m: size2.unwrap_or(w),
                }
            }
            other => return Err(err(n, format!("unknown shape `{other}`"))),
        };
        let shot = Shot {
            x_m: num(cx)?,
            y_m: num(cy)?,
            weight: num(cw)?,
            blur_m: opt(cb)?.unwrap_or(0.0),
            shape,
        };
        shot.validate().map_err(|e| err(n, e.to_string()))?;
        shots.push(shot);
    }
    Ok(shots)
}
