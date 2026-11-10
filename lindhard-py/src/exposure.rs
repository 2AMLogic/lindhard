//! Exposure functions: the dose of a shot list from a Gaussian-sum PSF, and
//! iso-dose contours of a gridded dose. Thin wrappers over
//! [`lindhard::exposure`] and [`lindhard::contour`], the same code the
//! `lindhard exposure` command runs, so the arrays equal its CSV values bit
//! for bit at any thread count.

use lindhard::contour::{contours_of, Contour};
use lindhard::exposure::{
    parse_shots_csv, uniform_area_dose, Exposure, Grid, Shot, DEFAULT_SUBSAMPLES,
};
use lindhard::tally::GaussianPsf;
use numpy::{IntoPyArray, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::errors;

/// The PSF from keyword values; `gamma_m` and `nu` give the triple form.
fn psf(
    total_ev: f64,
    alpha_m: f64,
    beta_m: f64,
    eta: f64,
    gamma_m: Option<f64>,
    nu: Option<f64>,
) -> PyResult<GaussianPsf> {
    match (gamma_m, nu) {
        (Some(g), Some(n)) => GaussianPsf::triple(total_ev, alpha_m, beta_m, eta, g, n),
        (None, None) => GaussianPsf::double(total_ev, alpha_m, beta_m, eta),
        _ => {
            return Err(errors::input(
                "gamma_m and nu must both be given or both omitted",
            ))
        }
    }
    .map_err(errors::input)
}

/// Shots from a CSV string (see `lindhard exposure`) or an `(n, 3)` array of
/// `x_m, y_m, weight` point shots.
fn shots(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Shot>> {
    if let Ok(text) = obj.extract::<String>() {
        return parse_shots_csv(&text).map_err(errors::input);
    }
    let a: PyReadonlyArray2<f64> = obj.extract().map_err(|_| {
        errors::input("shots must be a CSV string or a float array of shape (n, 3)")
    })?;
    let a = a.as_array();
    if a.ncols() != 3 {
        return Err(errors::input(format!(
            "shots must have shape (n, 3): x_m, y_m, weight; got {} columns",
            a.ncols()
        )));
    }
    Ok(a.rows()
        .into_iter()
        .map(|r| Shot::point(r[0], r[1], r[2]))
        .collect())
}

fn with_threads<T: Send>(threads: Option<usize>, f: impl FnOnce() -> T + Send) -> PyResult<T> {
    match threads {
        None => Ok(f()),
        Some(0) => Err(errors::input("threads must be at least 1")),
        Some(n) => Ok(rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .map_err(|e| errors::run(&e))?
            .install(f)),
    }
}

/// Dose, eV/m², at the points `(x_m[i], y_m[i])` of the exposure of `shots`.
///
/// `shots` is a CSV string (columns `x_m,y_m,weight` and optionally
/// `blur_m,shape,size_m,size2_m`) or an `(n, 3)` array of point shots. The PSF
/// is the double Gaussian `(total_ev, alpha_m, beta_m, eta)`, or the triple
/// with `gamma_m` and `nu`. `blur_m` is the beam width `b` of `exp(-r²/b²)`.
/// Raises `InputError` for an invalid value.
#[pyfunction]
#[pyo3(signature = (shots_, x_m, y_m, *, total_ev, alpha_m, beta_m, eta, gamma_m=None, nu=None,
                    blur_m=0.0, subsamples=DEFAULT_SUBSAMPLES, threads=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn exposure_dose<'py>(
    py: Python<'py>,
    shots_: &Bound<'py, PyAny>,
    x_m: PyReadonlyArray1<'py, f64>,
    y_m: PyReadonlyArray1<'py, f64>,
    total_ev: f64,
    alpha_m: f64,
    beta_m: f64,
    eta: f64,
    gamma_m: Option<f64>,
    nu: Option<f64>,
    blur_m: f64,
    subsamples: usize,
    threads: Option<usize>,
) -> PyResult<Bound<'py, numpy::PyArray1<f64>>> {
    let (x, y) = (x_m.as_array(), y_m.as_array());
    if x.len() != y.len() {
        return Err(errors::input("x_m and y_m must have the same length"));
    }
    let points: Vec<[f64; 2]> = x.iter().zip(y.iter()).map(|(a, b)| [*a, *b]).collect();
    let e = Exposure::with_subsamples(
        psf(total_ev, alpha_m, beta_m, eta, gamma_m, nu)?,
        shots(shots_)?,
        blur_m,
        subsamples,
    )
    .map_err(errors::input)?;
    let dose = py.detach(|| with_threads(threads, || e.dose_at_points(&points)))?;
    Ok(dose.into_pyarray(py))
}

/// Dose, eV/m², on a regular grid: an array of shape `(ny, nx)` whose element
/// `[j, i]` is at `(x0_m + i dx_m, y0_m + j dy_m)`. Arguments as for
/// `exposure_dose`.
#[pyfunction]
#[pyo3(signature = (shots_, x0_m, y0_m, dx_m, dy_m, nx, ny, *, total_ev, alpha_m, beta_m, eta,
                    gamma_m=None, nu=None, blur_m=0.0, subsamples=DEFAULT_SUBSAMPLES,
                    threads=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn exposure_dose_grid<'py>(
    py: Python<'py>,
    shots_: &Bound<'py, PyAny>,
    x0_m: f64,
    y0_m: f64,
    dx_m: f64,
    dy_m: f64,
    nx: usize,
    ny: usize,
    total_ev: f64,
    alpha_m: f64,
    beta_m: f64,
    eta: f64,
    gamma_m: Option<f64>,
    nu: Option<f64>,
    blur_m: f64,
    subsamples: usize,
    threads: Option<usize>,
) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
    let e = Exposure::with_subsamples(
        psf(total_ev, alpha_m, beta_m, eta, gamma_m, nu)?,
        shots(shots_)?,
        blur_m,
        subsamples,
    )
    .map_err(errors::input)?;
    let grid = Grid {
        x0_m,
        y0_m,
        dx_m,
        dy_m,
        nx,
        ny,
    };
    grid.validate().map_err(errors::input)?;
    let dose = py
        .detach(|| with_threads(threads, || e.dose_on_grid(&grid)))?
        .map_err(errors::input)?;
    dose.dose_ev_per_m2.into_pyarray(py).reshape([ny, nx])
}

/// The dose far inside a large uniform exposure of `weight_per_m2` primaries
/// per m², as a dict `total_ev_per_m2`, `forward_ev_per_m2` (the `alpha` term)
/// and `long_range_ev_per_m2` (the backscatter terms).
#[pyfunction]
#[pyo3(signature = (*, total_ev, alpha_m, beta_m, eta, gamma_m=None, nu=None, weight_per_m2=1.0))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn area_dose<'py>(
    py: Python<'py>,
    total_ev: f64,
    alpha_m: f64,
    beta_m: f64,
    eta: f64,
    gamma_m: Option<f64>,
    nu: Option<f64>,
    weight_per_m2: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let a = uniform_area_dose(
        &psf(total_ev, alpha_m, beta_m, eta, gamma_m, nu)?,
        weight_per_m2,
    );
    let d = PyDict::new(py);
    d.set_item("total_ev_per_m2", a.total_ev_per_m2)?;
    d.set_item("forward_ev_per_m2", a.forward_ev_per_m2)?;
    d.set_item("long_range_ev_per_m2", a.long_range_ev_per_m2)?;
    Ok(d)
}

/// Iso-dose contours of `dose` (shape `(ny, nx)`, as `exposure_dose_grid`
/// returns) at `threshold`: a list of dicts with `points` (an `(n, 2)` array,
/// m), `closed`, `area_m2` (signed: positive for the outer boundary of a
/// feature, negative for a hole) and `centroid_m`. A contour that leaves the
/// grid is open and has `area_m2` and `centroid_m` of `None`.
#[pyfunction]
pub(crate) fn iso_contours<'py>(
    py: Python<'py>,
    dose: PyReadonlyArray2<'py, f64>,
    x0_m: f64,
    y0_m: f64,
    dx_m: f64,
    dy_m: f64,
    threshold: f64,
) -> PyResult<Bound<'py, PyList>> {
    let a = dose.as_array();
    let grid = Grid {
        x0_m,
        y0_m,
        dx_m,
        dy_m,
        nx: a.ncols(),
        ny: a.nrows(),
    };
    grid.validate().map_err(errors::input)?;
    let values: Vec<f64> = a.iter().copied().collect();
    let list = PyList::empty(py);
    for c in contours_of(&grid, &values, threshold) {
        list.append(contour_dict(py, &c)?)?;
    }
    Ok(list)
}

fn contour_dict<'py>(py: Python<'py>, c: &Contour) -> PyResult<Bound<'py, PyDict>> {
    let n = c.points.len();
    let flat: Vec<f64> = c.points.iter().flatten().copied().collect();
    let d = PyDict::new(py);
    d.set_item("points", flat.into_pyarray(py).reshape([n, 2])?)?;
    d.set_item("closed", c.closed)?;
    d.set_item("area_m2", c.signed_area_m2)?;
    d.set_item("centroid_m", c.centroid_m)?;
    Ok(d)
}
