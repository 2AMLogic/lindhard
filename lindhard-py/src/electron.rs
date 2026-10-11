//! [`ElectronRun`] and [`ElectronRunResult`]: an electron run (an input with
//! an `[electron]` table) and its tallies.
//!
//! The configuration is the library's own [`ElectronInput`] and the run is the
//! CLI's driver ([`lindhard_cli::electron::simulate_electron`]); the summary
//! is the CLI's `electron_summary.json` and the arrays are computed with the
//! expressions of the CLI's CSV writers
//! ([`lindhard_cli::electron::spectra_csv`] and friends), so a result equals
//! the command's for the same input and seed. Nothing about the physics, the
//! schema or the table construction is duplicated here.

use std::path::PathBuf;

use lindhard::input::electron::{ElectronInput, ResolvedElectron};
use lindhard::tally::{ElectronReport, FaceEmission};
use lindhard_cli::electron::{self, ElectronSimulation};
use lindhard_cli::output::NM;
use numpy::{IntoPyArray, PyArrayMethods};
use pyo3::exceptions::PyUserWarning;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::errors;
use crate::result::Histogram;
use crate::run::io_err;

/// An electron run: the whole input document of an electron `lindhard run`
/// (the `[electron]` table, the target and `[run]`), read from the CLI's TOML
/// schema.
///
/// The first increment is TOML-first: `ElectronRun.from_toml(text)` and
/// `ElectronRun.from_toml_file(path)` parse the document, `to_toml()` writes
/// it back, and `histories`, `seed` and `threads` can be changed. Relative
/// paths to optical-ELF and other data files resolve against `base_dir`,
/// which `from_toml_file` sets to the file's directory (the working directory
/// otherwise), as the command resolves them.
#[pyclass(module = "lindhard")]
#[derive(Debug)]
pub struct ElectronRun {
    input: ElectronInput,
    #[pyo3(get, set)]
    base_dir: Option<PathBuf>,
}

#[pymethods]
impl ElectronRun {
    /// Parse a TOML document in the CLI's electron schema.
    ///
    /// Raises `InputError` naming the key for an unknown key, a missing key or
    /// a wrong type.
    #[staticmethod]
    #[pyo3(signature = (text, base_dir=None))]
    fn from_toml(text: &str, base_dir: Option<PathBuf>) -> PyResult<Self> {
        let input = ElectronInput::from_toml_str(text).map_err(errors::input)?;
        Ok(Self { input, base_dir })
    }

    /// Read and parse a TOML file; relative data-file paths resolve against
    /// its directory. Raises `FileNotFoundError` or `OSError` naming the path
    /// if it cannot be read.
    #[staticmethod]
    fn from_toml_file(path: PathBuf) -> PyResult<Self> {
        let text = std::fs::read_to_string(&path).map_err(|e| io_err(&path, e))?;
        let input = ElectronInput::from_toml_str(&text)
            .map_err(|e| errors::input(format!("{}: invalid input: {e}", path.display())))?;
        let base_dir = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .map(PathBuf::from);
        Ok(Self { input, base_dir })
    }

    /// The run as a TOML document in the CLI's schema. Parsing it with
    /// `ElectronRun.from_toml` gives back an equal run.
    fn to_toml(&self) -> PyResult<String> {
        self.input.to_toml_string().map_err(errors::input)
    }

    /// Check the run, including reading its data files, without building
    /// tables or transporting. Returns the list of warnings (the CLI prints
    /// them on stderr); raises `InputError` naming the field or file if the
    /// run is invalid.
    fn validate(&self) -> PyResult<Vec<String>> {
        Ok(self.resolve(&self.input)?.warnings)
    }

    /// Build the cross-section tables, run the histories and return an
    /// `ElectronRunResult`.
    ///
    /// `histories`, `seed` and `threads` override the run's own values for
    /// this call only. The GIL is released while the tables are built and
    /// transport runs. Results depend only on the input and the seed, never on
    /// `threads`. No table cache is used.
    ///
    /// Raises `InputError` for an invalid configuration (resolution warnings
    /// are issued as `UserWarning`) and `RunError` if table construction or
    /// transport fails.
    #[pyo3(signature = (*, histories=None, seed=None, threads=None))]
    fn run(
        &self,
        py: Python<'_>,
        histories: Option<u64>,
        seed: Option<u64>,
        threads: Option<usize>,
    ) -> PyResult<ElectronRunResult> {
        let mut input = self.input.clone();
        if let Some(n) = histories {
            input.run.histories = n;
        }
        if let Some(s) = seed {
            input.run.seed = s;
        }
        if threads.is_some() {
            input.run.threads = threads;
        }
        let resolved = self.resolve(&input)?;
        for w in &resolved.warnings {
            PyErr::warn(
                py,
                &py.get_type::<PyUserWarning>(),
                &std::ffi::CString::new(w.as_str()).unwrap_or_default(),
                1,
            )?;
        }
        let sim = py
            .detach(|| electron::simulate_electron(&resolved, input.run.threads, None))
            .map_err(|e| errors::run(&e))?;
        Ok(ElectronRunResult { resolved, sim })
    }

    /// Number of primary electron histories.
    #[getter]
    fn histories(&self) -> u64 {
        self.input.run.histories
    }

    #[setter]
    fn set_histories(&mut self, histories: u64) {
        self.input.run.histories = histories;
    }

    /// Run seed.
    #[getter]
    fn seed(&self) -> u64 {
        self.input.run.seed
    }

    #[setter]
    fn set_seed(&mut self, seed: u64) {
        self.input.run.seed = seed;
    }

    /// Worker threads (`None`: all available cores). Never changes results.
    #[getter]
    fn threads(&self) -> Option<usize> {
        self.input.run.threads
    }

    #[setter]
    fn set_threads(&mut self, threads: Option<usize>) {
        self.input.run.threads = threads;
    }

    fn __repr__(&self) -> String {
        format!(
            "ElectronRun(energy_ev={}, histories={}, seed={})",
            self.input.electron.beam.energy_ev, self.input.run.histories, self.input.run.seed
        )
    }
}

impl ElectronRun {
    fn resolve(&self, input: &ElectronInput) -> PyResult<ResolvedElectron> {
        let base = self
            .base_dir
            .as_deref()
            .unwrap_or_else(|| std::path::Path::new("."));
        input.resolve_in(base).map_err(errors::input)
    }
}

/// The result of `ElectronRun.run()`.
///
/// * `summary()`: the command's `electron_summary.json` as a dict (input
///   echo, models and provenance, transport metadata, table diagnostics,
///   energy budget, yields, per-layer deposition, PSF totals and fits, and a
///   `run` block of thread count and timings, the only part that depends on
///   `threads`).
/// * `escape_spectra`: per face, `Histogram`s of the escaping electrons.
/// * `cartesian`, `cylindrical`: the configured energy-deposition grids as
///   dicts of NumPy arrays, or `None` if the input does not configure them.
///
/// The arrays are copies; the result is read-only.
#[pyclass(module = "lindhard", frozen)]
#[derive(Debug)]
pub struct ElectronRunResult {
    resolved: ResolvedElectron,
    sim: ElectronSimulation,
}

fn face_dict<'py>(
    py: Python<'py>,
    f: &FaceEmission,
    histories: u64,
) -> PyResult<Bound<'py, PyDict>> {
    let deg = 180.0 / std::f64::consts::PI;
    let d = PyDict::new(py);
    d.set_item(
        "energy",
        Histogram::from_lib(&f.energy_histogram, histories, 1.0, 1.0, "eV"),
    )?;
    d.set_item(
        "polar_slow",
        Histogram::from_lib(&f.slow.polar_histogram, histories, 1.0, deg, "deg"),
    )?;
    d.set_item(
        "polar_fast",
        Histogram::from_lib(&f.fast.polar_histogram, histories, 1.0, deg, "deg"),
    )?;
    Ok(d)
}

fn edges_nm(b: &lindhard::tally::Binning) -> Vec<f64> {
    (0..=b.bins).map(|i| b.edge(i) / NM).collect()
}

fn report(r: &ElectronSimulation) -> &ElectronReport {
    &r.report
}

#[pymethods]
impl ElectronRunResult {
    /// Primary histories run.
    #[getter]
    fn histories(&self) -> u64 {
        report(&self.sim).histories
    }

    /// Escape spectra of the escaping electrons, as the CLI's
    /// `electron_escape_spectra.csv` writes them: a dict with keys `"front"`
    /// and `"back"` (the back face of a finite stack; all counts are zero for a substrate),
    /// each a dict of three `Histogram`s: `energy` (all escaping electrons,
    /// eV) and `polar_slow` / `polar_fast` (polar angle in degrees from the
    /// outward normal of the slow, secondary-class and fast,
    /// backscattered-class electrons). `density` is per primary per eV or per
    /// degree.
    #[getter]
    fn escape_spectra<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let r = report(&self.sim);
        let d = PyDict::new(py);
        d.set_item("front", face_dict(py, &r.front, r.histories)?)?;
        d.set_item("back", face_dict(py, &r.back, r.histories)?)?;
        Ok(d)
    }

    /// The Cartesian energy-deposition grid (`[electron.tally.cartesian]`) or
    /// `None` if it is not configured. A dict with `x_edges_nm` (depth),
    /// `y_edges_nm`, `z_edges_nm` (lengths `n + 1`), `energy_ev` of shape
    /// `(nx, ny, nz)` (energy deposited in each voxel, summed over all
    /// histories, eV; index order `(ix, iy, iz)`, `z` fastest),
    /// `ev_per_primary_per_nm3` of the same shape (`energy_ev` divided by the
    /// histories and the voxel volume) and `outside_ev` (energy deposited
    /// outside the grid).
    #[getter]
    fn cartesian<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let r = report(&self.sim);
        let Some(c) = r.deposition.cartesian.as_ref() else {
            return Ok(None);
        };
        let g = &c.grid;
        let shape = [g.x.bins, g.y.bins, g.z.bins];
        let n = r.histories.max(1) as f64;
        let v = g.voxel_volume_m3() / (NM * NM * NM);
        let dens: Vec<f64> = c.energy_ev.iter().map(|e| e / (n * v)).collect();
        let d = PyDict::new(py);
        d.set_item("x_edges_nm", edges_nm(&g.x).into_pyarray(py))?;
        d.set_item("y_edges_nm", edges_nm(&g.y).into_pyarray(py))?;
        d.set_item("z_edges_nm", edges_nm(&g.z).into_pyarray(py))?;
        d.set_item(
            "energy_ev",
            c.energy_ev.clone().into_pyarray(py).reshape(shape)?,
        )?;
        d.set_item(
            "ev_per_primary_per_nm3",
            dens.into_pyarray(py).reshape(shape)?,
        )?;
        d.set_item("outside_ev", c.outside_ev)?;
        Ok(Some(d))
    }

    /// The cylindrical energy-deposition grid (`[electron.tally.cylindrical]`,
    /// about the beam axis) or `None` if it is not configured. A dict with
    /// `r_edges_nm` and `depth_edges_nm` (lengths `n + 1`), `energy_ev` of
    /// shape `(nr, ndepth)` (energy deposited in each cell, summed over all
    /// histories, eV; index order `(ir, ix)`, depth fastest),
    /// `ev_per_primary_per_nm3` of the same shape (`energy_ev` divided by the
    /// histories and the annulus volume) and `outside_ev`.
    #[getter]
    fn cylindrical<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let r = report(&self.sim);
        let Some(c) = r.deposition.cylindrical.as_ref() else {
            return Ok(None);
        };
        let g = &c.grid;
        let shape = [g.r.bins, g.depth.bins];
        let n = r.histories.max(1) as f64;
        let mut dens = Vec::with_capacity(c.energy_ev.len());
        for ir in 0..g.r.bins {
            let v = g.cell_volume_m3(ir) / (NM * NM * NM);
            for ix in 0..g.depth.bins {
                dens.push(c.at(ir, ix) / (n * v));
            }
        }
        let d = PyDict::new(py);
        d.set_item("r_edges_nm", edges_nm(&g.r).into_pyarray(py))?;
        d.set_item("depth_edges_nm", edges_nm(&g.depth).into_pyarray(py))?;
        d.set_item(
            "energy_ev",
            c.energy_ev.clone().into_pyarray(py).reshape(shape)?,
        )?;
        d.set_item(
            "ev_per_primary_per_nm3",
            dens.into_pyarray(py).reshape(shape)?,
        )?;
        d.set_item("outside_ev", c.outside_ev)?;
        Ok(Some(d))
    }

    /// The command's `electron_summary.json` as a dict. Everything except the
    /// `run` block (thread count and timings) and the `origin` of the tables
    /// is a pure function of the input and the seed.
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let text = self.summary_json()?;
        py.import("json")?.call_method1("loads", (text,))
    }

    /// The command's `electron_summary.json` text, exactly as it writes it.
    fn summary_json(&self) -> PyResult<String> {
        electron::summary_json(&self.resolved, &self.sim).map_err(|e| errors::run(&e))
    }

    fn __repr__(&self) -> String {
        let r = report(&self.sim);
        format!(
            "ElectronRunResult({} histories: eta {:.4}, delta {:.4})",
            r.histories, r.yields.backscatter_eta, r.yields.secondary_delta
        )
    }
}
