//! [`RunResult`] and [`Histogram`]: the tallies of a finished run as NumPy
//! arrays, and the summary as a dict.
//!
//! Every number here is computed with the same expressions as the CLI's CSV
//! writers ([`lindhard_cli::output`]), so the arrays equal what the command
//! writes for the same input and seed, bit for bit.

use std::path::PathBuf;

use lindhard::input::Resolved;
use lindhard::tally::Histogram as LibHistogram;
use lindhard_cli::output::{self, density, symbol, NM};
use lindhard_cli::sim::Simulation;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::errors;

/// Counts over uniform bins, with the bin edges and per-ion densities.
///
/// `edges` has one more entry than `counts`. `underflow` and `overflow` count
/// samples outside the binned range (the stopped-primary `depth` histogram
/// has none: its last bin collects everything deeper and its upper edge is
/// `inf`). `density` is `counts / (histories * bin width)`: per incident ion
/// per unit of `unit` (`NaN` where the CSV leaves it empty). Lengths are in
/// nm, energies in eV and polar angles in degrees.
#[pyclass(module = "lindhard", frozen)]
pub struct Histogram {
    edges: Vec<f64>,
    counts: Vec<u64>,
    density: Vec<f64>,
    underflow: u64,
    overflow: u64,
    unit: &'static str,
}

#[pymethods]
impl Histogram {
    /// Bin edges (length `len(counts) + 1`).
    #[getter]
    fn edges<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.edges.clone().into_pyarray(py)
    }

    /// Bin centres (the last bin of an open-ended histogram uses its lower
    /// edge plus half the previous width).
    #[getter]
    fn centers<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let n = self.counts.len();
        let w = if n >= 2 {
            self.edges[1] - self.edges[0]
        } else {
            0.0
        };
        let c: Vec<f64> = (0..n)
            .map(|i| {
                let hi = self.edges[i + 1];
                if hi.is_finite() {
                    0.5 * (self.edges[i] + hi)
                } else {
                    self.edges[i] + 0.5 * w
                }
            })
            .collect();
        c.into_pyarray(py)
    }

    /// Count per bin (`uint64`).
    #[getter]
    fn counts<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        self.counts.clone().into_pyarray(py)
    }

    /// Count per incident ion per `unit`.
    #[getter]
    fn density<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.density.clone().into_pyarray(py)
    }

    #[getter]
    fn underflow(&self) -> u64 {
        self.underflow
    }

    #[getter]
    fn overflow(&self) -> u64 {
        self.overflow
    }

    /// Unit of the edges: `"nm"`, `"eV"` or `"deg"`.
    #[getter]
    fn unit(&self) -> &'static str {
        self.unit
    }

    fn __len__(&self) -> usize {
        self.counts.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Histogram({} bins in {}, underflow={}, overflow={})",
            self.counts.len(),
            self.unit,
            self.underflow,
            self.overflow
        )
    }
}

impl Histogram {
    /// A library histogram, edges and densities as the CSV writers do:
    /// edges `edge(i) / unit * scale`, widths `width / unit * scale`.
    fn from_lib(
        h: &LibHistogram,
        histories: u64,
        unit: f64,
        scale: f64,
        label: &'static str,
    ) -> Self {
        let b = &h.binning;
        let w = b.width() / unit * scale;
        Self {
            edges: (0..=b.bins).map(|i| b.edge(i) / unit * scale).collect(),
            counts: h.counts.clone(),
            density: h.counts.iter().map(|&c| density(c, histories, w)).collect(),
            underflow: h.underflow,
            overflow: h.overflow,
            unit: label,
        }
    }
}

/// The result of `Run.run()`.
///
/// Histograms are `Histogram` objects: `depth` (stopped primaries by depth),
/// `lateral_y`, `lateral_z` and `radial` (stopped primaries by lateral
/// position and radial distance), `vacancies`, `interstitials` and
/// `replacements` (cascade defects by depth) and `escapes` (energy and
/// polar-angle spectra of every escaping species). `summary()` is the
/// command's `summary.json` as a dict; `ions` holds the final state of every
/// primary as arrays (or `None` if `Tally.per_ion` is false).
///
/// The arrays are copies; the result is read-only.
#[pyclass(module = "lindhard", frozen)]
pub struct RunResult {
    resolved: Resolved,
    sim: Simulation,
}

impl RunResult {
    pub fn new(resolved: Resolved, sim: Simulation) -> Self {
        Self { resolved, sim }
    }

    fn lib_hist(&self, h: &LibHistogram, py: Python<'_>) -> PyResult<Py<Histogram>> {
        Py::new(
            py,
            Histogram::from_lib(h, self.sim.report.histories, NM, 1.0, "nm"),
        )
    }
}

#[pymethods]
impl RunResult {
    /// Primary histories run.
    #[getter]
    fn histories(&self) -> u64 {
        self.sim.tally.summary.histories
    }

    /// Stopped-primary depth histogram, nm (the CLI's `depth_profile.csv`).
    #[getter]
    fn depth(&self) -> Histogram {
        let s = &self.sim.tally.summary;
        let w = s.bin_width_m / NM;
        let n = s.histories as f64;
        let last = s.depth_hist.len() - 1;
        let edges = (0..=s.depth_hist.len())
            .map(|i| {
                if i > last {
                    f64::INFINITY
                } else {
                    i as f64 * w
                }
            })
            .collect();
        let density = s
            .depth_hist
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                if i == last {
                    f64::NAN
                } else {
                    c as f64 / (n * w)
                }
            })
            .collect();
        Histogram {
            edges,
            counts: s.depth_hist.clone(),
            density,
            underflow: 0,
            overflow: 0,
            unit: "nm",
        }
    }

    /// Stopped primaries by lateral position `y`, nm.
    #[getter]
    fn lateral_y(&self, py: Python<'_>) -> PyResult<Py<Histogram>> {
        self.lib_hist(&self.sim.report.range.lateral_y_histogram, py)
    }

    /// Stopped primaries by lateral position `z`, nm.
    #[getter]
    fn lateral_z(&self, py: Python<'_>) -> PyResult<Py<Histogram>> {
        self.lib_hist(&self.sim.report.range.lateral_z_histogram, py)
    }

    /// Stopped primaries by radial distance from the beam axis, nm.
    #[getter]
    fn radial(&self, py: Python<'_>) -> PyResult<Py<Histogram>> {
        self.lib_hist(&self.sim.report.range.radial_histogram, py)
    }

    /// Vacancies (displacements minus replacements) by depth, nm.
    #[getter]
    fn vacancies(&self, py: Python<'_>) -> PyResult<Py<Histogram>> {
        self.lib_hist(&self.sim.report.damage.vacancy_histogram, py)
    }

    /// Interstitials by depth, nm.
    #[getter]
    fn interstitials(&self, py: Python<'_>) -> PyResult<Py<Histogram>> {
        self.lib_hist(&self.sim.report.damage.interstitial_histogram, py)
    }

    /// Replacements by depth, nm.
    #[getter]
    fn replacements(&self, py: Python<'_>) -> PyResult<Py<Histogram>> {
        self.lib_hist(&self.sim.report.damage.replacement_histogram, py)
    }

    /// Escape spectra: a list of dicts with keys `z`, `symbol`, `beam`
    /// (whether the species is the beam species), `face` (`"front"` or
    /// `"back"`), `energy` (eV) and `polar` (degrees from the outward
    /// normal), the last two `Histogram`s.
    #[getter]
    fn escapes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let out = PyList::empty(py);
        let n = self.sim.report.histories;
        for s in &self.sim.report.escapes.species {
            for (face, f) in [("front", &s.front), ("back", &s.back)] {
                let d = PyDict::new(py);
                d.set_item("z", s.z)?;
                d.set_item("symbol", symbol(s.z))?;
                d.set_item("beam", s.beam)?;
                d.set_item("face", face)?;
                d.set_item(
                    "energy",
                    Histogram::from_lib(&f.energy_histogram, n, 1.0, 1.0, "eV"),
                )?;
                d.set_item(
                    "polar",
                    Histogram::from_lib(
                        &f.polar_histogram,
                        n,
                        1.0,
                        180.0 / std::f64::consts::PI,
                        "deg",
                    ),
                )?;
                out.append(d)?;
            }
        }
        Ok(out)
    }

    /// Final state of every primary, by history index, as a dict of arrays:
    /// `index` (`uint64`), `fate` (`uint8` codes into `lindhard.FATE_NAMES`),
    /// `position_nm` and `direction` (`(n, 3)`), `energy_ev` and `layer`
    /// (`uint64`). `None` if the tally was configured with `per_ion=False`.
    #[getter]
    fn ions<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let t = &self.sim.tally;
        if !t.per_ion {
            return Ok(None);
        }
        let n = t.finals.len();
        let flat = |f: &dyn Fn(&lindhard_cli::tally::FinalState) -> [f64; 3]| {
            let v: Vec<f64> = t.finals.iter().flat_map(f).collect();
            v.into_pyarray(py).reshape([n, 3])
        };
        let d = PyDict::new(py);
        d.set_item(
            "index",
            t.finals
                .iter()
                .map(|f| f.index)
                .collect::<Vec<_>>()
                .into_pyarray(py),
        )?;
        d.set_item(
            "fate",
            t.finals
                .iter()
                .map(|f| f.fate.code())
                .collect::<Vec<_>>()
                .into_pyarray(py),
        )?;
        d.set_item("position_nm", flat(&|f| f.pos.map(|x| x / NM))?)?;
        d.set_item(
            "energy_ev",
            t.finals
                .iter()
                .map(|f| f.energy_ev)
                .collect::<Vec<_>>()
                .into_pyarray(py),
        )?;
        d.set_item("direction", flat(&|f| f.dir)?)?;
        d.set_item(
            "layer",
            t.finals
                .iter()
                .map(|f| f.layer as u64)
                .collect::<Vec<_>>()
                .into_pyarray(py),
        )?;
        Ok(Some(d))
    }

    /// The command's `summary.json` as a dict: counts, yields, energy budget,
    /// range moments and Pearson fits, damage, sputtering, escapes and the
    /// models in use. Everything except the `run` block (thread count and
    /// timings) is a pure function of the input.
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let text = self.summary_json()?;
        py.import("json")?.call_method1("loads", (text,))
    }

    /// The command's `summary.json` text, exactly as it writes it.
    fn summary_json(&self) -> PyResult<String> {
        let s = &self.sim;
        output::summary_json(
            &self.resolved,
            &s.table,
            &s.tally,
            &s.report,
            &s.crystals,
            s.info,
        )
        .map_err(|e| errors::run(&e))
    }

    /// Write `summary.json` and the CSV profiles into `out_dir` (created if
    /// missing), the same files and bytes the `lindhard run` command writes
    /// (apart from the timing block of `summary.json`).
    fn write(&self, out_dir: PathBuf) -> PyResult<()> {
        let s = &self.sim;
        std::fs::create_dir_all(&out_dir).map_err(|e| {
            pyo3::exceptions::PyOSError::new_err(format!("{}: {e}", out_dir.display()))
        })?;
        let put = |name: &str, text: String| {
            let p = out_dir.join(name);
            std::fs::write(&p, text)
                .map_err(|e| pyo3::exceptions::PyOSError::new_err(format!("{}: {e}", p.display())))
        };
        put(output::SUMMARY_FILE, self.summary_json()?)?;
        put(output::DEPTH_FILE, output::depth_csv(&s.tally))?;
        put(output::LATERAL_FILE, output::lateral_csv(&s.report))?;
        put(output::DAMAGE_FILE, output::damage_csv(&s.report))?;
        put(output::ESCAPES_FILE, output::escapes_csv(&s.report))?;
        if s.tally.per_ion {
            put(output::IONS_FILE, output::ions_csv(&s.tally))?;
        }
        Ok(())
    }

    fn __repr__(&self) -> String {
        let s = &self.sim.tally.summary;
        format!(
            "RunResult({} ions: {} stopped, {} backscattered, {} transmitted, {} sputtered)",
            s.histories, s.primaries_stopped, s.backscattered, s.transmitted, s.sputtered
        )
    }
}

// The 2-D reshape above needs this trait in scope on some numpy versions.
#[allow(unused_imports)]
use PyArray2 as _;
