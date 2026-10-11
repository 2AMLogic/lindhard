//! [`DynamicRunResult`]: the fluence-step time series and slab composition
//! snapshots of a finished dynamic run, as NumPy arrays, and the summary as a
//! dict.
//!
//! Every number is computed with the same expressions as the CLI's dynamic
//! CSV writers ([`lindhard_cli::output::dynamic_steps_csv`] and
//! [`lindhard_cli::output::dynamic_composition_csv`]), so the arrays equal
//! what the command writes for the same input and seed, bit for bit. Units
//! are the CSVs' units: lengths in nm, fluence and areal inventory per cm^2.

use std::path::PathBuf;

use lindhard::input::Resolved;
use lindhard_cli::dynamic::DynamicSimulation;
use lindhard_cli::output::{self, symbol, NM};
use numpy::{IntoPyArray, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::errors;

/// The result of `Run.run_dynamic()`.
///
/// `species` lists the symbols of every element the run can contain (the beam
/// species and every target element, sorted by atomic number; `species_z`
/// holds the atomic numbers). `steps` is a dict of per-step arrays whose
/// index 0 is the initial target (step 0) followed by every accepted step;
/// `composition` is a list with one dict per entry of `steps`, holding the
/// slab profile after that step. The number of slabs can change between
/// steps (erosion removes emptied slabs), so each snapshot has its own
/// length: the snapshots are not stacked into one rectangular array.
///
/// Units match the dynamic CSVs: lengths in nm, `fluence_cm2` in ions/cm^2,
/// areal inventory in atoms/cm^2. `summary()` is `dynamic_summary.json` as a
/// dict. The arrays are copies; the result is read-only.
#[pyclass(module = "lindhard", frozen)]
#[derive(Debug)]
pub struct DynamicRunResult {
    resolved: Resolved,
    sim: DynamicSimulation,
}

impl DynamicRunResult {
    pub fn new(resolved: Resolved, sim: DynamicSimulation) -> Self {
        Self { resolved, sim }
    }
}

#[pymethods]
impl DynamicRunResult {
    /// Symbols of the species the run can contain, sorted by atomic number.
    #[getter]
    fn species(&self) -> Vec<&'static str> {
        self.sim.species.iter().map(|&z| symbol(z)).collect()
    }

    /// Atomic numbers of `species` (`uint8`).
    #[getter]
    fn species_z<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<u8>> {
        self.sim.species.clone().into_pyarray(py)
    }

    /// Number of entries of `steps` and `composition`: the accepted steps
    /// plus step 0.
    #[getter]
    fn n_steps(&self) -> usize {
        self.sim.steps.len()
    }

    /// Attempts that were rejected and retried by the adaptive step bound,
    /// over the whole run.
    #[getter]
    fn rejected_attempts(&self) -> u64 {
        self.sim.rejected_attempts
    }

    /// Per-step time series as a dict of arrays, one entry per step with
    /// step 0 (the initial target) first.
    ///
    /// Keys: `step`, `first_index`, `ions`, `ions_done`, `attempts`,
    /// `clamped`, `removed_slabs`, `n_slabs` (`uint64`; `attempts`, `clamped`
    /// are small integers), `fluence_cm2` (ions/cm^2 delivered so far),
    /// `max_change`, `surface_nm` (cumulative surface recession, 0 with
    /// erosion off; depths of a step plus this are depths in the original
    /// frame), `thickness_nm` (summed slab thickness),
    /// `cum_backscattered`, `cum_transmitted`, `cum_stopped_in_target`,
    /// `cum_stopped_in_substrate`, `cum_sputtered`, `cum_recoils_transmitted`
    /// (`uint64`, summed over steps 1..=this one), `sputter_yield`
    /// (cumulative sputtered atoms per ion, `NaN` at step 0) and
    /// `cum_sputtered_by_species` (`uint64`, shape `(n_steps, len(species))`).
    #[getter]
    fn steps<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let st = &self.sim.steps;
        let d = PyDict::new(py);
        let u = |f: &dyn Fn(&lindhard_cli::dynamic::StepRow) -> u64| -> Vec<u64> {
            st.iter().map(f).collect()
        };
        let x = |f: &dyn Fn(&lindhard_cli::dynamic::StepRow) -> f64| -> Vec<f64> {
            st.iter().map(f).collect()
        };
        d.set_item("step", u(&|s| s.step).into_pyarray(py))?;
        d.set_item("first_index", u(&|s| s.first_index).into_pyarray(py))?;
        d.set_item("ions", u(&|s| s.ions).into_pyarray(py))?;
        d.set_item("ions_done", u(&|s| s.ions_done).into_pyarray(py))?;
        d.set_item("attempts", u(&|s| u64::from(s.attempts)).into_pyarray(py))?;
        d.set_item("clamped", u(&|s| u64::from(s.clamped)).into_pyarray(py))?;
        d.set_item(
            "removed_slabs",
            u(&|s| s.removed_slabs as u64).into_pyarray(py),
        )?;
        d.set_item("n_slabs", u(&|s| s.slabs.len() as u64).into_pyarray(py))?;
        d.set_item("fluence_cm2", x(&|s| s.fluence_m2 * 1e-4).into_pyarray(py))?;
        d.set_item("max_change", x(&|s| s.max_change).into_pyarray(py))?;
        d.set_item("surface_nm", x(&|s| s.recession_m / NM).into_pyarray(py))?;
        d.set_item(
            "thickness_nm",
            x(&|s| s.slabs.iter().map(|b| b.back_m - b.front_m).sum::<f64>() / NM).into_pyarray(py),
        )?;
        d.set_item(
            "cum_backscattered",
            u(&|s| s.cumulative.backscattered).into_pyarray(py),
        )?;
        d.set_item(
            "cum_transmitted",
            u(&|s| s.cumulative.transmitted).into_pyarray(py),
        )?;
        d.set_item(
            "cum_stopped_in_target",
            u(&|s| s.cumulative.primaries_in_slabs).into_pyarray(py),
        )?;
        d.set_item(
            "cum_stopped_in_substrate",
            u(&|s| s.cumulative.primaries_in_substrate).into_pyarray(py),
        )?;
        d.set_item(
            "cum_sputtered",
            u(&|s| s.cumulative.sputtered_total()).into_pyarray(py),
        )?;
        d.set_item(
            "sputter_yield",
            x(&|s| {
                let y = &s.cumulative;
                if y.histories > 0 {
                    y.sputtered_total() as f64 / y.histories as f64
                } else {
                    f64::NAN
                }
            })
            .into_pyarray(py),
        )?;
        d.set_item(
            "cum_recoils_transmitted",
            u(&|s| s.cumulative.recoils_transmitted).into_pyarray(py),
        )?;
        let by_species: Vec<u64> = st
            .iter()
            .flat_map(|s| {
                self.sim
                    .species
                    .iter()
                    .map(|z| s.cumulative.sputtered.get(z).copied().unwrap_or(0))
            })
            .collect();
        d.set_item(
            "cum_sputtered_by_species",
            by_species
                .into_pyarray(py)
                .reshape([st.len(), self.sim.species.len()])?,
        )?;
        Ok(d)
    }

    /// Slab profile after each step, a list parallel to the arrays of
    /// `steps` (index 0: the initial target). Each entry is a dict with
    /// `step` (int), `front_nm`, `back_nm` and `thickness_nm` (shape
    /// `(n_slabs,)`, depths from the front surface of that step), and
    /// `atoms_per_cm2` and `fraction` (shape `(n_slabs, len(species))`,
    /// columns in `species` order; `fraction` is the atom fraction of the
    /// slab, `NaN` for an empty slab). `n_slabs` varies between entries when
    /// erosion removes slabs.
    #[getter]
    fn composition<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let out = PyList::empty(py);
        let ns = self.sim.species.len();
        for s in &self.sim.steps {
            let n = s.slabs.len();
            let d = PyDict::new(py);
            d.set_item("step", s.step)?;
            d.set_item(
                "front_nm",
                s.slabs
                    .iter()
                    .map(|b| b.front_m / NM)
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
            d.set_item(
                "back_nm",
                s.slabs
                    .iter()
                    .map(|b| b.back_m / NM)
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
            d.set_item(
                "thickness_nm",
                s.slabs
                    .iter()
                    .map(|b| (b.back_m - b.front_m) / NM)
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
            let mut atoms = Vec::with_capacity(n * ns);
            let mut frac = Vec::with_capacity(n * ns);
            for b in &s.slabs {
                let get = |z: u8| {
                    b.inventory
                        .iter()
                        .find(|&&(zz, _)| zz == z)
                        .map_or(0.0, |&(_, a)| a)
                };
                let total: f64 = b.inventory.iter().map(|&(_, a)| a).sum();
                atoms.extend(self.sim.species.iter().map(|&z| get(z) * 1e-4));
                frac.extend(self.sim.species.iter().map(|&z| get(z) / total));
            }
            d.set_item("atoms_per_cm2", atoms.into_pyarray(py).reshape([n, ns])?)?;
            d.set_item("fraction", frac.into_pyarray(py).reshape([n, ns])?)?;
            out.append(d)?;
        }
        Ok(out)
    }

    /// The command's `dynamic_summary.json` as a dict: the echoed input, the
    /// models, run totals and the file list. Everything except the `run`
    /// block (thread count and timings) is a pure function of the input.
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let text = self.summary_json()?;
        py.import("json")?.call_method1("loads", (text,))
    }

    /// The command's `dynamic_summary.json` text, exactly as it writes it.
    fn summary_json(&self) -> PyResult<String> {
        output::dynamic_summary_json(&self.resolved, &self.sim).map_err(|e| errors::run(&e))
    }

    /// Write `dynamic_summary.json`, `dynamic_steps.csv` and
    /// `dynamic_composition.csv` into `out_dir` (created if missing), the same
    /// files and bytes the `lindhard run` command writes (apart from the
    /// timing block of the summary). Like the command, it removes the
    /// previous `dynamic_summary.json` before replacing any CSV and writes
    /// the new summary last. Raises `OSError` naming the path on failure.
    fn write(&self, out_dir: PathBuf) -> PyResult<()> {
        std::fs::create_dir_all(&out_dir).map_err(|e| {
            pyo3::exceptions::PyOSError::new_err(format!("{}: {e}", out_dir.display()))
        })?;
        let summary = self.summary_json()?;
        lindhard_cli::publish::publish_dynamic(&out_dir, summary, &self.sim)
            .map_err(|e| pyo3::exceptions::PyOSError::new_err(format!("{e:#}")))
    }

    fn __repr__(&self) -> String {
        let last = self.sim.steps.last().expect("step 0");
        format!(
            "DynamicRunResult({} ions in {} steps: {} sputtered, {} slabs left)",
            last.ions_done,
            last.step,
            last.cumulative.sputtered_total(),
            last.slabs.len()
        )
    }
}
