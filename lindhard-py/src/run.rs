//! [`Run`]: a complete run description (the whole input document) and the
//! call that executes it.

use std::collections::BTreeMap;
use std::path::PathBuf;

use lindhard::input::{Input, Resolved, RunSpec, StoppingSpec};
use pyo3::exceptions::PyUserWarning;
use pyo3::prelude::*;

use crate::errors;
use crate::result::RunResult;
use crate::spec::{Beam, Material, Physics, Tally, Target};

/// A complete run description: beam, target, physics, tally options, run size
/// and seed; the whole input document of the `lindhard` command.
///
/// `materials` names inline materials that layers may refer to by name;
/// `stopping_tables` lists user stopping-table files (relative paths resolve
/// against `base_dir`, which `Run.from_toml_file` sets to the file's
/// directory). `beam`, `target`, `physics` and `tally` return copies: assign a
/// modified copy back to change the run.
///
/// `Run.from_toml(text)` and `run.to_toml()` use the same schema as the
/// command line, so a configuration can go either way.
#[pyclass(module = "lindhard")]
pub struct Run {
    input: Input,
    #[pyo3(get, set)]
    base_dir: Option<PathBuf>,
}

fn io_err(path: &std::path::Path, e: std::io::Error) -> PyErr {
    let msg = format!("{}: {e}", path.display());
    if e.kind() == std::io::ErrorKind::NotFound {
        pyo3::exceptions::PyFileNotFoundError::new_err(msg)
    } else {
        pyo3::exceptions::PyOSError::new_err(msg)
    }
}

#[pymethods]
impl Run {
    #[new]
    #[pyo3(signature = (
        beam, target, physics, ions, seed, *, threads=None, tally=None,
        materials=None, stopping_tables=None, base_dir=None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        beam: PyRef<'_, Beam>,
        target: PyRef<'_, Target>,
        physics: PyRef<'_, Physics>,
        ions: u64,
        seed: u64,
        threads: Option<usize>,
        tally: Option<PyRef<'_, Tally>>,
        materials: Option<BTreeMap<String, PyRef<'_, Material>>>,
        stopping_tables: Option<Vec<String>>,
        base_dir: Option<PathBuf>,
    ) -> PyResult<Self> {
        let input = Input {
            beam: beam.to_spec(),
            materials: materials
                .unwrap_or_default()
                .into_iter()
                .map(|(k, m)| (k, m.to_spec()))
                .collect(),
            target: target.to_spec(),
            physics: physics.to_spec()?,
            stopping: StoppingSpec {
                tables: stopping_tables.unwrap_or_default(),
            },
            run: RunSpec {
                ions,
                seed,
                threads,
            },
            tally: tally.map(|t| t.to_spec()).unwrap_or_default(),
            dynamic: None,
            crystal: Vec::new(),
        };
        Ok(Self { input, base_dir })
    }

    /// Parse a TOML document in the CLI's schema.
    ///
    /// Raises `InputError` naming the key and line for an unknown key, a
    /// missing key or a wrong type.
    #[staticmethod]
    #[pyo3(signature = (text, base_dir=None))]
    fn from_toml(text: &str, base_dir: Option<PathBuf>) -> PyResult<Self> {
        let input = Input::from_toml_str(text).map_err(errors::input)?;
        Ok(Self { input, base_dir })
    }

    /// Read and parse a TOML file; relative stopping-table paths resolve
    /// against its directory.
    #[staticmethod]
    fn from_toml_file(path: PathBuf) -> PyResult<Self> {
        let text = std::fs::read_to_string(&path).map_err(|e| io_err(&path, e))?;
        let input = Input::from_toml_str(&text)
            .map_err(|e| errors::input(format!("{}: invalid input: {e}", path.display())))?;
        let base_dir = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .map(PathBuf::from);
        Ok(Self { input, base_dir })
    }

    /// The run as a TOML document in the CLI's schema. Parsing it with
    /// `Run.from_toml` gives back an equal run.
    fn to_toml(&self) -> PyResult<String> {
        self.input.to_toml_string().map_err(errors::input)
    }

    /// Check the run without executing it. Returns the list of warnings (the
    /// CLI prints them on stderr); raises `InputError` if the run is invalid.
    fn validate(&self) -> PyResult<Vec<String>> {
        Ok(self.resolve(&self.input)?.warnings)
    }

    /// Run the simulation and return a `RunResult`.
    ///
    /// `ions`, `seed` and `threads` override the run's own values for this
    /// call only. The GIL is released while transport runs. Results depend
    /// only on the input and the seed, never on `threads`.
    ///
    /// Raises `InputError` for an invalid configuration and `RunError` if
    /// transport fails.
    #[pyo3(signature = (*, ions=None, seed=None, threads=None))]
    fn run(
        &self,
        py: Python<'_>,
        ions: Option<u64>,
        seed: Option<u64>,
        threads: Option<usize>,
    ) -> PyResult<RunResult> {
        let mut input = self.input.clone();
        if let Some(n) = ions {
            input.run.ions = n;
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
        if input.run.threads == Some(0) {
            return Err(errors::input("threads must be at least 1"));
        }
        let sim = py
            .detach(|| lindhard_cli::sim::simulate(&resolved, input.run.threads))
            .map_err(|e| errors::run(&e))?;
        Ok(RunResult::new(resolved, sim))
    }

    #[getter]
    fn beam(&self) -> Beam {
        Beam::from_spec(&self.input.beam)
    }

    #[setter]
    fn set_beam(&mut self, beam: PyRef<'_, Beam>) {
        self.input.beam = beam.to_spec();
    }

    #[getter]
    fn target(&self) -> Target {
        Target::from_spec(&self.input.target)
    }

    #[setter]
    fn set_target(&mut self, target: PyRef<'_, Target>) {
        self.input.target = target.to_spec();
    }

    #[getter]
    fn physics(&self) -> Physics {
        Physics::from_spec(&self.input.physics)
    }

    #[setter]
    fn set_physics(&mut self, physics: PyRef<'_, Physics>) -> PyResult<()> {
        self.input.physics = physics.to_spec()?;
        Ok(())
    }

    #[getter]
    fn tally(&self) -> Tally {
        Tally::from_spec(&self.input.tally)
    }

    #[setter]
    fn set_tally(&mut self, tally: PyRef<'_, Tally>) {
        self.input.tally = tally.to_spec();
    }

    /// Named inline materials (a copy).
    #[getter]
    fn materials(&self) -> BTreeMap<String, Material> {
        self.input
            .materials
            .iter()
            .map(|(k, m)| (k.clone(), Material::from_spec(m)))
            .collect()
    }

    #[setter]
    fn set_materials(&mut self, materials: BTreeMap<String, PyRef<'_, Material>>) {
        self.input.materials = materials
            .into_iter()
            .map(|(k, m)| (k, m.to_spec()))
            .collect();
    }

    #[getter]
    fn stopping_tables(&self) -> Vec<String> {
        self.input.stopping.tables.clone()
    }

    #[setter]
    fn set_stopping_tables(&mut self, tables: Vec<String>) {
        self.input.stopping.tables = tables;
    }

    /// Number of primary ions (histories).
    #[getter]
    fn ions(&self) -> u64 {
        self.input.run.ions
    }

    #[setter]
    fn set_ions(&mut self, ions: u64) {
        self.input.run.ions = ions;
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
            "Run(ion={}, energy_ev={}, ions={}, seed={})",
            self.input.beam.ion,
            self.input.beam.energy_ev,
            self.input.run.ions,
            self.input.run.seed
        )
    }
}

impl Run {
    fn resolve(&self, input: &Input) -> PyResult<Resolved> {
        let base = self
            .base_dir
            .as_deref()
            .unwrap_or_else(|| std::path::Path::new("."));
        input.resolve_in(base).map_err(errors::input)
    }
}
