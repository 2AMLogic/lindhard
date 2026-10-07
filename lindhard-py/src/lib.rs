//! Python bindings for `lindhard`.
//!
//! The classes here are thin value wrappers over the library's own input
//! schema ([`lindhard::input`]); a run is the CLI's run
//! ([`lindhard_cli::sim::simulate`]) and the summary and CSV text are the
//! CLI's ([`lindhard_cli::output`]). Nothing about the physics or the TOML
//! schema is duplicated here, so a configuration can go to the command line
//! and back, and a result matches the command's for the same input and seed.

mod errors;
mod result;
mod run;
mod spec;

use pyo3::prelude::*;

/// The compiled module `lindhard._lindhard`; the `lindhard` package
/// re-exports its contents.
#[pymodule]
fn _lindhard(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", lindhard::VERSION)?;

    errors::register(m)?;
    m.add_class::<spec::Element>()?;
    m.add_class::<spec::Material>()?;
    m.add_class::<spec::Layer>()?;
    m.add_class::<spec::Target>()?;
    m.add_class::<spec::Beam>()?;
    m.add_class::<spec::Physics>()?;
    m.add_class::<spec::Tally>()?;
    m.add_class::<run::Run>()?;
    m.add_class::<result::Histogram>()?;
    m.add_class::<result::RunResult>()?;
    m.add(
        "FATE_NAMES",
        lindhard_cli::tally::Fate::ALL.map(|f| f.as_str()),
    )?;
    Ok(())
}
