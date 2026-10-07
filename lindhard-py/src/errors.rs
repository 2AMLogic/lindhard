//! Python exception types and the mapping from Rust errors.
//!
//! * `InputError` (a `LindhardError` and a `ValueError`): the configuration is
//!   malformed or has an invalid value; the message names the offending field.
//! * `RunError` (a `LindhardError` and a `RuntimeError`): the configuration was
//!   accepted but the run failed, for example a stopping table queried outside
//!   its range.
//!
//! `except lindhard.LindhardError` catches everything this package raises on
//! purpose. The classes are created at import time with two bases each, which
//! `create_exception!` cannot express, and looked up when an error is raised.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple, PyType};

const MODULE: &str = "lindhard._lindhard";

/// Create the exception classes and add them to the module.
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let base = PyErr::new_type(
        py,
        c"lindhard.LindhardError",
        Some(c"Base class of every exception raised by lindhard."),
        Some(&py.get_type::<pyo3::exceptions::PyException>()),
        None,
    )?;
    let derived = |name: &str, doc: &str, other: Bound<'_, PyType>| -> PyResult<Py<PyAny>> {
        let bases = PyTuple::new(py, [base.bind(py).clone(), other])?;
        let ns = PyDict::new(py);
        ns.set_item("__doc__", doc)?;
        ns.set_item("__module__", "lindhard")?;
        Ok(py.get_type::<PyType>().call1((name, bases, ns))?.unbind())
    };
    let input = derived(
        "InputError",
        "The configuration is malformed or has an invalid value.",
        py.get_type::<PyValueError>(),
    )?;
    let run = derived(
        "RunError",
        "The configuration was accepted but the run failed.",
        py.get_type::<PyRuntimeError>(),
    )?;
    m.add("LindhardError", base)?;
    m.add("InputError", input)?;
    m.add("RunError", run)?;
    Ok(())
}

fn raise(name: &str, msg: String) -> PyErr {
    Python::attach(|py| {
        let cls = py
            .import(MODULE)
            .and_then(|m| m.getattr(name))
            .and_then(|c| c.cast_into::<PyType>().map_err(PyErr::from));
        match cls {
            Ok(cls) => PyErr::from_type(cls, msg),
            Err(e) => e,
        }
    })
}

/// An invalid configuration.
pub fn input(msg: impl std::fmt::Display) -> PyErr {
    raise("InputError", msg.to_string())
}

/// A failed run, with the full context chain (outermost first) as the CLI
/// prints it.
pub fn run(err: &impl std::fmt::Display) -> PyErr {
    raise("RunError", format!("{err:#}"))
}
