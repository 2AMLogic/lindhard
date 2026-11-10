//! The damage-export file as NumPy arrays: [`RunResult::damage_export`]
//! (ion runs) and [`read_damage_export`] (any engine's file). Both go through
//! the library's [`VoxelReport`], so the arrays of a written file equal the
//! arrays of the run that wrote it, bit for bit.

use std::path::PathBuf;

use lindhard::tally::{ExportProvenance, VoxelReport};
use numpy::{IntoPyArray, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::errors;

/// Dense `(nx, ny, nz)` array of one per-voxel value; unoccupied voxels are 0.
fn dense<'py, T: numpy::Element + Copy + Default>(
    py: Python<'py>,
    r: &VoxelReport,
    shape: &[usize],
    get: impl Fn(&lindhard::tally::VoxelRecord) -> T,
) -> PyResult<Bound<'py, PyAny>> {
    let [nx, ny, nz] = r.grid.dims();
    let mut v = vec![T::default(); nx * ny * nz];
    for rec in &r.voxels {
        let [i, j, k] = rec.index;
        v[(i * ny + j) * nz + k] = get(rec);
    }
    Ok(v.into_pyarray(py).reshape(shape.to_vec())?.into_any())
}

/// `[nspecies, nx, ny, nz]` array of one defect kind.
fn species_array<'py>(
    py: Python<'py>,
    r: &VoxelReport,
    get: impl Fn(&lindhard::tally::SpeciesDefects) -> u64,
) -> PyResult<Bound<'py, PyAny>> {
    let [nx, ny, nz] = r.grid.dims();
    let ns = r.species_z.len();
    let mut v = vec![0u64; ns * nx * ny * nz];
    for rec in &r.voxels {
        let [i, j, k] = rec.index;
        for (s, d) in rec.defects.iter().enumerate() {
            v[((s * nx + i) * ny + j) * nz + k] = get(d);
        }
    }
    Ok(v.into_pyarray(py).reshape(vec![ns, nx, ny, nz])?.into_any())
}

/// The dict documented in `docs/python.md` ("Damage export").
pub(crate) fn to_dict<'py>(
    py: Python<'py>,
    r: &VoxelReport,
    p: &ExportProvenance,
) -> PyResult<Bound<'py, PyDict>> {
    let [nx, ny, nz] = r.grid.dims();
    let shape = [nx, ny, nz];
    let d = PyDict::new(py);
    d.set_item("origin_m", r.grid.origin_m().to_vec().into_pyarray(py))?;
    d.set_item("spacing_m", r.grid.spacing_m().to_vec().into_pyarray(py))?;
    d.set_item("dims", (nx, ny, nz))?;
    d.set_item("histories", r.histories)?;
    d.set_item("species_z", PyList::new(py, r.species_z.iter().copied())?)?;
    d.set_item("implanted", dense(py, r, &shape, |v| v.cell.implanted)?)?;
    d.set_item(
        "displacements",
        dense(py, r, &shape, |v| v.cell.displacements)?,
    )?;
    d.set_item(
        "electronic_ev",
        dense(py, r, &shape, |v| v.cell.electronic_ev)?,
    )?;
    d.set_item("lattice_ev", dense(py, r, &shape, |v| v.cell.lattice_ev)?)?;
    d.set_item("rest_ev", dense(py, r, &shape, |v| v.cell.rest_ev)?)?;
    d.set_item(
        "deposited_ev",
        dense(py, r, &shape, |v| v.cell.deposited_ev())?,
    )?;
    d.set_item("vacancies", species_array(py, r, |s| s.vacancies)?)?;
    d.set_item("interstitials", species_array(py, r, |s| s.interstitials)?)?;
    d.set_item("replacements", species_array(py, r, |s| s.replacements)?)?;
    let o = &r.outside;
    let outside = PyDict::new(py);
    outside.set_item("implanted", o.implanted)?;
    outside.set_item("displacements", o.displacements)?;
    outside.set_item("electronic_ev", o.electronic_ev)?;
    outside.set_item("lattice_ev", o.lattice_ev)?;
    outside.set_item("rest_ev", o.rest_ev)?;
    outside.set_item("deposited_ev", o.deposited_ev())?;
    let od: Vec<u64> = r
        .outside_defects
        .iter()
        .flat_map(|s| [s.vacancies, s.interstitials, s.replacements])
        .collect();
    outside.set_item(
        "defects",
        od.into_pyarray(py).reshape(vec![r.species_z.len(), 3])?,
    )?;
    d.set_item("outside", outside)?;
    let prov = PyDict::new(py);
    prov.set_item("format_version", lindhard::tally::EXPORT_FORMAT_VERSION)?;
    prov.set_item("engine", &p.engine)?;
    prov.set_item("lindhard_version", &p.lindhard_version)?;
    prov.set_item("seed", p.seed)?;
    prov.set_item("input_sha256", &p.input_sha256)?;
    d.set_item("provenance", prov)?;
    Ok(d)
}

/// Read a damage-export file (written by `lindhard run` or by
/// `RunResult.write`) into a dict of NumPy arrays.
#[pyfunction]
pub(crate) fn read_damage_export(py: Python<'_>, path: PathBuf) -> PyResult<Bound<'_, PyDict>> {
    let text = std::fs::read_to_string(&path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            pyo3::exceptions::PyFileNotFoundError::new_err(format!("{}: {e}", path.display()))
        } else {
            pyo3::exceptions::PyOSError::new_err(format!("{}: {e}", path.display()))
        }
    })?;
    let (r, p) = VoxelReport::read_export_csv(&text)
        .map_err(|e| errors::input(format!("{}: {e}", path.display())))?;
    to_dict(py, &r, &p)
}
