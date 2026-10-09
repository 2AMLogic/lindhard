//! The configuration classes: thin, mutable value wrappers over the library's
//! input schema ([`lindhard::input`]). Each converts to and from its spec
//! type; validation of values is the library's ([`Input::resolve`]), run when
//! a [`crate::run::Run`] is resolved, so the messages are the CLI's.

use std::collections::BTreeMap;

use lindhard::input::{
    BeamSpec, DivergenceSpec, EnergyOverride, LayerSpec, MaterialRef, PhysicsSpec, TallySpec,
    TargetSpec,
};
use lindhard::material::{ElementSpec, MaterialSpec};
use pyo3::prelude::*;
use serde::de::DeserializeOwned;
use serde::Serialize;

use crate::errors;

/// The kebab-case name of a unit-variant enum, as written in the TOML.
fn enum_name<T: Serialize>(v: &T) -> String {
    match toml::Value::try_from(v) {
        Ok(toml::Value::String(s)) => s,
        _ => unreachable!("choice enums serialize to strings"),
    }
}

/// Parse a choice from its TOML name, naming the field on failure.
fn enum_parse<T: DeserializeOwned>(field: &str, s: &str) -> PyResult<T> {
    toml::Value::String(s.to_owned())
        .try_into()
        .map_err(|e: toml::de::Error| errors::input(format!("{field}: {e}")))
}

/// One element of a material.
///
/// Give exactly one of `symbol` or `z`, and exactly one of `atom_fraction`
/// (or a stoichiometric count) or `mass_fraction`; all elements of a material
/// use the same kind. `e_d_ev`, `e_b_ev` and `e_s_ev` override the displacement,
/// lattice-binding and surface-binding energies (eV).
#[pyclass(module = "lindhard", get_all, set_all, from_py_object)]
#[derive(Debug, Clone, Default)]
pub struct Element {
    pub symbol: Option<String>,
    pub z: Option<u8>,
    pub atom_fraction: Option<f64>,
    pub mass_fraction: Option<f64>,
    pub e_d_ev: Option<f64>,
    pub e_b_ev: Option<f64>,
    pub e_s_ev: Option<f64>,
}

#[pymethods]
impl Element {
    #[new]
    #[pyo3(signature = (symbol=None, *, z=None, atom_fraction=None, mass_fraction=None, e_d_ev=None, e_b_ev=None, e_s_ev=None))]
    fn new(
        symbol: Option<String>,
        z: Option<u8>,
        atom_fraction: Option<f64>,
        mass_fraction: Option<f64>,
        e_d_ev: Option<f64>,
        e_b_ev: Option<f64>,
        e_s_ev: Option<f64>,
    ) -> Self {
        Self {
            symbol,
            z,
            atom_fraction,
            mass_fraction,
            e_d_ev,
            e_b_ev,
            e_s_ev,
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.to_spec())
    }
}

impl Element {
    pub fn to_spec(&self) -> ElementSpec {
        ElementSpec {
            symbol: self.symbol.clone(),
            z: self.z,
            atom_fraction: self.atom_fraction,
            mass_fraction: self.mass_fraction,
            e_d_ev: self.e_d_ev,
            e_b_ev: self.e_b_ev,
            e_s_ev: self.e_s_ev,
        }
    }

    pub fn from_spec(s: &ElementSpec) -> Self {
        Self {
            symbol: s.symbol.clone(),
            z: s.z,
            atom_fraction: s.atom_fraction,
            mass_fraction: s.mass_fraction,
            e_d_ev: s.e_d_ev,
            e_b_ev: s.e_b_ev,
            e_s_ev: s.e_s_ev,
        }
    }
}

/// A material: elements with fractions, and a mass density (required for
/// compounds; optional for a single element, which then takes its tabulated
/// density).
#[pyclass(module = "lindhard", get_all, set_all, skip_from_py_object)]
#[derive(Debug, Clone, Default)]
pub struct Material {
    pub elements: Vec<Element>,
    pub density_g_cm3: Option<f64>,
    pub name: Option<String>,
}

#[pymethods]
impl Material {
    #[new]
    #[pyo3(signature = (elements, density_g_cm3=None, name=None))]
    fn new(elements: Vec<Element>, density_g_cm3: Option<f64>, name: Option<String>) -> Self {
        Self {
            elements,
            density_g_cm3,
            name,
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.to_spec())
    }
}

impl Material {
    pub fn to_spec(&self) -> MaterialSpec {
        MaterialSpec {
            name: self.name.clone(),
            density_g_cm3: self.density_g_cm3,
            elements: self.elements.iter().map(Element::to_spec).collect(),
        }
    }

    pub fn from_spec(s: &MaterialSpec) -> Self {
        Self {
            elements: s.elements.iter().map(Element::from_spec).collect(),
            density_g_cm3: s.density_g_cm3,
            name: s.name.clone(),
        }
    }
}

/// Convert a Python `str` (a name from the run's materials, or an element
/// symbol) or [`Material`] to a [`MaterialRef`].
pub(crate) fn material_ref(field: &str, obj: &Bound<'_, PyAny>) -> PyResult<MaterialRef> {
    if let Ok(s) = obj.extract::<String>() {
        Ok(MaterialRef::Name(s))
    } else if let Ok(m) = obj.extract::<PyRef<'_, Material>>() {
        Ok(MaterialRef::Inline(m.to_spec()))
    } else {
        Err(errors::input(format!(
            "{field}: expected a material name (str) or a Material, got {}",
            obj.get_type().name()?
        )))
    }
}

fn material_obj(py: Python<'_>, r: &MaterialRef) -> PyResult<Py<PyAny>> {
    Ok(match r {
        MaterialRef::Name(s) => s.into_pyobject(py)?.into_any().unbind(),
        MaterialRef::Inline(m) => Bound::new(py, Material::from_spec(m))?.into_any().unbind(),
    })
}

/// One finite layer: a material (a name, or an inline `Material`) and a
/// thickness in nm.
#[pyclass(module = "lindhard")]
#[derive(Debug)]
pub struct Layer {
    material: MaterialRef,
    /// Thickness, nm.
    #[pyo3(get, set)]
    pub thickness_nm: f64,
}

#[pymethods]
impl Layer {
    #[new]
    fn new(material: &Bound<'_, PyAny>, thickness_nm: f64) -> PyResult<Self> {
        Ok(Self {
            material: material_ref("layer.material", material)?,
            thickness_nm,
        })
    }

    /// The layer material: a `str` (a name or an element symbol) or a `Material`.
    #[getter]
    fn material(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        material_obj(py, &self.material)
    }

    #[setter]
    fn set_material(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        self.material = material_ref("layer.material", value)?;
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!(
            "Layer({:?}, thickness_nm={})",
            self.material, self.thickness_nm
        )
    }
}

impl Layer {
    pub fn to_spec(&self) -> LayerSpec {
        LayerSpec {
            material: self.material.clone(),
            thickness_nm: self.thickness_nm,
        }
    }

    pub fn from_spec(s: &LayerSpec) -> Self {
        Self {
            material: s.material.clone(),
            thickness_nm: s.thickness_nm,
        }
    }
}

/// The layered target: finite layers front to back, then an optional
/// semi-infinite substrate. At least one of the two is required; without a
/// substrate the target has a back face and particles can be transmitted.
#[pyclass(module = "lindhard")]
#[derive(Debug, Default)]
pub struct Target {
    layers: Vec<LayerSpec>,
    substrate: Option<MaterialRef>,
}

#[pymethods]
impl Target {
    #[new]
    #[pyo3(signature = (layers=None, substrate=None))]
    fn new(
        layers: Option<Vec<PyRef<'_, Layer>>>,
        substrate: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let mut t = Self {
            layers: layers
                .unwrap_or_default()
                .iter()
                .map(|l| l.to_spec())
                .collect(),
            substrate: None,
        };
        t.set_substrate(substrate)?;
        Ok(t)
    }

    /// The finite layers, front first. Returns copies: assign a new list to
    /// change them.
    #[getter]
    fn layers(&self, py: Python<'_>) -> PyResult<Vec<Py<Layer>>> {
        self.layers
            .iter()
            .map(|l| Py::new(py, Layer::from_spec(l)))
            .collect()
    }

    #[setter]
    fn set_layers(&mut self, layers: Vec<PyRef<'_, Layer>>) {
        self.layers = layers.iter().map(|l| l.to_spec()).collect();
    }

    /// The semi-infinite substrate (a `str` or `Material`), or `None`.
    #[getter]
    fn substrate(&self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        self.substrate
            .as_ref()
            .map(|m| material_obj(py, m))
            .transpose()
    }

    #[setter]
    fn set_substrate(&mut self, value: Option<&Bound<'_, PyAny>>) -> PyResult<()> {
        self.substrate = match value {
            None => None,
            Some(v) if v.is_none() => None,
            Some(v) => Some(material_ref("target.substrate", v)?),
        };
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.to_spec())
    }
}

impl Target {
    pub fn to_spec(&self) -> TargetSpec {
        TargetSpec {
            layers: self.layers.clone(),
            substrate: self.substrate.clone(),
        }
    }

    pub fn from_spec(s: &TargetSpec) -> Self {
        Self {
            layers: s.layers.clone(),
            substrate: s.substrate.clone(),
        }
    }
}

/// The incident beam: element symbol, energy (eV), mass (u; default the
/// standard atomic weight), polar angle of incidence from the surface normal
/// in `[0, 90)` degrees, and azimuth of the incidence plane (degrees).
#[pyclass(module = "lindhard", get_all, set_all, skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct Beam {
    pub ion: String,
    pub energy_ev: f64,
    pub mass_amu: Option<f64>,
    pub tilt_deg: f64,
    pub azimuth_deg: f64,
    pub divergence_model: Option<String>,
    pub divergence_deg: Option<f64>,
}

#[pymethods]
impl Beam {
    #[new]
    #[pyo3(signature = (ion, energy_ev, mass_amu=None, tilt_deg=0.0, azimuth_deg=0.0, divergence_model=None, divergence_deg=None))]
    fn new(
        ion: String,
        energy_ev: f64,
        mass_amu: Option<f64>,
        tilt_deg: f64,
        azimuth_deg: f64,
        divergence_model: Option<String>,
        divergence_deg: Option<f64>,
    ) -> Self {
        Self {
            ion,
            energy_ev,
            mass_amu,
            tilt_deg,
            azimuth_deg,
            divergence_model,
            divergence_deg,
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.to_spec())
    }
}

impl Beam {
    pub fn to_spec(&self) -> BeamSpec {
        BeamSpec {
            ion: self.ion.clone(),
            mass_amu: self.mass_amu,
            energy_ev: self.energy_ev,
            tilt_deg: self.tilt_deg,
            azimuth_deg: self.azimuth_deg,
            divergence: self.divergence_spec(),
        }
    }

    /// The `[beam.divergence]` table. An unknown model name or a model with no
    /// width cannot be represented; it becomes a Gaussian / cone of a NaN
    /// width so that resolution reports `beam.divergence.*` as invalid.
    fn divergence_spec(&self) -> Option<DivergenceSpec> {
        let model = self.divergence_model.as_deref()?;
        let w = self.divergence_deg.unwrap_or(f64::NAN);
        Some(match model {
            "uniform-cone" => DivergenceSpec::UniformCone { half_angle_deg: w },
            // "gaussian" and anything else: resolved as Gaussian.
            _ => DivergenceSpec::Gaussian { sigma_deg: w },
        })
    }

    pub fn from_spec(s: &BeamSpec) -> Self {
        Self {
            ion: s.ion.clone(),
            energy_ev: s.energy_ev,
            mass_amu: s.mass_amu,
            tilt_deg: s.tilt_deg,
            azimuth_deg: s.azimuth_deg,
            divergence_model: s.divergence.map(|d| {
                match d {
                    DivergenceSpec::Gaussian { .. } => "gaussian",
                    DivergenceSpec::UniformCone { .. } => "uniform-cone",
                }
                .to_string()
            }),
            divergence_deg: s.divergence.map(|d| match d {
                DivergenceSpec::Gaussian { sigma_deg } => sigma_deg,
                DivergenceSpec::UniformCone { half_angle_deg } => half_angle_deg,
            }),
        }
    }
}

const ENERGY_KEYS: [&str; 3] = ["e_d_ev", "e_b_ev", "e_s_ev"];

/// Model choices and cutoffs; the `[physics]` table of the input.
///
/// `potential`, `screening_length`, `stopping` and `free_path` are the TOML
/// names (`"zbl"`, `"kr-c"`, `"lindhard-scharff"`, `"constant"`, ...). Both
/// cutoffs (eV) are required. `energies` maps an element symbol to a dict
/// with any of `e_d_ev`, `e_b_ev`, `e_s_ev`, overriding that element's
/// energies in every layer.
#[pyclass(module = "lindhard", get_all, set_all, skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct Physics {
    pub primary_cutoff_ev: f64,
    pub recoil_cutoff_ev: f64,
    pub potential: String,
    pub screening_length: Option<String>,
    pub stopping: String,
    pub free_path: String,
    pub min_cm_angle_deg: Option<f64>,
    pub weak_collisions: u8,
    pub follow_recoils: bool,
    pub primary_surface_binding_ev: f64,
    pub energies: BTreeMap<String, BTreeMap<String, f64>>,
}

#[pymethods]
impl Physics {
    #[new]
    #[pyo3(signature = (
        primary_cutoff_ev, recoil_cutoff_ev, *, potential="zbl".to_string(),
        screening_length=None, stopping="lindhard-scharff".to_string(),
        free_path="constant".to_string(), min_cm_angle_deg=None, weak_collisions=0,
        follow_recoils=true, primary_surface_binding_ev=0.0, energies=None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        primary_cutoff_ev: f64,
        recoil_cutoff_ev: f64,
        potential: String,
        screening_length: Option<String>,
        stopping: String,
        free_path: String,
        min_cm_angle_deg: Option<f64>,
        weak_collisions: u8,
        follow_recoils: bool,
        primary_surface_binding_ev: f64,
        energies: Option<BTreeMap<String, BTreeMap<String, f64>>>,
    ) -> Self {
        Self {
            primary_cutoff_ev,
            recoil_cutoff_ev,
            potential,
            screening_length,
            stopping,
            free_path,
            min_cm_angle_deg,
            weak_collisions,
            follow_recoils,
            primary_surface_binding_ev,
            energies: energies.unwrap_or_default(),
        }
    }

    fn __repr__(&self) -> String {
        match self.to_spec() {
            Ok(s) => format!("{s:?}"),
            Err(_) => "Physics(<invalid>)".to_string(),
        }
    }
}

impl Physics {
    pub fn to_spec(&self) -> PyResult<PhysicsSpec> {
        let mut energies = BTreeMap::new();
        for (el, m) in &self.energies {
            for k in m.keys() {
                if !ENERGY_KEYS.contains(&k.as_str()) {
                    return Err(errors::input(format!(
                        "physics.energies.{el}: unknown key `{k}` (expected one of {})",
                        ENERGY_KEYS.join(", ")
                    )));
                }
            }
            energies.insert(
                el.clone(),
                EnergyOverride {
                    e_d_ev: m.get("e_d_ev").copied(),
                    e_b_ev: m.get("e_b_ev").copied(),
                    e_s_ev: m.get("e_s_ev").copied(),
                },
            );
        }
        Ok(PhysicsSpec {
            potential: enum_parse("physics.potential", &self.potential)?,
            screening_length: self
                .screening_length
                .as_deref()
                .map(|s| enum_parse("physics.screening_length", s))
                .transpose()?,
            stopping: enum_parse("physics.stopping", &self.stopping)?,
            free_path: enum_parse("physics.free_path", &self.free_path)?,
            min_cm_angle_deg: self.min_cm_angle_deg,
            weak_collisions: self.weak_collisions,
            primary_cutoff_ev: self.primary_cutoff_ev,
            recoil_cutoff_ev: self.recoil_cutoff_ev,
            follow_recoils: self.follow_recoils,
            primary_surface_binding_ev: self.primary_surface_binding_ev,
            energies,
            // Tuning is CLI/TOML-only in the pilot; Python runs are untuned.
            tuning: lindhard::input::NO_TUNING.to_string(),
        })
    }

    pub fn from_spec(s: &PhysicsSpec) -> Self {
        let mut energies = BTreeMap::new();
        for (el, o) in &s.energies {
            let mut m = BTreeMap::new();
            for (k, v) in ENERGY_KEYS.iter().zip([o.e_d_ev, o.e_b_ev, o.e_s_ev]) {
                if let Some(v) = v {
                    m.insert((*k).to_string(), v);
                }
            }
            energies.insert(el.clone(), m);
        }
        Self {
            primary_cutoff_ev: s.primary_cutoff_ev,
            recoil_cutoff_ev: s.recoil_cutoff_ev,
            potential: enum_name(&s.potential),
            screening_length: s.screening_length.as_ref().map(enum_name),
            stopping: enum_name(&s.stopping),
            free_path: enum_name(&s.free_path),
            min_cm_angle_deg: s.min_cm_angle_deg,
            weak_collisions: s.weak_collisions,
            follow_recoils: s.follow_recoils,
            primary_surface_binding_ev: s.primary_surface_binding_ev,
            energies,
        }
    }
}

/// What the run records; the `[tally]` table of the input. Defaults are the
/// CLI's.
#[pyclass(module = "lindhard", get_all, set_all, skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct Tally {
    pub depth_bin_nm: f64,
    pub depth_bins: usize,
    pub per_ion: bool,
    pub lateral_bin_nm: f64,
    pub lateral_bins: usize,
    pub escape_energy_max_ev: Option<f64>,
    pub escape_energy_bins: usize,
    pub escape_polar_bins: usize,
    pub dual_pearson: bool,
}

#[pymethods]
impl Tally {
    #[new]
    #[pyo3(signature = (
        depth_bin_nm=None, depth_bins=None, per_ion=None, lateral_bin_nm=None,
        lateral_bins=None, escape_energy_max_ev=None, escape_energy_bins=None,
        escape_polar_bins=None, dual_pearson=None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        depth_bin_nm: Option<f64>,
        depth_bins: Option<usize>,
        per_ion: Option<bool>,
        lateral_bin_nm: Option<f64>,
        lateral_bins: Option<usize>,
        escape_energy_max_ev: Option<f64>,
        escape_energy_bins: Option<usize>,
        escape_polar_bins: Option<usize>,
        dual_pearson: Option<bool>,
    ) -> Self {
        let d = TallySpec::default();
        Self {
            depth_bin_nm: depth_bin_nm.unwrap_or(d.depth_bin_nm),
            depth_bins: depth_bins.unwrap_or(d.depth_bins),
            per_ion: per_ion.unwrap_or(d.per_ion),
            lateral_bin_nm: lateral_bin_nm.unwrap_or(d.lateral_bin_nm),
            lateral_bins: lateral_bins.unwrap_or(d.lateral_bins),
            escape_energy_max_ev,
            escape_energy_bins: escape_energy_bins.unwrap_or(d.escape_energy_bins),
            escape_polar_bins: escape_polar_bins.unwrap_or(d.escape_polar_bins),
            dual_pearson: dual_pearson.unwrap_or(d.dual_pearson),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.to_spec())
    }
}

impl Tally {
    pub fn to_spec(&self) -> TallySpec {
        TallySpec {
            depth_bin_nm: self.depth_bin_nm,
            depth_bins: self.depth_bins,
            per_ion: self.per_ion,
            lateral_bin_nm: self.lateral_bin_nm,
            lateral_bins: self.lateral_bins,
            escape_energy_max_ev: self.escape_energy_max_ev,
            escape_energy_bins: self.escape_energy_bins,
            escape_polar_bins: self.escape_polar_bins,
            dual_pearson: self.dual_pearson,
        }
    }

    pub fn from_spec(s: &TallySpec) -> Self {
        Self {
            depth_bin_nm: s.depth_bin_nm,
            depth_bins: s.depth_bins,
            per_ion: s.per_ion,
            lateral_bin_nm: s.lateral_bin_nm,
            lateral_bins: s.lateral_bins,
            escape_energy_max_ev: s.escape_energy_max_ev,
            escape_energy_bins: s.escape_energy_bins,
            escape_polar_bins: s.escape_polar_bins,
            dual_pearson: s.dual_pearson,
        }
    }
}
