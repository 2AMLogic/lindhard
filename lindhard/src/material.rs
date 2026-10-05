//! Materials: elements combined by atom or mass fraction, with a mass density,
//! a derived atom number density, and per-element displacement (`E_d`), lattice
//! binding (`E_b`) and surface binding (`E_s`) energies.
//!
//! # Conventions (read this before trusting a default)
//!
//! * Fractions passed in are **relative weights**: they are normalised by
//!   dividing by their sum, so stoichiometric counts (`Si 1, O 2`) are valid.
//!   The sum must be finite and strictly positive and every entry finite and
//!   non-negative; otherwise construction fails. Normalisation is exact
//!   division, with no tolerance window.
//! * `E_d`, `E_b`, `E_s` are **model parameters, not measured properties of
//!   your sample**. Defaults come from [`crate::elements`] and are conventions
//!   with citations there: `E_s` defaults to the elemental cohesive
//!   (sublimation) energy where one is tabulated; `E_d` defaults to the ASTM
//!   E521 convention value where one is recorded; `E_b` defaults to 0 eV, the
//!   common choice for amorphous-target binary-collision runs. Where there is no
//!   default (`None`) the getter returns [`MaterialError::EnergyNotSet`] so the
//!   user has to decide. Energies are stored in eV (suffix `_ev`); see
//!   [`crate::units`].
//! * Elements keep the order they were given in, so every output derived from a
//!   material is deterministic.
//! * Atom number density follows `n_i = x_i * rho * N_A / sum_j(x_j * M_j)`
//!   with `x` the atom fractions, `rho` the mass density and `M` the atomic
//!   weights (CODATA 2022 `N_A`).

use serde::{Deserialize, Serialize};

use crate::constants::AVOGADRO;
use crate::elements::{element, element_by_symbol};
use crate::units::{g_cm3_to_kg_m3, kg_m3_to_g_cm3};

/// Which per-element energy a message is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyKind {
    /// Displacement energy `E_d`.
    Displacement,
    /// Lattice binding energy `E_b`.
    LatticeBinding,
    /// Surface binding energy `E_s`.
    SurfaceBinding,
}

impl std::fmt::Display for EnergyKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Displacement => "displacement energy E_d",
            Self::LatticeBinding => "lattice binding energy E_b",
            Self::SurfaceBinding => "surface binding energy E_s",
        })
    }
}

/// Errors from building or querying a [`Material`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MaterialError {
    /// No elements were given.
    #[error("material has no elements")]
    Empty,
    /// Atomic number outside 1..=92.
    #[error("unknown atomic number {0}")]
    UnknownAtomicNumber(u8),
    /// Symbol not in the element table.
    #[error("unknown element symbol {0:?}")]
    UnknownSymbol(String),
    /// The same element appeared twice.
    #[error("element Z={0} listed more than once")]
    DuplicateElement(u8),
    /// A fraction was negative, NaN or infinite.
    #[error("invalid fraction {value} for element Z={z}: must be finite and non-negative")]
    InvalidFraction {
        /// Atomic number.
        z: u8,
        /// Offending value.
        value: f64,
    },
    /// Fractions sum to zero or to a non-finite value, so cannot be normalised.
    #[error("fractions sum to {0}: cannot be normalised")]
    NotNormalizable(f64),
    /// A density was zero, negative or not finite.
    #[error("invalid density {0} kg/m^3: must be finite and positive")]
    InvalidDensity(f64),
    /// No density given and none can be inferred (compound, or element with no
    /// tabulated solid density).
    #[error("no density given and none can be inferred for this composition")]
    MissingDensity,
    /// An energy was negative, NaN or infinite.
    #[error("invalid {kind} {value} eV for element Z={z}: must be finite and non-negative")]
    InvalidEnergy {
        /// Atomic number.
        z: u8,
        /// Which energy.
        kind: EnergyKind,
        /// Offending value.
        value: f64,
    },
    /// The element is not part of the material.
    #[error("element Z={0} is not in this material")]
    ElementNotInMaterial(u8),
    /// No default exists for this energy and the user has not set one.
    #[error("{kind} for element Z={z} has no default; set it explicitly")]
    EnergyNotSet {
        /// Atomic number.
        z: u8,
        /// Which energy.
        kind: EnergyKind,
    },
    /// A serialized element entry was malformed (spec level problems).
    #[error("invalid material specification: {0}")]
    InvalidSpec(String),
}

/// One element within a material.
#[derive(Debug, Clone, PartialEq)]
pub struct Component {
    z: u8,
    atom_fraction: f64,
    e_d_ev: Option<f64>,
    e_b_ev: Option<f64>,
    e_s_ev: Option<f64>,
}

impl Component {
    /// Atomic number.
    pub fn z(&self) -> u8 {
        self.z
    }
    /// Normalised atom fraction.
    pub fn atom_fraction(&self) -> f64 {
        self.atom_fraction
    }
}

/// A material: a list of elements, a mass density and per-element energies.
///
/// Serializes (serde) to the canonical form of [`MaterialSpec`] with atom
/// fractions and explicit densities and energies; deserializes from any valid
/// spec, so a TOML file may give atom or mass fractions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "MaterialSpec", into = "MaterialSpec")]
pub struct Material {
    name: Option<String>,
    components: Vec<Component>,
    density_kg_m3: f64,
}

fn check_fraction(z: u8, v: f64) -> Result<(), MaterialError> {
    if v.is_finite() && v >= 0.0 {
        Ok(())
    } else {
        Err(MaterialError::InvalidFraction { z, value: v })
    }
}

fn check_energy(z: u8, kind: EnergyKind, v: f64) -> Result<f64, MaterialError> {
    if v.is_finite() && v >= 0.0 {
        Ok(v)
    } else {
        Err(MaterialError::InvalidEnergy { z, kind, value: v })
    }
}

/// Validate, normalise and (if mass fractions) convert to atom fractions.
fn normalise(parts: &[(u8, f64)], by_mass: bool) -> Result<Vec<(u8, f64)>, MaterialError> {
    if parts.is_empty() {
        return Err(MaterialError::Empty);
    }
    let mut out: Vec<(u8, f64)> = Vec::with_capacity(parts.len());
    for &(z, f) in parts {
        let el = element(z).ok_or(MaterialError::UnknownAtomicNumber(z))?;
        check_fraction(z, f)?;
        if out.iter().any(|&(zz, _)| zz == z) {
            return Err(MaterialError::DuplicateElement(z));
        }
        out.push((z, if by_mass { f / el.atomic_weight } else { f }));
    }
    let sum: f64 = out.iter().map(|&(_, f)| f).sum();
    if !(sum.is_finite() && sum > 0.0) {
        return Err(MaterialError::NotNormalizable(sum));
    }
    for p in &mut out {
        p.1 /= sum;
    }
    Ok(out)
}

impl Material {
    fn build(
        parts: &[(u8, f64)],
        by_mass: bool,
        density_kg_m3: Option<f64>,
    ) -> Result<Self, MaterialError> {
        let norm = normalise(parts, by_mass)?;
        let density = match density_kg_m3 {
            Some(d) if d.is_finite() && d > 0.0 => d,
            Some(d) => return Err(MaterialError::InvalidDensity(d)),
            None if norm.len() == 1 => element(norm[0].0)
                .and_then(|e| e.density_kg_m3())
                .ok_or(MaterialError::MissingDensity)?,
            None => return Err(MaterialError::MissingDensity),
        };
        let components = norm
            .into_iter()
            .map(|(z, x)| {
                // `normalise` already checked the element exists.
                let el = element(z).expect("validated");
                Component {
                    z,
                    atom_fraction: x,
                    e_d_ev: el.default_displacement_ev,
                    e_b_ev: Some(0.0),
                    e_s_ev: el.default_surface_binding_ev,
                }
            })
            .collect();
        Ok(Self {
            name: None,
            components,
            density_kg_m3: density,
        })
    }

    /// Build from `(Z, atom fraction)` pairs. Fractions are relative weights
    /// (see the module docs). `density_kg_m3` may be `None` only for a single
    /// element with a tabulated density.
    pub fn from_atom_fractions(
        parts: &[(u8, f64)],
        density_kg_m3: Option<f64>,
    ) -> Result<Self, MaterialError> {
        Self::build(parts, false, density_kg_m3)
    }

    /// Build from `(Z, mass fraction)` pairs. Same rules as
    /// [`Material::from_atom_fractions`].
    pub fn from_mass_fractions(
        parts: &[(u8, f64)],
        density_kg_m3: Option<f64>,
    ) -> Result<Self, MaterialError> {
        Self::build(parts, true, density_kg_m3)
    }

    /// Set a display name (does not affect physics).
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Display name, if set.
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Components in the order they were given.
    pub fn components(&self) -> &[Component] {
        &self.components
    }

    /// Mass density, kg/m³.
    pub fn mass_density(&self) -> f64 {
        self.density_kg_m3
    }

    /// Override the mass density, kg/m³.
    pub fn set_mass_density(&mut self, density_kg_m3: f64) -> Result<(), MaterialError> {
        if !(density_kg_m3.is_finite() && density_kg_m3 > 0.0) {
            return Err(MaterialError::InvalidDensity(density_kg_m3));
        }
        self.density_kg_m3 = density_kg_m3;
        Ok(())
    }

    /// Mean atomic weight `sum(x_j M_j)`, g/mol.
    pub fn mean_atomic_weight(&self) -> f64 {
        self.components
            .iter()
            .map(|c| c.atom_fraction * element(c.z).expect("validated").atomic_weight)
            .sum()
    }

    /// Total atom number density, atoms per m³.
    pub fn atom_number_density(&self) -> f64 {
        // rho [kg/m^3] -> g/m^3 is 1e3; divide by g/mol for mol/m^3.
        self.density_kg_m3 * 1.0e3 * AVOGADRO / self.mean_atomic_weight()
    }

    /// Number density of one element, atoms per m³.
    pub fn number_density_of(&self, z: u8) -> Result<f64, MaterialError> {
        let c = self.component(z)?;
        Ok(c.atom_fraction * self.atom_number_density())
    }

    /// Atom fractions in component order.
    pub fn atom_fractions(&self) -> Vec<(u8, f64)> {
        self.components
            .iter()
            .map(|c| (c.z, c.atom_fraction))
            .collect()
    }

    /// Mass fractions in component order.
    pub fn mass_fractions(&self) -> Vec<(u8, f64)> {
        let mean = self.mean_atomic_weight();
        self.components
            .iter()
            .map(|c| {
                let m = element(c.z).expect("validated").atomic_weight;
                (c.z, c.atom_fraction * m / mean)
            })
            .collect()
    }

    fn component(&self, z: u8) -> Result<&Component, MaterialError> {
        self.components
            .iter()
            .find(|c| c.z == z)
            .ok_or(MaterialError::ElementNotInMaterial(z))
    }

    fn component_mut(&mut self, z: u8) -> Result<&mut Component, MaterialError> {
        self.components
            .iter_mut()
            .find(|c| c.z == z)
            .ok_or(MaterialError::ElementNotInMaterial(z))
    }

    /// Displacement energy `E_d` of element `z` in this material, eV.
    pub fn displacement_energy_ev(&self, z: u8) -> Result<f64, MaterialError> {
        self.component(z)?
            .e_d_ev
            .ok_or(MaterialError::EnergyNotSet {
                z,
                kind: EnergyKind::Displacement,
            })
    }

    /// Lattice binding energy `E_b` of element `z`, eV.
    pub fn lattice_binding_energy_ev(&self, z: u8) -> Result<f64, MaterialError> {
        self.component(z)?
            .e_b_ev
            .ok_or(MaterialError::EnergyNotSet {
                z,
                kind: EnergyKind::LatticeBinding,
            })
    }

    /// Surface binding energy `E_s` of element `z`, eV.
    pub fn surface_binding_energy_ev(&self, z: u8) -> Result<f64, MaterialError> {
        self.component(z)?
            .e_s_ev
            .ok_or(MaterialError::EnergyNotSet {
                z,
                kind: EnergyKind::SurfaceBinding,
            })
    }

    /// Set `E_d` for element `z`, eV.
    pub fn set_displacement_energy_ev(&mut self, z: u8, ev: f64) -> Result<(), MaterialError> {
        let v = check_energy(z, EnergyKind::Displacement, ev)?;
        self.component_mut(z)?.e_d_ev = Some(v);
        Ok(())
    }

    /// Set `E_b` for element `z`, eV.
    pub fn set_lattice_binding_energy_ev(&mut self, z: u8, ev: f64) -> Result<(), MaterialError> {
        let v = check_energy(z, EnergyKind::LatticeBinding, ev)?;
        self.component_mut(z)?.e_b_ev = Some(v);
        Ok(())
    }

    /// Set `E_s` for element `z`, eV.
    pub fn set_surface_binding_energy_ev(&mut self, z: u8, ev: f64) -> Result<(), MaterialError> {
        let v = check_energy(z, EnergyKind::SurfaceBinding, ev)?;
        self.component_mut(z)?.e_s_ev = Some(v);
        Ok(())
    }

    /// Energies that still have no value, in component order: a run should
    /// refuse to start while this is non-empty.
    pub fn unset_energies(&self) -> Vec<(u8, EnergyKind)> {
        let mut v = Vec::new();
        for c in &self.components {
            if c.e_d_ev.is_none() {
                v.push((c.z, EnergyKind::Displacement));
            }
            if c.e_b_ev.is_none() {
                v.push((c.z, EnergyKind::LatticeBinding));
            }
            if c.e_s_ev.is_none() {
                v.push((c.z, EnergyKind::SurfaceBinding));
            }
        }
        v
    }
}

/// Serialized form of one element entry. Give exactly one of `symbol` / `z`
/// and exactly one of `atom_fraction` / `mass_fraction`; every element of a
/// material must use the same fraction kind.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElementSpec {
    /// Element symbol, e.g. `"Si"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symbol: Option<String>,
    /// Atomic number (alternative to `symbol`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub z: Option<u8>,
    /// Relative atom fraction (or stoichiometric count).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub atom_fraction: Option<f64>,
    /// Relative mass fraction.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mass_fraction: Option<f64>,
    /// Displacement energy override, eV.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub e_d_ev: Option<f64>,
    /// Lattice binding energy override, eV.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub e_b_ev: Option<f64>,
    /// Surface binding energy override, eV.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub e_s_ev: Option<f64>,
}

/// Serialized form of a [`Material`] (what the CLI's TOML describes).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialSpec {
    /// Display name.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// Mass density override, g/cm³. Required for compounds.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub density_g_cm3: Option<f64>,
    /// Elements, in order.
    pub elements: Vec<ElementSpec>,
}

impl TryFrom<MaterialSpec> for Material {
    type Error = MaterialError;

    fn try_from(spec: MaterialSpec) -> Result<Self, MaterialError> {
        let bad = |s: &str| MaterialError::InvalidSpec(s.to_string());
        let mut parts = Vec::with_capacity(spec.elements.len());
        let mut kind: Option<bool> = None; // Some(true) = by mass
        for e in &spec.elements {
            let z = match (&e.symbol, e.z) {
                (Some(s), None) => {
                    element_by_symbol(s)
                        .ok_or_else(|| MaterialError::UnknownSymbol(s.clone()))?
                        .z
                }
                (None, Some(z)) => z,
                (Some(s), Some(z)) => {
                    let el = element_by_symbol(s)
                        .ok_or_else(|| MaterialError::UnknownSymbol(s.clone()))?;
                    if el.z != z {
                        return Err(bad("symbol and z disagree"));
                    }
                    z
                }
                (None, None) => return Err(bad("element needs `symbol` or `z`")),
            };
            let (by_mass, f) = match (e.atom_fraction, e.mass_fraction) {
                (Some(f), None) => (false, f),
                (None, Some(f)) => (true, f),
                _ => return Err(bad("give exactly one of atom_fraction / mass_fraction")),
            };
            match kind {
                None => kind = Some(by_mass),
                Some(k) if k != by_mass => {
                    return Err(bad("cannot mix atom_fraction and mass_fraction"))
                }
                _ => {}
            }
            parts.push((z, f));
        }
        let density = spec.density_g_cm3.map(g_cm3_to_kg_m3);
        let mut m = Material::build(&parts, kind.unwrap_or(false), density)?;
        m.name = spec.name;
        for (e, &(z, _)) in spec.elements.iter().zip(&parts) {
            if let Some(v) = e.e_d_ev {
                m.set_displacement_energy_ev(z, v)?;
            }
            if let Some(v) = e.e_b_ev {
                m.set_lattice_binding_energy_ev(z, v)?;
            }
            if let Some(v) = e.e_s_ev {
                m.set_surface_binding_energy_ev(z, v)?;
            }
        }
        Ok(m)
    }
}

impl From<Material> for MaterialSpec {
    fn from(m: Material) -> Self {
        MaterialSpec {
            name: m.name,
            density_g_cm3: Some(kg_m3_to_g_cm3(m.density_kg_m3)),
            elements: m
                .components
                .iter()
                .map(|c| ElementSpec {
                    symbol: Some(element(c.z).expect("validated").symbol.to_string()),
                    z: None,
                    atom_fraction: Some(c.atom_fraction),
                    mass_fraction: None,
                    e_d_ev: c.e_d_ev,
                    e_b_ev: c.e_b_ev,
                    e_s_ev: c.e_s_ev,
                })
                .collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rel(a: f64, b: f64) -> f64 {
        ((a - b) / b).abs()
    }

    // Hand values: n = rho[kg/m^3] * 1e3 * N_A * (atoms per formula unit)
    //                  / (sum of atomic weights over the formula unit)
    // with N_A = 6.02214076e23, weights Si 28.085, O 15.999, Ga 69.723,
    // As 74.921595, C 12.011, and these densities (g/cm^3):
    //   Si   2.329  (table value, CRC)
    //   SiO2 2.648  (alpha-quartz, 20 C; crystallographic density)
    //   GaAs 5.3176 (zinc blende, 300 K; crystallographic density)
    //   SiC  3.21   (3C/6H SiC, room temperature)
    // Results, per m^3 (rho in g/m^3 = 1e6 * the g/cm^3 value):
    //   Si   4.99397e28 m^-3   (4.994e22 cm^-3)
    //   SiO2 3 * 2.648e6*N_A/60.083   = 7.9623e28 m^-3
    //   GaAs 2 * 5.3176e6*N_A/144.644595 = 4.4279e28 m^-3
    //   SiC  2 * 3.21e6*N_A/40.096    = 9.6424e28 m^-3
    #[test]
    fn si_number_density() {
        let m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        assert!(rel(m.atom_number_density(), 4.99397e28) < 1e-5);
    }

    #[test]
    fn sio2_number_density() {
        let m = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2648.0)).unwrap();
        assert!(rel(m.atom_number_density(), 7.9623e28) < 1e-4);
        assert!(rel(m.number_density_of(8).unwrap(), 7.9623e28 * 2.0 / 3.0) < 1e-4);
    }

    #[test]
    fn gaas_number_density() {
        let m = Material::from_atom_fractions(&[(31, 1.0), (33, 1.0)], Some(5317.6)).unwrap();
        assert!(rel(m.atom_number_density(), 4.4279e28) < 1e-4);
    }

    #[test]
    fn sic_number_density() {
        let m = Material::from_atom_fractions(&[(14, 1.0), (6, 1.0)], Some(3210.0)).unwrap();
        assert!(rel(m.atom_number_density(), 9.6424e28) < 1e-4);
    }

    #[test]
    fn mass_atom_fraction_round_trip() {
        let a = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2648.0)).unwrap();
        let w = a.mass_fractions();
        let sum: f64 = w.iter().map(|p| p.1).sum();
        assert!((sum - 1.0).abs() < 1e-14);
        // Hand value: w_Si = 28.085 / 60.083 = 0.467435...
        assert!((w[0].1 - 28.085 / 60.083).abs() < 1e-12);
        let b = Material::from_mass_fractions(&w, Some(2648.0)).unwrap();
        for (x, y) in a.atom_fractions().iter().zip(b.atom_fractions()) {
            assert_eq!(x.0, y.0);
            assert!((x.1 - y.1).abs() < 1e-14);
        }
        assert!(rel(a.atom_number_density(), b.atom_number_density()) < 1e-14);
    }

    #[test]
    fn error_paths() {
        use MaterialError::*;
        assert_eq!(Material::from_atom_fractions(&[], Some(1.0)), Err(Empty));
        assert_eq!(
            Material::from_atom_fractions(&[(0, 1.0)], Some(1.0)),
            Err(UnknownAtomicNumber(0))
        );
        assert_eq!(
            Material::from_atom_fractions(&[(93, 1.0)], Some(1.0)),
            Err(UnknownAtomicNumber(93))
        );
        assert!(matches!(
            Material::from_atom_fractions(&[(14, -1.0)], None),
            Err(InvalidFraction { z: 14, .. })
        ));
        assert!(matches!(
            Material::from_atom_fractions(&[(14, f64::NAN)], None),
            Err(InvalidFraction { .. })
        ));
        assert!(matches!(
            Material::from_atom_fractions(&[(14, 0.0), (8, 0.0)], Some(1.0)),
            Err(NotNormalizable(_))
        ));
        assert_eq!(
            Material::from_atom_fractions(&[(14, 1.0), (14, 1.0)], Some(1.0)),
            Err(DuplicateElement(14))
        );
        assert!(matches!(
            Material::from_atom_fractions(&[(14, 1.0)], Some(-5.0)),
            Err(InvalidDensity(_))
        ));
        assert_eq!(
            Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], None),
            Err(MissingDensity)
        );
        // Gas: no tabulated solid density.
        assert_eq!(
            Material::from_atom_fractions(&[(18, 1.0)], None),
            Err(MissingDensity)
        );
        let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        assert!(matches!(
            m.set_surface_binding_energy_ev(14, -1.0),
            Err(InvalidEnergy { .. })
        ));
        assert_eq!(
            m.set_displacement_energy_ev(8, 10.0),
            Err(ElementNotInMaterial(8))
        );
        // Oxygen has no default E_s: must be set by the user.
        let mut o = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2648.0)).unwrap();
        assert_eq!(
            o.surface_binding_energy_ev(8),
            Err(EnergyNotSet {
                z: 8,
                kind: EnergyKind::SurfaceBinding
            })
        );
        assert!(!o.unset_energies().is_empty());
        o.set_surface_binding_energy_ev(8, 2.0).unwrap();
        assert_eq!(o.surface_binding_energy_ev(8), Ok(2.0));
    }

    #[test]
    fn defaults_and_overrides() {
        let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        assert_eq!(m.surface_binding_energy_ev(14), Ok(4.63));
        assert_eq!(m.lattice_binding_energy_ev(14), Ok(0.0));
        m.set_surface_binding_energy_ev(14, 3.0).unwrap();
        assert_eq!(m.surface_binding_energy_ev(14), Ok(3.0));
    }

    const TOML_FIXTURE: &str = r#"
name = "SiO2"
density_g_cm3 = 2.648

[[elements]]
symbol = "Si"
mass_fraction = 28.085
e_d_ev = 21.0

[[elements]]
z = 8
mass_fraction = 31.998
e_s_ev = 2.0
"#;

    #[test]
    fn toml_fixture_and_serde_round_trip() {
        let m: Material = toml::from_str(TOML_FIXTURE).unwrap();
        assert_eq!(m.name(), Some("SiO2"));
        assert!(rel(m.atom_number_density(), 7.9623e28) < 1e-4);
        assert_eq!(m.displacement_energy_ev(14), Ok(21.0));
        assert_eq!(m.surface_binding_energy_ev(8), Ok(2.0));
        let s = toml::to_string(&m).unwrap();
        let back: Material = toml::from_str(&s).unwrap();
        assert_eq!(m.components().len(), back.components().len());
        assert!(rel(m.atom_number_density(), back.atom_number_density()) < 1e-14);
        assert_eq!(back.surface_binding_energy_ev(8), Ok(2.0));
        // Deterministic output.
        assert_eq!(s, toml::to_string(&back).unwrap());
        let j = serde_json::to_string(&m).unwrap();
        let _: Material = serde_json::from_str(&j).unwrap();
    }

    #[test]
    fn bad_specs_are_rejected() {
        let mixed = r#"
[[elements]]
symbol = "Si"
atom_fraction = 1.0
[[elements]]
symbol = "O"
mass_fraction = 1.0
"#;
        assert!(toml::from_str::<Material>(mixed).is_err());
        let unknown = "[[elements]]\nsymbol = \"Xx\"\natom_fraction = 1.0\n";
        assert!(toml::from_str::<Material>(unknown).is_err());
        assert!(toml::from_str::<Material>("elements = []").is_err());
    }
}
