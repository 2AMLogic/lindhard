//! Dynamic composition: a finite-slab target whose per-slab composition can
//! change, with volume relaxation.
//!
//! This is the stateful target model for fluence-dependent runs. The
//! [`CompositionGrid`] itself is only the bookkeeping: callers supply explicit
//! per-slab, per-element inventory deltas. Two submodule layers sit on top of
//! it: [`InventoryTally`] turns the transport events of a block of primaries
//! into those deltas (it is a [`BcaTally`](crate::ion::bca::BcaTally)), and
//! [`DynamicRun`] is the fluence stepping loop (fixed or adaptive steps,
//! deterministic at any thread count; see its module docs).
//!
//! # Inventory and units
//!
//! Each slab stores an **areal inventory** `A_i` in atoms/m² per atomic number
//! (kept sorted by `Z`, so every derived quantity is deterministic). Thickness
//! is in metres. A slab built from an existing layer gets
//! `A_i = n_i * thickness`, with `n_i` the layer's atom number density
//! ([`Material::number_density_of`]).
//!
//! # Coordinates
//!
//! As in [`crate::geometry`]: depth `x = 0` is the front surface and grows into
//! the target. After every relaxation the cumulative boundaries are rebuilt
//! from the front surface by [`Stack::new`], so layers are contiguous. The
//! front surface is fixed at `x = 0`: swelling moves the interior interfaces
//! and the back face, never the front surface. With sputter erosion on
//! ([`DynamicConfig::erosion`]) the lost material is removed from the front
//! and the grid is re-anchored, so the new surface is again `x = 0`; the
//! cumulative recession `R` (m, reported per step) maps a depth back to the
//! original frame as `x + R`. With erosion off the surface never moves. Only
//! finite slabs are mutable. An optional semi-infinite substrate is an
//! immutable backing material and is never part of the inventory.
//!
//! # Volume relaxation
//!
//! After an update, each slab's thickness is recomputed from its inventory so
//! atomic densities stay consistent with a stated convention
//! ([`Relaxation`]). The convention is applied to **every** slab on every
//! update, including untouched ones, so a zero update is the identity if and
//! only if the grid is already consistent with the convention.
//!
//! * [`Relaxation::IdealMixing`]: additive atomic volumes (Vegard-type
//!   additivity applied to atoms). With reference atomic volumes `v_i`
//!   (m³/atom), `thickness = sum_i A_i v_i`. Atom fractions are
//!   `A_i / sum_j A_j` and the mass density is the total mass per area divided
//!   by the thickness. Choose `v_i` as the atomic volume in the phase you wish
//!   to represent; for an element in its own solid,
//!   [`atomic_volume_from_density`] gives `M / (N_A rho)`. The caller must
//!   supply the volume of every species that occurs: none is inferred, in
//!   particular not for gases with no elemental solid density. A compound's
//!   real density generally differs from the ideal-mixture prediction, and so
//!   does the first relaxation of a grid built from a layer of such a compound;
//!   the grid's initial thicknesses are kept as given until the first update.
//!   Tradeoff: simple and parameter-light, but ignores chemistry, voids and
//!   amorphization swelling.
//! * [`Relaxation::FixedNumberDensity`]: the caller specifies one total atom
//!   number density `n_mix` (atoms/m³) for every slab, and
//!   `thickness = sum_i A_i / n_mix`. Needs no per-species data and reproduces
//!   a known compound density, but the density cannot respond to composition:
//!   a slab that goes from Si to mostly implanted species keeps the same
//!   atoms per volume, and a single value serves all slabs, which suits
//!   compositions that stay near one phase.
//!
//! Neither option is an equation of state with a literature-backed pressure or
//! phase model; both are stated conventions.
//!
//! # Conservation and update semantics
//!
//! A slab inventory changes only by the deltas passed to
//! [`CompositionGrid::apply`]. How a caller derives them:
//!
//! * A retained primary adds one atom of its species at its rest slab. A
//!   primary that is transmitted or backscattered without retention removes no
//!   target atom and adds nothing.
//! * A recoil moves an existing atom: subtract at recoil creation, add where it
//!   comes to rest. An internal relocation therefore cancels in the global
//!   total of that species ([`CompositionGrid::total_inventory`]); a recoil
//!   that escapes the target is a loss, with no addition after escape.
//! * Beam and target atoms of the same element are the same species in the
//!   inventory; the roles differ only in how the deltas are produced.
//! * Aggregate front/back escape counts cannot say which slab an escaped atom
//!   came from; loss deltas need the origin slab from the transport events,
//!   not the exit position. The default tally gets it by subtracting at the
//!   recoil event; [`InventoryTally`] also counts sputtered atoms per origin
//!   slab (`Particle::origin_layer`).
//! * With erosion on, the loss of a sputtered atom is not taken at its origin
//!   slab but from the front of the target, element by element, the origin
//!   loss being cancelled so nothing is removed twice (see
//!   [`DynamicRun`], "Erosion").
//!
//! # Transactions
//!
//! [`CompositionGrid::apply`] either succeeds completely or returns an error
//! and leaves the grid untouched. It rejects non-finite deltas, unknown slabs
//! or elements, results below zero, missing or invalid reference volumes and
//! overflow to non-finite thickness or density. A result a few ulps below zero
//! (within `4 eps` of the larger of the inventory and the delta) from removing
//! exactly what was there is rounded to zero. A species at exactly zero is
//! dropped from its slab; a slab with no species left is removed, and its
//! index is reported in [`UpdateOutcome::removed_slabs`] (indices of the grid
//! before the update). If every slab disappears the update is an error unless
//! a substrate backs the grid, in which case the grid is the bare substrate.
//!
//! # Energies
//!
//! Displacement, lattice-binding and surface-binding energies given on the
//! initial materials are remembered per slab and per element, including
//! explicitly overridden values, and re-applied when materials are rebuilt.
//! An element new to a slab gets the defaults of [`crate::elements`].

mod adapter;
mod driver;

pub use adapter::{InventoryTally, Yields};
pub use driver::{DynamicConfig, DynamicRun, DynamicRunError, StepPolicy, StepRecord};

use std::collections::BTreeMap;

use crate::constants::AVOGADRO;
use crate::elements::element;
use crate::geometry::{GeometryError, Stack};
use crate::material::{Material, MaterialError};

/// Errors from building or updating a [`CompositionGrid`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum DynamicError {
    /// A delta was NaN or infinite.
    #[error("non-finite inventory delta {delta} atoms/m^2 for slab {slab}, Z={z}")]
    NonFiniteDelta {
        /// Slab index.
        slab: usize,
        /// Atomic number.
        z: u8,
        /// Offending value.
        delta: f64,
    },
    /// A delta referred to a slab that does not exist.
    #[error("slab {slab} out of range ({n_slabs} slabs)")]
    SlabOutOfRange {
        /// Slab index.
        slab: usize,
        /// Number of slabs.
        n_slabs: usize,
    },
    /// Unknown atomic number.
    #[error("unknown atomic number {0}")]
    UnknownAtomicNumber(u8),
    /// An update would make an inventory negative.
    #[error("slab {slab}, Z={z}: removing more than present ({result} atoms/m^2 would remain)")]
    NegativeInventory {
        /// Slab index.
        slab: usize,
        /// Atomic number.
        z: u8,
        /// Resulting (negative) value.
        result: f64,
    },
    /// The relaxation convention has no volume for a species present.
    #[error("no reference atomic volume for Z={0}")]
    MissingVolume(u8),
    /// A reference volume or number density was not finite and positive.
    #[error("invalid relaxation parameter {0}: must be finite and positive")]
    InvalidParameter(f64),
    /// An element has no tabulated solid density to derive a volume from.
    #[error("element Z={0} has no elemental solid density; supply its volume explicitly")]
    NoElementalDensity(u8),
    /// A relaxed thickness or density was not finite and positive.
    #[error("slab {slab}: relaxed thickness or density is not finite and positive (overflow)")]
    Overflow {
        /// Slab index.
        slab: usize,
    },
    /// Every slab was depleted and there is no substrate.
    #[error("update would leave an empty target")]
    EmptyTarget,
    /// A slab of the initial stack had no usable composition.
    #[error("initial composition invalid: {0}")]
    InvalidInitial(String),
    /// Stack construction failed.
    #[error(transparent)]
    Geometry(#[from] GeometryError),
    /// Material construction failed.
    #[error(transparent)]
    Material(#[from] MaterialError),
}

/// Atomic volume `M / (N_A rho)` in m³/atom of element `z` in its elemental
/// solid, from the default density of [`crate::elements`] (conventional
/// atomic weight, CODATA `N_A`). Errors if `z` is unknown or has no tabulated
/// density (for example a gas): no volume is invented for those.
pub fn atomic_volume_from_density(z: u8) -> Result<f64, DynamicError> {
    let el = element(z).ok_or(DynamicError::UnknownAtomicNumber(z))?;
    let rho = el
        .density_kg_m3()
        .ok_or(DynamicError::NoElementalDensity(z))?;
    Ok(el.atomic_weight * 1.0e-3 / (AVOGADRO * rho))
}

/// Volume-relaxation convention. See the [module docs](self).
#[derive(Debug, Clone, PartialEq)]
pub enum Relaxation {
    /// Ideal mixing of atomic volumes, `thickness = sum A_i v_i`.
    IdealMixing {
        /// Reference atomic volume per atomic number, m³/atom.
        atomic_volumes_m3: BTreeMap<u8, f64>,
    },
    /// Caller-specified total atom number density, `thickness = sum A_i / n`.
    FixedNumberDensity {
        /// Atoms per m³.
        number_density_m3: f64,
    },
}

impl Relaxation {
    /// Ideal mixing with explicit volumes (m³/atom) per atomic number.
    pub fn ideal_mixing(volumes: &[(u8, f64)]) -> Result<Self, DynamicError> {
        let mut map = BTreeMap::new();
        for &(z, v) in volumes {
            if element(z).is_none() {
                return Err(DynamicError::UnknownAtomicNumber(z));
            }
            if !(v.is_finite() && v > 0.0) {
                return Err(DynamicError::InvalidParameter(v));
            }
            map.insert(z, v);
        }
        Ok(Self::IdealMixing {
            atomic_volumes_m3: map,
        })
    }

    /// Ideal mixing with each listed element at its elemental solid volume
    /// ([`atomic_volume_from_density`]).
    pub fn ideal_mixing_elemental(zs: &[u8]) -> Result<Self, DynamicError> {
        let v: Result<Vec<_>, _> = zs
            .iter()
            .map(|&z| atomic_volume_from_density(z).map(|v| (z, v)))
            .collect();
        Self::ideal_mixing(&v?)
    }

    /// Fixed total number density, atoms/m³.
    pub fn fixed_number_density(number_density_m3: f64) -> Result<Self, DynamicError> {
        if !(number_density_m3.is_finite() && number_density_m3 > 0.0) {
            return Err(DynamicError::InvalidParameter(number_density_m3));
        }
        Ok(Self::FixedNumberDensity { number_density_m3 })
    }

    /// Thickness in m of an inventory (all entries positive).
    fn thickness_m(&self, inv: &BTreeMap<u8, f64>) -> Result<f64, DynamicError> {
        match self {
            Self::IdealMixing { atomic_volumes_m3 } => {
                let mut t = 0.0;
                for (&z, &a) in inv {
                    let v = atomic_volumes_m3
                        .get(&z)
                        .ok_or(DynamicError::MissingVolume(z))?;
                    t += a * v;
                }
                Ok(t)
            }
            Self::FixedNumberDensity { number_density_m3 } => {
                Ok(inv.values().sum::<f64>() / number_density_m3)
            }
        }
    }
}

/// An inventory change for one element in one slab, atoms/m² (positive adds).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InventoryDelta {
    /// Slab index (0 = front), in the grid before the update.
    pub slab: usize,
    /// Atomic number.
    pub z: u8,
    /// Change in areal inventory, atoms/m².
    pub delta_atoms_m2: f64,
}

/// Result of a successful [`CompositionGrid::apply`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UpdateOutcome {
    /// Indices (before the update) of slabs removed because they emptied.
    pub removed_slabs: Vec<usize>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct Energies {
    e_d: Option<f64>,
    e_b: Option<f64>,
    e_s: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
struct Slab {
    inventory: BTreeMap<u8, f64>,
    thickness_m: f64,
    energies: BTreeMap<u8, Energies>,
    name: Option<String>,
}

impl Slab {
    fn from_layer(material: &Material, thickness_m: f64) -> Result<Self, DynamicError> {
        let mut inventory = BTreeMap::new();
        let mut energies = BTreeMap::new();
        for c in material.components() {
            let z = c.z();
            let a = material.number_density_of(z)? * thickness_m;
            if a > 0.0 {
                inventory.insert(z, a);
            }
            energies.insert(
                z,
                Energies {
                    e_d: material.displacement_energy_ev(z).ok(),
                    e_b: material.lattice_binding_energy_ev(z).ok(),
                    e_s: material.surface_binding_energy_ev(z).ok(),
                },
            );
        }
        if inventory.is_empty() {
            return Err(DynamicError::InvalidInitial(
                "a layer has no atoms".to_string(),
            ));
        }
        Ok(Self {
            inventory,
            thickness_m,
            energies,
            name: material.name().map(str::to_owned),
        })
    }

    /// Material for this slab's current inventory and thickness.
    fn material(&self) -> Result<Material, DynamicError> {
        let parts: Vec<(u8, f64)> = self.inventory.iter().map(|(&z, &a)| (z, a)).collect();
        let atoms_per_area: f64 = self.inventory.values().sum();
        let mass_per_area: f64 = self
            .inventory
            .iter()
            .map(|(&z, &a)| a * element(z).expect("validated").atomic_weight * 1.0e-3 / AVOGADRO)
            .sum();
        let density = mass_per_area / self.thickness_m;
        if !(atoms_per_area.is_finite() && density.is_finite() && density > 0.0) {
            return Err(DynamicError::Overflow { slab: usize::MAX });
        }
        let mut m = Material::from_atom_fractions(&parts, Some(density))?;
        if let Some(n) = &self.name {
            m = m.with_name(n.clone());
        }
        for &(z, _) in &parts {
            if let Some(e) = self.energies.get(&z) {
                if let Some(v) = e.e_d {
                    m.set_displacement_energy_ev(z, v)?;
                }
                if let Some(v) = e.e_b {
                    m.set_lattice_binding_energy_ev(z, v)?;
                }
                if let Some(v) = e.e_s {
                    m.set_surface_binding_energy_ev(z, v)?;
                }
            }
        }
        Ok(m)
    }
}

/// A finite-slab target with per-slab species inventories, an optional
/// immutable substrate, and a [`Relaxation`] convention.
#[derive(Debug, Clone, PartialEq)]
pub struct CompositionGrid {
    slabs: Vec<Slab>,
    substrate: Option<Material>,
    relaxation: Relaxation,
}

impl CompositionGrid {
    /// Build from a [`Stack`]: finite layers become slabs (inventory
    /// `n_i * thickness`, thicknesses kept as given), and a semi-infinite last
    /// layer becomes the immutable substrate.
    pub fn from_stack(stack: &Stack, relaxation: Relaxation) -> Result<Self, DynamicError> {
        let mut slabs = Vec::new();
        let mut substrate = None;
        for layer in stack.layers() {
            if layer.thickness_m().is_finite() {
                slabs.push(Slab::from_layer(layer.material(), layer.thickness_m())?);
            } else {
                substrate = Some(layer.material().clone());
            }
        }
        Ok(Self {
            slabs,
            substrate,
            relaxation,
        })
    }

    /// Number of mutable slabs.
    pub fn n_slabs(&self) -> usize {
        self.slabs.len()
    }

    /// The relaxation convention in use.
    pub fn relaxation(&self) -> &Relaxation {
        &self.relaxation
    }

    /// Slab thicknesses in m, front first.
    pub fn thicknesses_m(&self) -> Vec<f64> {
        self.slabs.iter().map(|s| s.thickness_m).collect()
    }

    /// Inventory of one slab as `(Z, atoms/m^2)` sorted by `Z`.
    pub fn inventory(&self, slab: usize) -> Option<Vec<(u8, f64)>> {
        self.slabs
            .get(slab)
            .map(|s| s.inventory.iter().map(|(&z, &a)| (z, a)).collect())
    }

    /// Record energies for element `z` in every slab, for the case that `z` is
    /// (or later becomes) present: only the values that are still unset are
    /// filled, so explicit values and defaults already in place are kept.
    /// This is how an element that only enters through implantation or
    /// recoil mixing (the beam species, a substrate element) gets the
    /// displacement, lattice-binding and surface-binding energies the engine
    /// requires of every element in a layer.
    pub fn seed_energies(
        &mut self,
        z: u8,
        e_d_ev: Option<f64>,
        e_b_ev: Option<f64>,
        e_s_ev: Option<f64>,
    ) {
        let defaults = Material::from_atom_fractions(&[(z, 1.0)], None).ok();
        for s in &mut self.slabs {
            let e = s.energies.entry(z).or_default();
            e.e_d = e.e_d.or(e_d_ev).or_else(|| {
                defaults
                    .as_ref()
                    .and_then(|m| m.displacement_energy_ev(z).ok())
            });
            e.e_b = e.e_b.or(e_b_ev).or_else(|| {
                defaults
                    .as_ref()
                    .and_then(|m| m.lattice_binding_energy_ev(z).ok())
            });
            e.e_s = e.e_s.or(e_s_ev).or_else(|| {
                defaults
                    .as_ref()
                    .and_then(|m| m.surface_binding_energy_ev(z).ok())
            });
        }
    }

    /// Total areal inventory of element `z` over all finite slabs, atoms/m².
    /// The substrate is not counted. Summed front to back.
    pub fn total_inventory(&self, z: u8) -> f64 {
        self.slabs
            .iter()
            .map(|s| s.inventory.get(&z).copied().unwrap_or(0.0))
            .sum()
    }

    /// Re-relax without changing any inventory; same as `apply(&[])`.
    pub fn relax(&mut self) -> Result<UpdateOutcome, DynamicError> {
        self.apply(&[])
    }

    /// Apply inventory deltas, then relax every slab. Transactional: on error
    /// the grid is unchanged. See the [module docs](self).
    pub fn apply(&mut self, deltas: &[InventoryDelta]) -> Result<UpdateOutcome, DynamicError> {
        let n = self.slabs.len();
        let mut inv: Vec<BTreeMap<u8, f64>> =
            self.slabs.iter().map(|s| s.inventory.clone()).collect();
        // Largest magnitude seen per (slab, Z), for the rounding tolerance.
        let mut scale: BTreeMap<(usize, u8), f64> = BTreeMap::new();
        for d in deltas {
            if d.slab >= n {
                return Err(DynamicError::SlabOutOfRange {
                    slab: d.slab,
                    n_slabs: n,
                });
            }
            if element(d.z).is_none() {
                return Err(DynamicError::UnknownAtomicNumber(d.z));
            }
            if !d.delta_atoms_m2.is_finite() {
                return Err(DynamicError::NonFiniteDelta {
                    slab: d.slab,
                    z: d.z,
                    delta: d.delta_atoms_m2,
                });
            }
            let e = inv[d.slab].entry(d.z).or_insert(0.0);
            let s = scale.entry((d.slab, d.z)).or_insert(e.abs());
            *s = s.max(e.abs()).max(d.delta_atoms_m2.abs());
            *e += d.delta_atoms_m2;
            if !e.is_finite() {
                return Err(DynamicError::Overflow { slab: d.slab });
            }
        }
        // Clean up: tolerance, negatives, zeros.
        for (i, m) in inv.iter_mut().enumerate() {
            for (&z, a) in m.iter_mut() {
                if *a < 0.0 {
                    let tol = 4.0 * f64::EPSILON * scale.get(&(i, z)).copied().unwrap_or(0.0);
                    if -*a <= tol {
                        *a = 0.0;
                    } else {
                        return Err(DynamicError::NegativeInventory {
                            slab: i,
                            z,
                            result: *a,
                        });
                    }
                }
            }
            m.retain(|_, a| *a > 0.0);
        }
        let removed_slabs: Vec<usize> = (0..n).filter(|&i| inv[i].is_empty()).collect();
        if removed_slabs.len() == n && self.substrate.is_none() {
            return Err(DynamicError::EmptyTarget);
        }
        // Relax and rebuild into a new slab list.
        let mut new_slabs = Vec::with_capacity(n);
        for (i, m) in inv.into_iter().enumerate() {
            if m.is_empty() {
                continue;
            }
            let t = self.relaxation.thickness_m(&m)?;
            if !(t.is_finite() && t > 0.0) {
                return Err(DynamicError::Overflow { slab: i });
            }
            let energies = self.slabs[i].energies.clone();
            let slab = Slab {
                inventory: m,
                thickness_m: t,
                energies,
                name: self.slabs[i].name.clone(),
            };
            // Validate materials now so failure cannot leave a partial grid.
            slab.material().map_err(|e| match e {
                DynamicError::Overflow { .. } => DynamicError::Overflow { slab: i },
                other => other,
            })?;
            new_slabs.push(slab);
        }
        self.slabs = new_slabs;
        Ok(UpdateOutcome { removed_slabs })
    }

    /// Convert to the engine's layer model: contiguous finite layers from the
    /// front surface, then the substrate if any.
    pub fn to_stack(&self) -> Result<Stack, DynamicError> {
        let layers = self
            .slabs
            .iter()
            .map(|s| Ok((s.material()?, s.thickness_m)))
            .collect::<Result<Vec<_>, DynamicError>>()?;
        Ok(Stack::new(layers, self.substrate.clone())?)
    }
}
