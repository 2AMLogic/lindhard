//! Assembling a resolved electron input into cross-section tables and a
//! transport engine: the one implementation every caller uses.
//!
//! [`ElectronInput::resolve_in`](super::ElectronInput::resolve_in) validates a
//! run description and returns [`ResolvedElectron`]; this module turns that
//! into the elastic and inelastic tables of every material, the optional
//! inner-shell channels, and the [`Transport`] built on them. Nothing here is
//! physics: every model is [`crate::electron`]'s, chosen and parameterised by
//! the input.
//!
//! It lives in the library rather than in `lindhard-cli` so that the library's
//! own benches, examples and tests build the same tables `lindhard run` builds
//! (#750). Before, `lindhard-cli/src/electron.rs` held the only copy and
//! `lindhard/benches/electron_common/mod.rs` held a hand copy that had drifted
//! from it: the bench copy set no rate refinement (#339), fitted Mermin
//! oscillators with the default options instead of the material's (#306),
//! built no inner-shell channels, and ignored the input's elastic corrections
//! on the DHFS potential.
//!
//! What stays in the CLI is what is CLI-specific: the table cache, the run
//! timings and the `anyhow` context strings.
//!
//! # Example
//!
//! ```no_run
//! use lindhard::input::electron::assemble;
//! use lindhard::input::electron::ElectronInput;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let input = ElectronInput::from_toml_str(&std::fs::read_to_string("run.toml")?)?;
//! let resolved = input.resolve_in(std::path::Path::new("."))?;
//! let tables = assemble::material_tables(&resolved)?;
//! let transport = assemble::transport(&resolved, &tables)?;
//! # let _ = transport;
//! # Ok(())
//! # }
//! ```

use crate::electron::data::{CrossSectionTable, ElectronDataError, ShellChannelTable};
use crate::electron::elastic::table::{
    combine, default_probability_grid, AtomicElastic, ElasticTableError, PotentialSource,
    SalvatDhfsTable, ThomasFermiYukawa, DEFAULT_REFINE_TOLERANCE,
};
use crate::electron::elastic::{ElasticError, SalvatDhfs, SolverOptions};
use crate::electron::inelastic::table::{
    build_inelastic_table_for_model, build_shell_channel_tables, EnergyAxis, InelasticTableError,
    InelasticTableOptions, RateRefinement,
};
use crate::electron::inelastic::{
    hydrogenic_shell_elfs, PennAlgorithm, PennInelastic, ShellElfGrid, ShellResolvedChannels,
    SinglePolePenn,
};
use crate::electron::transport::{LayerTables, Transport, TransportError};
use crate::material::{Material, MaterialError};

use super::{PotentialChoice, ResolvedElectron, ResolvedElectronMaterial};

/// Errors from assembling the tables or the transport of a resolved electron
/// input.
///
/// The model and table builders have their own error types; this wraps them,
/// adding the material name and the layer index where a message would
/// otherwise not say which part of the input failed.
#[derive(Debug, thiserror::Error)]
pub enum ElectronAssemblyError {
    /// No screened potential for an element of a material.
    #[error(transparent)]
    Elastic(#[from] ElasticError),
    /// The elastic table builder failed.
    #[error(transparent)]
    ElasticTable(#[from] ElasticTableError),
    /// A model or a table failed its validation.
    #[error(transparent)]
    Data(#[from] ElectronDataError),
    /// The inelastic table builder failed.
    #[error(transparent)]
    InelasticTable(#[from] InelasticTableError),
    /// A material lookup failed (an element's number density).
    #[error(transparent)]
    Material(#[from] MaterialError),
    /// The stack, the tables and the configuration are not a valid transport.
    #[error("setting up the electron transport: {source}")]
    Transport {
        /// The engine's error.
        #[source]
        source: TransportError,
    },
    /// Attaching the insulator channels of a layer failed.
    #[error("electron.materials.{material}.channels: layer {layer}: {source}")]
    InsulatorChannels {
        /// The material of the layer.
        material: String,
        /// Index of the layer in the stack.
        layer: usize,
        /// The engine's error.
        #[source]
        source: TransportError,
    },
    /// Attaching the inner-shell channels of a layer failed.
    #[error("electron.materials.{material}.inner_shells: layer {layer}: {source}")]
    InnerShells {
        /// The material of the layer.
        material: String,
        /// Index of the layer in the stack.
        layer: usize,
        /// The engine's error.
        #[source]
        source: TransportError,
    },
    /// The number of table pairs or inner-shell channel lists differs from the
    /// number of materials of the input.
    #[error("the input has {materials} materials but {given} {what} were given")]
    MaterialCount {
        /// Materials of the resolved input.
        materials: usize,
        /// What was given too few or too many of.
        what: &'static str,
        /// How many were given.
        given: usize,
    },
}

/// The per-material tables of an input, in [`ResolvedElectron::materials`]
/// order: the two cross-section tables and the inner-shell channels (empty
/// without `inner_shells`) of every material.
#[derive(Debug, Clone)]
pub struct MaterialTables {
    /// The elastic and inelastic table of each material.
    pub tables: Vec<LayerTables>,
    /// The inner-shell channel tables of each material.
    pub shells: Vec<Vec<ShellChannelTable>>,
}

impl MaterialTables {
    /// Pair up already-built tables and channels, checking both against the
    /// material count of `r`.
    ///
    /// [`material_tables`] builds them; this is for a caller that builds or
    /// caches the tables itself (`lindhard-cli`'s table cache) and still wants
    /// [`transport`] to assemble them.
    pub fn new(
        r: &ResolvedElectron,
        tables: Vec<LayerTables>,
        shells: Vec<Vec<ShellChannelTable>>,
    ) -> Result<Self, ElectronAssemblyError> {
        let materials = r.materials.len();
        if tables.len() != materials {
            return Err(ElectronAssemblyError::MaterialCount {
                materials,
                what: "table pairs",
                given: tables.len(),
            });
        }
        if shells.len() != materials {
            return Err(ElectronAssemblyError::MaterialCount {
                materials,
                what: "inner-shell channel lists",
                given: shells.len(),
            });
        }
        Ok(Self { tables, shells })
    }
}

/// The elastic table of one material.
///
/// Each element's potential doubles as the electron density the optional
/// exchange and correlation-polarization corrections need
/// ([`AtomicElastic::compute_corrected`]): the Yukawa stand-in's own Poisson
/// density, or the DHFS Poisson density of Salvat et al. (1987) Eq. (12)
/// ([`crate::electron::elastic::corrections`]). Without corrections this is
/// the plain partial-wave table of
/// [`build_elastic_table`](crate::electron::elastic::table::build_elastic_table).
pub fn elastic_table(
    r: &ResolvedElectron,
    m: &Material,
) -> Result<CrossSectionTable, ElectronAssemblyError> {
    let grid = &r.table_energy_ev;
    let solver = SolverOptions::default();
    let desc = match r.elastic.potential {
        PotentialChoice::SalvatDhfs => SalvatDhfsTable.description(),
        PotentialChoice::ThomasFermiYukawa => ThomasFermiYukawa.description(),
    };
    let desc = format!("{desc}{}", r.elastic.corrections_description());
    let mut atoms = Vec::new();
    for c in m.components() {
        if c.atom_fraction() <= 0.0 {
            continue;
        }
        let z = c.z();
        let corrections = r.elastic.corrections(z);
        atoms.push(match r.elastic.potential {
            PotentialChoice::SalvatDhfs => {
                let p = SalvatDhfs::for_element(u32::from(z))?;
                AtomicElastic::compute_corrected(z, &p, &p, &desc, grid, &corrections, solver)?
            }
            PotentialChoice::ThomasFermiYukawa => {
                let y = ThomasFermiYukawa::yukawa(z)?;
                AtomicElastic::compute_corrected(z, &y, &y, &desc, grid, &corrections, solver)?
            }
        });
    }
    Ok(combine(
        m,
        &atoms,
        &default_probability_grid(),
        Some(DEFAULT_REFINE_TOLERANCE),
    )?)
}

/// The inelastic table of one material, on the axis and with the Fermi
/// energy of [`inelastic_axis`] and the rate refinement of
/// [`rate_refinement`].
pub fn inelastic_table(
    r: &ResolvedElectron,
    m: &ResolvedElectronMaterial,
) -> Result<CrossSectionTable, ElectronAssemblyError> {
    let (fermi_ev, axis) = inelastic_axis(r, m);
    // The Mermin oscillators are fitted with the material's options
    // (`[electron.materials.<name>.mermin_fit]`, default
    // `MerminFitOptions::default()`, #306).
    let model = match r.inelastic {
        PennAlgorithm::Mermin => PennInelastic::mermin(m.optical_elf.clone(), &m.mermin_fit)?,
        other => PennInelastic::try_new(other, m.optical_elf.clone())?,
    }
    .with_fermi_energy_ev(fermi_ev)?;
    let mut options = InelasticTableOptions::new(r.table_energy_ev.clone()).with_axis(axis);
    options.rate_refinement = rate_refinement(r, m);
    Ok(build_inelastic_table_for_model(
        &model,
        &m.material,
        &options,
    )?)
}

/// The inner-shell channel tables of one material (`inner_shells` given):
/// one hydrogenic channel per requested subshell of every element, each
/// tabulated with the valence channel's single-pole model on the axis of
/// the inelastic table. The valence table of the result is not used: the
/// run keeps the material's own inelastic table, which is the same build.
///
/// Empty for a material without `inner_shells`.
pub fn shell_tables(
    r: &ResolvedElectron,
    m: &ResolvedElectronMaterial,
) -> Result<Vec<ShellChannelTable>, ElectronAssemblyError> {
    let Some(i) = &m.inner_shells else {
        return Ok(Vec::new());
    };
    let (fermi_ev, axis) = inelastic_axis(r, m);
    let valence = SinglePolePenn::new(m.optical_elf.clone()).with_fermi_energy_ev(fermi_ev)?;
    let mut wanted = Vec::new();
    for c in m.material.components() {
        if c.atom_fraction() <= 0.0 {
            continue;
        }
        let n = m.material.number_density_of(c.z())?;
        for &s in &i.subshells {
            wanted.push((c.z(), s, n));
        }
    }
    let shells = hydrogenic_shell_elfs(
        m.optical_elf.material(),
        &i.binding,
        &wanted,
        ShellElfGrid::default(),
    )?;
    let channels = ShellResolvedChannels::new(valence, shells)?;
    let mut options = InelasticTableOptions::new(r.table_energy_ev.clone()).with_axis(axis);
    options.rate_refinement = rate_refinement(r, m);
    let built =
        build_shell_channel_tables(&channels, &m.material, &options, i.binding.provenance())?;
    Ok(built.shells)
}

/// The rate refinement of the inelastic table of `m` (#339): the library's
/// default ([`RateRefinement::default`]) on the band-bottom axis of a
/// material with a band, where the rate climbs from zero at the Fermi
/// level; `None` on the model's own axis, whose table keeps exactly the
/// rows of `electron.tables`. The refined table has the grid's rows and the
/// rows the refinement adds between 2 and 20 eV above the Fermi level
/// ([`crate::electron::inelastic::table`], "Rate refinement"); the elastic
/// table keeps the grid. Part of `lindhard-cli`'s table cache key.
pub fn rate_refinement(
    r: &ResolvedElectron,
    m: &ResolvedElectronMaterial,
) -> Option<RateRefinement> {
    (inelastic_axis(r, m).1 == EnergyAxis::BandBottom).then(RateRefinement::default)
}

/// The model's Fermi energy and the table axis of the inelastic table of `m`.
///
/// For a material with a band, the table is built on the band-bottom axis
/// the transport reads it on, with the model's Fermi energy set to the
/// band's minimum excitation energy
/// ([`crate::electron::boundary::BandStructure::min_excitation_ev`]: the
/// Fermi energy of a metal, the conduction-band bottom of an insulator):
/// the row at `E` is the model's at `T = E - E_F`, kinematics on `E`, losses
/// below `E - E_F`. This is the convention of cstool's
/// `compile_full_imfp_icdf` (Nebula-simulator/cstool commit `0c739eb`,
/// `cstool/dielectric_function/compile.py`; module docs of
/// [`crate::electron::inelastic::table`], "Energy axis"). Without a band
/// the table is on the model's own axis with
/// `[electron.inelastic] fermi_energy_ev`, as before (#241).
pub fn inelastic_axis(r: &ResolvedElectron, m: &ResolvedElectronMaterial) -> (f64, EnergyAxis) {
    match &m.band {
        Some(b) => (b.min_excitation_ev(), EnergyAxis::BandBottom),
        None => (r.inelastic_fermi_ev, EnergyAxis::ModelFermiLevel),
    }
}

/// Both cross-section tables and the inner-shell channels of every material
/// of `r`, in [`ResolvedElectron::materials`] order.
///
/// Run this inside a [`rayon::ThreadPool::install`] to choose the threads the
/// table build uses; the result does not depend on the thread count.
pub fn material_tables(r: &ResolvedElectron) -> Result<MaterialTables, ElectronAssemblyError> {
    let tables = r
        .materials
        .iter()
        .map(|m| {
            Ok(LayerTables {
                elastic: elastic_table(r, &m.material)?,
                inelastic: inelastic_table(r, m)?,
            })
        })
        .collect::<Result<Vec<_>, ElectronAssemblyError>>()?;
    let shells = r
        .materials
        .iter()
        .map(|m| shell_tables(r, m))
        .collect::<Result<Vec<_>, ElectronAssemblyError>>()?;
    MaterialTables::new(r, tables, shells)
}

/// The transport engine of `r` on `tables`, with the band parameters, the
/// insulator channels and the inner-shell channels the input asks for.
///
/// `tables` is indexed by material; the engine is built per layer through
/// [`ResolvedElectron::layer_material`].
pub fn transport(
    r: &ResolvedElectron,
    tables: &MaterialTables,
) -> Result<Transport, ElectronAssemblyError> {
    if tables.tables.len() != r.materials.len() {
        return Err(ElectronAssemblyError::MaterialCount {
            materials: r.materials.len(),
            what: "table pairs",
            given: tables.tables.len(),
        });
    }
    if tables.shells.len() != r.materials.len() {
        return Err(ElectronAssemblyError::MaterialCount {
            materials: r.materials.len(),
            what: "inner-shell channel lists",
            given: tables.shells.len(),
        });
    }
    let layer_tables: Vec<LayerTables> = r
        .layer_material
        .iter()
        .map(|&i| tables.tables[i].clone())
        .collect();
    let bands: Option<Vec<_>> = r
        .layer_material
        .iter()
        .map(|&i| r.materials[i].band.clone())
        .collect();
    let mut t = match bands {
        Some(b) => Transport::with_band_structures(r.stack.clone(), layer_tables, b, r.config),
        None => Transport::new(r.stack.clone(), layer_tables, r.config),
    }
    .map_err(|source| ElectronAssemblyError::Transport { source })?;
    for (layer, &i) in r.layer_material.iter().enumerate() {
        let c = &r.materials[i].channels;
        if c.any() {
            t = t
                .with_insulator_channels(layer, c.clone())
                .map_err(|source| ElectronAssemblyError::InsulatorChannels {
                    material: r.materials[i].name.clone(),
                    layer,
                    source,
                })?;
        }
        if !tables.shells[i].is_empty() {
            t = t
                .with_inner_shells(layer, tables.shells[i].clone())
                .map_err(|source| ElectronAssemblyError::InnerShells {
                    material: r.materials[i].name.clone(),
                    layer,
                    source,
                })?;
        }
    }
    Ok(t)
}

#[cfg(test)]
mod elastic_tests {
    use super::*;
    use crate::electron::data::CrossSectionTable;
    use crate::input::electron::ElectronInput;
    use std::path::Path;

    /// The CLI example on a two-point grid with the given potential and
    /// exchange setting.
    fn resolved(potential: PotentialChoice, exchange: bool) -> ResolvedElectron {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../examples/electron");
        let text = std::fs::read_to_string(dir.join("e_10keV_si.toml")).unwrap();
        let mut r = ElectronInput::from_toml_str(&text)
            .unwrap()
            .resolve_in(&dir)
            .unwrap();
        r.elastic.potential = potential;
        r.elastic.exchange = exchange;
        r.table_energy_ev = vec![100.0, 1000.0];
        r
    }

    fn table(potential: PotentialChoice, exchange: bool) -> CrossSectionTable {
        let r = resolved(potential, exchange);
        elastic_table(&r, &r.materials[0].material).unwrap()
    }

    /// The corrections reach the DHFS table as they reach the stand-in's
    /// (#169: they were dropped for `salvat-dhfs`).
    #[test]
    fn exchange_applies_to_both_potentials() {
        for p in [
            PotentialChoice::SalvatDhfs,
            PotentialChoice::ThomasFermiYukawa,
        ] {
            let (off, on) = (table(p, false), table(p, true));
            assert!(!off.model().contains("exchange"), "{p:?}: {}", off.model());
            assert!(on.model().contains("exchange"), "{p:?}: {}", on.model());
            assert_ne!(
                off.inverse_mfp_per_m(),
                on.inverse_mfp_per_m(),
                "{p:?}: exchange left the table unchanged"
            );
        }
        let dhfs = table(PotentialChoice::SalvatDhfs, false);
        assert!(dhfs.model().contains("Salvat"), "{}", dhfs.model());
    }
}

#[cfg(test)]
mod mermin_fit_tests {
    use super::*;
    use crate::input::electron::ElectronInput;
    use std::path::Path;

    /// A Cu input with `oscillators = 6` (#306) builds its inelastic table
    /// (on a two-energy grid) with a six-oscillator fit, recorded in the
    /// table's `model`.
    #[test]
    fn cu_six_oscillators_builds_its_table() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../validation/data/optical");
        if !dir.join("cu_elf_hagemann1975.toml").is_file() {
            eprintln!("{} not found; skipped", dir.display());
            return;
        }
        let text = "[electron.beam]\nenergy_ev = 1000.0\n[electron.transport]\ncutoff_ev = 20.0\n\
                    [electron.elastic]\npotential = \"thomas-fermi-yukawa\"\n\
                    [electron.inelastic]\nmodel = \"mermin-melf\"\n\
                    [electron.materials.Cu]\noptical_elf = \"cu_elf_hagemann1975.toml\"\n\
                    [electron.materials.Cu.mermin_fit]\noscillators = 6\n\
                    [target]\nsubstrate = \"Cu\"\n[run]\nhistories = 1\nseed = 1\n";
        let mut r = ElectronInput::from_toml_str(text)
            .unwrap()
            .resolve_in(&dir)
            .unwrap();
        r.table_energy_ev = vec![100.0, 150.0];
        let t = inelastic_table(&r, &r.materials[0]).unwrap();
        assert!(t.model().contains("(6 oscillators "), "{}", t.model());
        assert_eq!(t.energy_ev(), r.table_energy_ev.as_slice());
        assert!(t
            .inverse_mfp_per_m()
            .iter()
            .all(|&x| x.is_finite() && x > 0.0));
    }
}
