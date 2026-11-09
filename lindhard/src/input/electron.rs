//! The electron run description: one TOML document whose `[electron]` table
//! describes the beam, the transport choices, the elastic and inelastic
//! models, the per-material electron data and the tallies, next to the same
//! `[materials]` and `[target]` tables as an ion run.
//!
//! This module only **exposes** what the library already does
//! ([`crate::electron`], [`crate::tally::electron`]); it adds no physics. The
//! schema is documented, with the rules for extending it, in `docs/cli.md`
//! (section "Electron runs"). In short:
//!
//! ```toml
//! [electron.beam]
//! energy_ev = 10000.0
//!
//! [electron.transport]
//! cutoff_ev = 50.0
//!
//! [electron.elastic]
//! potential = "thomas-fermi-yukawa"
//!
//! [electron.materials.Si]
//! optical_elf = "si_elf.toml"
//!
//! [target]
//! substrate = "Si"
//!
//! [run]
//! histories = 1000
//! seed = 1
//! ```
//!
//! [`ElectronInput::resolve_in`] validates the description, reads the data
//! files it names and checks their provenance with the loaders of
//! [`crate::electron::data`] (an optical ELF without a provenance is refused;
//! so are band, phonon and polaron parameters without one), and returns the
//! engine's types ([`ResolvedElectron`]). Building the cross-section tables
//! and running the transport is the caller's (`lindhard-cli`).
//!
//! Like the ion input, every table rejects unknown keys and the echoed input
//! ([`ResolvedElectron::input`]) has every default filled in, so a result can be
//! reproduced from its own header.

mod resolve;
mod resolved;
mod schema;
#[cfg(test)]
mod tests;

pub use resolved::{
    DataFile, ElasticChoice, ResolvedElectron, ResolvedElectronMaterial, ResolvedInnerShells,
};
pub use schema::*;
