//! `lindhard` — clean-room Monte Carlo transport of ions and electrons in matter.
//!
//! The roadmap and module plan are in `docs/architecture.md`. Implemented so
//! far: the shared core (units, constants, the element table, materials, the 1D
//! layered geometry and deterministic random streams), the ion screening
//! potentials and scattering integral, electronic stopping, and the amorphous
//! BCA engine with full recoil cascades (`ion::bca`), damage models
//! (`ion::damage`), and the tallies: moments, Pearson IV / dual-Pearson, range,
//! damage and escape statistics (`tally`). For electrons: the validated
//! data boundary (`electron::data`) and the event loop on cross-section tables
//! (`electron::transport`).

#![forbid(unsafe_code)]

pub mod constants;
pub mod electron;
pub mod elements;
pub mod geometry;
pub mod input;
pub mod ion;
pub mod material;
pub mod rng;
pub mod tally;
pub mod units;

/// Crate version, as stamped by Cargo.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
