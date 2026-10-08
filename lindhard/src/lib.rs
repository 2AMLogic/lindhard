//! `lindhard` — clean-room Monte Carlo transport of ions and electrons in matter.
//!
//! Roadmap, milestone status and the module plan are in `docs/architecture.md`.
//! Each module's own documentation is the maintained description of what it
//! provides:
//!
//! * [`constants`]: physical constants in SI units, each with its source.
//! * [`electron`]: electron transport, from validated data to tallies.
//! * [`elements`]: the per-element data table, Z = 1..=92.
//! * [`geometry`]: the `Geometry` trait and the layered, voxel and mesh targets.
//! * [`input`]: the TOML run description (beam, target, physics, run size).
//! * [`ion`]: ion transport, screening, stopping, damage and crystal targets.
//! * [`material`]: elements combined into materials with their energies.
//! * [`rng`]: deterministic per-particle random streams and the parallel driver.
//! * [`tally`]: histograms, moments, distributions and other result statistics.
//! * [`units`]: unit conventions and conversions.

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
