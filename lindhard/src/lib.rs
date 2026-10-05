//! `lindhard` — clean-room Monte Carlo transport of ions and electrons in matter.
//!
//! The roadmap and module plan are in `docs/architecture.md`. Nothing is
//! implemented yet beyond the shared core: units, constants, the element table,
//! materials and deterministic random streams.

#![forbid(unsafe_code)]

pub mod constants;
pub mod elements;
pub mod ion;
pub mod material;
pub mod rng;
pub mod units;

/// Crate version, as stamped by Cargo.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
