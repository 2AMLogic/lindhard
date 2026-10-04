//! `lindhard` — clean-room Monte Carlo transport of ions and electrons in matter.
//!
//! The roadmap and module plan are in `docs/architecture.md`. Nothing is
//! implemented yet; this crate is a scaffold.

#![forbid(unsafe_code)]

/// Crate version, as stamped by Cargo.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
