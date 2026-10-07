//! Electron transport (milestone M1, `docs/architecture.md`).
//!
//! [`data`] holds the validated types the electron engine consumes (optical
//! energy-loss functions, subshell binding energies and precomputed
//! cross-section tables). [`elastic`] solves the radial Dirac equation for one
//! screened potential at one energy (phase shifts, Mott cross sections).
//! [`inelastic`] holds the dielectric-function models of inelastic
//! scattering; so far the single-pole Penn algorithm, with its IMFP, stopping
//! power and the sum rules of the optical input. [`transport`] is the
//! event-by-event loop over layered stacks on those tables; [`secondary`] makes
//! secondary electrons at its inelastic events and [`boundary`] holds the
//! per-layer band parameters and the potential step at faces. The models that
//! fill the tables from first principles are later work.

pub mod boundary;
pub mod data;
pub mod elastic;
pub mod inelastic;
pub mod secondary;
pub mod transport;
