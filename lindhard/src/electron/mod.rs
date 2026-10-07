//! Electron transport (milestone M1, `docs/architecture.md`).
//!
//! [`data`] holds the validated types the electron engine consumes (optical
//! energy-loss functions, subshell binding energies and precomputed
//! cross-section tables). [`elastic`] solves the radial Dirac equation for one
//! screened potential at one energy (phase shifts, Mott cross sections).
//! [`transport`] is the event-by-event loop over layered stacks on those
//! tables. The models that fill the tables from first principles are later
//! work.

pub mod data;
pub mod elastic;
pub mod transport;
