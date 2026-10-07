//! Electron transport (milestone M1, `docs/architecture.md`).
//!
//! [`data`] holds the validated types the electron engine consumes (optical
//! energy-loss functions, subshell binding energies and precomputed
//! cross-section tables). [`transport`] is the event-by-event loop over
//! layered stacks on those tables. The scattering models that fill the tables
//! are later work.

pub mod data;
pub mod transport;
