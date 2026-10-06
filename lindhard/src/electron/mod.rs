//! Electron transport (milestone M1, `docs/architecture.md`).
//!
//! Only the data boundary exists so far: [`data`] holds the validated types the
//! electron engine will consume (optical energy-loss functions, subshell
//! binding energies and precomputed cross-section tables). The scattering
//! models and the transport loop are later work.

pub mod data;
