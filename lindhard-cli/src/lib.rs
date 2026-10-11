//! The pieces of the `lindhard` command that other front ends reuse: the run
//! tally ([`tally`]), the `summary.json` and CSV writers ([`output`]), the
//! one-call simulation driver ([`sim`]), the fluence-stepping driver for
//! `[dynamic]` inputs ([`dynamic`]) and the electron run mode for inputs with
//! an `[electron]` table ([`electron`]), with its cross-section table cache
//! ([`table_cache`]). The `lindhard` binary and the Python
//! bindings (`lindhard-py`) are both thin layers over these, which is what
//! makes their results identical for the same input and seed.

pub mod dynamic;
pub mod electron;
pub mod output;
pub mod publish;
pub mod sim;
pub mod table_cache;
pub mod tally;
