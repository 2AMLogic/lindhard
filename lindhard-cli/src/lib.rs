//! The pieces of the `lindhard` command that other front ends reuse: the run
//! tally ([`tally`]), the `summary.json` and CSV writers ([`output`]), the
//! one-call simulation driver ([`sim`]) and the fluence-stepping driver for
//! `[dynamic]` inputs ([`dynamic`]). The `lindhard` binary and the Python
//! bindings (`lindhard-py`) are both thin layers over these, which is what
//! makes their results identical for the same input and seed.

pub mod dynamic;
pub mod output;
pub mod sim;
pub mod tally;
