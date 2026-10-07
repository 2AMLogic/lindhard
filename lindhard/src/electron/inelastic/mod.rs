//! Inelastic electron scattering from dielectric-function models.
//!
//! The first model is the **single-pole Penn algorithm** ([`penn`]): from a
//! user-supplied optical energy-loss function
//! ([`crate::electron::data::OpticalElf`]) it builds the momentum-dependent
//! loss function, the differential inverse inelastic mean free path (DIIMFP),
//! the inelastic mean free path and the stopping power, with nonrelativistic
//! kinematics. [`sum_rules`] checks an optical ELF against the f-sum and
//! perfect-screening sum rules, [`bethe`] gives the nonrelativistic Bethe
//! stopping the model must approach at high energy, and [`drude`] is an
//! analytic Drude-Lorentz ELF with closed-form sum rules, used as a synthetic
//! fixture by the tests, examples and benchmarks.
//!
//! No optical data is committed (`docs/data-provenance.md`); every ELF comes
//! from the caller.
//!
//! The full Penn algorithm, Mermin fits, sampling tables, inner shells with
//! the exchange correction and relativistic kinematics are later work.

pub mod bethe;
pub mod drude;
pub mod penn;
mod quadrature;
pub mod sum_rules;

pub use drude::{DrudeLorentz, DrudeLorentzOscillator};
pub use penn::{InelasticPoint, SinglePolePenn, DEFAULT_RELATIVE_TOLERANCE};
pub use sum_rules::SumRuleReport;
