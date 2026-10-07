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
//! [`inner_shell`] resolves the losses into a valence channel and
//! inner-shell ionization channels (secondary energy `ω - B`), each from its
//! own caller-supplied optical ELF, and the optional Born-Ochkur
//! [`ExchangeCorrection`] makes the primary and the struck electron
//! indistinguishable at low energy.
//!
//! The full Penn algorithm, Mermin fits, sampling tables and relativistic
//! kinematics are later work.

pub mod bethe;
pub mod drude;
pub mod inner_shell;
pub mod penn;
mod quadrature;
pub mod sum_rules;

pub use drude::{DrudeLorentz, DrudeLorentzOscillator};
pub use inner_shell::{
    Channel, ChannelDiimfp, ChannelInverseImfp, InnerShell, ShellResolvedChannels,
};
pub use penn::{
    born_ochkur_factor, ExchangeCorrection, InelasticPoint, SinglePolePenn,
    DEFAULT_RELATIVE_TOLERANCE,
};
pub use sum_rules::SumRuleReport;
