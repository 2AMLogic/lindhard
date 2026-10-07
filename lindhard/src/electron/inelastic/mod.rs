//! Inelastic electron scattering from dielectric-function models.
//!
//! The simplest model is the **single-pole Penn algorithm** ([`penn`]): from a
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
//! [`table`] builds the inelastic energy-loss `CrossSectionTable` and a
//! momentum-transfer sampler from the model.
//!
//! The **full Penn algorithm** ([`full_penn`]) expands the same optical ELF
//! over the Lindhard dielectric functions of free-electron gases
//! ([`lindhard_gas`]) and costs about a second per energy where the single
//! pole costs milliseconds; [`model`] selects either per material and gives
//! the string to record in the metadata of derived tables.
//!
//! Mermin fits and relativistic kinematics are later work.

pub mod bethe;
pub mod drude;
pub mod full_penn;
pub mod inner_shell;
pub mod lindhard_gas;
pub mod model;
pub mod penn;
mod quadrature;
pub mod sum_rules;
pub mod table;

pub use drude::{DrudeLorentz, DrudeLorentzOscillator};
pub use full_penn::{FullPenn, DEFAULT_FULL_TOLERANCE};
pub use inner_shell::{
    Channel, ChannelDiimfp, ChannelInverseImfp, InnerShell, ShellResolvedChannels,
};
pub use lindhard_gas::{LindhardGas, LindhardPlasmon};
pub use model::{PennAlgorithm, PennInelastic};
pub use penn::{
    born_ochkur_factor, ExchangeCorrection, InelasticPoint, SinglePolePenn,
    DEFAULT_RELATIVE_TOLERANCE,
};
pub use sum_rules::SumRuleReport;
