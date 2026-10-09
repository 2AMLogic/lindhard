//! The run description: one TOML document describing the beam, the layered
//! target and its materials, the physics choices and the run size.
//!
//! The same types are serialized back into the output metadata
//! (`summary.json`), so a result carries the input that produced it. Every
//! table uses `deny_unknown_fields`, so a misspelt key is an error naming the
//! key rather than a silently ignored setting. Defaults are filled in on
//! parsing, so the echoed input shows every choice explicitly.
//!
//! The schema is documented with an example in `docs/cli.md`. In short:
//!
//! ```toml
//! [beam]
//! ion = "B"              # element symbol; mass_amu optional
//! energy_ev = 5000.0
//! tilt_deg = 7.0          # polar angle from the surface normal
//!
//! [target]
//! substrate = "Si"        # a name from [materials], or an element symbol
//!
//! [physics]
//! primary_cutoff_ev = 5.0
//! recoil_cutoff_ev = 2.0
//! [physics.energies.Si]
//! e_d_ev = 15.0
//!
//! [run]
//! ions = 1000
//! seed = 1
//! ```
//!
//! [`Input::resolve`] validates the whole description and turns it into the
//! engine's types ([`Resolved`]); [`Resolved::models`] lists every model in
//! use with its citation, for the output metadata.
//!
//! An electron run is described by a different document, with an
//! `[electron]` table and the same `[materials]` and `[target]`; its schema is
//! [`electron`].

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::material::MaterialSpec;

pub mod electron;
mod resolve;
mod resolved;
mod schema;
#[cfg(test)]
mod tests;

pub(crate) use resolved::resolve_material;
pub use resolved::{stopping_model, ModelInfo, Resolved, ResolvedLayer, TABLE_SPEC};
pub use schema::*;

/// Errors from reading or validating an [`Input`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum InputError {
    /// The document is not valid TOML or does not match the schema (unknown
    /// key, wrong type, missing required key). The message from the parser
    /// names the key and the line.
    #[error("{0}")]
    Parse(String),
    /// A value is out of range or inconsistent. `field` is the dotted path of
    /// the offending key, e.g. `target.layers[1].thickness_nm`.
    #[error("{field}: {message}")]
    Invalid {
        /// Dotted path of the key.
        field: String,
        /// What is wrong with it.
        message: String,
    },
}

fn invalid(field: impl Into<String>, message: impl Into<String>) -> InputError {
    InputError::Invalid {
        field: field.into(),
        message: message.into(),
    }
}

/// The whole run description.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    /// The incident beam.
    pub beam: BeamSpec,
    /// Named materials, referenced by name from the target layers.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub materials: BTreeMap<String, MaterialSpec>,
    /// The layered target.
    pub target: TargetSpec,
    /// Physics model choices and cutoffs.
    pub physics: PhysicsSpec,
    /// User-supplied electronic stopping tables (optional).
    #[serde(default, skip_serializing_if = "StoppingSpec::is_empty")]
    pub stopping: StoppingSpec,
    /// Run size, seed and threads.
    pub run: RunSpec,
    /// What the run records.
    #[serde(default)]
    pub tally: TallySpec,
    /// Fluence-dependent target (optional). Absent: the static run.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dynamic: Option<DynamicSpec>,
}

impl Input {
    /// Parse a TOML document. Schema errors (unknown or missing keys, wrong
    /// types) are reported with the key and line; values are checked by
    /// [`Input::resolve`].
    pub fn from_toml_str(text: &str) -> Result<Self, InputError> {
        toml::from_str(text).map_err(|e| InputError::Parse(e.to_string()))
    }

    /// The input as a TOML document in the same schema [`Input::from_toml_str`]
    /// reads: parsing the text gives back an equal [`Input`].
    pub fn to_toml_string(&self) -> Result<String, InputError> {
        toml::to_string(self).map_err(|e| InputError::Parse(e.to_string()))
    }

    /// The input as echoed into output metadata: defaults filled in, and
    /// `run.threads` removed (it does not affect results).
    pub fn echo(&self) -> Input {
        let mut e = self.clone();
        e.run.threads = None;
        if e.physics.screening_length.is_none() {
            e.physics.screening_length = Some(LengthChoice::from_length(
                e.physics.potential.screening().default_length(),
            ));
        }
        e
    }
}
