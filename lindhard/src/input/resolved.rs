//! The resolved run: [`Resolved`], its parts and the shared material resolver.

use crate::elements::element_by_symbol;
use crate::geometry::Stack;
use crate::ion::bca::{BcaConfig, Beam, CrystalTarget, MeanFreePath};
use crate::ion::potential::{Screening, ScreeningLength};
use crate::ion::scattering::TableSpec;
use crate::ion::stopping::bethe::BetheBloch;
use crate::ion::stopping::lindhard_scharff::LindhardScharff;
use crate::ion::stopping::table::TableOverride;
use crate::ion::stopping::ElectronicStopping;
use crate::material::{Material, MaterialSpec};
use serde::Serialize;
use std::borrow::Cow;
use std::collections::BTreeMap;

use super::schema::*;
use super::{invalid, Input, InputError};

/// One model in use, for the output metadata.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct ModelInfo {
    /// What the model is for, e.g. `"screening function"`.
    pub role: &'static str,
    /// Stable model name.
    pub name: &'static str,
    /// Published source (for a user table: its path and provenance).
    pub citation: Cow<'static, str>,
}

/// One target layer after resolution.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedLayer {
    /// How the input referred to the material (`[materials]` key, element
    /// symbol, or `"inline"`).
    pub source: String,
    /// The material with every energy set.
    pub material: Material,
}

/// One crystal assignment after resolution.
#[derive(Debug, Clone)]
pub struct ResolvedCrystal {
    /// Stack layers (engine regions) the crystal fills.
    pub regions: Vec<usize>,
    /// Lattice, orientation (with the beam's tilt and azimuth) and search
    /// parameters, ready for [`crate::ion::bca::Bca::with_crystal`].
    pub target: CrystalTarget,
}

/// A validated input, in the engine's types.
#[derive(Debug, Clone)]
pub struct Resolved {
    /// The input with defaults filled in (including the screening length).
    pub input: Input,
    /// The beam.
    pub beam: Beam,
    /// Beam divergence about the nominal direction ([`crate::ion::crystal::Divergence::None`]
    /// without `[beam.divergence]`).
    pub divergence: crate::ion::crystal::Divergence,
    /// The target.
    pub stack: Stack,
    /// Layers in stack order, with where each material came from.
    pub layers: Vec<ResolvedLayer>,
    /// Engine configuration (including the seed).
    pub config: BcaConfig,
    /// Screening function.
    pub screening: Screening,
    /// Screening length.
    pub screening_length: ScreeningLength,
    /// Scattering-table grid.
    pub table_spec: TableSpec,
    /// User stopping tables, in input order (empty without `[stopping]`).
    pub stopping_tables: Vec<LoadedStoppingTable>,
    /// Crystal assignments, in input order (empty: amorphous run).
    pub crystals: Vec<ResolvedCrystal>,
    /// Non-fatal advice (e.g. beam energy outside a model's validity range).
    pub warnings: Vec<String>,
    /// The tuning applied, if `[physics] tuning` named a set.
    pub tuning: Option<TuningReport>,
}

/// Scattering-table grid used for every run: `per_decade` points per decade
/// of reduced energy and impact parameter. Angles outside the energy range
/// fall back to direct quadrature in the engine; impact parameters above
/// `beta_max` give no deflection (the constant-free-path `p_max` is about 20
/// screening lengths or less for solid densities). The measured interpolation
/// error is reported in the output metadata.
pub const TABLE_SPEC: TableSpec = TableSpec {
    eps_min: 1e-6,
    eps_max: 1e4,
    beta_min: 1e-5,
    beta_max: 1e2,
    per_decade: 32,
};

/// Resolve a [`MaterialRef`] against the `[materials]` table: a key of
/// `materials`, else an element symbol (the pure element at its tabulated
/// density), or an inline table. `field` names the key in errors. Shared by
/// the ion input and the electron input ([`electron`]).
pub(crate) fn resolve_material(
    materials: &BTreeMap<String, MaterialSpec>,
    field: &str,
    r: &MaterialRef,
) -> Result<ResolvedLayer, InputError> {
    match r {
        MaterialRef::Inline(spec) => Material::try_from(spec.clone())
            .map(|material| ResolvedLayer {
                source: "inline".into(),
                material,
            })
            .map_err(|e| invalid(field, e.to_string())),
        MaterialRef::Name(name) => {
            if let Some(spec) = materials.get(name) {
                let mut material = Material::try_from(spec.clone())
                    .map_err(|e| invalid(format!("materials.{name}"), e.to_string()))?;
                if material.name().is_none() {
                    material = material.with_name(name.clone());
                }
                return Ok(ResolvedLayer {
                    source: name.clone(),
                    material,
                });
            }
            let el = element_by_symbol(name).ok_or_else(|| {
                invalid(
                    field,
                    format!(
                        "unknown material {name:?}: not a key of [materials] \
                         and not an element symbol"
                    ),
                )
            })?;
            let material = Material::from_atom_fractions(&[(el.z, 1.0)], None)
                .map_err(|e| {
                    invalid(
                        field,
                        format!("element {name}: {e}; define it in [materials] with a density"),
                    )
                })?
                .with_name(name.clone());
            Ok(ResolvedLayer {
                source: name.clone(),
                material,
            })
        }
    }
}

/// The electronic stopping model for a choice. For
/// [`StoppingChoice::EquipartitionLsOr`] the engine does not use this model
/// (it carries its own Lindhard-Scharff/Oen-Robinson mix); Lindhard-Scharff
/// is returned so the validity range can still be checked.
pub fn stopping_model(choice: StoppingChoice) -> Box<dyn ElectronicStopping + Send + Sync> {
    match choice {
        StoppingChoice::LindhardScharff | StoppingChoice::EquipartitionLsOr => {
            Box::new(LindhardScharff::new())
        }
        StoppingChoice::BetheBloch => Box::new(BetheBloch::new()),
    }
}

impl Resolved {
    /// The electronic stopping model of the run: the `[physics] stopping`
    /// choice, with any `[stopping]` tables layered over it for the pairs
    /// they declare. This is what to pass to the engine.
    pub fn stopping_model(&self) -> Box<dyn ElectronicStopping + Send + Sync> {
        let base = stopping_model(self.input.physics.stopping);
        if self.stopping_tables.is_empty() {
            base
        } else {
            Box::new(TableOverride::new(
                self.stopping_tables
                    .iter()
                    .map(|l| l.table.clone())
                    .collect(),
                base,
            ))
        }
    }

    /// Every model in use, with its published source, in a fixed order.
    pub fn models(&self) -> Vec<ModelInfo> {
        let mut v = vec![ModelInfo {
            role: "transport",
            name: "amorphous-bca",
            citation: Cow::Borrowed(
                "J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174 (1980) 257; \
                       M. T. Robinson and I. M. Torrens, Phys. Rev. B 9 (1974) 5008; \
                       W. Eckstein, Computer Simulation of Ion-Solid Interactions (Springer, 1991)",
            ),
        }];
        if !self.crystals.is_empty() {
            v.push(ModelInfo {
                role: "crystal transport",
                name: "lattice-site-bca",
                citation: Cow::Borrowed(
                    "collision partners from explicit lattice sites after DISPLATH \
                     (https://github.com/permissionx/DISPLATH, MIT); cubic lattices: \
                     M. J. Mehl et al., Comput. Mater. Sci. 136 (2017) S1; thermal \
                     displacements from the Debye model (lindhard::ion::crystal::debye); \
                     the electronic-constants-unverified flag and the known channeling \
                     deviations are stated in docs/cli.md, [[crystal]]",
                ),
            });
        }
        v.push(match self.screening {
            Screening::ZblUniversal => ModelInfo {
                role: "screening function",
                name: "zbl-universal",
                citation: Cow::Borrowed(
                    "J. F. Ziegler, J. P. Biersack, U. Littmark, The Stopping and Range \
                           of Ions in Solids (Pergamon, 1985), ch. 2",
                ),
            },
            Screening::KrC => ModelInfo {
                role: "screening function",
                name: "kr-c",
                citation: Cow::Borrowed(
                    "W. D. Wilson, L. G. Haggmark, J. P. Biersack, Phys. Rev. B 15 (1977) 2458",
                ),
            },
            Screening::Moliere => ModelInfo {
                role: "screening function",
                name: "moliere",
                citation: Cow::Borrowed("G. Moliere, Z. Naturforsch. A 2 (1947) 133"),
            },
            Screening::LenzJensen => ModelInfo {
                role: "screening function",
                name: "lenz-jensen",
                citation: Cow::Borrowed(
                    "W. Lenz, Z. Phys. 77 (1932) 713; H. Jensen, Z. Phys. 77 (1932) 722",
                ),
            },
        });
        v.push(match self.screening_length {
            ScreeningLength::Universal => ModelInfo {
                role: "screening length",
                name: "universal",
                citation: Cow::Borrowed(
                    "J. F. Ziegler, J. P. Biersack, U. Littmark, The Stopping and Range \
                           of Ions in Solids (Pergamon, 1985), ch. 2",
                ),
            },
            ScreeningLength::Firsov => ModelInfo {
                role: "screening length",
                name: "firsov",
                citation: Cow::Borrowed("O. B. Firsov, Sov. Phys. JETP 6 (1958) 534"),
            },
            ScreeningLength::Lindhard => ModelInfo {
                role: "screening length",
                name: "lindhard",
                citation: Cow::Borrowed(
                    "J. Lindhard, M. Scharff, H. E. Schiott, Mat. Fys. Medd. Dan. Vid. \
                           Selsk. 33 (14) (1963)",
                ),
            },
        });
        v.push(ModelInfo {
            role: "scattering angle",
            name: "gauss-mehler-quadrature-table",
            citation: Cow::Borrowed(
                "scattering integral by Gauss-Mehler quadrature, tabulated in \
                       (reduced energy, reduced impact parameter); see lindhard::ion::scattering",
            ),
        });
        v.extend(match self.input.physics.stopping {
            StoppingChoice::LindhardScharff => vec![ModelInfo {
                role: "electronic stopping",
                name: "lindhard-scharff",
                citation: Cow::Borrowed("J. Lindhard and M. Scharff, Phys. Rev. 124 (1961) 128"),
            }],
            StoppingChoice::BetheBloch => vec![ModelInfo {
                role: "electronic stopping",
                name: "bethe-bloch",
                citation: Cow::Borrowed(
                    "H. Bethe, Ann. Phys. 5 (1930) 325; F. Bloch, Ann. Phys. 408 (1933) 285; \
                           mean excitation energy by the Bloch rule I = 10 eV Z2",
                ),
            }],
            StoppingChoice::EquipartitionLsOr => vec![
                ModelInfo {
                    role: "electronic stopping",
                    name: "lindhard-scharff",
                    citation: Cow::Borrowed(
                        "J. Lindhard and M. Scharff, Phys. Rev. 124 (1961) 128",
                    ),
                },
                ModelInfo {
                    role: "electronic loss partition",
                    name: "equipartition-ls-or",
                    citation: Cow::Borrowed(
                        "half nonlocal Lindhard-Scharff, half local Oen-Robinson: \
                               O. S. Oen and M. T. Robinson, Nucl. Instrum. Methods 132 (1976) 647",
                    ),
                },
            ],
        });
        for l in &self.stopping_tables {
            v.push(ModelInfo {
                role: "electronic stopping (user table)",
                name: "user-table",
                citation: Cow::Owned(format!("{}: {}", l.path, l.table.provenance())),
            });
        }
        v.push(ModelInfo {
            role: "compound stopping",
            name: "bragg-additivity",
            citation: Cow::Borrowed(
                "W. H. Bragg and R. Kleeman, Phil. Mag. 10 (1905) 318; no compound correction",
            ),
        });
        v.push(match self.config.mean_free_path {
            MeanFreePath::Constant => ModelInfo {
                role: "free path",
                name: "constant",
                citation: Cow::Borrowed(
                    "J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174 (1980) 257",
                ),
            },
            MeanFreePath::EnergyDependent { .. } => ModelInfo {
                role: "free path",
                name: "energy-dependent",
                citation: Cow::Borrowed(
                    "W. Eckstein, Computer Simulation of Ion-Solid Interactions \
                           (Springer, 1991)",
                ),
            },
        });
        v.push(ModelInfo {
            role: "displacement criterion",
            name: "e_d-e_b",
            citation: Cow::Borrowed(
                "J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174 (1980) 257; \
                       W. Eckstein (1991)",
            ),
        });
        v.push(ModelInfo {
            role: "surface barrier",
            name: "planar",
            citation: Cow::Borrowed(
                "W. Eckstein, Computer Simulation of Ion-Solid Interactions (Springer, 1991)",
            ),
        });
        v.push(ModelInfo {
            role: "random numbers",
            name: "chacha8-per-history-stream",
            citation: Cow::Borrowed(
                "D. J. Bernstein, ChaCha, a variant of Salsa20 (2008); \
                       J. K. Salmon et al., Proc. SC'11 (2011)",
            ),
        });
        v
    }
}
