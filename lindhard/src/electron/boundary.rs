//! Surfaces and interfaces for the electron loop: per-material band
//! parameters (Fermi energy and inner potential) and the quantum-mechanical
//! potential step.
//!
//! # Energies
//!
//! Inside a material, the electron's kinetic energy is measured from the
//! bottom of that material's band, the convention of Kieft and Bosch,
//! J. Phys. D: Appl. Phys. 41, 215310 (2008), doi:10.1088/0022-3727/41/21/215310,
//! as written out by T. Verduin, *Quantum Noise Effects in e-Beam Lithography
//! and Metrology*, PhD thesis, Delft University of Technology (2017),
//! doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c (cited below as
//! "Verduin" with the thesis's printed page and equation numbers). The
//! **inner potential** `U` is the depth of that band bottom below the vacuum
//! level, so an electron of kinetic energy `E` inside has `E - U` in vacuum.
//! For a metal `U = E_F + Φ` (Fermi energy above the band bottom plus work
//! function; Verduin p. 98 and Fig. 3.28). For a semiconductor or insulator
//! the band bottom is the bottom of the valence band and
//! `U = W_v + E_g + χ` (valence band width, band gap, electron affinity), with
//! the Fermi level placed mid-gap, `E_F = W_v + E_g / 2`. That band model is
//! cstool's (`cstool/input_data/band_structure.py`, see [`BandStructure`]).
//!
//! # The step
//!
//! An electron crossing from a region of inner potential `U` into one of `U'`
//! sees the step `ΔU = U' - U` (Verduin Eq. 3.136, p. 99; vacuum has `U = 0`).
//! With `θ` the angle between its direction and the face normal and
//! `E_n = E cos²θ` its normal energy:
//!
//! - if `E_n + ΔU <= 0` it is always reflected (specularly);
//! - otherwise it is transmitted with probability
//!   `T = 4 s / (1 + s)²`, `s = sqrt(1 + ΔU / E_n)` (Verduin Eq. 3.145, p. 102,
//!   the plane-wave step of Eqs. 3.140-3.144), and reflected otherwise;
//! - on transmission its kinetic energy becomes `E + ΔU` and it is refracted,
//!   conserving the momentum parallel to the face:
//!   `sin θ' = sin θ / sqrt(1 + ΔU / E)` (Verduin Eq. 3.139, p. 100).
//!
//! The step and the refraction are ported from Nebula,
//! `source/physics/boundary_intersect.h` at commit
//! `a50a8e83207980ed221ac83cfe874751d39b86a9`
//! (<https://github.com/Nebula-simulator/nebula>), BSD-3-Clause,
//! Copyright (c) 2020, Nebula-simulator; the notice is in
//! `THIRD_PARTY_LICENSES.md`. Nebula's optional empirical interface absorption
//! (off by default there) is not ported.
//!
//! # Defaults
//!
//! A [`BandStructure`] is made from values the caller supplies, and it
//! refuses an empty provenance. The one computed helper is
//! [`BandStructure::free_electron_metal`], which takes the Fermi energy from
//! the free-electron relation (Verduin Eq. 3.133) and the layer's atom
//! density.
//!
//! A small table of **cited** per-material parameters, [`BAND_DEFAULTS`],
//! covers Al, Cu, Au, W, Si and SiO2. Look an entry up with
//! [`band_defaults`]. Each value in it names the document, table or figure,
//! page, and the date the document was opened. A parameter with no source
//! that could be opened has **no number**: the field is `None` and the
//! entry's [`BandDefaults::gap`] note says what was tried. As of 2026-10-08
//! the table holds:
//!
//! - Si: band gap 1.1 eV and electron affinity 4.05 eV (Robertson and
//!   Wallace 2015). No valence band width.
//! - SiO2: band gap 9 eV (Robertson and Wallace 2015). No electron affinity
//!   and no valence band width.
//! - Al, Cu, Au, W: nothing. No work function and no measured Fermi energy
//!   could be cited to a source that was opened.
//!
//! So no entry is complete yet. [`BandStructure::from_defaults`] builds a
//! band from an entry alone and says which parameters are missing.
//! [`BandDefaults::complete`] takes the missing ones from the caller (with
//! their provenance) and keeps the cited ones, and its provenance string
//! names both. The full list of sources tried, opened or not, is in
//! `docs/data-provenance.md` (rows "Per-material barrier parameters" and
//! "Cited band defaults").

use serde::Serialize;

use crate::constants::{ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR};
use crate::material::Material;

/// Errors from building a [`BandStructure`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum BandError {
    /// A parameter is out of range.
    #[error("invalid {what}: {why}")]
    Invalid {
        /// The parameter.
        what: &'static str,
        /// Why it was rejected.
        why: String,
    },
}

fn invalid<T>(what: &'static str, why: impl Into<String>) -> Result<T, BandError> {
    Err(BandError::Invalid {
        what,
        why: why.into(),
    })
}

/// The band model of one material, after cstool's `band_structure` class
/// (`cstool/input_data/band_structure.py` at commit
/// `0c739eb3fcc3fe5297e74c601ac4a9546db596cf`,
/// <https://github.com/Nebula-simulator/cstool>, BSD-3-Clause, Copyright (c)
/// 2020, Nebula-simulator; notice in `THIRD_PARTY_LICENSES.md`). cstool treats
/// semiconductors and insulators alike, and so does this type.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum BandModel {
    /// A metal: the conduction band is filled up to the Fermi level.
    Metal {
        /// Fermi energy above the bottom of the band, eV.
        fermi_ev: f64,
        /// Work function, eV.
        work_function_ev: f64,
    },
    /// A semiconductor or insulator.
    Insulator {
        /// Width of the valence band, eV (its bottom is the energy origin).
        valence_band_width_ev: f64,
        /// Band gap, eV.
        band_gap_ev: f64,
        /// Electron affinity (vacuum level above the conduction band
        /// bottom), eV. May be negative, as long as the inner potential is
        /// not.
        affinity_ev: f64,
    },
}

/// The band parameters of one layer, with where they came from.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BandStructure {
    model: BandModel,
    /// Binding energy below the Fermi level of the electron that a valence
    /// event liberates in a metal, eV; `None` is the default, the Fermi
    /// level (see [`BandStructure::with_valence_binding_ev`]).
    #[serde(skip_serializing_if = "Option::is_none")]
    valence_binding_ev: Option<f64>,
    provenance: String,
}

impl BandStructure {
    /// Check and wrap. Every parameter must be finite; the Fermi energy, the
    /// work function, the valence band width and the inner potential must be
    /// non-negative and the band gap positive. `provenance` (the source of
    /// the numbers) must not be blank.
    pub fn new(model: BandModel, provenance: impl Into<String>) -> Result<Self, BandError> {
        let provenance = provenance.into();
        if provenance.trim().is_empty() {
            return invalid("provenance", "must name the source of the parameters");
        }
        let nonneg = |what: &'static str, v: f64| -> Result<(), BandError> {
            if v.is_finite() && v >= 0.0 {
                Ok(())
            } else {
                invalid(what, format!("{v} must be finite and non-negative"))
            }
        };
        match model {
            BandModel::Metal {
                fermi_ev,
                work_function_ev,
            } => {
                nonneg("Fermi energy", fermi_ev)?;
                nonneg("work function", work_function_ev)?;
            }
            BandModel::Insulator {
                valence_band_width_ev,
                band_gap_ev,
                affinity_ev,
            } => {
                nonneg("valence band width", valence_band_width_ev)?;
                if !(band_gap_ev.is_finite() && band_gap_ev > 0.0) {
                    return invalid(
                        "band gap",
                        format!("{band_gap_ev} must be finite and positive"),
                    );
                }
                if !affinity_ev.is_finite() {
                    return invalid("electron affinity", "must be finite");
                }
                nonneg(
                    "inner potential",
                    valence_band_width_ev + band_gap_ev + affinity_ev,
                )?;
            }
        }
        Ok(Self {
            model,
            valence_binding_ev: None,
            provenance,
        })
    }

    /// The same metal band, with the electron that a **valence** inelastic
    /// event liberates bound `binding_ev` below the Fermi level instead of
    /// at it (the default, binding 0; [`crate::electron::secondary`],
    /// "Energy"). A valence loss `W > B` then lifts it to `E_F + W - B`
    /// (Verduin Eq. 3.86, p. 78, with this `B`), and a loss `W <= B`
    /// liberates nothing and leaves `W` in the solid.
    ///
    /// This is the rule of M. Azzolini et al., "Secondary electron emission
    /// and yield spectra of metals from Monte Carlo simulations and
    /// experiments", arXiv:1809.00859v1 (2018), p. 6: "should the energy
    /// loss be larger than the first ionization energy B ..., a secondary
    /// electron is emitted with kinetic energy equal to W - B", with `B` of
    /// their Table I, p. 7 ("the mean ionization energy characteristic of
    /// each sample"). Their kinetic energies inside the solid are counted
    /// from the Fermi level (p. 3: the incident energy is increased by the
    /// work function χ, and an electron escapes if `E cos²θ >= χ`), so their
    /// `W - B` is `E_F + W - B` on the band-bottom axis of this module.
    ///
    /// It is a switch for a sensitivity run, off by default; inner-shell
    /// events keep their own shell binding. The source of `binding_ev`
    /// belongs in the band's provenance. Refused: an insulator (its valence
    /// events already use the band gap) and a binding that is not finite and
    /// non-negative.
    pub fn with_valence_binding_ev(mut self, binding_ev: f64) -> Result<Self, BandError> {
        if !matches!(self.model, BandModel::Metal { .. }) {
            return invalid(
                "valence binding energy",
                "only a metal takes one; an insulator's valence events use its band gap",
            );
        }
        if !(binding_ev.is_finite() && binding_ev >= 0.0) {
            return invalid(
                "valence binding energy",
                format!("{binding_ev} must be finite and non-negative"),
            );
        }
        self.valence_binding_ev = Some(binding_ev);
        Ok(self)
    }

    /// The binding energy below the Fermi level of the electron that a
    /// valence event liberates in a metal, eV: the value set by
    /// [`BandStructure::with_valence_binding_ev`], or `None` (the Fermi
    /// level, binding 0).
    pub fn valence_binding_ev(&self) -> Option<f64> {
        self.valence_binding_ev
    }

    /// A free-electron metal: the Fermi energy is computed from the
    /// conduction-electron density `n = z N` (`z` valence electrons per atom,
    /// `N` the material's atom density) by [`free_electron_fermi_energy_ev`];
    /// the work function is the caller's. `z` is a model choice the caller
    /// makes and should name in `provenance` together with the source of the
    /// work function.
    pub fn free_electron_metal(
        material: &Material,
        valence_electrons_per_atom: f64,
        work_function_ev: f64,
        provenance: impl Into<String>,
    ) -> Result<Self, BandError> {
        let z = valence_electrons_per_atom;
        if !(z.is_finite() && z > 0.0) {
            return invalid(
                "valence electrons per atom",
                format!("{z} must be finite and positive"),
            );
        }
        let n = z * material.atom_number_density();
        Self::new(
            BandModel::Metal {
                fermi_ev: free_electron_fermi_energy_ev(n),
                work_function_ev,
            },
            provenance,
        )
    }

    /// The band model.
    pub fn model(&self) -> BandModel {
        self.model
    }

    /// Where the parameters came from.
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// Fermi energy above the band bottom, eV: the metal's own, or mid-gap
    /// (`W_v + E_g / 2`) for an insulator (cstool `get_fermi`).
    pub fn fermi_ev(&self) -> f64 {
        match self.model {
            BandModel::Metal { fermi_ev, .. } => fermi_ev,
            BandModel::Insulator {
                valence_band_width_ev,
                band_gap_ev,
                ..
            } => valence_band_width_ev + 0.5 * band_gap_ev,
        }
    }

    /// The lowest energy above the band bottom at which an electron can
    /// excite a secondary, eV: the Fermi energy of a metal, the
    /// conduction-band bottom `W_v + E_g` of an insulator, where both the
    /// primary and the secondary must end in the conduction band (cstool
    /// `get_min_excitation`, `cstool/input_data/band_structure.py` at commit
    /// `0c739eb3fcc3fe5297e74c601ac4a9546db596cf`). cstool's table compiler caps the loss of an electron of
    /// energy `K` at `K` minus this energy
    /// ([`crate::electron::inelastic::table`], "Energy axis").
    pub fn min_excitation_ev(&self) -> f64 {
        match self.model {
            BandModel::Metal { fermi_ev, .. } => fermi_ev,
            BandModel::Insulator {
                valence_band_width_ev,
                band_gap_ev,
                ..
            } => valence_band_width_ev + band_gap_ev,
        }
    }

    /// Inner potential `U` (vacuum level above the band bottom), eV:
    /// `E_F + Φ` for a metal, `W_v + E_g + χ` for an insulator (cstool
    /// `get_barrier`; Verduin p. 98 and Eq. 3.136).
    pub fn inner_potential_ev(&self) -> f64 {
        match self.model {
            BandModel::Metal {
                fermi_ev,
                work_function_ev,
            } => fermi_ev + work_function_ev,
            BandModel::Insulator {
                valence_band_width_ev,
                band_gap_ev,
                affinity_ev,
            } => valence_band_width_ev + band_gap_ev + affinity_ev,
        }
    }

    /// The band gap, eV, or `None` for a metal.
    pub fn band_gap_ev(&self) -> Option<f64> {
        match self.model {
            BandModel::Metal { .. } => None,
            BandModel::Insulator { band_gap_ev, .. } => Some(band_gap_ev),
        }
    }
}

/// Whether a [`BandDefaults`] entry describes a metal or an insulator, i.e.
/// which [`BandModel`] it builds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BandKind {
    /// [`BandModel::Metal`]: Fermi energy and work function.
    Metal,
    /// [`BandModel::Insulator`]: valence band width, band gap and electron
    /// affinity.
    Insulator,
}

/// One cited band parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CitedValue {
    /// The value, eV.
    pub value_ev: f64,
    /// Where it was read: document, table or figure, page, and the date the
    /// document was opened.
    pub source: &'static str,
}

/// The cited band parameters of one material. A field is `None` where no
/// source that could be opened gives the value; [`BandDefaults::gap`] says
/// what was tried. Fields that do not belong to the entry's [`BandKind`] are
/// always `None`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandDefaults {
    /// The material, as a chemical formula (`"Si"`, `"SiO2"`, ...).
    pub material: &'static str,
    /// Metal or insulator.
    pub kind: BandKind,
    /// Fermi energy above the band bottom, eV (metals).
    pub fermi_ev: Option<CitedValue>,
    /// Work function, eV (metals).
    pub work_function_ev: Option<CitedValue>,
    /// Valence band width, eV (insulators).
    pub valence_band_width_ev: Option<CitedValue>,
    /// Band gap, eV (insulators).
    pub band_gap_ev: Option<CitedValue>,
    /// Electron affinity, eV (insulators).
    pub affinity_ev: Option<CitedValue>,
    /// Why the missing parameters are missing.
    pub gap: &'static str,
}

/// Caller-supplied values for the parameters a [`BandDefaults`] entry lacks.
/// Leave a field `None` to use the cited value; a field set for a parameter
/// the entry already cites is refused (build with [`BandStructure::new`] to
/// replace a cited value).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BandFill {
    /// Fermi energy above the band bottom, eV (metals).
    pub fermi_ev: Option<f64>,
    /// Work function, eV (metals).
    pub work_function_ev: Option<f64>,
    /// Valence band width, eV (insulators).
    pub valence_band_width_ev: Option<f64>,
    /// Band gap, eV (insulators).
    pub band_gap_ev: Option<f64>,
    /// Electron affinity, eV (insulators).
    pub affinity_ev: Option<f64>,
}

/// J. Robertson and R. M. Wallace, read in the authors' manuscript.
macro_rules! robertson_wallace {
    () => {
        "J. Robertson and R. M. Wallace, \"High-K materials and metal gates for CMOS \
         applications\", Mater. Sci. Eng. R 88, 1-41 (2015), doi:10.1016/j.mser.2014.11.001, \
         read in the authors' manuscript in the University of Cambridge repository \
         (hdl:1810/246441; page numbers are the manuscript's), opened 2026-10-08"
    };
}

/// What was tried for the metals' parameters (the full list is in
/// `docs/data-provenance.md`).
macro_rules! metal_gap {
    () => {
        "no work function and no measured Fermi energy could be cited to a source that was \
         opened (searched 2026-10-07 and 2026-10-08): the compilations of Michaelson (1977) and \
         Derry et al. (2015) and the primary photoemission papers are closed access; the open \
         NASA report of J. R. Smith (1969) only repeats Fomenko's 1966 recommended values; see \
         docs/data-provenance.md, 'Per-material barrier parameters'. \
         BandStructure::free_electron_metal computes a Fermi energy from a valence count"
    };
}

const fn metal(material: &'static str) -> BandDefaults {
    BandDefaults {
        material,
        kind: BandKind::Metal,
        fermi_ev: None,
        work_function_ev: None,
        valence_band_width_ev: None,
        band_gap_ev: None,
        affinity_ev: None,
        gap: metal_gap!(),
    }
}

/// The cited band parameters, one entry per material (see the module docs,
/// "Defaults"). Every value cites its table or figure and page.
pub const BAND_DEFAULTS: &[BandDefaults] = &[
    metal("Al"),
    metal("Cu"),
    metal("Au"),
    metal("W"),
    BandDefaults {
        material: "Si",
        kind: BandKind::Insulator,
        fermi_ev: None,
        work_function_ev: None,
        valence_band_width_ev: None,
        // Robertson and Wallace, Table 1, manuscript p. 43: "experimental
        // band gap", Si 1.1 eV.
        band_gap_ev: Some(CitedValue {
            value_ev: 1.1,
            source: concat!(
                robertson_wallace!(),
                ", Table 1, p. 43 (experimental band gap of Si, 1.1 eV; Table 3, p. 45, also \
                 prints 1.1 eV); cross-check: S. Tanuma, C. J. Powell and D. R. Penn, Surf. \
                 Interface Anal. 43, 689 (2011), doi:10.1002/sia.3522, Table 1 (authors' \
                 manuscript at NIMS MDR, p. 32), Eg = 1.1 eV, opened 2026-10-08"
            ),
        }),
        // Robertson and Wallace, Fig. 42, manuscript p. 63: chi = 4.05 eV
        // drawn from the vacuum level to the Si conduction band; Section 6.1,
        // p. 20: n+ and p+ poly-Si gates have work functions of 4.05 and
        // 5.15 eV (the two band edges), a difference of the 1.1 eV gap.
        affinity_ev: Some(CitedValue {
            value_ev: 4.05,
            source: concat!(
                robertson_wallace!(),
                ", Fig. 42, p. 63 (electron affinity of Si, chi = 4.05 eV, vacuum level to \
                 conduction band edge; the valence band edge is at 5.15 eV), and Section 6.1, \
                 p. 20 (n+ and p+ poly-Si gate work functions 4.05 and 5.15 eV)"
            ),
        }),
        gap: "valence band width: no source that was opened gives the bottom of the Si valence \
              band relative to its top (L. Ley et al., LBL-688 (1972), the preprint of Phys. Rev. \
              Lett. 29, 1088, Table 1, place the Gamma1 feature at 14.7 eV below the Fermi level \
              only); see docs/data-provenance.md",
    },
    BandDefaults {
        material: "SiO2",
        kind: BandKind::Insulator,
        fermi_ev: None,
        work_function_ev: None,
        valence_band_width_ev: None,
        // Robertson and Wallace, Table 1, manuscript p. 43: SiO2 gap 9 eV.
        band_gap_ev: Some(CitedValue {
            value_ev: 9.0,
            source: concat!(
                robertson_wallace!(),
                ", Table 1, p. 43 (experimental band gap of SiO2, 9 eV; Table 4, p. 45, prints \
                 9.0 eV); cross-check: S. Tanuma, C. J. Powell and D. R. Penn, Surf. Interface \
                 Anal. 17, 927 (1991), doi:10.1002/sia.740171305, Table 5 (authors' manuscript \
                 at NIMS MDR, p. 21), Eg = 9.1 eV, opened 2026-10-08"
            ),
        }),
        affinity_ev: None,
        gap: "electron affinity and valence band width: no source that was opened states them. \
              Robertson and Wallace (2015), Table 1, give only the Si/SiO2 conduction band offset \
              (3.2 eV), which equals the affinity difference only without an interface dipole \
              (their Section 4.4), so no affinity is derived from it; see docs/data-provenance.md",
    },
];

/// The cited band parameters of `material` (a chemical formula such as
/// `"Si"` or `"SiO2"`, matched exactly), or `None` if the table has no entry
/// for it.
pub fn band_defaults(material: &str) -> Option<&'static BandDefaults> {
    BAND_DEFAULTS.iter().find(|d| d.material == material)
}

impl BandDefaults {
    /// Build a [`BandStructure`] from the cited values, taking each missing
    /// parameter from `fill`. `provenance` names the source of the filled
    /// values and must not be blank when any is filled. Refused: a parameter
    /// with neither a cited value nor a fill (the error quotes
    /// [`BandDefaults::gap`]), a fill for a parameter that is already cited,
    /// and a fill for a parameter of the other [`BandKind`]. The returned
    /// band's provenance lists every cited value with its source, then the
    /// filled values with `provenance`.
    pub fn complete(&self, fill: BandFill, provenance: &str) -> Result<BandStructure, BandError> {
        let (needed, foreign): (&[_], &[_]) = match self.kind {
            BandKind::Metal => (
                &[
                    ("Fermi energy", self.fermi_ev, fill.fermi_ev),
                    (
                        "work function",
                        self.work_function_ev,
                        fill.work_function_ev,
                    ),
                ],
                &[
                    ("valence band width", fill.valence_band_width_ev),
                    ("band gap", fill.band_gap_ev),
                    ("electron affinity", fill.affinity_ev),
                ],
            ),
            BandKind::Insulator => (
                &[
                    (
                        "valence band width",
                        self.valence_band_width_ev,
                        fill.valence_band_width_ev,
                    ),
                    ("band gap", self.band_gap_ev, fill.band_gap_ev),
                    ("electron affinity", self.affinity_ev, fill.affinity_ev),
                ],
                &[
                    ("Fermi energy", fill.fermi_ev),
                    ("work function", fill.work_function_ev),
                ],
            ),
        };
        let material = self.material;
        if let Some((what, _)) = foreign.iter().find(|(_, v)| v.is_some()) {
            let kind = match self.kind {
                BandKind::Metal => "a metal",
                BandKind::Insulator => "a semiconductor or insulator",
            };
            return invalid(what, format!("{material} is {kind}; it has no {what}"));
        }
        let mut values = Vec::with_capacity(needed.len());
        let mut cited = Vec::new();
        let mut filled = Vec::new();
        for &(what, default, given) in needed {
            match (default, given) {
                (Some(c), None) => {
                    values.push(c.value_ev);
                    cited.push(format!("{what} {} eV from {}", c.value_ev, c.source));
                }
                (Some(c), Some(_)) => {
                    return invalid(
                        what,
                        format!(
                            "{material} has a cited default ({} eV); build with \
                             BandStructure::new to replace it",
                            c.value_ev
                        ),
                    )
                }
                (None, Some(v)) => {
                    values.push(v);
                    filled.push(format!("{what} {v} eV"));
                }
                (None, None) => {
                    return invalid(
                        what,
                        format!(
                            "no cited default for {material} ({}); supply it with \
                             BandDefaults::complete",
                            self.gap
                        ),
                    )
                }
            }
        }
        if !filled.is_empty() && provenance.trim().is_empty() {
            return invalid(
                "provenance",
                "must name the source of the caller-supplied parameters",
            );
        }
        let mut text = format!(
            "{material}: cited defaults (lindhard::electron::boundary::BAND_DEFAULTS): {}",
            if cited.is_empty() {
                "none".to_string()
            } else {
                cited.join("; ")
            }
        );
        if !filled.is_empty() {
            text.push_str(&format!(
                ". Caller-supplied: {}, from: {}",
                filled.join(", "),
                provenance.trim()
            ));
        }
        let model = match self.kind {
            BandKind::Metal => BandModel::Metal {
                fermi_ev: values[0],
                work_function_ev: values[1],
            },
            BandKind::Insulator => BandModel::Insulator {
                valence_band_width_ev: values[0],
                band_gap_ev: values[1],
                affinity_ev: values[2],
            },
        };
        BandStructure::new(model, text)
    }
}

impl BandStructure {
    /// The band of `material` from [`BAND_DEFAULTS`] alone. Fails if the
    /// table has no entry for `material` or the entry lacks a parameter (the
    /// error names the first missing one and why it is missing); as of
    /// 2026-10-08 every entry lacks at least one, so use
    /// [`BandDefaults::complete`] to supply the rest.
    pub fn from_defaults(material: &str) -> Result<Self, BandError> {
        let Some(d) = band_defaults(material) else {
            let known: Vec<&str> = BAND_DEFAULTS.iter().map(|d| d.material).collect();
            return invalid(
                "material",
                format!(
                    "no band defaults for {material:?}; the table covers {}",
                    known.join(", ")
                ),
            );
        };
        d.complete(BandFill::default(), "")
    }
}

/// Free-electron Fermi energy, eV, for `n` conduction electrons per m³:
/// Verduin Eq. 3.133 (p. 95), `n = (2 sqrt 2 / 3π²) (m_e E_F / ħ²)^(3/2)`,
/// solved for `E_F`, which gives `E_F = (ħ² / 2 m_e) (3π² n)^(2/3)`. CODATA
/// constants from [`crate::constants`].
pub fn free_electron_fermi_energy_ev(n_per_m3: f64) -> f64 {
    let k_f = (3.0 * std::f64::consts::PI.powi(2) * n_per_m3).cbrt();
    HBAR * HBAR * k_f * k_f / (2.0 * ELECTRON_MASS) / ELEMENTARY_CHARGE
}

/// Conduction-electron density, per m³, at Fermi energy `fermi_ev`: Verduin
/// Eq. 3.133 (p. 95) as printed. The inverse of
/// [`free_electron_fermi_energy_ev`].
pub fn free_electron_density_per_m3(fermi_ev: f64) -> f64 {
    let e = fermi_ev * ELEMENTARY_CHARGE;
    2.0 * std::f64::consts::SQRT_2 / (3.0 * std::f64::consts::PI.powi(2))
        * (ELECTRON_MASS * e / (HBAR * HBAR)).powf(1.5)
}

/// Transmission probability through a potential step `delta_u_ev` for an
/// electron of normal energy `normal_energy_ev = E cos²θ` (Verduin Eq. 3.145,
/// p. 102): zero when `E_n + ΔU <= 0` (the electron cannot get over), else
/// `4 s / (1 + s)²` with `s = sqrt(1 + ΔU / E_n)`. Symmetric under exchanging
/// the two sides (`s -> 1/s`).
pub fn step_transmission(normal_energy_ev: f64, delta_u_ev: f64) -> f64 {
    if normal_energy_ev.is_nan() || normal_energy_ev <= 0.0 || normal_energy_ev + delta_u_ev <= 0.0
    {
        return 0.0;
    }
    let s = (1.0 + delta_u_ev / normal_energy_ev).sqrt();
    4.0 * s / ((1.0 + s) * (1.0 + s))
}

/// The outcome of [`cross_step`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StepOutcome {
    /// The electron got through: its new direction and kinetic energy.
    Transmitted {
        /// Unit direction after the face.
        dir: [f64; 3],
        /// Kinetic energy after the face, eV (`E + ΔU`).
        energy_ev: f64,
    },
    /// The electron was reflected specularly: its new direction.
    Reflected {
        /// Unit direction after the reflection (the normal component
        /// reversed).
        dir: [f64; 3],
    },
}

/// One electron at a face normal to `x`, ported from Nebula's
/// `boundary_intersect::execute` (file and commit in the module docs).
///
/// `dir` is the unit direction (its `x` component must be nonzero; its sign
/// says which way the electron crosses), `energy_ev` the kinetic energy on
/// the near side, `delta_u_ev` the step `U' - U`. `u` is a uniform draw on
/// `[0, 1)`; it is only consulted when the electron has the normal energy to
/// get over, and then transmission happens when `u < T`. With
/// `quantum_transmission` false, `T = 1` whenever the electron can get over;
/// with `refraction` false the direction is kept on transmission.
pub fn cross_step(
    dir: [f64; 3],
    energy_ev: f64,
    delta_u_ev: f64,
    quantum_transmission: bool,
    refraction: bool,
    u: f64,
) -> StepOutcome {
    let mu = dir[0];
    let normal = energy_ev * mu * mu;
    if normal + delta_u_ev > 0.0 {
        let s = (1.0 + delta_u_ev / normal).sqrt();
        let t = if quantum_transmission {
            4.0 * s / ((1.0 + s) * (1.0 + s))
        } else {
            1.0
        };
        if u < t {
            // Verduin Eq. 3.139 in Nebula's vector form: the normal component
            // is scaled by `s`, the parallel ones are kept.
            let dir = if refraction {
                let v = [s * mu, dir[1], dir[2]];
                let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                [v[0] / n, v[1] / n, v[2] / n]
            } else {
                dir
            };
            return StepOutcome::Transmitted {
                dir,
                energy_ev: energy_ev + delta_u_ev,
            };
        }
    }
    StepOutcome::Reflected {
        dir: [-mu, dir[1], dir[2]],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transmission_is_one_without_a_step_and_symmetric() {
        assert_eq!(step_transmission(3.0, 0.0), 1.0);
        // Crossing up by dU at normal energy E_n is the reverse of crossing
        // down by -dU at E_n + dU: the same s, inverted.
        for (en, du) in [(2.0, 5.0), (10.0, 0.5), (0.3, 12.0)] {
            let up = step_transmission(en, du);
            let down = step_transmission(en + du, -du);
            assert!((up - down).abs() < 1e-14, "{up} vs {down}");
        }
    }

    #[test]
    fn below_the_barrier_is_reflected_for_every_draw() {
        for u in [0.0, 1e-12, 0.5, 0.999_999] {
            let d = [-0.6, 0.8, 0.0];
            // E cos^2 = 10 * 0.36 = 3.6 < 4.
            match cross_step(d, 10.0, -4.0, true, true, u) {
                StepOutcome::Reflected { dir } => assert_eq!(dir, [0.6, 0.8, 0.0]),
                o => panic!("transmitted below the barrier: {o:?}"),
            }
        }
    }

    #[test]
    fn refraction_conserves_parallel_momentum() {
        let d = {
            let v: [f64; 3] = [0.5, 0.7, -0.3];
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            [v[0] / n, v[1] / n, v[2] / n]
        };
        for du in [-1.5, 2.0, 7.0] {
            let e = 9.0;
            match cross_step(d, e, du, false, true, 0.0) {
                StepOutcome::Transmitted { dir, energy_ev } => {
                    assert_eq!(energy_ev, e + du);
                    let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
                    assert!((len - 1.0).abs() < 1e-14);
                    assert!(dir[0] * d[0] > 0.0, "kept its side of travel");
                    // Verduin Eq. 3.138: sqrt(E) sin(theta) = sqrt(E') sin(theta').
                    let sin = (1.0 - d[0] * d[0]).sqrt();
                    let sin2 = (1.0 - dir[0] * dir[0]).sqrt();
                    assert!((e.sqrt() * sin - energy_ev.sqrt() * sin2).abs() < 1e-13);
                    // The azimuth is unchanged.
                    assert!((dir[1] * d[2] - dir[2] * d[1]).abs() < 1e-14);
                }
                o => panic!("expected transmission: {o:?}"),
            }
        }
    }

    #[test]
    fn fermi_energy_inverts_eq_3_133() {
        // Arbitrary densities: a round trip, not a material check.
        for n in [1e27, 4e28, 2e29] {
            let ef = free_electron_fermi_energy_ev(n);
            let back = free_electron_density_per_m3(ef);
            assert!((back - n).abs() <= 1e-12 * n, "{n} -> {ef} eV -> {back}");
        }
        // E_F scales as n^(2/3).
        let r = free_electron_fermi_energy_ev(8e28) / free_electron_fermi_energy_ev(1e28);
        assert!((r - 4.0).abs() < 1e-12);
    }

    #[test]
    fn band_structure_checks_its_inputs() {
        let metal = |f, w| {
            BandStructure::new(
                BandModel::Metal {
                    fermi_ev: f,
                    work_function_ev: w,
                },
                "test",
            )
        };
        assert!(metal(5.0, 4.0).is_ok());
        assert!(metal(-1.0, 4.0).is_err());
        assert!(metal(5.0, f64::NAN).is_err());
        assert!(BandStructure::new(
            BandModel::Metal {
                fermi_ev: 5.0,
                work_function_ev: 4.0
            },
            "  "
        )
        .is_err());
        let ins = |v, g, a| {
            BandStructure::new(
                BandModel::Insulator {
                    valence_band_width_ev: v,
                    band_gap_ev: g,
                    affinity_ev: a,
                },
                "test",
            )
        };
        let b = ins(6.0, 2.0, 1.0).unwrap();
        assert_eq!(b.fermi_ev(), 7.0);
        assert_eq!(b.inner_potential_ev(), 9.0);
        assert_eq!(b.band_gap_ev(), Some(2.0));
        assert!(ins(6.0, 0.0, 1.0).is_err());
        assert!(ins(0.0, 1.0, -2.0).is_err(), "negative inner potential");
        assert!(ins(6.0, 1.0, -0.5).is_ok(), "negative affinity is allowed");
    }

    /// Valence band width used only to complete the Si and SiO2 entries in
    /// these tests: a synthetic round number, not a material value.
    const SYNTHETIC_WV: f64 = 12.0;
    const SYNTHETIC: &str = "synthetic test value, not material data";

    #[test]
    fn silicon_defaults_match_robertson_wallace() {
        let si = band_defaults("Si").expect("Si entry");
        assert_eq!(si.kind, BandKind::Insulator);
        let gap = si.band_gap_ev.expect("Si band gap");
        assert_eq!(gap.value_ev, 1.1);
        for s in [
            "Robertson and R. M. Wallace",
            "Table 1, p. 43",
            "2026-10-08",
        ] {
            assert!(gap.source.contains(s), "{s:?} missing from {}", gap.source);
        }
        let chi = si.affinity_ev.expect("Si electron affinity");
        assert_eq!(chi.value_ev, 4.05);
        for s in ["Fig. 42, p. 63", "Section 6.1", "2026-10-08"] {
            assert!(chi.source.contains(s), "{s:?} missing from {}", chi.source);
        }
        assert!(si.valence_band_width_ev.is_none(), "no cited W_v");
        assert!(si.gap.contains("valence band width"));
    }

    #[test]
    fn silicon_dioxide_defaults_match_robertson_wallace() {
        let d = band_defaults("SiO2").expect("SiO2 entry");
        assert_eq!(d.kind, BandKind::Insulator);
        let gap = d.band_gap_ev.expect("SiO2 band gap");
        assert_eq!(gap.value_ev, 9.0);
        for s in ["Table 1, p. 43", "Table 5", "2026-10-08"] {
            assert!(gap.source.contains(s), "{s:?} missing from {}", gap.source);
        }
        assert!(d.affinity_ev.is_none() && d.valence_band_width_ev.is_none());
        assert!(d.gap.contains("electron affinity"));
    }

    #[test]
    fn metals_have_no_numbers() {
        for m in ["Al", "Cu", "Au", "W"] {
            let d = band_defaults(m).expect("metal entry");
            assert_eq!(d.kind, BandKind::Metal);
            assert!(d.fermi_ev.is_none() && d.work_function_ev.is_none(), "{m}");
            assert!(d.gap.contains("work function"), "{m}");
        }
    }

    #[test]
    fn entries_are_unique_and_hold_only_their_kind() {
        for (i, d) in BAND_DEFAULTS.iter().enumerate() {
            assert!(
                BAND_DEFAULTS[..i].iter().all(|e| e.material != d.material),
                "duplicate {}",
                d.material
            );
            assert!(!d.gap.trim().is_empty());
            let (own, other) = match d.kind {
                BandKind::Metal => (
                    vec![d.fermi_ev, d.work_function_ev],
                    vec![d.valence_band_width_ev, d.band_gap_ev, d.affinity_ev],
                ),
                BandKind::Insulator => (
                    vec![d.valence_band_width_ev, d.band_gap_ev, d.affinity_ev],
                    vec![d.fermi_ev, d.work_function_ev],
                ),
            };
            assert!(other.iter().all(Option::is_none), "{}", d.material);
            for c in own.into_iter().flatten() {
                assert!(c.value_ev.is_finite(), "{}", d.material);
                assert!(c.source.contains("opened 2026-"), "{}", c.source);
            }
        }
        assert!(band_defaults("si").is_none(), "matched exactly");
        assert!(band_defaults("Ge").is_none());
    }

    #[test]
    fn no_entry_is_complete_on_its_own() {
        for d in BAND_DEFAULTS {
            let e = BandStructure::from_defaults(d.material).unwrap_err();
            assert!(e.to_string().contains("no cited default"), "{e}");
        }
        let e = BandStructure::from_defaults("Ge").unwrap_err();
        assert!(e.to_string().contains("Si, SiO2"), "{e}");
    }

    #[test]
    fn completed_silicon_passes_the_checks_and_cites_both_sources() {
        let fill = BandFill {
            valence_band_width_ev: Some(SYNTHETIC_WV),
            ..BandFill::default()
        };
        let b = band_defaults("Si")
            .unwrap()
            .complete(fill, SYNTHETIC)
            .unwrap();
        assert_eq!(
            b.model(),
            BandModel::Insulator {
                valence_band_width_ev: SYNTHETIC_WV,
                band_gap_ev: 1.1,
                affinity_ev: 4.05,
            }
        );
        assert_eq!(b.band_gap_ev(), Some(1.1));
        assert!((b.inner_potential_ev() - (SYNTHETIC_WV + 1.1 + 4.05)).abs() < 1e-12);
        // The same numbers pass `BandStructure::new` directly.
        assert!(BandStructure::new(b.model(), SYNTHETIC).is_ok());
        let p = b.provenance();
        for s in [
            "Fig. 42, p. 63",
            "Table 1, p. 43",
            "valence band width 12 eV",
            SYNTHETIC,
        ] {
            assert!(p.contains(s), "{s:?} missing from {p}");
        }
    }

    #[test]
    fn completion_refuses_bad_fills() {
        let si = band_defaults("Si").unwrap();
        let wv = BandFill {
            valence_band_width_ev: Some(SYNTHETIC_WV),
            ..BandFill::default()
        };
        // Replacing a cited value.
        let e = si
            .complete(
                BandFill {
                    band_gap_ev: Some(1.2),
                    ..wv
                },
                SYNTHETIC,
            )
            .unwrap_err();
        assert!(e.to_string().contains("cited default"), "{e}");
        // A metal parameter on an insulator.
        let e = si
            .complete(
                BandFill {
                    work_function_ev: Some(4.0),
                    ..wv
                },
                SYNTHETIC,
            )
            .unwrap_err();
        assert!(e.to_string().contains("work function"), "{e}");
        // Filled values need a provenance.
        assert!(si.complete(wv, "  ").is_err());
        // Out-of-range fills still go through `BandStructure::new`.
        let neg = BandFill {
            valence_band_width_ev: Some(-1.0),
            ..BandFill::default()
        };
        assert!(si.complete(neg, SYNTHETIC).is_err());
        // SiO2 still lacks its affinity with W_v filled.
        let e = band_defaults("SiO2")
            .unwrap()
            .complete(wv, SYNTHETIC)
            .unwrap_err();
        assert!(e.to_string().contains("electron affinity"), "{e}");
        // A metal completed entirely from the caller.
        let cu = band_defaults("Cu").unwrap();
        let b = cu
            .complete(
                BandFill {
                    fermi_ev: Some(5.0),
                    work_function_ev: Some(4.0),
                    ..BandFill::default()
                },
                SYNTHETIC,
            )
            .unwrap();
        assert_eq!(b.inner_potential_ev(), 9.0);
        assert!(b.provenance().contains("cited defaults"));
        assert!(cu
            .complete(
                BandFill {
                    fermi_ev: Some(5.0),
                    work_function_ev: Some(4.0),
                    band_gap_ev: Some(1.0),
                    ..BandFill::default()
                },
                SYNTHETIC,
            )
            .is_err());
    }
}
