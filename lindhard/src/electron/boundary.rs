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
//! No material's band parameters are built in. A [`BandStructure`] is made
//! from values the caller supplies, and it refuses an empty provenance. The
//! one computed helper is [`BandStructure::free_electron_metal`], which takes
//! the Fermi energy from the free-electron relation (Verduin Eq. 3.133) and
//! the layer's atom density. Why there are no tabulated defaults is recorded
//! in `docs/data-provenance.md`.

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
        Ok(Self { model, provenance })
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
    if !(normal_energy_ev > 0.0) || normal_energy_ev + delta_u_ev <= 0.0 {
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
}
