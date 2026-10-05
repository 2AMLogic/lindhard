//! Bragg additivity for compounds and mixtures.
//!
//! W. H. Bragg and R. Kleeman, Phil. Mag. 10, 318 (1905): the stopping cross
//! section of a compound is the atom-fraction-weighted sum of the element
//! cross sections,
//!
//! `S_compound (per average atom) = Σ_j x_j S_j`,
//!
//! so `-dE/dx = N Σ_j x_j S_j · f`, with `N` the total atom density and `f` a
//! per-compound correction factor ([`CompoundCorrection`]; default 1) as a hook
//! for chemical/phase corrections (e.g. cores-and-bonds) supplied by the caller
//! from a cited source.

use super::{ElectronicStopping, Ion, StoppingError};
use crate::material::Material;

/// Per-compound correction factor applied to the Bragg sum.
pub trait CompoundCorrection {
    /// Multiplicative factor at the given ion energy. Must be finite and
    /// non-negative (zero is allowed); [`bragg_cross_section_per_atom`] rejects
    /// anything else with [`StoppingError::InvalidParameter`].
    fn factor(&self, material: &Material, ion: &Ion, energy_ev: f64) -> f64;
}

/// No correction: factor exactly 1.
#[derive(Debug, Clone, Copy, Default)]
pub struct NoCorrection;

impl CompoundCorrection for NoCorrection {
    fn factor(&self, _: &Material, _: &Ion, _: f64) -> f64 {
        1.0
    }
}

/// An energy-independent constant factor.
#[derive(Debug, Clone, Copy)]
pub struct ConstantCorrection(pub f64);

impl CompoundCorrection for ConstantCorrection {
    fn factor(&self, _: &Material, _: &Ion, _: f64) -> f64 {
        self.0
    }
}

/// Stopping cross section per average atom of `material`, J m²
/// (`f · Σ_j x_j S_j`).
pub fn bragg_cross_section_per_atom(
    model: &dyn ElectronicStopping,
    correction: &dyn CompoundCorrection,
    ion: &Ion,
    material: &Material,
    energy_ev: f64,
) -> Result<f64, StoppingError> {
    let mut sum = 0.0;
    for c in material.components() {
        sum += c.atom_fraction() * model.stopping(ion, c.z(), energy_ev)?;
    }
    let factor = correction.factor(material, ion, energy_ev);
    if !(factor.is_finite() && factor >= 0.0) {
        return Err(StoppingError::InvalidParameter {
            name: "compound_correction",
            value: factor,
        });
    }
    Ok(factor * sum)
}

/// Linear stopping power `-dE/dx`, J/m: `N` times
/// [`bragg_cross_section_per_atom`].
pub fn bragg_stopping_power(
    model: &dyn ElectronicStopping,
    correction: &dyn CompoundCorrection,
    ion: &Ion,
    material: &Material,
    energy_ev: f64,
) -> Result<f64, StoppingError> {
    Ok(material.atom_number_density()
        * bragg_cross_section_per_atom(model, correction, ion, material, energy_ev)?)
}

#[cfg(test)]
mod tests {
    use super::super::lindhard_scharff::LindhardScharff;
    use super::*;

    #[test]
    fn single_element_equals_element_stopping_exactly() {
        let m = LindhardScharff::new();
        let ion = Ion::new(15).unwrap();
        let si = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        let want = m.stopping(&ion, 14, 3.0e4).unwrap();
        let got = bragg_cross_section_per_atom(&m, &NoCorrection, &ion, &si, 3.0e4).unwrap();
        assert_eq!(got.to_bits(), want.to_bits());
    }

    #[test]
    fn sio2_is_hand_weighted_sum() {
        let m = LindhardScharff::new();
        let ion = Ion::new(5).unwrap();
        let sio2 = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
        let e = 2.0e4;
        let (s_si, s_o) = (
            m.stopping(&ion, 14, e).unwrap(),
            m.stopping(&ion, 8, e).unwrap(),
        );
        let hand = (s_si + 2.0 * s_o) / 3.0;
        let got = bragg_cross_section_per_atom(&m, &NoCorrection, &ion, &sio2, e).unwrap();
        assert!((got / hand - 1.0).abs() < 1e-14);
        let p = bragg_stopping_power(&m, &NoCorrection, &ion, &sio2, e).unwrap();
        assert!((p / (sio2.atom_number_density() * hand) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn correction_hook_scales() {
        let m = LindhardScharff::new();
        let ion = Ion::new(5).unwrap();
        let sio2 = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
        let a = bragg_cross_section_per_atom(&m, &NoCorrection, &ion, &sio2, 1e4).unwrap();
        let b =
            bragg_cross_section_per_atom(&m, &ConstantCorrection(0.9), &ion, &sio2, 1e4).unwrap();
        assert!((b / a - 0.9).abs() < 1e-14);
    }

    #[test]
    fn rejects_bad_compound_correction() {
        let m = LindhardScharff::new();
        let ion = Ion::new(5).unwrap();
        let sio2 = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let r = bragg_cross_section_per_atom(&m, &ConstantCorrection(bad), &ion, &sio2, 1e4);
            assert!(matches!(
                r,
                Err(StoppingError::InvalidParameter {
                    name: "compound_correction",
                    ..
                })
            ));
            let r = bragg_stopping_power(&m, &ConstantCorrection(bad), &ion, &sio2, 1e4);
            assert!(r.is_err());
        }
        let z = bragg_cross_section_per_atom(&m, &ConstantCorrection(0.0), &ion, &sio2, 1e4);
        assert_eq!(z.unwrap(), 0.0);
    }

    #[test]
    fn bad_element_correction_propagates_through_bragg() {
        let m = LindhardScharff::new().with_correction(8, -1.0);
        let ion = Ion::new(5).unwrap();
        let sio2 = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
        assert!(matches!(
            bragg_cross_section_per_atom(&m, &NoCorrection, &ion, &sio2, 1e4),
            Err(StoppingError::InvalidParameter { .. })
        ));
    }
}
