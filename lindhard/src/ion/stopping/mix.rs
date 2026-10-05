//! Equipartition mix of Lindhard-Scharff (nonlocal) and Oen-Robinson (local).
//!
//! Each collision loses half of its electronic energy by the nonlocal LS
//! (continuous, along the free flight) channel and half by the local OR
//! channel (at the collision, depending on `r_min`). This "equipartition"
//! split is the combination described in the BCA literature that builds on
//! Oen and Robinson (1976). Averaged over impact parameter the mix equals the
//! LS stopping.

use super::lindhard_scharff::LindhardScharff;
use super::oen_robinson::OenRobinson;
use super::{ElectronicStopping, Ion, StoppingError, ValidityRange};

/// Equipartition LS/OR model.
#[derive(Debug, Clone, Default)]
pub struct EquipartitionMix {
    ls: LindhardScharff,
    or: OenRobinson,
}

impl EquipartitionMix {
    /// New mix with plain LS coefficient.
    pub fn new() -> Self {
        Self::default()
    }

    /// Nonlocal (free-flight) loss per unit areal density, J m²: half of LS.
    pub fn nonlocal_stopping(
        &self,
        ion: &Ion,
        target_z: u8,
        energy_ev: f64,
    ) -> Result<f64, StoppingError> {
        Ok(0.5 * self.ls.stopping(ion, target_z, energy_ev)?)
    }

    /// Local loss at one collision, J: half of the OR local loss.
    pub fn local_loss(
        &self,
        ion: &Ion,
        target_z: u8,
        energy_ev: f64,
        r_min_m: f64,
    ) -> Result<f64, StoppingError> {
        Ok(0.5 * self.or.local_loss(ion, target_z, energy_ev, r_min_m)?)
    }
}

impl ElectronicStopping for EquipartitionMix {
    fn name(&self) -> &'static str {
        "equipartition-ls-or"
    }

    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        // Half nonlocal + half local averaged over impact parameter.
        Ok(self.nonlocal_stopping(ion, target_z, energy_ev)?
            + 0.5 * self.or.stopping(ion, target_z, energy_ev)?)
    }

    fn validity(&self, ion: &Ion) -> ValidityRange {
        self.ls.validity(ion)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn average_equals_ls() {
        let ion = Ion::new(31).unwrap();
        let a = EquipartitionMix::new().stopping(&ion, 14, 5e4).unwrap();
        let b = LindhardScharff::new().stopping(&ion, 14, 5e4).unwrap();
        assert!((a / b - 1.0).abs() < 1e-14);
    }
}
