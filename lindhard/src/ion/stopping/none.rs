//! No electronic stopping: the nuclear-only limit.
//!
//! [`NoStopping`] returns zero for every ion, target and energy, so a run
//! loses energy only in nuclear collisions. This is not a physical model of
//! any target. It is the `k = 0` limit of the classical range theory, where
//! the electronic stopping coefficient vanishes and only nuclear stopping
//! remains (J. Lindhard, M. Scharff and H. E. Schiott, Mat. Fys. Medd. Dan.
//! Vid. Selsk. 33 (14) (1963)), and it serves two purposes here:
//!
//! * like-for-like transport comparisons against a code run with its
//!   electronic stopping switched off, and
//! * nuclear-only checks (range and damage limits) that isolate the
//!   scattering and cascade machinery from any stopping model.
//!
//! It has no energy range of its own: [`ElectronicStopping::validity`] is
//! unbounded, so no out-of-range warning is ever raised for it.

use super::{ElectronicStopping, Ion, StoppingError, ValidityRange};

/// Zero electronic stopping (name `"none"`). Stateless, so bitwise
/// reproducible across threads like every model in this module.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NoStopping;

impl NoStopping {
    /// The (only) zero-stopping model.
    pub fn new() -> Self {
        Self
    }
}

impl ElectronicStopping for NoStopping {
    fn name(&self) -> &'static str {
        "none"
    }

    /// Always `0.0` J m². The arguments are not checked: there is no formula
    /// to be outside of.
    fn stopping(&self, _ion: &Ion, _target_z: u8, _energy_ev: f64) -> Result<f64, StoppingError> {
        Ok(0.0)
    }

    /// Unbounded: `[0, inf)` eV.
    fn validity(&self, _ion: &Ion) -> ValidityRange {
        ValidityRange {
            min_energy_ev: 0.0,
            max_energy_ev: f64::INFINITY,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_everywhere_and_unbounded() {
        let m = NoStopping::new();
        for z1 in [1, 5, 14, 33, 92] {
            let ion = Ion::new(z1).unwrap();
            for e in [1.0e-3, 1.0, 5.0e3, 1.0e9] {
                assert_eq!(m.stopping(&ion, 14, e).unwrap().to_bits(), 0.0f64.to_bits());
                assert!(m.validity(&ion).contains(e));
            }
        }
        assert_eq!(m.name(), "none");
    }
}
