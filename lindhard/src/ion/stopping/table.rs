//! User-supplied `S_e(E)` tables, loaded from TOML, with mandatory provenance.
//!
//! ```toml
//! provenance = "Author, Journal vol, page (year), Table N"   # required
//! ion_z = 14
//! ion_mass_amu = 28.0855            # optional; default: standard weight
//! target_z = 14
//! energy_ev = [1.0e3, 1.0e4, 1.0e5]
//! stopping_ev_1e15_cm2 = [10.0, 30.0, 60.0]   # eV 1e-15 cm^2 per atom
//! ```
//!
//! A table without a non-empty `provenance` is rejected at load time
//! (see `docs/data-provenance.md`: data without an origin is not admitted).
//! Do not load SRIM- or ICRU-derived tables for anything committed to the
//! repository.
//!
//! Interpolation is piecewise linear in `ln S` versus `ln E`, which preserves
//! the monotonicity of each segment of the data. **Out of range** queries are
//! an error ([`StoppingError::OutOfTableRange`]); nothing is extrapolated or
//! clamped silently.

use super::{from_ev_1e15_cm2, ElectronicStopping, Ion, StoppingError, ValidityRange};
use serde::Deserialize;
use std::path::Path;

#[derive(Deserialize)]
struct RawTable {
    provenance: Option<String>,
    ion_z: u8,
    ion_mass_amu: Option<f64>,
    target_z: u8,
    energy_ev: Vec<f64>,
    stopping_ev_1e15_cm2: Vec<f64>,
}

/// A validated user stopping table for one ion/target pair.
#[derive(Debug, Clone, PartialEq)]
pub struct StoppingTable {
    provenance: String,
    ion_z: u8,
    target_z: u8,
    /// ln(E / eV)
    ln_e: Vec<f64>,
    /// ln(S / (J m^2))
    ln_s: Vec<f64>,
    energy_ev: Vec<f64>,
}

impl StoppingTable {
    /// Parse and validate a TOML table.
    pub fn from_toml_str(text: &str) -> Result<Self, StoppingError> {
        let raw: RawTable =
            toml::from_str(text).map_err(|e| StoppingError::InvalidTable(e.to_string()))?;
        let provenance = raw.provenance.unwrap_or_default().trim().to_string();
        if provenance.is_empty() {
            return Err(StoppingError::MissingProvenance);
        }
        let _ = raw.ion_mass_amu; // identity is by Z; mass is not needed for S_e(E) lookup
        super::target(raw.ion_z)?;
        super::target(raw.target_z)?;
        let (e, s) = (raw.energy_ev, raw.stopping_ev_1e15_cm2);
        if e.len() != s.len() || e.len() < 2 {
            return Err(StoppingError::InvalidTable(
                "energy and stopping need the same length, at least 2 points".into(),
            ));
        }
        if e.iter().any(|x| !(x.is_finite() && *x > 0.0))
            || s.iter().any(|x| !(x.is_finite() && *x > 0.0))
        {
            return Err(StoppingError::InvalidTable(
                "energies and stopping values must be finite and positive".into(),
            ));
        }
        if e.windows(2).any(|w| w[1] <= w[0]) {
            return Err(StoppingError::InvalidTable(
                "energies must be strictly increasing".into(),
            ));
        }
        Ok(Self {
            provenance,
            ion_z: raw.ion_z,
            target_z: raw.target_z,
            ln_e: e.iter().map(|x| x.ln()).collect(),
            ln_s: s.iter().map(|x| from_ev_1e15_cm2(*x).ln()).collect(),
            energy_ev: e,
        })
    }

    /// Load from a TOML file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self, StoppingError> {
        let text = std::fs::read_to_string(path.as_ref())
            .map_err(|e| StoppingError::InvalidTable(e.to_string()))?;
        Self::from_toml_str(&text)
    }

    /// The provenance string (to be copied into run metadata).
    pub fn provenance(&self) -> &str {
        &self.provenance
    }
}

impl ElectronicStopping for StoppingTable {
    fn name(&self) -> &'static str {
        "user-table"
    }

    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        super::check_energy(energy_ev)?;
        if ion.z != self.ion_z || target_z != self.target_z {
            return Err(StoppingError::TableMismatch {
                table_ion: self.ion_z,
                table_target: self.target_z,
                ion: ion.z,
                target: target_z,
            });
        }
        let (lo, hi) = (self.energy_ev[0], self.energy_ev[self.energy_ev.len() - 1]);
        if energy_ev < lo || energy_ev > hi {
            return Err(StoppingError::OutOfTableRange {
                energy_ev,
                min_ev: lo,
                max_ev: hi,
            });
        }
        let x = energy_ev.ln();
        let i = self
            .ln_e
            .partition_point(|&v| v <= x)
            .clamp(1, self.ln_e.len() - 1);
        let t = (x - self.ln_e[i - 1]) / (self.ln_e[i] - self.ln_e[i - 1]);
        Ok((self.ln_s[i - 1] + t * (self.ln_s[i] - self.ln_s[i - 1])).exp())
    }

    fn validity(&self, _ion: &Ion) -> ValidityRange {
        ValidityRange {
            min_energy_ev: self.energy_ev[0],
            max_energy_ev: self.energy_ev[self.energy_ev.len() - 1],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GOOD: &str = r#"
provenance = "synthetic test data, made up for unit tests"
ion_z = 14
target_z = 14
energy_ev = [1.0e3, 1.0e4, 1.0e5]
stopping_ev_1e15_cm2 = [10.0, 40.0, 90.0]
"#;

    #[test]
    fn rejects_missing_or_blank_provenance() {
        let no = GOOD.replace(
            "provenance = \"synthetic test data, made up for unit tests\"\n",
            "",
        );
        let msg = StoppingTable::from_toml_str(&no).unwrap_err().to_string();
        assert!(msg.contains("provenance"), "{msg}");
        let blank = GOOD.replace("synthetic test data, made up for unit tests", "  ");
        assert_eq!(
            StoppingTable::from_toml_str(&blank).unwrap_err(),
            StoppingError::MissingProvenance
        );
    }

    #[test]
    fn interpolates_exactly_at_nodes_and_log_log_between() {
        let t = StoppingTable::from_toml_str(GOOD).unwrap();
        let ion = Ion::new(14).unwrap();
        let at = |e| super::super::to_ev_1e15_cm2(t.stopping(&ion, 14, e).unwrap());
        assert!((at(1.0e4) / 40.0 - 1.0).abs() < 1e-12);
        assert!((at(1.0e5) / 90.0 - 1.0).abs() < 1e-12);
        let mid = at(10f64.powf(3.5));
        assert!((mid / (10.0f64 * 40.0).sqrt() - 1.0).abs() < 1e-12);
        let mut prev = 0.0;
        for i in 0..200 {
            let v = at(1.0e3 * 10f64.powf(f64::from(i) / 100.0));
            assert!(v > prev);
            prev = v;
        }
    }

    #[test]
    fn out_of_range_and_mismatch_are_errors() {
        let t = StoppingTable::from_toml_str(GOOD).unwrap();
        let ion = Ion::new(14).unwrap();
        assert!(matches!(
            t.stopping(&ion, 14, 10.0),
            Err(StoppingError::OutOfTableRange { .. })
        ));
        assert!(matches!(
            t.stopping(&ion, 8, 1.0e4),
            Err(StoppingError::TableMismatch { .. })
        ));
    }

    #[test]
    fn rejects_unsorted() {
        let bad = GOOD.replace("[1.0e3, 1.0e4, 1.0e5]", "[1.0e4, 1.0e3, 1.0e5]");
        assert!(StoppingTable::from_toml_str(&bad).is_err());
    }
}
