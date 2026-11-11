//! SRIM-data-sourced electronic stopping from a versioned coefficient file.
//!
//! The model evaluates the published five-coefficient fit of the stopping
//! cross section per atom, with the coefficients read from a data file rather
//! than written into the source:
//!
//! * `E < 10 keV/u`: `S = A1 E^0.45`
//! * `E >= 10 keV/u`: `S_low = A2 E^0.45`,
//!   `S_high = (A3 / E) ln(1 + A4 / E + A5 E)`,
//!   `S = S_low S_high / (S_low + S_high)`
//!
//! with `E` the projectile energy per atomic mass unit in keV/u and `S` in
//! eV 1e-15 cm² per atom. Form: H. H. Andersen and J. F. Ziegler, *Hydrogen
//! Stopping Powers and Ranges in All Elements* (Pergamon, 1977); J. F.
//! Ziegler, J. P. Biersack and U. Littmark, *The Stopping and Range of Ions in
//! Solids* (Pergamon, 1985), ch. 3. The coefficient values are data from that
//! book (see `docs/data-provenance.md`); no SRIM code is read or ported
//! (`CONTRIBUTING.md`).
//!
//! **Verification note:** the functional form above is written from the
//! published description and must be checked against the book page by whoever
//! enters the full tables (issue #633) before those values are trusted; the
//! loader's format version is bumped if the form changes.
//!
//! # File format
//!
//! ```toml
//! format = "lindhard-srim-coefficients/1"   # required
//! name = "free text"                         # optional
//!
//! [[row]]
//! ion_z = 1
//! target_z = 14
//! a1 = 1.0
//! a2 = 1.0
//! a3 = 1.0
//! a4 = 1.0
//! a5 = 1.0
//! e_min_kev_per_u = 1.0     # validity range of the row, inclusive
//! e_max_kev_per_u = 1.0e4
//! citation = "book, edition, table, page"   # required, non-blank, per row
//! ```
//!
//! A row without a non-blank `citation` is rejected. A pair appears once.
//! A query outside a row's energy range, or for a pair with no row, is an
//! error ([`StoppingError::OutOfTableRange`],
//! [`StoppingError::NoSrimCoefficients`]); nothing is extrapolated.
//!
//! Output from this model is SRIM-data-sourced: agreement with SRIM is
//! consistency, not independent validation (`CLAUDE.md`).

use super::{
    check_energy, from_ev_1e15_cm2, ElectronicStopping, Ion, StoppingError, ValidityRange,
};
use serde::Deserialize;
use std::path::Path;

/// Value of the `format` key this loader reads.
pub const FORMAT_VERSION: &str = "lindhard-srim-coefficients/1";

/// Energy per amu below which the single-term low-energy form applies, keV/u.
pub const LOW_ENERGY_SPLIT_KEV_PER_U: f64 = 10.0;

const LOW_ENERGY_EXPONENT: f64 = 0.45;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawFile {
    format: Option<String>,
    #[serde(default)]
    name: Option<String>,
    #[serde(default)]
    row: Vec<RawRow>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawRow {
    ion_z: u8,
    target_z: u8,
    a1: f64,
    a2: f64,
    a3: f64,
    a4: f64,
    a5: f64,
    e_min_kev_per_u: f64,
    e_max_kev_per_u: f64,
    citation: Option<String>,
}

/// One validated coefficient row for an (ion, target element) pair.
#[derive(Debug, Clone, PartialEq)]
pub struct SrimRow {
    ion_z: u8,
    target_z: u8,
    a: [f64; 5],
    e_min_kev_per_u: f64,
    e_max_kev_per_u: f64,
    citation: String,
}

impl SrimRow {
    /// Projectile atomic number.
    pub fn ion_z(&self) -> u8 {
        self.ion_z
    }

    /// Target atomic number.
    pub fn target_z(&self) -> u8 {
        self.target_z
    }

    /// The row's citation (edition, table, page).
    pub fn citation(&self) -> &str {
        &self.citation
    }

    /// Stopping at `e` keV/u, eV 1e-15 cm² per atom (no range check).
    fn eval(&self, e: f64) -> f64 {
        let [a1, a2, a3, a4, a5] = self.a;
        if e < LOW_ENERGY_SPLIT_KEV_PER_U {
            a1 * e.powf(LOW_ENERGY_EXPONENT)
        } else {
            let low = a2 * e.powf(LOW_ENERGY_EXPONENT);
            let high = a3 / e * (1.0 + a4 / e + a5 * e).ln();
            low * high / (low + high)
        }
    }
}

/// The SRIM-data-sourced stopping model: a set of validated coefficient rows.
///
/// `Default` is an empty model, for which every query is an error.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SrimStopping {
    name: String,
    rows: Vec<SrimRow>,
}

impl SrimStopping {
    /// Parse and validate a coefficient file.
    pub fn from_toml_str(text: &str) -> Result<Self, StoppingError> {
        let raw: RawFile =
            toml::from_str(text).map_err(|e| StoppingError::InvalidTable(e.to_string()))?;
        match raw.format.as_deref() {
            Some(FORMAT_VERSION) => {}
            Some(other) => {
                return Err(StoppingError::InvalidTable(format!(
                    "unsupported format {other:?}; this loader reads {FORMAT_VERSION:?}"
                )))
            }
            None => {
                return Err(StoppingError::InvalidTable(format!(
                    "missing `format` key (expected {FORMAT_VERSION:?})"
                )))
            }
        }
        let mut rows: Vec<SrimRow> = Vec::with_capacity(raw.row.len());
        for (i, r) in raw.row.into_iter().enumerate() {
            let bad = |why: String| StoppingError::InvalidTable(format!("row {}: {why}", i + 1));
            let citation = r.citation.unwrap_or_default().trim().to_string();
            if citation.is_empty() {
                return Err(StoppingError::MissingProvenance);
            }
            super::target(r.ion_z)?;
            super::target(r.target_z)?;
            let a = [r.a1, r.a2, r.a3, r.a4, r.a5];
            if a.iter().any(|x| !x.is_finite()) {
                return Err(bad("coefficients must be finite".into()));
            }
            if !(r.e_min_kev_per_u.is_finite()
                && r.e_min_kev_per_u > 0.0
                && r.e_max_kev_per_u.is_finite()
                && r.e_max_kev_per_u > r.e_min_kev_per_u)
            {
                return Err(bad(
                    "energy range must satisfy 0 < e_min_kev_per_u < e_max_kev_per_u".into(),
                ));
            }
            if rows
                .iter()
                .any(|p| p.ion_z == r.ion_z && p.target_z == r.target_z)
            {
                return Err(bad(format!(
                    "duplicate pair Z1={}, Z2={}",
                    r.ion_z, r.target_z
                )));
            }
            let row = SrimRow {
                ion_z: r.ion_z,
                target_z: r.target_z,
                a,
                e_min_kev_per_u: r.e_min_kev_per_u,
                e_max_kev_per_u: r.e_max_kev_per_u,
                citation,
            };
            // The fit must give a finite, positive stopping everywhere in the
            // row's range; check it on a log grid including both ends.
            const N: i32 = 64;
            let (lo, hi) = (row.e_min_kev_per_u.ln(), row.e_max_kev_per_u.ln());
            for k in 0..=N {
                let e = (lo + (hi - lo) * f64::from(k) / f64::from(N)).exp();
                let s = row.eval(e);
                if !(s.is_finite() && s > 0.0) {
                    return Err(bad(format!(
                        "coefficients give a non-positive or non-finite stopping at {e} keV/u"
                    )));
                }
            }
            rows.push(row);
        }
        Ok(Self {
            name: raw.name.unwrap_or_default(),
            rows,
        })
    }

    /// Load from a TOML file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self, StoppingError> {
        let text = std::fs::read_to_string(path.as_ref())
            .map_err(|e| StoppingError::InvalidTable(e.to_string()))?;
        Self::from_toml_str(&text)
    }

    /// The file's `name`, empty if absent.
    pub fn dataset_name(&self) -> &str {
        &self.name
    }

    /// The validated rows, in file order.
    pub fn rows(&self) -> &[SrimRow] {
        &self.rows
    }

    /// Distinct citations in file order, for run metadata.
    pub fn citations(&self) -> Vec<&str> {
        let mut v: Vec<&str> = Vec::new();
        for r in &self.rows {
            if !v.contains(&r.citation.as_str()) {
                v.push(&r.citation);
            }
        }
        v
    }

    fn row(&self, ion_z: u8, target_z: u8) -> Option<&SrimRow> {
        self.rows
            .iter()
            .find(|r| r.ion_z == ion_z && r.target_z == target_z)
    }

    /// Whether a row exists for this pair.
    pub fn covers_pair(&self, ion_z: u8, target_z: u8) -> bool {
        self.row(ion_z, target_z).is_some()
    }
}

impl ElectronicStopping for SrimStopping {
    fn name(&self) -> &'static str {
        "srim"
    }

    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        check_energy(energy_ev)?;
        let row = self
            .row(ion.z(), target_z)
            .ok_or(StoppingError::NoSrimCoefficients {
                ion: ion.z(),
                target: target_z,
            })?;
        let mass = ion.mass_amu();
        let e = energy_ev / mass / 1.0e3;
        if e < row.e_min_kev_per_u || e > row.e_max_kev_per_u {
            return Err(StoppingError::OutOfTableRange {
                energy_ev,
                min_ev: row.e_min_kev_per_u * mass * 1.0e3,
                max_ev: row.e_max_kev_per_u * mass * 1.0e3,
            });
        }
        Ok(from_ev_1e15_cm2(row.eval(e)))
    }

    /// The union of the row ranges for this ion; `(0, inf)` when the file has
    /// no row for it (the query itself then fails).
    fn validity(&self, ion: &Ion) -> ValidityRange {
        let k = ion.mass_amu() * 1.0e3;
        let mut it = self.rows.iter().filter(|r| r.ion_z == ion.z());
        match it.next() {
            None => ValidityRange {
                min_energy_ev: 0.0,
                max_energy_ev: f64::INFINITY,
            },
            Some(first) => {
                let (lo, hi) = it.fold(
                    (first.e_min_kev_per_u, first.e_max_kev_per_u),
                    |(lo, hi), r| (lo.min(r.e_min_kev_per_u), hi.max(r.e_max_kev_per_u)),
                );
                ValidityRange {
                    min_energy_ev: lo * k,
                    max_energy_ev: hi * k,
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ion::stopping::to_ev_1e15_cm2;

    // Synthetic numbers made up for unit tests. They are NOT book values.
    const FIXTURE: &str = r#"
format = "lindhard-srim-coefficients/1"
name = "synthetic test fixture"

[[row]]
ion_z = 1
target_z = 14
a1 = 2.0
a2 = 3.0
a3 = 400.0
a4 = 5.0
a5 = 0.01
e_min_kev_per_u = 1.0
e_max_kev_per_u = 1.0e4
citation = "synthetic test data, made up for unit tests"

[[row]]
ion_z = 1
target_z = 8
a1 = 1.5
a2 = 2.5
a3 = 300.0
a4 = 4.0
a5 = 0.02
e_min_kev_per_u = 5.0
e_max_kev_per_u = 1.0e3
citation = "synthetic test data, made up for unit tests (second row)"
"#;

    fn proton() -> Ion {
        Ion::with_mass(1, 1.0).unwrap()
    }

    #[test]
    fn low_and_high_energy_branches_match_the_formulas() {
        let m = SrimStopping::from_toml_str(FIXTURE).unwrap();
        let at = |ev: f64| to_ev_1e15_cm2(m.stopping(&proton(), 14, ev).unwrap());
        // 4 keV/u: A1 E^0.45
        let e: f64 = 4.0;
        assert!((at(4.0e3) / (2.0 * e.powf(0.45)) - 1.0).abs() < 1e-12);
        // 100 keV/u: harmonic combination
        let e: f64 = 100.0;
        let low = 3.0 * e.powf(0.45);
        let high = 400.0 / e * (1.0 + 5.0 / e + 0.01 * e).ln();
        assert!((at(1.0e5) / (low * high / (low + high)) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn energy_is_per_amu() {
        let m = SrimStopping::from_toml_str(FIXTURE).unwrap();
        // A mass-2 ion at twice the energy is the same speed: same S.
        let d = Ion::with_mass(1, 2.0).unwrap();
        let a = m.stopping(&proton(), 14, 5.0e4).unwrap();
        let b = m.stopping(&d, 14, 1.0e5).unwrap();
        assert!((a / b - 1.0).abs() < 1e-12);
    }

    #[test]
    fn out_of_range_is_an_error_not_an_extrapolation() {
        let m = SrimStopping::from_toml_str(FIXTURE).unwrap();
        assert!(matches!(
            m.stopping(&proton(), 14, 500.0),
            Err(StoppingError::OutOfTableRange { .. })
        ));
        assert!(matches!(
            m.stopping(&proton(), 14, 2.0e7),
            Err(StoppingError::OutOfTableRange { .. })
        ));
        // inclusive ends
        assert!(m.stopping(&proton(), 14, 1.0e3).is_ok());
        assert!(m.stopping(&proton(), 14, 1.0e7).is_ok());
    }

    #[test]
    fn unknown_pair_is_an_error() {
        let m = SrimStopping::from_toml_str(FIXTURE).unwrap();
        assert!(matches!(
            m.stopping(&proton(), 29, 1.0e5),
            Err(StoppingError::NoSrimCoefficients { ion: 1, target: 29 })
        ));
        let si = Ion::new(14).unwrap();
        assert!(matches!(
            m.stopping(&si, 14, 1.0e5),
            Err(StoppingError::NoSrimCoefficients { .. })
        ));
        assert!(matches!(
            SrimStopping::default().stopping(&proton(), 14, 1.0e5),
            Err(StoppingError::NoSrimCoefficients { .. })
        ));
    }

    #[test]
    fn rejects_missing_format_blank_citation_duplicates_and_bad_range() {
        let no_format = FIXTURE.replace("format = \"lindhard-srim-coefficients/1\"\n", "");
        assert!(SrimStopping::from_toml_str(&no_format).is_err());
        let v2 = FIXTURE.replace("coefficients/1", "coefficients/2");
        assert!(SrimStopping::from_toml_str(&v2).is_err());
        let blank = FIXTURE.replace("synthetic test data, made up for unit tests\"", "  \"");
        assert_eq!(
            SrimStopping::from_toml_str(&blank).unwrap_err(),
            StoppingError::MissingProvenance
        );
        let dup = FIXTURE.replace("target_z = 8", "target_z = 14");
        assert!(SrimStopping::from_toml_str(&dup)
            .unwrap_err()
            .to_string()
            .contains("duplicate"));
        let range = FIXTURE.replace("e_min_kev_per_u = 5.0", "e_min_kev_per_u = 5000.0");
        assert!(SrimStopping::from_toml_str(&range).is_err());
        let neg = FIXTURE.replace("a1 = 2.0", "a1 = -2.0");
        assert!(SrimStopping::from_toml_str(&neg).is_err());
    }

    #[test]
    fn validity_is_the_union_of_the_ion_rows() {
        let m = SrimStopping::from_toml_str(FIXTURE).unwrap();
        let v = m.validity(&proton());
        assert!((v.min_energy_ev - 1.0e3).abs() < 1e-6);
        assert!((v.max_energy_ev - 1.0e7).abs() < 1e-3);
        assert_eq!(m.citations().len(), 2);
    }
}
