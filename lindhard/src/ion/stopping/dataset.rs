//! Stopping datasets: measured `S_e` points with a citation on every point, and
//! the comparison of an [`ElectronicStopping`] model against them.
//!
//! The format is specified in `docs/stopping-data-format.md`. In short: a
//! `#`-prefixed `key: value` header, one CSV header row with the units in the
//! column names, then one point per line:
//!
//! ```text
//! # format: lindhard-stopping-data/1
//! # name: example
//! # description: free text
//! ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method,citation
//! 5,,Si,1.0e4,12.5,0.3,transmission,"Author, Journal vol, page (year), Table N"
//! 5,11.009,SiO2,2.0e4,10.1,0.2,backscattering,"Author, Journal vol, page (year)"
//! ```
//!
//! * Stopping and its uncertainty are cross sections per target atom in
//!   eV·10⁻¹⁵ cm² (the unit of [`to_ev_1e15_cm2`]). For a compound target the
//!   value is per **average atom** (per molecule divided by the number of atoms
//!   in the formula), the quantity [`bragg_cross_section_per_atom`] returns.
//! * Every row needs a non-empty `method` and `citation`; a row without one is
//!   rejected with its line number (see `docs/data-provenance.md`: data without
//!   an origin is not admitted).
//! * The parser is a small hand-written CSV reader: comma-separated fields,
//!   double quotes around a field that contains a comma, `""` for a literal
//!   quote inside one. Quoted fields do not span lines.
//!
//! No real measurement is shipped with this module: which published data may
//! enter the tree is decided separately (issue #34). Tests use synthetic points
//! generated from this repository's own models.
//!
//! Every result is a pure function of the model and the dataset, computed in
//! dataset order with sequential summation, so it is bitwise reproducible
//! whatever the thread count of the caller.

use super::bragg::{bragg_cross_section_per_atom, CompoundCorrection};
use super::{to_ev_1e15_cm2, ElectronicStopping, Ion, StoppingError};
use crate::elements::{element, element_by_symbol};
use crate::material::Material;
use std::path::Path;

/// Value of the `format` header key this loader reads.
pub const FORMAT_VERSION: &str = "lindhard-stopping-data/1";

/// The column header row, in order. The units are part of the names
/// (`_ev`, `_amu`, `_ev_1e15_cm2`) so a file in another unit cannot be read
/// silently.
pub const COLUMNS: [&str; 8] = [
    "ion_z",
    "ion_mass_amu",
    "target",
    "energy_ev",
    "s_e_ev_1e15_cm2",
    "sigma_ev_1e15_cm2",
    "method",
    "citation",
];

/// Mass density given to the [`Material`] built for a compound target,
/// kg/m³. It is a placeholder: under Bragg's rule (W. H. Bragg and R. Kleeman,
/// Phil. Mag. 10, 318 (1905)) the stopping cross section per atom does not
/// depend on the density, and a dataset does not record one. A
/// [`CompoundCorrection`] used with [`compare`] must therefore not depend on
/// [`Material::mass_density`]; it can identify the compound by
/// [`Material::name`], which is set to the formula as written in the file.
pub const PLACEHOLDER_DENSITY_KG_M3: f64 = 1000.0;

/// The target of a point: one element, or a compound given by its formula.
#[derive(Debug, Clone, PartialEq)]
pub enum Target {
    /// An element, by atomic number. A formula with a single element (`Si`,
    /// `O2`) is an element target.
    Element(u8),
    /// A compound: the formula as written and a material built from it
    /// (atom fractions from the formula counts, named with the formula,
    /// density [`PLACEHOLDER_DENSITY_KG_M3`]).
    Compound {
        /// The formula as written in the file, e.g. `SiO2`.
        formula: String,
        /// The material handed to [`bragg_cross_section_per_atom`].
        material: Material,
    },
}

impl Target {
    /// Parse an element symbol or a formula such as `SiO2`, `Al2O3` or
    /// `Si0.5Ge0.5`: element symbols (case-sensitive), each followed by an
    /// optional positive count. Parentheses and repeated elements are not
    /// accepted.
    pub fn parse(text: &str) -> Result<Self, String> {
        let parts = parse_formula(text)?;
        if parts.len() == 1 {
            return Ok(Self::Element(parts[0].0));
        }
        let material = Material::from_atom_fractions(&parts, Some(PLACEHOLDER_DENSITY_KG_M3))
            .map_err(|e| format!("target {text:?}: {e}"))?
            .with_name(text);
        Ok(Self::Compound {
            formula: text.to_string(),
            material,
        })
    }

    /// The target as written for a compound, or the element symbol.
    pub fn label(&self) -> &str {
        match self {
            // `parse` only builds elements that exist.
            Self::Element(z) => element(*z).map_or("?", |e| e.symbol),
            Self::Compound { formula, .. } => formula,
        }
    }
}

fn parse_formula(text: &str) -> Result<Vec<(u8, f64)>, String> {
    if text.is_empty() {
        return Err("empty target".into());
    }
    let b = text.as_bytes();
    let mut i = 0;
    let mut parts: Vec<(u8, f64)> = Vec::new();
    while i < b.len() {
        if !b[i].is_ascii_uppercase() {
            return Err(format!(
                "target {text:?}: expected an element symbol at {:?}",
                &text[i..]
            ));
        }
        let start = i;
        i += 1;
        while i < b.len() && b[i].is_ascii_lowercase() {
            i += 1;
        }
        let symbol = &text[start..i];
        let el = element_by_symbol(symbol)
            .ok_or_else(|| format!("target {text:?}: unknown element {symbol:?}"))?;
        let cstart = i;
        while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.') {
            i += 1;
        }
        let count = if cstart == i {
            1.0
        } else {
            text[cstart..i]
                .parse::<f64>()
                .map_err(|_| format!("target {text:?}: bad count {:?}", &text[cstart..i]))?
        };
        if !(count.is_finite() && count > 0.0) {
            return Err(format!(
                "target {text:?}: count of {symbol} must be positive"
            ));
        }
        if parts.iter().any(|&(z, _)| z == el.z) {
            return Err(format!("target {text:?}: {symbol} appears more than once"));
        }
        parts.push((el.z, count));
    }
    Ok(parts)
}

/// One measured stopping point.
#[derive(Debug, Clone, PartialEq)]
pub struct StoppingPoint {
    line: usize,
    ion: Ion,
    target: Target,
    energy_ev: f64,
    s_e_ev_1e15_cm2: f64,
    sigma_ev_1e15_cm2: f64,
    method: String,
    citation: String,
}

impl StoppingPoint {
    /// 1-based line number in the source text.
    pub fn line(&self) -> usize {
        self.line
    }

    /// The projectile (with the mass given in the row, or the standard atomic
    /// weight when the column is empty).
    pub fn ion(&self) -> &Ion {
        &self.ion
    }

    /// The target.
    pub fn target(&self) -> &Target {
        &self.target
    }

    /// Ion kinetic energy, eV.
    pub fn energy_ev(&self) -> f64 {
        self.energy_ev
    }

    /// Measured stopping cross section per (average) target atom,
    /// eV·10⁻¹⁵ cm².
    pub fn s_e_ev_1e15_cm2(&self) -> f64 {
        self.s_e_ev_1e15_cm2
    }

    /// Stated one-standard-deviation uncertainty of
    /// [`StoppingPoint::s_e_ev_1e15_cm2`], same unit; 0 if none was stated.
    pub fn sigma_ev_1e15_cm2(&self) -> f64 {
        self.sigma_ev_1e15_cm2
    }

    /// Measurement method (free text).
    pub fn method(&self) -> &str {
        &self.method
    }

    /// Citation of the point's source (free text: a paper, or a database name
    /// and record identifier).
    pub fn citation(&self) -> &str {
        &self.citation
    }
}

/// A validated stopping dataset.
#[derive(Debug, Clone, PartialEq)]
pub struct StoppingDataset {
    name: String,
    description: Option<String>,
    metadata: Vec<(String, String)>,
    points: Vec<StoppingPoint>,
}

fn bad(line: usize, reason: impl Into<String>) -> StoppingError {
    StoppingError::InvalidDataset {
        line,
        reason: reason.into(),
    }
}

/// A header line `# key: value` with a lowercase identifier as the key.
fn metadata_entry(comment: &str) -> Option<(&str, &str)> {
    let (key, value) = comment.split_once(':')?;
    let key = key.trim();
    let is_key = !key.is_empty()
        && key
            .bytes()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == b'_');
    is_key.then(|| (key, value.trim()))
}

/// Split one CSV line into trimmed fields (see the module docs for the
/// quoting rules).
fn split_csv(line: &str) -> Result<Vec<String>, String> {
    let mut fields = Vec::new();
    let mut chars = line.chars().peekable();
    loop {
        while chars.peek().is_some_and(|c| *c == ' ' || *c == '\t') {
            chars.next();
        }
        let mut field = String::new();
        if chars.peek() == Some(&'"') {
            chars.next();
            loop {
                match chars.next() {
                    Some('"') if chars.peek() == Some(&'"') => {
                        chars.next();
                        field.push('"');
                    }
                    Some('"') => break,
                    Some(c) => field.push(c),
                    None => return Err("unterminated quoted field".into()),
                }
            }
            while chars.peek().is_some_and(|c| *c == ' ' || *c == '\t') {
                chars.next();
            }
            match chars.next() {
                None => {
                    fields.push(field);
                    return Ok(fields);
                }
                Some(',') => fields.push(field),
                Some(c) => return Err(format!("unexpected {c:?} after a quoted field")),
            }
        } else {
            loop {
                match chars.next() {
                    None => {
                        fields.push(field.trim().to_string());
                        return Ok(fields);
                    }
                    Some(',') => break,
                    Some('"') => return Err("quote inside an unquoted field".into()),
                    Some(c) => field.push(c),
                }
            }
            fields.push(field.trim().to_string());
        }
    }
}

fn number(line: usize, column: &str, text: &str) -> Result<f64, StoppingError> {
    text.parse::<f64>()
        .map_err(|_| bad(line, format!("{column}: not a number: {text:?}")))
}

fn finite_positive(line: usize, column: &str, text: &str) -> Result<f64, StoppingError> {
    let v = number(line, column, text)?;
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(bad(
            line,
            format!("{column} must be finite and positive, got {v}"),
        ))
    }
}

fn parse_point(line: usize, f: &[String]) -> Result<StoppingPoint, StoppingError> {
    if f.len() != COLUMNS.len() {
        return Err(bad(
            line,
            format!("expected {} fields, found {}", COLUMNS.len(), f.len()),
        ));
    }
    let ion_z: u8 = f[0]
        .parse()
        .map_err(|_| bad(line, format!("ion_z: not an atomic number: {:?}", f[0])))?;
    if element(ion_z).is_none() {
        return Err(bad(line, format!("unknown ion Z = {ion_z}")));
    }
    let ion = if f[1].is_empty() {
        Ion::new(ion_z)
    } else {
        Ion::with_mass(ion_z, number(line, "ion_mass_amu", &f[1])?)
    }
    .map_err(|e| bad(line, e.to_string()))?;
    let target = Target::parse(&f[2]).map_err(|e| bad(line, e))?;
    let energy_ev = finite_positive(line, "energy_ev", &f[3])?;
    let s_e = finite_positive(line, "s_e_ev_1e15_cm2", &f[4])?;
    let sigma = number(line, "sigma_ev_1e15_cm2", &f[5])?;
    if !(sigma.is_finite() && sigma >= 0.0) {
        return Err(bad(
            line,
            format!("sigma_ev_1e15_cm2 must be finite and non-negative, got {sigma}"),
        ));
    }
    if f[6].trim().is_empty() {
        return Err(bad(line, "empty method; every point needs one"));
    }
    if f[7].trim().is_empty() {
        return Err(bad(line, "empty citation; every point needs one"));
    }
    Ok(StoppingPoint {
        line,
        ion,
        target,
        energy_ev,
        s_e_ev_1e15_cm2: s_e,
        sigma_ev_1e15_cm2: sigma,
        method: f[6].trim().to_string(),
        citation: f[7].trim().to_string(),
    })
}

impl StoppingDataset {
    /// Parse and validate a dataset. A byte-order mark, CRLF line endings,
    /// blank lines and `#` comment lines in the data block are accepted.
    ///
    /// Errors name the offending line ([`StoppingError::InvalidDataset`]);
    /// a file with no column header or no points is
    /// [`StoppingError::InvalidTable`].
    pub fn from_csv_str(text: &str) -> Result<Self, StoppingError> {
        let text = text.strip_prefix('\u{feff}').unwrap_or(text);
        let mut metadata: Vec<(String, String)> = Vec::new();
        let mut in_header = true;
        let mut points = Vec::new();
        for (i, raw) in text.lines().enumerate() {
            let line = i + 1;
            let l = raw.strip_suffix('\r').unwrap_or(raw).trim();
            if l.is_empty() {
                continue;
            }
            if let Some(comment) = l.strip_prefix('#') {
                if in_header {
                    if let Some((key, value)) = metadata_entry(comment) {
                        if metadata.iter().any(|(k, _)| k == key) {
                            return Err(bad(line, format!("header key {key:?} given twice")));
                        }
                        metadata.push((key.to_string(), value.to_string()));
                    }
                }
                continue;
            }
            let fields = split_csv(l).map_err(|e| bad(line, e))?;
            if in_header {
                check_columns(line, &fields)?;
                check_metadata(line, &metadata)?;
                in_header = false;
            } else {
                points.push(parse_point(line, &fields)?);
            }
        }
        if in_header {
            return Err(StoppingError::InvalidTable(
                "stopping dataset has no column header row".into(),
            ));
        }
        if points.is_empty() {
            return Err(StoppingError::InvalidTable(
                "stopping dataset has no points".into(),
            ));
        }
        let get = |key: &str| {
            metadata
                .iter()
                .find(|(k, _)| k == key)
                .map(|(_, v)| v.clone())
        };
        Ok(Self {
            name: get("name").unwrap_or_default(),
            description: get("description").filter(|d| !d.is_empty()),
            metadata,
            points,
        })
    }

    /// Load from a file.
    pub fn from_csv_file(path: impl AsRef<Path>) -> Result<Self, StoppingError> {
        let text = std::fs::read_to_string(path.as_ref())
            .map_err(|e| StoppingError::InvalidTable(e.to_string()))?;
        Self::from_csv_str(&text)
    }

    /// The `name` header value.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The `description` header value, if given and non-empty.
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Every header `key: value` pair, in file order.
    pub fn metadata(&self) -> &[(String, String)] {
        &self.metadata
    }

    /// The points, in file order.
    pub fn points(&self) -> &[StoppingPoint] {
        &self.points
    }
}

fn check_columns(line: usize, fields: &[String]) -> Result<(), StoppingError> {
    for (k, want) in COLUMNS.iter().enumerate() {
        match fields.get(k) {
            Some(got) if got == want => {}
            Some(got) => {
                return Err(bad(
                    line,
                    format!("column {} must be {want:?}, found {got:?}", k + 1),
                ))
            }
            None => return Err(bad(line, format!("missing column {want:?}"))),
        }
    }
    if fields.len() > COLUMNS.len() {
        return Err(bad(
            line,
            format!("unexpected extra column {:?}", fields[COLUMNS.len()]),
        ));
    }
    Ok(())
}

fn check_metadata(line: usize, metadata: &[(String, String)]) -> Result<(), StoppingError> {
    let get = |key: &str| metadata.iter().find(|(k, _)| k == key).map(|(_, v)| v);
    match get("format") {
        Some(v) if v == FORMAT_VERSION => {}
        Some(v) => {
            return Err(bad(
                line,
                format!("unsupported format {v:?}; expected {FORMAT_VERSION:?}"),
            ))
        }
        None => return Err(bad(line, "header has no `format` key before the columns")),
    }
    match get("name") {
        Some(v) if !v.is_empty() => Ok(()),
        _ => Err(bad(
            line,
            "header has no non-empty `name` key before the columns",
        )),
    }
}

/// Model-versus-measurement residual at one point.
///
/// Sign convention: positive means the model is **above** the measurement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointResidual {
    /// Index of the point in [`StoppingDataset::points`].
    pub index: usize,
    /// Source line of the point.
    pub line: usize,
    /// Model stopping per (average) atom, eV·10⁻¹⁵ cm².
    pub model_ev_1e15_cm2: f64,
    /// Measured stopping, eV·10⁻¹⁵ cm².
    pub measured_ev_1e15_cm2: f64,
    /// `model / measured − 1`.
    pub relative: f64,
    /// `model − measured`, eV·10⁻¹⁵ cm².
    pub absolute_ev_1e15_cm2: f64,
    /// `(model − measured) / sigma`, or `None` when no uncertainty was stated
    /// (sigma 0).
    pub pull: Option<f64>,
}

/// A point the model could not evaluate, with the model's error. Such points
/// are not in any statistic.
#[derive(Debug, Clone, PartialEq)]
pub struct SkippedPoint {
    /// Index of the point in [`StoppingDataset::points`].
    pub index: usize,
    /// Source line of the point.
    pub line: usize,
    /// Why the model returned no value (for example
    /// [`StoppingError::NotApplicable`]).
    pub reason: StoppingError,
}

/// Summary statistics of the relative residuals of a set of points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResidualSummary {
    /// Number of points evaluated.
    pub n: usize,
    /// Number of points skipped (model error).
    pub n_skipped: usize,
    /// Mean of [`PointResidual::relative`]; NaN when `n == 0`.
    pub mean_relative: f64,
    /// Root mean square of [`PointResidual::relative`]; NaN when `n == 0`.
    pub rms_relative: f64,
    /// Largest `|relative|`; NaN when `n == 0`.
    pub max_abs_relative: f64,
}

impl ResidualSummary {
    /// Statistics over `relative` in the order given (sequential sums, so the
    /// result does not depend on how the caller was parallelised).
    fn of(relative: &[f64], n_skipped: usize) -> Self {
        let n = relative.len();
        if n == 0 {
            return Self {
                n,
                n_skipped,
                mean_relative: f64::NAN,
                rms_relative: f64::NAN,
                max_abs_relative: f64::NAN,
            };
        }
        let (mut sum, mut sum_sq, mut max_abs) = (0.0, 0.0, 0.0_f64);
        for &r in relative {
            sum += r;
            sum_sq += r * r;
            max_abs = max_abs.max(r.abs());
        }
        let nf = n as f64;
        Self {
            n,
            n_skipped,
            mean_relative: sum / nf,
            rms_relative: (sum_sq / nf).sqrt(),
            max_abs_relative: max_abs,
        }
    }
}

/// Statistics for one ion–target system: same ion Z, same ion mass, same
/// target as written.
#[derive(Debug, Clone, PartialEq)]
pub struct SystemSummary {
    /// Ion atomic number.
    pub ion_z: u8,
    /// Ion mass, u.
    pub ion_mass_amu: f64,
    /// Target label ([`Target::label`]).
    pub target: String,
    /// Statistics over this system's points.
    pub summary: ResidualSummary,
}

/// Result of [`compare`].
#[derive(Debug, Clone, PartialEq)]
pub struct Comparison {
    /// The model's [`ElectronicStopping::name`].
    pub model: &'static str,
    /// One entry per evaluated point, in dataset order.
    pub residuals: Vec<PointResidual>,
    /// Points the model could not evaluate, in dataset order.
    pub skipped: Vec<SkippedPoint>,
    /// Statistics over all evaluated points.
    pub overall: ResidualSummary,
    /// Statistics per ion–target system, in order of first appearance.
    pub systems: Vec<SystemSummary>,
}

/// Ion Z, ion mass bits, target label.
type SystemKey = (u8, u64, String);

/// Compare `model` against every point of `dataset`.
///
/// An element target is evaluated with [`ElectronicStopping::stopping`]; a
/// compound with [`bragg_cross_section_per_atom`] and `correction` (which is
/// not applied to element targets). Residuals follow the sign convention of
/// [`PointResidual`]: positive when the model is above the measurement. A
/// point where the model returns an error (such as
/// [`StoppingError::NotApplicable`]) is listed in [`Comparison::skipped`]
/// with that error, never dropped silently.
pub fn compare(
    model: &dyn ElectronicStopping,
    correction: &dyn CompoundCorrection,
    dataset: &StoppingDataset,
) -> Comparison {
    let mut residuals = Vec::with_capacity(dataset.points.len());
    let mut skipped = Vec::new();
    // Per system: key, relative residuals, skipped count.
    let mut systems: Vec<(SystemKey, Vec<f64>, usize)> = Vec::new();
    for (index, p) in dataset.points.iter().enumerate() {
        let key = (
            p.ion.z(),
            p.ion.mass_amu().to_bits(),
            p.target.label().to_string(),
        );
        let slot = match systems.iter().position(|s| s.0 == key) {
            Some(k) => k,
            None => {
                systems.push((key, Vec::new(), 0));
                systems.len() - 1
            }
        };
        let s = match &p.target {
            Target::Element(z) => model.stopping(&p.ion, *z, p.energy_ev),
            Target::Compound { material, .. } => {
                bragg_cross_section_per_atom(model, correction, &p.ion, material, p.energy_ev)
            }
        };
        match s {
            Ok(s) => {
                let m = to_ev_1e15_cm2(s);
                let meas = p.s_e_ev_1e15_cm2;
                let relative = m / meas - 1.0;
                let absolute = m - meas;
                residuals.push(PointResidual {
                    index,
                    line: p.line,
                    model_ev_1e15_cm2: m,
                    measured_ev_1e15_cm2: meas,
                    relative,
                    absolute_ev_1e15_cm2: absolute,
                    pull: (p.sigma_ev_1e15_cm2 > 0.0).then(|| absolute / p.sigma_ev_1e15_cm2),
                });
                systems[slot].1.push(relative);
            }
            Err(reason) => {
                skipped.push(SkippedPoint {
                    index,
                    line: p.line,
                    reason,
                });
                systems[slot].2 += 1;
            }
        }
    }
    let all: Vec<f64> = residuals.iter().map(|r| r.relative).collect();
    Comparison {
        model: model.name(),
        overall: ResidualSummary::of(&all, skipped.len()),
        systems: systems
            .into_iter()
            .map(
                |((ion_z, mass_bits, target), rel, n_skipped)| SystemSummary {
                    ion_z,
                    ion_mass_amu: f64::from_bits(mass_bits),
                    target,
                    summary: ResidualSummary::of(&rel, n_skipped),
                },
            )
            .collect(),
        residuals,
        skipped,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_quoting() {
        let f = split_csv(r#"1, "a, b" ,"say ""hi""",,x"#).unwrap();
        assert_eq!(f, ["1", "a, b", "say \"hi\"", "", "x"]);
        assert!(split_csv(r#"1,"open"#).is_err());
        assert!(split_csv(r#"1,"a"b"#).is_err());
        assert!(split_csv(r#"1,a"b"#).is_err());
    }

    #[test]
    fn formulas() {
        assert_eq!(Target::parse("Si").unwrap(), Target::Element(14));
        assert_eq!(Target::parse("O2").unwrap(), Target::Element(8));
        let Target::Compound { material, formula } = Target::parse("SiO2").unwrap() else {
            panic!("compound expected")
        };
        assert_eq!(formula, "SiO2");
        assert_eq!(material.name(), Some("SiO2"));
        let x = material.atom_fractions();
        assert_eq!(x[0].0, 14);
        assert!((x[0].1 - 1.0 / 3.0).abs() < 1e-15);
        assert!(Target::parse("Si0.5Ge0.5").is_ok());
        for bad in ["", "si", "Xx", "SiSi", "Si0", "Si(O)2", "Si-1", "Si1.2.3"] {
            assert!(Target::parse(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn metadata_keys() {
        assert_eq!(metadata_entry(" name: x: y"), Some(("name", "x: y")));
        assert_eq!(metadata_entry(" Note: free"), None);
        assert_eq!(metadata_entry(" a comment"), None);
    }
}
