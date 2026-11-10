//! Photon interaction data: the EPDL2017 tables and the reader that makes them.
//!
//! The committed table (`epdl2017.toml`, embedded in the library) holds, for
//! each element `Z = 1..=92`, five curves read from the Evaluated Photon Data
//! Library (EPDL97, as distributed in EPICS2017; D. E. Cullen, IAEA-NDS-225
//! Rev. 1 (February 2018)) file `EPDL2017.ALL`, in ENDF-6 format:
//!
//! | Curve | ENDF-6 section | Axes |
//! |---|---|---|
//! | `photoelectric` | MF=23, MT=522 (total photoionization) | photon energy, eV; barns per atom |
//! | `incoherent` | MF=23, MT=504 | photon energy, eV; barns per atom |
//! | `coherent` | MF=23, MT=502 | photon energy, eV; barns per atom |
//! | `form_factor` | MF=27, MT=502 (atomic form factor `F(x, Z)`) | `x`, 1/Å; dimensionless |
//! | `scattering_function` | MF=27, MT=504 (incoherent scattering function `S(x, Z)`) | `x`, 1/Å; dimensionless |
//!
//! `x = sin(θ/2)/λ` with `λ` the photon wavelength in Å (checked against the
//! hydrogenic form factor `1/(1 + (2π a₀ x)²)²`, which the hydrogen table
//! reproduces to 4 digits).
//!
//! Format read from the ENDF-6 Formats Manual, CSEWG Document ENDF-102,
//! BNL-203218-2018-INRE (2018): a `TAB1` record is a HEAD card `[ZA, AWR, ..]`,
//! a card `[0, 0, 0, 0, NR, NP]`, `NR` interpolation pairs, then `NP` pairs
//! `(x, y)`, three per card. The EPDL2017 file tabulates every curve with one
//! interpolation region of law 2 (linear-linear) on a dense grid; the reader
//! refuses any other record.
//!
//! # What is kept
//!
//! Each curve is restricted to a window (energy `1 keV ..= 1 MeV`, below the
//! 1.022 MeV pair-production threshold, which is not modelled; `x` from
//! `1e-3` to `100` 1/Å, which covers `E ≤ 1 MeV` at every angle) and then
//! thinned: a source point is dropped when log-log interpolation between its
//! kept neighbours reproduces it within [`THIN_TOLERANCE`] (relative). Pairs of
//! points at equal energy (an absorption edge) are always both kept. The
//! window ends are the source's own linear interpolation there. The tolerance
//! bounds the effect on each tabulated source point, so the committed values
//! are the source's values within 0.1 % (tested against a local copy of the
//! file when `LINDHARD_EPDL2017` is set). Nothing is typed by hand.

use serde::{Deserialize, Serialize};
use std::path::Path;

/// Lowest photon energy of the committed tables, eV.
pub const ENERGY_MIN_EV: f64 = 1.0e3;
/// Highest photon energy of the committed tables, eV (below the 1.022 MeV
/// pair-production threshold).
pub const ENERGY_MAX_EV: f64 = 1.0e6;
/// Lowest `x` of the committed form-factor and scattering-function curves, 1/Å.
pub const X_MIN: f64 = 1.0e-3;
/// Highest `x` of the committed curves, 1/Å (`E/hc` at 1 MeV is 80.7).
pub const X_MAX: f64 = 100.0;
/// Relative error allowed when a source point is dropped by the thinning.
pub const THIN_TOLERANCE: f64 = 1.0e-3;
/// Highest element number in the table.
pub const Z_MAX: u8 = 92;

const EPDL2017_TOML: &str = include_str!("epdl2017.toml");

/// Errors from reading or validating photon data.
#[derive(Debug, thiserror::Error)]
pub enum PhotonDataError {
    /// The text is not a valid table.
    #[error("photon data parse error: {0}")]
    Parse(String),
    /// A value is out of range or inconsistent.
    #[error("invalid photon data: {0}")]
    Invalid(String),
    /// A file could not be read.
    #[error("cannot read {path}: {reason}")]
    Io {
        /// The file.
        path: String,
        /// The OS error.
        reason: String,
    },
}

fn invalid(what: impl std::fmt::Display) -> PhotonDataError {
    PhotonDataError::Invalid(what.to_string())
}

/// A positive curve `y(x)` on strictly increasing `x`, except that a pair of
/// equal `x` is allowed (a discontinuity, as at an absorption edge). Evaluated
/// by log-log interpolation.
///
/// Below the first point the curve continues as the power law of its first
/// segment (for a scattering function, `S ∝ x²`; for a flat start it stays
/// constant); above the last point it stays at the last value.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "RawCurve", into = "RawCurve")]
pub struct Curve {
    x: Vec<f64>,
    y: Vec<f64>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawCurve {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl TryFrom<RawCurve> for Curve {
    type Error = PhotonDataError;
    fn try_from(r: RawCurve) -> Result<Self, PhotonDataError> {
        Curve::new(r.x, r.y)
    }
}

impl From<Curve> for RawCurve {
    fn from(c: Curve) -> Self {
        RawCurve { x: c.x, y: c.y }
    }
}

impl Curve {
    /// Validate and build a curve: equal lengths, at least two points, all
    /// values finite, `x > 0` and `y > 0`, `x` non-decreasing with the first
    /// two points distinct and no value of `x` repeated three times.
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Result<Self, PhotonDataError> {
        if x.len() != y.len() || x.len() < 2 {
            return Err(invalid("a curve needs at least two (x, y) pairs"));
        }
        if x.iter().chain(&y).any(|v| !(v.is_finite() && *v > 0.0)) {
            return Err(invalid("curve values must be finite and positive"));
        }
        if x.windows(2).any(|w| w[1] < w[0]) || x[1] == x[0] {
            return Err(invalid(
                "curve x must not decrease and the first two points must differ",
            ));
        }
        if x.windows(3).any(|w| w[0] == w[2]) {
            return Err(invalid("curve x repeated more than twice"));
        }
        Ok(Self { x, y })
    }

    /// The abscissae.
    pub fn x(&self) -> &[f64] {
        &self.x
    }

    /// The values.
    pub fn y(&self) -> &[f64] {
        &self.y
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Always false: a curve has at least two points.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// The value at `x > 0`: log-log interpolation; at an `x` equal to a
    /// repeated abscissa, the upper value.
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.x.len();
        let i = self.x.partition_point(|&v| v <= x);
        let (a, b) = match i {
            0 => (0, 1),
            _ if i == n => return self.y[n - 1],
            _ => (i - 1, i),
        };
        loglog(self.x[a], self.y[a], self.x[b], self.y[b], x)
    }
}

/// Power law through `(xa, ya)` and `(xb, yb)` evaluated at `x`.
fn loglog(xa: f64, ya: f64, xb: f64, yb: f64, x: f64) -> f64 {
    let slope = (yb / ya).ln() / (xb / xa).ln();
    ya * (x / xa).powf(slope)
}

/// The photon data of one element.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtomPhotonData {
    /// Atomic number.
    pub z: u8,
    /// Photoelectric (total photoionization) cross section, barns per atom
    /// against photon energy in eV.
    pub photoelectric: Curve,
    /// Incoherent (Compton) cross section, barns per atom against eV.
    pub incoherent: Curve,
    /// Coherent (Rayleigh) cross section, barns per atom against eV.
    pub coherent: Curve,
    /// Atomic form factor `F(x, Z)` against `x` in 1/Å.
    pub form_factor: Curve,
    /// Incoherent scattering function `S(x, Z)` against `x` in 1/Å.
    pub scattering_function: Curve,
}

impl AtomPhotonData {
    /// Sum of the three cross sections at `energy_ev`, barns per atom.
    pub fn total_barns(&self, energy_ev: f64) -> f64 {
        self.photoelectric.eval(energy_ev)
            + self.incoherent.eval(energy_ev)
            + self.coherent.eval(energy_ev)
    }
}

/// Photon data for a set of elements, with its provenance.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhotonTable {
    provenance: String,
    atoms: Vec<AtomPhotonData>,
}

impl PhotonTable {
    /// Build a table; atoms in strictly increasing `Z`, `1..=92`, and a
    /// non-blank provenance.
    pub fn new(
        provenance: impl Into<String>,
        atoms: Vec<AtomPhotonData>,
    ) -> Result<Self, PhotonDataError> {
        let provenance = provenance.into();
        if provenance.trim().is_empty() {
            return Err(invalid("provenance must not be blank"));
        }
        if atoms.is_empty() {
            return Err(invalid("photon table has no atoms"));
        }
        if atoms.iter().any(|a| a.z == 0 || a.z > Z_MAX)
            || atoms.windows(2).any(|w| w[1].z <= w[0].z)
        {
            return Err(invalid(
                "atoms must have strictly increasing Z within 1..=92",
            ));
        }
        Ok(Self { provenance, atoms })
    }

    /// Provenance string (copied into run metadata).
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    /// All atoms, increasing `Z`.
    pub fn atoms(&self) -> &[AtomPhotonData] {
        &self.atoms
    }

    /// The data of element `z`, if tabulated.
    pub fn atom(&self, z: u8) -> Option<&AtomPhotonData> {
        self.atoms
            .binary_search_by_key(&z, |a| a.z)
            .ok()
            .map(|i| &self.atoms[i])
    }

    /// Parse the TOML form.
    pub fn from_toml_str(text: &str) -> Result<Self, PhotonDataError> {
        let t: PhotonTable =
            toml::from_str(text).map_err(|e| PhotonDataError::Parse(e.to_string()))?;
        Self::new(t.provenance, t.atoms)
    }

    /// Serialize to the TOML form.
    pub fn to_toml_string(&self) -> Result<String, PhotonDataError> {
        toml::to_string(self).map_err(|e| PhotonDataError::Parse(e.to_string()))
    }

    /// The committed EPDL2017 table, `Z = 1..=92`, 1 keV to 1 MeV (see the
    /// module docs; the provenance row is in `docs/data-provenance.md`).
    ///
    /// # Panics
    ///
    /// Never in a correct build: the embedded file is validated by the tests.
    pub fn epdl2017() -> Self {
        Self::from_toml_str(EPDL2017_TOML).expect("the embedded EPDL2017 table is valid")
    }

    /// Read the five curves of every element `Z = 1..=92` from an
    /// `EPDL2017.ALL` file's text (ENDF-6; see the module docs) and thin them.
    /// Materials beyond `Z = 92` are skipped. `provenance` names the file,
    /// edition and retrieval.
    pub fn from_epdl_endf6_str(text: &str, provenance: &str) -> Result<Self, PhotonDataError> {
        let sections = read_sections(text)?;
        let mut atoms = Vec::new();
        for z in 1..=Z_MAX {
            let get = |mf: u32, mt: u32, lo: f64, hi: f64| -> Result<Curve, PhotonDataError> {
                let pts = sections.get(&(z, mf, mt)).ok_or_else(|| {
                    invalid(format!("Z = {z}: no section MF={mf} MT={mt} in the file"))
                })?;
                windowed_thinned(pts, lo, hi)
                    .map_err(|e| invalid(format!("Z = {z}, MF={mf} MT={mt}: {e}")))
            };
            atoms.push(AtomPhotonData {
                z,
                photoelectric: get(23, 522, ENERGY_MIN_EV, ENERGY_MAX_EV)?,
                incoherent: get(23, 504, ENERGY_MIN_EV, ENERGY_MAX_EV)?,
                coherent: get(23, 502, ENERGY_MIN_EV, ENERGY_MAX_EV)?,
                form_factor: get(27, 502, X_MIN, X_MAX)?,
                scattering_function: get(27, 504, X_MIN, X_MAX)?,
            });
        }
        Self::new(provenance, atoms)
    }

    /// [`PhotonTable::from_epdl_endf6_str`] on a file.
    pub fn from_epdl_endf6_file(
        path: impl AsRef<Path>,
        provenance: &str,
    ) -> Result<Self, PhotonDataError> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| PhotonDataError::Io {
            path: path.display().to_string(),
            reason: e.to_string(),
        })?;
        Self::from_epdl_endf6_str(&text, provenance)
    }
}

type SectionKey = (u8, u32, u32);
type Sections = std::collections::BTreeMap<SectionKey, Vec<(f64, f64)>>;
/// The cards of one section being collected: (line number, first 66 columns).
type OpenSection<'a> = Option<(SectionKey, Vec<(usize, &'a str)>)>;

fn parse_float(raw: &str) -> Option<f64> {
    let s = raw.trim().replace(['D', 'd'], "E");
    if s.is_empty() {
        return Some(0.0);
    }
    if let Ok(v) = s.parse::<f64>() {
        return Some(v);
    }
    // ENDF omits the `E` of the exponent: `1.234567+2`.
    let p = s.get(1..)?.rfind(['+', '-'])? + 1;
    format!("{}E{}", &s[..p], &s[p..]).parse::<f64>().ok()
}

/// Every `TAB1` section of MF 23 (MT 502, 504, 522) and MF 27 (MT 502, 504),
/// keyed by `(Z, MF, MT)`, for `Z <= 92`.
fn read_sections(text: &str) -> Result<Sections, PhotonDataError> {
    const WANTED: [(u32, u32); 5] = [(23, 502), (23, 504), (23, 522), (27, 502), (27, 504)];
    let mut out = Sections::new();
    let mut cur: OpenSection<'_> = None;
    let mut flush = |cur: &mut OpenSection<'_>| match cur.take() {
        Some((key, cards)) => parse_tab1(&cards).map(|pts| {
            out.insert(key, pts);
        }),
        None => Ok(()),
    };
    for (i, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.len() < 75 || !line.is_ascii() {
            continue;
        }
        let num = |a: usize, b: usize| line[a..b].trim().parse::<u32>().ok();
        let (Some(mat), Some(mf), Some(mt)) = (num(66, 70), num(70, 72), num(72, 75)) else {
            continue;
        };
        let key = (mat / 100, mf, mt);
        let wanted = WANTED.contains(&(mf, mt)) && (1..=u32::from(Z_MAX)).contains(&key.0);
        if !wanted {
            flush(&mut cur)?;
            continue;
        }
        let key = (key.0 as u8, mf, mt);
        match &mut cur {
            Some((k, cards)) if *k == key => cards.push((i + 1, &line[..66])),
            _ => {
                flush(&mut cur)?;
                cur = Some((key, vec![(i + 1, &line[..66])]));
            }
        }
    }
    flush(&mut cur)?;
    Ok(out)
}

fn parse_tab1(cards: &[(usize, &str)]) -> Result<Vec<(f64, f64)>, PhotonDataError> {
    let perr =
        |line: usize, what: &str| PhotonDataError::Parse(format!("ENDF line {line}: {what}"));
    if cards.len() < 4 {
        return Err(perr(cards[0].0, "TAB1 record too short"));
    }
    let field = |k: usize, f: usize| -> Result<f64, PhotonDataError> {
        let (ln, c) = cards[k];
        let raw = c.get(11 * f..11 * (f + 1)).unwrap_or("");
        parse_float(raw).ok_or_else(|| perr(ln, &format!("field {} is not a number", f + 1)))
    };
    let (nr, np) = (field(1, 4)?, field(1, 5)?);
    if nr != 1.0 {
        return Err(perr(cards[1].0, "more than one interpolation region"));
    }
    if field(2, 1)? != 2.0 || field(2, 0)? != np {
        return Err(perr(
            cards[2].0,
            "interpolation law is not 2 (linear-linear) over the whole table",
        ));
    }
    let np = np as usize;
    if cards.len() != 3 + np.div_ceil(3) {
        return Err(perr(
            cards[0].0,
            "TAB1 point count does not match its cards",
        ));
    }
    let mut pts = Vec::with_capacity(np);
    for k in 0..np {
        let (card, f) = (3 + k / 3, 2 * (k % 3));
        pts.push((field(card, f)?, field(card, f + 1)?));
    }
    Ok(pts)
}

/// Source linear-linear value at `x` between points (the law of the file).
fn lin_at(pts: &[(f64, f64)], x: f64) -> f64 {
    let i = pts.partition_point(|p| p.0 <= x).clamp(1, pts.len() - 1);
    let ((xa, ya), (xb, yb)) = (pts[i - 1], pts[i]);
    ya + (yb - ya) * (x - xa) / (xb - xa)
}

/// Clip to `[lo, hi]` (ends from the source's own interpolation) and thin.
fn windowed_thinned(pts: &[(f64, f64)], lo: f64, hi: f64) -> Result<Curve, PhotonDataError> {
    if pts.first().is_none_or(|p| p.0 > lo) || pts.last().is_none_or(|p| p.0 < hi) {
        return Err(invalid("source does not cover the window"));
    }
    let mut w = vec![(lo, lin_at(pts, lo))];
    w.extend(pts.iter().copied().filter(|p| p.0 > lo && p.0 < hi));
    w.push((hi, lin_at(pts, hi)));
    let kept = thin(&w, THIN_TOLERANCE);
    Curve::new(
        kept.iter().map(|p| p.0).collect(),
        kept.iter().map(|p| p.1).collect(),
    )
}

/// Greedy thinning: keep a point only when log-log interpolation between the
/// kept neighbours misses it by more than `tol` (relative). Equal-`x` pairs
/// and the ends are always kept.
fn thin(w: &[(f64, f64)], tol: f64) -> Vec<(f64, f64)> {
    let n = w.len();
    let mut must = vec![false; n];
    must[0] = true;
    must[n - 1] = true;
    for k in 0..n - 1 {
        if w[k].0 == w[k + 1].0 {
            must[k] = true;
            must[k + 1] = true;
        }
    }
    let mut out = vec![w[0]];
    let mut a = 0;
    while a < n - 1 {
        let mut b = a + 1;
        while b + 1 < n && !must[b] && segment_ok(w, a, b + 1, tol) {
            b += 1;
        }
        out.push(w[b]);
        a = b;
    }
    out
}

/// Whether the power law through `w[a]` and `w[b]` reproduces `w[a+1..b]`.
fn segment_ok(w: &[(f64, f64)], a: usize, b: usize, tol: f64) -> bool {
    if w[b].0 == w[a].0 {
        return false;
    }
    (a + 1..b).all(|k| {
        let v = loglog(w[a].0, w[a].1, w[b].0, w[b].1, w[k].0);
        (v - w[k].1).abs() <= tol * w[k].1
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curve_interpolates_log_log_and_keeps_edges() {
        // y = x^2 exactly, then a jump at x = 4 (equal-x pair).
        let c = Curve::new(
            vec![1.0, 2.0, 4.0, 4.0, 8.0],
            vec![1.0, 4.0, 16.0, 32.0, 128.0],
        )
        .unwrap();
        assert!((c.eval(3.0) - 9.0).abs() < 1e-12);
        assert_eq!(c.eval(4.0), 32.0); // upper value at the edge
        assert!((c.eval(3.999_999) / 16.0 - 1.0).abs() < 1e-5); // lower side of the edge
        assert!((c.eval(0.5) - 0.25).abs() < 1e-12); // first-segment power law
        assert_eq!(c.eval(100.0), 128.0); // clamped above
    }

    #[test]
    fn curve_validation() {
        assert!(Curve::new(vec![1.0], vec![1.0]).is_err());
        assert!(Curve::new(vec![1.0, 2.0], vec![1.0, 0.0]).is_err());
        assert!(Curve::new(vec![2.0, 1.0], vec![1.0, 1.0]).is_err());
        assert!(Curve::new(vec![1.0, 1.0], vec![1.0, 1.0]).is_err());
        assert!(Curve::new(vec![1.0, 2.0, 2.0, 2.0], vec![1.0; 4]).is_err());
    }

    #[test]
    fn thinning_keeps_edges_and_bounds_the_error() {
        let mut w = Vec::new();
        for i in 0..400 {
            let e = 1.0e3 * 10f64.powf(3.0 * f64::from(i) / 399.0);
            // A K-like edge at 1e4 with a different power law on each side.
            let y = if e < 1.0e4 {
                e.powf(-2.5)
            } else {
                20.0 * e.powf(-2.8)
            };
            w.push((e, y));
            if i == 133 {
                w.push((e, 7.0 * y)); // equal-x pair
            }
        }
        let kept = thin(&w, 1e-3);
        assert!(kept.len() < w.len() / 10, "kept {}", kept.len());
        assert!(kept.windows(2).any(|p| p[0].0 == p[1].0), "edge pair lost");
        let c = Curve::new(
            kept.iter().map(|p| p.0).collect(),
            kept.iter().map(|p| p.1).collect(),
        )
        .unwrap();
        for p in w.windows(2).filter(|p| p[0].0 != p[1].0).map(|p| p[0]) {
            assert!((c.eval(p.0) / p.1 - 1.0).abs() < 2e-3, "at {}", p.0);
        }
    }

    #[test]
    fn tab1_parser_reads_a_record() {
        let card = |s: &str, ln: u32| format!("{s:<66} 10023502{ln:>5}\n");
        let mut text = String::new();
        text += &card(
            " 1000.00000 .999242000          0          0          0          0",
            1,
        );
        text += &card(
            " 0.0        0.0                 0          0          1          2",
            2,
        );
        text += &card("          2          2", 3);
        text += &card(" 1.00000000 4.52522E-6 2.00000000 7.39734-5", 4);
        let s = read_sections(&text).unwrap();
        assert_eq!(s[&(1, 23, 502)], vec![(1.0, 4.52522e-6), (2.0, 7.39734e-5)]);
    }

    #[test]
    fn tab1_parser_refuses_other_interpolation() {
        let card = |s: &str, ln: u32| format!("{s:<66} 10023502{ln:>5}\n");
        let mut text = String::new();
        text += &card(
            " 1000.00000 .999242000          0          0          0          0",
            1,
        );
        text += &card(
            " 0.0        0.0                 0          0          1          2",
            2,
        );
        text += &card("          2          5", 3);
        text += &card(" 1.00000000 4.52522E-6 2.00000000 7.39734-5", 4);
        assert!(read_sections(&text).is_err());
    }
}
