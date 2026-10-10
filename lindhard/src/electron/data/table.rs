use super::{
    check_energy_grid, check_in_range, check_nonnegative, checked_text, invalid, lerp_on_grid,
    parse_error, read_text, ElectronDataError, Result, Subshell, CACHE_FORMAT_VERSION,
};
use serde::de::IgnoredAny;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

/// What a [`CrossSectionTable`] samples. Elastic and inelastic tables are
/// distinct types of data and are never interchangeable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SamplingAxis {
    /// Elastic scattering: the polar deflection angle `θ` in radians,
    /// `0 <= θ <= π`.
    ElasticPolarAngle,
    /// Inelastic scattering: the energy loss `W` in eV, `0 <= W <= E` for an
    /// incident energy `E`.
    InelasticEnergyLoss,
}

impl SamplingAxis {
    /// Unit of the sampled variable, as written in the cache file.
    pub fn unit(self) -> &'static str {
        match self {
            SamplingAxis::ElasticPolarAngle => "rad",
            SamplingAxis::InelasticEnergyLoss => "eV",
        }
    }
}

/// The unvalidated contents of a [`CrossSectionTable`], for
/// [`CrossSectionTable::new`]. See [`CrossSectionTable`] for the meaning and
/// invariants of each field.
#[derive(Debug, Clone, PartialEq)]
pub struct CrossSectionTableParts {
    /// Identity of the model that produced the table (name and version).
    pub model: String,
    /// Material identity.
    pub material: String,
    /// Provenance: the inputs of the model and their origin.
    pub provenance: String,
    /// The sampled variable.
    pub axis: SamplingAxis,
    /// Incident electron energies, eV.
    pub energy_ev: Vec<f64>,
    /// Inverse mean free path at each energy, m⁻¹.
    pub inverse_mfp_per_m: Vec<f64>,
    /// Cumulative-probability grid of the inverse CDFs.
    pub probability: Vec<f64>,
    /// Inverse-CDF ordinates, one row per energy (empty for zero-rate rows).
    pub quantiles: Vec<Vec<f64>>,
}

/// A precomputed cross-section table for one material and one process,
/// cacheable on disk (format version [`CACHE_FORMAT_VERSION`]).
///
/// # Contents and dimensions
///
/// With `N` incident energies and `M` probability points:
///
/// - `energy_ev[N]`: incident electron energies, eV; finite, positive,
///   strictly increasing, `N >= 2`.
/// - `inverse_mfp_per_m[N]`: inverse mean free path `λ⁻¹(E) = n σ(E)` for the
///   process, in m⁻¹; finite and non-negative. Zero means no interaction at
///   that energy (for example below a threshold).
/// - `probability[M]`: the cumulative-probability axis shared by every row;
///   finite, strictly increasing, starting at exactly 0 and ending at exactly
///   1, `M >= 2`.
/// - `quantiles[N][..]`: **inverse-CDF ordinates** (the cache stores the
///   inverse CDF, not the forward CDF). `quantiles[i][j]` is the value `x` of
///   the sampled variable at which the cumulative distribution at energy
///   `energy_ev[i]` reaches `probability[j]`. A row is non-decreasing and lies
///   in the axis's domain (`[0, π]` rad for the elastic polar angle, `[0,
///   energy_ev[i]]` eV for the inelastic energy loss). A row whose inverse
///   mean free path is zero must be **empty**: there is no distribution to
///   store, and [`CrossSectionTable::inverse_cdf`] returns
///   [`ElectronDataError::ZeroRate`] for it. Every other row has exactly `M`
///   entries.
///
/// The file also records the model identity, the material, the provenance,
/// the sampling axis and its unit (`ordinate_unit`, which must match the
/// axis), and the format version.
///
/// # Sampling
///
/// [`CrossSectionTable::inverse_cdf`] evaluates a stored row at a
/// probability `u` by linear interpolation between neighbouring probability
/// points (exact at the stored points). Interpolation between incident
/// energies is left to the transport code, which owns that choice.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "RawCacheIn", into = "RawCacheOut")]
pub struct CrossSectionTable {
    parts: CrossSectionTableParts,
}

/// Deserialization form. Every field is optional so that the format version
/// is checked first, whatever else the file holds; unknown keys are collected
/// and rejected after the version check.
#[derive(Debug, Deserialize)]
struct RawCacheIn {
    format_version: Option<u32>,
    model: Option<String>,
    material: Option<String>,
    provenance: Option<String>,
    axis: Option<SamplingAxis>,
    ordinate_unit: Option<String>,
    energy_ev: Option<Vec<f64>>,
    inverse_mfp_per_m: Option<Vec<f64>>,
    probability: Option<Vec<f64>>,
    quantiles: Option<Vec<Vec<f64>>>,
    #[serde(flatten)]
    unknown: BTreeMap<String, IgnoredAny>,
}

#[derive(Debug, Clone, Serialize)]
struct RawCacheOut {
    format_version: u32,
    model: String,
    material: String,
    provenance: String,
    axis: SamplingAxis,
    ordinate_unit: String,
    energy_ev: Vec<f64>,
    inverse_mfp_per_m: Vec<f64>,
    probability: Vec<f64>,
    quantiles: Vec<Vec<f64>>,
}

fn required<T>(field: &'static str, v: Option<T>) -> Result<T> {
    v.ok_or_else(|| ElectronDataError::Parse(format!("missing field `{field}`")))
}

impl TryFrom<RawCacheIn> for CrossSectionTable {
    type Error = ElectronDataError;
    fn try_from(r: RawCacheIn) -> Result<Self> {
        if r.format_version != Some(CACHE_FORMAT_VERSION) {
            return Err(ElectronDataError::UnsupportedVersion {
                found: r.format_version,
                supported: CACHE_FORMAT_VERSION,
            });
        }
        if let Some(k) = r.unknown.keys().next() {
            return Err(ElectronDataError::Parse(format!(
                "unknown field `{k}` in a version {CACHE_FORMAT_VERSION} cache"
            )));
        }
        let provenance = checked_text("provenance", r.provenance.as_deref())?;
        let model = checked_text("model", r.model.as_deref())?;
        let material = checked_text("material", r.material.as_deref())?;
        let axis = required("axis", r.axis)?;
        let unit = checked_text("ordinate_unit", r.ordinate_unit.as_deref())?;
        if unit != axis.unit() {
            return Err(ElectronDataError::MalformedText {
                field: "ordinate_unit",
                reason: format!(
                    "{unit:?} does not match axis {axis:?} (unit {:?})",
                    axis.unit()
                ),
            });
        }
        Self::new(CrossSectionTableParts {
            model,
            material,
            provenance,
            axis,
            energy_ev: required("energy_ev", r.energy_ev)?,
            inverse_mfp_per_m: required("inverse_mfp_per_m", r.inverse_mfp_per_m)?,
            probability: required("probability", r.probability)?,
            quantiles: required("quantiles", r.quantiles)?,
        })
    }
}

impl From<CrossSectionTable> for RawCacheOut {
    fn from(t: CrossSectionTable) -> Self {
        let p = t.parts;
        Self {
            format_version: CACHE_FORMAT_VERSION,
            ordinate_unit: p.axis.unit().to_string(),
            model: p.model,
            material: p.material,
            provenance: p.provenance,
            axis: p.axis,
            energy_ev: p.energy_ev,
            inverse_mfp_per_m: p.inverse_mfp_per_m,
            probability: p.probability,
            quantiles: p.quantiles,
        }
    }
}

impl CrossSectionTable {
    /// Build a validated table. See the type docs for the invariants.
    pub fn new(parts: CrossSectionTableParts) -> Result<Self> {
        let provenance = checked_text("provenance", Some(&parts.provenance))?;
        let model = checked_text("model", Some(&parts.model))?;
        let material = checked_text("material", Some(&parts.material))?;
        let n = parts.energy_ev.len();
        check_energy_grid("incident energy grid", &parts.energy_ev)?;
        if parts.inverse_mfp_per_m.len() != n {
            return Err(invalid(
                "inverse mean free path",
                format!(
                    "need one value per energy: {n} energies, {} values",
                    parts.inverse_mfp_per_m.len()
                ),
            ));
        }
        check_nonnegative("inverse mean free path", &parts.inverse_mfp_per_m)?;
        let p = &parts.probability;
        if p.len() < 2 {
            return Err(invalid("probability axis", "needs at least 2 points"));
        }
        if p.iter().any(|x| !x.is_finite()) {
            return Err(invalid("probability axis", "values must be finite"));
        }
        if p[0] != 0.0 || p[p.len() - 1] != 1.0 {
            return Err(invalid(
                "probability axis",
                format!(
                    "must start at 0 and end at 1, got {} .. {}",
                    p[0],
                    p[p.len() - 1]
                ),
            ));
        }
        if let Some(w) = p.windows(2).find(|w| w[1] <= w[0]) {
            return Err(invalid(
                "probability axis",
                format!("must be strictly increasing, got {} then {}", w[0], w[1]),
            ));
        }
        if parts.quantiles.len() != n {
            return Err(invalid(
                "inverse-CDF rows",
                format!(
                    "need one row per energy: {n} energies, {} rows",
                    parts.quantiles.len()
                ),
            ));
        }
        for (i, row) in parts.quantiles.iter().enumerate() {
            let e = parts.energy_ev[i];
            if parts.inverse_mfp_per_m[i] == 0.0 {
                if !row.is_empty() {
                    return Err(invalid(
                        "inverse-CDF rows",
                        format!(
                            "row {i} ({e} eV) has zero inverse mean free path and must be empty"
                        ),
                    ));
                }
                continue;
            }
            if row.len() != p.len() {
                return Err(invalid(
                    "inverse-CDF rows",
                    format!(
                        "row {i} ({e} eV) has {} entries; the probability axis has {}",
                        row.len(),
                        p.len()
                    ),
                ));
            }
            let max = match parts.axis {
                SamplingAxis::ElasticPolarAngle => std::f64::consts::PI,
                SamplingAxis::InelasticEnergyLoss => e,
            };
            if let Some(x) = row
                .iter()
                .find(|x| !(x.is_finite() && **x >= 0.0 && **x <= max))
            {
                return Err(invalid(
                    "inverse-CDF rows",
                    format!(
                        "row {i} ({e} eV): {x} is outside [0, {max}] {}",
                        parts.axis.unit()
                    ),
                ));
            }
            if let Some(w) = row.windows(2).find(|w| w[1] < w[0]) {
                return Err(invalid(
                    "inverse-CDF rows",
                    format!(
                        "row {i} ({e} eV) must be non-decreasing, got {} then {}",
                        w[0], w[1]
                    ),
                ));
            }
        }
        Ok(Self {
            parts: CrossSectionTableParts {
                model,
                material,
                provenance,
                ..parts
            },
        })
    }

    /// Parse and validate the TOML cache form.
    pub fn from_toml_str(text: &str) -> Result<Self> {
        // Parse the raw form, then convert, so the caller gets the typed
        // error (an unsupported version, a missing provenance) rather than
        // the message toml wraps a conversion error in.
        let raw: RawCacheIn = toml::from_str(text).map_err(parse_error)?;
        Self::try_from(raw)
    }

    /// Serialize to the TOML cache form.
    pub fn to_toml_string(&self) -> Result<String> {
        toml::to_string(self).map_err(|e| ElectronDataError::Parse(e.to_string()))
    }

    /// Load a TOML cache file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_toml_str(&read_text(path.as_ref())?)
    }

    /// Write a TOML cache file.
    pub fn write_toml_file(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        std::fs::write(path, self.to_toml_string()?)
            .map_err(|e| ElectronDataError::Io(format!("{}: {e}", path.display())))
    }

    /// Format version of this table ([`CACHE_FORMAT_VERSION`]).
    pub fn format_version(&self) -> u32 {
        CACHE_FORMAT_VERSION
    }

    /// Model identity.
    pub fn model(&self) -> &str {
        &self.parts.model
    }

    /// Material identity.
    pub fn material(&self) -> &str {
        &self.parts.material
    }

    /// Provenance string.
    pub fn provenance(&self) -> &str {
        &self.parts.provenance
    }

    /// The sampled variable.
    pub fn axis(&self) -> SamplingAxis {
        self.parts.axis
    }

    /// Incident energies, eV.
    pub fn energy_ev(&self) -> &[f64] {
        &self.parts.energy_ev
    }

    /// Inverse mean free paths, m⁻¹, one per energy.
    pub fn inverse_mfp_per_m(&self) -> &[f64] {
        &self.parts.inverse_mfp_per_m
    }

    /// The cumulative-probability axis.
    pub fn probability(&self) -> &[f64] {
        &self.parts.probability
    }

    /// The inverse-CDF ordinates at grid energy `i`, or `None` for a
    /// zero-rate row (or an index past the end).
    pub fn quantiles(&self, i: usize) -> Option<&[f64]> {
        self.parts
            .quantiles
            .get(i)
            .filter(|r| !r.is_empty())
            .map(Vec::as_slice)
    }

    /// The inverse CDF of grid row `i` at cumulative probability `u` in
    /// `[0, 1]`, in the axis's unit. A zero-rate row is
    /// [`ElectronDataError::ZeroRate`]: the caller must not sample an
    /// interaction that cannot happen.
    pub fn inverse_cdf(&self, i: usize, u: f64) -> Result<f64> {
        let n = self.parts.energy_ev.len();
        if i >= n {
            return Err(ElectronDataError::IndexOutOfRange { index: i, len: n });
        }
        check_in_range("cumulative probability", u, 0.0, 1.0)?;
        let Some(row) = self.quantiles(i) else {
            return Err(ElectronDataError::ZeroRate {
                energy_ev: self.parts.energy_ev[i],
            });
        };
        Ok(lerp_on_grid(&self.parts.probability, row, u))
    }

    /// The validated contents.
    pub fn parts(&self) -> &CrossSectionTableParts {
        &self.parts
    }
}

/// The format of a [`ShellChannelTable`] file this build reads and writes.
///
/// It versions only the wrapper (shell identity, binding energy, binding
/// provenance); the wrapped [`CrossSectionTable`] keeps its own
/// [`CACHE_FORMAT_VERSION`], which this type does not change. Bump it whenever
/// the meaning or layout of the wrapper changes.
pub const SHELL_CHANNEL_FORMAT_VERSION: u32 = 1;

/// The energy-loss table of one inner-shell ionisation channel, with the
/// identity of its shell.
///
/// The wrapped [`CrossSectionTable`] holds the channel's inverse mean free
/// path and the inverse CDF of its loss `ω`, built from that shell's own
/// optical ELF (`crate::electron::inelastic::inner_shell`); the wrapper adds
/// the shell (`Z`, subshell), its binding energy `B` and the provenance of
/// `B` (for example [`SubshellBindingTable::provenance`](super::SubshellBindingTable::provenance)). Built by
/// [`crate::electron::inelastic::table::build_shell_channel_tables`] and
/// consumed by [`crate::electron::transport::Transport::with_inner_shells`].
///
/// # Invariants
///
/// - `Z >= 1`; `B` finite and positive; a non-blank one-line binding
///   provenance.
/// - The table samples [`SamplingAxis::InelasticEnergyLoss`].
/// - **Every stored loss is at least `B`**: a shell contributes only above its
///   edge (de Vera et al. 2022, eq. (2), the step `Θ(E - B)`), so the
///   secondary energy `ω - B` is never negative. Since a loss is also at most
///   the row energy, every row with a nonzero rate lies at an energy of at
///   least `B`.
///
/// # File form
///
/// TOML: `format_version` ([`SHELL_CHANNEL_FORMAT_VERSION`]), `z`,
/// `subshell` (the ENDF label, e.g. `"K"`, `"L3"`), `binding_energy_ev`,
/// `binding_provenance`, and a `[table]` holding the [`CrossSectionTable`]
/// cache form (with its own `format_version`). The wrapper version is checked
/// first, unknown keys are rejected, and the same checks as
/// [`ShellChannelTable::new`] run on read.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "RawShellChannelIn", into = "RawShellChannelOut")]
pub struct ShellChannelTable {
    z: u8,
    subshell: Subshell,
    binding_energy_ev: f64,
    binding_provenance: String,
    table: CrossSectionTable,
}

#[derive(Debug, Deserialize)]
struct RawShellChannelIn {
    format_version: Option<u32>,
    z: Option<u8>,
    subshell: Option<String>,
    binding_energy_ev: Option<f64>,
    binding_provenance: Option<String>,
    table: Option<RawCacheIn>,
    #[serde(flatten)]
    unknown: BTreeMap<String, IgnoredAny>,
}

#[derive(Debug, Clone, Serialize)]
struct RawShellChannelOut {
    format_version: u32,
    z: u8,
    subshell: String,
    binding_energy_ev: f64,
    binding_provenance: String,
    table: CrossSectionTable,
}

impl TryFrom<RawShellChannelIn> for ShellChannelTable {
    type Error = ElectronDataError;
    fn try_from(r: RawShellChannelIn) -> Result<Self> {
        if r.format_version != Some(SHELL_CHANNEL_FORMAT_VERSION) {
            return Err(ElectronDataError::UnsupportedVersion {
                found: r.format_version,
                supported: SHELL_CHANNEL_FORMAT_VERSION,
            });
        }
        if let Some(k) = r.unknown.keys().next() {
            return Err(ElectronDataError::Parse(format!(
                "unknown field `{k}` in a version {SHELL_CHANNEL_FORMAT_VERSION} shell-channel table"
            )));
        }
        let label = checked_text("subshell", r.subshell.as_deref())?;
        let subshell =
            Subshell::from_label(&label).ok_or_else(|| ElectronDataError::MalformedText {
                field: "subshell",
                reason: format!("{label:?} is not an ENDF subshell label"),
            })?;
        Self::new(
            required("z", r.z)?,
            subshell,
            required("binding_energy_ev", r.binding_energy_ev)?,
            r.binding_provenance.as_deref().unwrap_or_default(),
            CrossSectionTable::try_from(required("table", r.table)?)?,
        )
    }
}

impl From<ShellChannelTable> for RawShellChannelOut {
    fn from(t: ShellChannelTable) -> Self {
        Self {
            format_version: SHELL_CHANNEL_FORMAT_VERSION,
            z: t.z,
            subshell: t.subshell.label().to_string(),
            binding_energy_ev: t.binding_energy_ev,
            binding_provenance: t.binding_provenance,
            table: t.table,
        }
    }
}

impl ShellChannelTable {
    /// A validated channel table: shell `subshell` of element `z` with
    /// binding energy `binding_energy_ev` (whose origin `binding_provenance`
    /// states) and its loss `table`. See the type docs for the invariants.
    pub fn new(
        z: u8,
        subshell: Subshell,
        binding_energy_ev: f64,
        binding_provenance: &str,
        table: CrossSectionTable,
    ) -> Result<Self> {
        let name = format!("Z = {z} {}", subshell.label());
        if z == 0 {
            return Err(invalid("inner-shell channel", "Z must be at least 1"));
        }
        if !(binding_energy_ev.is_finite() && binding_energy_ev > 0.0) {
            return Err(invalid(
                "inner-shell channel",
                format!(
                    "{name}: binding energy must be finite and positive, got {binding_energy_ev} eV"
                ),
            ));
        }
        let binding_provenance = checked_text("provenance", Some(binding_provenance))?;
        if table.axis() != SamplingAxis::InelasticEnergyLoss {
            return Err(invalid(
                "inner-shell channel",
                format!(
                    "{name}: the table samples {:?}, not the energy loss",
                    table.axis()
                ),
            ));
        }
        for i in 0..table.energy_ev().len() {
            if let Some(q) = table.quantiles(i) {
                if let Some(w) = q.iter().find(|&&w| w < binding_energy_ev) {
                    return Err(invalid(
                        "inner-shell channel",
                        format!(
                            "{name}: row {i} ({} eV) has the loss {w} eV below the binding \
                             energy {binding_energy_ev} eV (a shell contributes only above \
                             its edge)",
                            table.energy_ev()[i]
                        ),
                    ));
                }
            }
        }
        Ok(Self {
            z,
            subshell,
            binding_energy_ev,
            binding_provenance,
            table,
        })
    }

    /// Atomic number of the shell's element.
    pub fn z(&self) -> u8 {
        self.z
    }

    /// The subshell.
    pub fn subshell(&self) -> Subshell {
        self.subshell
    }

    /// Binding energy `B`, eV.
    pub fn binding_energy_ev(&self) -> f64 {
        self.binding_energy_ev
    }

    /// Where `B` came from.
    pub fn binding_provenance(&self) -> &str {
        &self.binding_provenance
    }

    /// The channel's loss table.
    pub fn table(&self) -> &CrossSectionTable {
        &self.table
    }

    /// Format version of the wrapper ([`SHELL_CHANNEL_FORMAT_VERSION`]).
    pub fn format_version(&self) -> u32 {
        SHELL_CHANNEL_FORMAT_VERSION
    }

    /// Parse and validate the TOML file form.
    pub fn from_toml_str(text: &str) -> Result<Self> {
        let raw: RawShellChannelIn = toml::from_str(text).map_err(parse_error)?;
        Self::try_from(raw)
    }

    /// Serialize to the TOML file form.
    pub fn to_toml_string(&self) -> Result<String> {
        toml::to_string(self).map_err(|e| ElectronDataError::Parse(e.to_string()))
    }

    /// Load a TOML file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_toml_str(&read_text(path.as_ref())?)
    }

    /// Write a TOML file.
    pub fn write_toml_file(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        std::fs::write(path, self.to_toml_string()?)
            .map_err(|e| ElectronDataError::Io(format!("{}: {e}", path.display())))
    }
}
