//! The `lindhard stopping-table` command: electronic and nuclear stopping and
//! the CSDA range of ions in targets over an energy grid. Models and
//! conventions: [`lindhard::ion::stopping`], [`lindhard::ion::range`]; input
//! and files: `docs/cli.md`.
//!
//! No transport is run. The ion, target, screening, stopping model and user
//! tables are resolved by the same code as `lindhard run` (each (ion, target)
//! pair is resolved as a one-layer [`Input`]), so a model or a compound rule
//! has one implementation.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use clap::Args;
use lindhard::elements::element;
use lindhard::input::{
    resolve_material, stopping_model, Input, LengthChoice, MaterialRef, PotentialChoice, Resolved,
    StoppingChoice, StoppingSpec,
};
use lindhard::ion::potential::Potential;
use lindhard::ion::range::{csda_ranges, nuclear_stopping_cross_section};
use lindhard::ion::stopping::bragg::{bragg_cross_section_per_atom, NoCorrection};
use lindhard::ion::stopping::Ion;
use lindhard::material::{Material, MaterialSpec};
use lindhard::units::J_PER_EV;
use serde::Deserialize;
use serde_json::{json, Value};

use crate::publish::OutputSet;

/// The table of every (ion, target) pair, one CSV.
pub const CSV_FILE: &str = "stopping_table.csv";
/// The same table with the inputs, models and provenance echoed. Published
/// last, so it certifies the completed output set.
pub const JSON_FILE: &str = "stopping_table.json";

/// Arguments of `lindhard stopping-table`.
#[derive(Args, Debug)]
pub struct StoppingTableArgs {
    /// Stopping-table input TOML file (`docs/cli.md`, "Stopping tables").
    pub input: PathBuf,
    /// Output directory (created if missing; files in it are overwritten).
    #[arg(long, short)]
    pub out: PathBuf,
}

/// One projectile of the batch.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IonSpec {
    /// Element symbol.
    pub ion: String,
    /// Mass, u. Default: the standard atomic weight.
    pub mass_amu: Option<f64>,
}

/// `[physics]` of a stopping-table input: the model choices of `run`'s
/// `[physics]` that a stopping table depends on, with the same defaults.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhysicsOpts {
    /// Screening function (nuclear stopping).
    #[serde(default)]
    pub potential: PotentialChoice,
    /// Screening length. Default: the one paired with the screening function.
    pub screening_length: Option<LengthChoice>,
    /// Electronic stopping model.
    #[serde(default)]
    pub stopping: StoppingChoice,
}

/// Spacing of a generated energy grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Scale {
    /// Constant ratio between points.
    #[default]
    Log,
    /// Constant difference between points.
    Linear,
}

/// `[energy]`: an explicit list, or a generated grid.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnergySpec {
    /// Explicit energies, eV, strictly ascending.
    pub values_ev: Option<Vec<f64>>,
    /// First energy of a generated grid, eV.
    pub start_ev: Option<f64>,
    /// Last energy of a generated grid, eV (omit for one point).
    pub stop_ev: Option<f64>,
    /// Number of points of a generated grid.
    pub points: Option<usize>,
    /// Spacing of a generated grid. Default `"log"`.
    pub scale: Option<Scale>,
}

fn finite_pos(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

impl EnergySpec {
    /// The energies, eV, strictly ascending.
    pub fn grid(&self) -> Result<Vec<f64>> {
        if let Some(v) = &self.values_ev {
            if self.start_ev.is_some()
                || self.stop_ev.is_some()
                || self.points.is_some()
                || self.scale.is_some()
            {
                bail!("energy: give either values_ev or start_ev/stop_ev/points/scale, not both");
            }
            if v.is_empty() {
                bail!("energy.values_ev: must not be empty");
            }
            for (i, &e) in v.iter().enumerate() {
                if !finite_pos(e) {
                    bail!("energy.values_ev[{i}]: {e} eV must be finite and positive");
                }
                if i > 0 && e <= v[i - 1] {
                    bail!("energy.values_ev[{i}]: {e} eV is not above the previous energy");
                }
            }
            return Ok(v.clone());
        }
        let start = self
            .start_ev
            .ok_or_else(|| anyhow!("energy: needs values_ev, or start_ev and points"))?;
        let n = self
            .points
            .ok_or_else(|| anyhow!("energy.points: required with start_ev"))?;
        if !finite_pos(start) {
            bail!("energy.start_ev: {start} eV must be finite and positive");
        }
        if n == 0 {
            bail!("energy.points: must be at least 1");
        }
        if n == 1 {
            if self.stop_ev.is_some_and(|s| s != start) {
                bail!("energy.stop_ev: with points = 1 omit it or make it equal to start_ev");
            }
            return Ok(vec![start]);
        }
        let stop = self
            .stop_ev
            .ok_or_else(|| anyhow!("energy.stop_ev: required with points > 1"))?;
        if !(stop.is_finite() && stop > start) {
            bail!("energy.stop_ev: {stop} eV must be finite and above start_ev ({start} eV)");
        }
        let last = (n - 1) as f64;
        let mut v: Vec<f64> = (0..n)
            .map(|i| match self.scale.unwrap_or_default() {
                Scale::Log => start * (stop / start).powf(i as f64 / last),
                Scale::Linear => start + (stop - start) * i as f64 / last,
            })
            .collect();
        v[n - 1] = stop;
        Ok(v)
    }
}

/// `[range]`: the CSDA range column.
#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct RangeSpec {
    /// Write the range column. Default `true`.
    pub enabled: bool,
    /// Energy the range is integrated from, eV: the range is the path length
    /// to slow from the row's energy down to this one, and the path below it
    /// is not counted. Fixed (not the grid's first energy), so a row does not
    /// depend on the other rows. Default 1 eV.
    pub from_ev: f64,
}

impl Default for RangeSpec {
    fn default() -> Self {
        Self {
            enabled: true,
            from_ev: 1.0,
        }
    }
}

/// Unit of the stopping columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
pub enum StoppingUnit {
    /// Linear stopping power, eV/nm.
    #[default]
    #[serde(rename = "eV/nm")]
    EvPerNm,
    /// Linear stopping power, keV/um.
    #[serde(rename = "keV/um")]
    KevPerUm,
    /// Mass stopping power, MeV/(mg/cm^2).
    #[serde(rename = "MeV/(mg/cm^2)")]
    MevPerMgCm2,
    /// Stopping cross section per atom, eV/(1e15 atoms/cm^2).
    #[serde(rename = "eV/1e15 atoms/cm^2")]
    EvPer1e15AtomsCm2,
}

impl StoppingUnit {
    /// The name as written in the input.
    pub fn label(self) -> &'static str {
        match self {
            Self::EvPerNm => "eV/nm",
            Self::KevPerUm => "keV/um",
            Self::MevPerMgCm2 => "MeV/(mg/cm^2)",
            Self::EvPer1e15AtomsCm2 => "eV/1e15 atoms/cm^2",
        }
    }

    /// Factor from a stopping cross section per atom in J m² to this unit,
    /// for atom density `n` (m⁻³) and mass density `rho` (kg/m³).
    fn factor(self, n: f64, rho: f64) -> f64 {
        match self {
            // N S in J/m, to eV/m, to eV/nm.
            Self::EvPerNm => n / J_PER_EV * 1e-9,
            // eV/m to keV/m to keV/um.
            Self::KevPerUm => n / J_PER_EV * 1e-3 * 1e-6,
            // N S / rho in J m^2/kg; 1 mg/cm^2 = 1e-2 kg/m^2.
            Self::MevPerMgCm2 => n / rho / J_PER_EV * 1e-2 * 1e-6,
            // 1e-15 cm^2 = 1e-19 m^2.
            Self::EvPer1e15AtomsCm2 => 1.0 / J_PER_EV / 1e-19,
        }
    }

    /// Convert a stopping cross section per atom, J m², to this unit.
    pub fn from_si(self, s_j_m2: f64, n: f64, rho: f64) -> f64 {
        s_j_m2 * self.factor(n, rho)
    }

    /// Inverse of [`from_si`](Self::from_si).
    pub fn to_si(self, value: f64, n: f64, rho: f64) -> f64 {
        value / self.factor(n, rho)
    }
}

/// Unit of the range column.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
pub enum RangeUnit {
    /// Nanometres.
    #[default]
    #[serde(rename = "nm")]
    Nm,
    /// Micrometres.
    #[serde(rename = "um")]
    Um,
    /// Areal mass density, mg/cm^2.
    #[serde(rename = "mg/cm^2")]
    MgCm2,
}

impl RangeUnit {
    /// The name as written in the input.
    pub fn label(self) -> &'static str {
        match self {
            Self::Nm => "nm",
            Self::Um => "um",
            Self::MgCm2 => "mg/cm^2",
        }
    }

    fn factor(self, rho: f64) -> f64 {
        match self {
            Self::Nm => 1e9,
            Self::Um => 1e6,
            // kg/m^2 over 1e-2 kg/m^2 per mg/cm^2.
            Self::MgCm2 => rho * 1e2,
        }
    }

    /// Convert a length in metres to this unit, for mass density `rho` kg/m³.
    pub fn from_si(self, metres: f64, rho: f64) -> f64 {
        metres * self.factor(rho)
    }

    /// Inverse of [`from_si`](Self::from_si).
    pub fn to_si(self, value: f64, rho: f64) -> f64 {
        value / self.factor(rho)
    }
}

/// `[units]`: the units of the output columns.
#[derive(Debug, Clone, Copy, Default, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct UnitsSpec {
    /// Stopping columns. Default `"eV/nm"`.
    pub stopping: StoppingUnit,
    /// Range column. Default `"nm"`.
    pub range: RangeUnit,
}

/// A stopping-table input.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TableInput {
    /// Projectiles: each is combined with every target.
    pub ions: Vec<IonSpec>,
    /// Targets: a `[materials]` key, an element symbol, or an inline material.
    pub targets: Vec<MaterialRef>,
    /// Named materials (the same form as in a run input).
    #[serde(default)]
    pub materials: BTreeMap<String, MaterialSpec>,
    /// Model choices.
    #[serde(default)]
    pub physics: PhysicsOpts,
    /// User electronic stopping tables.
    #[serde(default)]
    pub stopping: StoppingSpec,
    /// The energies.
    pub energy: EnergySpec,
    /// The range column.
    #[serde(default)]
    pub range: RangeSpec,
    /// Output units.
    #[serde(default)]
    pub units: UnitsSpec,
}

impl TableInput {
    /// Parse a TOML document (unknown keys are errors).
    pub fn from_toml_str(text: &str) -> Result<Self> {
        toml::from_str(text).map_err(|e| anyhow!("{e}"))
    }
}

/// One energy of a table, SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Row {
    /// Ion energy, eV.
    pub energy_ev: f64,
    /// Electronic stopping cross section per atom, J m².
    pub electronic: f64,
    /// Nuclear stopping cross section per atom, J m².
    pub nuclear: f64,
    /// CSDA range, m (`None` with `range.enabled = false`).
    pub range_m: Option<f64>,
}

/// The table of one (ion, target) pair.
#[derive(Debug, Clone)]
pub struct Table {
    /// The pair as resolved (for the echoed models and tables).
    pub resolved: Resolved,
    /// Name of the target in the output.
    pub target: String,
    /// One row per energy.
    pub rows: Vec<Row>,
    /// Non-fatal advice.
    pub warnings: Vec<String>,
}

impl Table {
    fn ion(&self) -> Ion {
        self.resolved.beam.ion
    }

    fn material(&self) -> &Material {
        &self.resolved.layers[0].material
    }
}

/// Resolve one (ion, target) pair through the run resolver.
fn resolve_pair(
    input: &TableInput,
    (i, ion): (usize, &IonSpec),
    (j, target): (usize, &MaterialRef),
    energies: &[f64],
    base_dir: &Path,
) -> Result<Resolved> {
    let e_min = energies[0];
    let e_max = energies[energies.len() - 1];
    let mut physics = json!({
        "potential": input.physics.potential,
        "stopping": input.physics.stopping,
        // Unused by a stopping table, but required and validated by the
        // run resolver; both are below every grid energy.
        "primary_cutoff_ev": 0.5 * e_min,
        "recoil_cutoff_ev": 0.5 * e_min,
        "follow_recoils": false,
    });
    // Displacement and binding energies are required by the run resolver and
    // unused here: a placeholder for every element, so the target needs none.
    let layer = resolve_material(&input.materials, &format!("targets[{j}]"), target)
        .map_err(|e| anyhow!("{e}"))?;
    physics["energies"] = layer
        .material
        .components()
        .iter()
        .filter_map(|c| element(c.z()))
        .map(|e| {
            (
                e.symbol.to_string(),
                json!({ "e_d_ev": 1.0, "e_b_ev": 1.0, "e_s_ev": 1.0 }),
            )
        })
        .collect::<serde_json::Map<String, Value>>()
        .into();
    if let Some(l) = input.physics.screening_length {
        physics["screening_length"] = serde_json::to_value(l)?;
    }
    let mut beam = json!({ "ion": ion.ion, "energy_ev": e_max });
    if let Some(m) = ion.mass_amu {
        beam["mass_amu"] = json!(m);
    }
    let doc = json!({
        "beam": beam,
        "materials": input.materials,
        "target": { "substrate": target },
        "physics": physics,
        "stopping": { "tables": input.stopping.tables },
        "run": { "ions": 1, "seed": 0 },
    });
    let as_input: Input = serde_json::from_value(doc).map_err(|e| anyhow!("{e}"))?;
    as_input.resolve_in(base_dir).map_err(|e| {
        // Name the keys of this input, not the synthetic run input.
        anyhow!(
            "{}",
            e.to_string()
                .replace("target.substrate", &format!("targets[{j}]"))
                .replace("beam.mass_amu", &format!("ions[{i}].mass_amu"))
                .replace("beam.ion", &format!("ions[{i}].ion"))
        )
    })
}

/// Name of a target in the output.
fn target_name(r: &Resolved, j: usize) -> String {
    let l = &r.layers[0];
    if l.source != "inline" {
        l.source.clone()
    } else {
        l.material
            .name()
            .map_or_else(|| format!("targets[{j}]"), str::to_string)
    }
}

/// Where energy `i` came from, for error messages.
fn energy_key(spec: &EnergySpec, i: usize, e: f64) -> String {
    if spec.values_ev.is_some() {
        format!("energy.values_ev[{i}] ({e} eV)")
    } else {
        format!("energy point {i} ({e} eV)")
    }
}

/// The stopping table of every (ion, target) pair of `input`, in input order
/// (ions outer, targets inner). Relative `[stopping]` table paths resolve
/// against `base_dir`.
pub fn compute(input: &TableInput, base_dir: &Path) -> Result<Vec<Table>> {
    if input.ions.is_empty() {
        bail!("ions: needs at least one ion");
    }
    if input.targets.is_empty() {
        bail!("targets: needs at least one target");
    }
    let energies = input.energy.grid()?;
    let range = input.range;
    if range.enabled {
        if !finite_pos(range.from_ev) {
            bail!(
                "range.from_ev: {} eV must be finite and positive",
                range.from_ev
            );
        }
        if energies[0] < range.from_ev {
            bail!(
                "range.from_ev: {} eV is above the lowest energy ({} eV); lower it, raise \
                 the energies, or set range.enabled = false",
                range.from_ev,
                energies[0]
            );
        }
    }
    let mut tables = Vec::new();
    for (i, ion) in input.ions.iter().enumerate() {
        for (j, target) in input.targets.iter().enumerate() {
            let resolved = resolve_pair(input, (i, ion), (j, target), &energies, base_dir)?;
            let name = target_name(&resolved, j);
            let t = pair_table(input, &energies, resolved, name)
                .with_context(|| format!("ions[{i}] ({}) in targets[{j}]", ion.ion))?;
            tables.push(t);
        }
    }
    Ok(tables)
}

fn pair_table(
    input: &TableInput,
    energies: &[f64],
    resolved: Resolved,
    target: String,
) -> Result<Table> {
    let ion = resolved.beam.ion;
    let material = resolved.layers[0].material.clone();
    let model = resolved.stopping_model();
    let (screening, length) = (resolved.screening, resolved.screening_length);
    let z1 = f64::from(ion.z());
    let electronic = |e: f64| -> Result<f64, lindhard::ion::stopping::StoppingError> {
        bragg_cross_section_per_atom(&*model, &NoCorrection, &ion, &material, e)
    };
    // Bragg additivity for the nuclear term too: the atom-fraction-weighted
    // sum over the elements of the target.
    let nuclear = |e: f64| -> Result<f64> {
        let mut sum = 0.0;
        for c in material.components() {
            let m2 = element(c.z()).ok_or_else(|| anyhow!("unknown element Z = {}", c.z()))?;
            let pot = Potential::new(screening, z1, f64::from(c.z())).with_length(length);
            let s = nuclear_stopping_cross_section(&pot, ion.mass_amu(), m2.atomic_weight, e)
                .ok_or_else(|| anyhow!("nuclear stopping undefined at {e} eV"))?;
            sum += c.atom_fraction() * s;
        }
        Ok(sum)
    };
    let mut rows = Vec::with_capacity(energies.len());
    for (k, &e) in energies.iter().enumerate() {
        let se =
            electronic(e).map_err(|err| anyhow!("{}: {err}", energy_key(&input.energy, k, e)))?;
        rows.push(Row {
            energy_ev: e,
            electronic: se,
            nuclear: nuclear(e)?,
            range_m: None,
        });
    }
    if input.range.enabled {
        let from = input.range.from_ev;
        let total = |e: f64| -> Result<f64> {
            let s = electronic(e).map_err(|err| {
                anyhow!(
                    "range: {err} (the range integral reaches {e} eV, between range.from_ev = \
                     {from} eV and the energy; raise range.from_ev or set range.enabled = false)"
                )
            })? + nuclear(e)?;
            if finite_pos(s) {
                Ok(s)
            } else {
                bail!("range: total stopping {s} J m^2 at {e} eV is not finite and positive")
            }
        };
        let r = csda_ranges(total, material.atom_number_density(), from, energies)?;
        for (row, r) in rows.iter_mut().zip(r) {
            row.range_m = Some(r);
        }
    }

    let mut warnings: Vec<String> = resolved
        .warnings
        .iter()
        // The resolver judged the placeholder beam energy; the grid is
        // checked below.
        .filter(|w| !w.starts_with("beam energy "))
        .cloned()
        .collect();
    let uncovered = material.components().iter().any(|c| {
        !resolved
            .stopping_tables
            .iter()
            .any(|t| t.table.covers_pair(ion.z(), c.z()))
    });
    let fallback = stopping_model(input.physics.stopping);
    let v = fallback.validity(&ion);
    if uncovered && energies.iter().any(|&e| !v.contains(e)) {
        warnings.push(format!(
            "some energies are outside the advisory validity range of the {} stopping model \
             for {} ([{}, {}] eV)",
            fallback.name(),
            element(ion.z()).map_or("?", |e| e.symbol),
            v.min_energy_ev,
            v.max_energy_ev
        ));
    }
    Ok(Table {
        resolved,
        target,
        rows,
        warnings,
    })
}

fn csv_field(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Mass density and atom density of a table's target.
fn densities(t: &Table) -> (f64, f64) {
    (
        t.material().mass_density(),
        t.material().atom_number_density(),
    )
}

/// The CSV of every table: one row per (ion, target, energy). Columns with
/// the stopping and range units in their names; the range column is empty
/// when disabled.
pub fn table_csv(tables: &[Table], units: &UnitsSpec) -> String {
    let (su, ru) = (units.stopping.label(), units.range.label());
    let mut s = format!(
        "ion,target,energy_ev,electronic_stopping [{su}],nuclear_stopping [{su}],\
         total_stopping [{su}],csda_range [{ru}]\n"
    );
    for t in tables {
        let (rho, n) = densities(t);
        let symbol = element(t.ion().z()).map_or("?", |e| e.symbol);
        for r in &t.rows {
            let st = |x: f64| units.stopping.from_si(x, n, rho);
            let range = r.range_m.map_or(String::new(), |m| {
                format!("{:e}", units.range.from_si(m, rho))
            });
            s.push_str(&format!(
                "{symbol},{},{:e},{:e},{:e},{:e},{range}\n",
                csv_field(&t.target),
                r.energy_ev,
                st(r.electronic),
                st(r.nuclear),
                st(r.electronic + r.nuclear),
            ));
        }
    }
    s
}

/// The JSON: units, grid and range settings, and per table the ion, target,
/// models, user-table provenance, warnings and rows.
pub fn table_json(input: &TableInput, tables: &[Table]) -> Result<String> {
    let units = &input.units;
    let out: Vec<Value> = tables
        .iter()
        .map(|t| {
            let (rho, n) = densities(t);
            let ion = t.ion();
            let symbol = element(ion.z()).map_or("?", |e| e.symbol);
            let st = |x: f64| units.stopping.from_si(x, n, rho);
            let mut models: Vec<Value> = t
                .resolved
                .models()
                .iter()
                .filter(|m| {
                    m.role.starts_with("screening")
                        || m.role.starts_with("electronic")
                        || m.role == "compound stopping"
                })
                .map(|m| json!({ "role": m.role, "name": m.name, "citation": m.citation }))
                .collect();
            models.push(json!({
                "role": "nuclear stopping",
                "name": "reduced-sn-quadrature",
                "citation": "J. Lindhard, M. Scharff, H. E. Schiott, Mat. Fys. Medd. Dan. Vid. \
                             Selsk. 33 (14) (1963); J. F. Ziegler, J. P. Biersack, U. Littmark, \
                             The Stopping and Range of Ions in Solids (Pergamon, 1985), ch. 2; \
                             summed over the elements of a compound by atom fraction",
            }));
            if input.range.enabled {
                models.push(json!({
                    "role": "range",
                    "name": "csda",
                    "citation": "R(E) = Int dE'/(N (S_e + S_n)) from range.from_ev; mean path \
                                 length of continuous slowing down, not the projected range",
                }));
            }
            let composition: Vec<Value> = t
                .material()
                .components()
                .iter()
                .map(|c| {
                    json!({
                        "symbol": element(c.z()).map_or("?", |e| e.symbol),
                        "z": c.z(),
                        "atom_fraction": c.atom_fraction(),
                    })
                })
                .collect();
            let user_tables: Vec<Value> = t
                .resolved
                .stopping_tables
                .iter()
                .map(|l| {
                    json!({
                        "path": l.path,
                        "sha256": l.sha256,
                        "provenance": l.table.provenance(),
                    })
                })
                .collect();
            let rows: Vec<Value> = t
                .rows
                .iter()
                .map(|r| {
                    json!({
                        "energy_ev": r.energy_ev,
                        "electronic_stopping": st(r.electronic),
                        "nuclear_stopping": st(r.nuclear),
                        "total_stopping": st(r.electronic + r.nuclear),
                        "csda_range": r.range_m.map(|m| units.range.from_si(m, rho)),
                    })
                })
                .collect();
            json!({
                "ion": { "symbol": symbol, "z": ion.z(), "mass_amu": ion.mass_amu() },
                "target": {
                    "name": t.target,
                    "mass_density_kg_m3": rho,
                    "atom_number_density_cm3": n * 1e-6,
                    "composition": composition,
                },
                "models": models,
                "stopping_tables": user_tables,
                "warnings": t.warnings,
                "rows": rows,
            })
        })
        .collect();
    let doc = json!({
        "lindhard_version": env!("CARGO_PKG_VERSION"),
        "units": {
            "energy": "eV",
            "stopping": units.stopping.label(),
            "range": units.range.label(),
        },
        "range": { "enabled": input.range.enabled, "from_ev": input.range.from_ev },
        "tables": out,
    });
    serde_json::to_string_pretty(&doc)
        .map(|mut s| {
            s.push('\n');
            s
        })
        .map_err(Into::into)
}

/// Run `lindhard stopping-table`. Nothing is written unless every table was
/// computed.
pub fn run(a: &StoppingTableArgs) -> Result<()> {
    let text = std::fs::read_to_string(&a.input)
        .with_context(|| format!("reading {}", a.input.display()))?;
    let base = a
        .input
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let input = TableInput::from_toml_str(&text)
        .with_context(|| format!("{}: invalid input", a.input.display()))?;
    let tables =
        compute(&input, base).with_context(|| format!("{}: invalid input", a.input.display()))?;
    for t in &tables {
        for w in &t.warnings {
            eprintln!("warning: {w}");
        }
    }
    let json = table_json(&input, &tables)?;
    std::fs::create_dir_all(&a.out).with_context(|| format!("creating {}", a.out.display()))?;
    let set = OutputSet::begin(&a.out, JSON_FILE, json)?;
    set.write(CSV_FILE, table_csv(&tables, &input.units))?;
    set.publish()?;
    eprintln!(
        "{} tables, {} energies each; wrote {}",
        tables.len(),
        tables[0].rows.len(),
        a.out.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stopping_units_round_trip() {
        let (n, rho) = (4.9e28, 2329.0);
        for u in [
            StoppingUnit::EvPerNm,
            StoppingUnit::KevPerUm,
            StoppingUnit::MevPerMgCm2,
            StoppingUnit::EvPer1e15AtomsCm2,
        ] {
            let s = 7.3e-35;
            let back = u.to_si(u.from_si(s, n, rho), n, rho);
            assert!((back / s - 1.0).abs() < 1e-13, "{u:?}");
        }
        for u in [RangeUnit::Nm, RangeUnit::Um, RangeUnit::MgCm2] {
            let back = u.to_si(u.from_si(1.5e-7, rho), rho);
            assert!((back / 1.5e-7 - 1.0).abs() < 1e-13, "{u:?}");
        }
    }

    #[test]
    fn units_agree_with_each_other() {
        // 100 eV/nm in Si (n = 5e28 m^-3, rho = 2329 kg/m^3).
        let (n, rho) = (5.0e28, 2329.0);
        let s = StoppingUnit::EvPerNm.to_si(100.0, n, rho);
        assert!((StoppingUnit::KevPerUm.from_si(s, n, rho) - 100.0).abs() < 1e-9);
        // 100 eV/nm = 100e7 eV/cm = 1e9 eV/cm; per g/cm^2 /2.329; MeV/(mg/cm^2): 1e9/2.329 eV
        // per g/cm^2 = 1e9/2.329/1e6 MeV per 1e3 mg/cm^2.
        let want = 1e9 / 2.329 / 1e6 / 1e3;
        assert!((StoppingUnit::MevPerMgCm2.from_si(s, n, rho) / want - 1.0).abs() < 1e-12);
        // 1e15 atoms/cm^2 of Si is 1e15/5e22 cm = 2e-8 cm = 0.2 nm: 20 eV.
        assert!((StoppingUnit::EvPer1e15AtomsCm2.from_si(s, n, rho) - 20.0).abs() < 1e-9);
        // 1 um of Si is 2329e-6 * 100 mg/cm^2... 1e-4 cm * 2.329 g/cm^3 = 0.2329 mg/cm^2.
        assert!((RangeUnit::MgCm2.from_si(1e-6, rho) - 0.2329).abs() < 1e-12);
    }

    #[test]
    fn grid_forms() {
        let log = EnergySpec {
            values_ev: None,
            start_ev: Some(10.0),
            stop_ev: Some(1e4),
            points: Some(4),
            scale: None,
        };
        let g = log.grid().unwrap();
        assert_eq!(g.len(), 4);
        assert!((g[1] / 100.0 - 1.0).abs() < 1e-12 && g[3] == 1e4);
        let one = EnergySpec {
            points: Some(1),
            stop_ev: None,
            ..log.clone()
        };
        assert_eq!(one.grid().unwrap(), vec![10.0]);
        let bad = EnergySpec {
            values_ev: Some(vec![5.0, 5.0]),
            start_ev: None,
            stop_ev: None,
            points: None,
            scale: None,
        };
        assert!(bad.grid().unwrap_err().to_string().contains("values_ev[1]"));
    }
}
