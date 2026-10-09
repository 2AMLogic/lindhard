//! The measurements of issue #173: why the default (single-pole Penn)
//! inelastic model overestimates the secondary-electron yield δ of Al and Au.
//!
//! ```text
//! cargo run --release -p lindhard-cli --example inelastic_low_energy -- INPUT MODE
//! ```
//!
//! `INPUT` is an electron input with one material, a band and the
//! Kieft-Bosch secondary model, such as the inputs that
//! `validation/experiments/se_yield.py` writes. `MODE` is one of:
//!
//! - `imfp`: the inelastic mean free path of the single-pole, full Penn and
//!   Mermin models of the input's optical ELF at a few energies from 12 eV to
//!   1 keV (model Fermi energy 0, as the tables of a default input). Full
//!   Penn costs minutes per material.
//! - `clamp`: run the input with its own tables and report, per bin of the
//!   energy above the Fermi level `E - E_F` before the event, the inelastic
//!   events per primary, the fraction that hit the transport's `E - E_F`
//!   clamp, the mean loss, and δ (front-face escapes below 50 eV per primary).
//! - `fermi-reference`: as `clamp`, with the inelastic table rebuilt in the
//!   convention of the cstool table compiler (rows at `T = E - E_F` above
//!   the Fermi level with the model's Fermi energy set to the band's, so the
//!   kinematics use `T' = E` and the losses stop at `E - E_F`; stored on the
//!   band-bottom energy axis `E` the transport reads). No event then reaches
//!   the clamp.
//! - `splice:MODEL:CUT`: as `clamp`, with the rows of the inelastic table
//!   below the band-bottom energy `CUT` (eV) replaced by the rows of `MODEL`
//!   (`mermin` or `single-pole`) on the same energy grid. This attributes a
//!   δ difference between two models to an energy range.
//!
//! Runs use two threads and seed 1. Nothing here is a model: the spliced and
//! rebuilt tables are diagnostics, not options of the library.

use std::path::Path;

use anyhow::{bail, Context, Result};
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts};
use lindhard::electron::inelastic::table::{
    build_inelastic_table_for_model, InelasticTableOptions,
};
use lindhard::electron::inelastic::{PennAlgorithm, PennInelastic};
use lindhard::electron::transport::{
    ElectronState, ElectronTally, Face, Fate, LayerTables, Transport,
};
use lindhard::input::electron::{ElectronInput, ResolvedElectron};

/// Bins of `E - E_F` before an inelastic event, eV.
const BINS: [f64; 9] = [0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, f64::INFINITY];

#[derive(Clone, Default)]
struct ClampTally {
    fermi_ev: f64,
    events: [u64; 8],
    clamped: [u64; 8],
    loss_ev: [f64; 8],
    slow: u64,
    histories: u64,
}

impl ElectronTally for ClampTally {
    fn inelastic(&mut self, after: &ElectronState, w_ev: f64) {
        let above = after.energy_ev + w_ev - self.fermi_ev;
        let b = BINS
            .windows(2)
            .position(|x| above >= x[0] && above < x[1])
            .unwrap_or(0);
        self.events[b] += 1;
        self.loss_ev[b] += w_ev;
        // A clamped loss leaves the primary exactly at the Fermi level.
        if (after.energy_ev - self.fermi_ev).abs() <= 1e-9 * self.fermi_ev.max(1.0) {
            self.clamped[b] += 1;
        }
    }
    fn escaped(&mut self, at: &ElectronState, face: Face) {
        if face == Face::Front && at.energy_ev < 50.0 {
            self.slow += 1;
        }
    }
    fn end_history(&mut self, _index: u64, _fate: Fate) {
        self.histories += 1;
    }
    fn merge(&mut self, o: Self) {
        for i in 0..self.events.len() {
            self.events[i] += o.events[i];
            self.clamped[i] += o.clamped[i];
            self.loss_ev[i] += o.loss_ev[i];
        }
        self.slow += o.slow;
        self.histories += o.histories;
    }
}

fn table_of(r: &ResolvedElectron, alg: PennAlgorithm, fermi_ev: f64) -> Result<CrossSectionTable> {
    let m = &r.materials[0];
    let model =
        PennInelastic::try_new(alg, m.optical_elf.clone())?.with_fermi_energy_ev(fermi_ev)?;
    Ok(build_inelastic_table_for_model(
        &model,
        &m.material,
        &InelasticTableOptions::new(r.table_energy_ev.clone()),
    )?)
}

fn with_rows(
    t: &CrossSectionTable,
    energy_ev: Vec<f64>,
    inverse_mfp_per_m: Vec<f64>,
    probability: Vec<f64>,
    quantiles: Vec<Vec<f64>>,
    what: &str,
) -> Result<CrossSectionTable> {
    Ok(CrossSectionTable::new(CrossSectionTableParts {
        model: format!("{} [diagnostic: {what}]", t.model()),
        material: t.material().to_string(),
        provenance: format!("{} [diagnostic: {what}]", t.provenance()),
        axis: t.axis(),
        energy_ev,
        inverse_mfp_per_m,
        probability,
        quantiles,
    })?)
}

/// The table in the cstool convention (module docs, `fermi-reference`).
fn fermi_reference(r: &ResolvedElectron, fermi_ev: f64) -> Result<CrossSectionTable> {
    let t = table_of(r, r.inelastic, fermi_ev)?;
    let n = t.energy_ev().len();
    let q = (0..n)
        .map(|i| t.quantiles(i).map(<[f64]>::to_vec).unwrap_or_default())
        .collect();
    with_rows(
        &t,
        t.energy_ev().iter().map(|e| e + fermi_ev).collect(),
        t.inverse_mfp_per_m().to_vec(),
        t.probability().to_vec(),
        q,
        "rows at E - E_F with the band Fermi energy, on the band-bottom axis",
    )
}

/// Rows below `cut_ev` from `low`, the rest from `high`, on the union of
/// their probability grids (each row resampled with the table's own linear
/// interpolation, so sampling is unchanged).
fn splice(
    high: &CrossSectionTable,
    low: &CrossSectionTable,
    cut_ev: f64,
) -> Result<CrossSectionTable> {
    let mut p: Vec<f64> = high
        .probability()
        .iter()
        .chain(low.probability())
        .copied()
        .collect();
    p.sort_by(f64::total_cmp);
    p.dedup();
    let grid = high.energy_ev().to_vec();
    let (mut inv, mut q) = (Vec::new(), Vec::new());
    for (i, &e) in grid.iter().enumerate() {
        let t = if e < cut_ev { low } else { high };
        inv.push(t.inverse_mfp_per_m()[i]);
        q.push(match t.quantiles(i) {
            Some(_) => p
                .iter()
                .map(|&u| t.inverse_cdf(i, u))
                .collect::<Result<_, _>>()?,
            None => Vec::new(),
        });
    }
    with_rows(
        high,
        grid,
        inv,
        p,
        q,
        &format!("rows below {cut_ev} eV spliced in"),
    )
}

fn imfp(r: &ResolvedElectron) -> Result<()> {
    let elf = &r.materials[0].optical_elf;
    let models = [
        PennAlgorithm::SinglePole,
        PennAlgorithm::Full,
        PennAlgorithm::Mermin,
    ]
    .into_iter()
    .map(|a| PennInelastic::try_new(a, elf.clone()))
    .collect::<Result<Vec<_>, _>>()?;
    println!(
        "{:>8} {:>12} {:>12} {:>12}",
        "E / eV", "spa / nm", "full / nm", "mermin / nm"
    );
    for e in [12.0, 15.0, 19.0, 24.0, 30.0, 50.0, 100.0, 300.0, 1000.0] {
        let mut row = format!("{e:8.1}");
        for m in &models {
            row += &format!(" {:12.3}", m.imfp_and_stopping(e)?.imfp_m() * 1e9);
        }
        println!("{row}");
    }
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        bail!("usage: inelastic_low_energy INPUT MODE (module docs)");
    }
    let path = Path::new(&args[1]);
    let mode = args[2].as_str();
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let r =
        ElectronInput::from_toml_str(&text)?.resolve_in(path.parent().unwrap_or(Path::new(".")))?;
    if r.materials.len() != 1 || r.layer_material.len() != 1 {
        bail!("the input must have one material in one layer");
    }
    if mode == "imfp" {
        return imfp(&r);
    }
    let band = r.materials[0]
        .band
        .clone()
        .context("the input's material needs a band")?;
    let fermi_ev = band.fermi_ev();
    // The input's own tables, built as `lindhard run` builds them.
    let mut tables_only = r.clone();
    tables_only.input.run.histories = 1;
    let sim = lindhard_cli::electron::simulate_electron(&tables_only, Some(2))?;
    let own = sim.tables[0].inelastic.clone();
    let inelastic = match mode {
        "clamp" => own,
        "fermi-reference" => fermi_reference(&r, fermi_ev)?,
        _ => match mode
            .strip_prefix("splice:")
            .map(|s| s.split(':').collect::<Vec<_>>())
        {
            Some(v) if v.len() == 2 => {
                let alg = match v[0] {
                    "mermin" => PennAlgorithm::Mermin,
                    "single-pole" => PennAlgorithm::SinglePole,
                    other => bail!("unknown splice model {other}"),
                };
                let cut: f64 = v[1].parse().context("splice cut, eV")?;
                splice(&own, &table_of(&r, alg, r.inelastic_fermi_ev)?, cut)?
            }
            _ => bail!("unknown mode {mode} (module docs)"),
        },
    };
    let transport = Transport::with_band_structures(
        r.stack.clone(),
        vec![LayerTables {
            elastic: sim.tables[0].elastic.clone(),
            inelastic,
        }],
        vec![band],
        r.config,
    )?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build()?;
    let d = pool
        .install(|| {
            transport.run(1, r.input.run.histories, 16, &r.primary, || ClampTally {
                fermi_ev,
                ..Default::default()
            })
        })?
        .tally;
    let n = d.histories as f64;
    println!(
        "{mode}: E_F {fermi_ev:.4} eV, {} histories, delta {:.3}",
        d.histories,
        d.slow as f64 / n
    );
    println!("  E - E_F / eV   events/primary   clamped   mean loss / eV");
    for i in 0..d.events.len() {
        if d.events[i] > 0 {
            println!(
                "  [{:>4}, {:>4})   {:14.2}   {:7.3}   {:14.2}",
                BINS[i],
                BINS[i + 1],
                d.events[i] as f64 / n,
                d.clamped[i] as f64 / d.events[i] as f64,
                d.loss_ev[i] / d.events[i] as f64
            );
        }
    }
    Ok(())
}
