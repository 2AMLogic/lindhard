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
//!   Since #241 `lindhard run` builds the inelastic table of a material with
//!   a band in the convention of the cstool table compiler (rows at
//!   `T = E - E_F` with the model's Fermi energy set to the band's, so the
//!   kinematics use `T' = E` and the losses stop at `E - E_F`; stored on the
//!   band-bottom energy axis `E` the transport reads), so no event reaches
//!   the clamp.
//! - `legacy`: as `clamp`, with the inelastic table built as `lindhard run`
//!   built it before #241: on the model's own axis with the model's Fermi
//!   energy `[electron.inelastic] fermi_energy_ev` (0), read at the
//!   band-bottom energy, so the losses reach `E` and the clamp acts.
//! - `splice:MODEL:CUT`: as `clamp`, with the rows of the inelastic table
//!   below the band-bottom energy `CUT` (eV) replaced by the rows of `MODEL`
//!   (`mermin` or `single-pole`) on the same energy grid and axis. This
//!   attributes a δ difference between two models to an energy range.
//! - `rate` (#291): no transport. The linear-interpolation error of the
//!   inelastic rate `1/λ` on the band-bottom axis just above the Fermi level,
//!   for the single-pole and Mermin models of the input's ELF, on the grid
//!   `lindhard run` builds (`log_energy_grid` from the input's grid ends,
//!   and from the library's default `min_energy_ev` 10 eV) at 10, 20, 40 and
//!   80 points per decade. The reference is the model's row at each
//!   evaluation energy itself (a table built at exactly those energies, the
//!   limit of a dense grid), so it carries no interpolation error. Only the
//!   coarse rows that bracket `E - E_F` = 0.5 to 20 eV are built. The
//!   interpolation is the transport's (`electron::transport`, "Free
//!   flight": linear in energy, held beyond the grid ends). Printed: the
//!   signed relative error `(interpolated - exact) / exact` at
//!   `E - E_F` = 0.5, 1, 2, 3, 5, 10 and 20 eV, and its largest magnitude
//!   over `E - E_F` = 2 to 20 eV in steps of 0.1 eV (with the `E - E_F` at
//!   which it occurs), gated against [`RATE_TOLERANCE`]. `-` marks a point
//!   whose exact rate is zero or whose row the table builder refuses (the
//!   refusal is printed); the gate covers the points with a positive exact
//!   rate. A coarse table the builder refuses is printed as refused.
//!
//! Runs use two threads and seed 1. Nothing here is a model: the spliced and
//! legacy tables are diagnostics, not options of the library.

use std::path::Path;

use anyhow::{bail, Context, Result};
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts};
use lindhard::electron::elastic::table::{log_energy_grid, DEFAULT_MIN_ENERGY_EV};
use lindhard::electron::inelastic::table::{
    build_inelastic_table_for_model, EnergyAxis, InelasticTableOptions,
};
use lindhard::electron::inelastic::{PennAlgorithm, PennInelastic};
use lindhard::electron::transport::{
    ElectronState, ElectronTally, Face, Fate, LayerTables, Transport,
};
use lindhard::input::electron::{ElectronInput, ResolvedElectron};
use rayon::prelude::*;

/// Bins of `E - E_F` before an inelastic event, eV.
const BINS: [f64; 9] = [0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, f64::INFINITY];

/// The tolerance of the `rate` mode, fixed before any value was measured
/// (#291): the largest relative error of the interpolated `1/λ` over
/// `E - E_F` = 2 to 20 eV must not exceed 1 %. It is the statistical error of
/// the most precise committed δ entry, rounded up: the Poisson floor of Al
/// `default` at 400 eV in `validation/experiments/se_yield_results.json` is
/// 0.0571 on δ = 6.52 (0.88 %; the other entries at 400 eV are 1.4 to 3.9 %).
/// A grid error below it cannot be told from the noise of those tables if δ
/// responds at most in proportion to the rate in that window; that response
/// is not measured here. The window starts at 2 eV, below the transport
/// cutoff of every `default` row (1 eV above the vacuum level, at least
/// 5 eV above `E_F`). Below 2 eV the error is printed, not gated: there
/// the rate starts from zero within about one coarse cell.
const RATE_TOLERANCE: f64 = 0.01;
/// `E - E_F` at which the `rate` mode prints the error, eV.
const RATE_POINTS_EV: [f64; 7] = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0];
/// Points per decade of the coarse grids of the `rate` mode.
const RATE_GRIDS: [f64; 4] = [10.0, 20.0, 40.0, 80.0];

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

fn table_of(
    r: &ResolvedElectron,
    alg: PennAlgorithm,
    fermi_ev: f64,
    axis: EnergyAxis,
) -> Result<CrossSectionTable> {
    let m = &r.materials[0];
    let model =
        PennInelastic::try_new(alg, m.optical_elf.clone())?.with_fermi_energy_ev(fermi_ev)?;
    Ok(build_inelastic_table_for_model(
        &model,
        &m.material,
        &InelasticTableOptions::new(r.table_energy_ev.clone()).with_axis(axis),
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

/// The rate of `t` at `e`, interpolated as the transport interpolates it:
/// linear in energy, held beyond the grid ends.
fn interpolated_rate(t: &CrossSectionTable, e: f64) -> f64 {
    let (g, r) = (t.energy_ev(), t.inverse_mfp_per_m());
    let n = g.len();
    if e <= g[0] {
        return r[0];
    }
    if e >= g[n - 1] {
        return r[n - 1];
    }
    let i = g.partition_point(|&x| x <= e) - 1;
    r[i] + (e - g[i]) / (g[i + 1] - g[i]) * (r[i + 1] - r[i])
}

fn rate(r: &ResolvedElectron) -> Result<()> {
    let m = &r.materials[0];
    let (fermi_ev, axis) = lindhard_cli::electron::inelastic_axis(r, m);
    if axis != EnergyAxis::BandBottom {
        bail!("the rate mode needs a material with a band");
    }
    let (lo_input, hi) = (
        r.table_energy_ev[0],
        r.table_energy_ev[r.table_energy_ev.len() - 1],
    );
    let mut lows = vec![lo_input];
    if DEFAULT_MIN_ENERGY_EV != lo_input {
        lows.push(DEFAULT_MIN_ENERGY_EV);
    }
    // E - E_F of the gate (2 to 20 eV in 0.1 eV steps) and of the printed points.
    let gate: Vec<f64> = (20..=200).map(|k| f64::from(k) * 0.1).collect();
    let mut above: Vec<f64> = gate.iter().chain(&RATE_POINTS_EV).copied().collect();
    above.sort_by(f64::total_cmp);
    above.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    let exact_e: Vec<f64> = above.iter().map(|t| fermi_ev + t).collect();
    let options = |e: Vec<f64>| InelasticTableOptions::new(e).with_axis(EnergyAxis::BandBottom);
    println!(
        "rate: {}, E_F {fermi_ev:.4} eV (band bottom), grid top {hi} eV, tolerance {:.0} % \
         over E - E_F = 2 to 20 eV",
        m.name,
        RATE_TOLERANCE * 100.0
    );
    for (name, alg) in [
        ("single-pole", PennAlgorithm::SinglePole),
        ("mermin", PennAlgorithm::Mermin),
    ] {
        let t0 = std::time::Instant::now();
        let model =
            PennInelastic::try_new(alg, m.optical_elf.clone())?.with_fermi_energy_ev(fermi_ev)?;
        // One table at every evaluation energy. If the builder refuses a
        // row, each energy is built again on its own (paired with the row at
        // E_F + 20 eV, as a table needs two rows), so that a refused row is
        // reported and left out instead of ending the mode.
        let rows: Vec<Result<f64, String>> =
            match build_inelastic_table_for_model(&model, &m.material, &options(exact_e.clone())) {
                Ok(t) => t.inverse_mfp_per_m().iter().map(|&x| Ok(x)).collect(),
                Err(_) => exact_e
                    .par_iter()
                    .map(|&e| {
                        let pair = if e < fermi_ev + 20.0 {
                            vec![e, fermi_ev + 20.0]
                        } else {
                            vec![fermi_ev + 19.0, e]
                        };
                        let k = usize::from(e >= fermi_ev + 20.0);
                        build_inelastic_table_for_model(&model, &m.material, &options(pair))
                            .map(|t| t.inverse_mfp_per_m()[k])
                            .map_err(|err| err.to_string())
                    })
                    .collect(),
            };
        for r in rows.iter().filter_map(|r| r.as_ref().err()) {
            println!("  {name}: row refused: {r}");
        }
        let exact: Vec<Option<f64>> = rows.into_iter().map(Result::ok).collect();
        let exact_at = |t: f64| {
            let i = above
                .iter()
                .position(|&a| (a - t).abs() < 1e-9)
                .expect("an evaluation point");
            exact[i]
        };
        let refused: Vec<f64> = above
            .iter()
            .zip(&exact)
            .filter(|(_, x)| x.is_none())
            .map(|(&t, _)| t)
            .collect();
        println!(
            "  {name}: exact 1/λ at E - E_F = 2, 5, 20 eV: {:.4e}, {:.4e}, {:.4e} 1/m; \
             rows refused at E - E_F = {refused:?} eV",
            exact_at(2.0).unwrap_or(f64::NAN),
            exact_at(5.0).unwrap_or(f64::NAN),
            exact_at(20.0).unwrap_or(f64::NAN)
        );
        let mut head = format!("  {:>12} {:>5} {:>5}", "min_energy", "ppd", "rows");
        for t in RATE_POINTS_EV {
            head += &format!(" {:>8}", format!("{t} eV"));
        }
        println!("{head} {:>9} {:>6} {:>5}", "max 2-20", "at eV", "gate");
        for &lo in &lows {
            for ppd in RATE_GRIDS {
                let grid = log_energy_grid(lo, hi, ppd)?;
                // The rows that bracket E_F + 0.5 .. E_F + 20 eV (the
                // transport's interpolant there depends on no other row).
                let first = grid
                    .partition_point(|&x| x <= fermi_ev + RATE_POINTS_EV[0])
                    .saturating_sub(1);
                let last = grid
                    .partition_point(|&x| x < fermi_ev + 20.0)
                    .min(grid.len() - 1);
                let sub = grid[first..=last].to_vec();
                let mut row = format!("  {:>9} eV {:>5} {:>5}", lo, ppd, sub.len());
                let coarse =
                    match build_inelastic_table_for_model(&model, &m.material, &options(sub)) {
                        Ok(t) => t,
                        Err(e) => {
                            println!("{row}  table refused: {e}");
                            continue;
                        }
                    };
                // The relative error, or `None` where the exact row is
                // refused or zero.
                let err = |t: f64| {
                    exact_at(t)
                        .filter(|&x| x > 0.0)
                        .map(|x| (interpolated_rate(&coarse, fermi_ev + t) - x) / x)
                };
                for t in RATE_POINTS_EV {
                    row += &match err(t) {
                        Some(v) => format!(" {:>7.2}%", 100.0 * v),
                        None => format!(" {:>8}", "-"),
                    };
                }
                // The largest error and its E - E_F (the first on a tie).
                let (worst, at) = gate
                    .iter()
                    .filter_map(|&t| err(t).map(|v| (v.abs(), t)))
                    .fold((0.0, gate[0]), |a, b| if b.0 > a.0 { b } else { a });
                println!(
                    "{row} {:>8.2}% {:>6.1} {:>5}",
                    100.0 * worst,
                    at,
                    if worst <= RATE_TOLERANCE {
                        "pass"
                    } else {
                        "FAIL"
                    }
                );
            }
        }
        println!("  ({name}: {:.1} s)", t0.elapsed().as_secs_f64());
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
    if mode == "rate" {
        return rate(&r);
    }
    let band = r.materials[0]
        .band
        .clone()
        .context("the input's material needs a band")?;
    let fermi_ev = band.fermi_ev();
    // The input's own tables, built as `lindhard run` builds them.
    let mut tables_only = r.clone();
    tables_only.input.run.histories = 1;
    let sim = lindhard_cli::electron::simulate_electron(&tables_only, Some(2), None)?;
    let own = sim.tables[0].inelastic.clone();
    let inelastic = match mode {
        "clamp" => own,
        "legacy" => table_of(
            &r,
            r.inelastic,
            r.inelastic_fermi_ev,
            EnergyAxis::ModelFermiLevel,
        )?,
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
                let (f, axis) = lindhard_cli::electron::inelastic_axis(&r, &r.materials[0]);
                splice(&own, &table_of(&r, alg, f, axis)?, cut)?
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
