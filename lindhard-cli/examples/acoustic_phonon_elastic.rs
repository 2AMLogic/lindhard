//! A measurement of issue #242: what replacing the Mott elastic cross
//! section of slow electrons by the acoustic-phonon mean free path of
//! Verduin (2017) does to the secondary-electron yield δ.
//!
//! ```text
//! cargo run --release -p lindhard-cli --example acoustic_phonon_elastic -- INPUT MODE
//! ```
//!
//! `INPUT` is an electron input with one material, a band and the
//! Kieft-Bosch secondary model, such as the inputs that
//! `validation/experiments/se_yield.py` writes. `MODE` is one of:
//!
//! - `mott`: the control. The input's own tables, so δ equals that of
//!   `lindhard run` on the same input.
//! - `acoustic`: the elastic table is replaced below 200 eV as described
//!   under "The model". Only materials with a row in [`VERDUIN_TABLE_3_2`]
//!   (Al, Au) can be run.
//! - `imfp:E1,E2,...`: the inelastic mean free path of the input's
//!   inelastic model at the given energies (eV), with the model Fermi energy
//!   of the input. No run.
//!
//! One JSON object is printed on standard output. The run uses the input's
//! seed and history count and two threads.
//!
//! # The model
//!
//! T. Verduin, *Quantum Noise Effects in e-Beam Lithography and Metrology*,
//! PhD thesis, TU Delft (2017), doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c,
//! Section 3.4, pp. 87-98 (open access, read 2026-10-09). The thesis follows
//! Schreiber and Fitting; those papers were not opened.
//!
//! - **Mean free path.** Above the Brillouin-zone energy `E_BZ` the
//!   acoustic-phonon mean free path is the second line of Eq. 3.130 (p. 94),
//!   which is Eq. 3.127 (p. 92) written with the minimum of the curve:
//!   `λ(E) = λ_min x / (α (1 + α)) · [(1 + α) ln(1 + α) - α] / [ln(1 + x) - x / (1 + x)]`,
//!   `x = α E / E_min = E / A`, with `α = 2.162581587` (printed below
//!   Eq. 3.128) and `A` the screening parameter. `λ_min`, `E_min` and the
//!   cross-over energy `E_co = E_BZ / 4` are the entries of Table 3.2
//!   (p. 96), which has rows for Al, Si, Au and SiO2 only.
//! - **Angle.** Eq. 3.126 (p. 92), second line: the distribution in
//!   `s = (1 - cos θ) / 2` is proportional to `x s / (1 + x s)²`. Its
//!   integral, `[ln(1 + x s) - x s / (1 + x s)] / [ln(1 + x) - x / (1 + x)]`,
//!   is our own algebra on that equation; at `s = 1` it reproduces the
//!   energy dependence of Eq. 3.127.
//! - **Where it is used.** P. 97: the acoustic-phonon mean free path is used
//!   exclusively below 100 eV, the Mott cross section exclusively above
//!   200 eV, and "in the intermediate range, we linearly interpolate the
//!   scattering cross-section of both models". This is read here as
//!   `λ⁻¹ = (1 - w) λ_AC⁻¹ + w λ_Mott⁻¹` with `w = (E - 100 eV) / 100 eV`;
//!   the thesis does not print the formula. The angle between 100 and
//!   200 eV is drawn from the two distributions in proportion to their
//!   weighted rates; the thesis does not say what it does there.
//! - **Energy.** The thesis plots these curves against the kinetic energy
//!   in the material, which is the energy above the band bottom that this
//!   transport carries (p. 98, "the kinetic energy increases from E to
//!   E + E_F + Φ").
//!
//! # What is not in it
//!
//! - **No energy loss.** The thesis gives the net phonon loss per event
//!   (Eq. 3.116, p. 89) and evaluates it for Si only (12.3 meV); it prints
//!   no dispersion coefficients for Al or Au. The events here are elastic.
//!   The run reports the number of elastic events below 100 eV per primary,
//!   which is the number of events such a loss would apply to.
//! - **Nothing below `E_BZ`.** The first line of Eq. 3.130 (below `E_co`)
//!   and the interpolation between `E_co` and `4 E_co`, which the thesis
//!   describes only as "linearly interpolated", are not implemented. The
//!   run refuses an input whose stopping threshold lies below
//!   `E_BZ = 4 E_co`, so no transported electron reads such a row; the rows
//!   below it are filled with the high-energy line and never sampled.
//!
//! Nothing here is a model of the library: the replaced table is a
//! diagnostic, as in `inelastic_low_energy.rs`.

use std::path::Path;

use anyhow::{bail, Context, Result};
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts};
use lindhard::electron::inelastic::PennInelastic;
use lindhard::electron::transport::{
    ElectronState, ElectronTally, Face, Fate, LayerTables, Transport,
};
use lindhard::input::electron::ElectronInput;

/// One row of Verduin (2017), Table 3.2, p. 96 ("The four parameters in the
/// phenomenological AC phonon model for various materials").
struct AcousticRow {
    /// Element symbol, as the input names the material.
    material: &'static str,
    /// `λ_min`, the minimum mean free path, Å.
    lambda_min_angstrom: f64,
    /// The energy of that minimum, eV.
    e_min_ev: f64,
    /// The cross-over energy `E_co = E_BZ / 4`, eV.
    e_co_ev: f64,
}

/// The metal rows of Verduin (2017), Table 3.2, p. 96, as printed: aluminum
/// `λ_min` 4.57 Å at 99.2 eV, `E_co` 2.30 eV; gold 6.62 Å at 97.7 eV, `E_co`
/// 2.26 eV. The thesis derives them from the Fermi level and the electrical
/// resistivity (Eq. 3.135, p. 95) and a screening parameter of five times
/// the Brillouin-zone energy (p. 96); it does not print those inputs. The
/// table has no row for Cu.
const VERDUIN_TABLE_3_2: [AcousticRow; 2] = [
    AcousticRow {
        material: "Al",
        lambda_min_angstrom: 4.57,
        e_min_ev: 99.2,
        e_co_ev: 2.30,
    },
    AcousticRow {
        material: "Au",
        lambda_min_angstrom: 6.62,
        e_min_ev: 97.7,
        e_co_ev: 2.26,
    },
];

/// `α`, the position of the minimum of the mean free path in units of the
/// screening parameter: Verduin (2017), printed below Eq. 3.128, p. 93.
const ALPHA: f64 = 2.162581587;

/// Below this energy the acoustic-phonon mean free path is used alone, eV
/// (Verduin 2017, p. 97).
const ACOUSTIC_ONLY_BELOW_EV: f64 = 100.0;

/// Above this energy the Mott cross section is used alone, eV (Verduin 2017,
/// p. 97).
const MOTT_ONLY_ABOVE_EV: f64 = 200.0;

/// Escaping electrons below this energy are secondaries, eV (the split of
/// `validation/experiments/se_yield.py`).
const SLOW_BELOW_EV: f64 = 50.0;

/// `g(y) = ln(1 + y) - y / (1 + y)`, the angular integral of Eq. 3.126 and
/// the energy dependence of Eq. 3.127. The series avoids the cancellation
/// at small `y`.
fn g(y: f64) -> f64 {
    if y < 1e-4 {
        y * y * (0.5 - 2.0 * y / 3.0)
    } else {
        y.ln_1p() - y / (1.0 + y)
    }
}

impl AcousticRow {
    /// The screening parameter `A = E_min / α`, eV.
    fn screening_ev(&self) -> f64 {
        self.e_min_ev / ALPHA
    }

    /// The Brillouin-zone energy `E_BZ = 4 E_co`, eV.
    fn brillouin_zone_ev(&self) -> f64 {
        4.0 * self.e_co_ev
    }

    /// The inverse mean free path above `E_BZ`, m⁻¹ (Eq. 3.130, second
    /// line, with `(1 + α) ln(1 + α) - α = (1 + α) g(α)`).
    fn inverse_mfp_per_m(&self, energy_ev: f64) -> f64 {
        let x = energy_ev / self.screening_ev();
        let lambda_min_m = self.lambda_min_angstrom * 1e-10;
        ALPHA * g(x) / (x * g(ALPHA) * lambda_min_m)
    }

    /// The cumulative distribution of the polar angle above `E_BZ`
    /// (integral of Eq. 3.126, second line).
    fn cdf(&self, energy_ev: f64, theta: f64) -> f64 {
        let x = energy_ev / self.screening_ev();
        let s = (0.5 * theta).sin().powi(2);
        g(x * s) / g(x)
    }
}

/// The forward CDF of a stored inverse-CDF row at `theta`: the inverse of
/// the table's own linear interpolation.
fn table_cdf(probability: &[f64], quantiles: &[f64], theta: f64) -> f64 {
    let j = quantiles.partition_point(|&q| q <= theta);
    if j == 0 {
        return 0.0;
    }
    if j == quantiles.len() {
        return 1.0;
    }
    let (q0, q1) = (quantiles[j - 1], quantiles[j]);
    probability[j - 1] + (theta - q0) / (q1 - q0) * (probability[j] - probability[j - 1])
}

/// The elastic table with the acoustic-phonon rows spliced in (module docs).
fn acoustic_table(mott: &CrossSectionTable, row: &AcousticRow) -> Result<CrossSectionTable> {
    // The union of the Mott grid and a uniform one: a Mott row resampled on
    // it samples exactly as before, and the smooth acoustic distribution
    // gets a grid that does not depend on where Mott needed refinement.
    let mut p: Vec<f64> = mott
        .probability()
        .iter()
        .copied()
        .chain((0..=256).map(|k| f64::from(k) / 256.0))
        .collect();
    p.sort_by(f64::total_cmp);
    p.dedup();
    let grid = mott.energy_ev().to_vec();
    let (mut inv, mut q) = (Vec::new(), Vec::new());
    for (i, &e) in grid.iter().enumerate() {
        let w = ((e - ACOUSTIC_ONLY_BELOW_EV) / (MOTT_ONLY_ABOVE_EV - ACOUSTIC_ONLY_BELOW_EV))
            .clamp(0.0, 1.0);
        let mott_row = mott.quantiles(i);
        let mott_rate = if mott_row.is_some() {
            w * mott.inverse_mfp_per_m()[i]
        } else {
            0.0
        };
        let acoustic_rate = (1.0 - w) * row.inverse_mfp_per_m(e);
        let total = mott_rate + acoustic_rate;
        inv.push(total);
        if total == 0.0 {
            q.push(Vec::new());
            continue;
        }
        if acoustic_rate == 0.0 {
            q.push(
                p.iter()
                    .map(|&u| mott.inverse_cdf(i, u))
                    .collect::<Result<_, _>>()?,
            );
            continue;
        }
        let cdf = |theta: f64| {
            let m = mott_row.map_or(0.0, |r| table_cdf(mott.probability(), r, theta));
            (acoustic_rate * row.cdf(e, theta) + mott_rate * m) / total
        };
        let mut out = Vec::with_capacity(p.len());
        let mut lo = 0.0_f64;
        for &u in &p {
            let theta = if u <= 0.0 {
                0.0
            } else if u >= 1.0 {
                std::f64::consts::PI
            } else {
                // Bisection from the previous quantile: the CDF is
                // non-decreasing, so the row is too.
                let (mut a, mut b) = (lo, std::f64::consts::PI);
                for _ in 0..200 {
                    let mid = 0.5 * (a + b);
                    if mid <= a || mid >= b {
                        break;
                    }
                    if cdf(mid) < u {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                0.5 * (a + b)
            };
            lo = theta;
            out.push(theta);
        }
        q.push(out);
    }
    let what = format!(
        "acoustic-phonon mean free path and angle of Verduin (2017) Table 3.2 and Eqs. 3.126, \
         3.130 below {ACOUSTIC_ONLY_BELOW_EV} eV (lambda_min {} angstrom at {} eV), mixed \
         linearly with the table's own rows up to {MOTT_ONLY_ABOVE_EV} eV; no energy loss",
        row.lambda_min_angstrom, row.e_min_ev
    );
    Ok(CrossSectionTable::new(CrossSectionTableParts {
        model: format!("{} [diagnostic: {what}]", mott.model()),
        material: mott.material().to_string(),
        provenance: format!("{} [diagnostic: {what}]", mott.provenance()),
        axis: mott.axis(),
        energy_ev: grid,
        inverse_mfp_per_m: inv,
        probability: p,
        quantiles: q,
    })?)
}

#[derive(Clone, Default)]
struct YieldTally {
    histories: u64,
    slow: u64,
    fast: u64,
    elastic: u64,
    elastic_below_100_ev: u64,
}

impl ElectronTally for YieldTally {
    fn elastic(&mut self, after: &ElectronState, _theta: f64) {
        self.elastic += 1;
        if after.energy_ev < ACOUSTIC_ONLY_BELOW_EV {
            self.elastic_below_100_ev += 1;
        }
    }
    fn escaped(&mut self, at: &ElectronState, face: Face) {
        if face == Face::Front {
            if at.energy_ev < SLOW_BELOW_EV {
                self.slow += 1;
            } else {
                self.fast += 1;
            }
        }
    }
    fn end_history(&mut self, _index: u64, _fate: Fate) {
        self.histories += 1;
    }
    fn merge(&mut self, o: Self) {
        self.histories += o.histories;
        self.slow += o.slow;
        self.fast += o.fast;
        self.elastic += o.elastic;
        self.elastic_below_100_ev += o.elastic_below_100_ev;
    }
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        bail!("usage: acoustic_phonon_elastic INPUT MODE (module docs)");
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
    let material = &r.materials[0];
    if let Some(list) = mode.strip_prefix("imfp:") {
        let model = PennInelastic::try_new(r.inelastic, material.optical_elf.clone())?
            .with_fermi_energy_ev(r.inelastic_fermi_ev)?;
        let mut rows = Vec::new();
        for e in list.split(',') {
            let e: f64 = e.parse().context("imfp energy, eV")?;
            rows.push(serde_json::json!({
                "energy_ev": e,
                "imfp_angstrom": model.imfp_and_stopping(e)?.imfp_m() * 1e10,
            }));
        }
        println!(
            "{}",
            serde_json::json!({
                "mode": "imfp",
                "material": material.name,
                "model_fermi_ev": r.inelastic_fermi_ev,
                "rows": rows,
            })
        );
        return Ok(());
    }
    let band = material
        .band
        .clone()
        .context("the input's material needs a band")?;
    // The input's own tables, built as `lindhard run` builds them.
    let mut tables_only = r.clone();
    tables_only.input.run.histories = 1;
    let sim = lindhard_cli::electron::simulate_electron(&tables_only, Some(2), None)?;
    let own = sim.tables[0].elastic.clone();
    let row = VERDUIN_TABLE_3_2
        .iter()
        .find(|a| a.material == material.name);
    let elastic = match mode {
        "mott" => own,
        "acoustic" => {
            let row = row.with_context(|| {
                format!(
                    "Verduin (2017) Table 3.2 has no row for {} (Al and Au only)",
                    material.name
                )
            })?;
            acoustic_table(&own, row)?
        }
        _ => bail!("unknown mode {mode} (module docs)"),
    };
    let transport = Transport::with_band_structures(
        r.stack.clone(),
        vec![LayerTables {
            elastic,
            inelastic: sim.tables[0].inelastic.clone(),
        }],
        vec![band],
        r.config,
    )?;
    let threshold_ev = transport.stopping_thresholds_ev()[0];
    if let ("acoustic", Some(row)) = (mode, row) {
        if threshold_ev < row.brillouin_zone_ev() {
            bail!(
                "the stopping threshold {threshold_ev} eV lies below E_BZ = {} eV, where the \
                 acoustic-phonon model of this example is not implemented (module docs)",
                row.brillouin_zone_ev()
            );
        }
    }
    let acoustic_mfp = match (mode, row) {
        ("acoustic", Some(row)) => serde_json::json!({
            "at_e_min": 1e10 / row.inverse_mfp_per_m(row.e_min_ev),
            "at_20_ev": 1e10 / row.inverse_mfp_per_m(20.0),
        }),
        _ => serde_json::Value::Null,
    };
    let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build()?;
    let d = pool
        .install(|| {
            transport.run(
                r.input.run.seed,
                r.input.run.histories,
                16,
                &r.primary,
                YieldTally::default,
            )
        })?
        .tally;
    let n = d.histories as f64;
    println!(
        "{}",
        serde_json::json!({
            "mode": mode,
            "material": material.name,
            "seed": r.input.run.seed,
            "histories": d.histories,
            "stopping_threshold_ev": threshold_ev,
            "slow": d.slow,
            "fast": d.fast,
            "delta": d.slow as f64 / n,
            "eta": d.fast as f64 / n,
            "elastic_per_primary": d.elastic as f64 / n,
            "elastic_below_100_ev_per_primary": d.elastic_below_100_ev as f64 / n,
            // The spliced-in mean free path at two energies, to check the
            // formula against the table by eye (null in the control).
            "acoustic_mfp_angstrom": acoustic_mfp,
        })
    );
    Ok(())
}
