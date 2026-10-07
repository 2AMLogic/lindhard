//! Run a resolved input: build the scattering table, set up the engine and
//! tallies, and transport every history on a rayon pool of the requested size.
//! Shared by the CLI and the Python bindings so both produce the same numbers.

use std::time::Instant;

use anyhow::{bail, Context, Result};
use lindhard::input::Resolved;
use lindhard::ion::bca::Bca;
use lindhard::ion::potential::Potential;
use lindhard::ion::scattering::ScatteringTable;
use lindhard::tally::{IonReport, IonTally};

use crate::output::RunInfo;
use crate::tally::{ion_tally_config, CliTally};

/// Everything a finished run produced.
pub struct Simulation {
    /// The merged tallies.
    pub tally: CliTally,
    /// The scattering table the run used.
    pub table: ScatteringTable,
    /// The ion tally's report.
    pub report: IonReport,
    /// Thread count and timings (the only thread-dependent part of the output).
    pub info: RunInfo,
}

/// Run `r` on `threads` workers (`None`: all available cores). The thread
/// count never changes the results.
pub fn simulate(r: &Resolved, threads: Option<usize>) -> Result<Simulation> {
    if threads == Some(0) {
        bail!("threads must be at least 1");
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.unwrap_or(0))
        .build()
        .context("building the thread pool")?;

    // The angle table depends on the screening function and length only (the
    // engine works in reduced variables), so the Z pair here is immaterial.
    let t0 = Instant::now();
    let z2 = r.layers[0].material.components()[0].z();
    let pot = Potential::new(r.screening, f64::from(r.beam.ion.z()), f64::from(z2))
        .with_length(r.screening_length);
    let table = ScatteringTable::build(&pot, &r.table_spec);
    let table_build_s = t0.elapsed().as_secs_f64();

    let stopping = r.stopping_model();
    let bca = Bca::new(r.beam, &r.stack, r.config, &*stopping, &table)
        .context("setting up the transport engine")?;
    let tally_spec = &r.input.tally;
    let ion_proto = IonTally::new(&r.stack, &bca.species_z(), ion_tally_config(r)?)
        .context("setting up the ion tally")?;
    let t1 = Instant::now();
    let tally = pool
        .install(|| {
            bca.run(|| {
                CliTally::new(
                    tally_spec.depth_bin_nm * 1e-9,
                    tally_spec.depth_bins,
                    tally_spec.per_ion,
                    ion_proto.clone(),
                )
            })
        })
        .context("transport failed")?;
    let transport_s = t1.elapsed().as_secs_f64();
    let report = tally.ion.report(r.input.tally.dual_pearson);
    let info = RunInfo {
        threads: pool.current_num_threads(),
        table_build_s,
        transport_s,
        ions_per_s: tally.summary.histories as f64 / transport_s.max(f64::MIN_POSITIVE),
    };
    Ok(Simulation {
        tally,
        table,
        report,
        info,
    })
}
