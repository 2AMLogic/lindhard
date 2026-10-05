//! The `lindhard` command: one TOML input in, `summary.json` and CSV profiles
//! out. Input schema and output layout: `docs/cli.md`.

mod output;
mod tally;

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use lindhard::input::{stopping_model, Input, Resolved};
use lindhard::ion::bca::Bca;
use lindhard::ion::potential::Potential;
use lindhard::ion::scattering::ScatteringTable;

use lindhard::tally::IonTally;

use crate::tally::{ion_tally_config, CliTally};

const LONG_VERSION: &str = concat!(
    env!("CARGO_PKG_VERSION"),
    " (",
    env!("LINDHARD_GIT_DESCRIBE"),
    ")"
);

/// Monte Carlo transport of ions and electrons in matter.
#[derive(Parser)]
#[command(name = "lindhard", version = LONG_VERSION, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Parse and validate an input file without running it.
    Check {
        /// Input TOML file.
        input: PathBuf,
    },
    /// Run an input file and write summary.json and CSV profiles.
    Run {
        /// Input TOML file.
        input: PathBuf,
        /// Output directory (created if missing; files in it are overwritten).
        #[arg(long, short)]
        out: PathBuf,
        /// Override `run.ions` (echoed in the output).
        #[arg(long)]
        ions: Option<u64>,
        /// Override `run.seed` (echoed in the output).
        #[arg(long)]
        seed: Option<u64>,
        /// Override `run.threads`. Never changes the results.
        #[arg(long)]
        threads: Option<usize>,
    },
}

fn load(path: &Path) -> Result<Input> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    Input::from_toml_str(&text).with_context(|| format!("{}: invalid input", path.display()))
}

fn resolve(path: &Path, input: &Input) -> Result<Resolved> {
    let r = input
        .resolve()
        .with_context(|| format!("{}: invalid input", path.display()))?;
    for w in &r.warnings {
        eprintln!("warning: {w}");
    }
    Ok(r)
}

fn check(path: &Path) -> Result<()> {
    let input = load(path)?;
    let r = resolve(path, &input)?;
    let b = &r.input.beam;
    println!("{}: OK", path.display());
    println!(
        "  beam: {} at {} eV, tilt {} deg, azimuth {} deg; {} ions, seed {}",
        b.ion, b.energy_ev, b.tilt_deg, b.azimuth_deg, r.input.run.ions, r.input.run.seed
    );
    for (i, (g, l)) in r.stack.layers().iter().zip(&r.layers).enumerate() {
        let extent = if g.back_m().is_finite() {
            format!("{} nm", g.thickness_m() * 1e9)
        } else {
            "semi-infinite".to_string()
        };
        println!(
            "  layer {i}: {} ({}), {extent}, {:.4e} atoms/cm^3",
            l.source,
            l.material
                .components()
                .iter()
                .map(|c| format!(
                    "{} {:.4}",
                    lindhard::elements::element(c.z()).expect("valid").symbol,
                    c.atom_fraction()
                ))
                .collect::<Vec<_>>()
                .join(", "),
            g.material().atom_number_density() * 1e-6
        );
    }
    for m in r.models() {
        println!("  {}: {}", m.role, m.name);
    }
    Ok(())
}

fn run(
    path: &Path,
    out: &Path,
    ions: Option<u64>,
    seed: Option<u64>,
    threads: Option<usize>,
) -> Result<()> {
    let mut input = load(path)?;
    if let Some(n) = ions {
        input.run.ions = n;
    }
    if let Some(s) = seed {
        input.run.seed = s;
    }
    if threads.is_some() {
        input.run.threads = threads;
    }
    let r = resolve(path, &input)?;
    if input.run.threads == Some(0) {
        bail!("--threads must be at least 1");
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(input.run.threads.unwrap_or(0))
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

    let stopping = stopping_model(r.input.physics.stopping);
    let bca = Bca::new(r.beam, &r.stack, r.config, &*stopping, &table)
        .context("setting up the transport engine")?;
    let tally_spec = &r.input.tally;
    let ion_proto = IonTally::new(&r.stack, &bca.species_z(), ion_tally_config(&r)?)
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

    std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    let write = |name: &str, text: String| -> Result<()> {
        let p = out.join(name);
        std::fs::write(&p, text).with_context(|| format!("writing {}", p.display()))
    };
    let info = output::RunInfo {
        threads: pool.current_num_threads(),
        table_build_s,
        transport_s,
        ions_per_s: tally.summary.histories as f64 / transport_s.max(f64::MIN_POSITIVE),
    };
    let report = tally.ion.report(r.input.tally.dual_pearson);
    write(
        output::SUMMARY_FILE,
        output::summary_json(&r, &table, &tally, &report, info)?,
    )?;
    write(output::DEPTH_FILE, output::depth_csv(&tally))?;
    write(output::LATERAL_FILE, output::lateral_csv(&report))?;
    write(output::DAMAGE_FILE, output::damage_csv(&report))?;
    write(output::ESCAPES_FILE, output::escapes_csv(&report))?;
    if tally.per_ion {
        write(output::IONS_FILE, output::ions_csv(&tally))?;
    }
    let s = &tally.summary;
    eprintln!(
        "{} ions: {} stopped, {} backscattered, {} transmitted, {} sputtered atoms; wrote {}",
        s.histories,
        s.primaries_stopped,
        s.backscattered,
        s.transmitted,
        s.sputtered,
        out.display()
    );
    Ok(())
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Check { input } => check(&input),
        Command::Run {
            input,
            out,
            ions,
            seed,
            threads,
        } => run(&input, &out, ions, seed, threads),
    }
}
