//! The `lindhard` command: one TOML input in, `summary.json` and CSV profiles
//! out. Input schema and output layout: `docs/cli.md`.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use lindhard::input::{Input, Resolved};
use lindhard_cli::{dynamic, output, sim};

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
    // Relative [stopping] table paths are relative to the input file.
    let base = path.parent().filter(|p| !p.as_os_str().is_empty());
    let r = input
        .resolve_in(base.unwrap_or(Path::new(".")))
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
    if let Some(d) = &r.input.dynamic {
        println!(
            "  dynamic: {} ions/cm^2 in steps of {}{}",
            d.fluence_cm2,
            d.ions_per_step,
            d.max_change
                .map_or(String::new(), |c| format!(" (adaptive, max change {c})"))
        );
    }
    for m in r.models() {
        if m.name == "user-table" {
            // A user table is identified by its file: "path: provenance".
            println!("  {}: {} ({})", m.role, m.name, m.citation);
        } else {
            println!("  {}: {}", m.role, m.name);
        }
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

    if r.input.dynamic.is_some() {
        return run_dynamic(&r, input.run.threads, out);
    }

    let sim::Simulation {
        tally,
        table,
        report,
        info,
    } = sim::simulate(&r, input.run.threads)?;

    std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    let write = |name: &str, text: String| -> Result<()> {
        let p = out.join(name);
        std::fs::write(&p, text).with_context(|| format!("writing {}", p.display()))
    };
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

/// A run with a `[dynamic]` section: the fluence-step time series and
/// composition profiles instead of the static profiles.
fn run_dynamic(r: &Resolved, threads: Option<usize>, out: &Path) -> Result<()> {
    let d = dynamic::simulate_dynamic(r, threads)?;
    std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    let write = |name: &str, text: String| -> Result<()> {
        let p = out.join(name);
        std::fs::write(&p, text).with_context(|| format!("writing {}", p.display()))
    };
    write(
        output::DYNAMIC_SUMMARY_FILE,
        output::dynamic_summary_json(r, &d)?,
    )?;
    write(output::DYNAMIC_STEPS_FILE, output::dynamic_steps_csv(&d))?;
    write(
        output::DYNAMIC_COMPOSITION_FILE,
        output::dynamic_composition_csv(&d),
    )?;
    let last = d.steps.last().expect("step 0");
    eprintln!(
        "{} ions in {} steps ({} attempts rejected): {} sputtered atoms, {} slabs left; wrote {}",
        last.ions_done,
        last.step,
        d.rejected_attempts,
        last.cumulative.sputtered_total(),
        last.slabs.len(),
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
