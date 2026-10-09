//! The `lindhard` command: one TOML input in, `summary.json` and CSV profiles
//! out. Input schema and output layout: `docs/cli.md`.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use lindhard::input::electron::{ElectronInput, ResolvedElectron};
use lindhard::input::{Input, Resolved};
use lindhard_cli::table_cache::TableCache;
use lindhard_cli::{dynamic, electron, output, sim};

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
    /// Run an input file and write summary.json and CSV profiles
    /// (electron_summary.json and electron_*.csv for an `[electron]` input).
    Run {
        /// Input TOML file.
        input: PathBuf,
        /// Output directory (created if missing; files in it are overwritten).
        #[arg(long, short)]
        out: PathBuf,
        /// Override `run.ions`, or `run.histories` of an electron run
        /// (echoed in the output).
        #[arg(long, visible_alias = "histories")]
        ions: Option<u64>,
        /// Override `run.seed` (echoed in the output).
        #[arg(long)]
        seed: Option<u64>,
        /// Override `run.threads`. Never changes the results.
        #[arg(long)]
        threads: Option<usize>,
        /// Electron runs only: read the cross-section tables from this
        /// directory when it holds them for exactly this physics, grid and
        /// build, and store the tables built otherwise. Never changes the
        /// results (docs/cli.md, "Cross-section table cache").
        #[arg(long, value_name = "DIR")]
        table_cache: Option<PathBuf>,
    },
}

fn read(path: &Path) -> Result<String> {
    std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))
}

fn load(path: &Path, text: &str) -> Result<Input> {
    Input::from_toml_str(text).with_context(|| format!("{}: invalid input", path.display()))
}

/// The input file's directory: relative data paths resolve against it.
fn base_dir(path: &Path) -> &Path {
    path.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}

fn load_electron(path: &Path, text: &str) -> Result<ElectronInput> {
    ElectronInput::from_toml_str(text).with_context(|| format!("{}: invalid input", path.display()))
}

fn resolve_electron(path: &Path, input: &ElectronInput) -> Result<ResolvedElectron> {
    let r = input
        .resolve_in(base_dir(path))
        .with_context(|| format!("{}: invalid input", path.display()))?;
    for w in &r.warnings {
        eprintln!("warning: {w}");
    }
    Ok(r)
}

fn check_electron(path: &Path, text: &str) -> Result<()> {
    let input = load_electron(path, text)?;
    let r = resolve_electron(path, &input)?;
    let b = &r.input.electron.beam;
    println!("{}: OK (electron run)", path.display());
    println!(
        "  beam: electrons at {} eV, tilt {} deg, azimuth {} deg; {} histories, seed {}",
        b.energy_ev, b.tilt_deg, b.azimuth_deg, r.input.run.histories, r.input.run.seed
    );
    for (i, (g, &m)) in r.stack.layers().iter().zip(&r.layer_material).enumerate() {
        let extent = if g.back_m().is_finite() {
            format!("{} nm", g.thickness_m() * 1e9)
        } else {
            "semi-infinite".to_string()
        };
        println!("  layer {i}: {}, {extent}", r.materials[m].name);
    }
    for m in &r.materials {
        println!(
            "  {}: optical ELF {} ({})",
            m.name,
            m.optical_elf_file.path,
            m.optical_elf.provenance()
        );
    }
    let g = &r.table_energy_ev;
    println!(
        "  tables: {} energies, {} to {} eV",
        g.len(),
        g[0],
        g[g.len() - 1]
    );
    for m in r.models() {
        println!("  {}: {}", m.role, m.name);
    }
    Ok(())
}

/// Writes an optional output file when the current run produced it, and
/// otherwise removes the file left in a reused output directory by an earlier
/// run. A missing file is fine; any other failure is an error naming the path.
/// Only the one reserved file is touched, never the directory.
fn reconcile_optional(out: &Path, name: &str, text: Option<String>) -> Result<()> {
    let p = out.join(name);
    match text {
        Some(t) => std::fs::write(&p, t).with_context(|| format!("writing {}", p.display())),
        None => match std::fs::remove_file(&p) {
            Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
                Err(e).with_context(|| format!("removing stale {}", p.display()))
            }
            _ => Ok(()),
        },
    }
}

fn run_electron(
    path: &Path,
    text: &str,
    out: &Path,
    histories: Option<u64>,
    seed: Option<u64>,
    threads: Option<usize>,
    table_cache: Option<&Path>,
) -> Result<()> {
    let mut input = load_electron(path, text)?;
    if let Some(n) = histories {
        input.run.histories = n;
    }
    if let Some(s) = seed {
        input.run.seed = s;
    }
    if threads.is_some() {
        input.run.threads = threads;
    }
    let r = resolve_electron(path, &input)?;
    let cache = table_cache.map(TableCache::open).transpose()?;
    let sim = electron::simulate_electron(&r, input.run.threads, cache.as_ref())?;
    if let Some(c) = &cache {
        let read = sim
            .tables
            .iter()
            .flat_map(|t| [&t.elastic_origin, &t.inelastic_origin])
            .filter(|o| o.source == lindhard_cli::table_cache::TableSource::Cache)
            .count();
        eprintln!(
            "tables: {read} read from the cache, {} built ({:.1} s), cache {}",
            2 * sim.tables.len() - read,
            sim.info.table_build_s,
            c.dir().display()
        );
    }
    std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    let write = |name: &str, text: String| -> Result<()> {
        let p = out.join(name);
        std::fs::write(&p, text).with_context(|| format!("writing {}", p.display()))
    };
    let summary = electron::summary_json(&r, &sim)?;
    write(electron::SPECTRA_FILE, electron::spectra_csv(&sim.report))?;
    write(electron::TABLES_FILE, electron::tables_csv(&r, &sim))?;
    reconcile_optional(
        out,
        electron::CARTESIAN_FILE,
        electron::cartesian_csv(&sim.report),
    )?;
    reconcile_optional(
        out,
        electron::CYLINDRICAL_FILE,
        electron::cylindrical_csv(&sim.report),
    )?;
    reconcile_optional(
        out,
        electron::PSF_PROFILE_FILE,
        electron::psf_profile_csv(&sim),
    )?;
    reconcile_optional(
        out,
        electron::PSF_PARAMETERS_FILE,
        electron::psf_parameters_csv(&sim),
    )?;
    // Last, so the summary describes the completed output set.
    write(electron::SUMMARY_FILE, summary)?;
    let y = &sim.report.yields;
    eprintln!(
        "{} electrons: eta {:.4}, delta {:.4}, energy balance residual {:.2e}; wrote {}",
        sim.report.histories,
        y.backscatter_eta,
        y.secondary_delta,
        sim.report.budget.relative_imbalance,
        out.display()
    );
    Ok(())
}

fn resolve(path: &Path, input: &Input) -> Result<Resolved> {
    // Relative [stopping] table paths are relative to the input file.
    let r = input
        .resolve_in(base_dir(path))
        .with_context(|| format!("{}: invalid input", path.display()))?;
    for w in &r.warnings {
        eprintln!("warning: {w}");
    }
    Ok(r)
}

fn check(path: &Path) -> Result<()> {
    let text = read(path)?;
    if ElectronInput::is_electron_toml(&text) {
        return check_electron(path, &text);
    }
    let input = load(path, &text)?;
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
    for (i, c) in r.input.crystal.iter().enumerate() {
        println!(
            "  crystal[{i}]: {} on layers {:?}, normal {:?}, reference {:?}, wafer rotation {} deg{}",
            c.preset.name(),
            c.layers,
            c.normal,
            c.reference,
            c.wafer_rotation_deg,
            c.thermal
                .as_ref()
                .map_or(", static lattice".to_string(), |t| format!(
                    ", thermal {} K",
                    t.temperature_k
                ))
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
    table_cache: Option<&Path>,
) -> Result<()> {
    let text = read(path)?;
    if ElectronInput::is_electron_toml(&text) {
        return run_electron(path, &text, out, ions, seed, threads, table_cache);
    }
    if table_cache.is_some() {
        bail!("--table-cache applies to electron runs only (an input with an [electron] table)");
    }
    let mut input = load(path, &text)?;
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
        crystals,
        info,
    } = sim::simulate(&r, input.run.threads)?;

    std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    let write = |name: &str, text: String| -> Result<()> {
        let p = out.join(name);
        std::fs::write(&p, text).with_context(|| format!("writing {}", p.display()))
    };
    let summary = output::summary_json(&r, &table, &tally, &report, &crystals, info)?;
    write(output::DEPTH_FILE, output::depth_csv(&tally))?;
    write(output::LATERAL_FILE, output::lateral_csv(&report))?;
    write(output::DAMAGE_FILE, output::damage_csv(&report))?;
    write(output::ESCAPES_FILE, output::escapes_csv(&report))?;
    reconcile_optional(
        out,
        output::IONS_FILE,
        tally.per_ion.then(|| output::ions_csv(&tally)),
    )?;
    // Last, so the summary describes the completed output set.
    write(output::SUMMARY_FILE, summary)?;
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
            table_cache,
        } => run(&input, &out, ions, seed, threads, table_cache.as_deref()),
    }
}
