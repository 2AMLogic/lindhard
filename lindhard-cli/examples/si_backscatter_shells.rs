//! Does an inner-shell channel close the Si backscatter excess? The
//! library-level comparison of #449 (part of #148; `docs/validation.md`,
//! "Backscatter coefficient η(E, Z)", subsection "Inner-shell channels for
//! Si").
//!
//! ```text
//! cargo run --release -p lindhard-cli --example si_backscatter_shells -- \
//!     validation/experiments/backscatter/eta_si.toml \
//!     [--table-cache DIR] [--energies 10,20,30] [--seeds 1,2] \
//!     [--histories 100000] [--threads 2] [--determinism]
//! ```
//!
//! The input is the committed Si input of the η comparison; at each energy
//! only `electron.beam.energy_ev`, `run.histories` and `run.seed` change, as
//! in `validation/experiments/backscatter.py`. `lindhard run` builds no
//! inner-shell tables (#156), so the example builds them with the library
//! and runs the transport itself. Four variants, the same elastic table,
//! seeds, primaries and table grid in each:
//!
//! - **A**, the committed baseline: full Penn on the Yang et al. (2019)
//!   ELF, no shells. Built and run as `lindhard run` does; the example stops
//!   before the other variants if its η at 10, 20 or 30 keV (seed 1, 100000
//!   primaries) differs from the committed `fast_count` of
//!   `validation/experiments/backscatter_results.json`.
//! - **B**, the matched control: the single-pole Penn algorithm on the same
//!   ELF, no shells. The single-pole algorithm is what the shell channels
//!   support (`build_shell_channel_tables`), so A → B is the change of the
//!   inelastic model and B → C, D is what the shells add.
//! - **C**: B plus the K shell (EADL2017 binding energy, 2 electrons,
//!   hydrogenic ELF of `shell_elf`).
//! - **D**: C plus the L1, L2 and L3 hydrogenic ELFs **cut below the upper
//!   end of the valence ELF** (199 eV). The valence ELF already holds the L
//!   edges and everything up to 199 eV, so only what lies above is added.
//!   The cut is a modelling choice of this example: the hydrogenic ELF at
//!   the cut is kept as it is, and the valence ELF's own last point is not
//!   adjusted. Without the cut the L shells would count 99 to 199 eV twice.
//!
//! Secondaries stay off in every variant, so the shell effect is not mixed
//! with the secondary-generation effect (a shell channel then removes
//! `W` from the primary and deposits it: nothing is emitted).
//!
//! Printed per run: η, its binomial standard error `sqrt(η (1 - η) / N)`
//! (each primary escapes at most once with secondaries off), the residual
//! against the measured median of the committed groups, and the verdict
//! against the unchanged bound. Then, per energy, the seeds pooled and the
//! A → B and B → C, D differences in units of their standard error. With
//! `--determinism`, variant C at the first energy is also run on one and on
//! two threads and the example fails unless the two tallies are equal.
//!
//! Nothing here is fitted: no parameter is adjusted to the Si data.

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use lindhard::electron::data::{
    CrossSectionTable, OpticalElf, ShellChannelTable, Subshell, SubshellBindingTable,
};
use lindhard::electron::inelastic::table::{
    build_shell_channel_tables, stopping_power_ev_per_m, InelasticTableOptions,
};
use lindhard::electron::inelastic::{
    hydrogenic_shell_elf, InnerShell, ShellElfGrid, ShellResolvedChannels, SinglePolePenn,
};
use lindhard::electron::transport::{LayerTables, Transport};
use lindhard::input::electron::{ElectronInput, ResolvedElectron};
use lindhard::tally::{ElectronReport, FullElectronTally};
use lindhard_cli::electron::{inelastic_axis, rate_refinement, simulate_electron, CHUNK_SIZE};
use lindhard_cli::table_cache::TableCache;

/// The unchanged pass bound on |η − measured median| (docs/validation.md).
const TOLERANCE: f64 = 0.05;
/// The graded Z of this example.
const Z_SI: u8 = 14;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Variant {
    A,
    B,
    C,
    D,
}

impl Variant {
    const ALL: [Variant; 4] = [Variant::A, Variant::B, Variant::C, Variant::D];

    fn label(self) -> &'static str {
        match self {
            Variant::A => "A baseline (full Penn, no shells)",
            Variant::B => "B single-pole, no shells",
            Variant::C => "C single-pole + K",
            Variant::D => "D single-pole + K + L1-3 above the cut",
        }
    }
}

struct Args {
    input: PathBuf,
    cache: Option<PathBuf>,
    energies_kev: Vec<f64>,
    seeds: Vec<u64>,
    histories: u64,
    threads: usize,
    determinism: bool,
}

fn parse_args() -> Result<Args> {
    let usage = "usage: si_backscatter_shells INPUT [--table-cache DIR] [--energies 10,20,30] \
                 [--seeds 1,2] [--histories N] [--threads N] [--determinism]";
    let mut it = std::env::args().skip(1);
    let input = PathBuf::from(it.next().context(usage)?);
    let mut a = Args {
        input,
        cache: None,
        energies_kev: vec![10.0, 20.0, 30.0],
        seeds: vec![1, 2],
        histories: 100_000,
        threads: 2,
        determinism: false,
    };
    while let Some(flag) = it.next() {
        if flag == "--determinism" {
            a.determinism = true;
            continue;
        }
        let v = it
            .next()
            .with_context(|| format!("{flag} needs a value; {usage}"))?;
        match flag.as_str() {
            "--table-cache" => a.cache = Some(PathBuf::from(v)),
            "--energies" => a.energies_kev = list(&v).context("--energies")?,
            "--seeds" => a.seeds = list(&v).context("--seeds")?,
            "--histories" => a.histories = v.parse().context("--histories")?,
            "--threads" => a.threads = v.parse().context("--threads")?,
            _ => bail!("unknown flag {flag}; {usage}"),
        }
    }
    if a.threads == 0 || a.histories == 0 || a.energies_kev.is_empty() || a.seeds.is_empty() {
        bail!("threads, histories, energies and seeds must be non-empty and positive");
    }
    Ok(a)
}

fn list<T: std::str::FromStr>(s: &str) -> Result<Vec<T>>
where
    T::Err: std::error::Error + Send + Sync + 'static,
{
    s.split(',').map(|x| Ok(x.trim().parse()?)).collect()
}

/// The committed input with the beam energy, history count and seed set;
/// each line must occur exactly once.
fn variant_input(base: &str, energy_ev: f64, histories: u64, seed: u64) -> Result<String> {
    let mut t = base.to_string();
    for (key, value) in [
        ("energy_ev", format!("{energy_ev:?}")),
        ("histories", histories.to_string()),
        ("seed", seed.to_string()),
    ] {
        let lines: Vec<&str> = t
            .lines()
            .filter(|l| l.split('=').next().is_some_and(|k| k.trim() == key) && !l.starts_with('#'))
            .collect();
        let [line] = lines.as_slice() else {
            bail!(
                "expected one `{key} = ...` line in the input, found {}",
                lines.len()
            );
        };
        t = t.replacen(line, &format!("{key} = {value}"), 1);
    }
    Ok(t)
}

/// `elf` restricted to energies from `cut_ev` up: the value at the cut by
/// the ELF's own linear interpolation, then every later point. The cut must
/// lie inside the ELF's range.
fn cut_below(elf: &OpticalElf, cut_ev: f64) -> Result<OpticalElf> {
    let (lo, hi) = elf.energy_range_ev();
    if !(cut_ev >= lo && cut_ev < hi) {
        bail!("cut {cut_ev} eV is outside the ELF range {lo}..{hi} eV");
    }
    let mut energy = vec![cut_ev];
    let mut values = vec![elf.elf(cut_ev)?];
    for (&w, &e) in elf.energy_ev().iter().zip(elf.elf_values()) {
        if w > cut_ev {
            energy.push(w);
            values.push(e);
        }
    }
    Ok(OpticalElf::new(
        elf.material(),
        format!(
            "{}; kept from {cut_ev} eV up (below it the valence ELF is used instead, #449)",
            elf.provenance()
        ),
        energy,
        values,
    )?)
}

/// The shell channels of `variant` as (shell, ELF) pairs for target atom
/// density `n_si` (m⁻³).
fn shells(
    variant: Variant,
    binding: &SubshellBindingTable,
    material: &str,
    n_si: f64,
    cut_ev: f64,
) -> Result<Vec<(InnerShell, OpticalElf)>> {
    let grid = ShellElfGrid::default();
    let mut out = Vec::new();
    if matches!(variant, Variant::C | Variant::D) {
        out.push(hydrogenic_shell_elf(
            material,
            binding,
            Z_SI,
            Subshell::K,
            n_si,
            grid,
        )?);
    }
    if variant == Variant::D {
        for label in ["L1", "L2", "L3"] {
            let s = Subshell::from_label(label).context("subshell label")?;
            let (shell, elf) = hydrogenic_shell_elf(material, binding, Z_SI, s, n_si, grid)?;
            out.push((shell, cut_below(&elf, cut_ev)?));
        }
    }
    Ok(out)
}

/// The tallied η of one run, with the report.
struct Run {
    eta: f64,
    report: ElectronReport,
    seconds: f64,
}

fn run_transport(
    r: &ResolvedElectron,
    elastic: &CrossSectionTable,
    inelastic: &CrossSectionTable,
    shell_tables: Vec<ShellChannelTable>,
    threads: usize,
) -> Result<Run> {
    let mut transport = Transport::new(
        r.stack.clone(),
        vec![LayerTables {
            elastic: elastic.clone(),
            inelastic: inelastic.clone(),
        }],
        r.config,
    )?;
    if !shell_tables.is_empty() {
        transport = transport.with_inner_shells(0, shell_tables)?;
    }
    let proto = FullElectronTally::new(&transport, r.tally)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    let t0 = Instant::now();
    let run = pool.install(|| {
        transport.run(
            r.input.run.seed,
            r.input.run.histories,
            CHUNK_SIZE,
            &r.primary,
            || proto.clone(),
        )
    })?;
    let report = run.tally.report();
    Ok(Run {
        eta: report.yields.backscatter_eta,
        report,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

/// What the committed results hold for Si: the measured medians per energy
/// and the baseline `fast_count` per energy.
struct Committed {
    median: Vec<(f64, f64)>,
    fast_count: Vec<(f64, u64)>,
    histories: u64,
    seed: u64,
}

fn committed(path: &Path) -> Result<Committed> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let v: serde_json::Value = serde_json::from_str(&text)?;
    let si = v["targets"]
        .as_array()
        .and_then(|t| t.iter().find(|x| x["target"] == "Si"))
        .context("no Si target in the committed results")?;
    let median = si["groups"]
        .as_array()
        .context("no groups")?
        .iter()
        .filter_map(|g| Some((g["energy_kev"].as_f64()?, g["median"].as_f64()?)))
        .collect();
    let fast_count = si["runs"]
        .as_object()
        .context("no runs")?
        .iter()
        .filter_map(|(k, r)| Some((k.parse().ok()?, r["fast_count"].as_u64()?)))
        .collect();
    Ok(Committed {
        median,
        fast_count,
        histories: v["histories"].as_u64().context("histories")?,
        seed: v["seed"].as_u64().context("seed")?,
    })
}

fn lookup<T: Copy>(table: &[(f64, T)], e_kev: f64) -> Option<T> {
    table
        .iter()
        .find(|(k, _)| (k - e_kev).abs() < 1e-9)
        .map(|(_, v)| *v)
}

fn sigma(eta: f64, n: f64) -> f64 {
    (eta * (1.0 - eta) / n).sqrt()
}

fn main() -> Result<()> {
    let args = parse_args()?;
    let base = std::fs::read_to_string(&args.input)
        .with_context(|| format!("reading {}", args.input.display()))?;
    let dir = args.input.parent().unwrap_or(Path::new("."));
    let results = committed(&dir.join("../backscatter_results.json"))?;
    let resolve = |e_kev: f64, seed: u64| -> Result<ResolvedElectron> {
        Ok(ElectronInput::from_toml_str(&variant_input(
            &base,
            1000.0 * e_kev,
            args.histories,
            seed,
        )?)?
        .resolve_in(dir)?)
    };

    // The grid is the same for every energy, so one build serves all runs:
    // the baseline tables (elastic and full-Penn inelastic) as `lindhard run`
    // builds them, through the same table cache.
    let r0 = resolve(args.energies_kev[0], args.seeds[0])?;
    if r0.materials.len() != 1 || r0.layer_material.len() != 1 {
        bail!("the input must have one material in one layer");
    }
    let cache = args.cache.as_ref().map(TableCache::open).transpose()?;
    let mut tables_only = r0.clone();
    tables_only.input.run.histories = 1;
    let t0 = Instant::now();
    let sim = simulate_electron(&tables_only, Some(args.threads), cache.as_ref())?;
    println!(
        "baseline tables (elastic, full Penn): {:.0} s",
        t0.elapsed().as_secs_f64()
    );
    let elastic = sim.tables[0].elastic.clone();
    let full_penn = sim.tables[0].inelastic.clone();

    // The single-pole valence channel and the shell channels, built once.
    let m = &r0.materials[0];
    let (fermi_ev, axis) = inelastic_axis(&r0, m);
    let valence = SinglePolePenn::new(m.optical_elf.clone()).with_fermi_energy_ev(fermi_ev)?;
    let (_, valence_hi) = m.optical_elf.energy_range_ev();
    let binding = SubshellBindingTable::eadl2017();
    let n_si = m.material.atom_number_density();
    let mut options = InelasticTableOptions::new(r0.table_energy_ev.clone()).with_axis(axis);
    options.rate_refinement = rate_refinement(&r0, m);
    let name = m.optical_elf.material().to_string();
    let mut built: Vec<(Variant, CrossSectionTable, Vec<ShellChannelTable>)> = Vec::new();
    for v in Variant::ALL {
        if v == Variant::A {
            built.push((v, full_penn.clone(), Vec::new()));
            continue;
        }
        let t0 = Instant::now();
        let channels = ShellResolvedChannels::new(
            valence.clone(),
            shells(v, &binding, &name, n_si, valence_hi)?,
        )?;
        let tables =
            build_shell_channel_tables(&channels, &m.material, &options, binding.provenance())?;
        println!(
            "{}: tables built in {:.1} s ({} shell channels)",
            v.label(),
            t0.elapsed().as_secs_f64(),
            tables.shells.len()
        );
        built.push((v, tables.valence, tables.shells));
    }

    // Coverage and share of the stopping power, per channel.
    let grid = elastic.energy_ev();
    println!("\nstopping-power share of each channel (single-pole rows of variant D):");
    let d = built
        .iter()
        .find(|(v, ..)| *v == Variant::D)
        .context("variant D")?;
    for &e_kev in &args.energies_kev {
        let i = grid
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1 - 1000.0 * e_kev)
                    .abs()
                    .total_cmp(&(b.1 - 1000.0 * e_kev).abs())
            })
            .map(|(i, _)| i)
            .context("empty grid")?;
        let sv = stopping_power_ev_per_m(&d.1, i).context("valence row")?;
        let mut parts = vec![("valence".to_string(), sv)];
        for s in &d.2 {
            parts.push((
                format!(
                    "{} (B = {:.1} eV)",
                    s.subshell().label(),
                    s.binding_energy_ev()
                ),
                stopping_power_ev_per_m(s.table(), i).context("shell row")?,
            ));
        }
        let total: f64 = parts.iter().map(|p| p.1).sum();
        let line: Vec<String> = parts
            .iter()
            .map(|(n, s)| format!("{n} {:.3} eV/nm ({:.2} %)", s * 1e-9, 100.0 * s / total))
            .collect();
        println!(
            "  row {:.0} eV (grid ends {:.0} eV): {}",
            grid[i],
            grid[grid.len() - 1],
            line.join(", ")
        );
    }

    // The runs.
    let mut rows: Vec<(Variant, f64, u64, f64)> = Vec::new(); // variant, E, seed, eta
    println!(
        "\n{:<44} {:>6} {:>5} {:>9} {:>8} {:>9} {:>7}  {:>6}",
        "variant", "E/keV", "seed", "eta", "sigma", "resid", "<=0.05", "run s"
    );
    for &e in &args.energies_kev {
        for &seed in &args.seeds {
            let r = resolve(e, seed)?;
            let n = r.input.run.histories as f64;
            let median = lookup(&results.median, e);
            for (v, inelastic, shell_tables) in &built {
                let run =
                    run_transport(&r, &elastic, inelastic, shell_tables.clone(), args.threads)?;
                if *v == Variant::A && seed == results.seed && args.histories == results.histories {
                    if let Some(want) = lookup(&results.fast_count, e) {
                        let got = (run.eta * n).round() as u64;
                        if got != want {
                            bail!(
                                "baseline at {e} keV: {got} fast electrons, the committed run has {want}; \
                                 stopping before any variant is judged"
                            );
                        }
                    }
                }
                let resid = median.map(|m| run.eta - m);
                println!(
                    "{:<44} {:>6} {:>5} {:>9.5} {:>8.5} {:>9} {:>7}  {:>6.0}",
                    v.label(),
                    e,
                    seed,
                    run.eta,
                    sigma(run.eta, n),
                    resid.map_or("n/a".into(), |x| format!("{x:+.4}")),
                    resid.map_or("n/a", |x| if x.abs() <= TOLERANCE { "yes" } else { "NO" }),
                    run.seconds
                );
                let _ = &run.report;
                rows.push((*v, e, seed, run.eta));
            }
        }
    }

    // Pooled over the seeds, and the differences.
    let pooled = |v: Variant, e: f64| -> Option<(f64, f64)> {
        let etas: Vec<f64> = rows
            .iter()
            .filter(|r| r.0 == v && r.1 == e)
            .map(|r| r.3)
            .collect();
        if etas.is_empty() {
            return None;
        }
        let eta = etas.iter().sum::<f64>() / etas.len() as f64;
        Some((eta, sigma(eta, args.histories as f64 * etas.len() as f64)))
    };
    println!(
        "\npooled over {} seed(s); differences in units of their combined standard error:",
        args.seeds.len()
    );
    println!(
        "{:>6} {:<44} {:>9} {:>8} {:>10} {:>8}",
        "E/keV", "variant", "eta", "sigma", "vs ref", "sigmas"
    );
    for &e in &args.energies_kev {
        for (v, reference) in [
            (Variant::A, None),
            (Variant::B, Some(Variant::A)),
            (Variant::C, Some(Variant::B)),
            (Variant::D, Some(Variant::B)),
        ] {
            let (eta, s) = pooled(v, e).context("pooled")?;
            let diff = reference.and_then(|rv| {
                pooled(rv, e).map(|(re, rs)| (rv, eta - re, (s * s + rs * rs).sqrt()))
            });
            println!(
                "{:>6} {:<44} {:>9.5} {:>8.5} {:>10} {:>8}",
                e,
                v.label(),
                eta,
                s,
                diff.map_or("-".into(), |(rv, d, _)| format!("{d:+.5} ({rv:?})")),
                diff.map_or("-".into(), |(_, d, ds)| format!("{:+.1}", d / ds)),
            );
        }
    }

    if args.determinism {
        let e = args.energies_kev[0];
        let r = resolve(e, args.seeds[0])?;
        let (_, inelastic, shell_tables) =
            built.iter().find(|(v, ..)| *v == Variant::C).context("C")?;
        let one = run_transport(&r, &elastic, inelastic, shell_tables.clone(), 1)?;
        let two = run_transport(&r, &elastic, inelastic, shell_tables.clone(), 2)?;
        if one.report != two.report {
            bail!("variant C at {e} keV differs between 1 and 2 threads");
        }
        println!(
            "\ndeterminism: variant C at {e} keV is identical on 1 and 2 threads (eta {:.5})",
            one.eta
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "[electron.beam]\nenergy_ev = 10000.0\ntilt_deg = 0.0\n\n\
        [electron.tables]\nmin_energy_ev = 51.0\nmax_energy_ev = 31000.0\n\n\
        [run]\nhistories = 100000\nseed = 1\n";

    #[test]
    fn variant_input_changes_only_the_three_lines() {
        let t = variant_input(INPUT, 20000.0, 50, 7).unwrap();
        assert!(t.contains("energy_ev = 20000.0\n"));
        assert!(t.contains("histories = 50\n"));
        assert!(t.contains("seed = 7\n"));
        assert!(t.contains("min_energy_ev = 51.0") && t.contains("max_energy_ev = 31000.0"));
        assert_eq!(t.lines().count(), INPUT.lines().count());
    }

    #[test]
    fn variant_input_rejects_a_missing_line() {
        assert!(variant_input("[run]\nseed = 1\n", 1.0, 1, 1).is_err());
    }

    #[test]
    fn cut_below_keeps_the_elf_above_the_cut_and_nothing_below() {
        let elf = OpticalElf::new(
            "x",
            "test",
            vec![1.0, 2.0, 4.0, 8.0],
            vec![0.1, 0.2, 0.4, 0.8],
        )
        .unwrap();
        let c = cut_below(&elf, 3.0).unwrap();
        assert_eq!(c.energy_ev(), &[3.0, 4.0, 8.0]);
        assert!((c.elf_values()[0] - 0.3).abs() < 1e-12);
        assert_eq!(&c.elf_values()[1..], &[0.4, 0.8]);
        // an energy at an existing point is not duplicated
        assert_eq!(cut_below(&elf, 2.0).unwrap().energy_ev(), &[2.0, 4.0, 8.0]);
        assert!(cut_below(&elf, 8.0).is_err() && cut_below(&elf, 0.5).is_err());
    }

    /// The overlap guard: no L-shell ELF of variant D has a point inside the
    /// valence ELF's range, and the K shell starts far above it.
    #[test]
    fn variant_d_adds_nothing_below_the_valence_elf_end() {
        let binding = SubshellBindingTable::eadl2017();
        let valence = OpticalElf::from_toml_file(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../validation/data/optical/si_elf_yang2019.toml"),
        )
        .unwrap();
        let (_, hi) = valence.energy_range_ev();
        let n_si = 4.99e28;
        for v in [Variant::C, Variant::D] {
            for (shell, elf) in shells(v, &binding, "Si", n_si, hi).unwrap() {
                let first = elf.energy_range_ev().0;
                assert!(
                    first >= hi,
                    "{} starts at {first} eV, valence ends at {hi} eV",
                    shell.subshell.label()
                );
            }
        }
        assert_eq!(
            shells(Variant::B, &binding, "Si", n_si, hi).unwrap().len(),
            0
        );
        assert_eq!(
            shells(Variant::C, &binding, "Si", n_si, hi).unwrap().len(),
            1
        );
        assert_eq!(
            shells(Variant::D, &binding, "Si", n_si, hi).unwrap().len(),
            4
        );
    }
}
