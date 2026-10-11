//! Thread scaling and determinism of the electron engine on one benchmark
//! problem (#152; `docs/benchmarks.md`, "Electron engine").
//!
//! ```text
//! cargo run --release -p lindhard --example electron_scaling -- \
//!     --problem e_5keV_si --histories 2000 --threads 1,2,4 --repeat 3
//! cargo run --release -p lindhard --example electron_scaling -- \
//!     --input validation/oracle-runs/electron-bench/e_5keV_si/input.toml --threads 1,2
//! ```
//!
//! Builds the cross-section tables once (timed separately, on
//! `--table-threads` workers, default 1), then runs the same `--histories`
//! primaries with the same seed on an explicit rayon pool of each `--threads`
//! count, `--repeat` times each, and prints one JSON object: per thread count
//! the wall times of the transport, the median, electrons/s, the speedup over
//! the first count and the parallel efficiency (speedup times first count over
//! count), and the SHA-256 of the full tally report. The reports must be
//! bit-identical at every thread count (the engine's determinism invariant,
//! CONTRIBUTING.md); if any differs the program says so and exits with
//! status 1.
//!
//! `--problem` takes a problem id of `benches/electron_common` (the #150
//! matched problems) and runs its input: SYNTHETIC material data unless
//! `LINDHARD_BENCH_ELECTRON_INPUTS` names a directory of matched inputs.
//! `--input` runs any electron input file instead.

#[path = "../benches/electron_common/mod.rs"]
mod electron_common;

use std::path::PathBuf;
use std::time::Instant;

use electron_common as ec;

struct Args {
    problem: Option<String>,
    input: Option<PathBuf>,
    histories: u64,
    threads: Vec<usize>,
    table_threads: usize,
    repeat: usize,
}

fn usage() -> ! {
    eprintln!(
        "usage: electron_scaling (--problem ID | --input FILE) [--histories N] \
         [--threads 1,2,...] [--table-threads N] [--repeat R]"
    );
    std::process::exit(2);
}

fn parse() -> Args {
    let mut a = Args {
        problem: None,
        input: None,
        histories: 1000,
        threads: vec![1],
        table_threads: 1,
        repeat: 3,
    };
    let mut it = std::env::args().skip(1);
    while let Some(k) = it.next() {
        let mut v = || it.next().unwrap_or_else(|| usage());
        match k.as_str() {
            "--problem" => a.problem = Some(v()),
            "--input" => a.input = Some(PathBuf::from(v())),
            "--histories" => a.histories = v().parse().unwrap_or_else(|_| usage()),
            "--threads" => {
                a.threads = v()
                    .split(',')
                    .map(|s| s.trim().parse().unwrap_or_else(|_| usage()))
                    .collect()
            }
            "--table-threads" => a.table_threads = v().parse().unwrap_or_else(|_| usage()),
            "--repeat" => a.repeat = v().parse().unwrap_or_else(|_| usage()),
            _ => usage(),
        }
    }
    if a.problem.is_some() == a.input.is_some()
        || a.threads.is_empty()
        || a.threads.contains(&0)
        || a.table_threads == 0
        || a.repeat == 0
        || a.histories == 0
    {
        usage();
    }
    a
}

fn pool(n: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .expect("thread pool")
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

fn main() -> ec::Result<()> {
    let a = parse();
    let (r, source) = match (&a.problem, &a.input) {
        (Some(id), _) => {
            let p = ec::problem(id).ok_or_else(|| format!("unknown problem {id}"))?;
            let (r, matched) = ec::problem_input(&p, a.histories)?;
            (r, if matched { "matched" } else { "synthetic" })
        }
        (_, Some(path)) => (ec::resolve(path)?, "input file"),
        _ => unreachable!(),
    };

    let tp = pool(a.table_threads);
    let t0 = Instant::now();
    let mut elastic = Vec::new();
    for m in &r.materials {
        elastic.push(tp.install(|| ec::elastic_table(&r, &m.material))?);
    }
    let elastic_s = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let mut tables = Vec::new();
    for (m, el) in r.materials.iter().zip(elastic) {
        tables.push(lindhard::electron::transport::LayerTables {
            elastic: el,
            inelastic: tp.install(|| ec::inelastic_table(&r, m))?,
        });
    }
    let inelastic_s = t1.elapsed().as_secs_f64();
    let (transport, proto) = ec::transport(&r, &tables)?;

    let mut rows = Vec::new();
    let mut first: Option<(usize, f64, String)> = None;
    let mut all_identical = true;
    for &n in &a.threads {
        let p = pool(n);
        let mut walls = Vec::new();
        let mut digest = String::new();
        for _ in 0..a.repeat {
            let t = Instant::now();
            let report = p.install(|| ec::run(&r, &transport, &proto, a.histories))?;
            walls.push(t.elapsed().as_secs_f64());
            let d = ec::report_digest(&report)?;
            if !digest.is_empty() && d != digest {
                all_identical = false;
            }
            digest = d;
        }
        let med = median(&walls);
        let (n1, m1, d1) = first.get_or_insert((n, med, digest.clone())).clone();
        let identical = digest == d1;
        all_identical &= identical;
        let speedup = m1 / med;
        rows.push(serde_json::json!({
            "threads": n,
            "wall_s": walls,
            "wall_s_median": med,
            "electrons_per_s": a.histories as f64 / med,
            "speedup": speedup,
            "efficiency": speedup * n1 as f64 / n as f64,
            "report_sha256": digest,
            "identical_to_first": identical,
        }));
    }

    let out = serde_json::json!({
        "format": "lindhard-electron-scaling/1",
        "problem": a.problem,
        "input": a.input,
        "material_data": source,
        "energy_ev": r.primary.energy_ev,
        "materials": r.materials.iter().map(|m| m.name.clone()).collect::<Vec<_>>(),
        "histories": a.histories,
        "seed": ec::SEED,
        "chunk_size": ec::CHUNK_SIZE,
        "repeat": a.repeat,
        "tables": {
            "threads": a.table_threads,
            "elastic_s": elastic_s,
            "inelastic_s": inelastic_s,
            "energies": r.table_energy_ev.len(),
        },
        "runs": rows,
        "deterministic": all_identical,
    });
    println!("{}", serde_json::to_string_pretty(&out)?);
    if !all_identical {
        eprintln!("electron_scaling: the tally report differs between thread counts or repeats");
        std::process::exit(1);
    }
    Ok(())
}
