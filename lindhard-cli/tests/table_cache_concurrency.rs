//! Concurrent callers of the electron table cache (#305): one build per
//! cold entry across threads and across processes, no waiting between
//! different entries, and no dead holder blocking the next caller.
//!
//! The multi-process tests run this test binary again as the child process
//! (`child_process_entry`, selected with `--exact`, with its role in an
//! environment variable); run normally, that test does nothing. Builders
//! return a small synthetic table on the run's own energy grid, so the tests
//! exercise the cache's locking, lookup and publication, not the physics.

use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Barrier};
use std::time::{Duration, Instant};

use anyhow::Result;
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::input::electron::{ElectronInput, ResolvedElectron};
use lindhard_cli::table_cache::{
    key_text, sha256_hex, BuildId, CachedTable, TableCache, TableKind, TableSource,
};

/// Role of a child process (unset in the parent).
const CHILD_ROLE: &str = "LINDHARD_TABLE_CACHE_TEST_ROLE";
/// Working directory shared by the parent and its children.
const CHILD_DIR: &str = "LINDHARD_TABLE_CACHE_TEST_DIR";
/// Number of children of a `build-once` test.
const CHILD_COUNT: &str = "LINDHARD_TABLE_CACHE_TEST_COUNT";

/// Upper bound on every wait in these tests, so a broken lock fails a test
/// instead of hanging it.
const DEADLINE: Duration = Duration::from_secs(120);

fn example() -> ResolvedElectron {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../examples/electron");
    let text = std::fs::read_to_string(dir.join("e_10keV_si.toml")).unwrap();
    ElectronInput::from_toml_str(&text)
        .unwrap()
        .resolve_in(&dir)
        .unwrap()
}

fn scratch(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("lindhard-cli-table-cache-concurrency")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A fixed build identity, so every process (and thread) computes the same
/// keys without hashing the test executable.
fn cache_in(dir: &Path) -> TableCache {
    TableCache::with_build(
        dir.join("cache"),
        BuildId {
            lindhard_version: "0.0.0".into(),
            git_describe: "table-cache-concurrency-test".into(),
            executable_sha256: "00".into(),
        },
    )
    .unwrap()
}

/// A small valid table of `kind` on the run's energy grid.
fn synthetic_table(r: &ResolvedElectron, kind: TableKind) -> Result<CrossSectionTable> {
    let energy_ev = r.table_energy_ev.clone();
    // Each row runs from 0 to a value inside the axis's domain: 1 rad for
    // the polar angle, the incident energy for the energy loss.
    let (axis, quantiles) = match kind {
        TableKind::Elastic => (
            SamplingAxis::ElasticPolarAngle,
            energy_ev.iter().map(|_| vec![0.0, 1.0]).collect(),
        ),
        TableKind::Inelastic => (
            SamplingAxis::InelasticEnergyLoss,
            energy_ev.iter().map(|&e| vec![0.0, e]).collect(),
        ),
    };
    Ok(CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic (table cache concurrency test)".into(),
        material: "Si".into(),
        provenance: "synthetic test table, not physics".into(),
        axis,
        inverse_mfp_per_m: vec![1.0e9; energy_ev.len()],
        probability: vec![0.0, 1.0],
        quantiles,
        energy_ev,
    })?)
}

fn wait_until(what: &str, mut done: impl FnMut() -> bool) {
    let t0 = Instant::now();
    while !done() {
        assert!(t0.elapsed() < DEADLINE, "timed out waiting for {what}");
        std::thread::sleep(Duration::from_millis(10));
    }
}

fn files_in(dir: &Path) -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect();
    v.sort();
    v
}

fn count_prefixed(dir: &Path, prefix: &str) -> usize {
    files_in(dir)
        .iter()
        .filter(|f| f.starts_with(prefix))
        .count()
}

fn sha(t: &CachedTable) -> String {
    t.cache.as_ref().unwrap().sha256.clone()
}

/// Barrier-released threads asking for one cold entry: one builds, the rest
/// wait and read what it stored, and no rename collides.
#[test]
fn same_key_threads_build_once() {
    const N: usize = 8;
    let dir = scratch("threads");
    let cache = cache_in(&dir);
    let r = example();
    let m = &r.materials[0];
    let builds = AtomicUsize::new(0);
    let barrier = Barrier::new(N);
    let results: Vec<CachedTable> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..N)
            .map(|_| {
                s.spawn(|| {
                    barrier.wait();
                    cache.get_or_build(&r, m, TableKind::Inelastic, || {
                        builds.fetch_add(1, Ordering::SeqCst);
                        // Long enough for every other thread to reach the
                        // lock while this one holds it.
                        std::thread::sleep(Duration::from_millis(300));
                        synthetic_table(&r, TableKind::Inelastic)
                    })
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap().unwrap())
            .collect()
    });
    assert_eq!(builds.load(Ordering::SeqCst), 1, "the builder ran once");
    let built = results
        .iter()
        .filter(|t| t.source == TableSource::Built)
        .count();
    assert_eq!(built, 1);
    let first = &results[0];
    for t in &results {
        assert_eq!(t.table, first.table);
        assert_eq!(t.cache, first.cache);
    }
    let files = files_in(cache.dir());
    assert_eq!(files.len(), 4, "table, hash, key and lock: {files:?}");
    assert!(!files.iter().any(|f| f.contains(".tmp.")), "{files:?}");
}

/// Two different keys build at the same time: each builder waits (bounded)
/// until the other has started, which only happens if neither holds a lock
/// the other needs.
#[test]
fn different_keys_build_concurrently() {
    let dir = scratch("distinct-keys");
    let cache = cache_in(&dir);
    let r = example();
    let m = &r.materials[0];
    let started = AtomicUsize::new(0);
    let overlapped = AtomicUsize::new(0);
    std::thread::scope(|s| {
        for kind in [TableKind::Elastic, TableKind::Inelastic] {
            let (cache, r, started, overlapped) = (&cache, &r, &started, &overlapped);
            s.spawn(move || {
                let t = cache
                    .get_or_build(r, m, kind, || {
                        started.fetch_add(1, Ordering::SeqCst);
                        let t0 = Instant::now();
                        while started.load(Ordering::SeqCst) < 2 && t0.elapsed() < DEADLINE {
                            std::thread::sleep(Duration::from_millis(5));
                        }
                        if started.load(Ordering::SeqCst) == 2 {
                            overlapped.fetch_add(1, Ordering::SeqCst);
                        }
                        synthetic_table(r, kind)
                    })
                    .unwrap();
                assert_eq!(t.source, TableSource::Built);
            });
        }
    });
    assert_eq!(
        overlapped.load(Ordering::SeqCst),
        2,
        "both builders ran while the other was running"
    );
}

/// Starts this test binary as a child with `role` in `dir`.
fn spawn_child(role: &str, dir: &Path, count: usize) -> Child {
    Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "child_process_entry",
            "--nocapture",
            "--test-threads",
            "1",
        ])
        .env(CHILD_ROLE, role)
        .env(CHILD_DIR, dir)
        .env(CHILD_COUNT, count.to_string())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap()
}

/// The child side of the multi-process tests; a no-op unless started by
/// `spawn_child`.
#[test]
fn child_process_entry() {
    let Ok(role) = std::env::var(CHILD_ROLE) else {
        return;
    };
    let dir = PathBuf::from(std::env::var(CHILD_DIR).unwrap());
    let count: usize = std::env::var(CHILD_COUNT).unwrap().parse().unwrap();
    let pid = std::process::id();
    let cache = cache_in(&dir);
    let r = example();
    let m = &r.materials[0];
    match role.as_str() {
        // One of `count` processes asking for the same cold entry.
        "build-once" => {
            std::fs::write(dir.join(format!("calling-{pid}")), "").unwrap();
            let t = cache
                .get_or_build(&r, m, TableKind::Elastic, || {
                    // Record the build, then hold the lock until every child
                    // is at its door, and a little longer for each to block.
                    use std::io::Write;
                    let mut log = std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(dir.join("builds.log"))
                        .unwrap();
                    writeln!(log, "{pid}").unwrap();
                    drop(log);
                    wait_until("every child to call the cache", || {
                        count_prefixed(&dir, "calling-") == count
                    });
                    std::thread::sleep(Duration::from_millis(1000));
                    synthetic_table(&r, TableKind::Elastic)
                })
                .unwrap();
            let source = match t.source {
                TableSource::Built => "built",
                TableSource::Cache => "cache",
            };
            std::fs::write(
                dir.join(format!("result-{pid}")),
                format!("{source} {}", sha(&t)),
            )
            .unwrap();
        }
        // Takes the lock and never lets go: the parent kills it.
        "hold" => {
            let _ = cache.get_or_build(&r, m, TableKind::Elastic, || {
                std::fs::write(dir.join("holding"), "").unwrap();
                std::thread::sleep(Duration::from_secs(3600));
                unreachable!("the parent kills this process")
            });
        }
        other => panic!("unknown child role {other}"),
    }
}

/// Separate processes sharing one cold entry all succeed, and the build log
/// shows exactly one of them built the table; the rest waited for it (and
/// said so) and read it.
#[test]
fn separate_processes_build_once() {
    const N: usize = 4;
    let dir = scratch("processes");
    let children: Vec<Child> = (0..N).map(|_| spawn_child("build-once", &dir, N)).collect();
    let outputs: Vec<_> = children
        .into_iter()
        .map(|c| c.wait_with_output().unwrap())
        .collect();
    for o in &outputs {
        assert!(
            o.status.success(),
            "child failed: {:?}\nstdout:\n{}\nstderr:\n{}",
            o.status,
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        );
    }
    let log = std::fs::read_to_string(dir.join("builds.log")).unwrap();
    assert_eq!(log.lines().count(), 1, "one build across processes: {log}");
    let results: Vec<String> = files_in(&dir)
        .iter()
        .filter(|f| f.starts_with("result-"))
        .map(|f| std::fs::read_to_string(dir.join(f)).unwrap())
        .collect();
    assert_eq!(results.len(), N, "{results:?}");
    let built: Vec<_> = results.iter().filter(|s| s.starts_with("built ")).collect();
    assert_eq!(built.len(), 1, "{results:?}");
    let hash = built[0].trim_start_matches("built ");
    for s in &results {
        assert!(
            s.ends_with(hash),
            "every process got the same table: {results:?}"
        );
    }
    let waited = outputs
        .iter()
        .filter(|o| String::from_utf8_lossy(&o.stderr).contains("waiting for another process"))
        .count();
    assert_eq!(waited, N - 1, "every other process waited for the builder");
    let cache = dir.join("cache");
    assert!(!files_in(&cache).iter().any(|f| f.contains(".tmp.")));
}

/// A process killed while it holds an entry's lock does not block the next
/// caller, and the uncommitted files it (or an earlier crash) left are not
/// accepted: the next caller rebuilds the entry and clears the temporaries.
#[test]
fn killed_holder_does_not_block_the_next_caller() {
    let dir = scratch("killed-holder");
    let cache = cache_in(&dir);
    let r = example();
    let m = &r.materials[0];
    let key = key_text(cache.build(), &r, m, TableKind::Elastic);
    let paths = cache.paths(TableKind::Elastic, &sha256_hex(key.as_bytes()));

    let mut holder = spawn_child("hold", &dir, 1);
    wait_until("the child to hold the lock", || {
        dir.join("holding").exists()
    });
    // What a holder killed between renames leaves, written while the child
    // holds the lock: a table and hash without their key (never committed),
    // and a temporary file.
    std::fs::write(&paths.table, "not a table").unwrap();
    std::fs::write(&paths.hash, format!("{}\n", sha256_hex(b"not a table"))).unwrap();
    let stale = paths.table.with_file_name(format!(
        ".{}.tmp.4000000000.0",
        paths.table.file_name().unwrap().to_str().unwrap()
    ));
    std::fs::write(&stale, "partial").unwrap();
    holder.kill().unwrap();
    let status = holder.wait().unwrap();
    assert!(!status.success());

    // Bounded: a lock the dead child still held would time out here.
    let (tx, rx) = mpsc::channel();
    let (c, r2) = (cache.clone(), r.clone());
    std::thread::spawn(move || {
        let builds = AtomicUsize::new(0);
        let t = c.get_or_build(&r2, &r2.materials[0], TableKind::Elastic, || {
            builds.fetch_add(1, Ordering::SeqCst);
            synthetic_table(&r2, TableKind::Elastic)
        });
        let _ = tx.send((t.map_err(|e| format!("{e:#}")), builds.into_inner()));
    });
    let (t, builds) = rx
        .recv_timeout(DEADLINE)
        .expect("the next caller is not blocked by the killed holder");
    let t = t.unwrap();
    assert_eq!(builds, 1);
    assert_eq!(
        t.source,
        TableSource::Built,
        "the uncommitted table is rebuilt"
    );
    assert!(
        !stale.exists(),
        "the dead holder's temporary file is removed"
    );

    // The rebuilt entry is committed and read back.
    let again = cache
        .get_or_build(&r, m, TableKind::Elastic, || {
            panic!("a committed entry is not rebuilt")
        })
        .unwrap();
    assert_eq!(again.source, TableSource::Cache);
    assert_eq!(sha(&again), sha(&t));
    assert_eq!(again.table, t.table);
}

/// A committed entry that fails its checks is still refused (not taken for
/// a miss) with the lock in place, and the failed call releases the lock.
#[test]
fn corrupted_entry_is_refused_under_the_lock() {
    let dir = scratch("corrupted");
    let cache = cache_in(&dir);
    let r = example();
    let m = &r.materials[0];
    let t = cache
        .get_or_build(&r, m, TableKind::Inelastic, || {
            synthetic_table(&r, TableKind::Inelastic)
        })
        .unwrap();
    let path = t.cache.unwrap().path;
    let text = std::fs::read_to_string(&path).unwrap();
    let edited = text.replacen("1000000000", "2000000000", 1);
    assert_ne!(edited, text);
    std::fs::write(&path, edited).unwrap();
    for _ in 0..2 {
        let e = cache
            .get_or_build(&r, m, TableKind::Inelastic, || {
                panic!("a corrupted entry is not rebuilt")
            })
            .unwrap_err();
        assert!(format!("{e:#}").contains("SHA-256"), "{e:#}");
    }
}
