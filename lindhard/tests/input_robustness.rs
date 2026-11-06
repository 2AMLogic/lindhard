//! Malformed-input robustness: a deterministic, bounded mutation corpus over
//! the example inputs, run in child processes.
//!
//! Property: for any mutated input document, `from_toml_str` + `resolve_in`
//! returns (`Ok`, or an `Err` with a non-empty message); it never panics,
//! hangs or exhausts memory. Any other outcome is reported with a
//! classification (`Panic`, `Crash`, `Timeout`, `ResourceLimit`).
//!
//! Hostile cases are never resolved in the parent. The parent writes the cases
//! to a scratch directory and re-invokes this same test binary (the `worker`
//! test below) on ranges of them. A child that is silent for `CASE_TIMEOUT` is
//! killed; on unix the child's address space is capped with `ulimit -v`
//! (`MEM_CAP_KIB`), so an allocation blow-up aborts the child instead of the
//! host. The parent restarts a child on the remaining cases, so one bad case
//! does not hide later ones.
//!
//! Cases are enumerated in a fixed order (seed file, leaf in document order,
//! mutation), with no RNG, so a case id and index are stable and a failure
//! reproduces from its saved file with the command printed in the report.
//!
//! The memory cap needs `sh` and `ulimit`, so the whole file is `cfg(unix)`;
//! the repo has no `libc` dependency and denies `unsafe`.
#![cfg(unix)]

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::time::Duration;

use lindhard::input::electron::ElectronInput;
use lindhard::input::{Input, InputError};

/// Upper bound on the enumerated corpus. Exceeding it fails the test: the
/// corpus is never silently truncated or sampled.
const MAX_CASES: usize = 2_500;
/// Cases per child process before a fresh child is started.
const BATCH: usize = 250;
/// A child silent for this long is killed and its current case is a `Timeout`.
const CASE_TIMEOUT: Duration = Duration::from_secs(10);
/// Address-space cap for a child (`ulimit -v`), in KiB: 4 GiB.
const MEM_CAP_KIB: u64 = 4 * 1024 * 1024;
/// Known failing cases: (case id, issue number). Starts empty; may only shrink.
const KNOWN_FAILING: &[(&str, u32)] = &[];

const ENV_DIR: &str = "LINDHARD_ROBUSTNESS_DIR";
const ENV_RANGE: &str = "LINDHARD_ROBUSTNESS_RANGE";
const SYNTH_HANG: &str = "#!synthetic-hang";
const SYNTH_ALLOC: &str = "#!synthetic-alloc";

// ---------------------------------------------------------------- worker side

/// Resolve one document, as the CLI would. Returns the one-line outcome.
fn resolve_text(text: &str, base: &Path) -> String {
    let is_electron = text
        .parse::<toml::Table>()
        .map(|t| t.contains_key("electron"))
        .unwrap_or(false);
    let r: Result<(), InputError> = if is_electron {
        ElectronInput::from_toml_str(text).and_then(|i| i.resolve_in(base).map(|_| ()))
    } else {
        Input::from_toml_str(text).and_then(|i| i.resolve_in(base).map(|_| ()))
    };
    match r {
        Ok(()) => "OK".to_string(),
        Err(InputError::Parse(m)) => format!("ERR parse {}", one_line(&m)),
        Err(InputError::Invalid { field, message }) => {
            format!("ERR invalid {field}: {}", one_line(&message))
        }
    }
}

fn one_line(s: &str) -> String {
    s.replace(['\n', '\r'], " ")
}

/// Child entry point. A no-op unless the parent set the environment, so a
/// plain `cargo test` run passes it trivially.
#[test]
fn worker() {
    let Ok(dir) = std::env::var(ENV_DIR) else {
        return;
    };
    let dir = PathBuf::from(dir);
    let range = std::env::var(ENV_RANGE).expect("range");
    let (a, b) = range.split_once("..").expect("a..b");
    let (a, b): (usize, usize) = (a.parse().unwrap(), b.parse().unwrap());
    for idx in a..b {
        let text = std::fs::read_to_string(dir.join(format!("{idx}.toml"))).unwrap();
        let base = PathBuf::from(
            std::fs::read_to_string(dir.join(format!("{idx}.base")))
                .unwrap()
                .trim(),
        );
        println!("@@START {idx}");
        let outcome = if text.starts_with(SYNTH_HANG) {
            loop {
                std::thread::sleep(Duration::from_secs(3600));
            }
        } else if text.starts_with(SYNTH_ALLOC) {
            let v = vec![1u8; 16 << 30];
            format!("ERR synthetic {}", std::hint::black_box(v).len())
        } else {
            resolve_text(&text, &base)
        };
        println!("@@DONE {idx} {outcome}");
    }
}

// ---------------------------------------------------------------- parent side

#[derive(Debug, Clone, PartialEq)]
enum Class {
    ReturnedOk,
    ReturnedErr(String),
    Panic,
    Crash,
    Timeout,
    ResourceLimit,
}

impl Class {
    fn name(&self) -> &'static str {
        match self {
            Class::ReturnedOk => "ReturnedOk",
            Class::ReturnedErr(_) => "ReturnedErr",
            Class::Panic => "Panic",
            Class::Crash => "Crash",
            Class::Timeout => "Timeout",
            Class::ResourceLimit => "ResourceLimit",
        }
    }
    /// The property is "returns, does not panic/hang/blow up". A mutation may
    /// legitimately resolve (a deleted optional key, `seed = 0`), so `Ok`
    /// passes; an `Err` must carry a non-empty message.
    fn passes(&self) -> bool {
        match self {
            Class::ReturnedOk => true,
            Class::ReturnedErr(m) => {
                let m = m.trim();
                !(m.is_empty() || m == "parse" || m.ends_with(':'))
            }
            _ => false,
        }
    }
}

fn classify_dead_child(code: Option<i32>, stderr: &str) -> Class {
    if stderr.contains("memory allocation of") || stderr.contains("capacity overflow") {
        Class::ResourceLimit
    } else if code == Some(101) && stderr.contains("panicked at") {
        Class::Panic
    } else {
        Class::Crash
    }
}

/// Run cases `[start, end)` of the directory in children, restarting after a
/// timeout or death. Returns one class per case index.
fn run_range(dir: &Path, start: usize, end: usize, timeout: Duration) -> BTreeMap<usize, Class> {
    let exe = std::env::current_exe().unwrap();
    let mut done = BTreeMap::new();
    let mut cursor = start;
    while cursor < end {
        let mut child = Command::new("sh")
            .arg("-c")
            .arg("ulimit -v \"$1\"; shift; exec \"$@\"")
            .arg("sh")
            .arg(MEM_CAP_KIB.to_string())
            .arg(&exe)
            .args(["--exact", "worker", "--nocapture", "--test-threads=1"])
            .env(ENV_DIR, dir)
            .env(ENV_RANGE, format!("{cursor}..{end}"))
            .env("RAYON_NUM_THREADS", "2")
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("spawn worker");
        let mut stderr = child.stderr.take().unwrap();
        let err_thread = std::thread::spawn(move || {
            let mut s = String::new();
            let _ = stderr.read_to_string(&mut s);
            s
        });
        let stdout = child.stdout.take().unwrap();
        let (tx, rx) = mpsc::channel();
        let out_thread = std::thread::spawn(move || {
            for line in BufReader::new(stdout).lines().map_while(Result::ok) {
                if tx.send(line).is_err() {
                    break;
                }
            }
        });

        let mut open: Option<usize> = None;
        let mut timed_out = false;
        loop {
            match rx.recv_timeout(timeout) {
                Ok(line) => {
                    if let Some(rest) = line.split("@@START ").nth(1) {
                        open = rest.trim().parse().ok();
                    } else if let Some(rest) = line.split("@@DONE ").nth(1) {
                        let mut it = rest.splitn(3, ' ');
                        let idx: usize = it.next().unwrap().parse().unwrap();
                        let verdict = it.next().unwrap_or("");
                        let detail = it.next().unwrap_or("").to_string();
                        let class = if verdict == "OK" {
                            Class::ReturnedOk
                        } else {
                            Class::ReturnedErr(detail)
                        };
                        done.insert(idx, class);
                        open = None;
                    }
                }
                Err(mpsc::RecvTimeoutError::Timeout) => {
                    timed_out = true;
                    let _ = child.kill();
                    break;
                }
                Err(mpsc::RecvTimeoutError::Disconnected) => break,
            }
        }
        let status = child.wait().unwrap();
        let _ = out_thread.join();
        let stderr = err_thread.join().unwrap_or_default();

        // The case that was in flight when the child went silent or died.
        let next_unfinished = (cursor..end).find(|i| !done.contains_key(i));
        let Some(next) = next_unfinished else { break };
        let victim = open.unwrap_or(next);
        let class = if timed_out {
            Class::Timeout
        } else {
            classify_dead_child(status.code(), &stderr)
        };
        done.insert(victim, class);
        cursor = (victim + 1).max(cursor + 1);
        // Cases between `cursor` and `victim` were all reported before.
    }
    done
}

// ---------------------------------------------------------------- corpus

#[derive(Debug, Clone)]
struct Case {
    id: String,
    base: PathBuf,
    text: String,
}

#[derive(Debug, Clone)]
enum Seg {
    Key(String),
    Idx(usize),
}

fn path_string(p: &[Seg]) -> String {
    let mut s = String::new();
    for seg in p {
        match seg {
            Seg::Key(k) if s.is_empty() => s.push_str(k),
            Seg::Key(k) => {
                s.push('.');
                s.push_str(k);
            }
            Seg::Idx(i) => {
                let _ = write!(s, "[{i}]");
            }
        }
    }
    s
}

fn leaves(v: &toml::Value, path: &mut Vec<Seg>, out: &mut Vec<(Vec<Seg>, toml::Value)>) {
    match v {
        toml::Value::Table(t) => {
            for (k, c) in t {
                path.push(Seg::Key(k.clone()));
                leaves(c, path, out);
                path.pop();
            }
        }
        toml::Value::Array(a) => {
            for (i, c) in a.iter().enumerate() {
                path.push(Seg::Idx(i));
                leaves(c, path, out);
                path.pop();
            }
        }
        leaf => out.push((path.clone(), leaf.clone())),
    }
}

/// Apply `f` to the slot at `path` (`None` deletes it).
fn edit(root: &mut toml::Value, path: &[Seg], new: Option<toml::Value>) {
    let (last, parents) = path.split_last().unwrap();
    let mut cur = root;
    for seg in parents {
        cur = match seg {
            Seg::Key(k) => cur.as_table_mut().unwrap().get_mut(k).unwrap(),
            Seg::Idx(i) => &mut cur.as_array_mut().unwrap()[*i],
        };
    }
    match (last, new) {
        (Seg::Key(k), Some(v)) => {
            cur.as_table_mut().unwrap().insert(k.clone(), v);
        }
        (Seg::Key(k), None) => {
            cur.as_table_mut().unwrap().remove(k);
        }
        (Seg::Idx(i), Some(v)) => cur.as_array_mut().unwrap()[*i] = v,
        (Seg::Idx(i), None) => {
            cur.as_array_mut().unwrap().remove(*i);
        }
    }
}

fn seed_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let mut entries: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .collect();
    entries.sort();
    for p in entries {
        if p.is_dir() {
            seed_files(&p, out);
        } else if p.extension().is_some_and(|x| x == "toml") && is_input_document(&p) {
            out.push(p);
        }
    }
}

/// Input documents have a top-level `beam` (ion) or `electron` table; other
/// `.toml` files under `examples/` are data tables (e.g. a synthetic ELF).
fn is_input_document(p: &Path) -> bool {
    std::fs::read_to_string(p)
        .ok()
        .and_then(|t| t.parse::<toml::Table>().ok())
        .is_some_and(|t| t.contains_key("beam") || t.contains_key("electron"))
}

fn examples_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../examples")
}

fn seeds() -> Vec<PathBuf> {
    let mut v = Vec::new();
    seed_files(&examples_dir(), &mut v);
    assert!(v.len() >= 3, "expected the example inputs, found {v:?}");
    v
}

fn numeric_mutations(is_float: bool) -> Vec<(&'static str, toml::Value)> {
    use toml::Value::{Float, Integer};
    let mut m = vec![
        ("zero", if is_float { Float(0.0) } else { Integer(0) }),
        (
            "minus-one",
            if is_float { Float(-1.0) } else { Integer(-1) },
        ),
        ("tiny", Float(1e-300)),
        ("huge", Float(1e300)),
        ("nan", Float(f64::NAN)),
        ("inf", Float(f64::INFINITY)),
        ("neg-inf", Float(f64::NEG_INFINITY)),
        ("i64-max", Integer(i64::MAX)),
    ];
    if !is_float {
        // An integer leaf also gets the float forms; wrong-typed ones are
        // legitimate parse errors.
        m.retain(|(n, _)| *n != "tiny");
    }
    m
}

/// Enumerate the corpus in its fixed order.
fn corpus() -> Vec<Case> {
    let root = examples_dir();
    let mut cases = Vec::new();
    for seed in seeds() {
        let text = std::fs::read_to_string(&seed).unwrap();
        let doc: toml::Value = toml::from_str(&text).unwrap();
        let rel = seed
            .strip_prefix(&root)
            .unwrap()
            .to_string_lossy()
            .into_owned();
        let base = seed.parent().unwrap().to_path_buf();
        let mut ls = Vec::new();
        leaves(&doc, &mut Vec::new(), &mut ls);
        for (path, leaf) in ls {
            let p = path_string(&path);
            let mut muts: Vec<(String, Option<toml::Value>)> = Vec::new();
            let wrong = match &leaf {
                toml::Value::Integer(_) | toml::Value::Float(_) => {
                    for (n, v) in numeric_mutations(matches!(leaf, toml::Value::Float(_))) {
                        muts.push((n.to_string(), Some(v)));
                    }
                    toml::Value::String("x".into())
                }
                _ => toml::Value::Integer(7),
            };
            muts.push(("delete".into(), None));
            muts.push(("wrong-type".into(), Some(wrong)));
            for (name, newv) in muts {
                let mut d = doc.clone();
                edit(&mut d, &path, newv);
                cases.push(Case {
                    id: format!("{rel}|{p}|{name}"),
                    base: base.clone(),
                    text: toml::to_string(&d).unwrap(),
                });
            }
        }
    }
    cases
}

fn write_case(dir: &Path, idx: usize, text: &str, base: &Path) {
    std::fs::write(dir.join(format!("{idx}.toml")), text).unwrap();
    std::fs::write(
        dir.join(format!("{idx}.base")),
        base.to_string_lossy().as_bytes(),
    )
    .unwrap();
}

fn scratch_dir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("lindhard-robustness-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

// ---------------------------------------------------------------- tests

#[test]
fn seeds_resolve_ok() {
    for seed in seeds() {
        let text = std::fs::read_to_string(&seed).unwrap();
        let r = resolve_text(&text, seed.parent().unwrap());
        assert_eq!(r, "OK", "{}: {r}", seed.display());
    }
}

#[test]
fn corpus_is_bounded_and_deterministic() {
    let a = corpus();
    assert!(!a.is_empty());
    assert!(
        a.len() <= MAX_CASES,
        "corpus has {} cases, over MAX_CASES = {MAX_CASES}; raise it deliberately",
        a.len()
    );
    let b = corpus();
    let ids = |v: &[Case]| {
        v.iter()
            .map(|c| (c.id.clone(), c.text.clone()))
            .collect::<Vec<_>>()
    };
    assert_eq!(ids(&a), ids(&b));
    let mut seen: Vec<_> = a.iter().map(|c| &c.id).collect();
    seen.sort();
    seen.dedup();
    assert_eq!(seen.len(), a.len(), "case ids must be unique");
}

#[test]
fn containment_detects_timeout_and_resource_limit() {
    let dir = scratch_dir("synthetic");
    let base = Path::new(".");
    write_case(&dir, 0, SYNTH_HANG, base);
    write_case(&dir, 1, SYNTH_ALLOC, base);
    write_case(&dir, 2, "[nope]\n", base);
    let r = run_range(&dir, 0, 3, Duration::from_secs(3));
    assert_eq!(r[&0], Class::Timeout, "{r:?}");
    assert_eq!(r[&1], Class::ResourceLimit, "{r:?}");
    // The case after both failures still ran in a restarted child.
    assert!(matches!(r[&2], Class::ReturnedErr(_)), "{r:?}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn mutated_inputs_never_panic_hang_or_exhaust_memory() {
    let cases = corpus();
    assert!(cases.len() <= MAX_CASES, "{} > MAX_CASES", cases.len());
    let dir = scratch_dir("corpus");
    for (i, c) in cases.iter().enumerate() {
        write_case(&dir, i, &c.text, &c.base);
    }

    let started = std::time::Instant::now();
    let mut results = BTreeMap::new();
    let mut start = 0;
    while start < cases.len() {
        let end = (start + BATCH).min(cases.len());
        results.extend(run_range(&dir, start, end, CASE_TIMEOUT));
        start = end;
    }
    let ok = results
        .values()
        .filter(|c| **c == Class::ReturnedOk)
        .count();
    eprintln!(
        "input_robustness: {} cases ({ok} resolved Ok, {} Err) in {:.1?} \
         (BATCH={BATCH}, MEM_CAP_KIB={MEM_CAP_KIB})",
        cases.len(),
        cases.len() - ok,
        started.elapsed()
    );

    let mut report = String::new();
    let mut failures = 0;
    for (i, c) in cases.iter().enumerate() {
        let class = results.get(&i).cloned().unwrap_or(Class::Crash);
        if class.passes() {
            continue;
        }
        if KNOWN_FAILING.iter().any(|(id, _)| *id == c.id) {
            continue;
        }
        failures += 1;
        let _ = writeln!(
            report,
            "{} {} (batch {}, case {i})\n  file: {}\n  repro: {ENV_DIR}={} {ENV_RANGE}={i}..{} \
             cargo test -p lindhard --test input_robustness -- --exact worker --nocapture",
            class.name(),
            c.id,
            i / BATCH,
            dir.join(format!("{i}.toml")).display(),
            dir.display(),
            i + 1,
        );
        if let Class::ReturnedErr(m) = &class {
            let _ = writeln!(report, "  message: {m:?}");
        }
    }
    if failures == 0 {
        let _ = std::fs::remove_dir_all(&dir);
    } else {
        panic!(
            "{failures} of {} mutated inputs failed (scratch kept at {}):\n{report}",
            cases.len(),
            dir.display()
        );
    }
}
