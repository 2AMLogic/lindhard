//! End-to-end tests of the `lindhard` binary: every example runs at small N
//! and writes parseable output; output is byte-identical across thread counts
//! (apart from the timing block); invalid input fails naming the field.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_lindhard");

fn examples_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../examples")
}

fn examples() -> Vec<PathBuf> {
    let mut v: Vec<_> = std::fs::read_dir(examples_dir())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "toml"))
        .collect();
    v.sort();
    assert!(v.len() >= 3, "expected the example inputs, found {v:?}");
    v
}

fn scratch(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("lindhard-cli-tests")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn lindhard(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().unwrap()
}

fn ok(o: &Output) {
    assert!(
        o.status.success(),
        "exit {:?}\nstdout:\n{}\nstderr:\n{}",
        o.status,
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
}

fn run(input: &Path, out: &Path, extra: &[&str]) {
    let mut args = vec![
        "run",
        input.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
    ];
    args.extend_from_slice(extra);
    ok(&lindhard(&args));
}

fn json(path: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

#[test]
fn version_includes_git_describe() {
    let o = lindhard(&["--version"]);
    ok(&o);
    let s = String::from_utf8(o.stdout).unwrap();
    assert!(
        s.starts_with(&format!("lindhard {} (", lindhard::VERSION)),
        "{s}"
    );
}

#[test]
fn every_example_checks_and_runs_at_small_n() {
    const N: u64 = 200;
    for ex in examples() {
        ok(&lindhard(&["check", ex.to_str().unwrap()]));
        let stem = ex.file_stem().unwrap().to_str().unwrap();
        let out = scratch(&format!("example-{stem}"));
        run(&ex, &out, &["--ions", &N.to_string()]);

        let s = json(&out.join("summary.json"));
        assert_eq!(s["format"]["name"], "lindhard-summary");
        assert_eq!(s["software"]["version"], lindhard::VERSION);
        assert_eq!(s["input"]["run"]["ions"], N, "{stem}: override echoed");
        assert!(s["input"]["run"].get("threads").is_none());
        // The echoed input is a complete input: it deserializes and resolves.
        let echo: lindhard::input::Input = serde_json::from_value(s["input"].clone()).unwrap();
        assert_eq!(echo.run.ions, N);
        echo.resolve().unwrap();
        assert!(!s["physics"]["models"].as_array().unwrap().is_empty());
        let r = &s["results"];
        assert_eq!(r["histories"], N);
        let p = &r["primaries"];
        let total = p["stopped"].as_u64().unwrap()
            + p["backscattered"].as_u64().unwrap()
            + p["transmitted"].as_u64().unwrap();
        assert_eq!(total, N, "{stem}");
        assert!(
            r["energy_budget_ev_per_ion"]["max_relative_residual"]
                .as_f64()
                .unwrap()
                < 1e-9
        );

        let depth = std::fs::read_to_string(out.join("depth_profile.csv")).unwrap();
        let bins = s["input"]["tally"]["depth_bins"].as_u64().unwrap() as usize;
        assert_eq!(depth.lines().count(), bins + 1, "{stem}: header + bins");
        let hist: u64 = depth
            .lines()
            .skip(1)
            .map(|l| l.split(',').nth(2).unwrap().parse::<u64>().unwrap())
            .sum();
        assert_eq!(hist, p["stopped"].as_u64().unwrap());

        let ions = std::fs::read_to_string(out.join("ions.csv")).unwrap();
        let rows: Vec<_> = ions.lines().skip(1).collect();
        assert_eq!(rows.len() as u64, N, "{stem}: one row per primary");
        for (i, row) in rows.iter().enumerate() {
            let f: Vec<_> = row.split(',').collect();
            assert_eq!(f.len(), 10);
            assert_eq!(f[0].parse::<u64>().unwrap(), i as u64, "sorted by index");
            for x in &f[2..9] {
                x.parse::<f64>().unwrap();
            }
        }
    }
}

/// Everything before the trailing `run` (timing) block.
fn deterministic_part(summary: &str) -> &str {
    let i = summary
        .find("\n  \"run\": {")
        .expect("run block is the last key");
    &summary[..i]
}

#[test]
fn output_is_byte_identical_across_thread_counts() {
    // The layered example exercises an interface, cascades and sputtering;
    // 300 ions are five chunks of 64, so threads really interleave.
    for name in ["as_50keV_si_sio2.toml", "ar_1keV_cu.toml"] {
        let ex = examples_dir().join(name);
        let mut outs = Vec::new();
        for threads in ["1", "4"] {
            let out = scratch(&format!("det-{name}-{threads}"));
            run(&ex, &out, &["--ions", "300", "--threads", threads]);
            outs.push(out);
        }
        let read = |d: &Path, f: &str| std::fs::read_to_string(d.join(f)).unwrap();
        let (a, b) = (
            read(&outs[0], "summary.json"),
            read(&outs[1], "summary.json"),
        );
        assert_eq!(deterministic_part(&a), deterministic_part(&b), "{name}");
        assert_eq!(json(&outs[0].join("summary.json"))["run"]["threads"], 1);
        assert_eq!(json(&outs[1].join("summary.json"))["run"]["threads"], 4);
        for f in ["depth_profile.csv", "ions.csv"] {
            assert_eq!(read(&outs[0], f), read(&outs[1], f), "{name}: {f}");
        }
    }
}

#[test]
fn seed_override_changes_results() {
    let ex = examples_dir().join("b_5keV_si.toml");
    let (a, b) = (scratch("seed-a"), scratch("seed-b"));
    run(&ex, &a, &["--ions", "100", "--seed", "1"]);
    run(&ex, &b, &["--ions", "100", "--seed", "2"]);
    let ions = |d: &Path| std::fs::read_to_string(d.join("ions.csv")).unwrap();
    assert_ne!(ions(&a), ions(&b));
    assert_eq!(json(&b.join("summary.json"))["input"]["run"]["seed"], 2);
}

/// Write `text` as an input file and return the failing command's stderr.
fn fails(name: &str, text: &str) -> String {
    let dir = scratch(&format!("bad-{name}"));
    let input = dir.join("input.toml");
    std::fs::write(&input, text).unwrap();
    for cmd in ["check", "run"] {
        let mut args = vec![cmd, input.to_str().unwrap()];
        if cmd == "run" {
            args.extend(["--out", dir.to_str().unwrap()]);
        }
        let o = lindhard(&args);
        assert!(!o.status.success(), "{cmd} {name} should fail");
    }
    assert!(
        !dir.join("summary.json").exists(),
        "{name}: no output on error"
    );
    String::from_utf8(lindhard(&["check", input.to_str().unwrap()]).stderr).unwrap()
}

const GOOD: &str = r#"
[beam]
ion = "B"
energy_ev = 5000.0

[target]
substrate = "Si"

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 10
seed = 1
"#;

#[test]
fn invalid_input_fails_naming_the_field() {
    let e = fails("unknown-key", &GOOD.replace("energy_ev", "energy_kev"));
    assert!(e.contains("energy_kev"), "{e}");

    let e = fails(
        "missing-ed",
        &GOOD.replace("[physics.energies.Si]\ne_d_ev = 15.0\n", ""),
    );
    assert!(
        e.contains("target.substrate") && e.contains("e_d_ev"),
        "{e}"
    );

    let e = fails(
        "negative-thickness",
        &GOOD.replace(
            "substrate = \"Si\"",
            "[[target.layers]]\nmaterial = \"Si\"\nthickness_nm = -5.0",
        ),
    );
    assert!(e.contains("target.layers[0].thickness_nm"), "{e}");

    let e = fails(
        "bad-model",
        &GOOD.replace("[physics]", "[physics]\nstopping = \"srim\""),
    );
    assert!(e.contains("srim"), "{e}");
}
