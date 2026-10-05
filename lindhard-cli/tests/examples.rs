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

        check_tally_outputs(&out, &s, stem);

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

/// The range, damage, sputtering and escape sections and their CSVs agree
/// with the older summary counters and with each other.
fn check_tally_outputs(out: &Path, s: &serde_json::Value, stem: &str) {
    let r = &s["results"];
    let p = &r["primaries"];
    let n = r["histories"].as_u64().unwrap();
    for f in ["lateral_profile", "damage_profile", "escape_spectra"] {
        let name = s["files"][f].as_str().unwrap();
        assert!(out.join(name).exists(), "{stem}: {name}");
    }

    // Range: the legacy depth keys are the same quantity.
    let range = &r["range"];
    assert_eq!(range["stopped"], p["stopped"], "{stem}");
    if !range["depth"].is_null() {
        assert_eq!(range["depth"]["mean_nm"], p["stopped_depth_mean_nm"]);
        assert_eq!(range["depth"]["std_dev_nm"], p["stopped_depth_std_nm"]);
    }
    let by_layer: u64 = range["layers"]
        .as_array()
        .unwrap()
        .iter()
        .map(|l| l["stopped"].as_u64().unwrap())
        .sum();
    assert_eq!(by_layer, p["stopped"].as_u64().unwrap(), "{stem}");
    // The fit is reported (as a result or an error) only when requested.
    let want_fit = s["input"]["tally"]["dual_pearson"] == true;
    assert_eq!(
        range.get("dual_pearson").is_some() || range.get("dual_pearson_error").is_some(),
        want_fit,
        "{stem}"
    );

    // Sputtering and escapes against the summary counters.
    let sputtered = r["recoils"]["sputtered"].as_u64().unwrap();
    let by_element: u64 = r["sputtering"]["by_element"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| e["count"].as_u64().unwrap())
        .sum();
    assert_eq!(by_element, sputtered, "{stem}");
    let y = r["sputtering"]["yield_per_ion"].as_f64().unwrap();
    assert!((y - sputtered as f64 / n as f64).abs() < 1e-12, "{stem}");
    let e = &r["escapes"];
    let beam = &e["species"][0];
    assert_eq!(beam["beam"], true);
    assert_eq!(beam["front"]["count"], p["backscattered"], "{stem}");
    assert_eq!(beam["back"]["count"], p["transmitted"], "{stem}");

    // Damage: the cascade identity and per-layer sums.
    let d = &r["damage"];
    let c = &d["cascade"];
    assert_eq!(
        c["vacancies"].as_u64().unwrap(),
        c["displacements"].as_u64().unwrap() - c["replacements"].as_u64().unwrap()
    );
    assert_eq!(c["displacements"], r["recoils"]["displaced"], "{stem}");
    let layer_disp: u64 = d["layers"]
        .as_array()
        .unwrap()
        .iter()
        .map(|l| l["cascade"]["displacements"].as_u64().unwrap())
        .sum();
    assert_eq!(layer_disp, c["displacements"].as_u64().unwrap());

    // CSV totals match.
    let csv = |f: &str| std::fs::read_to_string(out.join(s["files"][f].as_str().unwrap())).unwrap();
    let col_sum = |text: &str, col: usize| -> u64 {
        text.lines()
            .skip(1)
            .map(|l| l.split(',').nth(col).unwrap().parse::<u64>().unwrap())
            .sum()
    };
    let dmg = csv("damage_profile");
    assert_eq!(col_sum(&dmg, 2), c["vacancies"].as_u64().unwrap(), "{stem}");
    assert_eq!(col_sum(&dmg, 3), c["interstitials"].as_u64().unwrap());
    assert_eq!(col_sum(&dmg, 4), c["replacements"].as_u64().unwrap());
    let lat = csv("lateral_profile");
    for q in ["y", "z", "radial"] {
        let rows: Vec<_> = lat
            .lines()
            .filter(|l| l.starts_with(&format!("{q},")))
            .collect();
        let total: u64 = rows
            .iter()
            .map(|l| l.split(',').nth(3).unwrap().parse::<u64>().unwrap())
            .sum();
        assert_eq!(total, p["stopped"].as_u64().unwrap(), "{stem}: {q}");
    }
    let esc = csv("escape_spectra");
    for (face, want) in [
        ("front", &beam["front"]["count"]),
        ("back", &beam["back"]["count"]),
    ] {
        for spec in ["energy_ev", "polar_deg"] {
            let tag = format!(",true,{face},{spec},");
            let total: u64 = esc
                .lines()
                .filter(|l| l.contains(&tag))
                .map(|l| l.split(',').nth(7).unwrap().parse::<u64>().unwrap())
                .sum();
            assert_eq!(total, want.as_u64().unwrap(), "{stem}: {face} {spec}");
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
        for f in [
            "depth_profile.csv",
            "ions.csv",
            "lateral_profile.csv",
            "damage_profile.csv",
            "escape_spectra.csv",
        ] {
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

#[test]
fn dual_pearson_toggle_and_tally_settings_are_honoured() {
    let dir = scratch("tally-settings");
    let input = dir.join("input.toml");
    let text = format!(
        "{GOOD}\n[tally]\ndual_pearson = true\nlateral_bin_nm = 2.0\nlateral_bins = 10\n\
         escape_energy_max_ev = 1000.0\nescape_energy_bins = 5\nescape_polar_bins = 3\n"
    );
    std::fs::write(&input, text).unwrap();
    let out = dir.join("out");
    run(&input, &out, &["--ions", "400"]);
    let s = json(&out.join("summary.json"));
    let range = &s["results"]["range"];
    assert!(
        !range["dual_pearson"].is_null() || !range["dual_pearson_error"].is_null(),
        "the requested fit is reported, as a result or an error"
    );
    assert_eq!(s["input"]["tally"]["lateral_bins"], 10);
    let lat = std::fs::read_to_string(out.join("lateral_profile.csv")).unwrap();
    // y, z: 20 bins each; radial: 10; each plus underflow and overflow rows.
    assert_eq!(lat.lines().count(), 1 + 2 * 22 + 12);
    let esc = std::fs::read_to_string(out.join("escape_spectra.csv")).unwrap();
    let beam_front_energy = esc
        .lines()
        .filter(|l| l.contains(",true,front,energy_ev,"))
        .count();
    assert_eq!(beam_front_energy, 5 + 2);
}
