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

/// A table for boron in silicon sampled from lindhard's own Lindhard-Scharff
/// model (no external data), as TOML text.
fn b_in_si_ls_table() -> String {
    use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
    use lindhard::ion::stopping::{to_ev_1e15_cm2, ElectronicStopping, Ion};
    let (ls, ion) = (LindhardScharff::new(), Ion::new(5).unwrap());
    let (lo, hi, n) = (1.0f64, 1.0e4f64, 80);
    let grid: Vec<f64> = (0..=n)
        .map(|i| lo * (hi / lo).powf(f64::from(i) / f64::from(n)))
        .collect();
    let list = |f: &dyn Fn(f64) -> f64| {
        grid.iter()
            .map(|&e| format!("{:e}", f(e)))
            .collect::<Vec<_>>()
            .join(", ")
    };
    format!(
        "provenance = \"Lindhard-Scharff model of this repository sampled on a grid (test)\"\n\
         ion_z = 5\ntarget_z = 14\nenergy_ev = [{}]\nstopping_ev_1e15_cm2 = [{}]\n",
        list(&|e| e),
        list(&|e| to_ev_1e15_cm2(ls.stopping(&ion, 14, e).unwrap()))
    )
}

#[test]
fn user_stopping_table_runs_deterministically_and_is_recorded() {
    let dir = scratch("table-run");
    std::fs::create_dir_all(dir.join("tables")).unwrap();
    let table_text = b_in_si_ls_table();
    std::fs::write(dir.join("tables/b_in_si.toml"), &table_text).unwrap();
    let input = dir.join("input.toml");
    // The table path is relative: it resolves against the input file's
    // directory, not the current directory.
    std::fs::write(
        &input,
        format!("{GOOD}\n[stopping]\ntables = [\"tables/b_in_si.toml\"]\n"),
    )
    .unwrap();
    ok(&lindhard(&["check", input.to_str().unwrap()]));

    let mut outs = Vec::new();
    for threads in ["1", "4"] {
        let out = dir.join(format!("out-{threads}"));
        run(&input, &out, &["--ions", "300", "--threads", threads]);
        outs.push(out);
    }
    let read = |d: &Path, f: &str| std::fs::read_to_string(d.join(f)).unwrap();
    let (a, b) = (
        read(&outs[0], "summary.json"),
        read(&outs[1], "summary.json"),
    );
    assert_eq!(deterministic_part(&a), deterministic_part(&b));
    for f in ["depth_profile.csv", "ions.csv", "damage_profile.csv"] {
        assert_eq!(read(&outs[0], f), read(&outs[1], f), "{f}");
    }

    let s = json(&outs[0].join("summary.json"));
    assert_eq!(
        s["input"]["stopping"]["tables"][0], "tables/b_in_si.toml",
        "echoed input carries the table path"
    );
    let rec = &s["physics"]["stopping_tables"][0];
    assert_eq!(rec["path"], "tables/b_in_si.toml");
    assert!(rec["resolved_path"]
        .as_str()
        .unwrap()
        .ends_with("tables/b_in_si.toml"));
    assert!(rec["provenance"]
        .as_str()
        .unwrap()
        .contains("Lindhard-Scharff"));
    let sha = rec["sha256"].as_str().unwrap();
    assert_eq!(sha.len(), 64);
    // Recompute independently: `shasum` ships with macOS and most Linux
    // systems; skip the cross-check where it is missing.
    if let Ok(o) = Command::new("shasum")
        .args([
            "-a",
            "256",
            dir.join("tables/b_in_si.toml").to_str().unwrap(),
        ])
        .output()
    {
        if o.status.success() {
            assert!(String::from_utf8_lossy(&o.stdout).starts_with(sha));
        }
    }
    let models = s["physics"]["models"].as_array().unwrap();
    assert!(models.iter().any(|m| m["name"] == "user-table"
        && m["citation"].as_str().unwrap().contains("Lindhard-Scharff")));
}

#[test]
fn ls_equivalent_table_matches_the_builtin_model() {
    let dir = scratch("table-equiv");
    std::fs::write(dir.join("t.toml"), b_in_si_ls_table()).unwrap();
    let plain = dir.join("plain.toml");
    let tabled = dir.join("tabled.toml");
    std::fs::write(&plain, GOOD).unwrap();
    std::fs::write(
        &tabled,
        format!("{GOOD}\n[stopping]\ntables = [\"t.toml\"]\n"),
    )
    .unwrap();
    let mean = |input: &Path, name: &str| {
        let out = dir.join(name);
        run(input, &out, &["--ions", "2000"]);
        json(&out.join("summary.json"))["results"]["range"]["depth"]["mean_nm"]
            .as_f64()
            .unwrap()
    };
    let (a, b) = (mean(&plain, "o-plain"), mean(&tabled, "o-tabled"));
    // Same seed and the same stopping up to interpolation: nearly identical.
    assert!((a / b - 1.0).abs() < 0.01, "plain {a} nm, table {b} nm");
}

#[test]
fn stopping_table_errors_fail_naming_the_field() {
    let e = fails(
        "table-missing",
        &format!("{GOOD}\n[stopping]\ntables = [\"no/such.toml\"]\n"),
    );
    assert!(
        e.contains("stopping.tables[0]") && e.contains("no/such.toml"),
        "{e}"
    );
    let e = fails(
        "table-unknown-key",
        &format!("{GOOD}\n[stopping]\ntable = \"x\"\n"),
    );
    assert!(e.contains("table"), "{e}");
}

/// A Si->Si table sampled from our own Lindhard-Scharff model on [lo, 1e4] eV.
fn si_in_si_ls_table(lo: f64) -> String {
    use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
    use lindhard::ion::stopping::{to_ev_1e15_cm2, ElectronicStopping, Ion};
    let (ls, ion) = (LindhardScharff::new(), Ion::new(14).unwrap());
    let (hi, n) = (1.0e4f64, 60);
    let grid: Vec<f64> = (0..=n)
        .map(|i| lo * (hi / lo).powf(f64::from(i) / f64::from(n)))
        .collect();
    let list = |f: &dyn Fn(f64) -> f64| {
        grid.iter()
            .map(|&e| format!("{:e}", f(e)))
            .collect::<Vec<_>>()
            .join(", ")
    };
    format!(
        "provenance = \"Lindhard-Scharff model of this repository sampled on a grid (test)\"\n\
         ion_z = 14\ntarget_z = 14\nenergy_ev = [{}]\nstopping_ev_1e15_cm2 = [{}]\n",
        list(&|e| e),
        list(&|e| to_ev_1e15_cm2(ls.stopping(&ion, 14, e).unwrap()))
    )
}

#[test]
fn recoil_species_tables_are_checked_before_the_run() {
    let dir = scratch("table-recoil");
    let input = dir.join("input.toml");
    std::fs::write(
        &input,
        format!("{GOOD}\n[stopping]\ntables = [\"si.toml\"]\n"),
    )
    .unwrap();

    // Starts above recoil_cutoff_ev: `check` now fails, naming the field.
    std::fs::write(dir.join("si.toml"), si_in_si_ls_table(100.0)).unwrap();
    let o = lindhard(&["check", input.to_str().unwrap()]);
    assert!(!o.status.success());
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(
        e.contains("stopping.tables[0]") && e.contains("recoil_cutoff_ev"),
        "{e}"
    );

    // Down to the cutoff: `check` lists the table by path, and a run in
    // which the table serves every Si recoil completes.
    std::fs::write(dir.join("si.toml"), si_in_si_ls_table(2.0)).unwrap();
    let o = lindhard(&["check", input.to_str().unwrap()]);
    ok(&o);
    let s = String::from_utf8_lossy(&o.stdout);
    assert!(s.contains("user-table (si.toml: Lindhard-Scharff"), "{s}");
    run(&input, &dir.join("out"), &["--ions", "50"]);
}

// ---- [dynamic] ---------------------------------------------------------------

fn dynamic_example() -> PathBuf {
    examples_dir().join("dynamic/as_1keV_si_film.toml")
}

#[test]
fn dynamic_example_checks_runs_and_writes_a_time_series() {
    let ex = dynamic_example();
    ok(&lindhard(&["check", ex.to_str().unwrap()]));
    let out = scratch("dynamic-example");
    run(&ex, &out, &["--ions", "600"]);

    let s = json(&out.join("dynamic_summary.json"));
    assert_eq!(s["format"]["name"], "lindhard-dynamic-summary");
    assert_eq!(s["input"]["run"]["ions"], 600);
    assert!(s["input"]["dynamic"]["fluence_cm2"].as_f64().unwrap() > 0.0);
    let echo: lindhard::input::Input = serde_json::from_value(s["input"].clone()).unwrap();
    echo.resolve().unwrap();
    assert!(!out.join("summary.json").exists(), "no static outputs");

    let steps = std::fs::read_to_string(out.join("dynamic_steps.csv")).unwrap();
    let header: Vec<_> = steps.lines().next().unwrap().split(',').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
    let rows: Vec<Vec<&str>> = steps
        .lines()
        .skip(1)
        .map(|l| l.split(',').collect())
        .collect();
    assert_eq!(rows[0][col("step")], "0");
    assert_eq!(
        rows.len() as u64 - 1,
        s["totals"]["steps"].as_u64().unwrap()
    );
    // Accepted steps tile the global primary indices 0..ions.
    let mut next = 0u64;
    for r in &rows[1..] {
        assert_eq!(r[col("first_index")].parse::<u64>().unwrap(), next);
        next += r[col("ions")].parse::<u64>().unwrap();
        assert_eq!(r[col("ions_done")].parse::<u64>().unwrap(), next);
        assert_eq!(r[col("surface_nm")], "0.0");
    }
    assert_eq!(next, 600);
    // Cumulative sputtering by element adds up to the total.
    let last = rows.last().unwrap();
    let by_el: u64 = header
        .iter()
        .enumerate()
        .filter(|(_, h)| h.starts_with("cum_sputtered_"))
        .map(|(i, _)| last[i].parse::<u64>().unwrap())
        .sum();
    assert_eq!(by_el, last[col("cum_sputtered")].parse::<u64>().unwrap());

    let comp = std::fs::read_to_string(out.join("dynamic_composition.csv")).unwrap();
    let n_steps = rows.len();
    assert!(comp.lines().count() > n_steps, "at least one slab per step");
    let h: Vec<_> = comp.lines().next().unwrap().split(',').collect();
    assert!(h.contains(&"fraction_As") && h.contains(&"atoms_per_cm2_Si"));
}

#[test]
fn dynamic_output_is_byte_identical_across_thread_counts() {
    let ex = dynamic_example();
    let read = |d: &Path, f: &str| std::fs::read_to_string(d.join(f)).unwrap();
    let mut outs = Vec::new();
    for threads in ["1", "2", "8"] {
        let out = scratch(&format!("dynamic-det-{threads}"));
        run(&ex, &out, &["--ions", "500", "--threads", threads]);
        outs.push(out);
    }
    for o in &outs[1..] {
        for f in ["dynamic_steps.csv", "dynamic_composition.csv"] {
            assert_eq!(read(&outs[0], f), read(o, f), "{f}");
        }
        assert_eq!(
            deterministic_part(&read(&outs[0], "dynamic_summary.json")),
            deterministic_part(&read(o, "dynamic_summary.json"))
        );
    }
}

#[test]
fn static_runs_are_unchanged_without_a_dynamic_section() {
    let out = scratch("static-no-dynamic");
    run(
        &examples_dir().join("b_5keV_si.toml"),
        &out,
        &["--ions", "100"],
    );
    let s = json(&out.join("summary.json"));
    assert!(s["input"].get("dynamic").is_none());
    assert!(!out.join("dynamic_steps.csv").exists());
}

#[test]
fn invalid_dynamic_settings_fail_naming_the_field() {
    let base = format!("{GOOD}\n[dynamic]\nfluence_cm2 = 1.0e15\nions_per_step = 5\n");
    let e = fails("dyn-fluence", &base.replace("1.0e15", "-1.0"));
    assert!(e.contains("dynamic.fluence_cm2"), "{e}");
    let e = fails("dyn-key", &base.replace("ions_per_step", "ions_each"));
    assert!(e.contains("ions_each"), "{e}");
    let e = fails("dyn-min", &format!("{base}min_ions_per_step = 2\n"));
    assert!(e.contains("dynamic.min_ions_per_step"), "{e}");
    let e = fails(
        "dyn-density",
        &format!("{base}relaxation = \"fixed-number-density\"\n"),
    );
    assert!(e.contains("dynamic.number_density_cm3"), "{e}");
    let e = fails("dyn-erosion", &format!("{base}erosion = 1\n"));
    assert!(e.contains("erosion"), "{e}");
}

// ---- electron runs -------------------------------------------------------------

/// The committed inputs of the level-3 backscatter comparison (#148,
/// `validation/experiments/backscatter/`) stay valid: `check` resolves them,
/// including the measured optical ELFs they name (no run: their tables take
/// minutes to build).
#[test]
fn backscatter_validation_inputs_check() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../validation/experiments/backscatter");
    let mut inputs: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "toml"))
        .collect();
    inputs.sort();
    assert_eq!(
        inputs.len(),
        5,
        "expected eta_{{al,au,c,cu,si}}.toml, found {inputs:?}"
    );
    for input in &inputs {
        let o = lindhard(&["check", input.to_str().unwrap()]);
        ok(&o);
        let stdout = String::from_utf8_lossy(&o.stdout);
        assert!(stdout.contains("electron run"), "{input:?}: {stdout}");
        assert!(stdout.contains("salvat-2003"), "{input:?}: corrections on");
        assert!(
            stdout.contains("elastic potential: salvat-dhfs"),
            "{input:?}: DHFS potential (#169)"
        );
        // Al, Au, C, Cu: Hagemann et al. (1975); Si: Yang et al. (2019).
        let si = input.file_name().is_some_and(|n| n == "eta_si.toml");
        let elf = if si { "Yang" } else { "Hagemann" };
        assert!(stdout.contains(elf), "{input:?}: measured ELF");
    }
}

fn electron_example() -> PathBuf {
    examples_dir().join("electron/e_10keV_si.toml")
}

/// Sum of column `col` of a CSV text's data rows whose first columns start
/// with `prefix`.
fn csv_sum(text: &str, prefix: &str, col: usize) -> f64 {
    text.lines()
        .skip(1)
        .filter(|l| l.starts_with(prefix))
        .map(|l| l.split(',').nth(col).unwrap().parse::<f64>().unwrap())
        .sum()
}

#[test]
fn electron_example_checks_runs_and_reports() {
    let ex = electron_example();
    let o = lindhard(&["check", ex.to_str().unwrap()]);
    ok(&o);
    assert!(String::from_utf8_lossy(&o.stdout).contains("electron run"));
    let out = scratch("electron-example");
    run(&ex, &out, &["--histories", "40"]);
    assert!(!out.join("summary.json").exists(), "no ion outputs");

    let s = json(&out.join("electron_summary.json"));
    assert_eq!(s["format"]["name"], "lindhard-electron-summary");
    assert_eq!(s["software"]["version"], lindhard::VERSION);
    assert_eq!(s["input"]["run"]["histories"], 40, "override echoed");
    assert!(s["input"]["run"].get("threads").is_none());
    // The echoed input is complete: it deserializes and resolves (data paths
    // against the example's directory).
    let echo: lindhard::input::electron::ElectronInput =
        serde_json::from_value(s["input"].clone()).unwrap();
    assert_eq!(echo.run.histories, 40);
    echo.resolve_in(&examples_dir().join("electron")).unwrap();
    assert!(
        echo.electron.tables.max_energy_ev.is_some(),
        "defaults filled"
    );

    // Every model choice and every data provenance is in the header.
    let p = &s["physics"];
    let models = p["models"].as_array().unwrap();
    for name in [
        "thomas-fermi-yukawa",
        "penn-single-pole",
        "kieft-bosch",
        "step-barrier",
    ] {
        assert!(models.iter().any(|m| m["name"] == name), "{name}");
    }
    let mat = &p["materials"][0];
    assert!(mat["optical_elf"]["provenance"]
        .as_str()
        .unwrap()
        .contains("synthetic"));
    assert_eq!(mat["optical_elf"]["sha256"].as_str().unwrap().len(), 64);
    assert!(mat["band"]["provenance"]
        .as_str()
        .unwrap()
        .contains("SYNTHETIC"));
    assert!(mat["elastic_table"]["model"]
        .as_str()
        .unwrap()
        .contains("STAND-IN"));
    assert_eq!(p["transport"]["n_histories"], 40);
    assert_eq!(p["transport"]["chunk_size"], 16);

    // Results: fates add up, the energy balance closes, yields agree.
    let r = &s["results"];
    assert_eq!(r["histories"], 40);
    let fates: u64 = r["fates"]
        .as_object()
        .unwrap()
        .values()
        .map(|v| v.as_u64().unwrap())
        .sum();
    assert_eq!(fates, 40);
    // Event-cap diagnostics: counts of tracks and histories, present even
    // when nothing was capped, and consistent with the primaries' caps.
    let caps = &r["event_caps"];
    let tracks = caps["secondary_tracks"].as_u64().unwrap();
    let affected = caps["affected_histories"].as_u64().unwrap();
    let capped = r["fates"]["event_capped"].as_u64().unwrap();
    assert!(affected >= capped && affected <= 40);
    assert!(affected <= capped + tracks);
    assert!(r["budget"]["relative_imbalance"].as_f64().unwrap() < 1e-9);
    let front = &r["front"];
    let eta = r["yields"]["backscatter_eta"].as_f64().unwrap();
    assert!((eta - front["fast"]["count"].as_f64().unwrap() / 40.0).abs() < 1e-12);

    // Table coverage: one entry per layer, the bounds of the tables the
    // layer used, and as many elastic as inelastic rate evaluations.
    let cov = r["table_coverage"].as_array().unwrap();
    assert_eq!(cov.len(), p["target"].as_array().unwrap().len());
    for (i, l) in cov.iter().enumerate() {
        assert_eq!(l["layer"], i);
        let total = |ch: &str| {
            ["below", "within", "above"]
                .iter()
                .map(|k| l[ch][k].as_u64().unwrap())
                .sum::<u64>()
        };
        assert!(total("elastic") > 0);
        assert_eq!(total("elastic"), total("inelastic"));
        for ch in ["elastic", "inelastic"] {
            let t = &mat[format!("{ch}_table")];
            assert_eq!(l[ch]["energy_min_ev"], t["energy_min_ev"], "{ch}");
            assert_eq!(l[ch]["energy_max_ev"], t["energy_max_ev"], "{ch}");
        }
    }

    // CSV files against the summary.
    let csv = |key: &str| {
        let name = s["files"][key].as_str().unwrap();
        std::fs::read_to_string(out.join(name)).unwrap()
    };
    let spectra = csv("escape_spectra");
    for face in ["front", "back"] {
        let n = r[face]["count"].as_f64().unwrap();
        assert_eq!(csv_sum(&spectra, &format!("{face},energy_ev,all,"), 5), n);
        for class in ["slow", "fast"] {
            assert_eq!(
                csv_sum(&spectra, &format!("{face},polar_deg,{class},"), 5),
                r[face][class]["count"].as_f64().unwrap(),
                "{face} {class}"
            );
        }
    }
    let cyl = csv("deposition_cylindrical");
    let inside = r["deposition"]["cylindrical"]["inside_ev"]
        .as_f64()
        .unwrap();
    assert!((csv_sum(&cyl, "", 6) / inside - 1.0).abs() < 1e-9);
    assert!(s["files"]["deposition_cartesian"].is_null());
    let tables = csv("tables");
    let n_grid = mat["inelastic_table"]["energies"].as_u64().unwrap() as usize;
    assert_eq!(tables.lines().count(), 1 + n_grid);
}

#[test]
fn electron_output_is_byte_identical_across_thread_counts() {
    // 80 histories are five chunks of 16, so the threads really interleave.
    let ex = electron_example();
    let mut outs = Vec::new();
    for threads in ["1", "4"] {
        let out = scratch(&format!("electron-det-{threads}"));
        run(&ex, &out, &["--histories", "80", "--threads", threads]);
        outs.push(out);
    }
    let read = |d: &Path, f: &str| std::fs::read_to_string(d.join(f)).unwrap();
    let (a, b) = (
        read(&outs[0], "electron_summary.json"),
        read(&outs[1], "electron_summary.json"),
    );
    assert_eq!(deterministic_part(&a), deterministic_part(&b));
    assert_eq!(
        json(&outs[0].join("electron_summary.json"))["run"]["threads"],
        1
    );
    assert_eq!(
        json(&outs[1].join("electron_summary.json"))["run"]["threads"],
        4
    );
    for f in [
        "electron_escape_spectra.csv",
        "electron_deposition_cylindrical.csv",
        "electron_tables.csv",
    ] {
        assert_eq!(read(&outs[0], f), read(&outs[1], f), "{f}");
    }
}

#[test]
fn electron_data_without_provenance_is_refused() {
    let dir = scratch("electron-noprov");
    let elf = std::fs::read_to_string(examples_dir().join("electron/synthetic_plasmon_elf.toml"))
        .unwrap();
    let stripped: String = elf
        .lines()
        .filter(|l| !l.starts_with("provenance"))
        .map(|l| format!("{l}\n"))
        .collect();
    std::fs::write(dir.join("synthetic_plasmon_elf.toml"), stripped).unwrap();
    let input = dir.join("input.toml");
    std::fs::copy(electron_example(), &input).unwrap();
    for cmd in ["check", "run"] {
        let mut args = vec![cmd, input.to_str().unwrap()];
        if cmd == "run" {
            args.extend(["--out", dir.to_str().unwrap()]);
        }
        let o = lindhard(&args);
        assert!(!o.status.success(), "{cmd} must refuse the ELF");
        let e = String::from_utf8_lossy(&o.stderr);
        assert!(
            e.contains("electron.materials.Si.optical_elf") && e.contains("provenance"),
            "{e}"
        );
    }
    assert!(!dir.join("electron_summary.json").exists());

    // Band parameters without a source are refused the same way.
    std::fs::write(dir.join("synthetic_plasmon_elf.toml"), elf).unwrap();
    let text = std::fs::read_to_string(electron_example()).unwrap();
    let start = text.find("provenance = \"SYNTHETIC").unwrap();
    let end = start + text[start..].find(" }").unwrap();
    std::fs::write(
        &input,
        format!("{}provenance = \"\"{}", &text[..start], &text[end..]),
    )
    .unwrap();
    let o = lindhard(&["check", input.to_str().unwrap()]);
    assert!(!o.status.success());
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(e.contains("electron.materials.Si.band"), "{e}");
}

#[test]
fn tuning_none_is_bit_identical_and_unknown_sets_are_rejected() {
    // `tuning = "none"` is the default: same physical output, same echo, and
    // no `physics.tuning` block.
    let text = std::fs::read_to_string(examples_dir().join("ar_1keV_cu.toml")).unwrap();
    let tuned = text.replace("[physics]", "[physics]\ntuning = \"none\"");
    assert_ne!(text, tuned, "example has a [physics] table");
    let mut summaries = Vec::new();
    for (name, body) in [("plain", &text), ("none", &tuned)] {
        let dir = scratch(&format!("tuning-{name}"));
        let input = dir.join("input.toml");
        std::fs::write(&input, body).unwrap();
        run(&input, &dir, &["--ions", "100", "--threads", "2"]);
        summaries.push(std::fs::read_to_string(dir.join("summary.json")).unwrap());
        assert!(json(&dir.join("summary.json"))["physics"]
            .get("tuning")
            .is_none());
    }
    assert_eq!(
        deterministic_part(&summaries[0]),
        deterministic_part(&summaries[1])
    );

    // A name the registry does not know is rejected.
    let e = fails(
        "unknown-tuning",
        &GOOD.replace("[physics]", "[physics]\ntuning = \"made-up\""),
    );
    assert!(e.contains("physics.tuning") && e.contains("made-up"), "{e}");
}

#[test]
fn shipped_tuning_set_is_deterministic_and_equals_an_explicit_override() {
    // The shipped pilot set on the Ar -> Cu example: identical physical output
    // on 1, 2 and 8 threads; the same results as an untuned run with the
    // effective E_s given explicitly (the set multiplies the resolved E_s once
    // and changes nothing else); and the metadata in summary.json.
    let set = lindhard::input::ES_SPUTTER_AR_V1;
    let text = std::fs::read_to_string(examples_dir().join("ar_1keV_cu.toml")).unwrap();
    let tuned = replaced(
        &text,
        "[physics]",
        &format!("[physics]\ntuning = \"{}\"", set.name),
    );
    let mut outs = Vec::new();
    for threads in ["1", "2", "8"] {
        let dir = scratch(&format!("tuning-set-{threads}"));
        let input = dir.join("input.toml");
        std::fs::write(&input, &tuned).unwrap();
        run(&input, &dir, &["--ions", "300", "--threads", threads]);
        outs.push(dir);
    }
    let read = |d: &Path, f: &str| std::fs::read_to_string(d.join(f)).unwrap();
    for o in &outs[1..] {
        assert_eq!(
            deterministic_part(&read(&outs[0], "summary.json")),
            deterministic_part(&read(o, "summary.json"))
        );
        for f in ["depth_profile.csv", "ions.csv", "escape_spectra.csv"] {
            assert_eq!(read(&outs[0], f), read(o, f), "{f}");
        }
    }
    let s = json(&outs[0].join("summary.json"));
    let t = &s["physics"]["tuning"];
    assert_eq!(t["set"], set.name);
    assert_eq!(t["version"], set.version);
    assert_eq!(t["quantity"], "surface-binding-energy");
    let c = &t["components"][0];
    let k = set.e_s_factors.iter().find(|p| p.0 == "Cu").unwrap().1;
    assert_eq!(c["element"], "Cu");
    assert_eq!(c["factor"].as_f64().unwrap(), k);
    let original = c["e_s_original_ev"].as_f64().unwrap();
    // summary.json rounds floats when written, so compare within that; the
    // explicit run below uses the unrounded product the engine applies.
    let effective = original * k;
    let reported = c["e_s_effective_ev"].as_f64().unwrap();
    assert!(
        (reported - effective).abs() < 1e-9 * effective,
        "{reported} vs {effective}"
    );
    assert_eq!(s["input"]["physics"]["tuning"], set.name);

    // The same run untuned, with E_s = the effective value set explicitly.
    let explicit = format!("{text}\n[physics.energies.Cu]\ne_s_ev = {effective:?}\n");
    let dir = scratch("tuning-set-explicit");
    let input = dir.join("input.toml");
    std::fs::write(&input, explicit).unwrap();
    run(&input, &dir, &["--ions", "300", "--threads", "2"]);
    let e = json(&dir.join("summary.json"));
    assert!(e["physics"].get("tuning").is_none());
    assert_eq!(s["results"], e["results"]);
    for f in ["depth_profile.csv", "ions.csv", "escape_spectra.csv"] {
        assert_eq!(read(&outs[0], f), read(&dir, f), "{f}");
    }

    // Other beams and unlisted elements are rejected.
    let e = fails(
        "tuning-other-ion",
        &replaced(&tuned, "ion = \"Ar\"", "ion = \"Xe\""),
    );
    assert!(
        e.contains("physics.tuning") && e.contains("beam ion"),
        "{e}"
    );
    let e = fails(
        "tuning-other-target",
        &replaced(&tuned, "substrate = \"Cu\"", "substrate = \"Ni\""),
    );
    assert!(
        e.contains("physics.tuning") && e.contains("no factor for Ni"),
        "{e}"
    );
}

/// Rewrites `text` with `from` replaced by `to`, asserting it was present.
fn replaced(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "{from:?} not in input");
    text.replacen(from, to, 1)
}

fn sentinel(out: &Path) -> PathBuf {
    let p = out.join("user_notes.txt");
    std::fs::write(&p, "keep").unwrap();
    p
}

#[test]
fn rerun_removes_stale_per_ion_csv_and_keeps_unrelated_files() {
    let dir = scratch("reuse-ions");
    let out = dir.join("out");
    let mk = |per_ion: bool| {
        let input = dir.join(format!("in-{per_ion}.toml"));
        std::fs::write(&input, format!("{GOOD}\n[tally]\nper_ion = {per_ion}\n")).unwrap();
        input
    };
    run(&mk(true), &out, &["--ions", "50"]);
    let note = sentinel(&out);
    let ions = out.join("ions.csv");
    assert!(ions.exists());
    run(&mk(false), &out, &["--ions", "50"]);
    assert!(!ions.exists(), "stale ions.csv survived");
    assert!(json(&out.join("summary.json"))["files"]["ions"].is_null());
    assert!(note.exists());
    // A missing optional file is harmless on a further disabled rerun.
    run(&mk(false), &out, &["--ions", "50"]);
    assert!(note.exists());
    run(&mk(true), &out, &["--ions", "50"]);
    let text = std::fs::read_to_string(&ions).unwrap();
    assert_eq!(text.lines().count(), 51);
    assert!(!json(&out.join("summary.json"))["files"]["ions"].is_null());
    assert!(note.exists());
}

#[test]
fn rerun_removes_stale_electron_deposition_csvs() {
    let dir = scratch("reuse-electron");
    let out = dir.join("out");
    let ex_dir = electron_example().parent().unwrap().to_path_buf();
    let src = std::fs::read_to_string(electron_example()).unwrap();
    let src = replaced(
        &src,
        "\"synthetic_plasmon_elf.toml\"",
        &format!("{:?}", ex_dir.join("synthetic_plasmon_elf.toml")),
    );
    let cyl = "[electron.tally.cylindrical]\nr = { lo_nm = 0.0, hi_nm = 1000.0, bins = 20 }\n\
               depth = { lo_nm = 0.0, hi_nm = 2000.0, bins = 40 }\n";
    let cart = "[electron.tally.cartesian]\nx = { lo_nm = 0.0, hi_nm = 2000.0, bins = 4 }\n\
                y = { lo_nm = -500.0, hi_nm = 500.0, bins = 2 }\n\
                z = { lo_nm = -500.0, hi_nm = 500.0, bins = 2 }\n";
    let without = replaced(&src, cyl, "");
    let input = |name: &str, with_cyl: bool, with_cart: bool| {
        let mut t = without.clone();
        if with_cyl {
            t = replaced(&t, "[target]", &format!("{cyl}\n[target]"));
        }
        if with_cart {
            t = replaced(&t, "[target]", &format!("{cart}\n[target]"));
        }
        let p = dir.join(format!("{name}.toml"));
        std::fs::write(&p, t).unwrap();
        p
    };
    let both = input("both", true, true);
    let only_cyl = input("cyl", true, false);
    let only_cart = input("cart", false, true);
    let none = input("none", false, false);
    let cy = out.join("electron_deposition_cylindrical.csv");
    let ca = out.join("electron_deposition_cartesian.csv");
    let summary = out.join("electron_summary.json");
    let args = ["--histories", "32"];
    let steps: [(&Path, bool, bool); 6] = [
        (&both, true, true),
        (&only_cyl, true, false),
        (&only_cart, false, true),
        (&none, false, false),
        (&none, false, false),
        (&both, true, true),
    ];
    let mut note = None;
    for (inp, want_cyl, want_cart) in steps {
        run(inp, &out, &args);
        let note = note.get_or_insert_with(|| sentinel(&out));
        assert_eq!(cy.exists(), want_cyl, "cylindrical for {inp:?}");
        assert_eq!(ca.exists(), want_cart, "cartesian for {inp:?}");
        let s = json(&summary);
        assert_eq!(!s["files"]["deposition_cylindrical"].is_null(), want_cyl);
        assert_eq!(!s["files"]["deposition_cartesian"].is_null(), want_cart);
        assert!(note.exists());
    }
}

#[test]
fn stale_file_that_cannot_be_removed_is_a_contextual_error() {
    let dir = scratch("reuse-undeletable");
    let out = dir.join("out");
    let input = dir.join("in.toml");
    std::fs::write(&input, format!("{GOOD}\n[tally]\nper_ion = false\n")).unwrap();
    // A directory under the reserved name cannot be removed as a file.
    std::fs::create_dir_all(out.join("ions.csv")).unwrap();
    let o = lindhard(&[
        "run",
        input.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--ions",
        "20",
    ]);
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("removing stale") && err.contains("ions.csv"),
        "{err}"
    );
}

/// Minimal RFC 4180 parser: records of fields, quoted fields may contain
/// commas, CR, LF and doubled quotes.
fn parse_csv(text: &str) -> Vec<Vec<String>> {
    let (mut rows, mut row, mut f) = (Vec::new(), Vec::new(), String::new());
    let (mut quoted, mut chars) = (false, text.chars().peekable());
    while let Some(c) = chars.next() {
        match (quoted, c) {
            (true, '"') if chars.peek() == Some(&'"') => {
                chars.next();
                f.push('"');
            }
            (true, '"') => quoted = false,
            (true, _) => f.push(c),
            (false, '"') => quoted = true,
            (false, ',') => row.push(std::mem::take(&mut f)),
            (false, '\n') => {
                row.push(std::mem::take(&mut f));
                rows.push(std::mem::take(&mut row));
            }
            (false, _) => f.push(c),
        }
    }
    assert!(!quoted, "unterminated quoted field");
    assert!(f.is_empty() && row.is_empty(), "missing final newline");
    rows
}

#[test]
fn electron_tables_csv_encodes_material_names() {
    let dir = scratch("electron-csv-names");
    let src = std::fs::read_to_string(electron_example()).unwrap();
    std::fs::copy(
        examples_dir().join("electron/synthetic_plasmon_elf.toml"),
        dir.join("synthetic_plasmon_elf.toml"),
    )
    .unwrap();
    // Two quoted TOML keys: a comma, and commas with embedded quotes. (Line
    // breaks are rejected by the material-identity check, so they cannot
    // reach the table; the encoder's unit test covers them.)
    let film = "Si,film";
    let sub = "Si \"q\", \"\"r\"\"";
    let mat = |name: &str| {
        let m = "[electron.materials.Si]";
        let body = src.split(m).nth(1).unwrap();
        let body = body.split("\n\n").next().unwrap();
        let key = name.replace('\\', "\\\\").replace('"', "\\\"");
        let key = key.replace('\r', "\\r").replace('\n', "\\n");
        format!("[electron.materials.\"{key}\"]{body}\n\n")
    };
    let key = |name: &str| {
        let k = name.replace('"', "\\\"").replace('\r', "\\r");
        k.replace('\n', "\\n")
    };
    let head = src.split("[electron.materials.Si]").next().unwrap();
    let tail = src.split("[electron.tally]").nth(1).unwrap();
    let tail = tail.split("[target]").next().unwrap();
    let text = format!(
        "{head}{}{}[electron.tally]{tail}\n\
         [materials.\"{}\"]\ndensity_g_cm3 = 2.33\nelements = [{{ symbol = \"Si\", atom_fraction = 1.0 }}]\n\n\
         [materials.\"{}\"]\ndensity_g_cm3 = 2.33\nelements = [{{ symbol = \"Si\", atom_fraction = 1.0 }}]\n\n\
         [target]\nsubstrate = \"{}\"\n\n[[target.layers]]\nmaterial = \"{}\"\nthickness_nm = 20.0\n\n\
         [run]\nhistories = 10\nseed = 1\n",
        mat(film),
        mat(sub),
        key(film),
        key(sub),
        key(sub),
        key(film),
    );
    let input = dir.join("in.toml");
    std::fs::write(&input, text).unwrap();
    let out = dir.join("out");
    run(&input, &out, &[]);
    let csv = std::fs::read_to_string(out.join("electron_tables.csv")).unwrap();
    let rows = parse_csv(&csv);
    assert_eq!(rows[0].len(), 6);
    assert_eq!(rows[0][0], "material");
    let mut seen = std::collections::BTreeSet::new();
    for r in &rows[1..] {
        assert_eq!(r.len(), 6, "{r:?}");
        seen.insert(r[0].clone());
        assert!(r[1].parse::<f64>().unwrap() > 0.0);
        for c in [3, 4, 5] {
            assert!(r[c].parse::<f64>().unwrap().is_finite(), "{r:?}");
        }
        assert!(r[2].is_empty() || r[2].parse::<f64>().unwrap().is_finite());
    }
    let want: std::collections::BTreeSet<String> =
        [film, sub].iter().map(|s| s.to_string()).collect();
    assert_eq!(seen, want);
    // Both materials share the Si data, so their numeric columns agree.
    let n = (rows.len() - 1) / 2;
    for i in 1..=n {
        assert_eq!(rows[i][1..], rows[i + n][1..]);
    }
}

// ---- cross-section table cache (`--table-cache`) ---------------------------

const ELECTRON_FILES: [&str; 3] = [
    "electron_escape_spectra.csv",
    "electron_deposition_cylindrical.csv",
    "electron_tables.csv",
];

/// The summary without the parts a table cache is allowed to change: the
/// `run` block (timings, threads) and each table's `source` and `cache`.
fn summary_without_cache_echo(dir: &Path) -> serde_json::Value {
    let mut s = json(&dir.join("electron_summary.json"));
    s.as_object_mut().unwrap().remove("run");
    for m in s["physics"]["materials"].as_array_mut().unwrap() {
        for t in ["elastic_table", "inelastic_table"] {
            let o = m[t].as_object_mut().unwrap();
            assert!(o.remove("source").is_some(), "{t}.source echoed");
            assert!(o.remove("cache").is_some(), "{t}.cache echoed");
        }
    }
    s
}

fn table_sources(dir: &Path) -> (String, String) {
    let s = json(&dir.join("electron_summary.json"));
    let m = &s["physics"]["materials"][0];
    (
        m["elastic_table"]["source"].as_str().unwrap().to_string(),
        m["inelastic_table"]["source"].as_str().unwrap().to_string(),
    )
}

fn cache_files(cache: &Path) -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_dir(cache)
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect();
    v.sort();
    v
}

/// A run that reads its tables from the cache gives the same results, bit
/// for bit, as one that builds them, at any thread count.
#[test]
fn electron_table_cache_reuse_is_bit_identical() {
    let ex = electron_example();
    let dir = scratch("electron-cache-reuse");
    let cache = dir.join("cache");
    let c = cache.to_str().unwrap();
    let out = |name: &str| dir.join(name);

    run(&ex, &out("plain"), &["--histories", "80", "--threads", "4"]);
    run(
        &ex,
        &out("miss"),
        &["--histories", "80", "--threads", "4", "--table-cache", c],
    );
    let files = cache_files(&cache);
    assert_eq!(
        files.len(),
        6,
        "two tables, two hashes and two keys: {files:?}"
    );
    assert!(files
        .iter()
        .any(|f| f.starts_with("elastic-") && f.ends_with(".toml")));
    assert!(files
        .iter()
        .any(|f| f.starts_with("inelastic-") && f.ends_with(".key.json")));
    assert!(files
        .iter()
        .any(|f| f.starts_with("elastic-") && f.ends_with(".sha256")));
    for threads in ["1", "8"] {
        let o = lindhard(&[
            "run",
            ex.to_str().unwrap(),
            "--out",
            out(&format!("hit-{threads}")).to_str().unwrap(),
            "--histories",
            "80",
            "--threads",
            threads,
            "--table-cache",
            c,
        ]);
        ok(&o);
        assert!(
            String::from_utf8_lossy(&o.stderr).contains("tables: 2 read from the cache, 0 built"),
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
    }
    assert_eq!(cache_files(&cache), files, "a hit writes nothing");

    assert_eq!(
        table_sources(&out("plain")),
        ("built".into(), "built".into())
    );
    assert_eq!(
        table_sources(&out("miss")),
        ("built".into(), "built".into())
    );
    assert_eq!(
        table_sources(&out("hit-1")),
        ("cache".into(), "cache".into())
    );
    let plain = json(&out("plain").join("electron_summary.json"));
    assert!(plain["physics"]["materials"][0]["inelastic_table"]["cache"].is_null());
    let miss = json(&out("miss").join("electron_summary.json"));
    let hit = json(&out("hit-1").join("electron_summary.json"));
    for t in ["elastic_table", "inelastic_table"] {
        let (a, b) = (
            &miss["physics"]["materials"][0][t]["cache"],
            &hit["physics"]["materials"][0][t]["cache"],
        );
        assert_eq!(a, b, "{t}: same file, same SHA-256");
        assert_eq!(a["sha256"].as_str().unwrap().len(), 64);
        let path = PathBuf::from(a["path"].as_str().unwrap());
        assert!(path.starts_with(std::fs::canonicalize(&cache).unwrap()));
        let key = a["key_sha256"].as_str().unwrap();
        assert!(path.file_name().unwrap().to_str().unwrap().contains(key));
    }

    let want = summary_without_cache_echo(&out("plain"));
    for name in ["miss", "hit-1", "hit-8"] {
        assert_eq!(summary_without_cache_echo(&out(name)), want, "{name}");
        for f in ELECTRON_FILES {
            assert_eq!(
                std::fs::read(out("plain").join(f)).unwrap(),
                std::fs::read(out(name).join(f)).unwrap(),
                "{name}: {f}"
            );
        }
    }
    // Two cached runs differ only in the run block.
    let read = |d: &str| std::fs::read_to_string(out(d).join("electron_summary.json")).unwrap();
    assert_eq!(
        deterministic_part(&read("hit-1")),
        deterministic_part(&read("hit-8"))
    );
}

/// Writes `text` as the cached table at `path` and records its SHA-256 in the
/// entry's hash file, as a store would: an edit the hash check cannot see,
/// to reach the checks behind it.
fn write_table_and_hash(path: &Path, text: &str) {
    std::fs::write(path, text).unwrap();
    std::fs::write(
        path.with_extension("sha256"),
        format!(
            "{}\n",
            lindhard_cli::table_cache::sha256_hex(text.as_bytes())
        ),
    )
    .unwrap();
}

/// A run whose physics differs is not served another run's table, and a
/// cache file that does not match its key or its recorded SHA-256 is refused
/// with the field or file named.
#[test]
fn electron_table_cache_misses_on_changed_physics_and_refuses_bad_files() {
    let dir = scratch("electron-cache-keys");
    let ex_dir = electron_example().parent().unwrap().to_path_buf();
    let src = std::fs::read_to_string(electron_example())
        .unwrap()
        .replace(
            "optical_elf = \"synthetic_plasmon_elf.toml\"",
            &format!(
                "optical_elf = {:?}",
                ex_dir.join("synthetic_plasmon_elf.toml").to_str().unwrap()
            ),
        );
    let input = dir.join("in.toml");
    let cache = dir.join("cache");
    let c = cache.to_str().unwrap();
    let go = |text: &str, name: &str| -> Output {
        std::fs::write(&input, text).unwrap();
        lindhard(&[
            "run",
            input.to_str().unwrap(),
            "--out",
            dir.join(name).to_str().unwrap(),
            "--histories",
            "16",
            "--table-cache",
            c,
        ])
    };
    ok(&go(&src, "a"));
    assert_eq!(cache_files(&cache).len(), 6);

    // Another band with the same inner potential (so the same table grid):
    // the inelastic table, built on the band-bottom axis with the band's
    // minimum excitation energy (#241), is built and stored beside the
    // first; the elastic table is reused.
    let fermi = src.replace(
        "valence_band_width_ev = 10.0, band_gap_ev = 2.0, affinity_ev = 3.0",
        "valence_band_width_ev = 10.5, band_gap_ev = 2.0, affinity_ev = 2.5",
    );
    assert_ne!(fermi, src);
    ok(&go(&fermi, "b"));
    assert_eq!(
        table_sources(&dir.join("b")),
        ("cache".into(), "built".into())
    );
    assert_eq!(cache_files(&cache).len(), 9);

    // A key file edited under its hash name is refused, naming the field.
    let key = cache_files(&cache)
        .into_iter()
        .find(|f| f.starts_with("elastic-") && f.ends_with(".key.json"))
        .unwrap();
    let key_path = cache.join(&key);
    let good = std::fs::read_to_string(&key_path).unwrap();
    let bad = good.replace("\"kind\":\"elastic\"", "\"kind\":\"inelastic\"");
    assert_ne!(bad, good);
    std::fs::write(&key_path, bad).unwrap();
    let o = go(&src, "c");
    assert!(!o.status.success());
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(e.contains("field `kind` differs"), "{e}");
    std::fs::write(&key_path, &good).unwrap();

    let table_path = cache.join(key.replace(".key.json", ".toml"));
    let hash_path = cache.join(key.replace(".key.json", ".sha256"));
    let table = std::fs::read_to_string(&table_path).unwrap();
    let hash = std::fs::read_to_string(&hash_path).unwrap();

    // One cross-section value edited so the table still parses and
    // validates: refused by its SHA-256, naming the file.
    let first = table.find("inverse_mfp_per_m = [").unwrap() + "inverse_mfp_per_m = [".len();
    let end = first + table[first..].find(',').unwrap();
    let value: f64 = table[first..end].trim().parse().unwrap();
    assert!(value > 0.0);
    let edited = format!("{}{:?}{}", &table[..first], value * 1.5, &table[end..]);
    lindhard::electron::data::CrossSectionTable::from_toml_str(&edited)
        .expect("the edited table is still a valid table");
    std::fs::write(&table_path, &edited).unwrap();
    let o = go(&src, "g");
    assert!(!o.status.success());
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(e.contains("SHA-256"), "{e}");
    assert!(
        e.contains(table_path.file_name().unwrap().to_str().unwrap()),
        "{e}"
    );
    std::fs::write(&table_path, &table).unwrap();

    // A missing hash file is an error naming it, not a silent reuse.
    std::fs::remove_file(&hash_path).unwrap();
    let o = go(&src, "h");
    assert!(!o.status.success());
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(
        e.contains(hash_path.file_name().unwrap().to_str().unwrap()),
        "{e}"
    );
    assert!(e.contains("missing"), "{e}");
    std::fs::write(&hash_path, &hash).unwrap();

    // A table of another cache format version under the key is an error,
    // not a silent rebuild (its hash rewritten to reach the loader).
    let line = format!(
        "format_version = {}",
        lindhard::electron::data::CACHE_FORMAT_VERSION
    );
    assert!(table.contains(&line));
    write_table_and_hash(&table_path, &table.replace(&line, "format_version = 0"));
    let o = go(&src, "d");
    assert!(!o.status.success());
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(e.contains("version"), "{e}");

    // A table whose grid was edited (and its hash with it) is refused.
    let first = table.find("energy_ev = [").unwrap() + "energy_ev = [".len();
    let end = first + table[first..].find(',').unwrap();
    let edited = format!("{}10.5{}", &table[..first], &table[end..]);
    write_table_and_hash(&table_path, &edited);
    let o = go(&src, "e");
    assert!(!o.status.success());
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(e.contains("field `energy_ev`"), "{e}");

    // Restored, it is read again.
    write_table_and_hash(&table_path, &table);
    assert_eq!(std::fs::read_to_string(&hash_path).unwrap(), hash);
    ok(&go(&src, "f"));
    assert_eq!(
        table_sources(&dir.join("f")),
        ("cache".into(), "cache".into())
    );
}

#[test]
fn table_cache_is_refused_for_ion_runs() {
    let dir = scratch("ion-table-cache");
    let o = lindhard(&[
        "run",
        examples_dir().join("ar_1keV_cu.toml").to_str().unwrap(),
        "--out",
        dir.to_str().unwrap(),
        "--table-cache",
        dir.join("cache").to_str().unwrap(),
    ]);
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("electron runs only"));
}

/// The electron example with a `[electron.tally.psf]` table appended and its
/// data files copied beside the input.
fn electron_psf_input(name: &str, psf: &str) -> (PathBuf, PathBuf) {
    let dir = scratch(name);
    let src = examples_dir().join("electron");
    for f in std::fs::read_dir(&src).unwrap() {
        let f = f.unwrap().path();
        if f.is_file() && f.extension().is_some_and(|e| e == "toml") {
            std::fs::copy(&f, dir.join(f.file_name().unwrap())).unwrap();
        }
    }
    let text = std::fs::read_to_string(electron_example()).unwrap();
    let input = dir.join("input.toml");
    std::fs::write(&input, format!("{text}\n[electron.tally.psf]\n{psf}\n")).unwrap();
    (dir, input)
}

const PSF_SLAB: &str =
    "depth_lo_nm = 0.0\ndepth_hi_nm = 300.0\nr_min_nm = 5.0\nr_max_nm = 2000.0\nbins = 12";

#[test]
fn electron_psf_writes_profile_parameters_and_summary() {
    let (dir, input) = electron_psf_input("electron-psf", PSF_SLAB);
    let out = dir.join("out");
    run(&input, &out, &["--histories", "160"]);
    let profile = std::fs::read_to_string(out.join("electron_psf_profile.csv")).unwrap();
    let params = std::fs::read_to_string(out.join("electron_psf_parameters.csv")).unwrap();
    assert!(params.starts_with("model,parameter,value,std_error\n"));
    // Central disc plus 12 log bins, and a header.
    assert_eq!(profile.lines().count(), 1 + 13);
    let s = json(&out.join("electron_summary.json"));
    assert_eq!(s["files"]["psf_profile"], "electron_psf_profile.csv");
    assert_eq!(s["files"]["psf_parameters"], "electron_psf_parameters.csv");
    let p = &s["results"]["psf"];
    assert_eq!(p["histories"], 160);
    assert_eq!(p["bins"], 13);
    assert!(s["results"]["deposition"].get("psf").is_none());
    // Totals: the profile energies plus what lies beyond r_max.
    let inside = csv_sum(&profile, "", 4);
    let total = p["total_ev"].as_f64().unwrap();
    let beyond = p["beyond_ev"].as_f64().unwrap();
    assert!(total > 0.0);
    assert!(
        ((inside + beyond) - total).abs() <= 1e-9 * total,
        "{inside} + {beyond} vs {total}"
    );
    // Echoed with defaults, and the echo resolves.
    assert_eq!(s["input"]["electron"]["tally"]["psf"]["fits"][0], "double");
    assert_eq!(
        s["input"]["electron"]["tally"]["psf"]["normalization"],
        "slab-total"
    );
    let echo: lindhard::input::electron::ElectronInput =
        serde_json::from_value(s["input"].clone()).unwrap();
    echo.resolve_in(&dir).unwrap();
    // Every requested fit is either reported or recorded as an error.
    let fits = p["fits"].as_array().unwrap().len() + p["fit_errors"].as_array().unwrap().len();
    assert_eq!(fits, 2);
    for f in p["fits"].as_array().unwrap() {
        assert!(f["converged"].is_boolean() && f["dof"].is_u64());
    }
}

#[test]
fn electron_psf_output_is_byte_identical_across_thread_counts() {
    let (dir, input) = electron_psf_input("electron-psf-det", PSF_SLAB);
    let mut outs = Vec::new();
    for threads in ["1", "4"] {
        let out = dir.join(format!("out-{threads}"));
        run(&input, &out, &["--histories", "80", "--threads", threads]);
        outs.push(out);
    }
    let read = |d: &Path, f: &str| std::fs::read_to_string(d.join(f)).unwrap();
    for f in ["electron_psf_profile.csv", "electron_psf_parameters.csv"] {
        assert_eq!(read(&outs[0], f), read(&outs[1], f), "{f}");
    }
    let (a, b) = (
        read(&outs[0], "electron_summary.json"),
        read(&outs[1], "electron_summary.json"),
    );
    assert_eq!(deterministic_part(&a), deterministic_part(&b));
}

#[test]
fn electron_psf_fit_error_is_reported_and_does_not_abort() {
    // A slab far below anything the electrons reach: an empty profile.
    let (dir, input) = electron_psf_input(
        "electron-psf-empty",
        "depth_lo_nm = 100000.0\ndepth_hi_nm = 200000.0\nr_min_nm = 5.0\nr_max_nm = 2000.0\nbins = 12",
    );
    let out = dir.join("out");
    run(&input, &out, &["--histories", "32"]);
    let s = json(&out.join("electron_summary.json"));
    let p = &s["results"]["psf"];
    assert_eq!(p["total_ev"], 0.0);
    assert!(p["fits"].as_array().unwrap().is_empty());
    let errs = p["fit_errors"].as_array().unwrap();
    assert_eq!(errs.len(), 2);
    assert!(errs
        .iter()
        .all(|e| !e["error"].as_str().unwrap().is_empty()));
    let profile = std::fs::read_to_string(out.join("electron_psf_profile.csv")).unwrap();
    assert_eq!(profile.lines().count(), 1 + 13);
}

#[test]
fn electron_psf_input_errors_name_the_field() {
    for (name, psf, field) in [
        (
            "slab",
            "depth_lo_nm = 5.0\ndepth_hi_nm = 5.0\nr_min_nm = 5.0\nr_max_nm = 50.0\nbins = 4",
            "electron.tally.psf.depth_hi_nm",
        ),
        (
            "radii",
            "depth_lo_nm = 0.0\ndepth_hi_nm = 5.0\nr_min_nm = 50.0\nr_max_nm = 5.0\nbins = 4",
            "electron.tally.psf.r_max_nm",
        ),
        (
            "bins",
            "depth_lo_nm = 0.0\ndepth_hi_nm = 5.0\nr_min_nm = 5.0\nr_max_nm = 50.0\nbins = 0",
            "electron.tally.psf.bins",
        ),
    ] {
        let (_, input) = electron_psf_input(&format!("electron-psf-bad-{name}"), psf);
        let o = lindhard(&["check", input.to_str().unwrap()]);
        assert!(!o.status.success(), "{name}");
        let e = String::from_utf8_lossy(&o.stderr);
        assert!(e.contains(field), "{name}: {e}");
    }
    let (_, input) = electron_psf_input(
        "electron-psf-bad-fit",
        &format!("{PSF_SLAB}\nfits = [\"quadruple\"]"),
    );
    let o = lindhard(&["check", input.to_str().unwrap()]);
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("quadruple"));
}

#[test]
fn rerun_removes_stale_electron_psf_csvs() {
    let (dir, with) = electron_psf_input("reuse-electron-psf", PSF_SLAB);
    let without = dir.join("without.toml");
    std::fs::write(
        &without,
        std::fs::read_to_string(electron_example()).unwrap(),
    )
    .unwrap();
    let out = dir.join("out");
    let files = ["electron_psf_profile.csv", "electron_psf_parameters.csv"];
    run(&with, &out, &["--histories", "32"]);
    let note = sentinel(&out);
    assert!(files.iter().all(|f| out.join(f).exists()));
    run(&without, &out, &["--histories", "32"]);
    assert!(files.iter().all(|f| !out.join(f).exists()));
    let s = json(&out.join("electron_summary.json"));
    assert!(s["files"]["psf_profile"].is_null() && s["files"]["psf_parameters"].is_null());
    assert!(s["results"].get("psf").is_none());
    // The existing key stays, as null, so the format version is unchanged.
    assert!(s["results"]["deposition"].get("psf").is_some());
    assert!(s["results"]["deposition"]["psf"].is_null());
    assert!(note.exists());
    run(&without, &out, &["--histories", "32"]);
    run(&with, &out, &["--histories", "32"]);
    assert!(files.iter().all(|f| out.join(f).exists()));
}
