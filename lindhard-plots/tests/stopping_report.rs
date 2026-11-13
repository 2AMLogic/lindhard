//! Tests of the stopping-residual report: the committed book page is current
//! (golden), the rendering is repeatable, a modified dataset makes the page
//! stale, skipped points are stated, a forbidden source is refused, and the
//! binary's `--check` mode reports and exits as it should.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use lindhard_plots::stopping_report::render;
use lindhard_plots::{compare, CheckOutcome};

const LABEL: &str = "validation/data/stopping";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

/// A scratch directory holding `files`, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str, files: &[(&str, String)]) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "lindhard-stopping-report-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        for (name, text) in files {
            fs::write(dir.join(name), text).unwrap();
        }
        Self(dir)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

const HEADER: &str = "# format: lindhard-stopping-data/1\n# name: synthetic\n\
    # description: invented numbers for the report tests, not measurements\n\
    ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method,citation\n";

/// Eight invented proton-in-Si points; the lowest energies are below the
/// Bethe-Bloch validity range, so that model skips them.
fn synthetic(scale: f64) -> String {
    let mut text = HEADER.to_string();
    for (i, e) in [1.0e3, 3.0e3, 1.0e4, 3.0e4, 1.0e5, 3.0e5, 1.0e6, 3.0e6]
        .iter()
        .enumerate()
    {
        let s = scale * (10.0 + i as f64);
        text.push_str(&format!(
            "1,,Si,{e},{s},0.5,invented,\"Test Author, Test J. 1, 1 (2000), Table {i}\"\n"
        ));
    }
    text
}

#[test]
fn committed_page_is_what_the_generator_renders() {
    let page = render(&root().join(LABEL), LABEL).unwrap();
    let committed = fs::read(root().join("book/src/models/stopping-residuals.md")).ok();
    assert_eq!(
        compare(committed.as_deref(), page.as_bytes()),
        CheckOutcome::Matches,
        "regenerate with `cargo run -p lindhard-plots --bin lindhard-stopping-report`"
    );
    // The committed dataset is in the page, with its per-system statistics.
    for needle in [
        "bader1956_p_cu.csv",
        "Points: 15. Energy span: 50\u{a0}keV to 600\u{a0}keV.",
        "## H on Cu, Bader 1956",
        "File: `bader1956_p_cu.csv`. Name: `bader1956-p-cu`.",
        "## Summary",
        "## H on Si, Grahmann 1976 \u{2014} srim.org",
        "|---|---|---:|---:|---:|",
        "Phys. Rev. 103, 32 (1956)",
        "| lindhard-scharff | 15 | 0 |",
        "**held out**",
        "lowers",
    ] {
        assert!(page.contains(needle), "page lacks {needle:?}");
    }
}

#[test]
fn summary_has_one_row_per_dataset_and_model() {
    let page = render(&root().join(LABEL), LABEL).unwrap();
    let summary: Vec<&str> = page
        .split("## Summary\n")
        .nth(1)
        .unwrap()
        .split("\n## ")
        .next()
        .unwrap()
        .lines()
        .filter(|l| l.starts_with("| ") && !l.starts_with("| Dataset"))
        .collect();
    let datasets = fs::read_dir(root().join(LABEL))
        .unwrap()
        .filter(|e| {
            let p = e.as_ref().unwrap().path();
            p.extension().is_some_and(|x| x == "csv")
        })
        .count();
    assert_eq!(summary.len(), 4 * datasets);
}

#[test]
fn rendering_twice_gives_identical_bytes() {
    let dir = Scratch::new("twice", &[("a.csv", synthetic(1.0))]);
    let a = render(&dir.0, "data").unwrap();
    let b = render(&dir.0, "data").unwrap();
    assert_eq!(a, b);
}

#[test]
fn a_modified_dataset_makes_the_page_stale() {
    let before = Scratch::new("stale-a", &[("a.csv", synthetic(1.0))]);
    let after = Scratch::new("stale-b", &[("a.csv", synthetic(1.01))]);
    let committed = render(&before.0, "data").unwrap();
    let fresh = render(&after.0, "data").unwrap();
    assert!(matches!(
        compare(Some(committed.as_bytes()), fresh.as_bytes()),
        CheckOutcome::Differs { .. }
    ));
}

#[test]
fn skipped_points_are_counted_listed_and_kept_out_of_the_statistics() {
    let dir = Scratch::new("skip", &[("a.csv", synthetic(1.0))]);
    let page = render(&dir.0, "data").unwrap();
    assert!(page.contains("Skipped points"), "{page}");
    let skipped_lines: Vec<&str> = page
        .lines()
        .filter(|l| l.starts_with("* bethe-bloch: line"))
        .collect();
    assert!(!skipped_lines.is_empty(), "{page}");
    assert!(skipped_lines.iter().all(|l| l.contains("skipped:")));
    let row = page
        .lines()
        .find(|l| l.starts_with("| bethe-bloch | "))
        .unwrap();
    let cols: Vec<&str> = row.split('|').map(str::trim).collect();
    let (n, skipped): (usize, usize) = (cols[2].parse().unwrap(), cols[3].parse().unwrap());
    assert_eq!(n + skipped, 8);
    assert_eq!(skipped, skipped_lines.len());
    assert!(skipped > 0);
    // A model that cannot evaluate every point is not fitted, and says why.
    assert!(page.contains("no correction fitted"), "{page}");
}

#[test]
fn a_system_with_too_few_points_states_that_no_correction_was_fitted() {
    let mut text = HEADER.to_string();
    for e in [1.0e4, 2.0e4, 4.0e4] {
        text.push_str(&format!(
            "1,,Si,{e},10.0,0,invented,\"Test, J. 1 (2000)\"\n"
        ));
    }
    let dir = Scratch::new("few", &[("a.csv", text)]);
    let page = render(&dir.0, "data").unwrap();
    assert!(page.contains("no correction fitted"), "{page}");
}

#[test]
fn a_dataset_citing_a_forbidden_source_is_refused() {
    for word in ["NIST", "IAEA", "SRIM", "PSTAR", "ASTAR"] {
        let text = synthetic(1.0).replace("Test Author", &format!("{word} database"));
        let dir = Scratch::new(&format!("forbid-{word}"), &[("a.csv", text)]);
        let err = render(&dir.0, "data").unwrap_err().to_string();
        assert!(err.contains("may not enter"), "{word}: {err}");
    }
}

#[test]
fn a_name_that_merely_contains_the_letters_is_not_refused() {
    let text = synthetic(1.0).replace("Test Author", "Srimson and Nistor");
    let dir = Scratch::new("substring", &[("a.csv", text)]);
    assert!(render(&dir.0, "data").is_ok());
}

const ORIGIN_HEADER: &str = "# format: lindhard-stopping-data/1\n# name: compilation\n\
    # description: invented points standing in for a digitized srim.org compilation\n\
    ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method,citation,\
    citation_key,origin\n";

fn compilation(origin: &str) -> String {
    let mut text = ORIGIN_HEADER.to_string();
    for (i, e) in [1.0e4, 3.0e4, 1.0e5, 3.0e5, 1.0e6, 3.0e6]
        .iter()
        .enumerate()
    {
        text.push_str(&format!(
            "1,,Si,{e},{},0.5,digitized plot,\"Digitized from the srim.org compilation: Test, J. 1 (2000)\",2000-Test,{origin}\n",
            10.0 + i as f64
        ));
    }
    text
}

#[test]
fn a_compilation_origin_file_may_name_its_source() {
    let dir = Scratch::new(
        "compilation",
        &[("a.csv", compilation("digitized-compilation"))],
    );
    assert!(render(&dir.0, "data").is_ok());
}

#[test]
fn an_unknown_origin_is_refused_by_the_loader() {
    let dir = Scratch::new("bad-origin", &[("a.csv", compilation("table"))]);
    let err = format!("{:#}", render(&dir.0, "data").unwrap_err());
    assert!(err.contains("unknown origin"), "{err}");
}

#[test]
fn an_eight_column_file_mentioning_srim_is_still_refused() {
    let text = synthetic(1.0).replace("Test Author", "SRIM compilation");
    let dir = Scratch::new("srim-eight", &[("a.csv", text)]);
    let err = render(&dir.0, "data").unwrap_err().to_string();
    assert!(err.contains("may not enter"), "{err}");
}

/// The binary's `--check` mode: the three `CheckOutcome` messages, the exit
/// codes, and the regenerate command named in each failure. The shared driver
/// in `lindhard_plots` (`parse_args` / `write_or_check` / `exit_code`) is what
/// produces these, so this is also the guard that it keeps saying the same
/// things for this binary as for `lindhard-plots`.
#[test]
fn check_mode_passes_fails_and_names_the_regenerate_command() {
    const REGENERATE: &str = "regenerate it with `cargo run -p lindhard-plots \
                              --bin lindhard-stopping-report` and commit it";
    let dir = Scratch::new("bin-check", &[("a.csv", synthetic(1.0))]);
    let out = dir.0.join("page.md");
    let run = |args: &[&std::ffi::OsStr]| {
        Command::new(env!("CARGO_BIN_EXE_lindhard-stopping-report"))
            .args(args)
            .output()
            .expect("binary runs")
    };
    let check = || {
        run(&[
            "--check".as_ref(),
            "--input".as_ref(),
            dir.0.as_os_str(),
            "--output".as_ref(),
            out.as_os_str(),
        ])
    };

    // Missing file: exit 1, and say which command regenerates it.
    let o = check();
    assert_eq!(o.status.code(), Some(1));
    let msg = String::from_utf8_lossy(&o.stderr);
    assert!(msg.contains("is missing"), "{msg}");
    assert!(msg.contains(REGENERATE), "{msg}");

    // Write, then check: exit 0, and the file is the library's rendering with
    // the given `--input` as the label the page records.
    let w = run(&[
        "--input".as_ref(),
        dir.0.as_os_str(),
        "--output".as_ref(),
        out.as_os_str(),
    ]);
    assert!(w.status.success(), "{}", String::from_utf8_lossy(&w.stderr));
    assert!(
        String::from_utf8_lossy(&w.stdout).contains("wrote "),
        "{}",
        String::from_utf8_lossy(&w.stdout)
    );
    assert_eq!(
        fs::read_to_string(&out).unwrap(),
        render(&dir.0, &dir.0.display().to_string()).unwrap()
    );
    let o = check();
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(
        String::from_utf8_lossy(&o.stdout).contains("is up to date"),
        "{}",
        String::from_utf8_lossy(&o.stdout)
    );

    // Altered file: exit 1 with the first differing byte, and the command.
    let mut bytes = fs::read(&out).unwrap();
    let at = bytes.len() / 2;
    bytes[at] ^= 1;
    fs::write(&out, &bytes).unwrap();
    let o = check();
    assert_eq!(o.status.code(), Some(1));
    let msg = String::from_utf8_lossy(&o.stderr);
    assert!(msg.contains("differs from a fresh render of"), "{msg}");
    assert!(
        msg.contains(&format!("first difference at byte {at}")),
        "{msg}"
    );
    assert!(msg.contains(REGENERATE), "{msg}");

    // An unexpected argument is an error (exit 2) quoting the usage line.
    let o = run(&["--nope".as_ref()]);
    assert_eq!(o.status.code(), Some(2));
    let msg = String::from_utf8_lossy(&o.stderr);
    assert!(msg.contains("unexpected argument"), "{msg}");
    assert!(msg.contains("usage: lindhard-stopping-report"), "{msg}");
}
