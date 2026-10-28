//! Tests of the stopping-residual report: the committed book page is current
//! (golden), the rendering is repeatable, a modified dataset makes the page
//! stale, skipped points are stated, and a forbidden source is refused.

use std::fs;
use std::path::{Path, PathBuf};

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
        "Points: 15. Energy span: 50 keV to 600 keV.",
        "Phys. Rev. 103, 32 (1956)",
        "| lindhard-scharff | 15 | 0 |",
        "**held out**",
        "lowers",
    ] {
        assert!(page.contains(needle), "page lacks {needle:?}");
    }
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
