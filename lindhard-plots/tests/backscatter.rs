//! Tests of the backscatter figure: the JSON-to-series mapping, byte-stable
//! rendering (against a committed golden, repeated, and across threads), the
//! committed `docs/img` figure, and the binary's `--check` mode.
//!
//! The fixture is a subset of `validation/experiments/backscatter_results.json`
//! (C and Au at 1, 5 and 10 keV, unneeded fields dropped). To regenerate the
//! golden after a deliberate change, run the tests with
//! `LINDHARD_PLOTS_BLESS=1`.

use std::path::{Path, PathBuf};
use std::process::Command;

use lindhard_plots::backscatter::{self, MeasuredPoint, RunPoint};
use lindhard_plots::{compare, CheckOutcome};

fn manifest_dir() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn fixture_path() -> PathBuf {
    manifest_dir().join("tests/fixtures/backscatter_results.json")
}

fn fixture() -> backscatter::Results {
    let json = std::fs::read_to_string(fixture_path()).expect("fixture readable");
    backscatter::parse(&json).expect("fixture parses")
}

#[test]
fn maps_json_to_sorted_series() {
    let r = fixture();
    assert_eq!(r.lindhard, "lindhard 0.0.1 (aa04511)");
    assert_eq!(r.histories, 100_000);
    assert_eq!(r.tolerance, 0.05);
    assert_eq!(r.pass_min_kev, 5.0);
    let names: Vec<&str> = r.targets.iter().map(|t| t.target.as_str()).collect();
    assert_eq!(names, ["C", "Au"], "targets keep file order");

    let c = &r.targets[0];
    // Run keys are "1", "10", "5" in lexical order; the series is numeric.
    assert_eq!(
        c.lindhard,
        [
            RunPoint {
                energy_kev: 1.0,
                eta: 0.0925,
                eta_se: 0.0009162082186926724
            },
            RunPoint {
                energy_kev: 5.0,
                eta: 0.05576,
                eta_se: 0.0007256088643339467
            },
            RunPoint {
                energy_kev: 10.0,
                eta: 0.04753,
                eta_se: 0.0006728365262082611
            },
        ]
    );
    assert_eq!(
        c.measured[1],
        MeasuredPoint {
            energy_kev: 5.0,
            median: 0.082,
            min: 0.058,
            max: 0.0857
        }
    );
    for t in &r.targets {
        for m in &t.measured {
            assert!(
                m.min <= m.median && m.median <= m.max,
                "{}: {m:?}",
                t.target
            );
        }
    }
}

#[test]
fn rejects_wrong_format_and_bad_run_key() {
    let json = std::fs::read_to_string(fixture_path()).unwrap();
    let wrong = json.replace(
        "lindhard-backscatter-results/1",
        "lindhard-backscatter-results/2",
    );
    let err = backscatter::parse(&wrong).unwrap_err();
    assert!(format!("{err:#}").contains("expected format"), "{err:#}");
    let bad_key = json.replacen("\"10\": {", "\"ten\": {", 1);
    let err = backscatter::parse(&bad_key).unwrap_err();
    assert!(format!("{err:#}").contains("is not an energy"), "{err:#}");
}

#[test]
fn fixture_matches_golden() {
    let svg = backscatter::render_svg(&fixture());
    let golden = manifest_dir().join("tests/golden/backscatter_fixture.svg");
    if std::env::var_os("LINDHARD_PLOTS_BLESS").is_some() {
        std::fs::write(&golden, &svg).unwrap();
    }
    let committed = std::fs::read(&golden).expect("golden present");
    assert_eq!(
        compare(Some(&committed), svg.as_bytes()),
        CheckOutcome::Matches,
        "rendering of the fixture changed; rerun with LINDHARD_PLOTS_BLESS=1 if deliberate"
    );
}

#[test]
fn rendering_is_byte_stable_across_calls_and_threads() {
    let results = fixture();
    let first = backscatter::render_svg(&results);
    assert_eq!(first, backscatter::render_svg(&results));
    let on_threads: Vec<String> = (0..4)
        .map(|_| {
            let r = results.clone();
            std::thread::spawn(move || backscatter::render_svg(&r))
        })
        .collect::<Vec<_>>()
        .into_iter()
        .map(|h| h.join().unwrap())
        .collect();
    for svg in on_threads {
        assert_eq!(svg, first);
    }
}

/// The figure embedded in `docs/validation.md` is what the committed results
/// render to. This runs on every OS/architecture in the CI test matrix, so it
/// is also the cross-machine determinism check.
#[test]
fn committed_docs_figure_is_current() {
    let root = manifest_dir().parent().unwrap();
    let json =
        std::fs::read_to_string(root.join("validation/experiments/backscatter_results.json"))
            .unwrap();
    let svg = backscatter::render_svg(&backscatter::parse(&json).unwrap());
    let committed = std::fs::read(root.join("docs/img/backscatter-eta.svg")).ok();
    assert_eq!(
        compare(committed.as_deref(), svg.as_bytes()),
        CheckOutcome::Matches,
        "docs/img/backscatter-eta.svg is stale or missing; run `cargo run -p lindhard-plots`"
    );
}

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("lindhard-plots-check");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name);
    let _ = std::fs::remove_file(&path);
    path
}

fn run_bin(args: &[&std::ffi::OsStr]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_lindhard-plots"))
        .args(args)
        .output()
        .expect("binary runs")
}

#[test]
fn check_mode_passes_fails_and_reports() {
    let input = fixture_path();
    let out = scratch("check.svg");
    let check = |out: &Path| {
        run_bin(&[
            "--check".as_ref(),
            "--input".as_ref(),
            input.as_os_str(),
            "--output".as_ref(),
            out.as_os_str(),
        ])
    };

    // Missing file: exit 1 and say so.
    let o = check(&out);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("is missing"));

    // Write, then check: exit 0, and the file is the library's rendering.
    let w = run_bin(&[
        "--input".as_ref(),
        input.as_os_str(),
        "--output".as_ref(),
        out.as_os_str(),
    ]);
    assert!(w.status.success(), "{}", String::from_utf8_lossy(&w.stderr));
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        backscatter::render_svg(&fixture())
    );
    let o = check(&out);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );

    // Altered file: exit 1 with the first differing byte.
    let mut bytes = std::fs::read(&out).unwrap();
    let at = bytes.len() / 2;
    bytes[at] ^= 1;
    std::fs::write(&out, &bytes).unwrap();
    let o = check(&out);
    assert_eq!(o.status.code(), Some(1));
    let msg = String::from_utf8_lossy(&o.stderr);
    assert!(
        msg.contains(&format!("first difference at byte {at}")),
        "{msg}"
    );
}

#[test]
fn compare_reports_prefix_difference() {
    assert_eq!(compare(None, b"abc"), CheckOutcome::Missing);
    assert_eq!(compare(Some(b"abc"), b"abc"), CheckOutcome::Matches);
    assert_eq!(
        compare(Some(b"ab"), b"abc"),
        CheckOutcome::Differs {
            first_difference: 2,
            committed_len: 2,
            rendered_len: 3
        }
    );
}
