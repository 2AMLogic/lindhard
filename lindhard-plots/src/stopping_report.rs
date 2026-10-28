//! The per-system residual report of the electronic-stopping models against
//! every committed measured stopping dataset (`lindhard-stopping-data/1`,
//! `docs/stopping-data-format.md`), as the Markdown page
//! `book/src/models/stopping-residuals.md`.
//!
//! Like the figures, the page is a pure function of committed files (the
//! datasets under `validation/data/stopping/` and the models in `lindhard`):
//! nothing is read from the environment and the text is identical from run to
//! run and whatever the thread count, so CI re-renders it and compares with
//! the committed copy (`--check` of the `lindhard-stopping-report` binary).
//!
//! Residuals come from [`lindhard::ion::stopping::dataset::compare`]. Where a
//! system is a single ion on a single element, a per-pair correction of each
//! model is also fitted with [`lindhard::ion::stopping::correction`] and the
//! uncorrected and corrected statistics on the held-out points are reported
//! side by side, with the split marked. The fit settings are the constants
//! [`FIT_SETTINGS`]; no bootstrap band is computed (it is not part of the
//! report).

use std::fmt::Write as _;
use std::path::Path;

use anyhow::{bail, Context, Result};
use lindhard::elements::element;
use lindhard::ion::stopping::bethe::BetheBloch;
use lindhard::ion::stopping::bragg::NoCorrection;
use lindhard::ion::stopping::correction::{fit_correction, FitSettings, ResidualStats};
use lindhard::ion::stopping::dataset::{compare, StoppingDataset, StoppingPoint, Target};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::mix::EquipartitionMix;
use lindhard::ion::stopping::oen_robinson::OenRobinson;
use lindhard::ion::stopping::ElectronicStopping;

/// Settings of the per-pair correction fits shown in the report.
pub const FIT_SETTINGS: FitSettings = FitSettings {
    degree: 2,
    held_out: 3,
    replicates: 0,
    seed: 20_261_015,
};

/// Words that mark a source the project may not use (CONTRIBUTING.md, "Tier
/// C", and `docs/data-provenance.md`). A dataset whose method, citation or
/// header mentions one is refused, so such values cannot reach the report.
const FORBIDDEN_SOURCES: [&str; 5] = ["nist", "iaea", "srim", "pstar", "astar"];

/// The models compared, in report order.
fn models() -> Vec<Box<dyn ElectronicStopping>> {
    vec![
        Box::new(LindhardScharff::new()),
        Box::new(OenRobinson::new()),
        Box::new(EquipartitionMix::new()),
        Box::new(BetheBloch::new()),
    ]
}

/// Render the report for every `*.csv` in `data_dir`, in file-name order.
/// `data_label` is how the directory is named in the page.
pub fn render(data_dir: &Path, data_label: &str) -> Result<String> {
    let mut files: Vec<_> = std::fs::read_dir(data_dir)
        .with_context(|| format!("reading {}", data_dir.display()))?
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "csv"))
        .collect();
    files.sort();
    let mut datasets = Vec::new();
    for path in &files {
        let text =
            std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .context("dataset file name is not UTF-8")?
            .to_string();
        refuse_forbidden_sources(&name, &text)?;
        let ds = StoppingDataset::from_csv_str(&text)
            .with_context(|| format!("in {}", path.display()))?;
        datasets.push(Loaded { name, text, ds });
    }
    Ok(page(&datasets, data_label))
}

fn refuse_forbidden_sources(file: &str, text: &str) -> Result<()> {
    let lower = text.to_lowercase();
    for word in FORBIDDEN_SOURCES {
        // Whole words only, so that an author's name containing the letters
        // does not trip the guard.
        let hit = lower
            .split(|c: char| !c.is_ascii_alphanumeric())
            .any(|w| w == word);
        if hit {
            bail!(
                "{file} mentions {word:?}; data from NIST, IAEA or SRIM-derived sources may not \
                 enter this repository (CONTRIBUTING.md, docs/data-provenance.md)"
            );
        }
    }
    Ok(())
}

/// A dataset file as read: its name, its text and its parse.
struct Loaded {
    name: String,
    text: String,
    ds: StoppingDataset,
}

/// One ion-target system of a dataset: indices into `points()`.
struct System {
    z: u8,
    mass_bits: u64,
    target: String,
    indices: Vec<usize>,
}

fn systems_of(ds: &StoppingDataset) -> Vec<System> {
    let mut out: Vec<System> = Vec::new();
    for (i, p) in ds.points().iter().enumerate() {
        let (z, bits, label) = (
            p.ion().z(),
            p.ion().mass_amu().to_bits(),
            p.target().label(),
        );
        match out
            .iter_mut()
            .find(|s| s.z == z && s.mass_bits == bits && s.target == label)
        {
            Some(s) => s.indices.push(i),
            None => out.push(System {
                z,
                mass_bits: bits,
                target: label.to_string(),
                indices: vec![i],
            }),
        }
    }
    out
}

fn symbol(z: u8) -> &'static str {
    element(z).map_or("?", |e| e.symbol)
}

/// A table cell must not break the row.
fn cell(text: &str) -> String {
    text.replace('|', "\\|")
}

fn pct(x: f64) -> String {
    if x.is_nan() {
        "n/a".into()
    } else {
        format!("{:+.2} %", 100.0 * x)
    }
}

fn pct_unsigned(x: f64) -> String {
    if x.is_nan() {
        "n/a".into()
    } else {
        format!("{:.2} %", 100.0 * x)
    }
}

fn val(x: f64) -> String {
    if x.is_nan() {
        "n/a".into()
    } else {
        format!("{x:+.3}")
    }
}

fn val_unsigned(x: f64) -> String {
    if x.is_nan() {
        "n/a".into()
    } else {
        format!("{x:.3}")
    }
}

fn energy(ev: f64) -> String {
    if ev >= 1.0e6 {
        format!("{} MeV", trim(ev / 1.0e6))
    } else if ev >= 1.0e3 {
        format!("{} keV", trim(ev / 1.0e3))
    } else {
        format!("{} eV", trim(ev))
    }
}

fn trim(x: f64) -> String {
    let s = format!("{x:.3}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

fn distinct(values: impl Iterator<Item = String>) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for v in values {
        if !out.contains(&v) {
            out.push(v);
        }
    }
    out
}

fn page(datasets: &[Loaded], data_label: &str) -> String {
    let mut o = String::new();
    o.push_str(
        "<!-- Generated by `cargo run -p lindhard-plots --bin lindhard-stopping-report`. \
         Do not edit by hand: CI regenerates this page and fails if it differs. -->\n",
    );
    o.push_str("# Stopping residuals against measurements\n\n");
    let _ = writeln!(
        o,
        "Each electronic-stopping model of the [stopping chapter](stopping.md) compared with every \
         measured stopping dataset committed under `{data_label}/` (format: \
         `lindhard-stopping-data/1`, `docs/stopping-data-format.md`). Only data admitted by \
         `docs/data-provenance.md` are in that directory; no NIST, IAEA or SRIM-derived value is \
         read or shown.\n"
    );
    o.push_str(
        "Regenerate this page with\n\n```text\ncargo run -p lindhard-plots --bin \
         lindhard-stopping-report\n```\n\nand check it is current with the same command and \
         `-- --check`. The output is a pure function of the committed datasets and the models, \
         so running it twice gives the same bytes.\n\n",
    );
    o.push_str("## How to read the tables\n\n");
    o.push_str(
        "* A residual is `model - measured`; positive means the model is **above** the \
         measurement. Relative residuals are `model / measured - 1`.\n\
         * Stopping is a cross section per target atom (per average atom for a compound), in \
         eV\u{b7}10\u{207b}\u{b9}\u{2075} cm\u{b2}. *Mean* and *rms* are the mean and root mean \
         square of the absolute residual in that unit; *rel. mean* and *rel. rms* are those of \
         the relative residual.\n\
         * A point the model cannot evaluate is **skipped**: it is counted in the *skipped* \
         column, listed with the model's reason, and is in no statistic.\n",
    );
    let _ = writeln!(
        o,
        "* Where a system is one ion on one element, a per-pair correction of each model \
         (a Chebyshev polynomial of degree {} in ln E for the log of the measured-to-model \
         ratio, `lindhard::ion::stopping::correction`) is fitted to the training points. \
         {} interior points chosen with seed {} are held out of the fit. The corrected and \
         uncorrected statistics on the **held-out** points are shown side by side; whether the \
         correction helps is stated either way. The correction applies only inside the \
         energy span of the training points.\n",
        FIT_SETTINGS.degree, FIT_SETTINGS.held_out, FIT_SETTINGS.seed
    );
    o.push_str("## Datasets\n\n| File | Name | Points | Systems |\n|---|---|---|---|\n");
    for Loaded { name: file, ds, .. } in datasets {
        let _ = writeln!(
            o,
            "| `{}` | {} | {} | {} |",
            file,
            cell(ds.name()),
            ds.points().len(),
            systems_of(ds).len()
        );
    }
    o.push('\n');
    if datasets.is_empty() {
        o.push_str("No dataset is committed.\n");
    }
    let models = models();
    for d in datasets {
        dataset_section(&mut o, d, &models);
    }
    o
}

fn dataset_section(o: &mut String, d: &Loaded, models: &[Box<dyn ElectronicStopping>]) {
    let Loaded { name: file, ds, .. } = d;
    let _ = writeln!(o, "## {}\n", cell(ds.name()));
    let _ = writeln!(o, "File: `{file}`.\n");
    if let Some(d) = ds.description() {
        let _ = writeln!(o, "{d}\n");
    }
    for key in ["source", "note"] {
        if let Some((_, v)) = ds.metadata().iter().find(|(k, _)| k == key) {
            let _ = writeln!(
                o,
                "{}: {v}\n",
                if key == "source" { "Source" } else { "Note" }
            );
        }
    }
    let comparisons: Vec<_> = models
        .iter()
        .map(|m| compare(m.as_ref(), &NoCorrection, ds))
        .collect();
    for sys in systems_of(ds) {
        let pts: Vec<&StoppingPoint> = sys.indices.iter().map(|&i| &ds.points()[i]).collect();
        let first = pts[0];
        let _ = writeln!(
            o,
            "### {} ({} u) on {}\n",
            symbol(sys.z),
            trim(first.ion().mass_amu()),
            sys.target
        );
        let (lo, hi) = pts.iter().fold((f64::INFINITY, 0.0_f64), |(l, h), p| {
            (l.min(p.energy_ev()), h.max(p.energy_ev()))
        });
        let _ = writeln!(
            o,
            "* Points: {}. Energy span: {} to {}.",
            pts.len(),
            energy(lo),
            energy(hi)
        );
        let sigmas = pts.iter().filter(|p| p.sigma_ev_1e15_cm2() > 0.0).count();
        let _ = writeln!(
            o,
            "* Stated uncertainty: {}.",
            match sigmas {
                0 => "none (sigma 0 on every point)".to_string(),
                n if n == pts.len() => "on every point".to_string(),
                n => format!("on {n} of {} points", pts.len()),
            }
        );
        for m in distinct(pts.iter().map(|p| p.method().to_string())) {
            let _ = writeln!(o, "* Method: {m}.");
        }
        for c in distinct(pts.iter().map(|p| p.citation().to_string())) {
            let _ = writeln!(o, "* Citation: {c}.");
        }
        o.push_str(
            "\n| Model | Points | Skipped | Mean | rms | Rel. mean | Rel. rms |\n\
             |---|---|---|---|---|---|---|\n",
        );
        let mut skipped_lines: Vec<String> = Vec::new();
        for (m, cmp) in models.iter().zip(&comparisons) {
            let pairs: Vec<(f64, f64)> = cmp
                .residuals
                .iter()
                .filter(|r| sys.indices.contains(&r.index))
                .map(|r| (r.model_ev_1e15_cm2, r.measured_ev_1e15_cm2))
                .collect();
            let st = ResidualStats::of(&pairs);
            let skipped: Vec<_> = cmp
                .skipped
                .iter()
                .filter(|s| sys.indices.contains(&s.index))
                .collect();
            let _ = writeln!(
                o,
                "| {} | {} | {} | {} | {} | {} | {} |",
                m.name(),
                st.n,
                skipped.len(),
                val(st.mean_absolute_ev_1e15_cm2),
                val_unsigned(st.rms_absolute_ev_1e15_cm2),
                pct(st.mean_relative),
                pct_unsigned(st.rms_relative),
            );
            for s in skipped {
                skipped_lines.push(format!(
                    "* {}: line {} ({}) skipped: {}.",
                    m.name(),
                    s.line,
                    energy(ds.points()[s.index].energy_ev()),
                    s.reason
                ));
            }
        }
        o.push('\n');
        if !skipped_lines.is_empty() {
            o.push_str("Skipped points (in no statistic above):\n\n");
            for l in skipped_lines {
                let _ = writeln!(o, "{l}");
            }
            o.push('\n');
        }
        corrections(o, d, &sys, &pts, models);
    }
}

fn corrections(
    o: &mut String,
    d: &Loaded,
    sys: &System,
    pts: &[&StoppingPoint],
    models: &[Box<dyn ElectronicStopping>],
) {
    if !matches!(pts[0].target(), Target::Element(_)) {
        o.push_str("No per-pair correction is fitted: the target is a compound.\n\n");
        return;
    }
    // The fit takes a dataset of one system, so it gets this system's rows.
    let Some(sub) = subset(d, &sys.indices) else {
        let _ = writeln!(
            o,
            "No per-pair correction is fitted: the system could not be isolated.\n"
        );
        return;
    };
    o.push_str("#### Held-out residuals, uncorrected and corrected\n\n");
    let mut rows = String::new();
    let mut notes = String::new();
    // (table text, models whose table it is): models that agree to the
    // printed digits share one table.
    let mut splits: Vec<(String, Vec<&'static str>)> = Vec::new();
    for m in models {
        match fit_correction(m.as_ref(), &sub, FIT_SETTINGS) {
            Err(e) => {
                let _ = writeln!(notes, "* {}: no correction fitted ({e}).", m.name());
            }
            Ok(fit) => {
                let rep = fit.held_out_report();
                let (a, b) = (rep.base, rep.corrected);
                let _ = writeln!(
                    rows,
                    "| {} | {} | {} | {} | {} | {} | {} | {} |",
                    m.name(),
                    rep.points.len(),
                    pct(a.mean_relative),
                    pct(b.mean_relative),
                    pct_unsigned(a.rms_relative),
                    pct_unsigned(b.rms_relative),
                    if rep.reduces_rms_relative() {
                        "lowers"
                    } else {
                        "does not lower"
                    },
                    {
                        let (l, h) = fit.span_ev();
                        format!("{} to {}", energy(l), energy(h))
                    },
                );
                let mut table = String::from(
                    "| Line | Energy | Measured | Uncorrected | Corrected | Split |\n\
                     |---|---|---|---|---|---|\n",
                );
                for (k, p) in sub.points().iter().enumerate() {
                    let held = fit.held_out_indices().contains(&k);
                    let (base, corr) = model_values(m.as_ref(), &fit, p);
                    let _ = writeln!(
                        table,
                        "| {} | {} | {} | {} | {} | {} |",
                        p.line(),
                        energy(p.energy_ev()),
                        val_unsigned(p.s_e_ev_1e15_cm2()),
                        val_unsigned(base),
                        val_unsigned(corr),
                        if held { "**held out**" } else { "train" },
                    );
                }
                match splits.iter_mut().find(|(t, _)| *t == table) {
                    Some((_, names)) => names.push(m.name()),
                    None => splits.push((table, vec![m.name()])),
                }
            }
        }
    }
    if !rows.is_empty() {
        o.push_str(
            "| Model | Held out | Rel. mean, uncorrected | Rel. mean, corrected | Rel. rms, \
             uncorrected | Rel. rms, corrected | Correction vs rel. rms | Correction span |\n\
             |---|---|---|---|---|---|---|---|\n",
        );
        o.push_str(&rows);
        o.push('\n');
    }
    if !notes.is_empty() {
        o.push_str(&notes);
        o.push('\n');
    }
    for (table, names) in &splits {
        let _ = writeln!(
            o,
            "Per point, for {}: held-out points are marked **held out**; the others trained the \
             fit.\n",
            names
                .iter()
                .map(|n| format!("`{n}`"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        o.push_str(table);
        o.push('\n');
    }
}

/// Base and corrected model values at a point (stopping in
/// eV\u{b7}10\u{207b}\u{b9}\u{2075} cm\u{b2}); NaN when the model declines.
fn model_values(
    base: &dyn ElectronicStopping,
    fit: &lindhard::ion::stopping::correction::CorrectionFit,
    p: &StoppingPoint,
) -> (f64, f64) {
    let Target::Element(z) = *p.target() else {
        return (f64::NAN, f64::NAN);
    };
    let Ok(s) = base.stopping(p.ion(), z, p.energy_ev()) else {
        return (f64::NAN, f64::NAN);
    };
    let s = lindhard::ion::stopping::to_ev_1e15_cm2(s);
    (s, s * fit.factor(p.energy_ev()).unwrap_or(1.0))
}

/// The rows `indices` of the dataset `d` as a dataset of their own: the
/// original text with the other points' lines removed.
fn subset(d: &Loaded, indices: &[usize]) -> Option<StoppingDataset> {
    let keep: Vec<usize> = indices.iter().map(|&i| d.ds.points()[i].line()).collect();
    let all: Vec<usize> = d.ds.points().iter().map(StoppingPoint::line).collect();
    let text: String = d
        .text
        .lines()
        .enumerate()
        .filter(|(i, _)| !all.contains(&(i + 1)) || keep.contains(&(i + 1)))
        .map(|(_, l)| format!("{l}\n"))
        .collect();
    StoppingDataset::from_csv_str(&text).ok()
}
