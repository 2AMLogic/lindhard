//! Stopping datasets: loader validation and `compare` against synthetic points.
//!
//! Every dataset here is synthetic. The "measured" values are sampled from this
//! repository's own Lindhard-Scharff model (and Bragg sums of it), or are
//! made-up numbers in rejection fixtures. No real measurement is used (see the
//! row in `docs/data-provenance.md`).

use lindhard::ion::stopping::bragg::{
    bragg_cross_section_per_atom, CompoundCorrection, ConstantCorrection, NoCorrection,
};
use lindhard::ion::stopping::dataset::{compare, StoppingDataset, Target};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::{
    to_ev_1e15_cm2, ElectronicStopping, Ion, StoppingError, ValidityRange,
};

const HEADER: &str = "# format: lindhard-stopping-data/1\n\
# name: synthetic-test\n\
# description: synthetic points from this repository's Lindhard-Scharff model\n";
const COLUMNS: &str =
    "ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method,citation\n";
const CITATION: &str = "\"Lindhard-Scharff model of this repository, sampled (synthetic, test)\"";

/// `(ion Z, target)` systems of the synthetic dataset: two element targets and
/// one compound (Bragg path).
const SYSTEMS: [(u8, &str); 3] = [(5, "Si"), (15, "Si"), (5, "SiO2")];

/// Model value for one point, eV 1e-15 cm^2 per (average) atom.
fn model_value(
    model: &dyn ElectronicStopping,
    correction: &dyn CompoundCorrection,
    ion_z: u8,
    target: &str,
    energy_ev: f64,
) -> f64 {
    let ion = Ion::new(ion_z).unwrap();
    let s = match Target::parse(target).unwrap() {
        Target::Element(z) => model.stopping(&ion, z, energy_ev).unwrap(),
        Target::Compound { material, .. } => {
            bragg_cross_section_per_atom(model, correction, &ion, &material, energy_ev).unwrap()
        }
    };
    to_ev_1e15_cm2(s)
}

/// A dataset sampled from plain Lindhard-Scharff with every value multiplied
/// by `scale`, on a log grid from 1 keV to 100 keV.
fn ls_dataset(scale: f64) -> String {
    let ls = LindhardScharff::new();
    let mut text = format!("{HEADER}{COLUMNS}");
    for (ion_z, target) in SYSTEMS {
        for i in 0..=8 {
            let e = 1.0e3 * 100f64.powf(f64::from(i) / 8.0);
            let s = scale * model_value(&ls, &NoCorrection, ion_z, target, e);
            text.push_str(&format!(
                "{ion_z},,{target},{e:e},{s:e},{:e},synthetic,{CITATION}\n",
                0.05 * s
            ));
        }
    }
    text
}

#[test]
fn self_comparison_has_zero_residual() {
    let ds = StoppingDataset::from_csv_str(&ls_dataset(1.0)).unwrap();
    assert_eq!(ds.name(), "synthetic-test");
    assert_eq!(ds.points().len(), 27);
    let c = compare(&LindhardScharff::new(), &NoCorrection, &ds);
    assert_eq!(c.model, "lindhard-scharff");
    assert!(c.skipped.is_empty());
    assert_eq!(c.residuals.len(), 27);
    assert_eq!(c.overall.n, 27);
    assert!(c.overall.max_abs_relative < 1e-14, "{:?}", c.overall);
    assert!(c.overall.rms_relative < 1e-14);
    for r in &c.residuals {
        assert!(r.absolute_ev_1e15_cm2.abs() < 1e-12 * r.measured_ev_1e15_cm2);
        assert!(r.pull.unwrap().abs() < 1e-10);
    }
    let labels: Vec<_> = c
        .systems
        .iter()
        .map(|s| (s.ion_z, s.target.as_str(), s.summary.n))
        .collect();
    assert_eq!(labels, [(5, "Si", 9), (15, "Si", 9), (5, "SiO2", 9)]);
}

#[test]
fn scaled_dataset_gives_the_scale_factor() {
    // measured = k * model, so model / measured - 1 = 1/k - 1 at every point,
    // positive because the model sits above the measurement.
    let k = 0.8;
    let ds = StoppingDataset::from_csv_str(&ls_dataset(k)).unwrap();
    let c = compare(&LindhardScharff::new(), &NoCorrection, &ds);
    let want = 1.0 / k - 1.0;
    assert!((c.overall.mean_relative - want).abs() < 1e-13);
    assert!((c.overall.rms_relative - want).abs() < 1e-13);
    assert!((c.overall.max_abs_relative - want).abs() < 1e-13);
    for s in &c.systems {
        assert!((s.summary.mean_relative - want).abs() < 1e-13, "{s:?}");
    }
    for r in &c.residuals {
        assert!(r.relative > 0.0);
        let measured = r.measured_ev_1e15_cm2;
        assert!((r.absolute_ev_1e15_cm2 - (1.0 / k - 1.0) * measured).abs() < 1e-12 * measured);
        // sigma = 5 % of measured: pull = (1/k - 1) / 0.05.
        assert!((r.pull.unwrap() - want / 0.05).abs() < 1e-10);
    }
}

#[test]
fn compound_correction_applies_to_compounds_only() {
    let ds = StoppingDataset::from_csv_str(&ls_dataset(1.0)).unwrap();
    let c = compare(&LindhardScharff::new(), &ConstantCorrection(1.1), &ds);
    for s in &c.systems {
        let want = if s.target == "SiO2" { 0.1 } else { 0.0 };
        assert!((s.summary.mean_relative - want).abs() < 1e-13, "{s:?}");
    }
}

#[test]
fn compound_value_is_per_average_atom() {
    let ls = LindhardScharff::new();
    let ion = Ion::new(5).unwrap();
    let e = 2.0e4;
    let hand = to_ev_1e15_cm2(
        (ls.stopping(&ion, 14, e).unwrap() + 2.0 * ls.stopping(&ion, 8, e).unwrap()) / 3.0,
    );
    let text = format!("{HEADER}{COLUMNS}5,,SiO2,{e:e},{hand:e},0,synthetic,{CITATION}\n");
    let c = compare(
        &ls,
        &NoCorrection,
        &StoppingDataset::from_csv_str(&text).unwrap(),
    );
    assert!(c.residuals[0].relative.abs() < 1e-14);
    assert_eq!(c.residuals[0].pull, None);
}

/// Lindhard-Scharff that refuses energies below 5 keV, standing in for a
/// model with a lower validity limit.
struct LimitedLs;

impl ElectronicStopping for LimitedLs {
    fn name(&self) -> &'static str {
        "limited-ls"
    }
    fn stopping(&self, ion: &Ion, target_z: u8, energy_ev: f64) -> Result<f64, StoppingError> {
        if energy_ev < 5.0e3 {
            return Err(StoppingError::NotApplicable {
                model: "limited-ls",
                energy_ev,
            });
        }
        LindhardScharff::new().stopping(ion, target_z, energy_ev)
    }
    fn validity(&self, _ion: &Ion) -> ValidityRange {
        ValidityRange {
            min_energy_ev: 5.0e3,
            max_energy_ev: f64::INFINITY,
        }
    }
}

#[test]
fn not_applicable_points_are_reported_as_skipped() {
    let ds = StoppingDataset::from_csv_str(&ls_dataset(1.0)).unwrap();
    let c = compare(&LimitedLs, &NoCorrection, &ds);
    // Grid energies 1e3 * 100^(i/8): i = 0, 1, 2 are below 5 keV.
    assert_eq!(c.skipped.len(), 9);
    assert_eq!(c.residuals.len(), 18);
    assert_eq!(c.overall.n, 18);
    assert_eq!(c.overall.n_skipped, 9);
    for s in &c.systems {
        assert_eq!((s.summary.n, s.summary.n_skipped), (6, 3));
    }
    for s in &c.skipped {
        assert!(matches!(s.reason, StoppingError::NotApplicable { .. }));
        assert_eq!(s.line, ds.points()[s.index].line());
    }
    assert!(c.overall.max_abs_relative < 1e-14);
}

#[test]
fn all_skipped_gives_nan_statistics() {
    let text = format!("{HEADER}{COLUMNS}5,,Si,1000,5,0.1,synthetic,{CITATION}\n");
    let c = compare(
        &LimitedLs,
        &NoCorrection,
        &StoppingDataset::from_csv_str(&text).unwrap(),
    );
    assert_eq!((c.overall.n, c.overall.n_skipped), (0, 1));
    assert!(c.overall.mean_relative.is_nan() && c.overall.rms_relative.is_nan());
}

#[test]
fn summary_statistics_are_bitwise_deterministic() {
    let text = ls_dataset(0.93);
    let run = || {
        let ds = StoppingDataset::from_csv_str(&text).unwrap();
        compare(&LindhardScharff::new(), &ConstantCorrection(1.05), &ds)
    };
    let c0 = run();
    let bits = |c: &lindhard::ion::stopping::dataset::Comparison| {
        let mut v = vec![
            c.overall.mean_relative.to_bits(),
            c.overall.rms_relative.to_bits(),
            c.overall.max_abs_relative.to_bits(),
        ];
        v.extend(c.residuals.iter().map(|r| r.relative.to_bits()));
        v.extend(c.systems.iter().map(|s| s.summary.rms_relative.to_bits()));
        v
    };
    let want = bits(&c0);
    assert_eq!(bits(&run()), want);
    std::thread::scope(|scope| {
        let hs: Vec<_> = (0..4).map(|_| scope.spawn(|| bits(&run()))).collect();
        for h in hs {
            assert_eq!(h.join().unwrap(), want);
        }
    });
}

#[test]
fn accepts_bom_crlf_comments_blank_lines_quotes_and_repeats() {
    let text = "\u{feff}# format: lindhard-stopping-data/1\r\n\
# name: synthetic edge cases\r\n\
# A free comment line in the header.\r\n\
# source_note: header keys other than format/name/description are kept\r\n\
\r\n\
ion_z, ion_mass_amu, target, energy_ev, s_e_ev_1e15_cm2, sigma_ev_1e15_cm2, method, citation\r\n\
# a comment inside the data block\r\n\
5,10.0129,Si,1.0e4,12.0,0,synthetic,\"Made-up, test only\"\r\n\
\r\n\
5,10.0129,Si,1.0e4,12.5,0.5,synthetic,\"Made-up \"\"repeat\"\", test only\"\r\n\
5,,Al2O3,2.0e4,9.0,0.2,synthetic,database X record 17 (made-up)\r\n\
5,,Si,1.0e4,12.0,0.5,synthetic,made-up\r\n\
\r\n\
\r\n";
    let ds = StoppingDataset::from_csv_str(text).unwrap();
    assert_eq!(ds.name(), "synthetic edge cases");
    assert_eq!(ds.description(), None);
    assert_eq!(ds.metadata().len(), 3);
    let p = ds.points();
    assert_eq!(p.len(), 4);
    assert_eq!(p[0].line(), 8);
    assert_eq!(p[0].citation(), "Made-up, test only");
    assert_eq!(p[1].citation(), "Made-up \"repeat\", test only");
    assert_eq!(p[1].line(), 10);
    assert_eq!(p[0].ion().mass_amu(), 10.0129);
    assert_eq!(p[0].sigma_ev_1e15_cm2(), 0.0);
    assert_eq!(p[0].energy_ev(), p[1].energy_ev());
    assert_eq!(p[2].target().label(), "Al2O3");
    assert!(matches!(p[2].target(), Target::Compound { .. }));
    let c = compare(&LindhardScharff::new(), &NoCorrection, &ds);
    assert_eq!(
        c.systems.len(),
        3,
        "different ion masses are different systems"
    );
    assert_eq!(c.systems[0].summary.n, 2);
}

/// Load `rows` under the standard header and return the error's line and
/// reason.
fn reject(header: &str, columns: &str, rows: &str) -> (usize, String) {
    match StoppingDataset::from_csv_str(&format!("{header}{columns}{rows}")) {
        Err(StoppingError::InvalidDataset { line, reason }) => (line, reason),
        other => panic!("expected InvalidDataset, got {other:?}"),
    }
}

fn reject_row(row: &str) -> (usize, String) {
    let good = "5,,Si,1.0e4,12.0,0.1,synthetic,made-up\n";
    reject(HEADER, COLUMNS, &format!("{good}{row}\n"))
}

#[test]
fn rejects_missing_or_empty_citation_and_method() {
    for (row, what) in [
        ("5,,Si,1e4,12,0.1,synthetic,", "citation"),
        ("5,,Si,1e4,12,0.1,synthetic,\"  \"", "citation"),
        ("5,,Si,1e4,12,0.1,synthetic,\"\"", "citation"),
        ("5,,Si,1e4,12,0.1,,made-up", "method"),
    ] {
        let (line, reason) = reject_row(row);
        assert_eq!(line, 6, "{row}");
        assert!(reason.contains(&format!("empty {what}")), "{row}: {reason}");
    }
    // A missing citation field altogether is a field-count error.
    let (line, reason) = reject_row("5,,Si,1e4,12,0.1,synthetic");
    assert_eq!(line, 6);
    assert!(reason.contains("expected 8 fields"), "{reason}");
}

#[test]
fn rejects_wrong_or_missing_unit_columns() {
    for cols in [
        "ion_z,ion_mass_amu,target,energy_kev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method,citation\n",
        "ion_z,ion_mass_amu,target,energy_ev,s_e_mev_cm2_mg,sigma_ev_1e15_cm2,method,citation\n",
        "ion_z,ion_mass_amu,target,energy_ev,s_e,sigma_ev_1e15_cm2,method,citation\n",
        "ion_z,ion_mass_amu,target,energy_ev,sigma_ev_1e15_cm2,method,citation\n",
        "ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method\n",
        "ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method,citation,x\n",
    ] {
        let (line, reason) = reject(HEADER, cols, "5,,Si,1e4,12,0.1,synthetic,made-up\n");
        assert_eq!(line, 4, "{cols}: {reason}");
        assert!(reason.contains("column"), "{reason}");
    }
}

#[test]
fn rejects_bad_energy_stopping_and_uncertainty() {
    for (row, col) in [
        ("5,,Si,0,12,0.1,synthetic,made-up", "energy_ev"),
        ("5,,Si,-1e4,12,0.1,synthetic,made-up", "energy_ev"),
        ("5,,Si,NaN,12,0.1,synthetic,made-up", "energy_ev"),
        ("5,,Si,inf,12,0.1,synthetic,made-up", "energy_ev"),
        ("5,,Si,1e4,0,0.1,synthetic,made-up", "s_e_ev_1e15_cm2"),
        ("5,,Si,1e4,-3,0.1,synthetic,made-up", "s_e_ev_1e15_cm2"),
        ("5,,Si,1e4,inf,0.1,synthetic,made-up", "s_e_ev_1e15_cm2"),
        ("5,,Si,1e4,12,-0.1,synthetic,made-up", "sigma_ev_1e15_cm2"),
        ("5,,Si,1e4,12,NaN,synthetic,made-up", "sigma_ev_1e15_cm2"),
        ("5,,Si,1e4,12,,synthetic,made-up", "sigma_ev_1e15_cm2"),
    ] {
        let (line, reason) = reject_row(row);
        assert_eq!(line, 6, "{row}");
        assert!(reason.starts_with(col), "{row}: {reason}");
    }
}

#[test]
fn rejects_unknown_ion_or_target_and_bad_mass() {
    for (row, needle) in [
        ("0,,Si,1e4,12,0.1,synthetic,made-up", "unknown ion"),
        ("93,,Si,1e4,12,0.1,synthetic,made-up", "unknown ion"),
        ("300,,Si,1e4,12,0.1,synthetic,made-up", "ion_z"),
        ("B,,Si,1e4,12,0.1,synthetic,made-up", "ion_z"),
        ("5,,Xx,1e4,12,0.1,synthetic,made-up", "unknown element"),
        ("5,,si,1e4,12,0.1,synthetic,made-up", "element symbol"),
        ("5,,,1e4,12,0.1,synthetic,made-up", "empty target"),
        ("5,,SiSi,1e4,12,0.1,synthetic,made-up", "more than once"),
        ("5,-1,Si,1e4,12,0.1,synthetic,made-up", "mass"),
    ] {
        let (line, reason) = reject_row(row);
        assert_eq!(line, 6, "{row}");
        assert!(reason.contains(needle), "{row}: {reason}");
    }
}

#[test]
fn rejects_malformed_rows() {
    for row in [
        "5,,Si,1e4,12,0.1,synthetic,\"unterminated",
        "5,,Si,1e4,12,0.1,synthetic,\"a\"b",
        "5,,Si,1e4,12,0.1,synthetic,a\"b",
        "5,,Si,1e4,twelve,0.1,synthetic,made-up",
        "5,,Si,1e4,12,0.1,synthetic,made-up,extra",
    ] {
        let (line, _) = reject_row(row);
        assert_eq!(line, 6, "{row}");
    }
}

#[test]
fn rejects_bad_headers() {
    let (line, reason) = reject("# name: x\n", COLUMNS, "");
    assert_eq!(line, 2);
    assert!(reason.contains("format"), "{reason}");
    let (line, reason) = reject(
        "# format: lindhard-stopping-data/2\n# name: x\n",
        COLUMNS,
        "",
    );
    assert_eq!(line, 3);
    assert!(reason.contains("unsupported format"), "{reason}");
    let (_, reason) = reject("# format: lindhard-stopping-data/1\n# name:\n", COLUMNS, "");
    assert!(reason.contains("name"), "{reason}");
    let (line, reason) = reject(
        "# format: lindhard-stopping-data/1\n# name: a\n# name: b\n",
        COLUMNS,
        "",
    );
    assert_eq!(line, 3);
    assert!(reason.contains("twice"), "{reason}");
    for text in [
        HEADER.to_string(),
        format!("{HEADER}{COLUMNS}\n# only a comment\n"),
    ] {
        assert!(matches!(
            StoppingDataset::from_csv_str(&text),
            Err(StoppingError::InvalidTable(_))
        ));
    }
}

#[test]
fn loads_from_file() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("lindhard-stopping-dataset");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("synthetic.csv");
    std::fs::write(&path, ls_dataset(1.0)).unwrap();
    let ds = StoppingDataset::from_csv_file(&path).unwrap();
    assert_eq!(ds.points().len(), 27);
    assert!(matches!(
        StoppingDataset::from_csv_file(dir.join("missing.csv")),
        Err(StoppingError::InvalidTable(_))
    ));
}
