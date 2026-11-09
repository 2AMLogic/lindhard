//! Per-pair stopping correction fits: recovery of a known correction from
//! synthetic noisy points, bit-reproducibility across seeds and thread counts,
//! the out-of-span fallback, held-out reporting, and a golden test on the
//! committed Bader et al. (1956) Cu dataset.
//!
//! The synthetic points are sampled from this repository's own
//! Lindhard-Scharff model times a made-up known factor, with made-up noise;
//! they are not measurements (`docs/data-provenance.md`).

use lindhard::ion::stopping::bragg::NoCorrection;
use lindhard::ion::stopping::correction::{fit_correction, Correction, CorrectionFit, FitSettings};
use lindhard::ion::stopping::dataset::{compare, StoppingDataset};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::{to_ev_1e15_cm2, ElectronicStopping, Ion, StoppingError};
use lindhard::rng::{stream, uniform};
use rand_core::Rng;
use std::path::PathBuf;

const COLUMNS: &str =
    "ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,method,citation\n";
const CITATION: &str = "\"Lindhard-Scharff model of this repository times a made-up factor, with made-up noise (synthetic, test)\"";

/// The known correction: a quadratic in `u = log10(E / 10 keV)`, so a degree-2
/// fit can represent it exactly.
fn true_factor(energy_ev: f64) -> f64 {
    let u = (energy_ev / 1.0e4).log10();
    (0.10 + 0.25 * u - 0.15 * u * u).exp()
}

/// Approximately standard normal (sum of 12 uniforms minus 6), made up.
fn noise(rng: &mut impl Rng) -> f64 {
    (0..12).map(|_| uniform(rng)).sum::<f64>() - 6.0
}

/// 40 points from 1 keV to 100 keV: B in Si, `LS * true_factor * (1 + 2 % N(0,1))`,
/// stated uncertainty 2 %.
fn synthetic() -> StoppingDataset {
    let ls = LindhardScharff::new();
    let ion = Ion::new(5).unwrap();
    let mut rng = stream(11, 0);
    let mut text =
        format!("# format: lindhard-stopping-data/1\n# name: synthetic-correction\n{COLUMNS}");
    for i in 0..40 {
        let e = 1.0e3 * 100f64.powf(f64::from(i) / 39.0);
        let s = to_ev_1e15_cm2(ls.stopping(&ion, 14, e).unwrap()) * true_factor(e);
        let meas = s * (1.0 + 0.02 * noise(&mut rng));
        text.push_str(&format!(
            "5,,Si,{e:e},{meas:e},{:e},synthetic,{CITATION}\n",
            0.02 * meas
        ));
    }
    StoppingDataset::from_csv_str(&text).unwrap()
}

const SETTINGS: FitSettings = FitSettings {
    degree: 2,
    held_out: 6,
    replicates: 400,
    seed: 20261015,
};

fn fit(ds: &StoppingDataset, settings: FitSettings) -> CorrectionFit {
    fit_correction(&LindhardScharff::new(), ds, settings).unwrap()
}

#[test]
fn recovers_a_known_correction_within_the_bootstrap_band() {
    let ds = synthetic();
    let f = fit(&ds, SETTINGS);
    assert!(f.weighted());
    assert_eq!(f.failed_replicates(), 0);
    let (lo, hi) = f.span_ev();
    for k in 0..=10 {
        let e = lo * (hi / lo).powf(f64::from(k) / 10.0);
        let truth = true_factor(e);
        let g = f.factor(e).unwrap();
        assert!(
            (g / truth - 1.0).abs() < 0.03,
            "g({e}) = {g}, truth {truth}"
        );
        let (b_lo, b_hi) = f.band(e, 0.99).unwrap();
        assert!(
            b_lo <= truth && truth <= b_hi,
            "truth {truth} outside [{b_lo}, {b_hi}] at {e}"
        );
        assert!(b_lo < g && g < b_hi);
    }
    // A wider level gives a wider band.
    let (a, b) = f.band(1.0e4, 0.68).unwrap();
    let (c, d) = f.band(1.0e4, 0.99).unwrap();
    assert!(c < a && b < d);
}

#[test]
fn held_out_points_are_disjoint_interior_and_reported_beside_the_base() {
    let ds = synthetic();
    let f = fit(&ds, SETTINGS);
    let pts = ds.points();
    assert_eq!(f.held_out_indices().len(), 6);
    assert_eq!(
        f.train_indices().len() + f.held_out_indices().len(),
        pts.len()
    );
    for i in f.held_out_indices() {
        assert!(!f.train_indices().contains(i));
        assert!(
            *i != 0 && *i != pts.len() - 1,
            "end points must stay in training"
        );
    }
    // The base statistics are exactly `compare`'s residuals on those points.
    let cmp = compare(&LindhardScharff::new(), &NoCorrection, &ds);
    let rel: Vec<f64> = f
        .held_out_indices()
        .iter()
        .map(|&i| {
            cmp.residuals
                .iter()
                .find(|r| r.index == i)
                .unwrap()
                .relative
        })
        .collect();
    let rms = (rel.iter().map(|r| r * r).sum::<f64>() / rel.len() as f64).sqrt();
    let rep = f.held_out_report();
    assert_eq!(rep.points.len(), 6);
    assert!((rep.base.rms_relative / rms - 1.0).abs() < 1e-12);
    // The synthetic correction is large (LS is off by tens of percent), so the
    // fit must reduce the held-out residuals, and be near the 2 % noise.
    assert!(rep.reduces_rms_relative());
    assert!(rep.corrected.rms_relative < 0.05, "{:?}", rep.corrected);
    assert!(rep.corrected.rms_absolute_ev_1e15_cm2 < rep.base.rms_absolute_ev_1e15_cm2);
}

#[test]
fn a_correction_that_does_not_help_is_reported_as_such() {
    // Points equal to the base model exactly: the held-out residual of the
    // base is zero, so a fit with any numerical noise cannot reduce it.
    let ls = LindhardScharff::new();
    let ion = Ion::new(5).unwrap();
    let mut text = format!("# format: lindhard-stopping-data/1\n# name: identity\n{COLUMNS}");
    for i in 0..12 {
        let e = 1.0e3 * 100f64.powf(f64::from(i) / 11.0);
        let s = to_ev_1e15_cm2(ls.stopping(&ion, 14, e).unwrap());
        text.push_str(&format!("5,,Si,{e:e},{s:e},0,synthetic,{CITATION}\n"));
    }
    let ds = StoppingDataset::from_csv_str(&text).unwrap();
    let f = fit(
        &ds,
        FitSettings {
            replicates: 20,
            ..SETTINGS
        },
    );
    let rep = f.held_out_report();
    assert!(!f.weighted());
    assert!(rep.base.rms_relative < 1e-12);
    assert!(!rep.reduces_rms_relative());
}

#[test]
fn same_seed_is_bit_identical_whatever_the_thread_count() {
    let ds = synthetic();
    let on = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| fit(&ds, SETTINGS))
    };
    let (a, b, c) = (on(1), on(2), fit(&ds, SETTINGS));
    assert_eq!(a, b);
    assert_eq!(a, c);
    for e in [2.0e3, 1.7e4, 8.0e4] {
        let (x, y) = (a.band(e, 0.68).unwrap(), b.band(e, 0.68).unwrap());
        assert_eq!(x.0.to_bits(), y.0.to_bits());
        assert_eq!(x.1.to_bits(), y.1.to_bits());
    }
    // A different seed gives a different split or band.
    let d = fit(
        &ds,
        FitSettings {
            seed: 7,
            ..SETTINGS
        },
    );
    assert_ne!(a.band(1.0e4, 0.68), d.band(1.0e4, 0.68));
}

#[test]
fn outside_the_span_or_for_another_pair_the_base_is_returned_unchanged() {
    let ds = synthetic();
    let f = fit(&ds, SETTINGS);
    let (lo, hi) = f.span_ev();
    let base = LindhardScharff::new();
    let model = Correction::new(LindhardScharff::new(), f);
    let ion = Ion::new(5).unwrap();
    let bits = |m: &dyn ElectronicStopping, ion: &Ion, z: u8, e: f64| {
        m.stopping(ion, z, e).unwrap().to_bits()
    };
    for e in [lo * 0.999, lo * 0.5, hi * 1.001, hi * 10.0] {
        assert_eq!(
            bits(&model, &ion, 14, e),
            bits(&base, &ion, 14, e),
            "E = {e}"
        );
        assert!(model.fit().factor(e).is_none() && model.fit().band(e, 0.68).is_none());
    }
    // Another target, another ion, another isotope: base.
    for (ion, z) in [
        (Ion::new(5).unwrap(), 6u8),
        (Ion::new(15).unwrap(), 14),
        (Ion::with_mass(5, 11.009).unwrap(), 14),
    ] {
        assert_eq!(bits(&model, &ion, z, 1.0e4), bits(&base, &ion, z, 1.0e4));
    }
    // Inside the span, at both ends included, the correction acts.
    for e in [lo, 1.0e4, hi] {
        let g = model.fit().factor(e).unwrap();
        let got = model.stopping(&ion, 14, e).unwrap();
        assert_eq!(
            got.to_bits(),
            (base.stopping(&ion, 14, e).unwrap() * g).to_bits()
        );
        assert!((g - 1.0).abs() > 1e-3);
    }
    assert_eq!(model.name(), "corrected");
}

#[test]
fn unusable_input_is_an_error_not_a_silent_drop() {
    let ls = LindhardScharff::new();
    let row = |ion: &str, target: &str, e: f64| format!("{ion},,{target},{e:e},10,0,m,c\n");
    let header = format!("# format: lindhard-stopping-data/1\n# name: x\n{COLUMNS}");
    let load = |body: String| StoppingDataset::from_csv_str(&format!("{header}{body}")).unwrap();
    let s = SETTINGS;
    // Two systems.
    let ds = load(
        [
            row("5", "Si", 1e4),
            row("5", "Si", 2e4),
            row("15", "Si", 3e4),
            row("5", "Si", 4e4),
        ]
        .concat(),
    );
    assert!(matches!(
        fit_correction(&ls, &ds, FitSettings { held_out: 0, ..s }),
        Err(StoppingError::InvalidFit(_))
    ));
    // Compound.
    let ds = load(
        [
            row("5", "SiO2", 1e4),
            row("5", "SiO2", 2e4),
            row("5", "SiO2", 3e4),
            row("5", "SiO2", 4e4),
        ]
        .concat(),
    );
    assert!(matches!(
        fit_correction(&ls, &ds, FitSettings { held_out: 0, ..s }),
        Err(StoppingError::InvalidFit(_))
    ));
    // Too few points for the degree and the held-out count.
    let ds = load(
        [
            row("5", "Si", 1e4),
            row("5", "Si", 2e4),
            row("5", "Si", 3e4),
        ]
        .concat(),
    );
    assert!(matches!(
        fit_correction(&ls, &ds, s),
        Err(StoppingError::InvalidFit(_))
    ));
    // One energy only.
    let ds = load(
        [
            row("5", "Si", 1e4),
            row("5", "Si", 1e4),
            row("5", "Si", 1e4),
            row("5", "Si", 1e4),
        ]
        .concat(),
    );
    assert!(matches!(
        fit_correction(
            &ls,
            &ds,
            FitSettings {
                held_out: 0,
                degree: 1,
                ..s
            }
        ),
        Err(StoppingError::InvalidFit(_))
    ));
}

fn bader() -> Option<StoppingDataset> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../validation/data/stopping/bader1956_p_cu.csv");
    if !p.is_file() {
        // A packaged crate does not ship validation/data.
        eprintln!("{} not found; Bader golden test skipped", p.display());
        return None;
    }
    Some(StoppingDataset::from_csv_file(p).unwrap())
}

#[test]
fn golden_fit_to_the_committed_bader_cu_dataset() {
    let Some(ds) = bader() else { return };
    assert_eq!(ds.points().len(), 15);
    let settings = FitSettings {
        degree: 2,
        held_out: 3,
        replicates: 500,
        seed: 20261015,
    };
    let f = fit(&ds, settings);
    assert!(!f.weighted(), "the paper states no uncertainty");
    assert_eq!(f.span_ev(), (5.0e4, 6.0e5));
    let rep = f.held_out_report();
    println!(
        "{:?}\n{:?}\n{:?}",
        f.coefficients(),
        f.held_out_indices(),
        rep
    );
    let close =
        |got: f64, want: f64| assert!((got / want - 1.0).abs() < 1e-6, "{got} vs pinned {want}");
    for (g, w) in f.coefficients().iter().zip(GOLDEN_COEFFICIENTS) {
        close(*g, w);
    }
    assert_eq!(f.held_out_indices(), GOLDEN_HELD_OUT);
    close(rep.base.rms_relative, GOLDEN_BASE_RMS);
    close(rep.corrected.rms_relative, GOLDEN_CORRECTED_RMS);
    // Reproducible, and the held-out result is not hidden whichever way it goes.
    assert_eq!(f, fit(&ds, settings));
    assert!(rep.reduces_rms_relative());
}

const GOLDEN_COEFFICIENTS: [f64; 3] = [
    -0.7167515384988218,
    -0.7486515905937396,
    -0.10258908539243347,
];
const GOLDEN_HELD_OUT: [usize; 3] = [5, 10, 12];
const GOLDEN_BASE_RMS: f64 = 2.286480730341239;
const GOLDEN_CORRECTED_RMS: f64 = 0.0021715284976015358;
