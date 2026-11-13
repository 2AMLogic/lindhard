//! Primary-energy spectrum in transport (issue #641): the opt-in
//! `Bca::with_spectrum` gives each history an energy drawn from a weighted
//! list of lines, on a copy of the history's own stream at
//! `ENERGY_STREAM_WORD`. Checks: the draw equals an independent sampler on
//! the documented segment; line frequencies follow the weights; one- and
//! multi-thread runs agree; the energy budget closes for every history; a
//! one-line spectrum reproduces the plain run bit for bit (the draw sits in a
//! disjoint stream segment, so no transport draw moves); a self-ion beam
//! (Si into Si, a recoil start) works; invalid spectra fail.

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, EnergyBudget, EnergySpectrum, ENERGY_STREAM_WORD,
};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::rng::{stream, uniform};

mod common;
use common::with_threads;
use common::zbl_si_table as table;

#[derive(Debug, Clone, PartialEq, Default)]
struct Record {
    index: u64,
    incident: f64,
    residual: f64,
    lattice: f64,
    backscattered: f64,
}

#[derive(Default)]
struct Recorder(Vec<Record>);

impl BcaTally for Recorder {
    fn end_history(&mut self, index: u64, b: &EnergyBudget) {
        self.0.push(Record {
            index,
            incident: b.incident,
            residual: b.residual(),
            lattice: b.lattice,
            backscattered: b.backscattered,
        });
    }
    fn merge(&mut self, other: Self) {
        self.0.extend(other.0);
    }
}

fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    m
}

fn config() -> BcaConfig {
    let mut c = BcaConfig::new(5.0, 2.0);
    c.seed = 41;
    c
}

/// Si into Si: the beam species is a target element.
fn run(spectrum: Option<&[(f64, f64)]>, threads: usize, n: u64) -> Vec<Record> {
    let stack = Stack::new(vec![], Some(si())).unwrap();
    let ls = LindhardScharff::new();
    let e_max = spectrum.map_or(1.0e3, |s| s.iter().map(|l| l.0).fold(0.0, f64::max));
    let beam = Beam::normal(Ion::new(14).unwrap(), e_max, n);
    let mut bca = Bca::new(beam, &stack, config(), &ls, table()).unwrap();
    if let Some(s) = spectrum {
        bca = bca.with_spectrum(EnergySpectrum::new(s).unwrap());
    }
    let mut recs = with_threads(threads, || bca.run(Recorder::default))
        .unwrap()
        .0;
    recs.sort_by_key(|r| r.index);
    recs
}

const LINES: [(f64, f64); 2] = [(1.0e3, 3.0), (2.0e3, 1.0)];

#[test]
fn energy_draw_matches_an_independent_sampler_on_the_documented_segment() {
    assert_eq!(ENERGY_STREAM_WORD, 1u128 << 64);
    let recs = run(Some(&LINES), 1, 200);
    for r in &recs {
        let mut rng = stream(41, r.index);
        rng.set_word_pos(ENERGY_STREAM_WORD);
        // Weights 3:1 normalise to 0.75 / 0.25.
        let e = if uniform(&mut rng) < 0.75 {
            1.0e3
        } else {
            2.0e3
        };
        assert_eq!(r.incident, e, "history {}", r.index);
    }
}

#[test]
fn line_frequencies_follow_the_weights_and_swap_with_them() {
    let n = 4000u64;
    let count = |s: &[(f64, f64)]| {
        run(Some(s), 2, n)
            .iter()
            .filter(|r| r.incident == 1.0e3)
            .count() as f64
    };
    // 3:1 over n draws: mean 3000, sd sqrt(n p q) = 27.4; 6 sd of slack.
    let sd = (n as f64 * 0.75 * 0.25).sqrt();
    let low = count(&LINES);
    assert!((low - 0.75 * n as f64).abs() < 6.0 * sd, "got {low}");
    let swapped = count(&[(1.0e3, 1.0), (2.0e3, 3.0)]);
    assert!(
        (swapped - 0.25 * n as f64).abs() < 6.0 * sd,
        "got {swapped}"
    );
}

#[test]
fn spectrum_runs_are_thread_count_independent() {
    assert_eq!(run(Some(&LINES), 1, 300), run(Some(&LINES), 3, 300));
}

#[test]
fn energy_budget_closes_for_every_history_of_a_self_ion_spectrum_run() {
    for r in run(Some(&LINES), 2, 300) {
        assert!(r.incident == 1.0e3 || r.incident == 2.0e3);
        assert!(
            r.residual.abs() <= 1e-9 * r.incident,
            "history {}: residual {} eV of {} eV",
            r.index,
            r.residual,
            r.incident
        );
    }
}

#[test]
fn a_one_line_spectrum_reproduces_the_plain_run_bit_for_bit() {
    // The draw consumes its own stream segment, so no transport draw moves.
    assert_eq!(run(None, 1, 200), run(Some(&[(1.0e3, 5.0)]), 1, 200));
}

#[test]
fn invalid_spectra_are_rejected() {
    assert!(EnergySpectrum::new(&[]).is_err());
    assert!(EnergySpectrum::new(&[(0.0, 1.0)]).is_err());
    assert!(EnergySpectrum::new(&[(f64::NAN, 1.0)]).is_err());
    assert!(EnergySpectrum::new(&[(1.0e3, 0.0)]).is_err());
    assert!(EnergySpectrum::new(&[(1.0e3, -1.0)]).is_err());
    assert!(EnergySpectrum::new(&[(1.0e3, f64::INFINITY)]).is_err());
    let ok = EnergySpectrum::new(&LINES).unwrap();
    let w: f64 = ok.lines().iter().map(|l| l.1).sum();
    assert!((w - 1.0).abs() < 1e-15);
    assert_eq!((ok.min_energy_ev(), ok.max_energy_ev()), (1.0e3, 2.0e3));
}
