//! Finite electron beam (`BeamSpread`): the sampled start positions and
//! directions follow the stated laws, a default beam leaves a run unchanged,
//! and the starts are identical at any thread count. The cross-section tables
//! are synthetic and empty (no interaction), since only the history start is
//! under test.

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::{
    BeamProfile, BeamSpread, ElectronState, ElectronTally, LayerTables, Primary, Transport,
    TransportConfig,
};
use lindhard::geometry::Stack;
use lindhard::ion::crystal::Divergence;
use lindhard::material::Material;

mod common;
use common::with_threads;

const NM: f64 = 1e-9;
const GRID: [f64; 3] = [1.0, 100.0, 10_000.0];

fn table(axis: SamplingAxis) -> CrossSectionTable {
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "computed in tests/electron_beam_profile.rs".into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: vec![0.0; GRID.len()],
        probability: vec![0.0, 1.0],
        quantiles: vec![vec![]; GRID.len()],
    })
    .unwrap()
}

fn transport() -> Transport {
    let m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    Transport::new(
        Stack::new(vec![(m, 5.0 * NM)], None).unwrap(),
        vec![LayerTables {
            elastic: table(SamplingAxis::ElasticPolarAngle),
            inelastic: table(SamplingAxis::InelasticEnergyLoss),
        }],
        TransportConfig::new(10.0),
    )
    .unwrap()
}

/// Records every primary's start state, in history order.
#[derive(Clone, Default)]
struct Starts(Vec<([f64; 3], [f64; 3])>);

impl ElectronTally for Starts {
    fn begin_history(&mut self, _index: u64, start: &ElectronState) {
        self.0.push((start.pos, start.dir));
    }
    fn merge(&mut self, other: Self) {
        self.0.extend(other.0);
    }
}

fn starts(primary: &Primary, beam: &BeamSpread, n: u64, seed: u64) -> Vec<([f64; 3], [f64; 3])> {
    transport()
        .run_beam(seed, n, 64, primary, beam, Starts::default)
        .unwrap()
        .tally
        .0
}

fn mean_var(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    (
        m,
        v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0),
    )
}

#[test]
fn default_beam_is_the_plain_run() {
    let p = Primary::normal(1_000.0);
    let t = transport();
    let plain = t.run(5, 300, 64, &p, Starts::default).unwrap();
    let beam = t
        .run_beam(5, 300, 64, &p, &BeamSpread::default(), Starts::default)
        .unwrap();
    assert_eq!(plain.tally.0, beam.tally.0);
    assert!(plain.tally.0.iter().all(|(pos, _)| *pos == [0.0; 3]));
    assert!(beam.metadata.beam.is_point());
    assert_eq!(plain.metadata, beam.metadata);
}

#[test]
fn point_beam_at_a_landing_point_draws_nothing() {
    let p = Primary::normal(1_000.0);
    let b = BeamSpread {
        landing_m: [3.0 * NM, -2.0 * NM],
        ..BeamSpread::default()
    };
    for (pos, dir) in starts(&p, &b, 50, 1) {
        assert_eq!(pos, [0.0, 3.0 * NM, -2.0 * NM]);
        assert_eq!(dir, [1.0, 0.0, 0.0]);
    }
}

#[test]
fn gaussian_profile_has_the_stated_mean_and_variance() {
    let sigma = 5.0 * NM;
    let (y0, z0) = (2.0 * NM, -3.0 * NM);
    let b = BeamSpread {
        landing_m: [y0, z0],
        profile: BeamProfile::Gaussian { sigma_m: sigma },
        divergence: Divergence::None,
    };
    let n = 40_000u64;
    let s = starts(&Primary::normal(1_000.0), &b, n, 7);
    let ys: Vec<f64> = s.iter().map(|(p, _)| p[1]).collect();
    let zs: Vec<f64> = s.iter().map(|(p, _)| p[2]).collect();
    for (v, m0) in [(&ys, y0), (&zs, z0)] {
        let (m, var) = mean_var(v);
        assert!((m - m0).abs() < 4.0 * sigma / (n as f64).sqrt(), "mean {m}");
        assert!((var / (sigma * sigma) - 1.0).abs() < 0.03, "var {var}");
    }
    // y and z are uncorrelated.
    let cov: f64 = ys
        .iter()
        .zip(&zs)
        .map(|(y, z)| (y - y0) * (z - z0))
        .sum::<f64>()
        / n as f64;
    assert!(cov.abs() < 0.03 * sigma * sigma, "cov {cov}");
}

#[test]
fn disc_is_uniform_in_area() {
    let r0 = 10.0 * NM;
    let b = BeamSpread {
        landing_m: [NM, NM],
        profile: BeamProfile::Disc { radius_m: r0 },
        divergence: Divergence::None,
    };
    let n = 40_000u64;
    let s = starts(&Primary::normal(1_000.0), &b, n, 11);
    let r: Vec<f64> = s
        .iter()
        .map(|(p, _)| ((p[1] - NM).powi(2) + (p[2] - NM).powi(2)).sqrt())
        .collect();
    assert!(r.iter().all(|&r| r <= r0));
    // Equal-area annuli: the fraction inside r0*sqrt(f) is f.
    for f in [0.25, 0.5, 0.75] {
        let inside = r.iter().filter(|&&r| r <= r0 * f64::sqrt(f)).count() as f64 / n as f64;
        assert!((inside - f).abs() < 0.015, "f={f}: {inside}");
    }
    // Uniform in azimuth: each quadrant holds a quarter.
    let q = s.iter().filter(|(p, _)| p[1] > NM && p[2] > NM).count() as f64 / n as f64;
    assert!((q - 0.25).abs() < 0.015, "quadrant {q}");
}

#[test]
fn rectangle_is_uniform_within_its_widths() {
    let (wy, wz) = (8.0 * NM, 4.0 * NM);
    let b = BeamSpread {
        landing_m: [0.0, 0.0],
        profile: BeamProfile::Rectangle {
            width_y_m: wy,
            width_z_m: wz,
        },
        divergence: Divergence::None,
    };
    let n = 40_000u64;
    let s = starts(&Primary::normal(1_000.0), &b, n, 13);
    assert!(s
        .iter()
        .all(|(p, _)| p[1].abs() <= 0.5 * wy && p[2].abs() <= 0.5 * wz));
    let (_, vy) = mean_var(&s.iter().map(|(p, _)| p[1]).collect::<Vec<_>>());
    assert!((vy / (wy * wy / 12.0) - 1.0).abs() < 0.03, "var {vy}");
    let left = s.iter().filter(|(p, _)| p[1] < -0.25 * wy).count() as f64 / n as f64;
    assert!((left - 0.25).abs() < 0.015, "left quarter {left}");
}

#[test]
fn divergence_cone_stays_inside_and_defaults_to_central() {
    let half = 2.0_f64.to_radians();
    let b = BeamSpread {
        divergence: Divergence::UniformCone {
            half_angle_rad: half,
        },
        ..BeamSpread::default()
    };
    let s = starts(&Primary::normal(1_000.0), &b, 2_000, 3);
    for (_, d) in &s {
        assert!((d.iter().map(|c| c * c).sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(d[0] >= half.cos() - 1e-12);
    }
    assert!(s.iter().any(|(_, d)| d[1] != 0.0 || d[2] != 0.0));
}

#[test]
fn divergence_is_conditioned_inward_at_grazing_tilt() {
    let tilt = 89.9_f64.to_radians();
    let p = Primary {
        energy_ev: 1_000.0,
        direction: [tilt.cos(), tilt.sin(), 0.0],
    };
    let b = BeamSpread {
        divergence: Divergence::Gaussian {
            sigma_rad: 5.0_f64.to_radians(),
        },
        ..BeamSpread::default()
    };
    let s = starts(&p, &b, 2_000, 9);
    assert!(s.iter().all(|(_, d)| d[0] > 0.0));
}

#[test]
fn starts_are_identical_at_any_thread_count() {
    let b = BeamSpread {
        landing_m: [NM, 0.0],
        profile: BeamProfile::Gaussian { sigma_m: 3.0 * NM },
        divergence: Divergence::Gaussian { sigma_rad: 0.02 },
    };
    let p = Primary::normal(1_000.0);
    let on = |threads: usize| with_threads(threads, || starts(&p, &b, 1_000, 21));
    assert_eq!(on(1), on(2));
}

#[test]
fn invalid_beams_are_refused() {
    let p = Primary::normal(1_000.0);
    let t = transport();
    let bad = [
        BeamProfile::Gaussian { sigma_m: 0.0 },
        BeamProfile::Gaussian { sigma_m: -1.0 },
        BeamProfile::Disc { radius_m: f64::NAN },
        BeamProfile::Rectangle {
            width_y_m: 1e-9,
            width_z_m: 0.0,
        },
    ];
    for profile in bad {
        let b = BeamSpread {
            profile,
            ..BeamSpread::default()
        };
        assert!(t.run_beam(1, 4, 2, &p, &b, Starts::default).is_err());
    }
    let b = BeamSpread {
        landing_m: [f64::INFINITY, 0.0],
        ..BeamSpread::default()
    };
    assert!(t.run_beam(1, 4, 2, &p, &b, Starts::default).is_err());
}
