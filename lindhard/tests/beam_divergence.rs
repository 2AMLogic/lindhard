//! Beam divergence in transport (issue #285): the opt-in
//! `Bca::with_divergence` samples each primary's direction about the nominal
//! beam direction on its own indexed stream, conditioned on pointing into the
//! target. Checks: no divergence is bit-identical to the engine without the
//! setter (including a mixed amorphous/crystal target); the recorded initial
//! directions are unit vectors that equal an independent rejection sampler on
//! the documented stream segment; narrow normal-incidence and near-grazing
//! nominal directions; one- and multi-thread reports agree; invalid settings
//! fail.

use std::f64::consts::PI;
use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, BcaError, Beam, CrystalTarget, DIVERGENCE_STREAM_WORD};
use lindhard::ion::crystal::{Divergence, Lattice, Orientation};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::rng::stream;
use lindhard::tally::{Binning, IonTally, IonTallyConfig};

const NM: f64 = 1e-9;

fn table() -> &'static ScatteringTable {
    static T: OnceLock<ScatteringTable> = OnceLock::new();
    T.get_or_init(|| {
        ScatteringTable::build(
            &Potential::new(Screening::ZblUniversal, 14.0, 14.0),
            &TableSpec {
                eps_min: 1e-6,
                eps_max: 1e4,
                beta_min: 1e-5,
                beta_max: 1e2,
                per_decade: 16,
            },
        )
    })
}

fn elemental(z: u8) -> Material {
    let mut m = Material::from_atom_fractions(&[(z, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(z, 15.0).unwrap();
    if m.surface_binding_energy_ev(z).is_err() {
        m.set_surface_binding_energy_ev(z, 2.0).unwrap();
    }
    m
}

/// A 20 nm amorphous Si screen over a crystalline Si substrate.
fn mixed_stack() -> Stack {
    Stack::new(vec![(elemental(14), 20.0 * NM)], Some(elemental(14))).unwrap()
}

fn beam(tilt_deg: f64, n: u64) -> Beam {
    Beam {
        ion: Ion::new(5).unwrap(),
        energy_ev: 5.0e3,
        polar_rad: tilt_deg.to_radians(),
        azimuth_rad: 22f64.to_radians(),
        count: n,
    }
}

fn config() -> BcaConfig {
    let mut c = BcaConfig::new(5.0, 2.0);
    c.seed = 11;
    c
}

fn crystal(tilt_deg: f64) -> CrystalTarget {
    let lat = Lattice::silicon();
    let o = Orientation::new(
        &lat,
        [0, 0, 1],
        [0, 1, 0],
        tilt_deg.to_radians(),
        22f64.to_radians(),
        0.0,
    )
    .unwrap();
    CrystalTarget::new(lat, o)
}

/// Serialised report of the mixed target on `threads` workers, with an
/// optional divergence setter.
fn report(div: Option<Divergence>, threads: usize, n: u64) -> String {
    let stack = mixed_stack();
    let ls = LindhardScharff::new();
    let mut bca = Bca::new(beam(7.0, n), &stack, config(), &ls, table())
        .unwrap()
        .with_crystal(crystal(7.0), &[1])
        .unwrap();
    if let Some(d) = div {
        bca = bca.with_divergence(d).unwrap();
    }
    let species = bca.species_z();
    let tc = IonTallyConfig {
        depth: Binning::new(0.0, 200.0 * NM, 100).unwrap(),
        lateral: Binning::new(-100.0 * NM, 100.0 * NM, 50).unwrap(),
        radial: Binning::new(0.0, 100.0 * NM, 25).unwrap(),
        escape_energy: Binning::new(0.0, 5.0e4, 25).unwrap(),
        escape_polar: Binning::new(0.0, PI / 2.0, 9).unwrap(),
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let tally = pool
        .install(|| bca.run(|| IonTally::new(&stack, &species, tc).unwrap()))
        .unwrap();
    serde_json::to_string(&tally.report(false)).unwrap()
}

#[test]
fn absent_divergence_preserves_fixed_seed_results_on_a_mixed_target() {
    let plain = report(None, 1, 300);
    assert_eq!(plain, report(Some(Divergence::None), 1, 300));
}

#[test]
fn divergent_reports_are_thread_count_independent_on_a_mixed_target() {
    for div in [
        Divergence::Gaussian { sigma_rad: 0.01 },
        Divergence::UniformCone {
            half_angle_rad: 0.02,
        },
    ] {
        let one = report(Some(div), 1, 300);
        assert_eq!(one, report(Some(div), 3, 300), "{div:?}");
        assert_ne!(one, report(None, 1, 300), "{div:?} must change transport");
    }
}

/// The documented law, written independently of the engine: the divergence
/// sample on the history's stream at word 2^65, redrawn until inward.
fn reference_direction(div: &Divergence, central: [f64; 3], seed: u64, index: u64) -> [f64; 3] {
    let mut rng = stream(seed, index);
    rng.set_word_pos(DIVERGENCE_STREAM_WORD);
    assert_eq!(DIVERGENCE_STREAM_WORD, 1u128 << 65);
    loop {
        let d = div.sample_direction(central, &mut rng);
        if d[0] > 0.0 {
            return d;
        }
    }
}

fn engine<'a>(stack: &'a Stack, ls: &'a LindhardScharff, tilt_deg: f64, d: Divergence) -> Bca<'a> {
    Bca::new(beam(tilt_deg, 1), stack, config(), ls, table())
        .unwrap()
        .with_divergence(d)
        .unwrap()
}

#[test]
fn initial_directions_follow_the_inward_conditioned_law() {
    let stack = Stack::semi_infinite(elemental(14));
    let ls = LindhardScharff::new();
    // Near-grazing nominal directions with a wide spread (many raw samples
    // point outward), a moderate tilt, and normal incidence.
    let cases = [
        (89.5, Divergence::Gaussian { sigma_rad: 0.15 }),
        (
            89.0,
            Divergence::UniformCone {
                half_angle_rad: 0.17,
            },
        ),
        (7.0, Divergence::Gaussian { sigma_rad: 0.01 }),
        (
            0.0,
            Divergence::UniformCone {
                half_angle_rad: 0.1,
            },
        ),
    ];
    for (tilt, div) in cases {
        let bca = engine(&stack, &ls, tilt, div);
        let central = bca.beam().direction();
        let mut rejected_first = 0;
        for i in 0..2000 {
            let d = bca.primary_direction(i).unwrap();
            let norm = d.iter().map(|c| c * c).sum::<f64>().sqrt();
            assert!((norm - 1.0).abs() < 1e-12, "unit vector, tilt {tilt}");
            assert!(d[0] > 0.0, "inward, tilt {tilt}");
            assert_eq!(d, reference_direction(&div, central, 11, i), "tilt {tilt}");
            let mut rng = stream(11, i);
            rng.set_word_pos(DIVERGENCE_STREAM_WORD);
            if div.sample_direction(central, &mut rng)[0] <= 0.0 {
                rejected_first += 1;
            }
        }
        if tilt > 80.0 {
            assert!(rejected_first > 100, "the near-surface case must reject");
        } else {
            assert_eq!(rejected_first, 0, "no conditioning at tilt {tilt}");
        }
    }
}

#[test]
fn narrow_normal_beam_keeps_the_gaussian_rayleigh_law() {
    let stack = Stack::semi_infinite(elemental(14));
    let ls = LindhardScharff::new();
    let sigma = 0.002;
    let bca = engine(&stack, &ls, 0.0, Divergence::Gaussian { sigma_rad: sigma });
    let n = 20_000;
    let mean: f64 = (0..n)
        .map(|i| bca.primary_direction(i).unwrap()[0].clamp(-1.0, 1.0).acos())
        .sum::<f64>()
        / n as f64;
    // Rayleigh mean: sigma sqrt(pi / 2).
    let want = sigma * (PI / 2.0).sqrt();
    assert!((mean - want).abs() < 0.03 * want, "mean {mean} vs {want}");
}

#[test]
fn absent_divergence_gives_the_nominal_direction_without_draws() {
    let stack = Stack::semi_infinite(elemental(14));
    let ls = LindhardScharff::new();
    let bca = engine(&stack, &ls, 7.0, Divergence::None);
    for i in 0..10 {
        assert_eq!(bca.primary_direction(i).unwrap(), bca.beam().direction());
    }
    assert!(bca.divergence_metadata().is_none());
}

#[test]
fn invalid_divergence_is_rejected() {
    let stack = Stack::semi_infinite(elemental(14));
    let ls = LindhardScharff::new();
    for bad in [
        Divergence::Gaussian { sigma_rad: -0.1 },
        Divergence::Gaussian {
            sigma_rad: f64::NAN,
        },
        Divergence::UniformCone {
            half_angle_rad: f64::INFINITY,
        },
        Divergence::UniformCone {
            half_angle_rad: -1.0,
        },
    ] {
        let r = Bca::new(beam(0.0, 1), &stack, config(), &ls, table())
            .unwrap()
            .with_divergence(bad);
        assert!(matches!(r, Err(BcaError::InvalidBeam(_))), "{bad:?}");
    }
}

#[test]
fn metadata_records_the_resolved_spread() {
    let stack = Stack::semi_infinite(elemental(14));
    let ls = LindhardScharff::new();
    let m = engine(&stack, &ls, 0.0, Divergence::Gaussian { sigma_rad: 0.01 })
        .divergence_metadata()
        .unwrap();
    assert_eq!(m.model, "gaussian");
    assert_eq!(m.width_kind, "sigma_per_plane");
    assert_eq!(m.width_rad, 0.01);
    assert_eq!(m.incidence, "inward-conditioned");
}
