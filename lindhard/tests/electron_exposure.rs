//! `tally::exposure`: the compound-Poisson mean, variance and threshold
//! probabilities computed from the per-electron histogram agree with direct
//! simulation of finite exposures (`N ~ Poisson(λ)` independent histories
//! summed, many exposures), and the exposure output is bit-identical across
//! thread counts. The cross-section tables are computed here from simple
//! formulas (constant mean free paths, isotropic angles), as in
//! `tests/psf.rs`; no measured or tabulated data.

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::{LayerTables, Primary, Transport, TransportConfig};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::rng::{stream, uniform};
use lindhard::tally::{
    BeamProfile, Binning, CompoundRequest, ElectronReport, ElectronTallyConfig, ExposureConfig,
    ExposureReport, ExposureVolume, FullElectronTally, VolumeShape,
};

const NM: f64 = 1e-9;
const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];
const ENERGY: f64 = 2_000.0;
const LAMBDA: f64 = 3.0;

fn table(axis: SamplingAxis, rate: f64, row: impl Fn(f64, f64) -> f64) -> CrossSectionTable {
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "computed in tests/electron_exposure.rs".into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: vec![rate; GRID.len()],
        probability: PROB.to_vec(),
        quantiles: GRID
            .iter()
            .map(|&e| PROB.iter().map(|&p| row(e, p)).collect())
            .collect(),
    })
    .unwrap()
}

fn transport() -> Transport {
    Transport::new(
        Stack::semi_infinite(Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()),
        vec![LayerTables {
            elastic: table(SamplingAxis::ElasticPolarAngle, 1.0 / (2.0 * NM), |_, p| {
                (1.0 - 2.0 * p).acos()
            }),
            inelastic: table(
                SamplingAxis::InelasticEnergyLoss,
                1.0 / (3.0 * NM),
                |e, p| 0.3 * e * p,
            ),
        }],
        TransportConfig::new(10.0),
    )
    .unwrap()
}

fn volume(name: &str, offset_y_nm: f64) -> ExposureVolume {
    ExposureVolume {
        name: name.into(),
        shape: VolumeShape::Cylinder { radius_m: 3.0 * NM },
        depth_lo_m: 0.0,
        depth_hi_m: 4.0 * NM,
        offset_y_m: offset_y_nm * NM,
        offset_z_m: 0.0,
    }
}

fn config(beam: BeamProfile, seed: u64) -> ElectronTallyConfig {
    let mut c = ElectronTallyConfig::new(
        Binning::new(0.0, ENERGY, 20).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    c.exposure = Some(ExposureConfig {
        volumes: vec![volume("target", 0.0), volume("neighbour", 6.0)],
        beam,
        seed,
        bins: 400,
        x_max_ev: ENERGY,
        lambda: LAMBDA,
        quantiles: vec![0.1, 0.5, 0.9],
        below_ev: vec![40.0],
        above_ev: vec![100.0],
        tail_tolerance: 1e-12,
    });
    c
}

fn run(beam: BeamProfile, threads: usize, seed: u64, n: u64) -> ElectronReport {
    let t = transport();
    let proto = FullElectronTally::new(&t, config(beam, 7)).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(|| {
            t.run(seed, n, 64, &Primary::normal(ENERGY), || proto.clone())
                .unwrap()
                .tally
                .report()
        })
}

fn exposure(r: &ElectronReport) -> &ExposureReport {
    r.deposition.exposure.as_ref().unwrap()
}

fn poisson(rng: &mut impl rand_core::Rng, mean: f64) -> u64 {
    let (l, mut k, mut p) = ((-mean).exp(), 0, 1.0);
    loop {
        p *= uniform(rng);
        if p <= l {
            return k;
        }
        k += 1;
    }
}

/// Direct finite exposures: `N ~ Poisson(λ)` fresh histories each, summed in
/// the target volume. The sum is read back as the histogram's `Σ X`.
fn direct_exposures(m: usize) -> Vec<f64> {
    (0..m)
        .map(|e| {
            let n = poisson(&mut stream(0xD1CE, e as u64), LAMBDA);
            if n == 0 {
                return 0.0;
            }
            let r = run(BeamProfile::Point, 2, 10_000 + e as u64, n);
            exposure(&r).volumes[0].histogram.sum_x
        })
        .collect()
}

#[test]
fn compound_poisson_matches_direct_simulation_of_finite_exposures() {
    let big = run(BeamProfile::Point, 2, 0xE5, 20_000);
    let v = &big.deposition.exposure.as_ref().unwrap().volumes[0];
    let c = v.compound.as_ref().unwrap();
    assert!(v.zero_fraction > 0.0 && v.zero_fraction < 1.0);

    let m = 3_000;
    let s = direct_exposures(m);
    let mf = m as f64;
    let mean = s.iter().sum::<f64>() / mf;
    let sq: Vec<f64> = s.iter().map(|x| (x - mean).powi(2)).collect();
    let var = sq.iter().sum::<f64>() / mf;
    let se_mean = (var / mf).sqrt();
    assert!(
        (mean - c.mean_ev).abs() < 5.0 * (se_mean.powi(2) + c.mean_std_err_ev.powi(2)).sqrt(),
        "mean {mean} vs {} (se {se_mean}, {})",
        c.mean_ev,
        c.mean_std_err_ev
    );
    let sq_mean = var;
    let se_var = (sq.iter().map(|x| (x - sq_mean).powi(2)).sum::<f64>() / mf / mf).sqrt();
    assert!(
        (var - c.variance_ev2).abs()
            < 5.0 * (se_var.powi(2) + c.variance_std_err_ev2.powi(2)).sqrt(),
        "variance {var} vs {} (se {se_var}, {})",
        c.variance_ev2,
        c.variance_std_err_ev2
    );

    // Threshold probabilities and quantiles lie in the binning bracket.
    let d = c.distribution.as_ref().unwrap();
    let tol = 5.0 * (0.25 / mf).sqrt() + 0.01;
    let p_above = s.iter().filter(|&&x| x > 100.0).count() as f64 / mf;
    let b = d.above[0].probability;
    assert!(
        p_above > b.low - tol && p_above < b.high + tol,
        "{p_above} {b:?}"
    );
    let p_below = s.iter().filter(|&&x| x < 40.0).count() as f64 / mf;
    let b = d.below[0].probability;
    assert!(
        p_below > b.high - tol && p_below < b.low + tol,
        "{p_below} {b:?}"
    );
    let p_zero = s.iter().filter(|&&x| x == 0.0).count() as f64 / mf;
    assert!((p_zero - c.p_zero).abs() < tol, "{p_zero} {}", c.p_zero);
    let mut sorted = s.clone();
    sorted.sort_by(f64::total_cmp);
    for q in &d.quantiles {
        let e = sorted[((q.p * mf) as usize).min(m - 1)];
        let slack = 2.0 * big_bin_width(&big) + 0.05 * e;
        assert!(
            e >= q.energy_ev.low - slack && e <= q.energy_ev.high + slack,
            "p {} empirical {e} vs {:?}",
            q.p,
            q.energy_ev
        );
    }
}

fn big_bin_width(r: &ElectronReport) -> f64 {
    exposure(r).volumes[0].histogram.bin_width_ev
}

#[test]
fn exposure_output_is_bit_identical_for_1_2_and_8_threads() {
    let beam = BeamProfile::Gaussian { sigma_m: 2.0 * NM };
    let r1 = run(beam, 1, 0x95F, 1_500);
    let e1 = exposure(&r1);
    assert_eq!(e1.volumes.len(), 2);
    assert!(e1.volumes[1].compound.is_some());
    let j1 = serde_json::to_string(e1).unwrap();
    for n in [2, 8] {
        let r = run(beam, n, 0x95F, 1_500);
        assert_eq!(r1, r, "report differs on {n} threads");
        assert_eq!(j1, serde_json::to_string(exposure(&r)).unwrap());
        assert_eq!(e1.histogram_csv(), exposure(&r).histogram_csv());
    }
}

#[test]
fn the_neighbour_receives_less_than_the_target_and_other_doses_need_no_rerun() {
    let r = run(BeamProfile::Point, 2, 0xAB, 6_000);
    let e = exposure(&r);
    let (t, n) = (&e.volumes[0], &e.volumes[1]);
    assert!(n.mean_ev < t.mean_ev);
    // The stored histogram gives any other dose.
    let req = CompoundRequest {
        quantiles: vec![0.5],
        below_ev: vec![],
        above_ev: vec![100.0],
        tail_tolerance: 1e-12,
    };
    let lo = t.histogram.compound_poisson(0.5, &req).unwrap();
    let hi = t.histogram.compound_poisson(5.0, &req).unwrap();
    assert!(
        lo.distribution.unwrap().above[0].probability.central
            < hi.distribution.unwrap().above[0].probability.central
    );
    assert!((hi.mean_ev / lo.mean_ev - 10.0).abs() < 1e-12);
}
