//! Integration tests for `lindhard::tally`: moments, merging, Pearson fits and
//! the ion tally wired to the BCA engine.

use std::f64::consts::PI;
use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, Beam, SummaryTally};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::rng::{run_particles, stream, ParticleRng};
use lindhard::tally::{
    Binning, DualPearson, Histogram, IonReport, IonTally, IonTallyConfig, Moments, PearsonError,
    PearsonIv,
};
use rand_core::Rng;

const NM: f64 = 1e-9;

fn uniform(rng: &mut ParticleRng) -> f64 {
    (rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// Box-Muller standard normal.
fn normal(rng: &mut ParticleRng) -> f64 {
    let u1 = 1.0 - uniform(rng);
    let u2 = uniform(rng);
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

fn with_threads<R: Send>(n: usize, f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .unwrap()
        .install(f)
}

#[test]
fn gaussian_moments_are_recovered() {
    let (mu, sigma) = (120.0 * NM, 35.0 * NM);
    let m = run_particles(
        2024,
        1_000_000,
        4096,
        Moments::new,
        |m, rng, _| m.push(mu + sigma * normal(rng)),
        |a, b| a.merge(&b),
    );
    let s = m.summary().unwrap();
    assert_eq!(s.n, 1_000_000);
    assert!((s.mean - mu).abs() < 4.0 * s.mean_std_err, "{s:?}");
    assert!((s.std_dev - sigma).abs() < 4.0 * s.std_dev_std_err, "{s:?}");
    assert!(s.skewness.abs() < 4.0 * s.skewness_std_err, "{s:?}");
    assert!((s.kurtosis - 3.0).abs() < 4.0 * s.kurtosis_std_err, "{s:?}");
    // Standard errors are close to the normal-theory values.
    let n = 1e6f64;
    assert!((s.skewness_std_err / (6.0 / n).sqrt() - 1.0).abs() < 0.05);
    assert!((s.kurtosis_std_err / (24.0 / n).sqrt() - 1.0).abs() < 0.1);
}

/// A skewed sample (exponential plus a constant).
fn skewed_sample(n: u64) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let mut rng = stream(7, i);
            5.0 - (1.0 - uniform(&mut rng)).ln() * 2.0
        })
        .collect()
}

fn chunked(xs: &[f64], sizes: &[usize]) -> Moments {
    let mut total = Moments::new();
    let mut at = 0;
    for &s in sizes {
        let mut c = Moments::new();
        for &x in &xs[at..at + s] {
            c.push(x);
        }
        total.merge(&c);
        at += s;
    }
    assert_eq!(at, xs.len());
    total
}

#[test]
fn merging_is_split_independent_and_deterministic() {
    let xs = skewed_sample(10_000);
    let single = chunked(&xs, &[10_000]);
    let s0 = single.summary().unwrap();
    for sizes in [
        vec![5_000, 5_000],
        vec![1, 9_998, 1],
        vec![3_333, 3_333, 3_334],
        vec![64; 156].into_iter().chain([16]).collect::<Vec<_>>(),
        (0..100).map(|_| 100).collect(),
    ] {
        let m = chunked(&xs, &sizes);
        let s = m.summary().unwrap();
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-11 * b.abs().max(1.0);
        assert_eq!(s.n, s0.n);
        assert!(close(s.mean, s0.mean), "{sizes:?}");
        assert!(close(s.std_dev, s0.std_dev));
        assert!(close(s.skewness, s0.skewness));
        assert!(close(s.kurtosis, s0.kurtosis));
        assert!(close(s.kurtosis_std_err, s0.kurtosis_std_err));
        // The same chunking gives the same bits.
        assert_eq!(chunked(&xs, &sizes), m);
    }
}

#[test]
fn moments_and_histograms_are_bit_identical_across_thread_counts() {
    let binning = Binning::new(0.0, 20.0, 40).unwrap();
    let run = || {
        run_particles(
            11,
            50_000,
            256,
            || (Moments::new(), Histogram::new(binning)),
            |(m, h), rng, _| {
                let x = 5.0 - (1.0 - uniform(rng)).ln() * 2.0;
                m.push(x);
                h.fill(x);
            },
            |a, b| {
                a.0.merge(&b.0);
                a.1.merge(&b.1);
            },
        )
    };
    let r1 = with_threads(1, run);
    for n in [2, 8] {
        let r = with_threads(n, run);
        assert_eq!(r1.0.mean.to_bits(), r.0.mean.to_bits());
        for (a, b) in r1.0.power_sums.iter().zip(&r.0.power_sums) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(r1.1, r.1);
    }
}

/// `∫ x^p f(x) dx` for `p = 0..=4` by the substitution `x = λ + a tan θ`
/// and composite Simpson in `θ` over `(-π/2, π/2)`.
fn numerical_raw_moments(p: &PearsonIv) -> [f64; 5] {
    let n = 200_000;
    let h = PI / n as f64;
    let mut out = [0.0; 5];
    for i in 1..n {
        let th = -0.5 * PI + i as f64 * h;
        let x = p.lambda + p.a * th.tan();
        let jac = p.a / (th.cos() * th.cos());
        let w = if i % 2 == 1 { 4.0 } else { 2.0 } * h / 3.0;
        let f = p.pdf(x) * jac * w;
        let mut xp = 1.0;
        for o in out.iter_mut() {
            *o += f * xp;
            xp *= x;
        }
    }
    out // the end points contribute zero
}

#[test]
fn pearson_iv_round_trips_its_moments() {
    for &(mean, sd, g, b) in &[
        (100.0, 30.0, -0.5, 3.8),
        (50.0, 10.0, 0.8, 4.5),
        (0.0, 1.0, 0.0, 4.0),
        (200.0, 40.0, -1.2, 6.0),
        (10.0, 2.0, 0.3, 3.4),
    ] {
        let p = PearsonIv::from_moments(mean, sd, g, b).unwrap();
        let r = numerical_raw_moments(&p);
        assert!((r[0] - 1.0).abs() < 1e-6, "norm {}", r[0]);
        let m = r[1] / r[0];
        let c2 = r[2] / r[0] - m * m;
        let c3 = r[3] / r[0] - 3.0 * m * r[2] / r[0] + 2.0 * m.powi(3);
        let c4 = r[4] / r[0] - 4.0 * m * r[3] / r[0] + 6.0 * m * m * r[2] / r[0] - 3.0 * m.powi(4);
        let (s, gg, bb) = (c2.sqrt(), c3 / c2.powf(1.5), c4 / (c2 * c2));
        let rel = |a: f64, b: f64, scale: f64| (a - b).abs() / scale;
        assert!(rel(m, mean, sd) < 1e-6, "mean {m} vs {mean}");
        assert!(rel(s, sd, sd) < 1e-6, "sd {s} vs {sd}");
        assert!(rel(gg, g, 1.0) < 1e-6, "skew {gg} vs {g}");
        assert!(rel(bb, b, b) < 1e-6, "kurt {bb} vs {b}");
    }
}

#[test]
fn moments_outside_type_iv_are_errors() {
    // Below the type IV / VI boundary for gamma = 1 (about 4.97).
    assert!(matches!(
        PearsonIv::from_moments(0.0, 1.0, 1.0, 4.5),
        Err(PearsonError::OutsideTypeIv {
            min_kurtosis: Some(_),
            ..
        })
    ));
    // Platykurtic (type II region) and the Gaussian point.
    assert!(PearsonIv::from_moments(0.0, 1.0, 0.0, 2.5).is_err());
    assert!(PearsonIv::from_moments(0.0, 1.0, 0.0, 3.0).is_err());
    // Skewness no type IV can reach.
    assert!(matches!(
        PearsonIv::from_moments(0.0, 1.0, 6.0, 100.0),
        Err(PearsonError::OutsideTypeIv {
            min_kurtosis: None,
            ..
        })
    ));
    assert!(matches!(
        PearsonIv::from_moments(0.0, -1.0, 0.0, 4.0),
        Err(PearsonError::InvalidMoments { .. })
    ));
    assert!(PearsonIv::from_moments(f64::NAN, 1.0, 0.0, 4.0).is_err());
}

/// Expected (rounded) counts of `n` samples of `pdf` over `binning`.
fn expected_histogram(binning: Binning, n: f64, pdf: impl Fn(f64) -> f64) -> Histogram {
    let mut h = Histogram::new(binning);
    let mut inside = 0;
    for i in 0..binning.bins {
        // Fine midpoint rule per bin.
        let (lo, hi) = (binning.edge(i), binning.edge(i + 1));
        let k = 64;
        let mass: f64 = (0..k)
            .map(|j| pdf(lo + (j as f64 + 0.5) * (hi - lo) / k as f64))
            .sum::<f64>()
            * (hi - lo)
            / k as f64;
        h.counts[i] = (n * mass).round() as u64;
        inside += h.counts[i];
    }
    h.overflow = (n as u64).saturating_sub(inside);
    h
}

#[test]
fn dual_pearson_fit_describes_a_two_component_profile() {
    let head = PearsonIv::from_moments(100.0, 25.0, -0.2, 3.5).unwrap();
    let tail = PearsonIv::from_moments(190.0, 45.0, 0.4, 4.5).unwrap();
    let truth = DualPearson {
        head_fraction: 0.8,
        head,
        tail,
    };
    let binning = Binning::new(0.0, 400.0, 100).unwrap();
    let hist = expected_histogram(binning, 1e6, |x| truth.pdf(x));
    let fit = DualPearson::fit(&hist).unwrap();
    let again = DualPearson::fit(&hist).unwrap();
    assert_eq!(fit, again, "the fit is deterministic");
    // A bimodal mixture's overall moments may lie outside the type IV
    // region; when they do not, the dual fit must beat the single Pearson.
    if let Some(single) = fit.single_chi_square {
        assert!(
            fit.chi_square < 0.1 * single,
            "{} vs single {single}",
            fit.chi_square
        );
    }
    // Expected counts, so the optimum is the truth up to rounding.
    assert!(fit.chi_square < 1.0, "{fit:?}");
    assert!((fit.profile.head_fraction - 0.8).abs() < 0.02, "{fit:?}");
    // The fitted density matches the true one bin by bin.
    let peak = (0..binning.bins)
        .map(|i| truth.pdf(binning.center(i)))
        .fold(0.0, f64::max);
    for i in 0..binning.bins {
        let x = binning.center(i);
        let d = (fit.profile.pdf(x) - truth.pdf(x)).abs();
        assert!(
            d < 0.01 * peak,
            "x = {x}: {} vs {}",
            fit.profile.pdf(x),
            truth.pdf(x)
        );
    }
    // Its moments are those of the profile.
    let (m, s, _, _) = fit.profile.moments();
    let (mt, st, _, _) = truth.moments();
    assert!((m - mt).abs() < 0.01 * st && (s / st - 1.0).abs() < 0.02);
}

#[test]
fn dual_pearson_rejects_thin_histograms() {
    let mut h = Histogram::new(Binning::new(0.0, 1.0, 50).unwrap());
    for _ in 0..30 {
        h.fill(0.5);
    }
    assert!(matches!(
        DualPearson::fit(&h),
        Err(PearsonError::BadHistogram(_))
    ));
}

// ---- The ion tally on the BCA engine ----

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

/// `E_d` = 15 eV for every element (a test parameter, not data).
fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m
}

fn config() -> IonTallyConfig {
    IonTallyConfig {
        depth: Binning::new(0.0, 60.0 * NM, 60).unwrap(),
        lateral: Binning::new(-30.0 * NM, 30.0 * NM, 60).unwrap(),
        radial: Binning::new(0.0, 30.0 * NM, 30).unwrap(),
        escape_energy: Binning::new(0.0, 3000.0, 30).unwrap(),
        escape_polar: Binning::new(0.0, 0.5 * PI, 18).unwrap(),
    }
}

fn run_ion(
    stack: &Stack,
    z1: u8,
    e_ev: f64,
    n: u64,
    recoil_cutoff: f64,
    threads: usize,
) -> (IonReport, SummaryTally) {
    let ls = LindhardScharff::new();
    let bca = Bca::new(
        Beam::normal(Ion::new(z1).unwrap(), e_ev, n),
        stack,
        BcaConfig::new(5.0, recoil_cutoff),
        &ls,
        table(),
    )
    .unwrap();
    let proto = IonTally::new(stack, &bca.species_z(), config()).unwrap();
    with_threads(threads, || {
        let t = bca.run(|| proto.clone()).unwrap();
        let s = bca.run(|| SummaryTally::new(NM, 60)).unwrap();
        (t.report(false), s)
    })
}

#[test]
fn ion_tally_is_bit_identical_across_thread_counts() {
    let stack = Stack::new(vec![(si(), 8.0 * NM)], Some(si())).unwrap();
    let (r1, _) = run_ion(&stack, 5, 3e3, 300, 2.0, 1);
    for n in [2, 8] {
        let (r, _) = run_ion(&stack, 5, 3e3, 300, 2.0, n);
        assert_eq!(r1, r, "report differs on {n} threads");
        let a = serde_json::to_string(&r1).unwrap();
        assert_eq!(a, serde_json::to_string(&r).unwrap());
    }
}

#[test]
fn ion_tally_agrees_with_the_engine_bookkeeping() {
    // B into a 20 nm Si film: backscattering, transmission and sputtering
    // all occur, and the beam (Z = 5) can never replace a Si atom.
    let stack = Stack::new(vec![(si(), 20.0 * NM)], None).unwrap();
    let (r, s) = run_ion(&stack, 5, 5e3, 400, 2.0, 4);
    assert_eq!(r.histories, 400);
    assert_eq!(r.range.stopped, s.primaries_stopped);
    let beam = &r.escapes.species[0];
    assert!(beam.beam && beam.z == 5);
    assert_eq!(beam.front.count, s.backscattered);
    assert_eq!(beam.back.count, s.transmitted);
    assert_eq!(r.range.stopped + beam.front.count + beam.back.count, 400);
    let si_esc = &r.escapes.species[1];
    assert_eq!(si_esc.front.count, s.sputtered);
    assert_eq!(si_esc.back.count, s.recoils_transmitted);
    assert!(s.sputtered > 0 && s.transmitted > 0);
    assert!((r.escapes.sputter_yield - s.sputtered as f64 / 400.0).abs() < 1e-15);
    assert!((r.escapes.backscatter_coefficient - s.backscattered as f64 / 400.0).abs() < 1e-15);

    // Cascade bookkeeping: every displaced atom stopped (interstitial or
    // replacement) or left the target.
    let c = r.damage.cascade;
    assert_eq!(c.displacements, s.recoils);
    assert_eq!(c.vacancies, c.displacements - c.replacements);
    assert_eq!(
        c.displacements,
        c.interstitials + c.replacements + s.sputtered + s.recoils_transmitted
    );
    assert_eq!(r.damage.vacancy_histogram.total(), c.vacancies);
    assert_eq!(r.damage.interstitial_histogram.total(), c.interstitials);

    // NRT: one PKA per primary-target displacement, damage energy below T.
    let nrt = r.damage.nrt;
    assert!(nrt.pka_count > 0 && nrt.pka_count <= c.displacements);
    assert!(nrt.damage_energy_ev < nrt.pka_energy_ev);
    assert!(
        nrt.nrt_displacements > 0.0 && nrt.kinchin_pease_displacements >= nrt.nrt_displacements
    );

    // Range moments and the depth histogram.
    let d = r.range.depth.unwrap();
    assert!((d.mean - s.mean_depth()).abs() < 1e-12 * d.mean);
    assert_eq!(r.range.depth_histogram.total(), r.range.stopped);
    assert_eq!(r.range.layers.len(), 1);
    assert_eq!(r.range.layers[0].stopped, r.range.stopped);
    // Normal incidence: lateral distributions centred on the axis.
    let y = r.range.lateral_y.unwrap();
    assert!(y.mean.abs() < 4.0 * y.mean_std_err);
    // Polar angles of escapes lie in [0, pi/2].
    assert_eq!(si_esc.front.polar_histogram.overflow, 0);
}

#[test]
fn replacements_are_counted_for_self_ions() {
    // Si into Si with the recoil cutoff at E_d: moving Si atoms that stop
    // where they displaced a Si atom fill the site.
    let stack = Stack::semi_infinite(si());
    let (r, _) = run_ion(&stack, 14, 2e3, 100, 15.0, 2);
    let c = r.damage.cascade;
    assert!(c.replacements > 0, "{c:?}");
    assert!(c.replacements < c.displacements);
    assert_eq!(c.vacancies + c.replacements, c.displacements);
}

#[test]
fn report_carries_a_pearson_fit_or_its_reason() {
    let stack = Stack::semi_infinite(si());
    let ls = LindhardScharff::new();
    let bca = Bca::new(
        Beam::normal(Ion::new(5).unwrap(), 3e3, 1000),
        &stack,
        BcaConfig::new(5.0, 2.0),
        &ls,
        table(),
    )
    .unwrap();
    let proto = IonTally::new(&stack, &bca.species_z(), config()).unwrap();
    let r = bca.run(|| proto.clone()).unwrap().report(true);
    assert!(r.range.pearson_iv.is_some() != r.range.pearson_iv_error.is_some());
    assert!(r.range.dual_pearson.is_some() != r.range.dual_pearson_error.is_some());
    if let Some(p) = r.range.pearson_iv {
        let d = r.range.depth.unwrap();
        assert_eq!((p.mean, p.std_dev), (d.mean, d.std_dev));
    }
}

#[test]
fn dual_pearson_is_never_worse_than_the_single_pearson() {
    let p = PearsonIv::from_moments(80.0, 20.0, 0.6, 4.2).unwrap();
    let binning = Binning::new(0.0, 250.0, 80).unwrap();
    let hist = expected_histogram(binning, 2e5, |x| p.pdf(x));
    let fit = DualPearson::fit(&hist).unwrap();
    let single = fit.single_chi_square.expect("type IV moments");
    assert!(fit.chi_square <= single);
}
