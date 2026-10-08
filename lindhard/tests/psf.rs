//! `tally::psf`: recovery of known double- and triple-Gaussian PSFs from
//! profiles with Poisson noise, normalisation of the fit against the energy
//! deposited in the slab of a transport run, bit-identical profiles and fits
//! across thread counts, and the CSV and JSON export.
//!
//! The synthetic profiles are expected bin energies of the PSF forms of
//! `tally::psf` (Mao et al. (2025) eq. (1); Rosa Figueiro (2015) eq. (72)) at
//! parameters chosen here, with Poisson counting noise drawn from
//! `lindhard::rng` streams. The transport runs use cross-section tables
//! computed here from simple formulas (constant mean free paths, isotropic
//! angles), as in `tests/electron_tally.rs`. No measured or tabulated data.

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::{LayerTables, Primary, Transport, TransportConfig};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::rng::stream;
use lindhard::tally::{
    fit_psf, Binning, CylindricalGrid, ElectronReport, ElectronTallyConfig, FullElectronTally,
    GaussianPsf, LogRadialBinning, PsfConfig, PsfError, PsfFit, PsfFitOptions, PsfModel,
    PsfNormalization, PsfReport, RadialProfile,
};
use rand_core::Rng;

const NM: f64 = 1e-9;

// ---------------------------------------------------------------------------
// Synthetic profiles with Poisson noise
// ---------------------------------------------------------------------------

fn uniform(rng: &mut impl Rng) -> f64 {
    (rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// Poisson(`mean`) by counting the arrivals of a unit-rate Poisson process
/// in `[0, mean]` (exact; cost proportional to `mean`).
fn poisson(rng: &mut impl Rng, mean: f64) -> u64 {
    let mut n = 0;
    let mut t = -(1.0 - uniform(rng)).ln();
    while t <= mean {
        n += 1;
        t -= (1.0 - uniform(rng)).ln();
    }
    n
}

/// The profile of `truth` on `bins`, as counts of quanta of 1 eV with
/// Poisson noise, each bin from its own stream; the error of a bin with `n`
/// counts is `n^(1/2)`.
fn noisy_profile(truth: &GaussianPsf, bins: LogRadialBinning, seed: u64) -> RadialProfile {
    let edges = bins.edges();
    let n = bins.len();
    let mut energy = Vec::with_capacity(n);
    let mut err = Vec::with_capacity(n);
    for i in 0..n {
        let mean = truth.annulus_energy_ev(edges[i], edges[i + 1]);
        let k = poisson(&mut stream(seed, i as u64), mean) as f64;
        energy.push(k);
        err.push(k.sqrt());
    }
    let beyond = poisson(
        &mut stream(seed, n as u64),
        truth.annulus_energy_ev(edges[n], 1.0),
    ) as f64;
    let total: f64 = energy.iter().sum::<f64>() + beyond;
    let p = RadialProfile {
        histories: total as u64,
        depth_lo_m: 0.0,
        depth_hi_m: 1.0 * NM,
        edges_m: edges[..=n].to_vec(),
        energy_ev: energy,
        std_err_ev: err,
        beyond_ev: beyond,
        beyond_std_err_ev: beyond.sqrt(),
        total_ev: total,
        total_std_err_ev: total.sqrt(),
    };
    p.validate().unwrap();
    p
}

fn bins() -> LogRadialBinning {
    // 0.2 nm to 30 µm: below every α and beyond 3β for every set.
    LogRadialBinning::new(0.2 * NM, 30_000.0 * NM, 100).unwrap()
}

/// Every fitted parameter within `k` standard errors of the truth.
fn within(fit: &PsfFit, truth: &GaussianPsf, k: f64) -> Result<(), String> {
    for ((name, v), (e, t)) in fit
        .parameter_names
        .iter()
        .zip(&fit.values)
        .zip(fit.std_errors.iter().zip(truth.values()))
    {
        let pull = (v - t) / e;
        if pull.abs() > k {
            return Err(format!(
                "{name}: fitted {v:e} +- {e:e}, truth {t:e}, pull {pull:.2}"
            ));
        }
    }
    Ok(())
}

/// Three double-Gaussian sets spanning β/α = 10, 100 and 1000.
fn double_sets() -> [GaussianPsf; 3] {
    [
        GaussianPsf::double(1.0e6, 20.0 * NM, 200.0 * NM, 0.6).unwrap(),
        GaussianPsf::double(1.0e6, 5.0 * NM, 500.0 * NM, 0.8).unwrap(),
        GaussianPsf::double(1.0e6, 3.0 * NM, 3_000.0 * NM, 0.4).unwrap(),
    ]
}

#[test]
fn double_gaussian_is_recovered_within_two_sigma_for_beta_over_alpha_10_to_1000() {
    for (k, truth) in double_sets().iter().enumerate() {
        let ratio = truth.beta_m / truth.alpha_m;
        let profile = noisy_profile(truth, bins(), 100 + k as u64);
        for normalization in [PsfNormalization::SlabTotal, PsfNormalization::Free] {
            let opts = PsfFitOptions {
                normalization,
                ..PsfFitOptions::default()
            };
            let fit = fit_psf(&profile, PsfModel::DoubleGaussian, &opts).unwrap();
            assert!(
                fit.converged,
                "beta/alpha = {ratio}, {normalization:?}: not converged"
            );
            within(&fit, truth, 2.0)
                .unwrap_or_else(|e| panic!("beta/alpha = {ratio}, {normalization:?}: {e}"));
        }
        let fit = fit_psf(
            &profile,
            PsfModel::DoubleGaussian,
            &PsfFitOptions::default(),
        )
        .unwrap();
        // A correct form with correct errors: reduced chi² near 1.
        let tol = 4.0 * (2.0 / fit.dof as f64).sqrt();
        assert!(
            (fit.reduced_chi2 - 1.0).abs() < tol,
            "beta/alpha = {ratio}: reduced chi2 {}",
            fit.reduced_chi2
        );
        assert_eq!(fit.residuals.len(), profile.len());
    }
}

#[test]
fn double_gaussian_errors_cover_the_truth_over_replicas() {
    // 40 noise replicas of the beta/alpha = 100 set: the fraction of
    // parameters within 2 sigma should be near 95 %. The bound allows for
    // the binomial spread of 160 draws (sd 1.7 %) and the correlation of the
    // four parameters of one replica.
    let truth = double_sets()[1];
    let (mut inside, mut all) = (0, 0);
    for seed in 0..40 {
        let p = noisy_profile(&truth, bins(), 1_000 + seed);
        let fit = fit_psf(&p, PsfModel::DoubleGaussian, &PsfFitOptions::default()).unwrap();
        for ((v, e), t) in fit.values.iter().zip(&fit.std_errors).zip(truth.values()) {
            all += 1;
            if ((v - t) / e).abs() <= 2.0 {
                inside += 1;
            }
        }
    }
    let frac = inside as f64 / all as f64;
    assert!(
        frac > 0.85,
        "only {inside} of {all} parameters within 2 sigma"
    );
}

#[test]
fn triple_gaussian_is_recovered_within_two_sigma() {
    let truth = GaussianPsf::triple(4.0e6, 4.0 * NM, 80.0 * NM, 0.5, 2_000.0 * NM, 0.3).unwrap();
    let profile = noisy_profile(&truth, bins(), 7);
    let fit = fit_psf(
        &profile,
        PsfModel::TripleGaussian,
        &PsfFitOptions::default(),
    )
    .unwrap();
    assert!(fit.converged);
    within(&fit, &truth, 2.0).unwrap();
    assert_eq!(fit.psf.model(), PsfModel::TripleGaussian);
    // The double form cannot describe three ranges.
    let d = fit_psf(
        &profile,
        PsfModel::DoubleGaussian,
        &PsfFitOptions::default(),
    )
    .unwrap();
    assert!(d.reduced_chi2 > 10.0 * fit.reduced_chi2);
}

#[test]
fn a_supplied_start_reaches_the_grid_start_minimum() {
    let truth = double_sets()[0];
    let profile = noisy_profile(&truth, bins(), 3);
    let grid = fit_psf(
        &profile,
        PsfModel::DoubleGaussian,
        &PsfFitOptions::default(),
    )
    .unwrap();
    let start = GaussianPsf::double(5.0e5, 10.0 * NM, 500.0 * NM, 1.5).unwrap();
    let opts = PsfFitOptions {
        start: Some(start),
        ..PsfFitOptions::default()
    };
    let s = fit_psf(&profile, PsfModel::DoubleGaussian, &opts).unwrap();
    for (a, b) in grid.values.iter().zip(&s.values) {
        assert!((a - b).abs() <= 1e-6 * a.abs(), "{a} vs {b}");
    }
    // A start of the wrong model is rejected.
    let wrong = PsfFitOptions {
        start: Some(GaussianPsf::triple(1.0, 1e-9, 1e-8, 1.0, 1e-7, 1.0).unwrap()),
        ..PsfFitOptions::default()
    };
    assert!(fit_psf(&profile, PsfModel::DoubleGaussian, &wrong).is_err());
}

// ---------------------------------------------------------------------------
// Transport runs
// ---------------------------------------------------------------------------

const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];

fn material() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

fn table(axis: SamplingAxis, rate: f64, row: impl Fn(f64, f64) -> f64) -> CrossSectionTable {
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "computed in tests/psf.rs".into(),
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

/// Isotropic elastic scattering and inelastic losses uniform on `[0, 0.3 E]`.
fn transport() -> Transport {
    Transport::new(
        Stack::semi_infinite(material()),
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

/// The PSF slab and a cylindrical grid whose single depth bin is the same
/// slab, for a cross-check of the slab total.
fn config() -> ElectronTallyConfig {
    let mut c = ElectronTallyConfig::new(
        Binning::new(0.0, 2_000.0, 20).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    c.psf = Some(
        PsfConfig::new(
            LogRadialBinning::new(0.05 * NM, 500.0 * NM, 60).unwrap(),
            0.0,
            4.0 * NM,
        )
        .unwrap(),
    );
    c.cylindrical = Some(CylindricalGrid {
        r: Binning::new(0.0, 1.0, 1).unwrap(),
        depth: Binning::new(0.0, 4.0 * NM, 1).unwrap(),
    });
    c
}

fn run(threads: usize, n: u64) -> ElectronReport {
    let t = transport();
    let proto = FullElectronTally::new(&t, config()).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(|| {
            t.run(0x95F, n, 64, &Primary::normal(2_000.0), || proto.clone())
                .unwrap()
                .tally
                .report()
        })
}

/// `∫₀^∞ E f(r) 2πr dr` by the trapezoidal rule on a fine log grid in `r`
/// (independent of the closed form the fit uses).
fn integral(psf: &GaussianPsf) -> f64 {
    let (lo, hi, n) = (1e-14_f64, 1e-3_f64, 200_000);
    let q = (hi / lo).ln() / n as f64;
    // dr = r d(ln r): integrand in ln r is E f(r) 2π r².
    let g = |r: f64| psf.density_ev_per_m2(r) * 2.0 * std::f64::consts::PI * r * r;
    let mut s = 0.5 * (g(lo) + g(hi));
    for i in 1..n {
        s += g(lo * (q * i as f64).exp());
    }
    // The disc below `lo` holds a negligible fraction.
    s * q
}

// Under `PsfNormalization::SlabTotal` the energy is fixed to the slab energy, so this
// checks the normalisation of the coded forms, not the quality of the fit. The
// fit-quality reading of the criterion is
// `a_free_scale_integrates_to_the_slab_energy_when_the_form_fits`.
#[test]
fn fitted_psf_integrates_to_the_slab_energy_within_one_percent() {
    let r = run(4, 4_000);
    let p = r.deposition.psf.as_ref().unwrap();
    // The slab total is the energy deposited in the slab: the one-cell
    // cylindrical grid over the same depths holds the same energy, and the
    // bins plus the part beyond r_max add up to it.
    let cyl = r.deposition.cylindrical.as_ref().unwrap();
    let slab = cyl.energy_ev[0];
    assert!(slab > 0.0);
    assert!(
        (p.total_ev - slab).abs() <= 1e-9 * slab,
        "{} vs {slab}",
        p.total_ev
    );
    let bins: f64 = p.energy_ev.iter().sum();
    assert!((bins + p.beyond_ev - p.total_ev).abs() <= 1e-9 * p.total_ev);
    assert!(p.total_std_err_ev > 0.0 && p.total_std_err_ev < 0.1 * p.total_ev);
    assert!(p.std_err_ev.iter().filter(|e| **e > 0.0).count() > 40);
    for model in [PsfModel::DoubleGaussian, PsfModel::TripleGaussian] {
        let fit = fit_psf(p, model, &PsfFitOptions::default()).unwrap();
        assert_eq!(fit.normalization, PsfNormalization::SlabTotal);
        // An independent quadrature of the fitted density over the plane.
        let int = integral(&fit.psf);
        let rel = (int - p.total_ev).abs() / p.total_ev;
        assert!(
            rel < 0.01,
            "{model:?}: integral {int:e} eV vs slab {:e} eV ({:.3} %)",
            p.total_ev,
            100.0 * rel
        );
        assert_eq!(fit.std_errors[0], p.total_std_err_ev);
        // The bins the fit covers hold the slab energy less what lies
        // beyond r_max, here nothing.
        let model_sum: f64 = fit.residuals.iter().map(|r| r.model_ev).sum();
        assert!((model_sum - p.total_ev).abs() < 0.01 * p.total_ev);
    }
    // Why the scale is fixed by default: the double Gaussian does not
    // describe this profile (a quarter of the slab energy lies on the axis,
    // r < 0.05 nm), and with E free its integral misses the slab energy by
    // far more than 1 % (about a third in this run).
    let free = PsfFitOptions {
        normalization: PsfNormalization::Free,
        ..PsfFitOptions::default()
    };
    let fit = fit_psf(p, PsfModel::DoubleGaussian, &free).unwrap();
    assert!(fit.reduced_chi2 > 10.0);
    assert!((integral(&fit.psf) - p.total_ev).abs() > 0.1 * p.total_ev);
}

#[test]
fn a_free_scale_integrates_to_the_slab_energy_when_the_form_fits() {
    // A profile the triple Gaussian describes, with E fitted.
    let truth = GaussianPsf::triple(2.0e6, 4.0 * NM, 80.0 * NM, 0.5, 2_000.0 * NM, 0.3).unwrap();
    let p = noisy_profile(&truth, bins(), 21);
    let opts = PsfFitOptions {
        normalization: PsfNormalization::Free,
        ..PsfFitOptions::default()
    };
    let fit = fit_psf(&p, PsfModel::TripleGaussian, &opts).unwrap();
    let int = integral(&fit.psf);
    assert!((int - fit.psf.total_ev).abs() <= 1e-6 * int);
    assert!(
        (int - p.total_ev).abs() < 0.01 * p.total_ev,
        "{int} vs {}",
        p.total_ev
    );
    assert!(fit.reduced_chi2 < 2.0);
}

#[test]
fn profile_and_fit_are_bit_identical_for_1_2_and_8_threads() {
    let fit_of = |r: &ElectronReport| {
        let p = r.deposition.psf.clone().unwrap();
        let fits = [PsfModel::DoubleGaussian, PsfModel::TripleGaussian]
            .into_iter()
            .map(|m| fit_psf(&p, m, &PsfFitOptions::default()).unwrap())
            .collect();
        PsfReport { profile: p, fits }
    };
    let r1 = run(1, 1_500);
    let f1 = fit_of(&r1);
    let j1 = serde_json::to_string_pretty(&f1).unwrap();
    for n in [2, 8] {
        let r = run(n, 1_500);
        assert_eq!(r1, r, "report differs on {n} threads");
        let f = fit_of(&r);
        assert_eq!(f1, f, "PSF fit differs on {n} threads");
        assert_eq!(j1, serde_json::to_string_pretty(&f).unwrap());
        assert_eq!(f1.profile_csv(), f.profile_csv());
        assert_eq!(f1.parameters_csv(), f.parameters_csv());
    }
}

// ---------------------------------------------------------------------------
// Export
// ---------------------------------------------------------------------------

#[test]
fn csv_and_json_export() {
    let truth = double_sets()[0];
    let profile = noisy_profile(&truth, bins(), 11);
    let fits: Vec<PsfFit> = [PsfModel::DoubleGaussian, PsfModel::TripleGaussian]
        .into_iter()
        .map(|m| fit_psf(&profile, m, &PsfFitOptions::default()).unwrap())
        .collect();
    let report = PsfReport {
        profile: profile.clone(),
        fits,
    };

    let csv = report.profile_csv();
    let lines: Vec<&str> = csv.lines().collect();
    assert_eq!(lines.len(), profile.len() + 1);
    let header: Vec<&str> = lines[0].split(',').collect();
    assert_eq!(header[0], "r_lo_m");
    assert_eq!(header[8], "model_double_gaussian_ev");
    assert_eq!(header[9], "model_triple_gaussian_ev");
    for (i, l) in lines[1..].iter().enumerate() {
        let f: Vec<f64> = l.split(',').map(|x| x.parse().unwrap()).collect();
        assert_eq!(f.len(), header.len());
        assert_eq!(f[0], profile.edges_m[i]);
        assert_eq!(f[4], profile.energy_ev[i]);
        assert_eq!(f[8], report.fits[0].residuals[i].model_ev);
    }

    let pcsv = report.parameters_csv();
    assert!(pcsv.starts_with("model,parameter,value,std_error\n"));
    assert_eq!(pcsv.lines().count(), 1 + (4 + 3) + (6 + 3));
    let alpha = pcsv
        .lines()
        .find(|l| l.starts_with("double_gaussian,alpha_m,"))
        .unwrap();
    let v: f64 = alpha.split(',').nth(2).unwrap().parse().unwrap();
    assert_eq!(v, report.fits[0].psf.alpha_m);

    let json = serde_json::to_string_pretty(&report).unwrap();
    let back: PsfReport = serde_json::from_str(&json).unwrap();
    assert_eq!(back.fits.len(), 2);
    assert_eq!(back.fits[0].model, PsfModel::DoubleGaussian);
    assert_eq!(back.fits[1].parameter_names.len(), 6);
    assert_eq!(back.profile.energy_ev, profile.energy_ev);
    let v: serde_json::Value = serde_json::from_str(&json).unwrap();
    assert!(v["fits"][0]["covariance"].as_array().unwrap().len() == 4);
    assert!(v["fits"][0]["source"].as_str().unwrap().contains("eq. (1)"));
    assert!(v["fits"][1]["source"]
        .as_str()
        .unwrap()
        .contains("eq. (72)"));
}

#[test]
fn too_few_bins_and_empty_profiles_are_rejected() {
    let cfg = PsfConfig::new(LogRadialBinning::new(1e-9, 1e-6, 3).unwrap(), 0.0, 1e-9).unwrap();
    let acc = lindhard::tally::RadialAccumulator::new(cfg).unwrap();
    let p = acc.profile();
    assert!(fit_psf(&p, PsfModel::DoubleGaussian, &PsfFitOptions::default()).is_err());
}

/// Each mutation of an otherwise valid profile's slab scalars or depth bounds
/// is rejected as `PsfError::Profile` by `validate` and by `fit_psf`, for
/// both normalisations, before the fit runs.
#[test]
fn invalid_slab_scalars_and_depth_bounds_are_rejected_as_profile_errors() {
    let truth = &double_sets()[0];
    let valid = noisy_profile(truth, bins(), 7);
    type Mutation = fn(&mut RadialProfile);
    let cases: [(&str, Mutation); 14] = [
        ("total_ev = NaN", |p| p.total_ev = f64::NAN),
        ("total_ev = inf", |p| p.total_ev = f64::INFINITY),
        ("total_ev < 0", |p| p.total_ev = -p.total_ev),
        ("total_ev inconsistent", |p| p.total_ev *= 1.01),
        ("total_std_err_ev = NaN", |p| p.total_std_err_ev = f64::NAN),
        ("total_std_err_ev < 0", |p| p.total_std_err_ev = -1.0),
        ("beyond_ev = NaN", |p| p.beyond_ev = f64::NAN),
        ("beyond_ev < 0", |p| {
            // Kept consistent with total_ev, so only the sign is wrong.
            p.beyond_ev -= 1.0;
            p.total_ev -= 1.0;
            assert!(p.beyond_ev < 0.0 && p.total_ev > 0.0);
        }),
        ("beyond_std_err_ev = NaN", |p| {
            p.beyond_std_err_ev = f64::NAN
        }),
        ("beyond_std_err_ev < 0", |p| p.beyond_std_err_ev = -1.0),
        ("depth bounds reversed", |p| {
            std::mem::swap(&mut p.depth_lo_m, &mut p.depth_hi_m)
        }),
        ("depth bounds equal", |p| p.depth_hi_m = p.depth_lo_m),
        ("depth_lo_m = NaN", |p| p.depth_lo_m = f64::NAN),
        ("depth_hi_m = inf", |p| p.depth_hi_m = f64::INFINITY),
    ];
    assert_eq!(valid.beyond_ev, 0.0, "the fixture changed");
    for (name, mutate) in cases {
        let mut p = valid.clone();
        mutate(&mut p);
        assert!(
            matches!(p.validate(), Err(PsfError::Profile(_))),
            "{name}: validate gave {:?}",
            p.validate()
        );
        for normalization in [PsfNormalization::SlabTotal, PsfNormalization::Free] {
            let opts = PsfFitOptions {
                normalization,
                ..PsfFitOptions::default()
            };
            let r = fit_psf(&p, PsfModel::DoubleGaussian, &opts);
            assert!(
                matches!(r, Err(PsfError::Profile(_))),
                "{name}, {normalization:?}: fit_psf gave {r:?}"
            );
        }
    }
    // The unmutated profile still fits, with finite errors throughout.
    let fit = fit_psf(&valid, PsfModel::DoubleGaussian, &PsfFitOptions::default()).unwrap();
    assert!(fit.std_errors.iter().all(|e| e.is_finite()));
    assert!(fit.covariance.iter().flatten().all(|c| c.is_finite()));
}
