//! Thermal vibration in the crystal flight model (issue #181, step 21b):
//! per-encounter Debye displacements of the lattice sites.
//!
//! The fast checks run in every `cargo test`: the static-lattice identity at
//! `T = 0` without the zero-point term, determinism across thread counts,
//! the run metadata, and that amorphous transport is untouched. The
//! statistical checks need thousands of ions and a release build, so they are
//! `#[ignore]`d; run them with
//! `cargo test --release -p lindhard --test crystal_thermal -- --ignored --nocapture`.
//! The measured numbers are in the module docs of `lindhard::ion::bca::crystal`.

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaError, Beam, CrystalTarget, ElectronicLoss, SummaryTally, Thermal,
};
use lindhard::ion::crystal::debye::{ThermalVibration, THETA_D_SI};
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
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

/// Silicon at the tabulated (crystal) density with `E_d` = 15 eV, as in
/// `tests/crystal_bca.rs`.
fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    if m.surface_binding_energy_ev(14).is_err() {
        m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    }
    m
}

fn beam(z: u8, e: f64, tilt_deg: f64, twist_deg: f64, count: u64) -> Beam {
    Beam {
        ion: Ion::new(z).unwrap(),
        energy_ev: e,
        polar_rad: tilt_deg.to_radians(),
        azimuth_rad: twist_deg.to_radians(),
        count,
    }
}

/// Si cut on (100), reference [010], matching `beam(.., tilt, twist, ..)`.
fn crystal_si(tilt_deg: f64, twist_deg: f64, thermal: Option<Thermal>) -> CrystalTarget {
    let lat = Lattice::silicon();
    let o = Orientation::new(
        &lat,
        [1, 0, 0],
        [0, 1, 0],
        tilt_deg.to_radians(),
        twist_deg.to_radians(),
        0.0,
    )
    .unwrap();
    let c = CrystalTarget::new(lat, o);
    match thermal {
        Some(t) => c.with_thermal(t),
        None => c,
    }
}

/// `T` K with the cited Si Debye temperature, zero-point term included.
fn at(t: f64) -> Option<Thermal> {
    Some(Thermal::new(t, THETA_D_SI))
}

fn run_summary(
    st: &Stack,
    b: Beam,
    seed: u64,
    crystal: Option<(CrystalTarget, Vec<usize>)>,
    tweak: impl Fn(&mut BcaConfig),
    threads: Option<usize>,
) -> SummaryTally {
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = seed;
    tweak(&mut cfg);
    let mut bca = Bca::new(b, st, cfg, &ls, table()).unwrap();
    if let Some((c, regions)) = crystal {
        bca = bca.with_crystal(c, &regions).unwrap();
    }
    let go = || bca.run(|| SummaryTally::new(0.5 * NM, 2000)).unwrap();
    match threads {
        None => go(),
        Some(n) => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(go),
    }
}

/// The complete serialised `IonReport` of a run.
fn report_json(
    st: &Stack,
    b: Beam,
    cfg: BcaConfig,
    crystal: Option<(CrystalTarget, Vec<usize>)>,
    threads: usize,
) -> String {
    let ls = LindhardScharff::new();
    let mut bca = Bca::new(b, st, cfg, &ls, table()).unwrap();
    if let Some((c, regions)) = crystal {
        bca = bca.with_crystal(c, &regions).unwrap();
    }
    let species = bca.species_z();
    let tc = IonTallyConfig {
        depth: Binning::new(0.0, 400.0 * NM, 400).unwrap(),
        lateral: Binning::new(-100.0 * NM, 100.0 * NM, 100).unwrap(),
        radial: Binning::new(0.0, 100.0 * NM, 50).unwrap(),
        escape_energy: Binning::new(0.0, 5.0e4, 50).unwrap(),
        escape_polar: Binning::new(0.0, std::f64::consts::FRAC_PI_2, 18).unwrap(),
    };
    let tally = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(|| {
            bca.run(|| IonTally::new(st, &species, tc).unwrap())
                .unwrap()
        });
    serde_json::to_string(&tally.report(false)).unwrap()
}

/// Fraction of stopped primaries deeper than `depth_m`.
fn fraction_deeper(t: &SummaryTally, depth_m: f64) -> f64 {
    let tot: u64 = t.depth_hist.iter().sum();
    let first = (depth_m / t.bin_width_m).ceil() as usize;
    t.depth_hist[first.min(t.depth_hist.len())..]
        .iter()
        .sum::<u64>() as f64
        / tot as f64
}

/// Depth below which a fraction `q` of the stopped primaries lie, nm.
fn percentile_nm(t: &SummaryTally, q: f64) -> f64 {
    let tot: u64 = t.depth_hist.iter().sum();
    let mut acc = 0u64;
    for (i, &h) in t.depth_hist.iter().enumerate() {
        acc += h;
        if acc as f64 >= q * tot as f64 {
            return (i as f64 + 0.5) * t.bin_width_m / NM;
        }
    }
    f64::NAN
}

/// At `T = 0` with the zero-point term switched off every amplitude is zero:
/// the thermal code path (margin search, per-encounter draws, re-sorting)
/// must then reproduce the static 20b lattice bit for bit, complete report
/// included, for a channeling and a random direction and with the optional
/// local electronic loss.
#[test]
fn zero_amplitude_reproduces_the_static_lattice_bit_for_bit() {
    let st = Stack::semi_infinite(si());
    let cold = Thermal {
        temperature_k: 0.0,
        debye_temperature_k: THETA_D_SI,
        include_zero_point: false,
    };
    for (label, z, e, tilt, twist, n, local) in [
        ("B 5 keV <110>", 5u8, 5.0e3, 45.0, 0.0, 120, false),
        ("As 20 keV 7/22, local loss", 33, 2.0e4, 7.0, 22.0, 60, true),
    ] {
        let mut cfg = BcaConfig::new(5.0, 2.0);
        cfg.seed = 5;
        if local {
            cfg.electronic = ElectronicLoss::EquipartitionLsOr;
        }
        let stat = report_json(
            &st,
            beam(z, e, tilt, twist, n),
            cfg,
            Some((crystal_si(tilt, twist, None), vec![0])),
            1,
        );
        let therm = report_json(
            &st,
            beam(z, e, tilt, twist, n),
            cfg,
            Some((crystal_si(tilt, twist, Some(cold)), vec![0])),
            1,
        );
        assert!(
            stat == therm,
            "{label}: T = 0 without zero point differs from static"
        );
    }
}

/// With vibration on, the result changes (the switch is not a no-op), and
/// the complete report is identical on 1, 2 and 8 threads.
#[test]
fn thermal_runs_are_identical_on_1_2_8_threads() {
    let st = Stack::semi_infinite(si());
    for (label, local) in [("nonlocal", false), ("local+nonlocal", true)] {
        let mut cfg = BcaConfig::new(5.0, 2.0);
        cfg.seed = 11;
        if local {
            cfg.electronic = ElectronicLoss::EquipartitionLsOr;
        }
        let run = |threads, thermal| {
            report_json(
                &st,
                beam(33, 1.0e4, 7.0, 22.0, 40),
                cfg,
                Some((crystal_si(7.0, 22.0, thermal), vec![0])),
                threads,
            )
        };
        let r1 = run(1, at(300.0));
        for n in [2, 8] {
            assert!(r1 == run(n, at(300.0)), "{label}: {n} threads differ");
        }
        assert!(r1 != run(1, None), "{label}: vibration changed nothing");
    }
    // Energy is conserved with vibration on (SummaryTally keeps the residual).
    let s = run_summary(
        &st,
        beam(5, 5.0e3, 45.0, 0.0, 100),
        3,
        Some((crystal_si(45.0, 0.0, at(600.0)), vec![0])),
        |_| {},
        Some(2),
    );
    assert!(s.max_relative_residual < 1e-9);
    assert!(s.recoils > 0);
}

/// An amorphous layer in front of a vibrating crystal that no particle
/// reaches: the report is that of the plain amorphous engine.
#[test]
fn unreachable_thermal_crystal_leaves_amorphous_run_identical() {
    let st = Stack::new(vec![(si(), 1000.0 * NM)], Some(si())).unwrap();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 1;
    let plain = report_json(&st, beam(5, 5.0e3, 7.0, 0.0, 400), cfg, None, 1);
    let with = report_json(
        &st,
        beam(5, 5.0e3, 7.0, 0.0, 400),
        cfg,
        Some((crystal_si(7.0, 0.0, at(300.0)), vec![1])),
        1,
    );
    assert!(plain == with);
}

/// The target temperature and the thermal model are recorded in the run
/// metadata, and serialise.
#[test]
fn temperature_is_recorded_in_the_run_metadata() {
    let st = Stack::semi_infinite(si());
    let ls = LindhardScharff::new();
    let cfg = BcaConfig::new(5.0, 2.0);
    let b = beam(5, 5.0e3, 7.0, 22.0, 1);
    let bca = Bca::new(b, &st, cfg, &ls, table())
        .unwrap()
        .with_crystal(crystal_si(7.0, 22.0, at(300.0)), &[0])
        .unwrap();
    let meta = bca.crystal_metadata();
    assert_eq!(meta.len(), 1);
    let m = &meta[0];
    assert_eq!(m.regions, vec![0]);
    let th = m.thermal.as_ref().expect("thermal metadata");
    assert_eq!(th.input.temperature_k, 300.0);
    assert_eq!(th.input.debye_temperature_k, THETA_D_SI);
    assert!(th.input.include_zero_point);
    assert_eq!(th.rms_1d_m.len(), 1);
    let (z, mass, u1) = th.rms_1d_m[0];
    assert_eq!(z, 14);
    let want = ThermalVibration::new(THETA_D_SI, mass, 300.0)
        .unwrap()
        .rms_1d();
    assert_eq!(u1, want);
    // 0.065 Angstrom for Si at 300 K with THETA_D_SI (the Debye formula).
    assert!((u1 / 0.0652e-10 - 1.0).abs() < 0.01, "u1 = {u1}");
    assert_eq!(th.search_margin_m, 6.0 * u1);
    let json = serde_json::to_string(&meta).unwrap();
    assert!(json.contains("\"temperature_k\":300.0"), "{json}");
    // A static crystal records no thermal model.
    let stat = Bca::new(b, &st, cfg, &ls, table())
        .unwrap()
        .with_crystal(crystal_si(7.0, 22.0, None), &[0])
        .unwrap();
    assert!(stat.crystal_metadata()[0].thermal.is_none());
}

#[test]
fn bad_thermal_input_is_reported() {
    let st = Stack::semi_infinite(si());
    let ls = LindhardScharff::new();
    let mk = || {
        Bca::new(
            beam(5, 5e3, 7.0, 22.0, 1),
            &st,
            BcaConfig::new(5.0, 2.0),
            &ls,
            table(),
        )
        .unwrap()
    };
    for th in [
        Thermal::new(-1.0, THETA_D_SI),
        Thermal::new(f64::NAN, THETA_D_SI),
        Thermal::new(300.0, 0.0),
    ] {
        match mk().with_crystal(crystal_si(7.0, 22.0, Some(th)), &[0]) {
            Err(BcaError::InvalidConfig(m)) => assert!(m.contains("thermal"), "{m}"),
            Err(e) => panic!("unexpected error {e}"),
            Ok(_) => panic!("{th:?} accepted"),
        }
    }
}

/// Binomial standard error of a fraction `f` of `n` ions.
fn binomial_sigma(f: f64, n: f64) -> f64 {
    (f * (1.0 - f) / n).sqrt()
}

/// The issue's channeling criterion: B 5 keV along <110>, fraction of ions
/// deeper than twice the amorphous Rp at T = 0 (zero-point motion only), 300
/// and 600 K. It must fall monotonically, each step by more than three
/// combined binomial standard errors.
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn channeling_tail_falls_with_temperature() {
    let n = 4000u64;
    let st = Stack::semi_infinite(si());
    let am = run_summary(&st, beam(5, 5.0e3, 45.0, 0.0, n), 1, None, |_| {}, None);
    let two_rp = 2.0 * am.mean_depth();
    let mut prev: Option<(f64, f64)> = None;
    let stat = run_summary(
        &st,
        beam(5, 5.0e3, 45.0, 0.0, n),
        1,
        Some((crystal_si(45.0, 0.0, None), vec![0])),
        |_| {},
        None,
    );
    println!(
        "amorphous Rp {:.1} nm; static lattice: Rp {:.1} nm, tail {:.4}",
        am.mean_depth() / NM,
        stat.mean_depth() / NM,
        fraction_deeper(&stat, two_rp)
    );
    for t in [0.0, 300.0, 600.0] {
        let cr = run_summary(
            &st,
            beam(5, 5.0e3, 45.0, 0.0, n),
            1,
            Some((crystal_si(45.0, 0.0, at(t)), vec![0])),
            |_| {},
            None,
        );
        let f = fraction_deeper(&cr, two_rp);
        let sf = binomial_sigma(f, cr.primaries_stopped as f64);
        println!(
            "T = {t} K: Rp {:.1} nm, tail fraction {f:.4} +- {sf:.4}",
            cr.mean_depth() / NM
        );
        if let Some((fp, sp)) = prev {
            let z = (fp - f) / (sp * sp + sf * sf).sqrt();
            println!("  drop {:.4} = {z:.1} sigma", fp - f);
            assert!(z > 3.0, "T = {t} K: tail {f} vs {fp} ({z:.1} sigma)");
        }
        prev = Some((f, sf));
    }
}

/// A random direction (30 degrees tilt, 17 degrees twist) at 300 K: Rp
/// within 10 % of amorphous and the 90th percentile within 10 %, for B 5 keV
/// and As 30 keV. The static lattice's 15 % dRp bound is not applied at
/// 300 K (#225): vibration feeds a few per cent of the ions into a tail
/// (`off_axis_tail_grows_with_vibration_under_both_losses`), which the
/// second moment weights heavily, so dRp is only held to its recorded value
/// (1.18 B, 1.17 As) within 0.12.
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn random_direction_rp_at_300_k_matches_amorphous() {
    let st = Stack::semi_infinite(si());
    for (z, e) in [(5u8, 5.0e3), (33, 3.0e4)] {
        let am = run_summary(&st, beam(z, e, 30.0, 17.0, 4000), 1, None, |_| {}, None);
        let cr = run_summary(
            &st,
            beam(z, e, 30.0, 17.0, 4000),
            1,
            Some((crystal_si(30.0, 17.0, at(300.0)), vec![0])),
            |_| {},
            None,
        );
        let rp = cr.mean_depth() / am.mean_depth();
        let drp = cr.depth_std() / am.depth_std();
        let p90 = percentile_nm(&cr, 0.9) / percentile_nm(&am, 0.9);
        println!("Z={z} 30/17 at 300 K: Rp ratio {rp:.3}, dRp ratio {drp:.3}, p90 ratio {p90:.3}");
        assert!((rp - 1.0).abs() < 0.10, "Z={z}: Rp ratio {rp}");
        assert!((p90 - 1.0).abs() < 0.10, "Z={z}: p90 ratio {p90}");
        assert!(drp < 1.18 + 0.12, "Z={z}: dRp ratio {drp}");
    }
}

/// The 7 degrees / 22 degrees orientation at 300 K (rescoped in #225).
///
/// That beam lies 2.6 to 2.7 degrees from a {100} and a {110} plane, so it is
/// not a random direction, and the #180 bound "dRp within 15 % of
/// amorphous" is not the criterion here (module docs of
/// `lindhard::ion::bca::crystal`, "The 7°/22° criterion"). No measured
/// profile at matched conditions has been found, so this is a
/// recorded-value regression check: the dRp ratio stays within
/// [`DRP_TOL_B`] / [`DRP_TOL_AS`] of the recorded value, the 90th percentile
/// stays within 30 % of the amorphous one (the excess width is a tail), and
/// the fraction of ions beyond twice the amorphous Rp is above the amorphous
/// fraction by more than five binomial standard errors (the tail exists).
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn near_planar_7_22_at_300_k() {
    let st = Stack::semi_infinite(si());
    for (z, e, drp_rec, tol) in [(5u8, 5.0e3, 1.43, DRP_TOL_B), (33, 3.0e4, 1.80, DRP_TOL_AS)] {
        let am = run_summary(&st, beam(z, e, 7.0, 22.0, 4000), 1, None, |_| {}, None);
        let cr = run_summary(
            &st,
            beam(z, e, 7.0, 22.0, 4000),
            1,
            Some((crystal_si(7.0, 22.0, at(300.0)), vec![0])),
            |_| {},
            None,
        );
        let rp = cr.mean_depth() / am.mean_depth();
        let drp = cr.depth_std() / am.depth_std();
        let two_rp = 2.0 * am.mean_depth();
        let (f, fa) = (fraction_deeper(&cr, two_rp), fraction_deeper(&am, two_rp));
        let sf = binomial_sigma(f, cr.primaries_stopped as f64);
        println!(
            "Z={z} 7/22 at 300 K: Rp ratio {rp:.3}, dRp ratio {drp:.3}, median {:.1} vs {:.1} nm, \
             p90 {:.1} vs {:.1} nm, p99 {:.1} vs {:.1} nm, tail {f:.4} +- {sf:.4} vs {fa:.4}",
            percentile_nm(&cr, 0.5),
            percentile_nm(&am, 0.5),
            percentile_nm(&cr, 0.9),
            percentile_nm(&am, 0.9),
            percentile_nm(&cr, 0.99),
            percentile_nm(&am, 0.99)
        );
        assert!(
            (drp - drp_rec).abs() <= tol,
            "Z={z}: dRp ratio {drp} moved from the recorded {drp_rec}; update the docs"
        );
        assert!(percentile_nm(&cr, 0.9) < 1.3 * percentile_nm(&am, 0.9));
        assert!(f - fa > 5.0 * sf, "Z={z}: tail {f} vs amorphous {fa}");
    }
}

/// Tolerance of the recorded dRp ratios at 7 degrees / 22 degrees (B 5 keV):
/// three times the largest seed-to-seed standard deviation over seeds 1-5
/// (0.018 static, 0.037 at 300 K), rounded up. It also covers the difference
/// between platforms of the crystal path (`docs/architecture.md`).
const DRP_TOL_B: f64 = 0.12;

/// The same for As 30 keV (seed-to-seed standard deviation 0.052 static,
/// 0.047 at 300 K).
const DRP_TOL_AS: f64 = 0.16;

/// The finding of #225: at the off-axis directions 7 degrees / 22 degrees
/// and 30 degrees / 17 degrees, the fraction of B 5 keV ions deeper than
/// twice the amorphous Rp is **larger** at 300 K than in the static lattice,
/// with the nonlocal loss and with the Oen-Robinson local half alike (the
/// amorphous reference uses the same loss). Each rise must exceed three
/// combined binomial standard errors (16000 ions; measured 5 to 9). This
/// records a property of the model, not a validated physical result: if it
/// flips, the module docs of `lindhard::ion::bca::crystal` must change.
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn off_axis_tail_grows_with_vibration_under_both_losses() {
    let n = 16000u64;
    let st = Stack::semi_infinite(si());
    for (tilt, twist) in [(7.0, 22.0), (30.0, 17.0)] {
        for loss in [ElectronicLoss::NonLocal, ElectronicLoss::EquipartitionLsOr] {
            let tw = move |c: &mut BcaConfig| c.electronic = loss;
            let am = run_summary(&st, beam(5, 5.0e3, tilt, twist, n), 1, None, tw, None);
            let two_rp = 2.0 * am.mean_depth();
            let tail = |th: Option<Thermal>| {
                let cr = run_summary(
                    &st,
                    beam(5, 5.0e3, tilt, twist, n),
                    1,
                    Some((crystal_si(tilt, twist, th), vec![0])),
                    tw,
                    None,
                );
                let f = fraction_deeper(&cr, two_rp);
                println!(
                    "{tilt}/{twist} {loss:?} T = {:?}: Rp ratio {:.3}, dRp ratio {:.3}, \
                     tail {f:.4} (amorphous {:.4})",
                    th.map(|t| t.temperature_k),
                    cr.mean_depth() / am.mean_depth(),
                    cr.depth_std() / am.depth_std(),
                    fraction_deeper(&am, two_rp)
                );
                (f, binomial_sigma(f, cr.primaries_stopped as f64))
            };
            let (f0, s0) = tail(None);
            let (f3, s3) = tail(at(300.0));
            let z = (f3 - f0) / (s0 * s0 + s3 * s3).sqrt();
            println!("  rise {:.4} = {z:.1} sigma", f3 - f0);
            assert!(z > 3.0, "{tilt}/{twist} {loss:?}: tail {f0} -> {f3}");
        }
    }
}
