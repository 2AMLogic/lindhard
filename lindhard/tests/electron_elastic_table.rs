//! Energy-grid elastic tables (`electron::elastic::table`, issue #90).
//!
//! The potentials are the Thomas-Fermi-length Yukawa **stand-in**
//! (`ThomasFermiYukawa`): the Salvat et al. (1987) DHFS table is not in the
//! tree (`docs/data-provenance.md`). To keep debug-mode test time modest the
//! energy grid here is 10 eV to 50 keV at 3 points per decade (13 points, both
//! ends of the default range included); the default grid is exercised by the
//! ignored `full_default_grid_recovers_sigma_tr1` test and the benchmark.

use std::sync::OnceLock;

use lindhard::electron::data::{ElectronDataError, SamplingAxis, CACHE_FORMAT_VERSION};
use lindhard::electron::elastic::table::{
    build_elastic_table, combine, default_energy_grid, default_probability_grid, log_energy_grid,
    mean_one_minus_cos, AtomicElastic, ElasticTableError, ElasticTableOptions, PotentialSource,
    SalvatDhfsTable, ThomasFermiYukawa, DEFAULT_REFINE_TOLERANCE,
};
use lindhard::electron::elastic::{ElasticError, ScreenedPotential, SolverOptions};
use lindhard::material::Material;

fn test_grid() -> Vec<f64> {
    log_energy_grid(10.0, 5.0e4, 3.0).unwrap()
}

fn atom(z: u8) -> AtomicElastic {
    let src = ThomasFermiYukawa;
    let pot = src.potential(z).unwrap();
    AtomicElastic::compute(
        z,
        pot.as_ref(),
        &src.description(),
        &test_grid(),
        SolverOptions::default(),
    )
    .unwrap()
}

fn si() -> &'static AtomicElastic {
    static A: OnceLock<AtomicElastic> = OnceLock::new();
    A.get_or_init(|| atom(14))
}

fn au() -> &'static AtomicElastic {
    static A: OnceLock<AtomicElastic> = OnceLock::new();
    A.get_or_init(|| atom(79))
}

fn pure(z: u8) -> Material {
    Material::from_atom_fractions(&[(z, 1.0)], None).unwrap()
}

fn table_of(
    material: &Material,
    atoms: &[AtomicElastic],
) -> lindhard::electron::data::CrossSectionTable {
    combine(
        material,
        atoms,
        &default_probability_grid(),
        Some(DEFAULT_REFINE_TOLERANCE),
    )
    .unwrap()
}

/// Worst relative error of `sigma_tr1` recovered from the stored CDF.
fn worst_sigma_tr1_error(
    atom: &AtomicElastic,
    table: &lindhard::electron::data::CrossSectionTable,
) -> (f64, f64) {
    let mut worst = (0.0_f64, 0.0);
    for (i, row) in atom.rows().iter().enumerate() {
        let q = table.quantiles(i).expect("non-zero rate");
        let recovered = row.sigma_el() * mean_one_minus_cos(table.probability(), q);
        let err = (recovered / row.sigma_tr1() - 1.0).abs();
        if err > worst.0 {
            worst = (err, row.energy_ev());
        }
    }
    worst
}

#[test]
fn default_energy_grid_spans_10ev_to_50kev_at_20_per_decade() {
    let g = default_energy_grid();
    assert_eq!(g[0], 10.0);
    assert_eq!(*g.last().unwrap(), 5.0e4);
    let decades = (5.0e4f64 / 10.0).log10();
    assert!((g.len() - 1) as f64 >= 20.0 * decades, "{} points", g.len());
    for w in g.windows(2) {
        assert!((w[1] / w[0]).log10() <= 1.0 / 20.0 + 1e-12);
    }
    assert_eq!(ElasticTableOptions::default().energy_ev, g);
}

#[test]
fn sigma_tr1_recovered_from_cdf_matches_solver() {
    for (atom, z) in [(si(), 14), (au(), 79)] {
        let t = table_of(&pure(z), std::slice::from_ref(atom));
        assert_eq!(t.axis(), SamplingAxis::ElasticPolarAngle);
        assert_eq!(t.energy_ev(), test_grid().as_slice());
        let (err, at) = worst_sigma_tr1_error(atom, &t);
        assert!(err < 1e-3, "Z={z}: sigma_tr1 off by {err:.3e} at {at} eV");
        for i in 0..t.energy_ev().len() {
            let q = t.quantiles(i).unwrap();
            assert_eq!(q[0], 0.0);
            assert_eq!(*q.last().unwrap(), std::f64::consts::PI);
        }
    }
}

#[test]
fn inverse_mean_free_path_is_n_sigma() {
    let m = pure(14);
    let t = table_of(&m, &[si().clone()]);
    let n = m.number_density_of(14).unwrap();
    for (row, rate) in si().rows().iter().zip(t.inverse_mfp_per_m()) {
        let want = n * row.sigma_el() * lindhard::electron::elastic::BOHR2_TO_M2;
        assert_eq!(*rate, want);
    }
}

#[test]
fn compound_is_the_number_density_weighted_sum_of_its_elements() {
    let mix = Material::from_atom_fractions(&[(14, 1.0), (79, 1.0)], Some(1.0e4)).unwrap();
    let (m_si, m_au) = (pure(14), pure(79));
    let atoms = [si().clone(), au().clone()];
    let t_mix = table_of(&mix, &atoms);
    let t_si = table_of(&m_si, &atoms[..1]);
    let t_au = table_of(&m_au, &atoms[1..]);
    let w_si = mix.number_density_of(14).unwrap() / m_si.number_density_of(14).unwrap();
    let w_au = mix.number_density_of(79).unwrap() / m_au.number_density_of(79).unwrap();

    // Rates add linearly.
    for i in 0..t_mix.energy_ev().len() {
        let want = w_si * t_si.inverse_mfp_per_m()[i] + w_au * t_au.inverse_mfp_per_m()[i];
        let got = t_mix.inverse_mfp_per_m()[i];
        assert!(
            (got - want).abs() <= 1e-12 * want,
            "row {i}: {got} vs {want}"
        );
    }

    // Angles: the mixture's forward CDF at its stored quantiles is the
    // number-density-weighted sum of the element CDFs, and equals the stored
    // probability.
    let (n_si, n_au) = (
        mix.number_density_of(14).unwrap(),
        mix.number_density_of(79).unwrap(),
    );
    for i in 0..t_mix.energy_ev().len() {
        let (rs, ra) = (&si().rows()[i], &au().rows()[i]);
        let total = n_si * rs.cumulative_total() + n_au * ra.cumulative_total();
        let q = t_mix.quantiles(i).unwrap();
        for (&u, &theta) in t_mix.probability().iter().zip(q) {
            let f = (n_si * rs.cumulative(theta) + n_au * ra.cumulative(theta)) / total;
            assert!((f - u).abs() <= 1e-12, "row {i}, u={u}: F={f}");
        }
    }

    // The element tables invert their own CDFs to the same precision.
    for (t, atom) in [(&t_si, si()), (&t_au, au())] {
        for (i, row) in atom.rows().iter().enumerate() {
            let q = t.quantiles(i).unwrap();
            for (&u, &theta) in t.probability().iter().zip(q) {
                let f = row.cumulative(theta) / row.cumulative_total();
                assert!(
                    (f - u).abs() <= 1e-12,
                    "Z={} row {i}, u={u}: F={f}",
                    atom.z()
                );
            }
        }
    }

    // And the mixture still recovers its own sigma_tr1.
    for i in 0..t_mix.energy_ev().len() {
        let (rs, ra) = (&si().rows()[i], &au().rows()[i]);
        let sel = n_si * rs.sigma_el() + n_au * ra.sigma_el();
        let str1 = n_si * rs.sigma_tr1() + n_au * ra.sigma_tr1();
        let rec = sel * mean_one_minus_cos(t_mix.probability(), t_mix.quantiles(i).unwrap());
        assert!((rec / str1 - 1.0).abs() < 1e-3, "row {i}");
    }
}

#[test]
fn cache_round_trip_checks_the_version() {
    let t = table_of(&pure(14), &[si().clone()]);
    assert_eq!(t.format_version(), CACHE_FORMAT_VERSION);
    let dir = std::env::temp_dir().join(format!("lindhard-elastic-table-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("si.toml");
    t.write_toml_file(&path).unwrap();
    let back = lindhard::electron::data::CrossSectionTable::from_toml_file(&path).unwrap();
    assert_eq!(back, t);

    let text = std::fs::read_to_string(&path).unwrap();
    let line = format!("format_version = {CACHE_FORMAT_VERSION}");
    assert!(text.contains(&line));
    let wrong = text.replace(
        &line,
        &format!("format_version = {}", CACHE_FORMAT_VERSION + 1),
    );
    match lindhard::electron::data::CrossSectionTable::from_toml_str(&wrong) {
        Err(ElectronDataError::UnsupportedVersion { found, supported }) => {
            assert_eq!(found, Some(CACHE_FORMAT_VERSION + 1));
            assert_eq!(supported, CACHE_FORMAT_VERSION);
        }
        other => panic!("expected a version error, got {other:?}"),
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn model_and_provenance_name_the_stand_in_potential() {
    let t = table_of(&pure(14), &[si().clone()]);
    assert!(t.model().contains("STAND-IN"), "{}", t.model());
    assert!(t.provenance().contains("STAND-IN"), "{}", t.provenance());
    assert!(!t.model().contains("DHFS atomic potential;"));
    assert!(t.material().contains("Si:1"), "{}", t.material());
}

#[test]
fn caller_chosen_two_point_grid() {
    let opts = ElasticTableOptions {
        energy_ev: vec![100.0, 1000.0],
        ..ElasticTableOptions::default()
    };
    let t = build_elastic_table(&pure(29), &ThomasFermiYukawa, &opts).unwrap();
    assert_eq!(t.energy_ev(), &[100.0, 1000.0]);
    assert!(t.inverse_mfp_per_m().iter().all(|r| *r > 0.0));

    // Without refinement the probability grid is stored as given.
    let p = vec![0.0, 0.25, 0.5, 0.75, 1.0];
    let opts = ElasticTableOptions {
        energy_ev: vec![100.0, 1000.0],
        probability: p.clone(),
        refine_tolerance: None,
        ..ElasticTableOptions::default()
    };
    let t = build_elastic_table(&pure(29), &ThomasFermiYukawa, &opts).unwrap();
    assert_eq!(t.probability(), p.as_slice());

    // One-point and unsorted grids are rejected.
    for bad in [vec![100.0], vec![1000.0, 100.0]] {
        let opts = ElasticTableOptions {
            energy_ev: bad,
            ..ElasticTableOptions::default()
        };
        assert!(matches!(
            build_elastic_table(&pure(29), &ThomasFermiYukawa, &opts),
            Err(ElasticTableError::Invalid { .. })
        ));
    }
}

/// A source that has no potential for one element.
struct MissingGold;

impl PotentialSource for MissingGold {
    fn potential(&self, z: u8) -> Result<Box<dyn ScreenedPotential>, ElasticError> {
        if z == 79 {
            Err(ElasticError::ScreeningCoefficientsUnavailable(79))
        } else {
            ThomasFermiYukawa.potential(z)
        }
    }
    fn description(&self) -> String {
        "test source without Au".into()
    }
}

#[test]
fn potential_errors_name_the_element() {
    let opts = ElasticTableOptions {
        energy_ev: vec![100.0, 200.0],
        ..ElasticTableOptions::default()
    };
    let mix = Material::from_atom_fractions(&[(14, 1.0), (79, 1.0)], Some(1.0e4)).unwrap();
    match build_elastic_table(&mix, &MissingGold, &opts) {
        Err(ElasticTableError::Potential { z: 79, .. }) => {}
        other => panic!("expected a potential error for Au, got {other:?}"),
    }
    // A zero-fraction component is skipped, so its missing potential is no error.
    let mix0 = Material::from_atom_fractions(&[(14, 1.0), (79, 0.0)], Some(2.3e3)).unwrap();
    assert!(build_elastic_table(&mix0, &MissingGold, &opts).is_ok());
    // The DHFS table is a documented gap.
    assert!(matches!(
        build_elastic_table(&pure(14), &SalvatDhfsTable, &opts),
        Err(ElasticTableError::Potential {
            z: 14,
            source: ElasticError::ScreeningCoefficientsUnavailable(14)
        })
    ));
}

#[test]
fn mixing_inconsistent_atomic_results_is_rejected() {
    let mix = Material::from_atom_fractions(&[(14, 1.0), (79, 1.0)], Some(1.0e4)).unwrap();
    let src = ThomasFermiYukawa;
    let pot = src.potential(79).unwrap();
    let other_grid = AtomicElastic::compute(
        79,
        pot.as_ref(),
        &src.description(),
        &[100.0, 200.0],
        SolverOptions::default(),
    )
    .unwrap();
    let r = combine(
        &mix,
        &[si().clone(), other_grid],
        &default_probability_grid(),
        None,
    );
    assert!(matches!(r, Err(ElasticTableError::Invalid { .. })));
    let r = combine(&mix, &[si().clone()], &default_probability_grid(), None);
    assert!(matches!(r, Err(ElasticTableError::Invalid { .. })));
}

#[test]
fn table_is_bit_identical_on_one_and_two_threads() {
    let mix = Material::from_atom_fractions(&[(6, 1.0), (79, 1.0)], Some(5.0e3)).unwrap();
    let opts = ElasticTableOptions {
        energy_ev: vec![20.0, 150.0, 1000.0],
        ..ElasticTableOptions::default()
    };
    let run = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| build_elastic_table(&mix, &ThomasFermiYukawa, &opts).unwrap())
    };
    let (a, b) = (run(1), run(2));
    assert_eq!(a, b);
    let bits = |t: &lindhard::electron::data::CrossSectionTable| -> Vec<u64> {
        (0..t.energy_ev().len())
            .flat_map(|i| t.quantiles(i).unwrap().to_vec())
            .chain(t.inverse_mfp_per_m().iter().copied())
            .chain(t.probability().iter().copied())
            .map(f64::to_bits)
            .collect()
    };
    assert_eq!(bits(&a), bits(&b));
}

/// The acceptance check on the full default grid (75 energies), for Si and
/// Au. Slow in a debug build; run with
/// `cargo test --release --test electron_elastic_table -- --ignored`.
#[test]
#[ignore]
fn full_default_grid_recovers_sigma_tr1() {
    let opts = ElasticTableOptions::default();
    let src = ThomasFermiYukawa;
    for z in [14u8, 79] {
        let pot = src.potential(z).unwrap();
        let atom = AtomicElastic::compute(
            z,
            pot.as_ref(),
            &src.description(),
            &opts.energy_ev,
            opts.solver,
        )
        .unwrap();
        let t = combine(
            &pure(z),
            std::slice::from_ref(&atom),
            &opts.probability,
            opts.refine_tolerance,
        )
        .unwrap();
        let (err, at) = worst_sigma_tr1_error(&atom, &t);
        println!(
            "Z={z}: {} energies, {} probability points, worst sigma_tr1 error {err:.3e} at {at} eV",
            t.energy_ev().len(),
            t.probability().len()
        );
        assert!(err < 1e-3, "Z={z}: {err:.3e} at {at} eV");
    }
}
