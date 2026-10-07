use super::*;
use std::f64::consts::PI;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

#[test]
fn bessel_matches_closed_forms_and_wronskian() {
    for &x in &[0.05, 0.7, 3.3, 12.5, 60.0] {
        let (j, y) = spherical_bessel(6, x).unwrap();
        let (s, c) = x.sin_cos();
        let j0 = s / x;
        let j1 = s / (x * x) - c / x;
        let j2 = (3.0 / (x * x) - 1.0) * s / x - 3.0 * c / (x * x);
        let y0 = -c / x;
        let y1 = -c / (x * x) - s / x;
        let y2 = (-3.0 / (x * x) + 1.0) * c / x - 3.0 * s / (x * x);
        for (a, b) in [
            (j[0], j0),
            (j[1], j1),
            (j[2], j2),
            (y[0], y0),
            (y[1], y1),
            (y[2], y2),
        ] {
            assert!((a - b).abs() < 1e-12 * (1.0 + b.abs()), "x={x}: {a} vs {b}");
        }
        for n in 0..6 {
            let w = j[n + 1] * y[n] - j[n] * y[n + 1];
            assert!(rel(w, 1.0 / (x * x)) < 1e-11, "wronskian x={x} n={n}");
        }
    }
    // deep forbidden region reports None instead of overflowing
    assert!(spherical_bessel(400, 1.0).is_none());
}

#[test]
fn gauss_legendre_integrates_polynomials_exactly() {
    let (x, w) = gauss_legendre(9);
    let int: f64 = x.iter().zip(&w).map(|(x, w)| w * x.powi(16)).sum();
    assert!(rel(int, 2.0 / 17.0) < 1e-13);
    assert!((w.iter().sum::<f64>() - 2.0).abs() < 1e-13);
}

/// Closed-form Dirac phase shift of a square well: the interior solution is the
/// free solution with `E -> E + V0`, matched at the wall to the free exterior.
fn well_exact(v0: f64, rw: f64, e: f64, kappa: i32, c: f64) -> f64 {
    let l = if kappa < 0 { -kappa - 1 } else { kappa } as usize;
    let lbar = if kappa < 0 { l + 1 } else { l - 1 };
    let sigma = if kappa < 0 { -1.0 } else { 1.0 };
    let ei = e + v0;
    let qi = (ei * (ei + 2.0 * c * c)).sqrt() / c;
    let si = c * qi / (ei + 2.0 * c * c);
    let k = (e * (e + 2.0 * c * c)).sqrt() / c;
    let s = c * k / (e + 2.0 * c * c);
    let (ji, _) = spherical_bessel(l.max(lbar), qi * rw).unwrap();
    let (jo, yo) = spherical_bessel(l.max(lbar), k * rw).unwrap();
    let (p, q) = (ji[l], sigma * si * ji[lbar]);
    let (jp, jq) = (jo[l], sigma * s * jo[lbar]);
    let (yp, yq) = (yo[l], sigma * s * yo[lbar]);
    ((q * jp - p * jq) / (q * yp - p * yq)).atan()
}

#[test]
fn square_well_matches_closed_form() {
    let c = 1.0 / FINE_STRUCTURE;
    let mut worst: f64 = 0.0;
    for &(v0, rw) in &[(0.5, 2.0), (3.0, 1.5), (20.0, 1.0)] {
        let well = SquareWell::new(v0, rw).unwrap();
        for &e_ev in &[50.0, 500.0, 5000.0] {
            let solver = ElasticSolver::new(&well, e_ev, SolverOptions::default()).unwrap();
            let e = e_ev / HARTREE_EV;
            for kappa in [-1, 1, -2, 2, -3, 3, -5, 6] {
                let num = solver.phase_shift(kappa).unwrap();
                let exact = well_exact(v0, rw, e, kappa, c);
                // phase shifts are defined mod pi; compare tan or the reduced value
                let d = (num - exact).abs();
                let scale = exact.abs().max(1e-3);
                worst = worst.max(d / scale);
                assert!(
                    d < 1e-6 * scale,
                    "V0={v0} R={rw} E={e_ev} kappa={kappa}: {num} vs {exact}"
                );
            }
        }
    }
    eprintln!("worst relative well error {worst:e}");
}

fn tf_yukawa(z: f64) -> Yukawa {
    // Thomas-Fermi length 0.8853 a0 / Z^(1/3): a fixture screening, not the
    // Salvat 1987 function (whose table is unavailable).
    Yukawa::new(z, 0.8853 / z.cbrt()).unwrap()
}

#[test]
fn optical_theorem_and_quadrature_agree() {
    for &z in &[6.0, 29.0, 79.0] {
        let pot = tf_yukawa(z);
        for &e in &[100.0, 1000.0, 10_000.0] {
            let s = ElasticSolver::new(&pot, e, SolverOptions::default()).unwrap();
            let pw = s.partial_waves().unwrap();
            let opt = pw.sigma_optical();
            let sum = pw.sigma_el();
            let quad = pw.sigma_el_quadrature();
            assert!(
                rel(quad, opt) < 1e-4,
                "Z={z} E={e}: quad {quad} optical {opt}"
            );
            assert!(
                rel(sum, opt) < 1e-10,
                "Z={z} E={e}: sum {sum} optical {opt}"
            );
            for t in [0.1, 1.0, 2.5, 3.1] {
                let sh = pw.sherman(t);
                assert!(sh.abs() <= 1.0 + 1e-12, "|S|<=1 violated: {sh}");
            }
        }
    }
}

#[test]
fn raising_l_max_by_20_percent_changes_nothing() {
    for &(z, e) in &[(6.0, 100.0), (29.0, 1000.0), (79.0, 10_000.0)] {
        let pot = tf_yukawa(z);
        let s = ElasticSolver::new(&pot, e, SolverOptions::default()).unwrap();
        let base = s.partial_waves().unwrap();
        let more = s
            .phase_shifts_up_to((base.l_max() as f64 * 1.2).ceil() as usize)
            .unwrap();
        assert!(more.l_max() > base.l_max());
        assert!(rel(more.sigma_el(), base.sigma_el()) < 1e-5);
        assert!(rel(more.sigma_tr1(), base.sigma_tr1()) < 1e-5);
    }
}

#[test]
fn step_size_convergence() {
    let pot = tf_yukawa(29.0);
    let coarse = SolverOptions::default();
    let fine = SolverOptions {
        step_scale: coarse.step_scale / 2.0,
        ..coarse
    };
    let a = ElasticSolver::new(&pot, 1000.0, coarse)
        .unwrap()
        .partial_waves()
        .unwrap();
    let b = ElasticSolver::new(&pot, 1000.0, fine)
        .unwrap()
        .partial_waves()
        .unwrap();
    assert!(rel(a.sigma_el(), b.sigma_el()) < 1e-7);
    assert!(rel(a.sigma_tr1(), b.sigma_tr1()) < 1e-7);
}

/// First-Born (screened-Rutherford) cross sections of `V = -Z e^{-r/a}/r` for a
/// Dirac electron: `DCS = gamma^2 Z^2 4/(q^2 + mu^2)^2 (1 - beta^2 sin^2(t/2))`
/// with `q = 2 k sin(t/2)`, `mu = 1/a` (McKinley-Feshbach lowest order),
/// integrated in `s = sin^2(t/2)`: `d Omega = 4 pi ds`, `1 - cos t = 2 s`.
fn born(z: f64, a: f64, e_ev: f64) -> (f64, f64) {
    let c = 1.0 / FINE_STRUCTURE;
    let e = e_ev / HARTREE_EV;
    let k = (e * (e + 2.0 * c * c)).sqrt() / c;
    let gamma = 1.0 + e / (c * c);
    let beta2 = 1.0 - 1.0 / (gamma * gamma);
    let (aa, bb) = (4.0 * k * k, 1.0 / (a * a));
    let pref = 16.0 * PI * gamma * gamma * z * z;
    let i0 = (1.0 / bb - 1.0 / (aa + bb)) / aa;
    let l = ((aa + bb) / bb).ln();
    let i1 = (l + bb / (aa + bb) - 1.0) / (aa * aa);
    let i2 = (1.0 - 2.0 * bb / aa * l + bb * bb / aa * (1.0 / bb - 1.0 / (aa + bb))) / (aa * aa);
    (pref * (i0 - beta2 * i1), 2.0 * pref * (i1 - beta2 * i2))
}

#[test]
fn born_closed_forms_match_brute_force() {
    let (z, a, e) = (1.0, 0.2, 5000.0);
    let c = 1.0 / FINE_STRUCTURE;
    let en = e / HARTREE_EV;
    let k = (en * (en + 2.0 * c * c)).sqrt() / c;
    let gamma = 1.0 + en / (c * c);
    let beta2 = 1.0 - 1.0 / (gamma * gamma);
    let n = 200_000;
    let (mut s0, mut s1) = (0.0, 0.0);
    for i in 0..n {
        // midpoint in u = ln-spaced s to resolve the forward peak
        let u0 = (1e-12_f64).ln();
        let u = u0 + (0.0 - u0) * (i as f64 + 0.5) / n as f64;
        let s = u.exp();
        let ds = s * (0.0 - u0) / n as f64;
        let d = gamma * gamma * z * z * 4.0 / (4.0 * k * k * s + 1.0 / (a * a)).powi(2)
            * (1.0 - beta2 * s);
        s0 += d * 4.0 * PI * ds;
        s1 += d * 2.0 * s * 4.0 * PI * ds;
    }
    let (b0, b1) = born(z, a, e);
    assert!(
        rel(s0, b0) < 1e-5 && rel(s1, b1) < 1e-4,
        "{s0} {b0} {s1} {b1}"
    );
}

/// High-energy limit. For `Z = 1` at 50 keV the first-Born parameter is
/// `Z alpha / beta = 0.0073/0.4127 = 0.018` (beta = 0.4127 at 50 keV), so the
/// second Born term, which is of relative order `pi Z alpha / beta` times an
/// angular factor, is below a few percent: hence the 2 % tolerance on
/// `sigma_el` and `sigma_tr1`. The gap to the first-Born value must also shrink
/// monotonically from 5 keV (`Z alpha/beta = 0.05`) to 50 keV.
#[test]
fn yukawa_high_energy_approaches_first_born() {
    let a = 0.2;
    let pot = Yukawa::new(1.0, a).unwrap();
    let mut prev = f64::INFINITY;
    for &e in &[5_000.0, 10_000.0, 20_000.0, 50_000.0] {
        let s = ElasticSolver::new(&pot, e, SolverOptions::default()).unwrap();
        let pw = s.partial_waves().unwrap();
        let (b_el, b_tr) = born(1.0, a, e);
        let g_el = rel(pw.sigma_el(), b_el);
        let g_tr = rel(pw.sigma_tr1(), b_tr);
        let gap = g_el.max(g_tr);
        eprintln!("E={e} l_max={} gap_el={g_el:e} gap_tr={g_tr:e}", pw.l_max());
        assert!(gap < prev, "gap did not shrink at {e} eV: {gap} >= {prev}");
        prev = gap;
        if e == 50_000.0 {
            assert!(g_el < 0.02 && g_tr < 0.02, "{g_el} {g_tr}");
        }
    }
}

#[test]
fn thread_count_does_not_change_bits() {
    let pot = tf_yukawa(29.0);
    let run = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| {
                ElasticSolver::new(&pot, 1000.0, SolverOptions::default())
                    .unwrap()
                    .partial_waves()
                    .unwrap()
            })
    };
    let (a, b, c) = (run(1), run(2), run(8));
    assert_eq!(a, b);
    assert_eq!(a, c);
    assert_eq!(a.sigma_tr1().to_bits(), c.sigma_tr1().to_bits());
}

#[test]
fn solve_reports_observables_on_a_grid() {
    let pot = tf_yukawa(6.0);
    let thetas = [0.0, 0.5, 1.0, 2.0, PI];
    let r = solve(&pot, 1000.0, &thetas, SolverOptions::default()).unwrap();
    assert_eq!(r.dcs.len(), 5);
    assert!(r.dcs.iter().all(|d| *d > 0.0));
    assert!(r.sigma_tr1 < r.sigma_el * 2.0 && r.sigma_tr1 > 0.0);
    assert!(solve(&pot, 1000.0, &[4.0], SolverOptions::default()).is_err());
    assert!(solve(&pot, -1.0, &thetas, SolverOptions::default()).is_err());
}

#[test]
fn salvat_table_is_a_documented_gap() {
    assert_eq!(
        SalvatDhfs::for_element(29),
        Err(ElasticError::ScreeningCoefficientsUnavailable(29))
    );
    assert!(SalvatDhfs::from_coefficients(29, [0.5, 0.3, 0.3], [1.0, 2.0, 3.0]).is_err());
    let p = SalvatDhfs::from_coefficients(1, [0.5, 0.3, 0.2], [1.0, 2.0, 3.0]).unwrap();
    assert!(p.energy(1.0) < 0.0);
}

// ---------------------------------------------------------------------------
// Exchange and correlation-polarization corrections (#91)

use super::corrections::{
    correlation_potential, exchange_potential, CorrelationPolarizationInfo,
    CORRELATION_POLARIZATION_MODEL, DEFAULT_OUTER_RADIUS, EXCHANGE_MODEL,
};

/// Synthetic polarizability for the fixture tests (bohr^3); not a property of
/// any element.
const FIXTURE_ALPHA: f64 = 40.0;

fn fixture_cp(outer_radius: f64) -> CorrelationPolarization {
    CorrelationPolarization {
        outer_radius,
        ..CorrelationPolarization::new(FIXTURE_ALPHA, "synthetic fixture value")
    }
}

fn both(outer_radius: f64) -> Corrections {
    Corrections {
        exchange: true,
        correlation_polarization: Some(fixture_cp(outer_radius)),
    }
}

#[test]
fn correlation_potential_matches_independent_forms() {
    // Salvat (2003) Eq. (7a), with the sign discussed in the module docs,
    // against the Perdew-Zunger A, B, C, D form (dsedu.org page, see docs).
    let (a, b, c, d) = (0.0311, -0.048, 0.002, -0.0116);
    let pz_high = |rs: f64| {
        a * rs.ln() + (b - a / 3.0) + 2.0 / 3.0 * c * rs * rs.ln() + (2.0 * d - c) / 3.0 * rs
    };
    let rho_of = |rs: f64| 3.0 / (4.0 * PI * rs.powi(3));
    for rs in [1e-3, 0.05, 0.3, 0.7, 0.99] {
        let v = correlation_potential(rho_of(rs));
        assert!(v < 0.0, "rs={rs}: {v}");
        // Salvat prints the coefficients rounded to 3-4 decimals.
        assert!(
            (v - pz_high(rs)).abs() < 2e-4,
            "rs={rs}: {v} vs {}",
            pz_high(rs)
        );
    }
    // Eq. (7b) at rs = 4, evaluated by hand from the printed constants:
    // -0.1423 (1 + 7/6 1.0529 * 2 + 4/3 0.3334 * 4)/(1 + 1.0529 * 2 + 0.3334 * 4)^2
    let num = 1.0 + 7.0 / 6.0 * 1.0529 * 2.0 + 4.0 / 3.0 * 0.3334 * 4.0;
    let den = (1.0 + 1.0529 * 2.0 + 0.3334 * 4.0_f64).powi(2);
    assert!((correlation_potential(rho_of(4.0)) + 0.1423 * num / den).abs() < 1e-15);
    // The two branches join at rs = 1 (to the rounding of the coefficients).
    let below = correlation_potential(rho_of(1.0 - 1e-12));
    let above = correlation_potential(rho_of(1.0 + 1e-12));
    assert!((below - above).abs() < 1e-4, "{below} vs {above}");
    assert!((above + 0.0668).abs() < 1e-4, "{above}");
    // Zero-density limit and invalid densities.
    assert_eq!(correlation_potential(0.0), 0.0);
    assert_eq!(correlation_potential(-1.0), 0.0);
    assert_eq!(correlation_potential(f64::NAN), 0.0);
    assert!(correlation_potential(1e-30).abs() < 1e-9);
}

#[test]
fn exchange_matches_the_printed_formula() {
    // Salvat (2003) Eq. (2) as printed, where it has no cancellation.
    let printed = |d: f64, rho: f64| 0.5 * d - 0.5 * (d * d + 4.0 * PI * rho).sqrt();
    for (d, rho) in [(0.5, 1.0), (3.0, 0.2), (30.0, 5.0), (0.0, 0.3), (-0.5, 0.1)] {
        let v = exchange_potential(d, rho);
        assert!(rel(v, printed(d, rho)) < 1e-13, "{d} {rho}");
        assert!(v <= 0.0);
    }
    // No density, no exchange.
    assert_eq!(exchange_potential(3.0, 0.0), 0.0);
    // High-energy limit -pi rho / D, where the printed form cancels.
    let (d, rho) = (367.0, 1e-9);
    assert!(rel(exchange_potential(d, rho), -PI * rho / d) < 1e-12);
}

#[test]
fn polarization_tail_join_and_outer_cut() {
    let pot = tf_yukawa(29.0);
    let e = 1000.0;
    let c = Corrections {
        exchange: false,
        correlation_polarization: Some(fixture_cp(30.0)),
    };
    let cp = CorrectedPotential::new(&pot, &pot, e, &c).unwrap();
    let info: &CorrelationPolarizationInfo =
        cp.metadata().correlation_polarization.as_ref().unwrap();
    // Eqs. (4)-(5): b^2 = (E - 50 eV)/16 eV, d^4 = alpha Z^(-1/3) b^2 / 2.
    let b2 = (e - 50.0) / 16.0;
    assert!(rel(info.b_pol_squared, b2) < 1e-15);
    let d4 = 0.5 * FIXTURE_ALPHA * 29.0_f64.powf(-1.0 / 3.0) * b2;
    assert!(rel(info.cutoff_d_bohr.powi(4), d4) < 1e-12);
    // Eq. (3) and its -alpha/(2 r^4) tail.
    for r in [0.5, 3.0, 12.0] {
        let want = -FIXTURE_ALPHA / (2.0 * (r * r + d4.sqrt()).powi(2));
        assert!(rel(cp.polarization_energy(r), want) < 1e-12);
    }
    let r = 1.0e3;
    assert!(
        rel(
            cp.polarization_energy(r),
            -FIXTURE_ALPHA / (2.0 * r.powi(4))
        ) < 1e-4
    );
    // Eq. (9): the outer crossing, V_pol beyond it, max below it, 0 past the cut.
    let r_cp = info
        .join_radius_bohr
        .expect("crossing above the scan floor");
    let v_co = |r: f64| correlation_potential(pot.density(r));
    assert!(rel(v_co(r_cp), cp.polarization_energy(r_cp)) < 1e-9);
    for r in [r_cp * 1.01, 0.5 * (r_cp + 30.0), 29.9] {
        assert_eq!(
            cp.correlation_polarization_energy(r),
            cp.polarization_energy(r)
        );
    }
    for r in [1e-3, 0.1, 1.0, r_cp * 0.99] {
        let want = v_co(r).max(cp.polarization_energy(r));
        assert_eq!(cp.correlation_polarization_energy(r), want);
    }
    assert_eq!(cp.correlation_polarization_energy(30.0), 0.0);
    assert_eq!(cp.correlation_polarization_energy(31.0), 0.0);
    // The cut is a breakpoint and the matching radius lies beyond it.
    assert!(cp.breakpoints().contains(&30.0));
    assert!(cp.matching_radius(1e-10) >= 30.0);
    // The total adds the pieces.
    let r = 0.7;
    let want = pot.energy(r) + cp.correlation_polarization_energy(r);
    assert_eq!(cp.energy(r), want);
}

#[test]
fn poisson_densities_integrate_to_z() {
    // Salvat et al. (1987) Eq. (12) and its Yukawa special case hold Z electrons.
    let dhfs = SalvatDhfs::from_coefficients(10, [0.2, 0.5, 0.3], [9.0, 2.0, 0.8]).unwrap();
    let yuk = tf_yukawa(10.0);
    for d in [&dhfs as &dyn ElectronDensity, &yuk] {
        // midpoint rule in u = ln r
        let (u0, u1, n) = ((1e-8_f64).ln(), (200.0_f64).ln(), 400_000);
        let h = (u1 - u0) / n as f64;
        let s: f64 = (0..n)
            .map(|i| {
                let r = (u0 + (i as f64 + 0.5) * h).exp();
                4.0 * PI * r * r * d.density(r) * r * h
            })
            .sum();
        assert!(rel(s, 10.0) < 1e-6, "{s}");
    }
}

#[test]
fn corrected_potential_rejects_a_mismatched_solver_energy() {
    let pot = tf_yukawa(6.0);
    let c = Corrections {
        exchange: true,
        correlation_polarization: Some(fixture_cp(20.0)),
    };
    let corrected = CorrectedPotential::new(&pot, &pot, 100.0, &c).unwrap();
    let opts = SolverOptions::default();
    assert!(matches!(
        ElasticSolver::new(&corrected, 10_000.0, opts),
        Err(ElasticError::Invalid { .. })
    ));
    assert!(solve(&corrected, 10_000.0, &[], opts).is_err());
    assert!(ElasticSolver::new(&corrected, 100.0, opts).is_ok());
    assert!(corrected.solver(opts).is_ok());
    // An unbound potential accepts any energy.
    assert!(ElasticSolver::new(&pot, 10_000.0, opts).is_ok());
}

#[test]
fn trait_based_solve_of_a_corrected_potential_records_its_corrections() {
    let pot = tf_yukawa(6.0);
    let c = Corrections {
        exchange: true,
        correlation_polarization: Some(fixture_cp(20.0)),
    };
    let corrected = CorrectedPotential::new(&pot, &pot, 500.0, &c).unwrap();
    let r = solve(
        &corrected,
        corrected.energy_ev(),
        &[],
        SolverOptions::default(),
    )
    .unwrap();
    assert_eq!(r.corrections, *corrected.metadata());
    assert!(r.corrections.exchange.is_some());
    assert!(r.corrections.correlation_polarization.is_some());
}

#[test]
fn every_option_combination_is_recorded() {
    let pot = tf_yukawa(6.0);
    let thetas = [0.0, 1.0, PI];
    let opts = SolverOptions::default();
    let combos = [
        (false, None),
        (true, None),
        (false, Some(fixture_cp(20.0))),
        (true, Some(fixture_cp(20.0))),
    ];
    let mut sigmas = Vec::new();
    for (exchange, cp) in combos {
        let c = Corrections {
            exchange,
            correlation_polarization: cp.clone(),
        };
        let r = solve_corrected(&pot, &pot, 500.0, &c, &thetas, opts).unwrap();
        let m = &r.corrections;
        assert_eq!(m.exchange.is_some(), exchange);
        if exchange {
            assert_eq!(m.exchange.as_ref().unwrap().model, EXCHANGE_MODEL);
        }
        assert_eq!(m.correlation_polarization.is_some(), cp.is_some());
        if let Some(info) = &m.correlation_polarization {
            assert_eq!(info.model, CORRELATION_POLARIZATION_MODEL);
            assert_eq!(info.polarizability_bohr3, FIXTURE_ALPHA);
            assert_eq!(info.polarizability_source, "synthetic fixture value");
            assert_eq!(info.outer_radius_bohr, 20.0);
        }
        assert_eq!(m.density_source.is_some(), exchange || cp.is_some());
        if !exchange && cp.is_none() {
            assert_eq!(*m, CorrectionMetadata::default());
        }
        sigmas.push(r.sigma_el);
    }
    // Each correction changes the answer.
    for i in 1..4 {
        assert!(rel(sigmas[i], sigmas[0]) > 1e-4, "{sigmas:?}");
    }
    // Plain solve reports both off.
    let r = solve(&pot, 500.0, &thetas, opts).unwrap();
    assert_eq!(r.corrections, CorrectionMetadata::default());
    assert!(Corrections::default().is_none());
}

#[test]
fn invalid_correction_inputs_are_rejected() {
    let pot = tf_yukawa(29.0);
    let mk = |cp: CorrelationPolarization, e: f64| {
        let c = Corrections {
            exchange: false,
            correlation_polarization: Some(cp),
        };
        CorrectedPotential::new(&pot, &pot, e, &c).map(|_| ())
    };
    assert!(mk(fixture_cp(30.0), 1000.0).is_ok());
    for alpha in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let cp = CorrelationPolarization {
            polarizability: alpha,
            ..fixture_cp(30.0)
        };
        assert!(mk(cp, 1000.0).is_err(), "alpha={alpha}");
    }
    let cp = CorrelationPolarization {
        polarizability_source: "  ".into(),
        ..fixture_cp(30.0)
    };
    assert!(mk(cp, 1000.0).is_err());
    // Seltzer's recipe needs E > 50 eV.
    assert!(mk(fixture_cp(30.0), 50.0).is_err());
    assert!(mk(fixture_cp(30.0), 40.0).is_err());
    for b2 in [-1.0, f64::NAN] {
        let cp = CorrelationPolarization {
            cutoff: PolarizationCutoff::BPolSquared(b2),
            ..fixture_cp(30.0)
        };
        assert!(mk(cp, 1000.0).is_err());
    }
    let cp = CorrelationPolarization {
        cutoff: PolarizationCutoff::BPolSquared(2.0),
        ..fixture_cp(30.0)
    };
    assert!(mk(cp, 40.0).is_ok());
    for r in [0.0, -3.0, f64::NAN] {
        assert!(mk(fixture_cp(r), 1000.0).is_err());
    }
    // An outer radius inside the region where V_co < V_pol has no outer crossing.
    assert!(mk(fixture_cp(0.05), 1000.0).is_err());
    for e in [0.0, -5.0, f64::NAN] {
        assert!(CorrectedPotential::new(&pot, &pot, e, &both(30.0)).is_err());
    }
    // Eq. (4) needs Z: a regular potential has none.
    let well = SquareWell::new(1.0, 2.0).unwrap();
    let c = Corrections {
        exchange: false,
        correlation_polarization: Some(fixture_cp(30.0)),
    };
    assert!(CorrectedPotential::new(&well, &pot, 1000.0, &c).is_err());
}

fn assert_same_bits(a: &ElasticResult, b: &ElasticResult) {
    assert_eq!(a.l_max, b.l_max);
    for l in 0..=a.l_max {
        for kappa in [-(l as i32) - 1, l as i32] {
            if kappa == 0 {
                continue;
            }
            let (x, y) = (
                a.partial_waves.phase_shift(kappa).unwrap(),
                b.partial_waves.phase_shift(kappa).unwrap(),
            );
            assert_eq!(x.to_bits(), y.to_bits(), "kappa={kappa}");
        }
    }
    for (x, y) in a
        .dcs
        .iter()
        .zip(&b.dcs)
        .chain(a.sherman.iter().zip(&b.sherman))
    {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    assert_eq!(a.sigma_el.to_bits(), b.sigma_el.to_bits());
    assert_eq!(a.sigma_tr1.to_bits(), b.sigma_tr1.to_bits());
    assert_eq!(
        a.partial_waves.sigma_optical().to_bits(),
        b.partial_waves.sigma_optical().to_bits()
    );
}

#[test]
fn corrections_off_is_bit_identical_to_the_plain_solver() {
    let thetas = [0.0, 0.2, 1.0, 2.5, PI];
    let opts = SolverOptions::default();
    let dhfs_fixture = SalvatDhfs::from_coefficients(29, [0.2, 0.5, 0.3], [6.0, 2.0, 0.8]).unwrap();
    for &e in &[100.0, 1000.0] {
        for pot in [&tf_yukawa(29.0) as &dyn ScreenedPotential, &dhfs_fixture] {
            let plain = solve(pot, e, &thetas, opts).unwrap();
            let off = solve_corrected(pot, &dhfs_fixture, e, &Corrections::none(), &thetas, opts)
                .unwrap();
            assert_same_bits(&plain, &off);
            assert_eq!(off.corrections, CorrectionMetadata::default());
            // A CorrectedPotential with nothing switched on is the static one.
            let wrapped =
                CorrectedPotential::new(pot, &dhfs_fixture, e, &Corrections::none()).unwrap();
            assert!(!wrapped.long_range());
            let via = solve(&wrapped, e, &thetas, opts).unwrap();
            assert_same_bits(&plain, &via);
        }
    }
}

#[test]
fn corrected_thread_count_does_not_change_bits() {
    let pot = tf_yukawa(29.0);
    let run = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| {
                solve_corrected(
                    &pot,
                    &pot,
                    1000.0,
                    &both(20.0),
                    &[0.0, 0.5, 2.0],
                    SolverOptions::default(),
                )
                .unwrap()
            })
    };
    let (a, b, c) = (run(1), run(2), run(8));
    assert_same_bits(&a, &b);
    assert_same_bits(&a, &c);
    assert_eq!(a, c);
}

#[test]
fn optical_theorem_holds_with_corrections() {
    for &z in &[6.0, 29.0, 79.0] {
        let pot = tf_yukawa(z);
        for &e in &[100.0, 1000.0, 10_000.0] {
            for c in [
                Corrections {
                    exchange: true,
                    correlation_polarization: None,
                },
                both(20.0),
            ] {
                let cp = CorrectedPotential::new(&pot, &pot, e, &c).unwrap();
                let pw = cp
                    .solver(SolverOptions::default())
                    .unwrap()
                    .partial_waves()
                    .unwrap();
                let (opt, sum, quad) =
                    (pw.sigma_optical(), pw.sigma_el(), pw.sigma_el_quadrature());
                assert!(
                    rel(quad, opt) < 1e-4,
                    "Z={z} E={e}: quad {quad} optical {opt}"
                );
                assert!(
                    rel(sum, opt) < 1e-10,
                    "Z={z} E={e}: sum {sum} optical {opt}"
                );
            }
        }
    }
}

/// Numerical convergence of the corrected potential, separately from the
/// optical theorem: halving the step, tightening the matching threshold and
/// the partial-wave tolerance, and moving the outer cut of the polarization
/// tail outward.
#[test]
fn corrected_numerics_converge() {
    let pot = tf_yukawa(29.0);
    let e = 1000.0;
    let base_opts = SolverOptions::default();
    let sig = |c: &Corrections, o: SolverOptions| {
        let r = solve_corrected(&pot, &pot, e, c, &[], o).unwrap();
        (r.sigma_el, r.sigma_tr1)
    };
    let base = sig(&both(20.0), base_opts);
    let off = sig(&Corrections::none(), base_opts);
    let effect = (base.0 - off.0).abs();
    for o in [
        SolverOptions {
            step_scale: base_opts.step_scale / 2.0,
            ..base_opts
        },
        SolverOptions {
            matching_threshold: 1e-14,
            ..base_opts
        },
        SolverOptions {
            phase_tolerance: 1e-10,
            ..base_opts
        },
    ] {
        let t = sig(&both(20.0), o);
        assert!(rel(t.0, base.0) < 1e-7, "{o:?}: {t:?} vs {base:?}");
        assert!(rel(t.1, base.1) < 1e-7, "{o:?}: {t:?} vs {base:?}");
    }
    // The tail beyond the cut: doubling it moves sigma_el by a small fraction
    // of the size of the correction.
    let far = sig(&both(40.0), base_opts);
    assert!(
        (far.0 - base.0).abs() < 0.01 * effect,
        "{far:?} vs {base:?}, effect {effect}"
    );
    assert!(rel(far.1, base.1) < 1e-5, "{far:?} vs {base:?}");
}

/// Algorithm check on the fixture (not Cu data): each correction changes
/// sigma_el at 10 keV by less than a fifth of its change at 100 eV.
#[test]
fn corrections_fade_with_energy_on_a_fixture() {
    let pot = tf_yukawa(29.0);
    let opts = SolverOptions::default();
    let delta = |c: &Corrections, e: f64| {
        let a = solve_corrected(&pot, &pot, e, c, &[], opts)
            .unwrap()
            .sigma_el;
        let b = solve(&pot, e, &[], opts).unwrap().sigma_el;
        (a - b).abs()
    };
    for c in [
        Corrections {
            exchange: true,
            correlation_polarization: None,
        },
        Corrections {
            exchange: false,
            correlation_polarization: Some(fixture_cp(20.0)),
        },
    ] {
        let (lo, hi) = (delta(&c, 100.0), delta(&c, 10_000.0));
        assert!(hi < lo / 5.0, "{c:?}: {lo} at 100 eV, {hi} at 10 keV");
    }
}

/// The Cu evidence of #91, with sourced inputs: Salvat et al. (1987) Table I
/// screening coefficients for Cu, its Eq. (12) density, and the Cu dipole
/// polarizability 47(1) a.u. of Schwerdtfeger & Nagle (2019), Fig. 1 (see
/// `docs/data-provenance.md`). Ignored by default (about a minute in release,
/// far longer in a debug build); run with
/// `cargo test -p lindhard --release --lib cu_corrections_fade_with_energy -- --ignored --nocapture`.
#[test]
#[ignore]
fn cu_corrections_fade_with_energy() {
    let cu = SalvatDhfs::from_coefficients(29, [0.0771, 0.7951, 0.1278], [25.326, 3.3928, 1.1426])
        .unwrap();
    let cp = |outer: f64| CorrelationPolarization {
        outer_radius: outer,
        ..CorrelationPolarization::new(
            47.0,
            "Schwerdtfeger & Nagle, Mol. Phys. 117, 1200 (2019), Fig. 1: Cu 47(1) a.u.",
        )
    };
    let opts = SolverOptions::default();
    eprintln!("Cu: Salvat 1987 Table I, Eq. (12) density, alpha_d = 47 bohr^3; {opts:?}");
    let models = [
        (
            "exchange",
            Corrections {
                exchange: true,
                correlation_polarization: None,
            },
        ),
        (
            "corr-pol",
            Corrections {
                exchange: false,
                correlation_polarization: Some(cp(DEFAULT_OUTER_RADIUS)),
            },
        ),
    ];
    let mut d = [[0.0; 2]; 2];
    for (j, &e) in [100.0, 10_000.0].iter().enumerate() {
        let base = solve(&cu, e, &[], opts).unwrap();
        for (i, (name, c)) in models.iter().enumerate() {
            let r = solve_corrected(&cu, &cu, e, c, &[], opts).unwrap();
            let diff = r.sigma_el - base.sigma_el;
            d[i][j] = diff;
            eprintln!(
                "E = {e:>7} eV  {name:<8}  sigma_el: off {:.8e}  on {:.8e} bohr^2  abs {:+.4e}  frac {:+.4e}  (l_max {} -> {}, {:?})",
                base.sigma_el,
                r.sigma_el,
                diff,
                diff / base.sigma_el,
                base.l_max,
                r.l_max,
                r.corrections.correlation_polarization.as_ref().map(|m| (m.b_pol_squared, m.cutoff_d_bohr, m.join_radius_bohr)),
            );
        }
        // Tail convergence of the polarization correction.
        let far = Corrections {
            exchange: false,
            correlation_polarization: Some(cp(2.0 * DEFAULT_OUTER_RADIUS)),
        };
        let r = solve_corrected(&cu, &cu, e, &far, &[], opts).unwrap();
        eprintln!(
            "E = {e:>7} eV  corr-pol, outer radius {} bohr: abs change {:+.4e}",
            2.0 * DEFAULT_OUTER_RADIUS,
            r.sigma_el - base.sigma_el
        );
    }
    for (i, (name, _)) in models.iter().enumerate() {
        let ratio = d[i][1].abs() / d[i][0].abs();
        eprintln!("{name}: |change at 10 keV| / |change at 100 eV| = {ratio:.4e}");
        assert!(ratio < 0.2, "{name}: {ratio}");
    }
}
