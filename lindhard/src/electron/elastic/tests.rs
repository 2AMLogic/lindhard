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
