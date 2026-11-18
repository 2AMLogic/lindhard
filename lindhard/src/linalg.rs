//! Small dense linear algebra shared by the fitting code.
//!
//! One Gaussian-elimination solver and one inverse, so every caller (the PSF
//! fit, the Mermin fit, the stopping-correction fit) agrees on when a system
//! is numerically singular.

/// Relative pivot floor: a pivot no larger than this times the matrix scale
/// counts as zero.
const PIVOT_RTOL: f64 = 1e-12;

/// Solves the `n x n` system `m x = r` (`m` row-major) by Gaussian
/// elimination with partial pivoting; `None` if the system is singular or the
/// solution is not finite.
///
/// **Singularity criterion.** The system is singular when a pivot is `NaN` or
/// its magnitude is `<= 1e-12 * scale`, where `scale` is the largest absolute
/// entry of the (equilibrated, see below) matrix, computed once before
/// elimination. The test is relative so that a numerically singular matrix
/// (for example `[[1e10, 1e10], [1e-8, -1e-8]]`, condition number about
/// `1e18`) is rejected instead of returning a finite answer that does not
/// satisfy the equations. An all-zero matrix has `scale = 0` and is rejected.
///
/// **Equilibration.** When every diagonal entry is positive and finite (the
/// symmetric positive-semidefinite normal matrices `JᵀJ` the callers pass),
/// the system is first scaled symmetrically, `B = D A D` with
/// `D = diag(1/sqrt(a_ii))`, solved for `y`, and `x = D y` returned. This makes
/// the criterion independent of the units of the parameters: the PSF fit's
/// covariance matrix has columns spanning many orders of magnitude, and a
/// floor relative to its largest entry would reject its legitimately small
/// pivots. Otherwise the matrix is used as given.
pub(crate) fn solve(mut m: Vec<f64>, mut r: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    debug_assert_eq!(m.len(), n * n);
    debug_assert_eq!(r.len(), n);
    let equilibrate = (0..n).all(|i| m[i * n + i] > 0.0 && m[i * n + i].is_finite());
    let d: Vec<f64> = (0..n)
        .map(|i| {
            if equilibrate {
                1.0 / m[i * n + i].sqrt()
            } else {
                1.0
            }
        })
        .collect();
    if equilibrate {
        for i in 0..n {
            for k in 0..n {
                m[i * n + k] *= d[i] * d[k];
            }
            r[i] *= d[i];
        }
    }
    let scale = m.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    for c in 0..n {
        let piv = (c..n).max_by(|&i, &j| m[i * n + c].abs().total_cmp(&m[j * n + c].abs()))?;
        let p = m[piv * n + c].abs();
        if p.is_nan() || !p.is_finite() || p <= PIVOT_RTOL * scale {
            return None;
        }
        if piv != c {
            for k in 0..n {
                m.swap(c * n + k, piv * n + k);
            }
            r.swap(c, piv);
        }
        for i in c + 1..n {
            let f = m[i * n + c] / m[c * n + c];
            for k in c..n {
                m[i * n + k] -= f * m[c * n + k];
            }
            r[i] -= f * r[c];
        }
    }
    let mut x = vec![0.0; n];
    for c in (0..n).rev() {
        let s: f64 = (c + 1..n).map(|k| m[c * n + k] * x[k]).sum();
        x[c] = (r[c] - s) / m[c * n + c];
    }
    for (xi, di) in x.iter_mut().zip(&d) {
        *xi *= di;
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Inverse of an `n x n` matrix (row-major), column by column with
/// [`solve`], so it uses the same singularity criterion; `None` if singular.
pub(crate) fn invert(m: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    let mut inv = vec![0.0; n * n];
    for c in 0..n {
        let mut e = vec![0.0; n];
        e[c] = 1.0;
        let x = solve(m.clone(), e, n)?;
        for r in 0..n {
            inv[r * n + c] = x[r];
        }
    }
    Some(inv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solver_and_inverse() {
        let m = vec![4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0];
        let x = solve(m.clone(), vec![1.0, 2.0, 3.0], 3).unwrap();
        let r: Vec<f64> = (0..3)
            .map(|i| (0..3).map(|k| m[i * 3 + k] * x[k]).sum())
            .collect();
        for (a, b) in r.iter().zip([1.0, 2.0, 3.0]) {
            assert!((a - b).abs() < 1e-12);
        }
        let inv = invert(m.clone(), 3).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let s: f64 = (0..3).map(|k| m[i * 3 + k] * inv[k * 3 + j]).sum();
                assert!((s - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn exactly_singular_is_none() {
        assert!(solve(vec![1.0, 2.0, 2.0, 4.0], vec![1.0, 1.0], 2).is_none());
        assert!(invert(vec![1.0, 2.0, 2.0, 4.0], 2).is_none());
    }

    #[test]
    fn near_singular_is_none() {
        // Condition number ~1e18: a finite but wrong answer without the
        // relative pivot test.
        assert!(solve(vec![1e10, 1e10, 1e-8, -1e-8], vec![1.0, 1.0], 2).is_none());
    }

    #[test]
    fn uniformly_scaled_but_conditioned_solves() {
        // [[2, 1], [1, 3]] x = [3, 4] has x = [1, 1]; scaling the whole
        // system by 1e-9 or 1e9 must not change that (the test is relative,
        // not an absolute pivot floor).
        for s in [1e-9, 1e9] {
            let x = solve(vec![2.0 * s, s, s, 3.0 * s], vec![3.0 * s, 4.0 * s], 2).unwrap();
            assert!((x[0] - 1.0).abs() < 1e-9 && (x[1] - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn mixed_unit_normal_matrix_solves_and_scaled_singular_does_not() {
        // JᵀJ for parameters in different units: D C D with D = diag(d), so
        // entries span 1e16 .. 1e-16 although C is well conditioned. A
        // matrix-wide pivot floor would reject it (the PSF covariance, which
        // failed with Singular); the equilibrated floor must not.
        let d = [1e8, 1.0, 1e-8];
        let c = [4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0];
        let m: Vec<f64> = (0..9).map(|i| c[i] * d[i / 3] * d[i % 3]).collect();
        let inv = invert(m.clone(), 3).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let s: f64 = (0..3).map(|k| m[i * 3 + k] * inv[k * 3 + j]).sum();
                assert!((s - if i == j { 1.0 } else { 0.0 }).abs() < 1e-6);
            }
        }
        // The same congruence applied to an exactly singular C stays singular.
        let cs = [1.0, 2.0, 0.0, 2.0, 4.0, 0.0, 0.0, 0.0, 1.0];
        let ms: Vec<f64> = (0..9).map(|i| cs[i] * d[i / 3] * d[i % 3]).collect();
        assert!(invert(ms, 3).is_none());
    }

    #[test]
    fn edge_cases() {
        assert!(solve(vec![0.0; 4], vec![1.0, 1.0], 2).is_none());
        assert!(solve(vec![f64::NAN, 0.0, 0.0, 1.0], vec![1.0, 1.0], 2).is_none());
        assert_eq!(solve(vec![4.0], vec![2.0], 1).unwrap(), vec![0.5]);
        assert_eq!(solve(vec![], vec![], 0).unwrap(), Vec::<f64>::new());
    }
}
