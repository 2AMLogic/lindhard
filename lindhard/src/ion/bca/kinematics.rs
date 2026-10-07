//! Two-body elastic kinematics, direction rotation and planar surface
//! refraction. Pure functions, no state.

use std::f64::consts::PI;

/// Laboratory deflection of projectile and recoil for a centre-of-mass
/// scattering angle `theta` and mass ratio `mu = M1 / M2`.
///
/// Classical elastic two-body kinematics (e.g. H. Goldstein, *Classical
/// Mechanics*, 2nd ed. (1980), Sec. 3.11; M. T. Robinson and I. M. Torrens,
/// Phys. Rev. B 9 (1974) 5008):
///
/// * projectile: `tan(psi) = sin(theta) / (cos(theta) + mu)`, `psi` in `[0, pi]`,
/// * recoil: `phi = (pi - theta) / 2`, on the opposite azimuth,
/// * energy transfer: `T = gamma E sin^2(theta/2)`,
///   `gamma = 4 M1 M2 / (M1 + M2)^2`.
///
/// Returns `(psi, phi)`.
pub fn lab_angles(theta: f64, mu: f64) -> (f64, f64) {
    let psi = theta.sin().atan2(theta.cos() + mu);
    (psi, 0.5 * (PI - theta))
}

/// Kinematic factor `gamma = 4 M1 M2 / (M1 + M2)^2`.
pub fn gamma(m1: f64, m2: f64) -> f64 {
    4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2))
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalise a vector to unit length.
pub fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

/// Rotate unit vector `u` by polar angle `polar` (from `u`) at azimuth
/// `azimuth` about `u`. The azimuth is measured in a frame built from `u` and
/// the coordinate axis least aligned with it; since the azimuth is drawn
/// uniformly the choice of frame does not matter, it only has to be fixed so
/// results are reproducible.
pub fn rotate(u: [f64; 3], polar: f64, azimuth: f64) -> [f64; 3] {
    rotate_sc(u, polar.sin_cos(), azimuth.sin_cos())
}

/// [`rotate`] with the polar and azimuthal angles given as `(sin, cos)`
/// pairs, so a caller that already has them (the collision code gets the
/// lab deflection from the half-angle tangent, without an inverse
/// trigonometric call) pays for no `sin_cos`.
///
/// `u` must be a unit vector to rounding and `polar` must have
/// `sin^2 + cos^2 = 1` to rounding. The result is renormalised to first
/// order, `v (3 - |v|^2) / 2`, which is exact to second order in
/// `|v|^2 - 1` (about 1e-16 here) and so costs a multiply instead of a square
/// root and a division; the rounding drift of a long chain of rotations
/// therefore does not accumulate.
#[inline]
pub fn rotate_sc(u: [f64; 3], polar: (f64, f64), azimuth: (f64, f64)) -> [f64; 3] {
    let k = if u[0].abs() < 0.5 {
        [1.0, 0.0, 0.0]
    } else if u[1].abs() < 0.5 {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let e1 = normalize(cross(u, k));
    let e2 = cross(u, e1);
    let (sp, cp) = polar;
    let (sa, ca) = azimuth;
    let v = [
        cp * u[0] + sp * (ca * e1[0] + sa * e2[0]),
        cp * u[1] + sp * (ca * e1[1] + sa * e2[1]),
        cp * u[2] + sp * (ca * e1[2] + sa * e2[2]),
    ];
    let f = 0.5 * (3.0 - (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
    [v[0] * f, v[1] * f, v[2] * f]
}

/// Lab deflection of the projectile, as `(sin psi, cos psi)`, from the
/// half-angle sine and cosine of the centre-of-mass angle (`s = sin(theta/2)`,
/// `c = cos(theta/2)`, `s^2 + c^2 = 1`) and `mu = M1 / M2`.
///
/// The same relation as in [`lab_angles`], `tan(psi) = sin(theta) /
/// (cos(theta) + mu)` (Goldstein, Sec. 3.11), written without the angle:
/// `sin theta = 2 s c`, `cos theta = 1 - 2 s^2`, and the vector
/// `(sin theta, cos theta + mu)` has length `sqrt(1 + 2 mu cos theta + mu^2)`.
/// The recoil deflection `phi = (pi - theta) / 2` has `(sin, cos) = (c, s)`.
#[inline]
pub fn lab_projectile_sc(s: f64, c: f64, mu: f64) -> (f64, f64) {
    let sin_t = 2.0 * s * c;
    let cos_t = 1.0 - 2.0 * s * s;
    let x = cos_t + mu;
    let inv = 1.0 / (sin_t * sin_t + x * x).sqrt();
    (sin_t * inv, x * inv)
}

/// Planar surface barrier at a face with outward normal along `x` (sign of
/// `dir[0]` gives which face).
///
/// A particle with energy `e` reaching the face with direction cosine
/// `c = |dir[0]|` to the normal escapes if `e c^2 > e_s` (the planar barrier
/// model: only the normal component of the momentum is reduced, by the
/// surface binding energy `E_s`). P. Sigmund, Phys. Rev. 184 (1969) 383;
/// W. Eckstein, *Computer Simulation of Ion-Solid Interactions* (Springer,
/// 1991). Outside, the energy is `e - e_s`, the parallel momentum is unchanged
/// and the normal component satisfies `e' c'^2 = e c^2 - e_s` (refraction away
/// from the normal).
///
/// Returns `Some((energy outside, direction outside))`, or `None` if the
/// particle cannot escape (the caller reflects it specularly).
pub fn refract_out(e: f64, dir: [f64; 3], e_s: f64) -> Option<(f64, [f64; 3])> {
    let c = dir[0].abs();
    let e_normal = e * c * c;
    if e_normal <= e_s {
        return None;
    }
    if e_s == 0.0 {
        return Some((e, dir));
    }
    let e_out = e - e_s;
    // Momentum components scale with sqrt(energy): parallel unchanged,
    // normal from the reduced normal energy.
    let f_par = (e / e_out).sqrt();
    let c_out = ((e_normal - e_s) / e_out).sqrt().copysign(dir[0]);
    Some((e_out, normalize([c_out, dir[1] * f_par, dir[2] * f_par])))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    #[test]
    fn momentum_and_energy_are_conserved() {
        // Projectile along +x; check p1 = p1' + p2 and E = E1' + T.
        for &(m1, m2) in &[(1.0, 28.0), (28.0, 28.0), (75.0, 12.0)] {
            for &theta in &[0.01f64, 0.5, 1.5, 2.5, 3.1] {
                let e = 1000.0;
                let mu = m1 / m2;
                let t = gamma(m1, m2) * e * (0.5 * theta).sin().powi(2);
                let (psi, phi) = lab_angles(theta, mu);
                let p1 = (2.0 * m1 * e).sqrt();
                let p1p = (2.0 * m1 * (e - t)).sqrt();
                let p2 = (2.0 * m2 * t).sqrt();
                // Projectile at +psi, recoil at -phi in the scattering plane.
                let px = p1p * psi.cos() + p2 * phi.cos();
                let py = p1p * psi.sin() - p2 * phi.sin();
                assert!((px - p1).abs() < 1e-9 * p1, "{m1} {m2} {theta}");
                assert!(py.abs() < 1e-9 * p1, "{m1} {m2} {theta}");
            }
        }
    }

    #[test]
    fn rotation_gives_requested_polar_angle() {
        let dirs = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            normalize([0.3, -0.4, 0.86]),
            normalize([-0.7, 0.7, 0.1]),
        ];
        for u in dirs {
            for &pol in &[0.0, 0.2, 1.0, 2.9] {
                for &az in &[0.0, 1.0, 4.0] {
                    let v = rotate(u, pol, az);
                    assert!((dot(v, v) - 1.0).abs() < 1e-14);
                    assert!((dot(u, v) - pol.cos()).abs() < 1e-12);
                }
            }
        }
        // Opposite azimuths lie on opposite sides of u.
        let u = normalize([0.6, 0.0, 0.8]);
        let a = rotate(u, 0.3, 0.7);
        let b = rotate(u, 0.3, 0.7 + PI);
        let s = normalize([a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
        assert!((dot(s, u) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sin_cos_forms_agree_with_the_angle_forms() {
        for &mu in &[0.05, 0.4, 1.0, 2.5, 6.0] {
            for &theta in &[1e-6f64, 0.01, 0.5, 1.5, 2.5, 3.1, PI - 1e-5] {
                let (psi, phi) = lab_angles(theta, mu);
                let (s, c) = (0.5 * theta).sin_cos();
                let (sp, cp) = lab_projectile_sc(s, c, mu);
                assert!((sp - psi.sin()).abs() < 1e-12, "{mu} {theta}");
                assert!((cp - psi.cos()).abs() < 1e-12, "{mu} {theta}");
                // Recoil: (sin phi, cos phi) = (c, s).
                assert!((c - phi.sin()).abs() < 1e-14);
                assert!((s - phi.cos()).abs() < 1e-14);
            }
        }
        let u = normalize([0.3, -0.4, 0.86]);
        let a = rotate(u, 1.1, 2.2);
        let b = rotate_sc(u, 1.1f64.sin_cos(), 2.2f64.sin_cos());
        assert_eq!(a, b);
    }

    #[test]
    fn long_rotation_chain_stays_unit() {
        let mut u = normalize([0.3, -0.4, 0.86]);
        for i in 0..200_000 {
            let x = i as f64;
            u = rotate(u, 0.1 + (x * 0.37).sin().abs(), x * 1.3);
        }
        let n = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        assert!((n - 1.0).abs() < 1e-14, "{n}");
    }

    #[test]
    fn surface_refraction() {
        // Normal incidence on the front face, E = 10, Es = 4.
        let (e, d) = refract_out(10.0, [-1.0, 0.0, 0.0], 4.0).unwrap();
        assert_eq!(e, 6.0);
        assert_eq!(d, [-1.0, 0.0, 0.0]);
        // Below the barrier in normal energy: trapped.
        let c = 0.5f64; // E cos^2 = 2.5 < 4
        let dir = [-c, (1.0 - c * c).sqrt(), 0.0];
        assert!(refract_out(10.0, dir, 4.0).is_none());
        // Oblique: parallel momentum conserved, normal energy reduced.
        let c = 0.8f64;
        let dir = [c, 0.0, (1.0 - c * c).sqrt()];
        let (e, d) = refract_out(10.0, dir, 2.0).unwrap();
        assert!((e - 8.0).abs() < 1e-14);
        assert!((e * d[0] * d[0] - (10.0 * c * c - 2.0)).abs() < 1e-12);
        assert!((e * d[2] * d[2] - 10.0 * (1.0 - c * c)).abs() < 1e-12);
        assert!(d[0] > 0.0);
    }
}
