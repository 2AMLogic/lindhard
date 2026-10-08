//! Wafer cut and beam angles: the lab-to-crystal rotation from
//! (tilt, twist, wafer rotation).
//!
//! # Conventions
//!
//! The figure and a worked example are in `docs/crystal-orientation.md`.
//!
//! * **Lab frame** (as in [`crate::geometry`]): `+x` is depth, pointing into
//!   the target along the surface normal; `y`, `z` are lateral.
//! * **Wafer cut.** The surface is the plane `(hkl)` (Miller indices on the
//!   conventional cubic cell). Its normal `n = g_hkl / |g_hkl|`
//!   ([`super::Lattice::plane_normal`]) is the crystal direction of lab
//!   `+x`, i.e. `(hkl)` names the **inward** normal: at zero tilt the beam
//!   travels along `n`. (To name the outward normal `(hkl)`, pass
//!   `(-h, -k, -l)`; it matters only for polar faces.) The **reference
//!   direction** `[uvw]` is a lattice direction in the surface plane (the
//!   zone law `h u + k v + l w = 0` must hold); its unit vector `r` is the
//!   crystal direction of lab `+y` at zero wafer rotation. The third axis is
//!   `t = n x r`, the crystal direction of lab `+z`, so `(n, r, t)` is
//!   right-handed like `(x, y, z)`.
//! * **Tilt** `θ` is the polar angle of the beam from lab `+x` (the inward
//!   normal), in `[0, π/2)`.
//! * **Twist** `φ` is the azimuth of the beam about lab `+x`, measured in the
//!   surface plane from lab `+y` towards lab `+z`. The beam direction in the
//!   lab is `(cos θ, sin θ cos φ, sin θ sin φ)`, the same as
//!   [`crate::ion::bca::Beam::direction`] with `polar_rad = θ` and
//!   `azimuth_rad = φ`.
//! * **Wafer rotation** `ω` turns the wafer (and the crystal) about lab `+x`
//!   by `ω`, right-handed: the reference direction `r` moves from lab `+y`
//!   towards lab `+z`.
//!
//! So the crystal-to-lab rotation is `M = R_x(ω) B`, where `B` has rows
//! `n`, `r`, `t` (it takes crystal components to wafer components) and
//! `R_x(ω)` is the right-handed rotation by `ω` about `x`. The lab-to-crystal
//! rotation is its transpose, and the beam in the crystal frame is
//!
//! ```text
//! d = cos θ n + sin θ cos(φ - ω) r + sin θ sin(φ - ω) t.
//! ```
//!
//! For the beam direction only the difference `φ - ω` matters; the two are
//! kept apart because the lab frame (where the beam and the target geometry
//! live) and the crystal frame are distinct.
//!
//! **Reference direction.** The meaning of "twist" differs between tools,
//! because the azimuth needs a zero. Here the zero is always the reference
//! direction the caller names, so there is no hidden default. As a cited
//! example of the polar/azimuth decomposition of an incidence direction on a
//! cubic crystal: K. Nordlund, F. Djurabekova and G. Hobler, "Large fraction
//! of crystal directions leads to ion channeling", Phys. Rev. B 94, 214109
//! (2016), doi:10.1103/PhysRevB.94.214109, Sec. II B, p. 214109-3, set up a
//! cell "with a [001] surface normal and tilting (θ) and twisting (ϕ) the
//! incoming ion direction", and reach `[011]` at `θ = 45°, ϕ = 0°` and
//! `[111]` at `θ = 54.73°, ϕ = 45°`: their twist is zero towards an in-plane
//! `<100>` axis. With `n = [001]` and `r = [010]` our convention gives the
//! same two channels (the second as `[-111]`, a member of `<111>`; the paper
//! does not fix the sense of `ϕ`, which we fix above). A wafer flat or notch
//! direction can be passed as the reference instead.

use super::{add, cross, dot, scale, CrystalError, Lattice};

/// A 3x3 rotation matrix, row-major: `v' = M v` with
/// `v'_i = sum_j M[i][j] v_j`.
pub type Mat3 = [[f64; 3]; 3];

fn mat_vec(m: &Mat3, v: [f64; 3]) -> [f64; 3] {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

fn transpose(m: &Mat3) -> Mat3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

fn check_angle(
    name: &'static str,
    value: f64,
    ok: bool,
    why: &'static str,
) -> Result<f64, CrystalError> {
    if value.is_finite() && ok {
        Ok(value)
    } else {
        Err(CrystalError::InvalidAngle { name, value, why })
    }
}

/// Unit beam direction in the lab frame for tilt `θ` and twist `φ`:
/// `(cos θ, sin θ cos φ, sin θ sin φ)`.
pub fn lab_direction(tilt_rad: f64, twist_rad: f64) -> [f64; 3] {
    let (s, c) = tilt_rad.sin_cos();
    let (sa, ca) = twist_rad.sin_cos();
    [c, s * ca, s * sa]
}

/// A wafer cut and the beam angles. See the module docs for the
/// conventions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orientation {
    normal_hkl: [i32; 3],
    reference_uvw: [i32; 3],
    /// Rows `n`, `r`, `t`: crystal components to wafer components.
    wafer: Mat3,
    tilt_rad: f64,
    twist_rad: f64,
    wafer_rotation_rad: f64,
}

impl Orientation {
    /// Orientation of a wafer of `lattice` cut on `(hkl)` (inward normal)
    /// with in-plane reference direction `[uvw]`, and a beam at `tilt_rad`
    /// (in `[0, π/2)`), `twist_rad` and `wafer_rotation_rad` (finite).
    pub fn new(
        lattice: &Lattice,
        normal_hkl: [i32; 3],
        reference_uvw: [i32; 3],
        tilt_rad: f64,
        twist_rad: f64,
        wafer_rotation_rad: f64,
    ) -> Result<Self, CrystalError> {
        let n = lattice.plane_normal(normal_hkl)?;
        let r = lattice.direction(reference_uvw)?;
        let zone: i64 = (0..3)
            .map(|i| i64::from(normal_hkl[i]) * i64::from(reference_uvw[i]))
            .sum();
        if zone != 0 {
            return Err(CrystalError::ReferenceNotInPlane {
                normal: normal_hkl,
                reference: reference_uvw,
                dot: zone,
            });
        }
        check_angle(
            "tilt",
            tilt_rad,
            (0.0..std::f64::consts::FRAC_PI_2).contains(&tilt_rad),
            "must be in [0, pi/2)",
        )?;
        check_angle("twist", twist_rad, true, "must be finite")?;
        check_angle("wafer rotation", wafer_rotation_rad, true, "must be finite")?;
        let t = cross(n, r);
        Ok(Self {
            normal_hkl,
            reference_uvw,
            wafer: [n, r, t],
            tilt_rad,
            twist_rad,
            wafer_rotation_rad,
        })
    }

    /// Miller indices of the surface plane (inward normal).
    pub fn normal_hkl(&self) -> [i32; 3] {
        self.normal_hkl
    }

    /// Indices of the in-plane reference direction.
    pub fn reference_uvw(&self) -> [i32; 3] {
        self.reference_uvw
    }

    /// Tilt, radians.
    pub fn tilt_rad(&self) -> f64 {
        self.tilt_rad
    }

    /// Twist, radians.
    pub fn twist_rad(&self) -> f64 {
        self.twist_rad
    }

    /// Wafer rotation, radians.
    pub fn wafer_rotation_rad(&self) -> f64 {
        self.wafer_rotation_rad
    }

    /// Crystal-to-lab rotation `M = R_x(ω) B`.
    pub fn crystal_to_lab(&self) -> Mat3 {
        let (s, c) = self.wafer_rotation_rad.sin_cos();
        let b = &self.wafer;
        // R_x(ω) = [[1, 0, 0], [0, c, -s], [0, s, c]] times B.
        [
            b[0],
            add(scale(b[1], c), scale(b[2], -s)),
            add(scale(b[1], s), scale(b[2], c)),
        ]
    }

    /// Lab-to-crystal rotation, the transpose of [`Self::crystal_to_lab`].
    pub fn lab_to_crystal(&self) -> Mat3 {
        transpose(&self.crystal_to_lab())
    }

    /// Rotate a lab-frame vector into the crystal frame.
    pub fn to_crystal(&self, v_lab: [f64; 3]) -> [f64; 3] {
        mat_vec(&self.lab_to_crystal(), v_lab)
    }

    /// Rotate a crystal-frame vector into the lab frame.
    pub fn to_lab(&self, v_crystal: [f64; 3]) -> [f64; 3] {
        mat_vec(&self.crystal_to_lab(), v_crystal)
    }

    /// Unit beam direction in the lab frame, [`lab_direction`] of the tilt
    /// and twist.
    pub fn beam_lab(&self) -> [f64; 3] {
        lab_direction(self.tilt_rad, self.twist_rad)
    }

    /// Unit beam direction in the crystal frame.
    pub fn beam_crystal(&self) -> [f64; 3] {
        self.to_crystal(self.beam_lab())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
        (0..3).all(|i| (a[i] - b[i]).abs() <= tol)
    }

    #[test]
    fn rotation_is_orthonormal_and_right_handed() {
        let si = Lattice::silicon();
        let o = Orientation::new(&si, [1, 1, 1], [1, -1, 0], 0.3, 1.1, -0.4).unwrap();
        let m = o.crystal_to_lab();
        for i in 0..3 {
            for j in 0..3 {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((dot(m[i], m[j]) - want).abs() < 1e-15);
            }
        }
        let det = dot(m[0], cross(m[1], m[2]));
        assert!((det - 1.0).abs() < 1e-15);
        let v = [0.2, -0.7, 0.4];
        assert!(close(o.to_lab(o.to_crystal(v)), v, 1e-15));
    }

    #[test]
    fn surface_normal_maps_to_lab_x() {
        let si = Lattice::silicon();
        let o = Orientation::new(&si, [1, 1, 0], [0, 0, 1], 0.2, 0.5, 0.9).unwrap();
        let n = si.plane_normal([1, 1, 0]).unwrap();
        assert!(close(o.to_lab(n), [1.0, 0.0, 0.0], 1e-15));
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        let si = Lattice::silicon();
        assert_eq!(
            Orientation::new(&si, [1, 0, 0], [1, 1, 0], 0.0, 0.0, 0.0),
            Err(CrystalError::ReferenceNotInPlane {
                normal: [1, 0, 0],
                reference: [1, 1, 0],
                dot: 1
            })
        );
        for tilt in [-0.1, std::f64::consts::FRAC_PI_2, f64::NAN] {
            assert!(matches!(
                Orientation::new(&si, [1, 0, 0], [0, 1, 0], tilt, 0.0, 0.0),
                Err(CrystalError::InvalidAngle { name: "tilt", .. })
            ));
        }
        assert!(matches!(
            Orientation::new(&si, [1, 0, 0], [0, 1, 0], 0.1, f64::INFINITY, 0.0),
            Err(CrystalError::InvalidAngle { name: "twist", .. })
        ));
        assert!(matches!(
            Orientation::new(&si, [1, 0, 0], [0, 1, 0], 0.1, 0.0, f64::NAN),
            Err(CrystalError::InvalidAngle {
                name: "wafer rotation",
                ..
            })
        ));
        assert_eq!(
            Orientation::new(&si, [0, 0, 0], [0, 1, 0], 0.0, 0.0, 0.0),
            Err(CrystalError::ZeroIndex([0; 3]))
        );
    }
}
