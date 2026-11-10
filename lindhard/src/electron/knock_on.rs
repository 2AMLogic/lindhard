//! Knock-on (displacement) threshold report for an electron beam.
//!
//! Pure kinematics on resolved input; it needs no transport. For a beam
//! electron of kinetic energy `E` hitting a free target atom of mass `M` at
//! rest head on, conservation of energy and momentum gives the largest energy
//! the atom can take (relativistic, electron rest energy `mc²`):
//!
//! ```text
//! T_max = 2 E (E + 2 m c²) / (M c²)
//! ```
//!
//! This follows from conservation of energy and momentum for a head-on
//! elastic collision with a heavy atom at rest (electron momentum
//! `p²c² = E (E + 2 m c²)`, recoil `T = 2 p²c² / (M c²)`); it is the form
//! used in electron-irradiation damage work, for example F. Seitz and J. S.
//! Koehler, Solid State Physics 2, 305 (1956) (not opened for this change;
//! the derivation above stands on its own). It treats the atom as free and at
//! rest and the recoil as non-relativistic (`T_max << Mc²`), so lattice
//! binding, thermal motion and the crystal direction are ignored: `E_d` is
//! the threshold displacement energy of the input, a damage-accounting
//! convention, not a measured property of one direction.
//!
//! An atom can be displaced when `T_max >= E_d`. Setting `T_max = E_d` and
//! solving the quadratic in `E` gives the threshold beam energy:
//!
//! ```text
//! E_th = -m c² + sqrt((m c²)² + E_d M c² / 2)
//! ```
//!
//! evaluated as `(E_d M c² / 2) / (m c² + sqrt(...))` so that it does not
//! cancel for small `E_d / Mc²`.
//!
//! The report says nothing about a displacement cross section or a dpa: it
//! is a threshold test only. An element with no `E_d` in the input is
//! reported with `e_d_ev = None` (printed "not set"); no default is invented
//! here.
//!
//! `M` is the atomic weight of the element (relative atomic mass) times the
//! atomic mass unit; the electron mass difference of the atom is neglected.

use serde::{Deserialize, Serialize};

use crate::constants::{ATOMIC_MASS_UNIT, ELECTRON_REST_ENERGY, ELEMENTARY_CHARGE, SPEED_OF_LIGHT};
use crate::elements::element;
use crate::material::Material;

/// Electron rest energy `m c²`, eV.
pub fn electron_rest_energy_ev() -> f64 {
    ELECTRON_REST_ENERGY / ELEMENTARY_CHARGE
}

/// Rest energy `M c²`, eV, of an atom of `atomic_weight` (relative atomic
/// mass, i.e. g/mol).
pub fn atom_rest_energy_ev(atomic_weight: f64) -> f64 {
    atomic_weight * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT * SPEED_OF_LIGHT / ELEMENTARY_CHARGE
}

/// Largest energy, eV, an electron of kinetic energy `beam_ev` can transfer
/// to a free atom at rest of rest energy `atom_rest_ev` (see the module docs).
pub fn t_max_ev(beam_ev: f64, atom_rest_ev: f64) -> f64 {
    2.0 * beam_ev * (beam_ev + 2.0 * electron_rest_energy_ev()) / atom_rest_ev
}

/// Beam energy, eV, at which `T_max` equals `e_d_ev` (see the module docs).
pub fn threshold_beam_energy_ev(e_d_ev: f64, atom_rest_ev: f64) -> f64 {
    let mc2 = electron_rest_energy_ev();
    let c = e_d_ev * atom_rest_ev / 2.0;
    c / (mc2 + (mc2 * mc2 + c).sqrt())
}

/// The knock-on result for one element of a layer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KnockOnElement {
    /// Atomic number.
    pub z: u8,
    /// Element symbol.
    pub symbol: String,
    /// Atomic weight used for `M`, g/mol.
    pub atomic_weight: f64,
    /// `T_max` at the beam energy, eV.
    pub t_max_ev: f64,
    /// The element's displacement energy `E_d`, eV; `None` if not set.
    pub e_d_ev: Option<f64>,
    /// Whether `T_max >= E_d`; `None` if `E_d` is not set.
    pub can_displace: Option<bool>,
    /// Beam energy at which `T_max = E_d`, eV; `None` if `E_d` is not set.
    pub threshold_beam_ev: Option<f64>,
}

/// The knock-on result for every element of one layer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KnockOnLayer {
    /// Layer index in stack order.
    pub layer: usize,
    /// The beam energy the test was made at, eV.
    pub beam_ev: f64,
    /// One entry per element, in the material's order.
    pub elements: Vec<KnockOnElement>,
}

/// Knock-on report of layer `layer`, made of `material`, for a beam of
/// `beam_ev`. `E_d` is read from the material; an element without one is
/// reported as not set.
pub fn layer_report(layer: usize, material: &Material, beam_ev: f64) -> KnockOnLayer {
    let elements = material
        .components()
        .iter()
        .map(|c| {
            let z = c.z();
            let el = element(z).expect("material elements are validated");
            let rest = atom_rest_energy_ev(el.atomic_weight);
            let t = t_max_ev(beam_ev, rest);
            let e_d = material.displacement_energy_ev(z).ok();
            KnockOnElement {
                z,
                symbol: el.symbol.to_string(),
                atomic_weight: el.atomic_weight,
                t_max_ev: t,
                e_d_ev: e_d,
                can_displace: e_d.map(|d| t >= d),
                threshold_beam_ev: e_d.map(|d| threshold_beam_energy_ev(d, rest)),
            }
        })
        .collect();
    KnockOnLayer {
        layer,
        beam_ev,
        elements,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Closed form at the threshold, written out independently: with
    /// `x = E_th`, `2 x (x + 2 mc²) = E_d Mc²`.
    fn check(z: u8, e_d: f64) {
        let el = element(z).unwrap();
        let m = atom_rest_energy_ev(el.atomic_weight);
        let mc2 = 510_998.950_0; // electron rest energy, eV (CODATA 2022)
        let th = threshold_beam_energy_ev(e_d, m);
        let direct = -mc2 + (mc2 * mc2 + e_d * m / 2.0).sqrt();
        assert!(
            (th - direct).abs() < 1e-6 * direct,
            "Z={z}: {th} vs {direct}"
        );
        // T_max at the threshold is E_d, and straddles it either side.
        assert!((t_max_ev(th, m) - e_d).abs() < 1e-9 * e_d);
        assert!(t_max_ev(th * 0.99, m) < e_d);
        assert!(t_max_ev(th * 1.01, m) > e_d);
    }

    #[test]
    fn threshold_matches_closed_form_for_three_elements() {
        check(6, 35.0); // C
        check(14, 21.0); // Si
        check(29, 25.0); // Cu
    }

    #[test]
    fn layer_report_flags_and_unset_ed() {
        let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        m.set_displacement_energy_ev(14, 21.0).unwrap();
        let el = element(14).unwrap();
        let th = threshold_beam_energy_ev(21.0, atom_rest_energy_ev(el.atomic_weight));
        let above = layer_report(0, &m, th * 1.1);
        assert_eq!(above.elements[0].can_displace, Some(true));
        let below = layer_report(0, &m, th * 0.9);
        assert_eq!(below.elements[0].can_displace, Some(false));
        assert!(below.elements[0].t_max_ev < 21.0);
        assert!((below.elements[0].threshold_beam_ev.unwrap() - th).abs() < 1e-9 * th);
    }

    #[test]
    fn element_without_ed_is_not_set_not_an_error() {
        // Silicon has no E_d convention value (src/elements.rs), so it is unset.
        let m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        let r = layer_report(2, &m, 1.0e6);
        assert_eq!(r.layer, 2);
        let e = &r.elements[0];
        assert_eq!(
            (e.e_d_ev, e.can_displace, e.threshold_beam_ev),
            (None, None, None)
        );
        assert!(e.t_max_ev > 0.0);
    }
}
