//! Displacement-damage models: the Lindhard partition, NRT and Kinchin-Pease.
//!
//! These are *model estimates* of the number of displacements made by a
//! primary knock-on atom (PKA) of a given energy. They are kept apart from the
//! *event counts* of a full-cascade simulation (vacancies, interstitials,
//! replacements, [`crate::tally::CascadeDefects`]), which the BCA engine
//! produces by following every recoil. The two answer different questions and
//! are reported side by side, never added or substituted for each other.
//!
//! # Damage energy (Lindhard partition)
//!
//! Of a PKA's energy `T`, the part eventually given to atomic motion (the
//! damage energy) is `T_dam = T / (1 + k g(ε))`, J. Lindhard, V. Nielsen,
//! M. Scharff and P. V. Thomsen, "Integral equations governing radiation
//! effects", Mat. Fys. Medd. Dan. Vid. Selsk. 33 (10) (1963), with the
//! analytic fit `g(ε) = 3.4008 ε^(1/6) + 0.40244 ε^(3/4) + ε` of M. T.
//! Robinson (in *Nuclear Fusion Reactors*, British Nuclear Energy Society,
//! London, 1970, p. 364), as adopted by the NRT standard below. Here `ε` is
//! the Lindhard reduced energy and `k` the Lindhard-Scharff electronic
//! stopping coefficient of the PKA (atomic number `Z1`, mass `A1`) in the
//! target (`Z2`, `A2`), as in [`crate::ion::stopping::lindhard_scharff`]:
//!
//! ```text
//! ε = T a A2 / (Z1 Z2 e^2 (A1 + A2)),   a = 0.8853 a0 (Z1^(2/3) + Z2^(2/3))^(-1/2),
//! k = 0.0793 Z1^(2/3) Z2^(1/2) (A1 + A2)^(3/2) / ((Z1^(2/3) + Z2^(2/3))^(3/4) A1^(3/2) A2^(1/2)).
//! ```
//!
//! For a self-ion (`Z1 = Z2 = Z`, `A1 = A2 = A`) these reduce to the forms
//! written in the NRT standard, `ε = T / (86.931 Z^(7/3))` (T in eV) and
//! `k = 0.1337 Z^(1/6) (Z/A)^(1/2)` (checked in a test).
//!
//! # NRT and Kinchin-Pease
//!
//! M. J. Norgett, M. T. Robinson and I. M. Torrens, "A proposed method of
//! calculating displacement dose rates", Nucl. Eng. Des. 33 (1975) 50:
//!
//! ```text
//! N_NRT = 0                      for T_dam < E_d,
//!       = 1                      for E_d <= T_dam < 2 E_d / 0.8,
//!       = 0.8 T_dam / (2 E_d)    for T_dam >= 2 E_d / 0.8,
//! ```
//!
//! with `E_d` the displacement threshold and 0.8 the displacement efficiency
//! the authors took from BCA simulations. The older G. H. Kinchin and R. S.
//! Pease, "The displacement of atoms in solids by radiation", Rep. Prog.
//! Phys. 18 (1955) 1 count is `0`, `1`, or `T / (2 E_d)` with thresholds `E_d`
//! and `2 E_d`; [`kinchin_pease`] is evaluated on whatever energy it is given
//! (the tallies pass the damage energy, the usual "modified" form).
//!
//! **NRT is not a cascade count.** R. E. Stoller, M. B. Toloczko, G. S. Was,
//! A. G. Certain, S. Dwaraknath and F. A. Garner, "On the use of SRIM for
//! computing radiation damage exposure", Nucl. Instrum. Methods B 310 (2013)
//! 75, discuss how BCA full-cascade vacancy counts and NRT displacement
//! estimates differ and recommend that, when a displacement dose comparable
//! with the NRT standard is wanted, it be computed from the damage energy with
//! the NRT formula rather than taken from the full-cascade vacancy count. NRT
//! is also known to overestimate the number of defects that survive
//! in-cascade recombination (as seen in molecular dynamics); it is a standard
//! exposure unit, not a prediction of surviving defects. This crate therefore
//! reports both and lets the user choose.
//!
//! # Compound and multi-element targets (an approximation)
//!
//! NRT (Norgett, Robinson and Torrens 1975) and the Lindhard partition it
//! uses are defined for a monatomic target. For a layer with more than one
//! element, the tally ([`crate::tally::IonTally`]) evaluates the partition
//! with:
//!
//! * `Z1`, `A1`: the PKA's own atomic number and mass;
//! * `Z2`, `A2`: the atom-fraction-weighted mean atomic number and mean
//!   mass of the layer the PKA starts in (non-integer in general);
//! * `E_d`: the PKA element's own displacement threshold in that layer.
//!
//! This is an approximation outside NRT's monatomic definition. It is
//! **this crate's own convention**; no published source is claimed for this
//! particular averaging. For a monatomic layer it reduces exactly to the
//! standard NRT evaluation. Treat NRT numbers for compounds as an exposure
//! index that is comparable within this crate, not as a value that follows
//! any compound-target standard.

use crate::constants::{BOHR_RADIUS, COULOMB_E2};
use crate::units::J_PER_EV;

/// NRT displacement efficiency (Norgett, Robinson and Torrens 1975).
pub const NRT_EFFICIENCY: f64 = 0.8;

/// Lindhard reduced energy `ε` of a PKA of energy `t_ev` (`Z1`, mass `a1`
/// u) in a target (`Z2`, `a2` u). `Z` and masses may be non-integer
/// (averaged) values.
pub fn lindhard_reduced_energy(t_ev: f64, z1: f64, a1: f64, z2: f64, a2: f64) -> f64 {
    let a = 0.8853 * BOHR_RADIUS / (z1.powf(2.0 / 3.0) + z2.powf(2.0 / 3.0)).sqrt();
    t_ev * J_PER_EV * a * a2 / (z1 * z2 * COULOMB_E2 * (a1 + a2))
}

/// Lindhard-Scharff electronic stopping coefficient `k` of the PKA in the
/// target (see the module docs).
pub fn lindhard_k(z1: f64, a1: f64, z2: f64, a2: f64) -> f64 {
    0.0793 * z1.powf(2.0 / 3.0) * z2.sqrt() * (a1 + a2).powf(1.5)
        / ((z1.powf(2.0 / 3.0) + z2.powf(2.0 / 3.0)).powf(0.75) * a1.powf(1.5) * a2.sqrt())
}

/// Robinson's fit `g(ε) = 3.4008 ε^(1/6) + 0.40244 ε^(3/4) + ε`.
pub fn robinson_g(eps: f64) -> f64 {
    3.4008 * eps.powf(1.0 / 6.0) + 0.40244 * eps.powf(0.75) + eps
}

/// Damage energy `T / (1 + k g(ε))`, eV, of a PKA of energy `t_ev`.
/// Returns 0 for non-positive `t_ev`.
pub fn damage_energy_ev(t_ev: f64, z1: f64, a1: f64, z2: f64, a2: f64) -> f64 {
    if t_ev.is_nan() || t_ev <= 0.0 {
        return 0.0;
    }
    let eps = lindhard_reduced_energy(t_ev, z1, a1, z2, a2);
    t_ev / (1.0 + lindhard_k(z1, a1, z2, a2) * robinson_g(eps))
}

/// NRT displacements for damage energy `t_dam_ev` and threshold `e_d_ev`.
pub fn nrt_displacements(t_dam_ev: f64, e_d_ev: f64) -> f64 {
    if t_dam_ev < e_d_ev {
        0.0
    } else if t_dam_ev < 2.0 * e_d_ev / NRT_EFFICIENCY {
        1.0
    } else {
        NRT_EFFICIENCY * t_dam_ev / (2.0 * e_d_ev)
    }
}

/// Kinchin-Pease displacements for energy `energy_ev` and threshold
/// `e_d_ev`.
pub fn kinchin_pease(energy_ev: f64, e_d_ev: f64) -> f64 {
    if energy_ev < e_d_ev {
        0.0
    } else if energy_ev < 2.0 * e_d_ev {
        1.0
    } else {
        energy_ev / (2.0 * e_d_ev)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn self_ion_forms_match_the_nrt_constants() {
        for (z, a) in [(14.0, 28.0855), (26.0, 55.845), (74.0, 183.84)] {
            let k = lindhard_k(z, a, z, a);
            let k_nrt = 0.1337 * f64::powf(z, 1.0 / 6.0) * (z / a).sqrt();
            assert!((k / k_nrt - 1.0).abs() < 3e-3, "k: {k} vs {k_nrt}");
            let eps = lindhard_reduced_energy(1e4, z, a, z, a);
            let eps_nrt = 1e4 / (86.931 * f64::powf(z, 7.0 / 3.0));
            assert!(
                (eps / eps_nrt - 1.0).abs() < 2e-4,
                "eps: {eps} vs {eps_nrt}"
            );
        }
    }

    #[test]
    fn damage_energy_is_below_t_and_tends_to_t_at_low_energy() {
        let (z, a) = (26.0, 55.845);
        let t = 1e5;
        let td = damage_energy_ev(t, z, a, z, a);
        assert!(td > 0.3 * t && td < t);
        let low = damage_energy_ev(10.0, z, a, z, a);
        assert!(low < 10.0 && low > 9.0);
        assert_eq!(damage_energy_ev(0.0, z, a, z, a), 0.0);
    }

    #[test]
    fn displacement_count_steps() {
        let ed = 40.0;
        assert_eq!(nrt_displacements(39.9, ed), 0.0);
        assert_eq!(nrt_displacements(40.0, ed), 1.0);
        assert_eq!(nrt_displacements(99.9, ed), 1.0);
        assert_eq!(nrt_displacements(100.0, ed), 1.0);
        assert!((nrt_displacements(1000.0, ed) - 10.0).abs() < 1e-12);
        assert_eq!(kinchin_pease(79.9, ed), 1.0);
        assert!((kinchin_pease(1000.0, ed) - 12.5).abs() < 1e-12);
    }
}
