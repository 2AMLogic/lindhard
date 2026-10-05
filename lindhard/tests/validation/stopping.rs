//! Electronic-stopping identities (issue #4). These check the formulas
//! against themselves and each other, not against measurement: comparing the
//! stopping magnitudes with experiment is level 3 (`docs/validation.md`).

use std::f64::consts::PI;

use lindhard::constants::{BOHR_RADIUS, COULOMB_E2};
use lindhard::ion::stopping::bragg::{bragg_cross_section_per_atom, NoCorrection};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::oen_robinson::{local_length, OenRobinson};
use lindhard::ion::stopping::{bohr_velocity, speed_nonrel, ElectronicStopping, Ion};
use lindhard::material::Material;

use crate::report::{pct, sci, Check};

pub fn checks() -> Vec<Check> {
    let ls = LindhardScharff::new();
    let mut out = Vec::new();

    // LS reduced form (k_L eps^(1/2), with xi_e = Z1^(1/6)) vs the
    // dimensional form S_e = Z1^(1/6) 8 pi e^2 a0 Z1 Z2 / (Z1^(2/3) +
    // Z2^(2/3))^(3/2) v/v0 (Lindhard and Scharff, Phys. Rev. 124 (1961) 128).
    let mut worst = 0.0f64;
    for (z1, z2) in [(1u8, 14u8), (5, 14), (15, 14), (33, 14), (11, 29), (5, 79)] {
        let ion = Ion::new(z1).unwrap();
        let e = 2.0e3;
        let s = ls.stopping(&ion, z2, e).unwrap();
        let (a, b) = (f64::from(z1), f64::from(z2));
        let v = speed_nonrel(&ion, e) / bohr_velocity();
        let dim = a.powf(1.0 / 6.0) * 8.0 * PI * COULOMB_E2 * BOHR_RADIUS * a * b
            / (a.powf(2.0 / 3.0) + b.powf(2.0 / 3.0)).powf(1.5)
            * v;
        worst = worst.max((s / dim - 1.0).abs());
    }
    out.push(Check::at_most(
        "stopping.ls_reduced_vs_dimensional",
        "Lindhard-Scharff: reduced form k_L eps^(1/2) vs dimensional form, H/B/P/As in Si, Na in Cu, B in Au",
        worst,
        1e-2,
        pct,
        "the rounded 0.0793 in k_L; LS 1961, LSS 1963",
    ));

    // Velocity proportionality: S_e / E^(1/2) constant over three decades.
    let ion = Ion::new(15).unwrap();
    let r0 = ls.stopping(&ion, 14, 1e2).unwrap() / 1e2f64.sqrt();
    let mut worst = 0.0f64;
    for e in [1e3, 1e4, 1e5] {
        worst = worst.max((ls.stopping(&ion, 14, e).unwrap() / e.sqrt() / r0 - 1.0).abs());
    }
    out.push(Check::at_most(
        "stopping.ls_velocity_proportional",
        "Lindhard-Scharff S_e / E^(1/2), P in Si, 100 eV to 100 keV, |rel. spread|",
        worst,
        1e-12,
        sci,
        "S_e proportional to v below v0 Z1^(2/3)",
    ));

    // Oen-Robinson: the impact-parameter average of the local loss is the LS
    // stopping (the normalisation of Oen and Robinson, NIM 132 (1976) 647).
    let or = OenRobinson::new();
    let s = or.stopping(&ion, 14, 1.0e4).unwrap();
    let a = local_length(15, 14);
    let (n, pmax) = (200_000u32, 60.0 * a);
    let dp = pmax / f64::from(n);
    let mut sum = 0.0;
    for i in 0..n {
        let p = (f64::from(i) + 0.5) * dp;
        sum += or.local_loss(&ion, 14, 1.0e4, p).unwrap() * 2.0 * PI * p * dp;
    }
    out.push(Check::at_most(
        "stopping.oen_robinson_normalisation",
        "Oen-Robinson local loss integrated over 2 pi p dp vs the LS stopping, P in Si at 10 keV",
        (sum / s - 1.0).abs(),
        1e-4,
        sci,
        "midpoint rule to 60 a; exp(-18) tail",
    ));

    // Bragg additivity: SiO2 per-atom cross section is the atom-fraction
    // weighted sum (W. H. Bragg and R. Kleeman, Phil. Mag. 10, 318 (1905)).
    let ion = Ion::new(5).unwrap();
    let sio2 = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
    let e = 2.0e4;
    let hand = (ls.stopping(&ion, 14, e).unwrap() + 2.0 * ls.stopping(&ion, 8, e).unwrap()) / 3.0;
    let got = bragg_cross_section_per_atom(&ls, &NoCorrection, &ion, &sio2, e).unwrap();
    out.push(Check::at_most(
        "stopping.bragg_additivity",
        "Bragg sum for SiO2 vs hand-weighted (S_Si + 2 S_O)/3, B at 20 keV",
        (got / hand - 1.0).abs(),
        1e-14,
        sci,
        "rounding only",
    ));
    out
}
