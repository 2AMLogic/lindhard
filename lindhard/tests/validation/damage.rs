//! Damage-model limits (issue #6): NRT and Kinchin-Pease, the Lindhard
//! partition, the engine's full-cascade displacement count against the
//! Kinchin-Pease estimate in the nuclear-only limit, and the cascade's
//! electronic share against the partition integral equation
//! (`partition.rs`, issue #64).

use std::f64::consts::PI;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, Beam};
use lindhard::ion::damage::{
    damage_energy_ev, kinchin_pease, lindhard_k, lindhard_reduced_energy, nrt_displacements,
    NRT_EFFICIENCY,
};
use lindhard::ion::scattering::ScatteringTable;
use lindhard::ion::stopping::none::NoStopping;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{Binning, IonTally, IonTallyConfig};

use crate::report::{num, pct, sci, Check};

pub fn checks(table: &ScatteringTable, quick: bool) -> Vec<Check> {
    let mut out = Vec::new();
    let e_d = 40.0;

    // NRT steps and its high-energy ratio to Kinchin-Pease (Norgett,
    // Robinson, Torrens, Nucl. Eng. Des. 33 (1975) 50; Kinchin and Pease,
    // Rep. Prog. Phys. 18 (1955) 1).
    let x = 2.0 * e_d / NRT_EFFICIENCY;
    let steps = nrt_displacements(e_d - 1e-9, e_d) == 0.0
        && nrt_displacements(e_d, e_d) == 1.0
        && nrt_displacements(x * (1.0 - 1e-12), e_d) == 1.0
        && (nrt_displacements(x, e_d) - 1.0).abs() < 1e-12;
    out.push(Check::holds(
        "damage.nrt_steps",
        "NRT: 0 below E_d, 1 from E_d, continuous at 2 E_d / 0.8 (E_d = 40 eV)",
        steps,
        if steps { "holds" } else { "violated" },
        "Norgett, Robinson and Torrens 1975",
    ));
    let mut worst = 0.0f64;
    for t in [1e3, 1e4, 1e6] {
        worst = worst.max((nrt_displacements(t, e_d) / kinchin_pease(t, e_d) - 0.8).abs());
    }
    out.push(Check::at_most(
        "damage.nrt_over_kp",
        "NRT / Kinchin-Pease above 2.5 E_d, |ratio - 0.8|",
        worst,
        1e-12,
        sci,
        "the NRT displacement efficiency 0.8",
    ));

    // Lindhard partition: the general k and eps reduce to the NRT self-ion
    // forms (k = 0.1337 Z^(1/6) (Z/A)^(1/2), eps = T / (86.931 Z^(7/3))).
    let (mut wk, mut we) = (0.0f64, 0.0f64);
    for (z, a) in [
        (6.0, 12.011),
        (14.0, 28.0855),
        (26.0, 55.845),
        (74.0, 183.84),
    ] {
        let k_nrt = 0.1337 * f64::powf(z, 1.0 / 6.0) * (z / a).sqrt();
        wk = wk.max((lindhard_k(z, a, z, a) / k_nrt - 1.0).abs());
        let eps_nrt = 1e4 / (86.931 * f64::powf(z, 7.0 / 3.0));
        we = we.max((lindhard_reduced_energy(1e4, z, a, z, a) / eps_nrt - 1.0).abs());
    }
    out.push(Check::at_most(
        "damage.self_ion_k",
        "Lindhard k (general form) vs NRT self-ion form, C/Si/Fe/W",
        wk,
        3e-3,
        pct,
        "rounding of 0.1337 and 0.0793 in the published forms",
    ));
    out.push(Check::at_most(
        "damage.self_ion_eps",
        "Lindhard reduced energy (general form) vs NRT self-ion form T/(86.931 Z^(7/3)), C/Si/Fe/W",
        we,
        2e-4,
        pct,
        "rounding of 86.931 and of the constants",
    ));

    // High-energy limit of the partition: eps_dam = eps / (1 + k g(eps))
    // tends to 1/k (g(eps) ~ eps), with a leading correction
    // 0.40244 eps^(-1/4) (Robinson's g). At eps = 1e10 that is 0.13 %.
    let (z, a) = (14.0, 28.0855);
    let k = lindhard_k(z, a, z, a);
    let eps1 = 1e10;
    let t = eps1 / lindhard_reduced_energy(1.0, z, a, z, a);
    let eps_dam = lindhard_reduced_energy(damage_energy_ev(t, z, a, z, a), z, a, z, a);
    out.push(Check::at_most(
        "damage.partition_saturates",
        "Lindhard partition, Si self-ion at eps = 1e10: |k eps_dam - 1|",
        (k * eps_dam - 1.0).abs(),
        3e-3,
        pct,
        "analytic limit of T/(1 + k g(eps)); expected 0.13 % from the eps^(3/4) term",
    ));

    // Engine: full-cascade defects vs Kinchin-Pease in the nuclear-only
    // limit (no electronic loss, E_b = 0). Reported: KP counts Frenkel pairs
    // for hard spheres; the BCA efficiency for screened potentials is lower,
    // which is where NRT's 0.8 came from (Robinson and Torrens, Phys. Rev. B
    // 9 (1974) 5008), but no published number applies to exactly this setup
    // (replacement rule, cutoffs, surface). Both cutoffs are set to E_d, the
    // KP assumption that an atom below E_d displaces nothing more and stays
    // where it is: then a displacer left below E_d stops at the site it
    // emptied and counts as a replacement. With a 1 eV recoil cutoff the
    // displacer moves on instead, and vacancies / KP rises to about 1.26
    // (measured at review time), a property of the replacement rule
    // documented on `CascadeDefects`.
    let count = if quick { 300 } else { 3000 };
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m.set_lattice_binding_energy_ev(14, 0.0).unwrap();
    let stack = Stack::semi_infinite(m);
    let mut cfg = BcaConfig::new(15.0, 15.0);
    cfg.seed = 0xD15;
    let bca = Bca::new(
        Beam::normal(Ion::new(14).unwrap(), 2.0e3, count),
        &stack,
        cfg,
        &NoStopping,
        table,
    )
    .unwrap();
    let tc = IonTallyConfig {
        depth: Binning::new(0.0, 20e-9, 20).unwrap(),
        lateral: Binning::new(-10e-9, 10e-9, 10).unwrap(),
        radial: Binning::new(0.0, 10e-9, 10).unwrap(),
        escape_energy: Binning::new(0.0, 2000.0, 10).unwrap(),
        escape_polar: Binning::new(0.0, 0.5 * PI, 9).unwrap(),
    };
    let proto = IonTally::new(&stack, &bca.species_z(), tc).unwrap();
    let r = bca.run(|| proto.clone()).unwrap().report(false);
    let e_avail = r.budget.lattice + r.budget.rest;
    let kp = e_avail / (2.0 * 15.0);
    let c = r.damage.cascade;
    out.push(Check::info(
        "damage.cascade_over_kp",
        format!(
            "Si 2 keV -> Si, no electronic loss, E_d = 15 eV, E_b = 0, cutoffs = E_d: engine vacancies / KP E/(2 E_d), E = energy kept in the target, {count} ions"
        ),
        format!(
            "{} (displacements {})",
            num(c.vacancies as f64 / kp),
            num(c.displacements as f64 / kp)
        ),
        "context: Robinson and Torrens (1974) found BCA efficiencies near 0.8 (NRT's factor); vacancies = displacements - replacements, so the ratio depends on the replacement rule and cutoffs (about 1.26 with a 1 eV recoil cutoff)",
    ));

    // The cascade's electronic share against the Lindhard partition (#64).
    out.extend(crate::partition::checks(table, quick));
    out
}
