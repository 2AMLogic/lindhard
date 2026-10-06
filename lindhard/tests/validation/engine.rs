//! Engine invariants (issue #5): energy conservation per history and
//! determinism across thread counts, plus the sputter-yield sensitivity to
//! `E_d` noted in the PR #46 review (reported, not asserted).

use std::f64::consts::PI;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, Beam, ElectronicLoss, MeanFreePath, SummaryTally};
use lindhard::ion::scattering::ScatteringTable;
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{Binning, IonTally, IonTallyConfig};

use crate::report::{num, sci, Check};

const NM: f64 = 1e-9;

/// `E_d` = 15 eV everywhere and `E_s` = 2 eV where the element table has
/// none (O). Test parameters, not data.
fn ready(mut m: Material) -> Material {
    for (z, _) in m.atom_fractions() {
        m.set_displacement_energy_ev(z, 15.0).unwrap();
        if m.surface_binding_energy_ev(z).is_err() {
            m.set_surface_binding_energy_ev(z, 2.0).unwrap();
        }
    }
    m
}

fn si() -> Material {
    ready(Material::from_atom_fractions(&[(14, 1.0)], None).unwrap())
}

fn sio2() -> Material {
    ready(Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap())
}

fn with_threads<R: Send>(n: usize, f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .unwrap()
        .install(f)
}

pub fn checks(table: &ScatteringTable, quick: bool) -> Vec<Check> {
    let ls = LindhardScharff::new();
    let mut out = Vec::new();

    // Energy conservation: every code path (two materials, front and back
    // faces, cascades, both free paths, both electronic modes).
    let stack = Stack::new(vec![(sio2(), 3.0 * NM), (si(), 4.0 * NM)], None).unwrap();
    let mut worst = 0.0f64;
    let mut histories = 0;
    for (mfp, el, weak) in [
        (MeanFreePath::Constant, ElectronicLoss::NonLocal, 0),
        (
            MeanFreePath::EnergyDependent {
                min_cm_angle_rad: 0.01,
            },
            ElectronicLoss::EquipartitionLsOr,
            0,
        ),
        (MeanFreePath::Constant, ElectronicLoss::EquipartitionLsOr, 3),
    ] {
        let mut cfg = BcaConfig::new(5.0, 1.0);
        cfg.mean_free_path = mfp;
        cfg.electronic = el;
        cfg.weak_collisions = weak;
        cfg.seed = 7;
        let beam = Beam {
            ion: Ion::new(18).unwrap(),
            energy_ev: 3.0e3,
            polar_rad: 0.6,
            azimuth_rad: 0.3,
            count: 300,
        };
        let t = Bca::new(beam, &stack, cfg, &ls, table)
            .unwrap()
            .run(|| SummaryTally::new(NM, 8))
            .unwrap();
        worst = worst.max(t.max_relative_residual);
        histories += t.histories;
    }
    out.push(Check::at_most(
        "engine.energy_conservation",
        format!(
            "Per-history energy budget residual / incident, worst of {histories} Ar 3 keV histories (SiO2/Si film, cascades, both free paths and electronic modes, with and without weak collisions)"
        ),
        worst,
        1e-9,
        sci,
        "deposited + escaped + bound + at rest = incident; rounding only",
    ));

    // Determinism: the full IonReport, bit-identical on 1, 2 and 8 threads.
    let stack = Stack::new(vec![(sio2(), 2.0 * NM)], Some(si())).unwrap();
    let mut cfg = BcaConfig::new(5.0, 1.0);
    cfg.seed = 0xC0FFEE;
    cfg.chunk_size = 16;
    let beam = Beam {
        ion: Ion::new(18).unwrap(),
        energy_ev: 2.0e3,
        polar_rad: 0.4,
        azimuth_rad: 0.0,
        count: 400,
    };
    let tc = IonTallyConfig {
        depth: Binning::new(0.0, 20.0 * NM, 40).unwrap(),
        lateral: Binning::new(-10.0 * NM, 10.0 * NM, 20).unwrap(),
        radial: Binning::new(0.0, 10.0 * NM, 10).unwrap(),
        escape_energy: Binning::new(0.0, 2000.0, 20).unwrap(),
        escape_polar: Binning::new(0.0, 0.5 * PI, 9).unwrap(),
    };
    let mut same = true;
    for weak in [0, 3] {
        cfg.weak_collisions = weak;
        let bca = Bca::new(beam, &stack, cfg, &ls, table).unwrap();
        let proto = IonTally::new(&stack, &bca.species_z(), tc).unwrap();
        let report = |n: usize| {
            with_threads(n, || {
                let r = bca.run(|| proto.clone()).unwrap().report(true);
                serde_json::to_string(&r).unwrap()
            })
        };
        let r1 = report(1);
        same &= [2, 8].iter().all(|&n| report(n) == r1);
    }
    out.push(Check::holds(
        "engine.determinism",
        "IonReport (moments, Pearson fits, damage, escapes) of a 400-ion Ar 2 keV cascade run, without and with K = 3 weak collisions, byte-identical JSON on 1, 2 and 8 threads",
        same,
        if same { "identical" } else { "DIFFERS" },
        "counter-based streams keyed on (seed, index), chunk-order merge (CONTRIBUTING.md)",
    ));

    // Sputter yield vs E_d (PR #46 review): a transfer between E_s and E_d
    // near the surface never moves, so it cannot sputter.
    let count = if quick { 600 } else { 3000 };
    let yield_at = |e_d: f64| {
        let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        m.set_displacement_energy_ev(14, e_d).unwrap();
        m.set_lattice_binding_energy_ev(14, 0.0).unwrap();
        let e_s = m.surface_binding_energy_ev(14).unwrap();
        let mut cfg = BcaConfig::new(2.0, 0.5);
        cfg.seed = 46;
        let t = Bca::new(
            Beam::normal(Ion::new(18).unwrap(), 1.0e3, count),
            &Stack::semi_infinite(m),
            cfg,
            &ls,
            table,
        )
        .unwrap()
        .run(|| SummaryTally::new(NM, 4))
        .unwrap();
        (t.sputtered as f64 / count as f64, e_s)
    };
    let (y15, e_s) = yield_at(15.0);
    let (ys, _) = yield_at(e_s);
    out.push(Check::info(
        "engine.sputter_ed_sensitivity",
        format!(
            "Ar 1 keV -> Si sputter yield with E_d = 15 eV vs E_d = E_s = {e_s} eV (E_b = 0), {count} ions"
        ),
        format!("{} vs {}", num(y15), num(ys)),
        "known model sensitivity (BH80 displacement criterion); set E_d <= E_s when yields matter (Eckstein 1991). Measured yields are compared at level 3",
    ));
    out
}
