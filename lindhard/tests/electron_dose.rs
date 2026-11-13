//! Dose per layer, depth-dose and interface transmission
//! (`tally::dose`): the by-construction sums, the dose arithmetic, thread-count
//! determinism, and that the tallies off leave the report unchanged. The
//! cross-section tables are synthetic (constant mean free paths, isotropic
//! angles), as in `electron_tally.rs`.

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::{LayerTables, Primary, Transport, TransportConfig};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::tally::{
    Binning, CartesianGrid, DepthDoseConfig, DoseConfig, ElectronReport, ElectronTallyConfig,
    ElectronTallyError, FullElectronTally, InterfaceRefine,
};
use lindhard::units::J_PER_EV;

mod common;
use common::close;
use common::with_threads;

const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];
const NM: f64 = 1e-9;

fn material() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

fn table(axis: SamplingAxis, rate: f64, row: impl Fn(f64, f64) -> f64) -> CrossSectionTable {
    let quantiles = GRID
        .iter()
        .map(|&e| {
            if rate == 0.0 {
                vec![]
            } else {
                PROB.iter().map(|&p| row(e, p)).collect()
            }
        })
        .collect();
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "computed in tests/electron_tally.rs".into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: vec![rate; GRID.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

/// Isotropic elastic scattering (`cos θ = 1 - 2p`) and inelastic losses
/// uniform on `[0, frac E]`.
fn scattering(lambda_el: f64, lambda_inel: f64, frac: f64) -> LayerTables {
    LayerTables {
        elastic: table(SamplingAxis::ElasticPolarAngle, 1.0 / lambda_el, |_, p| {
            (1.0 - 2.0 * p).acos()
        }),
        inelastic: table(
            SamplingAxis::InelasticEnergyLoss,
            1.0 / lambda_inel,
            move |e, p| frac * e * p,
        ),
    }
}

fn light() -> Material {
    Material::from_atom_fractions(&[(6, 1.0)], None).unwrap()
}

fn stack() -> Stack {
    Stack::new(
        vec![(material(), 6.0 * NM), (light(), 4.0 * NM)],
        Some(material()),
    )
    .unwrap()
}

fn transport() -> Transport {
    Transport::new(
        stack(),
        vec![
            scattering(3.0 * NM, 5.0 * NM, 0.3),
            scattering(4.0 * NM, 6.0 * NM, 0.3),
            scattering(2.0 * NM, 4.0 * NM, 0.3),
        ],
        TransportConfig::new(5.0),
    )
    .unwrap()
}

fn dose_config() -> DoseConfig {
    DoseConfig {
        depth: Some(DepthDoseConfig {
            bins_per_layer: 4,
            refine: Some(InterfaceRefine {
                width_m: 1.0 * NM,
                bins: 5,
            }),
            substrate_depth_m: 60.0 * NM,
        }),
        transmission: true,
    }
}

fn config(e0: f64, dose: Option<DoseConfig>) -> ElectronTallyConfig {
    let mut c = ElectronTallyConfig::new(
        Binning::new(0.0, e0, 20).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    // Covers the finite layers only, so the grid total differs from the
    // deposited total by `outside_ev`.
    c.cartesian = Some(CartesianGrid {
        x: Binning::new(0.0, 200.0 * NM, 40).unwrap(),
        y: Binning::new(-100.0 * NM, 100.0 * NM, 4).unwrap(),
        z: Binning::new(-100.0 * NM, 100.0 * NM, 4).unwrap(),
    });
    c.dose = dose;
    c
}

fn run(t: &Transport, cfg: ElectronTallyConfig, n: u64) -> ElectronReport {
    let proto = FullElectronTally::new(t, cfg).unwrap();
    t.run(0xD05E, n, 64, &Primary::normal(2_000.0), || proto.clone())
        .unwrap()
        .tally
        .report()
}

#[test]
fn layer_depth_and_cartesian_sums_agree_by_construction() {
    let r = run(&transport(), config(2_000.0, Some(dose_config())), 3_000);
    let dose = r.dose.as_ref().unwrap();
    let depth = dose.depth.as_ref().unwrap();
    let h = r.histories as f64;
    assert_eq!(depth.edges_m.len(), depth.energy_ev.len() + 1);
    for (i, l) in dose.layers.iter().enumerate() {
        // Energy per primary against the plain per-layer total.
        assert!(close(
            l.energy_ev.mean,
            r.deposition.per_layer_ev[i] / h,
            1e-12
        ));
        let bins: f64 = depth
            .layer_of_bin
            .iter()
            .zip(&depth.energy_ev)
            .filter(|(&b, _)| b == i)
            .map(|(_, e)| e.mean)
            .sum();
        let outside = if i == 2 { depth.outside_ev.mean } else { 0.0 };
        assert!(
            close(bins + outside, l.energy_ev.mean, 1e-12),
            "layer {i}: bins {bins} + outside {outside} vs {}",
            l.energy_ev.mean
        );
    }
    let total: f64 = dose.layers.iter().map(|l| l.energy_ev.mean).sum();
    let cart = r.deposition.cartesian.as_ref().unwrap();
    assert!(close(
        total,
        (cart.energy_ev.iter().sum::<f64>() + cart.outside_ev) / h,
        1e-12
    ));
    assert!(dose.layers.iter().all(|l| l.energy_ev.std_err > 0.0));
    // The energy balance is a property of the transport, not of this tally.
    assert!(r.budget.relative_imbalance < 1e-9);
}

#[test]
fn dose_is_energy_over_areal_mass_and_substrate_has_none() {
    let r = run(&transport(), config(2_000.0, Some(dose_config())), 2_000);
    let dose = r.dose.as_ref().unwrap();
    let rho = [material(), light(), material()].map(|m| m.mass_density());
    let t = [6.0 * NM, 4.0 * NM];
    for i in 0..2 {
        let l = &dose.layers[i];
        let m = rho[i] * t[i];
        assert!(close(l.areal_mass_kg_m2.unwrap(), m, 1e-12));
        let d = l.dose_gy_m2.unwrap();
        assert!(close(d.mean, l.energy_ev.mean * J_PER_EV / m, 1e-12));
        assert!(close(d.std_err, l.energy_ev.std_err * J_PER_EV / m, 1e-12));
    }
    assert!(dose.layers[2].dose_gy_m2.is_none());
    assert!(dose.layers[2].areal_mass_kg_m2.is_none());
    let depth = dose.depth.as_ref().unwrap();
    assert!(depth
        .dose_gy_m2
        .iter()
        .zip(&depth.layer_of_bin)
        .all(|(d, &l)| d.is_some() == (l < 2)));
    // Refined bands next to the interfaces: the first layer has 5 + 4 + 5 bins.
    assert_eq!(depth.layer_of_bin.iter().filter(|&&l| l == 0).count(), 14);
    // A band is a fifth of a nanometre wide.
    assert!(close(depth.edges_m[1] - depth.edges_m[0], 0.2 * NM, 1e-9));
}

#[test]
fn transmission_counts_crossings_in_both_directions() {
    let r = run(&transport(), config(2_000.0, Some(dose_config())), 2_000);
    let tr = &r.dose.as_ref().unwrap().transmission;
    // Two interfaces, forward and back, primary and secondary.
    assert_eq!(tr.len(), 8);
    let fwd0 = tr
        .iter()
        .find(|t| t.interface == 0 && t.from_layer == 0 && !t.secondary);
    let back0 = tr
        .iter()
        .find(|t| t.interface == 0 && t.from_layer == 1 && !t.secondary);
    let (fwd0, back0) = (fwd0.unwrap(), back0.unwrap());
    assert!(fwd0.count.mean > 0.0 && back0.count.mean > 0.0);
    assert!(fwd0.energy_ev.mean > 0.0);
    // No secondary model, so no secondaries.
    assert!(tr
        .iter()
        .filter(|t| t.secondary)
        .all(|t| t.count.mean == 0.0));
    let n: u64 =
        fwd0.spectrum.counts.iter().sum::<u64>() + fwd0.spectrum.underflow + fwd0.spectrum.overflow;
    assert!(close(n as f64, fwd0.count.mean * r.histories as f64, 1e-9));
    assert!(fwd0.energy_fraction > 0.0 && fwd0.energy_fraction < 1.0);
}

#[test]
fn tallies_off_leave_the_report_unchanged() {
    let t = transport();
    let off = run(&t, config(2_000.0, None), 1_500);
    let mut on = run(&t, config(2_000.0, Some(dose_config())), 1_500);
    assert!(off.dose.is_none());
    on.dose = None;
    on.metadata.config.dose = None;
    assert_eq!(off, on);
}

#[test]
fn report_is_bit_identical_across_thread_counts() {
    let t = transport();
    let go = |threads: usize| {
        with_threads(threads, || {
            run(&t, config(2_000.0, Some(dose_config())), 3_000)
        })
    };
    let r1 = go(1);
    let r2 = go(2);
    assert_eq!(r1.dose, r2.dose);
    assert_eq!(
        serde_json::to_string(&r1).unwrap(),
        serde_json::to_string(&r2).unwrap()
    );
}

#[test]
fn one_layer_zero_width_refinement_and_coarse_bins() {
    let t = Transport::new(
        Stack::new(vec![(material(), 10.0 * NM)], None).unwrap(),
        vec![scattering(3.0 * NM, 5.0 * NM, 0.3)],
        TransportConfig::new(5.0),
    )
    .unwrap();
    let cfg = DoseConfig {
        depth: Some(DepthDoseConfig {
            bins_per_layer: 1,
            refine: Some(InterfaceRefine {
                width_m: 0.0,
                bins: 3,
            }),
            substrate_depth_m: 0.0,
        }),
        transmission: true,
    };
    let r = run(&t, config(1_000.0, Some(cfg)), 500);
    let d = r.dose.unwrap();
    assert!(d.transmission.is_empty());
    let depth = d.depth.unwrap();
    assert_eq!(depth.energy_ev.len(), 1);
    assert!(close(
        depth.energy_ev[0].mean,
        d.layers[0].energy_ev.mean,
        1e-12
    ));
    // A band wider than half the layer is limited to half.
    let wide = DoseConfig {
        depth: Some(DepthDoseConfig {
            bins_per_layer: 2,
            refine: Some(InterfaceRefine {
                width_m: 1.0,
                bins: 2,
            }),
            substrate_depth_m: 0.0,
        }),
        transmission: false,
    };
    let r = run(&t, config(1_000.0, Some(wide)), 500);
    let depth = r.dose.unwrap().depth.unwrap();
    assert!(depth.edges_m.windows(2).all(|w| w[1] > w[0]));
}

#[test]
fn invalid_dose_configuration_is_refused() {
    let t = transport();
    let mut bad = dose_config();
    bad.depth.as_mut().unwrap().bins_per_layer = 0;
    assert!(matches!(
        FullElectronTally::new(&t, config(2_000.0, Some(bad))),
        Err(ElectronTallyError::Dose(_))
    ));
    let mut bad = dose_config();
    bad.depth.as_mut().unwrap().substrate_depth_m = 0.0;
    assert!(matches!(
        FullElectronTally::new(&t, config(2_000.0, Some(bad))),
        Err(ElectronTallyError::Dose(_))
    ));
}
