//! Run the electron transport loop with the full electron tally and print the
//! report.
//!
//! ```text
//! cargo run --release -p lindhard --example electron_tally            # summary
//! cargo run --release -p lindhard --example electron_tally -- --json  # full report
//! ```
//!
//! The cross-section tables are **synthetic**, computed here from simple
//! formulas: constant mean free paths, isotropic elastic scattering, and
//! inelastic losses uniform on `[0, 0.2 E]`. They are not physical data, so
//! the numbers only show what the tally reports, not what a real material
//! does. Real tables plug in through the same `LayerTables` once the elastic
//! and inelastic models land.

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::{LayerTables, Primary, Transport, TransportConfig};
use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::tally::{
    Binning, CartesianGrid, CylindricalGrid, ElectronTallyConfig, FullElectronTally,
};

const NM: f64 = 1e-9;

fn table(axis: SamplingAxis, lambda_m: f64, row: impl Fn(f64, f64) -> f64) -> CrossSectionTable {
    let energy_ev: Vec<f64> = (0..=40).map(|i| 10f64.powf(i as f64 / 10.0)).collect();
    let probability: Vec<f64> = (0..=20).map(|i| i as f64 / 20.0).collect();
    let quantiles = energy_ev
        .iter()
        .map(|&e| probability.iter().map(|&p| row(e, p)).collect())
        .collect();
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "computed in lindhard/examples/electron_tally.rs (not physical data)".into(),
        axis,
        inverse_mfp_per_m: vec![1.0 / lambda_m; energy_ev.len()],
        energy_ev,
        probability,
        quantiles,
    })
    .expect("valid synthetic table")
}

fn main() {
    let json = std::env::args().any(|a| a == "--json");
    let e0 = 1_000.0;

    let si = Material::from_atom_fractions(&[(14, 1.0)], None).expect("valid material");
    let tables = LayerTables {
        elastic: table(SamplingAxis::ElasticPolarAngle, 2.0 * NM, |_, p| {
            (1.0 - 2.0 * p).acos()
        }),
        inelastic: table(SamplingAxis::InelasticEnergyLoss, 3.0 * NM, |e, p| {
            0.2 * e * p
        }),
    };
    let transport = Transport::new(
        Stack::semi_infinite(si),
        vec![tables],
        TransportConfig::new(5.0),
    )
    .expect("valid transport");

    let mut config = ElectronTallyConfig::new(
        Binning::new(0.0, e0, 50).expect("valid binning"),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 18).expect("valid binning"),
    );
    config.cartesian = Some(CartesianGrid {
        x: Binning::new(0.0, 30.0 * NM, 30).expect("valid binning"),
        y: Binning::new(-15.0 * NM, 15.0 * NM, 30).expect("valid binning"),
        z: Binning::new(-15.0 * NM, 15.0 * NM, 30).expect("valid binning"),
    });
    config.cylindrical = Some(CylindricalGrid {
        r: Binning::new(0.0, 15.0 * NM, 30).expect("valid binning"),
        depth: Binning::new(0.0, 30.0 * NM, 30).expect("valid binning"),
    });
    let proto = FullElectronTally::new(&transport, config).expect("valid tally configuration");

    let run = transport
        .run(2026, 20_000, 256, &Primary::normal(e0), || proto.clone())
        .expect("valid primary");
    let report = run.tally.report();

    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&report).expect("serializable report")
        );
        return;
    }

    let b = &report.budget;
    println!("histories               {}", report.histories);
    println!(
        "SE/BSE split            {} eV ({})",
        report.metadata.se_bse_split_ev, report.metadata.se_bse_split_rule
    );
    println!(
        "backscatter yield eta   {:.4}",
        report.yields.backscatter_eta
    );
    println!(
        "secondary yield delta   {:.4}",
        report.yields.secondary_delta
    );
    println!("incident energy         {:.6e} eV", b.incident_ev);
    println!("  deposited             {:.6e} eV", b.deposited_ev);
    println!("  escaped               {:.6e} eV", b.escaped_ev);
    println!("  trapped               {:.6e} eV", b.trapped_ev);
    println!("  relative imbalance    {:.3e}", b.relative_imbalance);
    if let Some(g) = report.generation_volume {
        println!(
            "generation volume       mean depth {:.2} nm, depth sd {:.2} nm, rms radius {:.2} nm",
            g.mean_m[0] / NM,
            g.std_dev_m[0] / NM,
            g.rms_radius_m / NM
        );
    }
    if let Some(cyl) = &report.deposition.cylindrical {
        let inside: f64 = cyl.energy_ev.iter().sum();
        println!(
            "cylindrical grid        {:.1} % of the deposit inside",
            100.0 * inside / b.deposited_ev
        );
    }
}
