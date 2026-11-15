//! `[beam.divergence]` through the shared input, driver and summary: a TOML
//! run equals the engine set up by hand with `Bca::with_divergence`, is
//! independent of the thread count, leaves inputs without the table
//! untouched, records the resolved spread, and bad settings fail with the
//! field name (issue #285).

use lindhard::input::Input;
use lindhard::ion::bca::{Bca, CrystalTarget};
use lindhard::ion::crystal::{Divergence, Lattice, Orientation};
use lindhard::ion::potential::Potential;
use lindhard::ion::scattering::ScatteringTable;
use lindhard::tally::IonTally;
use lindhard_cli::output::RunInfo;
use lindhard_cli::sim::{simulate, Simulation};
use lindhard_cli::tally::CliTally;

mod common;
use common::{check, summary};

/// A 20 nm amorphous Si screen over a crystalline Si substrate.
const BASE: &str = r#"
[beam]
ion = "B"
energy_ev = 5000.0
tilt_deg = 7.0
azimuth_deg = 22.0

[target]
substrate = "Si"
[[target.layers]]
material = "Si"
thickness_nm = 20.0

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 150
seed = 7

[tally]
depth_bin_nm = 2.0
depth_bins = 100

[[crystal]]
layers = [1]
preset = "Si"
normal = [0, 0, 1]
reference = [0, 1, 0]
"#;

fn with_divergence(table: &str) -> String {
    BASE.replace(
        "azimuth_deg = 22.0\n",
        &format!("azimuth_deg = 22.0\n\n[beam.divergence]\n{table}\n"),
    )
}

const GAUSS: &str = "model = \"gaussian\"\nsigma_deg = 0.5";
const CONE: &str = "model = \"uniform-cone\"\nhalf_angle_deg = 1.0";

fn resolve(text: &str) -> lindhard::input::Resolved {
    Input::from_toml_str(text).unwrap().resolve().unwrap()
}

/// The engine set up by hand with the divergence given in radians.
fn direct(r: &lindhard::input::Resolved, div: Divergence, threads: usize) -> Simulation {
    let pot = Potential::new(r.screening, f64::from(r.beam.ion.z()), 14.0)
        .with_length(r.screening_length);
    let table = ScatteringTable::build(&pot, &r.table_spec);
    let lattice = Lattice::silicon();
    let orientation = Orientation::new(
        &lattice,
        [0, 0, 1],
        [0, 1, 0],
        7f64.to_radians(),
        22f64.to_radians(),
        0.0,
    )
    .unwrap();
    let target = CrystalTarget::new(lattice, orientation);
    let stopping = r.stopping_model();
    let bca = Bca::new(r.beam, &r.stack, r.config, &*stopping, &table)
        .unwrap()
        .with_divergence(div)
        .unwrap()
        .with_crystal(target, &[1])
        .unwrap();
    let crystals = bca.crystal_metadata();
    let proto = IonTally::new(&r.stack, &bca.species_z(), r.ion_tally_config().unwrap()).unwrap();
    let t = &r.input.tally;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let tally = pool
        .install(|| {
            bca.run(|| {
                CliTally::new(
                    t.depth_bin_nm * 1e-9,
                    t.depth_bins,
                    t.per_ion,
                    proto.clone(),
                )
            })
        })
        .unwrap();
    let report = tally.ion.report(t.dual_pearson);
    Simulation {
        tally,
        table,
        report,
        crystals,
        info: RunInfo {
            threads,
            table_build_s: 0.0,
            transport_s: 0.0,
            ions_per_s: 0.0,
        },
    }
}

#[test]
fn toml_divergence_run_equals_direct_engine_setup() {
    for (toml, div) in [
        (
            GAUSS,
            Divergence::Gaussian {
                sigma_rad: 0.5f64.to_radians(),
            },
        ),
        (
            CONE,
            Divergence::UniformCone {
                half_angle_rad: 1.0f64.to_radians(),
            },
        ),
    ] {
        let r = resolve(&with_divergence(toml));
        assert_eq!(r.divergence, div);
        let a = simulate(&r, Some(1)).unwrap();
        let b = direct(&r, div, 1);
        assert_eq!(a.tally.finals, b.tally.finals, "{toml}");
        assert_eq!(summary(&r, &a), summary(&r, &b), "{toml}");
        // Thread independence.
        let many = simulate(&r, Some(3)).unwrap();
        assert_eq!(a.tally.finals, many.tally.finals, "{toml}");
        assert_eq!(summary(&r, &a), summary(&r, &many), "{toml}");
    }
}

#[test]
fn absent_divergence_is_unchanged_and_present_one_is_recorded() {
    let base = resolve(BASE);
    let a = simulate(&base, Some(1)).unwrap();
    let sa = summary(&base, &a);
    assert!(sa["physics"].get("beam_divergence").is_none());
    assert!(sa["input"]["beam"].get("divergence").is_none());
    assert_eq!(sa["format"]["version"], 1);
    // A zero-width spread is a valid setting and still changes no scientific
    // output (every deflection is the nominal direction).
    let zero = resolve(&with_divergence("model = \"gaussian\"\nsigma_deg = 0.0"));
    let z = simulate(&zero, Some(1)).unwrap();
    assert_eq!(a.tally.finals, z.tally.finals);

    let r = resolve(&with_divergence(GAUSS));
    let s = simulate(&r, Some(1)).unwrap();
    assert_ne!(a.tally.finals, s.tally.finals, "divergence was applied");
    let v = summary(&r, &s);
    let m = &v["physics"]["beam_divergence"];
    assert_eq!(m["model"], "gaussian");
    assert_eq!(m["width_kind"], "sigma_per_plane");
    assert_eq!(
        m["width_deg"].as_f64().unwrap(),
        0.5f64.to_radians().to_degrees()
    );
    assert_eq!(m["incidence"], "inward-conditioned");
    assert_eq!(v["input"]["beam"]["divergence"]["model"], "gaussian");
    assert_eq!(v["input"]["beam"]["divergence"]["sigma_deg"], 0.5);
    assert_eq!(v["format"]["version"], 1);
}

#[test]
fn check_accepts_valid_and_names_the_field_of_invalid_settings() {
    let (ok, out) = check("div-ok", &with_divergence(GAUSS));
    assert!(ok, "{out}");
    assert!(out.contains("beam divergence: gaussian"), "{out}");

    let bad = |name: &str, text: String, want: &str| {
        let (ok, out) = check(name, &text);
        assert!(!ok, "{name}: should fail\n{out}");
        assert!(out.contains(want), "{name}: expected {want:?} in\n{out}");
    };
    for (name, w) in [
        ("neg", "-0.1"),
        ("nan", "nan"),
        ("inf", "inf"),
        ("wide", "45.0"),
    ] {
        bad(
            name,
            with_divergence(&format!("model = \"gaussian\"\nsigma_deg = {w}")),
            "beam.divergence.sigma_deg",
        );
        bad(
            name,
            with_divergence(&format!("model = \"uniform-cone\"\nhalf_angle_deg = {w}")),
            "beam.divergence.half_angle_deg",
        );
    }
    bad(
        "unknown-model",
        with_divergence("model = \"lorentzian\"\nsigma_deg = 0.1"),
        "lorentzian",
    );
    bad(
        "wrong-key",
        with_divergence("model = \"gaussian\"\nhalf_angle_deg = 0.1"),
        "half_angle_deg",
    );
    bad(
        "dynamic",
        with_divergence(GAUSS).replace(
            "[run]\n",
            "[dynamic]\nfluence_cm2 = 1e14\nions_per_step = 50\n[run]\n",
        ),
        "beam.divergence",
    );
}
