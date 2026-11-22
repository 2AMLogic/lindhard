//! `[beam.divergence]` through the shared input, driver and summary: a TOML
//! run equals the engine set up by hand with `Bca::with_divergence`, is
//! independent of the thread count, leaves inputs without the table
//! untouched, records the resolved spread, and bad settings fail with the
//! field name (issue #285).

use lindhard::ion::crystal::Divergence;
use lindhard_cli::sim::simulate;

mod common;
use common::{check, direct_simulation, resolve, si_001_target, summary};

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
        let b = direct_simulation(&r, 1, |bca| {
            bca.with_divergence(div)
                .unwrap()
                .with_crystal(si_001_target(), &[1])
                .unwrap()
        });
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
