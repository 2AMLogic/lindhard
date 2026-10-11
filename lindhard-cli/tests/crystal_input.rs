//! `[[crystal]]` through the shared input, driver and summary: a TOML run
//! equals the engine set up by hand, is independent of the thread count,
//! leaves amorphous inputs and amorphous layers untouched, and is diagnosed
//! by `check` before transport.

use std::path::{Path, PathBuf};
use std::process::Command;

use lindhard::input::Input;
use lindhard::ion::bca::{Bca, CrystalTarget, Thermal};
use lindhard::ion::crystal::debye::THETA_D_SI;
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::Potential;
use lindhard::ion::scattering::ScatteringTable;
use lindhard::tally::IonTally;
use lindhard_cli::output::{self, RunInfo};
use lindhard_cli::sim::{simulate, Simulation};
use lindhard_cli::tally::{ion_tally_config, CliTally};
use serde_json::Value;

const AMORPHOUS: &str = r#"
[beam]
ion = "B"
energy_ev = 5000.0
tilt_deg = 7.0
azimuth_deg = 22.0

[target]
substrate = "Si"

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
"#;

const CRYSTAL: &str = r#"
[[crystal]]
layers = [0]
preset = "Si"
normal = [0, 0, 1]
reference = [0, 1, 0]
[crystal.thermal]
temperature_k = 300.0
"#;

fn crystal_input() -> String {
    format!("{AMORPHOUS}{CRYSTAL}")
}

fn resolve(text: &str) -> lindhard::input::Resolved {
    Input::from_toml_str(text).unwrap().resolve().unwrap()
}

fn summary(r: &lindhard::input::Resolved, s: &Simulation) -> Value {
    let text = output::summary_json(r, &s.table, &s.tally, &s.report, &s.crystals, s.info).unwrap();
    let mut v: Value = serde_json::from_str(&text).unwrap();
    v.as_object_mut().unwrap().remove("run");
    v
}

/// The engine set up by hand: Si lattice, (001) cut, [010] reference, the
/// beam's tilt and azimuth, 300 K with the cited Debye temperature.
fn direct(r: &lindhard::input::Resolved, threads: usize) -> Simulation {
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
    let target =
        CrystalTarget::new(lattice, orientation).with_thermal(Thermal::new(300.0, THETA_D_SI));
    let stopping = r.stopping_model();
    let bca = Bca::new(r.beam, &r.stack, r.config, &*stopping, &table)
        .unwrap()
        .with_crystal(target, &[0])
        .unwrap();
    let crystals = bca.crystal_metadata();
    let proto = IonTally::new(&r.stack, &bca.species_z(), ion_tally_config(r).unwrap()).unwrap();
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
fn toml_crystal_run_equals_direct_engine_setup() {
    let r = resolve(&crystal_input());
    let a = simulate(&r, Some(1)).unwrap();
    let b = direct(&r, 1);
    assert_eq!(a.crystals, b.crystals, "resolved crystal metadata");
    assert_eq!(a.tally.finals, b.tally.finals, "per-ion final states");
    assert_eq!(summary(&r, &a), summary(&r, &b), "tallies and metadata");
    // The metadata is what the input asked for.
    let m = &a.crystals[0];
    assert_eq!(m.regions, vec![0]);
    assert_eq!(m.normal_hkl, [0, 0, 1]);
    assert_eq!(m.reference_uvw, [0, 1, 0]);
    assert_eq!(m.tilt_rad.to_bits(), 7f64.to_radians().to_bits());
    assert_eq!(m.twist_rad.to_bits(), 22f64.to_radians().to_bits());
    assert_eq!(
        m.thermal.as_ref().unwrap().input.debye_temperature_k,
        THETA_D_SI
    );
}

#[test]
fn crystal_summary_is_identical_across_thread_counts() {
    let r = resolve(&crystal_input());
    let one = simulate(&r, Some(1)).unwrap();
    let many = simulate(&r, Some(3)).unwrap();
    assert_eq!(one.tally.finals, many.tally.finals);
    assert_eq!(summary(&r, &one), summary(&r, &many));
    assert_eq!(many.info.threads, 3);
}

#[test]
fn crystal_changes_the_physics_but_amorphous_input_does_not_gain_the_field() {
    let ra = resolve(AMORPHOUS);
    let rc = resolve(&crystal_input());
    let a = simulate(&ra, Some(1)).unwrap();
    let c = simulate(&rc, Some(1)).unwrap();
    assert!(a.crystals.is_empty());
    let sa = summary(&ra, &a);
    assert!(sa["physics"].get("crystal").is_none());
    assert!(sa["input"].get("crystal").is_none());
    assert_eq!(sa["format"]["version"], 1);
    assert!(!sa["physics"]["models"]
        .as_array()
        .unwrap()
        .iter()
        .any(|m| m["role"] == "crystal transport"));
    let sc = summary(&rc, &c);
    assert!(sc["physics"]["crystal"].is_array());
    assert_eq!(sc["format"]["version"], 1);
    assert_ne!(
        sa["results"]["range"], sc["results"]["range"],
        "the crystal engine was actually used"
    );
    // The amorphous summary equals the one of the amorphous engine run by hand.
    let pot = Potential::new(ra.screening, 5.0, 14.0).with_length(ra.screening_length);
    let table = ScatteringTable::build(&pot, &ra.table_spec);
    let stopping = ra.stopping_model();
    let bca = Bca::new(ra.beam, &ra.stack, ra.config, &*stopping, &table).unwrap();
    let proto = IonTally::new(&ra.stack, &bca.species_z(), ion_tally_config(&ra).unwrap()).unwrap();
    let t = &ra.input.tally;
    let tally = bca
        .run(|| {
            CliTally::new(
                t.depth_bin_nm * 1e-9,
                t.depth_bins,
                t.per_ion,
                proto.clone(),
            )
        })
        .unwrap();
    assert_eq!(tally.finals, a.tally.finals);
}

const OVERLAYER: &str = r#"
[beam]
ion = "As"
energy_ev = 500.0
tilt_deg = 7.0
azimuth_deg = 22.0

[materials.SiO2]
density_g_cm3 = 2.2
elements = [
  { symbol = "Si", atom_fraction = 1.0 },
  { symbol = "O", atom_fraction = 2.0, e_d_ev = 20.0, e_s_ev = 2.0 },
]

[target]
substrate = "Si"
[[target.layers]]
material = "SiO2"
thickness_nm = 30.0

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 100
seed = 3
"#;

#[test]
fn amorphous_overlayer_stays_amorphous_while_substrate_is_crystal() {
    let crystal = format!(
        "{OVERLAYER}\n[[crystal]]\nlayers = [1]\npreset = \"Si\"\nnormal = [0,0,1]\nreference = [0,1,0]\n"
    );
    let rc = resolve(&crystal);
    let ra = resolve(OVERLAYER);
    let c = simulate(&rc, Some(1)).unwrap();
    let a = simulate(&ra, Some(1)).unwrap();
    // Only the substrate region is a crystal.
    assert_eq!(c.crystals.len(), 1);
    assert_eq!(c.crystals[0].regions, vec![1]);
    // 500 eV As and its recoils stop within the 30 nm oxide, so no history
    // ever reaches the crystal region: the run is bit for bit the amorphous
    // one, which it can only be if the oxide kept the amorphous transport.
    assert!(c.tally.finals.iter().all(|f| f.pos[0] < 30e-9));
    assert_eq!(c.tally.finals, a.tally.finals);
    assert_eq!(
        summary(&rc, &c)["results"],
        summary(&ra, &a)["results"],
        "overlayer results unchanged"
    );
}

#[test]
fn unverified_electronic_constants_flag_is_preserved() {
    let local = crystal_input().replace(
        "[physics]\n",
        "[physics]\nstopping = \"equipartition-ls-or\"\n",
    );
    let r = resolve(&local);
    assert!(r.warnings.iter().any(|w| w.contains("not verified")));
    let s = simulate(&r, Some(1)).unwrap();
    assert!(s.crystals[0].electronic_constants_unverified);
    let v = summary(&r, &s);
    assert_eq!(
        v["physics"]["crystal"][0]["electronic_constants_unverified"],
        true
    );
    // Nonlocal-only: the flag is absent.
    let s = simulate(&resolve(&crystal_input()), Some(1)).unwrap();
    assert!(!s.crystals[0].electronic_constants_unverified);
    assert!(
        summary(&resolve(&crystal_input()), &s)["physics"]["crystal"][0]
            .get("electronic_constants_unverified")
            .is_none()
    );
}

fn scratch(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("lindhard-cli-tests")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// `lindhard check` on `text`: (success, stderr).
fn check(name: &str, text: &str) -> (bool, String) {
    let p = scratch(name).join("in.toml");
    std::fs::write(&p, text).unwrap();
    let o = Command::new(env!("CARGO_BIN_EXE_lindhard"))
        .args(["check", p.to_str().unwrap()])
        .output()
        .unwrap();
    (
        o.status.success(),
        format!(
            "{}{}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        ),
    )
}

#[test]
fn check_accepts_valid_and_names_the_key_of_invalid_assignments() {
    let (ok, out) = check("crystal-ok", &crystal_input());
    assert!(ok, "{out}");
    assert!(out.contains("crystal[0]: Si"), "{out}");

    let bad = |name: &str, from: &str, to: &str, want: &str| {
        let text = crystal_input().replace(from, to);
        assert_ne!(text, crystal_input(), "{name}: replacement applied");
        let (ok, out) = check(name, &text);
        assert!(!ok, "{name}: should fail\n{out}");
        assert!(out.contains(want), "{name}: expected {want:?} in\n{out}");
    };
    bad(
        "zero-normal",
        "normal = [0, 0, 1]",
        "normal = [0, 0, 0]",
        "crystal[0].normal",
    );
    bad(
        "ref-not-in-plane",
        "reference = [0, 1, 0]",
        "reference = [0, 0, 1]",
        "crystal[0].reference",
    );
    bad(
        "zero-ref",
        "reference = [0, 1, 0]",
        "reference = [0, 0, 0]",
        "crystal[0].reference",
    );
    bad(
        "layer-range",
        "layers = [0]",
        "layers = [3]",
        "crystal[0].layers",
    );
    bad(
        "layer-empty",
        "layers = [0]",
        "layers = []",
        "crystal[0].layers",
    );
    bad(
        "layer-dup",
        "layers = [0]",
        "layers = [0, 0]",
        "already assigned",
    );
    bad(
        "preset-unknown",
        "preset = \"Si\"",
        "preset = \"Diamond\"",
        "preset",
    );
    // Preset/material mismatch: Ge lattice on a Si layer.
    bad(
        "preset-mismatch",
        "preset = \"Si\"",
        "preset = \"Ge\"",
        "crystal[0].preset",
    );
    bad(
        "bad-temperature",
        "temperature_k = 300.0",
        "temperature_k = -1.0",
        "crystal[0].thermal.temperature_k",
    );
    bad(
        "missing-reference",
        "reference = [0, 1, 0]\n",
        "",
        "reference",
    );
    bad(
        "weak",
        "[physics]\n",
        "[physics]\nweak_collisions = 1\n",
        "physics.weak_collisions",
    );
    bad(
        "free-path",
        "[physics]\n",
        "[physics]\nfree_path = \"energy-dependent\"\nmin_cm_angle_deg = 1.0\n",
        "physics.free_path",
    );
    bad(
        "dynamic",
        "[run]\n",
        "[dynamic]\nfluence_cm2 = 1e14\nions_per_step = 50\n[run]\n",
        "not supported with a [dynamic]",
    );
}

#[test]
fn crystal_with_an_electron_run_is_rejected() {
    let text = "[electron]\n[[crystal]]\nlayers = [0]\n";
    let (ok, out) = check("crystal-electron", text);
    assert!(!ok);
    assert!(out.contains("crystal") && out.contains("electron"), "{out}");
}
