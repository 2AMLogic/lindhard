//! A `[[target.solids]]` cuboid target against the equivalent
//! `[[target.layers]]` stack: the same physics on two descriptions of the
//! same solid.
//!
//! A stack is laterally infinite and (without a substrate) finite in depth;
//! a cuboid is finite in both, so the two agree only where the cuboid's
//! faces are out of reach: the cuboid is wide enough that lateral straggle
//! does not reach its sides, and the stack is a layers-only stack (no
//! substrate) whose thickness equals the cuboid's depth, so both have a back
//! face. Both runs use the same seed, so their differences are the
//! floating-point path divergence of the two geometry back ends, not
//! statistics; the tolerances below are far above that noise and far below
//! any systematic error. The sizing is asserted, not assumed:
//! `budget.lateral` must be zero (no side escape) and `budget.transmitted`
//! must match the stack's, so an undersized cuboid fails loudly instead of
//! biasing the depth profile quietly.

use lindhard::input::{Input, Resolved, Target};
use lindhard::tally::{IonReport, IonTally};

/// A 5 keV B run into a 30 nm Si target (a layers-only stack, or the same
/// slab as a cuboid): the projected range of 5 keV B in Si is inside the
/// slab, with a good fraction transmitted through the back face.
fn toml(target: &str) -> String {
    format!(
        r#"
[beam]
ion = "B"
energy_ev = 5000.0
tilt_deg = 7.0

[target]
{target}

[physics]
primary_cutoff_ev = 5.0
recoil_cutoff_ev = 2.0
[physics.energies.Si]
e_d_ev = 15.0

[run]
ions = 4000
seed = 42

[tally]
depth_bin_nm = 1.0
depth_bins = 30
"#
    )
}

const STACK_BODY: &str = r#"
[[target.layers]]
material = "Si"
thickness_nm = 30.0
"#;

const CUBOID_BODY: &str = r#"
[[target.solids]]
shape = "cuboid"
material = "Si"
lo_nm = [0.0, -5000.0, -5000.0]
hi_nm = [30.0, 5000.0, 5000.0]
"#;

fn run(text: &str) -> (IonReport, Resolved) {
    let input = Input::from_toml_str(text).unwrap();
    let r = input.resolve().unwrap();
    let table = r.scattering_table();
    let stopping = r.stopping_model();
    let bca = r.bca(&*stopping, &table).unwrap();
    let proto = IonTally::new(
        r.geometry(),
        &bca.species_z(),
        r.ion_tally_config().unwrap(),
    )
    .unwrap();
    let tally = bca.run(|| proto.clone()).unwrap();
    (tally.report(false), r)
}

#[test]
fn cuboid_target_matches_the_equivalent_stack() {
    let (stack, rs) = run(&toml(STACK_BODY));
    let (csg, rc) = run(&toml(CUBOID_BODY));
    assert!(matches!(rs.target, Target::Stack(_)));
    assert!(matches!(rc.target, Target::Csg(_)));

    // The sizing assertions: the cuboid is wide enough that nothing escapes
    // a side face, and the back face transmits like the stack's.
    assert_eq!(csg.budget.lateral, 0.0, "a side escape biases the profile");
    let t_stack = stack.budget.transmitted;
    let t_csg = csg.budget.transmitted;
    assert!(t_stack > 0.0, "the test needs transmitted primaries");
    assert!(
        (t_csg - t_stack).abs() <= 0.02 * t_stack,
        "transmitted {t_csg} vs stack {t_stack}"
    );

    // The depth profiles agree within the stated tolerance: each histogram
    // bin within 2 counts or 2% of the stack's, and the mean stopped depth
    // within 0.5 nm.
    let a = &stack.range.depth_histogram;
    let b = &csg.range.depth_histogram;
    assert_eq!(a.binning, b.binning, "same depth grid");
    assert_eq!(a.counts.len(), b.counts.len());
    for (i, (&ca, &cb)) in a.counts.iter().zip(&b.counts).enumerate() {
        assert!(
            (ca as f64 - cb as f64).abs() <= 2.0f64.max(0.02 * ca as f64),
            "depth bin {i}: cuboid {cb} vs stack {ca}"
        );
    }
    let mean = |r: &IonReport| {
        r.range
            .depth
            .as_ref()
            .map(|s| s.mean)
            .expect("stopped primaries")
    };
    assert!(
        (mean(&stack) - mean(&csg)).abs() <= 0.5e-9,
        "mean depth {} vs {}",
        mean(&stack),
        mean(&csg)
    );
    assert_eq!(stack.range.layers.len(), 1);
    assert_eq!(csg.range.layers.len(), 1, "one region: the single cuboid");
}
