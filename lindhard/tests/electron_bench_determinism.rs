//! The electron benchmark fixture (`benches/electron_common`, #152) gives a
//! bit-identical full tally report on 1 and 2 threads, with every physics
//! choice of the benchmarks on (Mott elastic with exchange, single-pole Penn
//! inelastic, Kieft-Bosch secondaries, step barrier, vacuum-level cutoff).
//! `examples/electron_scaling.rs` checks the same digest at every thread count
//! it measures; this keeps the fixture itself honest in CI. The table grid is
//! coarse (4 points per decade) so the test is quick in a debug build; the
//! determinism does not depend on the grid.
//!
//! The thread counts are 1 and 2 (the engine's 1, 2 and 8 thread checks are in
//! `tests/electron_secondaries.rs` and `tests/determinism.rs`).
//!
//! The second test guards the fixture against the drift of #750: the fixture
//! must build its tables with `lindhard::input::electron::assemble`, the
//! implementation `lindhard run` uses, and not with a copy of its own. The
//! fixture's materials all have a band, so the shared builder puts their
//! inelastic tables on the band-bottom axis and rate-refines them (#339),
//! which adds rows between 2 and 20 eV above the Fermi level. A local copy
//! that leaves `InelasticTableOptions::rate_refinement` at its `None` default
//! — as the bench copy did from 2026-10-10 until #750 — returns exactly the
//! rows of `electron.tables` instead, so the row count tells the two apart.

#[path = "../benches/electron_common/mod.rs"]
mod electron_common;

use electron_common as ec;
use lindhard::input::electron::assemble;
use lindhard::tally::FullElectronTally;

mod common;
use common::with_threads;

fn digest_on(threads: usize, p: &ec::ElectronProblem) -> String {
    let dir = std::env::temp_dir().join(format!(
        "lindhard-test-electron-bench-{}-{}",
        std::process::id(),
        p.id
    ));
    let settings = ec::Settings {
        points_per_decade: 4.0,
        rmax_nm: 1000.0,
        rbins: 50,
    };
    let histories = 48;
    let path = ec::write_synthetic(&dir, p, &settings, histories).unwrap();
    let r = ec::resolve(&path).unwrap();
    with_threads(threads, || {
        let tables = assemble::material_tables(&r).unwrap();
        let t = assemble::transport(&r, &tables).unwrap();
        let proto = FullElectronTally::new(&t, r.tally.clone()).unwrap();
        let report = ec::run(&r, &t, &proto, histories).unwrap();
        assert_eq!(report.histories, histories);
        ec::report_digest(&report).unwrap()
    })
}

#[test]
fn bench_fixture_report_is_bit_identical_on_1_and_2_threads() {
    for id in ["e_1keV_si", "e_1keV_cu"] {
        let p = ec::problem(id).unwrap();
        let one = digest_on(1, &p);
        let two = digest_on(2, &p);
        assert_eq!(one, two, "{id}: report differs between 1 and 2 threads");
    }
}

#[test]
fn fixture_inelastic_tables_are_rate_refined_like_a_cli_run() {
    for id in ["e_1keV_si", "e_1keV_cu"] {
        let p = ec::problem(id).unwrap();
        let dir = std::env::temp_dir().join(format!(
            "lindhard-test-electron-bench-refine-{}-{}",
            std::process::id(),
            p.id
        ));
        let settings = ec::Settings {
            points_per_decade: 4.0,
            rmax_nm: 1000.0,
            rbins: 50,
        };
        let path = ec::write_synthetic(&dir, &p, &settings, 1).unwrap();
        let r = ec::resolve(&path).unwrap();
        let m = &r.materials[0];

        // The fixture's materials have a band, so a CLI run refines them.
        assert!(
            assemble::rate_refinement(&r, m).is_some(),
            "{id}: the fixture material has no band-bottom axis to refine"
        );
        let t = assemble::inelastic_table(&r, m).unwrap();
        assert!(
            t.energy_ev().len() > r.table_energy_ev.len(),
            "{id}: the inelastic table has {} rows for a {}-row grid, so it was \
             built without the rate refinement a `lindhard run` applies (#750)",
            t.energy_ev().len(),
            r.table_energy_ev.len()
        );
    }
}
