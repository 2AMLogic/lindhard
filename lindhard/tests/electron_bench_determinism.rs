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

#[path = "../benches/electron_common/mod.rs"]
mod electron_common;

use electron_common as ec;

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
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    pool.install(|| {
        let tables = ec::tables(&r).unwrap();
        let (t, proto) = ec::transport(&r, &tables).unwrap();
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
