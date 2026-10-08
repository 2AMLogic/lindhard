//! Lattice neighbour search throughput: candidates found per second for a
//! typical search radius and segment, in Si and GaAs.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use std::hint::black_box;

use lindhard::ion::crystal::{Lattice, LatticeSearch};

fn bench(c: &mut Criterion) {
    let mut g = c.benchmark_group("crystal_search");
    for (name, lat) in [
        ("silicon", Lattice::silicon()),
        ("gaas", Lattice::gallium_arsenide()),
    ] {
        let a = lat.lattice_constant();
        let search = LatticeSearch::new(&lat);
        // Typical flight: p_max ~ 0.4 a, segment ~ 3 a; one off-axis
        // direction and one along <110>.
        for (dname, dir) in [("random", [0.31, 0.52, 0.79]), ("110", [1.0, 1.0, 0.0])] {
            let origin = [0.13 * a, 0.29 * a, 0.41 * a];
            let (p_max, len) = (0.4 * a, 3.0 * a);
            let n = search.search(origin, dir, p_max, len).unwrap().len();
            g.throughput(Throughput::Elements(n as u64));
            let mut buf = Vec::new();
            g.bench_function(format!("{name}_{dname}"), |b| {
                b.iter(|| {
                    search
                        .search_into(black_box(origin), black_box(dir), p_max, len, &mut buf)
                        .unwrap();
                    black_box(buf.len())
                })
            });
        }
    }
    g.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
