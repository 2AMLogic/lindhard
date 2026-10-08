//! Acceptance tests for the unit-cell lattice neighbour search
//! (`lindhard::ion::crystal::search`, issue #20): exact agreement with an
//! O(N) brute-force scan, edge cases, and determinism.

use lindhard::ion::crystal::search::{compare, path_metrics, site_position};
use lindhard::ion::crystal::{Candidate, CrystalError, Lattice, LatticeSearch};
use lindhard::rng::{stream, ParticleRng};
use rand_core::Rng;

fn uniform(rng: &mut ParticleRng) -> f64 {
    (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64
}

/// O(N) reference: scan every site of the cell block `lo..hi` per axis.
fn brute(
    lat: &Lattice,
    lo: [i64; 3],
    hi: [i64; 3],
    origin: [f64; 3],
    direction: [f64; 3],
    p_max: f64,
    length: f64,
) -> Vec<Candidate> {
    let a = lat.lattice_constant();
    let sites = lat.conventional_cell_sites();
    let n = direction.iter().map(|x| x * x).sum::<f64>().sqrt();
    let d = direction.map(|x| x / n);
    let mut out = Vec::new();
    for i in lo[0]..hi[0] {
        for j in lo[1]..hi[1] {
            for k in lo[2]..hi[2] {
                for (bi, &(z, f)) in sites.iter().enumerate() {
                    let cell = [i, j, k];
                    let pos = site_position(a, cell, f);
                    let (s, p2) = path_metrics(origin, d, pos);
                    if s >= 0.0 && s <= length && p2 <= p_max * p_max {
                        out.push(Candidate {
                            z,
                            cell,
                            basis: bi as u8,
                            position: pos,
                            s,
                            p: p2.sqrt(),
                        });
                    }
                }
            }
        }
    }
    out.sort_by(compare);
    out
}

fn in_block(c: &Candidate, lo: [i64; 3], hi: [i64; 3]) -> bool {
    (0..3).all(|k| c.cell[k] >= lo[k] && c.cell[k] < hi[k])
}

fn search_in_block(
    lat: &Lattice,
    lo: [i64; 3],
    hi: [i64; 3],
    ray: ([f64; 3], [f64; 3], f64, f64),
) -> Vec<Candidate> {
    let mut v = LatticeSearch::new(lat)
        .search(ray.0, ray.1, ray.2, ray.3)
        .unwrap();
    v.retain(|c| in_block(c, lo, hi));
    v
}

fn oracle(lat: &Lattice, seed: u64, n_rays: u64) -> usize {
    let a = lat.lattice_constant();
    let (lo, hi) = ([0i64; 3], [10i64; 3]);
    let search = LatticeSearch::new(lat);
    let mut total = 0;
    for idx in 0..n_rays {
        let mut rng = stream(seed, idx);
        let origin = [0, 1, 2].map(|_| a * (10.0 * uniform(&mut rng)));
        let dir = [0, 1, 2].map(|_| 2.0 * uniform(&mut rng) - 1.0);
        if dir.iter().all(|x| x.abs() < 1e-3) {
            continue;
        }
        let p_max = a * (0.02 + 1.5 * uniform(&mut rng));
        // Mix of sub-cell and many-cell segments.
        let length = a * (0.01 + 12.0 * uniform(&mut rng).powi(2));
        let want = brute(lat, lo, hi, origin, dir, p_max, length);
        let mut got = search.search(origin, dir, p_max, length).unwrap();
        got.retain(|c| in_block(c, lo, hi));
        assert_eq!(got, want, "ray {idx}: {origin:?} {dir:?} {p_max} {length}");
        total += want.len();
    }
    total
}

#[test]
fn matches_brute_force_silicon() {
    let n = oracle(&Lattice::silicon(), 20, 10_000);
    assert!(
        n > 10_000,
        "oracle must see non-trivial candidate sets: {n}"
    );
}

#[test]
fn matches_brute_force_gaas() {
    let lat = Lattice::gallium_arsenide();
    let n = oracle(&lat, 21, 10_000);
    assert!(n > 10_000, "{n}");
    // Both species appear.
    let a = lat.lattice_constant();
    let v = LatticeSearch::new(&lat)
        .search([0.3 * a, 0.2 * a, 0.1 * a], [1.0, 2.0, 3.0], a, 5.0 * a)
        .unwrap();
    assert!(v.iter().any(|c| c.z == 31) && v.iter().any(|c| c.z == 33));
}

fn sorted_by_s(v: &[Candidate]) {
    for w in v.windows(2) {
        assert!(compare(&w[0], &w[1]).is_le());
    }
}

#[test]
fn rays_along_110_and_100_with_many_equal_p() {
    let lat = Lattice::silicon();
    let a = lat.lattice_constant();
    let (lo, hi) = ([-6i64; 3], [12i64; 3]);
    for dir in [
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 1.0],
        [1.0, -1.0, 0.0],
    ] {
        // Start on a lattice site and at channel-centre-like offsets.
        for origin in [
            [0.0; 3],
            [0.25 * a, 0.25 * a, 0.25 * a],
            [0.5 * a, 0.0, 0.5 * a],
            [0.125 * a, 0.0, 0.0],
        ] {
            for (p_max, len) in [(0.3 * a, 4.0 * a), (a, 7.0 * a), (0.5 * a, 0.0)] {
                let ray = (origin, dir, p_max, len);
                let got = search_in_block(&lat, lo, hi, ray);
                let want = brute(&lat, lo, hi, origin, dir, p_max, len);
                assert_eq!(got, want, "{dir:?} {origin:?} {p_max} {len}");
                sorted_by_s(&got);
            }
        }
    }
    // Along [100] from the first atom of cell (0,0,0): the first candidate is
    // the start atom itself (s = 0, p = 0), and ties in s are ordered by cell
    // then basis.
    let p0 = site_position(a, [0, 0, 0], lat.conventional_cell_sites()[0].1);
    let v = LatticeSearch::new(&lat)
        .search(p0, [1.0, 0.0, 0.0], 0.4 * a, 2.0 * a)
        .unwrap();
    assert_eq!(
        (v[0].s, v[0].p, v[0].cell, v[0].basis),
        (0.0, 0.0, [0, 0, 0], 0)
    );
    let tied = v.windows(2).filter(|w| w[0].s == w[1].s).count();
    assert!(tied > 0, "expected equal-s sites");
}

#[test]
fn grazing_site_at_exactly_p_max_and_segment_end() {
    // a = 1 m makes every position a dyadic rational, so s and p are exact.
    let lat = Lattice::diamond(14, 1.0, 300.0).unwrap();
    let s = LatticeSearch::new(&lat);
    // Site (1/8, 1/8, 1/8) of cell (0,0,0). Ray along +x from (-3/8, -1/8,
    // 1/8): the site is at s = 1/2, p = 1/4 exactly.
    let origin = [-0.375, -0.125, 0.125];
    let has = |origin: [f64; 3], p_max: f64, len: f64| {
        s.search(origin, [3.0, 0.0, 0.0], p_max, len)
            .unwrap()
            .iter()
            .any(|c| c.cell == [0, 0, 0] && c.position == [0.125; 3])
    };
    assert!(has(origin, 0.25, 0.5), "closed at p_max and at s = L");
    assert!(
        !has(origin, 0.25 * (1.0 - 1e-12), 0.5),
        "inside p_max excludes"
    );
    assert!(
        !has(origin, 0.25, 0.5 * (1.0 - 1e-12)),
        "short segment excludes"
    );
    // Closed at the start (s = 0); a site behind the start is excluded.
    assert!(has([0.125; 3], 0.01, 0.1));
    assert!(!has([0.126, 0.125, 0.125], 0.01, 0.1));
    let v = s
        .search([0.126, 0.125, 0.125], [1.0, 0.0, 0.0], 0.3, 0.5)
        .unwrap();
    assert!(v.iter().all(|c| c.s >= 0.0));
}

#[test]
fn short_and_long_segments() {
    let lat = Lattice::gallium_arsenide();
    let a = lat.lattice_constant();
    // Shorter than a cell, and a zero-length segment on an atom.
    let (lo, hi) = ([-2i64; 3], [4i64; 3]);
    let ray = ([0.1 * a, 0.0, 0.0], [1.0, 0.3, 0.2], 0.4 * a, 0.2 * a);
    let got = search_in_block(&lat, lo, hi, ray);
    assert_eq!(got, brute(&lat, lo, hi, ray.0, ray.1, ray.2, ray.3));
    let one = LatticeSearch::new(&lat)
        .search([0.0; 3], [0.0, 0.0, 1.0], 1e-3 * a, 0.0)
        .unwrap();
    assert_eq!(one.len(), 1);
    // Many cells long.
    let (lo, hi) = ([0i64, -4, -4], [60i64, 5, 5]);
    for dir in [[1.0, 0.0, 0.0], [1.0, 0.05, -0.02], [1.0, 1.0, 0.0]] {
        let ray = ([0.0, 0.3 * a, 0.2 * a], dir, 0.7 * a, 50.0 * a);
        let got = search_in_block(&lat, lo, hi, ray);
        assert_eq!(got, brute(&lat, lo, hi, ray.0, ray.1, ray.2, ray.3));
        assert!(!got.is_empty());
        if dir[1] < 0.1 {
            assert!(got.len() > 100);
        }
    }
    // No candidates: empty, not an error.
    let none = LatticeSearch::new(&lat)
        .search(
            [0.125 * a, 0.125 * a, 0.125 * a],
            [1.0, 0.0, 0.0],
            1e-3 * a,
            a / 8.0,
        )
        .unwrap();
    assert!(none.is_empty());
}

#[test]
fn far_from_origin_works_with_constant_memory() {
    let lat = Lattice::silicon();
    let a = lat.lattice_constant();
    let off = 1.0e6 * a;
    let s = LatticeSearch::new(&lat);
    let near = s
        .search([0.0; 3], [1.0, 1.0, 0.0], 0.4 * a, 3.0 * a)
        .unwrap();
    let far = s
        .search([off, off, 0.0], [1.0, 1.0, 0.0], 0.4 * a, 3.0 * a)
        .unwrap();
    assert!(!far.is_empty() && (far.len() as i64 - near.len() as i64).abs() <= 2);
}

#[test]
fn rejects_bad_arguments() {
    let s = LatticeSearch::new(&Lattice::silicon());
    let d = [1.0, 0.0, 0.0];
    let o = [0.0; 3];
    let bad = |r: Result<Vec<Candidate>, CrystalError>| {
        assert!(matches!(r, Err(CrystalError::InvalidSearch { .. })))
    };
    bad(s.search(o, [0.0; 3], 1e-10, 1e-9));
    bad(s.search(o, [f64::NAN, 0.0, 0.0], 1e-10, 1e-9));
    bad(s.search(o, d, 0.0, 1e-9));
    bad(s.search(o, d, -1e-10, 1e-9));
    bad(s.search(o, d, 1e-10, -1e-9));
    bad(s.search(o, d, 1e-10, f64::INFINITY));
    bad(s.search([f64::NAN, 0.0, 0.0], d, 1e-10, 1e-9));
}

#[test]
fn identical_across_thread_counts() {
    let lat = Lattice::silicon();
    let a = lat.lattice_constant();
    let search = LatticeSearch::new(&lat);
    let run = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            use rayon::prelude::*;
            (0..200u64)
                .into_par_iter()
                .map(|i| {
                    let mut rng = stream(7, i);
                    let o = [0, 1, 2].map(|_| a * 5.0 * uniform(&mut rng));
                    let d = [0, 1, 2].map(|_| uniform(&mut rng) - 0.5);
                    search.search(o, d, 0.5 * a, 4.0 * a).unwrap()
                })
                .collect::<Vec<_>>()
        })
    };
    let one = run(1);
    assert_eq!(one, run(2));
    assert_eq!(one, run(2));
}
