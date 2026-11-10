//! A photon run must be bit-identical on 1, 2 and 8 threads.

use lindhard::geometry::Stack;
use lindhard::material::Material;
use lindhard::photon::{PhotonSource, PhotonTable, PhotonTransport};
use lindhard::tally::PhotonTally;

fn run(threads: usize) -> PhotonTally {
    let table = PhotonTable::epdl2017();
    let al = Material::from_atom_fractions(&[(13, 1.0)], Some(2699.0)).unwrap();
    let cu = Material::from_atom_fractions(&[(29, 1.0)], Some(8960.0)).unwrap();
    let stack = Stack::new(vec![(al.clone(), 2.0e-3), (cu, 5.0e-4)], Some(al)).unwrap();
    let tr = PhotonTransport::new(&stack, &table).unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    pool.install(|| {
        tr.run(&PhotonSource::beam(1.2e5), 0xC0FFEE, 20_000, 64)
            .unwrap()
    })
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn photon_tally_is_bit_identical_across_thread_counts() {
    let r1 = run(1);
    assert!(r1.total_deposited_ev() > 0.0);
    for n in [2, 8] {
        let r = run(n);
        assert_eq!(r1.n_histories, r.n_histories);
        assert_eq!(
            r1.incident_ev.to_bits(),
            r.incident_ev.to_bits(),
            "{n} threads"
        );
        assert_eq!(
            bits(&r1.photoelectric_ev),
            bits(&r.photoelectric_ev),
            "{n} threads"
        );
        assert_eq!(bits(&r1.compton_ev), bits(&r.compton_ev), "{n} threads");
        assert_eq!(bits(&r1.cutoff_ev), bits(&r.cutoff_ev), "{n} threads");
        for (a, b) in [
            (&r1.front, &r.front),
            (&r1.back, &r.back),
            (&r1.side, &r.side),
        ] {
            assert_eq!((a.n, a.n_uncollided), (b.n, b.n_uncollided));
            assert_eq!(a.energy_ev.to_bits(), b.energy_ev.to_bits());
        }
        assert_eq!(
            (
                r1.n_photoelectric,
                r1.n_incoherent,
                r1.n_coherent,
                r1.n_cutoff
            ),
            (r.n_photoelectric, r.n_incoherent, r.n_coherent, r.n_cutoff)
        );
    }
}
