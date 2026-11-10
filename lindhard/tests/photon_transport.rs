//! Photon transport: attenuation against the tabulated cross sections, energy
//! balance, edge cases, voxel geometry.

use lindhard::geometry::{Boundary, Geometry, Stack, VoxelGrid};
use lindhard::material::Material;
use lindhard::photon::{PhotonSource, PhotonTable, PhotonTransport, PhotonTransportError};
use lindhard::tally::PhotonTally;

const N_A: f64 = 6.022_140_76e23;

fn element(z: u8, rho_kg_m3: f64) -> Material {
    Material::from_atom_fractions(&[(z, 1.0)], Some(rho_kg_m3)).unwrap()
}

/// EPDL2017 `MF=23 MT=501` total cross sections, barns per atom, read from
/// the source file by linear interpolation in its own (linear-linear) law at
/// the energies below (the sum of partials equals the total there; pair
/// production is zero below 1.022 MeV). These are cited facts from the
/// operator-approved source (`docs/data-provenance.md`), independent of the
/// committed (thinned) table the transport reads.
const EPDL_TOTAL_BARNS: [(u8, f64, f64); 6] = [
    (13, 6.0e4, 12.4705),
    (13, 6.62e5, 3.34126),
    (29, 6.0e4, 168.638),
    (29, 6.62e5, 7.65147),
    (8, 6.0e4, 5.06699),
    (1, 6.0e4, 0.545205),
];

#[test]
fn committed_table_reproduces_the_source_totals() {
    let t = PhotonTable::epdl2017();
    for (z, e, barns) in EPDL_TOTAL_BARNS {
        let got = t.atom(z).unwrap().total_barns(e);
        assert!(
            (got / barns - 1.0).abs() < 2e-3,
            "Z = {z}, E = {e}: {got} vs {barns}"
        );
    }
}

/// Narrow beam through a slab of `mu t = 1`: the uncollided fraction leaving
/// the back face is `exp(-mu t)`, with `mu` from the EPDL totals.
fn narrow_beam(z: u8, rho: f64, e: f64, expected_barns: f64, a_gmol: f64, seed: u64) {
    let table = PhotonTable::epdl2017();
    let mat = element(z, rho);
    let n = rho * 1.0e3 * N_A / a_gmol; // atoms per m^3
    let mu = n * expected_barns * 1.0e-28;
    let thickness = 1.0 / mu;
    let stack = Stack::new(vec![(mat, thickness)], None).unwrap();
    let tr = PhotonTransport::new(&stack, &table).unwrap();
    // The coefficient the transport uses agrees with the EPDL total.
    let mu_t = tr.attenuation_coefficient(0, e);
    assert!((mu_t / mu - 1.0).abs() < 3e-3, "mu {mu_t} vs {mu}");

    let n_photons = 400_000u64;
    let tally = tr
        .run(&PhotonSource::beam(e), seed, n_photons, 1024)
        .unwrap();
    let p_expect = (-mu_t * thickness).exp();
    let p = tally.back.n_uncollided as f64 / n_photons as f64;
    let sigma = (p_expect * (1.0 - p_expect) / n_photons as f64).sqrt();
    assert!(
        (p - p_expect).abs() < 4.0 * sigma,
        "Z = {z}, E = {e}: {p} vs {p_expect} (sigma {sigma})"
    );
    // The EPDL-total expectation (exp(-1)) within the same error plus the
    // table's own 0.3 % in mu.
    assert!((p - (-1.0f64).exp()).abs() < 4.0 * sigma + 0.004);
    let bal = tally.energy_balance();
    assert!(bal.relative_residual() < 1e-12, "{bal:?}");
    assert_eq!(
        tally.front.n + tally.back.n + tally.n_photoelectric + tally.n_cutoff,
        n_photons
    );
}

#[test]
fn aluminium_narrow_beam_60_kev() {
    narrow_beam(13, 2699.0, 6.0e4, 12.4705, 26.981_538, 1);
}

#[test]
fn aluminium_narrow_beam_662_kev() {
    narrow_beam(13, 2699.0, 6.62e5, 3.34126, 26.981_538, 2);
}

#[test]
fn copper_narrow_beam_60_kev() {
    narrow_beam(29, 8960.0, 6.0e4, 168.638, 63.546, 3);
}

#[test]
fn copper_narrow_beam_662_kev() {
    narrow_beam(29, 8960.0, 6.62e5, 7.65147, 63.546, 4);
}

#[test]
fn energy_balance_closes_in_a_thick_multilayer_target() {
    let table = PhotonTable::epdl2017();
    let stack = Stack::new(
        vec![(element(13, 2699.0), 2.0e-3), (element(29, 8960.0), 1.0e-3)],
        Some(element(13, 2699.0)),
    )
    .unwrap();
    let tr = PhotonTransport::new(&stack, &table).unwrap();
    let t = tr.run(&PhotonSource::beam(1.5e5), 9, 50_000, 512).unwrap();
    let bal = t.energy_balance();
    assert!(bal.relative_residual() < 1e-12, "{bal:?}");
    assert_eq!(t.n_regions(), 3);
    assert!(t.photoelectric_ev.iter().all(|&e| e > 0.0));
    assert!(t.compton_ev.iter().all(|&e| e > 0.0));
    assert_eq!(t.side.n, 0, "a stack is laterally infinite");
    assert!(t.front.n > 0, "some photons are backscattered");
    // Every history ends exactly once.
    let ends = t.front.n + t.back.n + t.side.n + t.n_photoelectric + t.n_cutoff + t.n_missed;
    assert_eq!(ends, t.n_histories);
    assert!(t.n_coherent > 0 && t.n_incoherent > 0);
}

#[test]
fn photon_below_the_table_is_refused_and_one_scattered_below_is_absorbed() {
    let table = PhotonTable::epdl2017();
    let stack = Stack::semi_infinite(element(13, 2699.0));
    let tr = PhotonTransport::new(&stack, &table).unwrap();
    for e in [999.0, 1.0e6 + 1.0, f64::NAN] {
        assert!(matches!(
            tr.run(&PhotonSource::beam(e), 1, 1, 1),
            Err(PhotonTransportError::EnergyOutOfRange { .. })
        ));
    }
    // 1.02 keV in aluminium: any Compton loss goes below 1 keV or the photon
    // is absorbed; either way the balance closes.
    let t = tr.run(&PhotonSource::beam(1.0e3), 5, 10_000, 100).unwrap();
    assert!(t.energy_balance().relative_residual() < 1e-12);
    assert_eq!(t.n_photoelectric + t.n_cutoff + t.front.n, 10_000);
}

#[test]
fn source_validation() {
    let table = PhotonTable::epdl2017();
    let stack = Stack::semi_infinite(element(13, 2699.0));
    let tr = PhotonTransport::new(&stack, &table).unwrap();
    let mut s = PhotonSource::beam(1e5);
    s.direction = [0.0; 3];
    assert!(matches!(
        tr.run(&s, 1, 1, 1),
        Err(PhotonTransportError::BadDirection)
    ));
    s.direction = [f64::NAN, 0.0, 1.0];
    assert!(matches!(
        tr.run(&s, 1, 1, 1),
        Err(PhotonTransportError::BadDirection)
    ));
    s.direction = [1.0, 0.0, 0.0];
    s.position = Some([0.0, f64::INFINITY, 0.0]);
    assert!(matches!(
        tr.run(&s, 1, 1, 1),
        Err(PhotonTransportError::BadPosition)
    ));
}

#[test]
fn photon_starting_outside_the_target_is_counted_as_missed() {
    let table = PhotonTable::epdl2017();
    let stack = Stack::new(vec![(element(13, 2699.0), 1e-3)], None).unwrap();
    let tr = PhotonTransport::new(&stack, &table).unwrap();
    let mut s = PhotonSource::beam(1e5);
    s.position = Some([-1.0e-3, 0.0, 0.0]);
    let t = tr.run(&s, 1, 10, 4).unwrap();
    assert_eq!((t.n_missed, t.missed_ev), (10, 1.0e6));
    assert_eq!(t.energy_balance().residual_ev, 0.0);
}

#[test]
fn zero_thickness_layer_is_refused_by_the_geometry() {
    assert!(Stack::new(vec![(element(13, 2699.0), 0.0)], None).is_err());
}

#[test]
fn photon_born_on_an_interface_goes_either_way() {
    let table = PhotonTable::epdl2017();
    let stack = Stack::new(
        vec![(element(13, 2699.0), 1.0e-3), (element(29, 8960.0), 1.0e-3)],
        None,
    )
    .unwrap();
    let tr = PhotonTransport::new(&stack, &table).unwrap();
    for dir in [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]] {
        let mut s = PhotonSource::beam(8.0e4);
        s.position = Some([1.0e-3, 0.0, 0.0]);
        s.direction = dir;
        let t = tr.run(&s, 3, 20_000, 256).unwrap();
        assert!(t.energy_balance().relative_residual() < 1e-12);
        assert_eq!(t.n_missed, 0);
        if dir[0] > 0.0 {
            // Entering copper from aluminium's back face: the first flight is in layer 1.
            assert!(t.photoelectric_ev[1] > 0.0);
        } else {
            assert!(t.photoelectric_ev[0] > 0.0);
        }
    }
}

#[test]
fn photons_leave_through_each_face_of_a_voxel_grid() {
    let table = PhotonTable::epdl2017();
    // A 1 mm aluminium cube grid of 4x4x4 voxels, vacuum on every face.
    let grid = VoxelGrid::new(
        vec![element(13, 2699.0)],
        [4, 4, 4],
        [2.5e-4; 3],
        vec![0u32; 64],
        [Boundary::Vacuum; 3],
    )
    .unwrap();
    let tr = PhotonTransport::new(&grid, &table).unwrap();
    let mut s = PhotonSource::beam(1.0e5);
    s.position = Some([0.0, 5.0e-4, 5.0e-4]); // centre of the front face
    s.direction = [1.0, 0.0, 0.0];
    let t = tr.run(&s, 17, 40_000, 512).unwrap();
    assert_eq!(t.n_regions(), 64);
    assert!(
        t.energy_balance().relative_residual() < 1e-12,
        "{:?}",
        t.energy_balance()
    );
    assert!(t.front.n > 0 && t.back.n > 0 && t.side.n > 0, "{t:?}");
    assert!(t.back.n_uncollided > 0);
    assert_eq!(
        t.front.n + t.back.n + t.side.n + t.n_photoelectric + t.n_cutoff,
        t.n_histories
    );
    // Voxels along the beam line see more deposit than far corners.
    let line = grid.locate([1.2e-4, 5.0e-4, 5.0e-4]).unwrap();
    let corner = grid.locate([9.0e-4, 9.0e-4, 9.0e-4]).unwrap();
    assert!(t.deposited_ev(line) > t.deposited_ev(corner));
}

#[test]
fn missing_element_is_reported() {
    let table = PhotonTable::epdl2017();
    // Build a table lacking Z = 13 by keeping only hydrogen.
    let h_only = PhotonTable::new("test", vec![table.atom(1).unwrap().clone()]).unwrap();
    let stack = Stack::semi_infinite(element(13, 2699.0));
    assert!(matches!(
        PhotonTransport::new(&stack, &h_only),
        Err(PhotonTransportError::MissingElement { z: 13 })
    ));
}

#[test]
fn tally_type_is_reexported() {
    let _ = PhotonTally::new(1);
}
