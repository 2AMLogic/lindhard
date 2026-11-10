//! The committed EPDL2017 photon table: coverage, sanity and (with a local
//! copy of the source file) agreement with the source.

use lindhard::photon::data::{ENERGY_MAX_EV, ENERGY_MIN_EV, X_MAX, X_MIN};
use lindhard::photon::PhotonTable;

#[test]
fn epdl2017_table_covers_z_1_to_92() {
    let t = PhotonTable::epdl2017();
    assert!(t.provenance().contains("EPDL2017"));
    assert!(t.provenance().contains("IAEA-NDS-225"));
    for z in 1..=92u8 {
        let a = t.atom(z).unwrap_or_else(|| panic!("Z = {z} missing"));
        for c in [&a.photoelectric, &a.incoherent, &a.coherent] {
            assert_eq!(c.x()[0], ENERGY_MIN_EV);
            assert_eq!(*c.x().last().unwrap(), ENERGY_MAX_EV);
        }
        for c in [&a.form_factor, &a.scattering_function] {
            assert_eq!(c.x()[0], X_MIN);
            assert_eq!(*c.x().last().unwrap(), X_MAX);
        }
        // The form factor is Z at x -> 0 and the scattering function tends to Z.
        assert!(
            (a.form_factor.eval(X_MIN) / f64::from(z) - 1.0).abs() < 0.02,
            "Z = {z}"
        );
        assert!(
            (a.scattering_function.eval(X_MAX) / f64::from(z) - 1.0).abs() < 1e-3,
            "Z = {z}"
        );
    }
}

#[test]
fn hydrogen_form_factor_is_the_hydrogenic_one() {
    // F(x) = 1 / (1 + (2 pi a0 x)^2)^2, x = sin(theta/2)/lambda (a0 in Angstrom).
    let h = PhotonTable::epdl2017();
    let h = h.atom(1).unwrap();
    for x in [0.1, 0.3, 1.0, 3.0] {
        let exact = 1.0 / (1.0 + (2.0 * std::f64::consts::PI * 0.529_177 * x).powi(2)).powi(2);
        assert!(
            (h.form_factor.eval(x) / exact - 1.0).abs() < 2e-3,
            "x = {x}"
        );
    }
}

#[test]
fn cross_sections_are_physical() {
    let t = PhotonTable::epdl2017();
    // Compton per atom grows with Z roughly as Z (binding matters at 1 keV only).
    let at = |z: u8, e: f64| t.atom(z).unwrap().incoherent.eval(e);
    assert!((at(8, 1.0e6) / at(1, 1.0e6) / 8.0 - 1.0).abs() < 0.02);
    // Photoelectric falls steeply with energy away from edges; lead K edge at 88 keV.
    let pe = |e: f64| t.atom(82).unwrap().photoelectric.eval(e);
    assert!(pe(88.1e3) / pe(87.9e3) > 4.0, "no K edge jump");
    assert!(pe(200e3) < pe(100e3) / 4.0);
}

/// Needs a local `EPDL2017.ALL` (path in `LINDHARD_EPDL2017`): the committed
/// table equals what the reader makes of the file, and every source point in
/// the window is reproduced within the thinning tolerance.
#[test]
fn epdl2017_committed_table_matches_local_source_file() {
    let Ok(path) = std::env::var("LINDHARD_EPDL2017") else {
        eprintln!("LINDHARD_EPDL2017 not set; skipping");
        return;
    };
    let committed = PhotonTable::epdl2017();
    let fresh = PhotonTable::from_epdl_endf6_file(&path, committed.provenance()).unwrap();
    assert_eq!(committed, fresh);
}
