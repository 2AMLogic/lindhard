//! Golden results: small fixed-seed runs whose outputs are pinned to values
//! stored in this file, so that the CI OS matrix (x86-64 Linux, aarch64 Linux,
//! aarch64 macOS) compares each platform with a committed result instead of
//! only with itself. The thread-count tests (`determinism.rs`, `tally.rs`,
//! ...) never leave one process.
//!
//! What is exact where:
//!
//! - Amorphous BCA and electron runs: integer outputs (counts, histogram
//!   bins) are asserted exactly on every platform; floating point outputs to a
//!   relative tolerance of [`REL_TOL`], because `f64::exp`, `ln`, `sin`, `cos`,
//!   `powf`, `tan` and `atan2` defer to the platform libm and may differ in the
//!   last bits between platforms.
//! - Crystal BCA run: three tiers, chosen by platform.
//!   - Linux (x86-64 and aarch64, both verified in CI): exact (integers) and
//!     [`REL_TOL`] (floats) against [`EXPECTED_CRYSTAL`].
//!   - aarch64 macOS: exact in the same sense against the separately pinned
//!     [`EXPECTED_CRYSTAL_AARCH64_MACOS`].
//!   - Every other OS or architecture: only the statistical tolerances in
//!     [`check_crystal_statistical`] (conservation invariants exact, counts
//!     and sums to a stated relative tolerance, no per-bin histogram check).
//!
//!   The crystal set is platform-labelled because the platform libm differs in
//!   the last bit (glibc against macOS libm, for `tan`, `acos`, `sin`, `cos`,
//!   `exp`, `ln`, `powf` and others). The differences enter through the BCA
//!   path shared with the amorphous run, where they only move floats in the
//!   last bits and happen to leave the integer counts alone; the crystal
//!   path's lattice-site and simultaneous-collision comparisons amplify them
//!   into different histories (`recoils` = 67721 on aarch64 macOS against
//!   67786 on Linux, PR #216 CI). There is no single flipped comparison. The
//!   engine keeps its tolerance-only cross-platform promise (`docs/architecture.md`,
//!   "Reproducibility"); a pure-Rust libm that would make the two sets
//!   identical was considered and declined (#217). Pinning a second set keeps
//!   exact regression detection on the macOS runner. If the `macos-latest`
//!   image's libm changes, that set may need regenerating. Set
//!   `LINDHARD_GOLDEN_CRYSTAL_STATISTICAL=1` to force the statistical check on
//!   any platform.
//!
//! Regenerating the expected values (after an intended change to the physics,
//! the sampling order or the RNG streams):
//!
//! ```text
//! LINDHARD_GOLDEN_PRINT=1 cargo test -p lindhard --test golden -- --nocapture --test-threads=1
//! ```
//!
//! prints each pinned quantity as `name = value`; copy the values into the
//! `EXPECTED_*` constants below and say in the commit why they changed. All
//! values below (amorphous, crystal, electron) were generated on x86-64 Linux
//! except [`EXPECTED_CRYSTAL_AARCH64_MACOS`], which was generated on aarch64
//! macOS; `EXPECTED_CRYSTAL` is the exact reference for Linux only, and the
//! other platforms check the amorphous and electron sets exactly. Regenerate
//! the Linux sets on x86-64 Linux and the macOS crystal set separately, on
//! aarch64 macOS (or from the `macos-latest` CI log with the command above),
//! and state the platform in the commit.
//!
//! The electron run uses constant mean free paths and isotropic angles
//! computed in this file (not measured or tabulated data).

use std::sync::OnceLock;

use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::{LayerTables, Primary, Transport, TransportConfig};
use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, Beam, CrystalTarget, SummaryTally};
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{Binning, ElectronTallyConfig, FullElectronTally};

/// Relative tolerance for floating-point outputs.
const REL_TOL: f64 = 1e-9;
const NM: f64 = 1e-9;

fn printing() -> bool {
    std::env::var_os("LINDHARD_GOLDEN_PRINT").is_some()
}

/// Compare (or, when regenerating, print) one pinned quantity.
struct Golden {
    name: &'static str,
}

impl Golden {
    fn int(&self, what: &str, actual: u64, expected: u64) {
        if printing() {
            println!("{}.{what} = {actual}", self.name);
        } else {
            assert_eq!(actual, expected, "{}.{what}", self.name);
        }
    }

    fn ints(&self, what: &str, actual: &[u64], expected: &[u64]) {
        if printing() {
            println!("{}.{what} = {actual:?}", self.name);
        } else {
            assert_eq!(actual, expected, "{}.{what}", self.name);
        }
    }

    fn float(&self, what: &str, actual: f64, expected: f64) {
        if printing() {
            println!("{}.{what} = {actual:e}", self.name);
        } else {
            let tol = REL_TOL * expected.abs().max(actual.abs());
            assert!(
                (actual - expected).abs() <= tol,
                "{}.{what}: {actual:e} vs expected {expected:e} (relative tolerance {REL_TOL:e})",
                self.name
            );
        }
    }
}

/// Sum a histogram into `n` equal groups of bins so that it fits in a line.
fn fold(h: &[u64], n: usize) -> Vec<u64> {
    let per = h.len().div_ceil(n);
    h.chunks(per).map(|c| c.iter().sum()).collect()
}

// ---------------------------------------------------------------------------
// Ion BCA: amorphous and crystal
// ---------------------------------------------------------------------------

fn scattering_table() -> &'static ScatteringTable {
    static T: OnceLock<ScatteringTable> = OnceLock::new();
    T.get_or_init(|| {
        ScatteringTable::build(
            &Potential::new(Screening::ZblUniversal, 14.0, 14.0),
            &TableSpec {
                eps_min: 1e-6,
                eps_max: 1e4,
                beta_min: 1e-5,
                beta_max: 1e2,
                per_decade: 16,
            },
        )
    })
}

fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m
}

fn ion_run(crystal: bool) -> SummaryTally {
    let st = Stack::semi_infinite(si());
    let (tilt, twist): (f64, f64) = (7.0, 22.0);
    let beam = Beam {
        ion: Ion::new(33).unwrap(),
        energy_ev: 2.0e4,
        polar_rad: tilt.to_radians(),
        azimuth_rad: twist.to_radians(),
        count: 100,
    };
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 0x601D;
    let mut bca = Bca::new(beam, &st, cfg, &ls, scattering_table()).unwrap();
    if crystal {
        let lat = Lattice::silicon();
        let o = Orientation::new(
            &lat,
            [1, 0, 0],
            [0, 1, 0],
            tilt.to_radians(),
            twist.to_radians(),
            0.0,
        )
        .unwrap();
        bca = bca.with_crystal(CrystalTarget::new(lat, o), &[0]).unwrap();
    }
    bca.run(|| SummaryTally::new(0.5 * NM, 2000)).unwrap()
}

/// The pinned crystal set for this platform, or `None` where only the
/// statistical check applies (other platforms, or when forced off by
/// `LINDHARD_GOLDEN_CRYSTAL_STATISTICAL`).
fn crystal_exact() -> Option<&'static IonExpected> {
    if std::env::var_os("LINDHARD_GOLDEN_CRYSTAL_STATISTICAL").is_some() {
        None
    } else if cfg!(target_os = "linux") {
        Some(&EXPECTED_CRYSTAL)
    } else if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
        Some(&EXPECTED_CRYSTAL_AARCH64_MACOS)
    } else {
        None
    }
}

/// Relative tolerances for the crystal run on platforms where it is not exact.
/// Chosen to be well above the run-to-run spread of 100 histories (about 9%
/// Poisson on the sputter count, a few percent on cascade sizes) and well
/// below a real physics change; the observed macOS recoil count differs by
/// 0.1%.
const CRYSTAL_COUNT_TOL: f64 = 0.15;
const CRYSTAL_SUM_TOL: f64 = 0.15;

/// Documented relaxation for the crystal run on platforms without a pinned
/// set; see the module docs.
fn check_crystal_statistical(g: &Golden, t: &SummaryTally, e: &IonExpected) {
    let near = |what: &str, actual: f64, expected: f64, tol: f64| {
        if printing() {
            println!("{}.{what} = {actual:e}", g.name);
        } else {
            assert!(
                (actual - expected).abs() <= tol * expected.abs(),
                "{}.{what}: {actual:e} vs expected {expected:e} (statistical tolerance {tol})",
                g.name
            );
        }
    };
    // Conservation invariants stay exact.
    g.int("histories", t.histories, 100);
    g.int("primaries_stopped", t.primaries_stopped, e.stopped);
    g.int("backscattered", t.backscattered, e.backscattered);
    g.int("transmitted", t.transmitted, e.transmitted);
    let binned: u64 = t.depth_hist.iter().sum();
    g.int("depth_hist_total", binned, 100);
    near(
        "sputtered",
        t.sputtered as f64,
        e.sputtered as f64,
        CRYSTAL_COUNT_TOL,
    );
    near(
        "recoils",
        t.recoils as f64,
        e.recoils as f64,
        CRYSTAL_COUNT_TOL,
    );
    near("depth_sum_m", t.depth_sum, e.depth_sum, CRYSTAL_SUM_TOL);
    near(
        "depth_sq_sum_m2",
        t.depth_sq_sum,
        e.depth_sq_sum,
        2.0 * CRYSTAL_SUM_TOL,
    );
    near("lattice_ev", t.budget.lattice, e.lattice, CRYSTAL_SUM_TOL);
    near(
        "electronic_nonlocal_ev",
        t.budget.electronic_nonlocal,
        e.electronic_nonlocal,
        CRYSTAL_SUM_TOL,
    );
    assert!(t.max_relative_residual < 1e-9);
}

fn check_ion(g: &Golden, t: &SummaryTally, e: &IonExpected) {
    g.int("histories", t.histories, 100);
    g.int("primaries_stopped", t.primaries_stopped, e.stopped);
    g.int("backscattered", t.backscattered, e.backscattered);
    g.int("transmitted", t.transmitted, e.transmitted);
    g.int("sputtered", t.sputtered, e.sputtered);
    g.int("recoils", t.recoils, e.recoils);
    // Depth bins are 0.5 nm; the first 100 bins (50 nm) are folded into 20
    // groups of 2.5 nm, the rest (deeper stoppers) into one count.
    g.ints(
        "depth_hist_folded20",
        &fold(&t.depth_hist[..100], 20),
        e.depth_hist,
    );
    g.int(
        "depth_hist_beyond_50nm",
        t.depth_hist[100..].iter().sum(),
        e.beyond,
    );
    g.float("depth_sum_m", t.depth_sum, e.depth_sum);
    g.float("depth_sq_sum_m2", t.depth_sq_sum, e.depth_sq_sum);
    g.float("lattice_ev", t.budget.lattice, e.lattice);
    g.float(
        "electronic_nonlocal_ev",
        t.budget.electronic_nonlocal,
        e.electronic_nonlocal,
    );
    assert!(t.max_relative_residual < 1e-9);
}

struct IonExpected {
    stopped: u64,
    backscattered: u64,
    transmitted: u64,
    sputtered: u64,
    recoils: u64,
    depth_hist: &'static [u64],
    beyond: u64,
    depth_sum: f64,
    depth_sq_sum: f64,
    lattice: f64,
    electronic_nonlocal: f64,
}

const EXPECTED_AMORPHOUS: IonExpected = IonExpected {
    stopped: 100,
    backscattered: 0,
    transmitted: 0,
    sputtered: 155,
    recoils: 58821,
    depth_hist: &[
        0, 1, 6, 3, 7, 11, 16, 12, 14, 12, 4, 8, 2, 3, 1, 0, 0, 0, 0, 0,
    ],
    beyond: 0,
    depth_sum: 1.895393904897188e-6,
    depth_sq_sum: 4.0803189719301385e-14,
    lattice: 1.1947263335161665e6,
    electronic_nonlocal: 7.308933894123819e5,
};

const EXPECTED_CRYSTAL: IonExpected = IonExpected {
    stopped: 100,
    backscattered: 0,
    transmitted: 0,
    sputtered: 133,
    recoils: 67786,
    depth_hist: &[
        0, 0, 1, 6, 12, 15, 13, 13, 11, 6, 6, 6, 4, 2, 2, 0, 1, 1, 0, 0,
    ],
    beyond: 1,
    depth_sum: 1.972785412561225e-6,
    depth_sq_sum: 4.644724413036624e-14,
    lattice: 1.2876489887391336e6,
    electronic_nonlocal: 6.361895180662028e5,
};

/// Crystal reference for aarch64 macOS, generated on the `macos-latest` runner
/// (image `macos-26-arm64` 20261013.0424, macOS 26.6.2) from the fixed-seed
/// run above at 9e9697c, 2026-10-15. Differs from [`EXPECTED_CRYSTAL`] only
/// because the platform libm differs in the last bit; see the module docs.
/// Regenerate on aarch64 macOS, separately from the Linux sets.
const EXPECTED_CRYSTAL_AARCH64_MACOS: IonExpected = IonExpected {
    stopped: 100,
    backscattered: 0,
    transmitted: 0,
    sputtered: 133,
    recoils: 67721,
    depth_hist: &[
        0, 0, 1, 6, 12, 15, 13, 13, 11, 6, 6, 6, 4, 2, 2, 0, 1, 1, 0, 0,
    ],
    beyond: 1,
    depth_sum: 1.9714783626537694e-6,
    depth_sq_sum: 4.6289472492559004e-14,
    lattice: 1.2878137943367276e6,
    electronic_nonlocal: 6.360087276137662e5,
};

#[test]
fn golden_amorphous_bca_batch() {
    let g = Golden { name: "amorphous" };
    check_ion(&g, &ion_run(false), &EXPECTED_AMORPHOUS);
}

#[test]
fn golden_crystal_bca_batch() {
    let g = Golden { name: "crystal" };
    let t = ion_run(true);
    if let Some(expected) = crystal_exact() {
        check_ion(&g, &t, expected);
    } else {
        check_crystal_statistical(&g, &t, &EXPECTED_CRYSTAL);
    }
}

// ---------------------------------------------------------------------------
// Electron transport
// ---------------------------------------------------------------------------

const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];
const PROB: [f64; 3] = [0.0, 0.5, 1.0];

fn et_table(axis: SamplingAxis, rate: f64, row: impl Fn(f64, f64) -> f64) -> CrossSectionTable {
    let quantiles = GRID
        .iter()
        .map(|&e| PROB.iter().map(|&p| row(e, p)).collect())
        .collect();
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "computed in tests/golden.rs".into(),
        axis,
        energy_ev: GRID.to_vec(),
        inverse_mfp_per_m: vec![rate; GRID.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

/// Isotropic elastic scattering (`cos θ = 1 - 2p`) and inelastic losses
/// uniform on `[0, 0.3 E]`.
fn layer(lambda_el: f64, lambda_inel: f64) -> LayerTables {
    LayerTables {
        elastic: et_table(SamplingAxis::ElasticPolarAngle, 1.0 / lambda_el, |_, p| {
            (1.0 - 2.0 * p).acos()
        }),
        inelastic: et_table(
            SamplingAxis::InelasticEnergyLoss,
            1.0 / lambda_inel,
            |e, p| 0.3 * e * p,
        ),
    }
}

struct ElectronExpected {
    stopped: u64,
    escaped_front: u64,
    escaped_back: u64,
    front_hist: &'static [u64],
    deposited_ev: f64,
    escaped_ev: f64,
    backscatter_eta: f64,
    secondary_delta: f64,
}

const EXPECTED_ELECTRON: ElectronExpected = ElectronExpected {
    stopped: 172,
    escaped_front: 808,
    escaped_back: 20,
    front_hist: &[131, 48, 57, 50, 60, 60, 67, 73, 67, 39],
    deposited_ev: 1.0959295086172502e6,
    escaped_ev: 9.040704913827501e5,
    backscatter_eta: 7.48e-1,
    secondary_delta: 6e-2,
};

#[test]
fn golden_electron_batch() {
    let g = Golden { name: "electron" };
    let e0 = 2_000.0;
    let mat = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    let stack = Stack::new(vec![(mat.clone(), 5.0 * NM), (mat, 30.0 * NM)], None).unwrap();
    let t = Transport::new(
        stack,
        vec![layer(3.0 * NM, 5.0 * NM), layer(2.0 * NM, 4.0 * NM)],
        TransportConfig::new(5.0),
    )
    .unwrap();
    let cfg = ElectronTallyConfig::new(
        Binning::new(0.0, e0, 50).unwrap(),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 9).unwrap(),
    );
    let proto = FullElectronTally::new(&t, cfg).unwrap();
    let r = t
        .run(0xE1EC, 1_000, 64, &Primary::normal(e0), || proto.clone())
        .unwrap()
        .tally
        .report();
    let e = &EXPECTED_ELECTRON;
    g.int("histories", r.histories, 1_000);
    g.int("stopped", r.fates.stopped, e.stopped);
    g.int("escaped_front", r.fates.escaped_front, e.escaped_front);
    g.int("escaped_back", r.fates.escaped_back, e.escaped_back);
    g.ints(
        "front_energy_hist_folded10",
        &fold(&r.front.energy_histogram.counts, 10),
        e.front_hist,
    );
    g.float("deposited_ev", r.budget.deposited_ev, e.deposited_ev);
    g.float("escaped_ev", r.budget.escaped_ev, e.escaped_ev);
    g.float(
        "backscatter_eta",
        r.yields.backscatter_eta,
        e.backscatter_eta,
    );
    g.float(
        "secondary_delta",
        r.yields.secondary_delta,
        e.secondary_delta,
    );
}
