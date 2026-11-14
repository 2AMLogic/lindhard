//! Helpers shared by the integration tests. Cargo does not build this
//! directory as its own test target; each test binary pulls it in with
//! `mod common;`.

#![allow(dead_code)]

use std::sync::OnceLock;

use lindhard::electron::boundary::{BandModel, BandStructure};
use lindhard::electron::data::{CrossSectionTable, CrossSectionTableParts, SamplingAxis};
use lindhard::electron::transport::LayerTables;
use lindhard::geometry::Stack;
use lindhard::ion::bca::{Bca, BcaConfig, Beam, CrystalTarget, SummaryTally};
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;

/// One Si-Si ZBL table shared by every test in a binary (the angle depends
/// on the screening function only, so one table serves every pair). The grid
/// is deliberately coarser than `TableSpec::default()` for test speed.
pub(crate) fn zbl_si_table() -> &'static ScatteringTable {
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

// The Si and SiO2 targets below differ only in how the surface binding energy
// is set, which moves sputtering and any binding-dependent tally. Each test
// keeps the variant it always used so no recorded value changes. The 15 eV
// displacement energy and the 2 eV surface binding are test parameters, not
// data.

/// Si with `E_d` = 15 eV and no surface binding energy set.
pub(crate) fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m
}

/// Si with `E_d` = 15 eV and `E_s` = 2 eV set unconditionally.
pub(crate) fn si_with_surface_binding() -> Material {
    let mut m = si();
    m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    m
}

/// Si with `E_d` = 15 eV and `E_s` = 2 eV only where it is not already set.
pub(crate) fn si_surface_binding_if_unset() -> Material {
    let mut m = si();
    if m.surface_binding_energy_ev(14).is_err() {
        m.set_surface_binding_energy_ev(14, 2.0).unwrap();
    }
    m
}

/// Amorphous SiO2 at 2200 kg/m3 with `E_d` = 15 eV for both elements and
/// `E_s` = 2 eV where it is not already set.
pub(crate) fn sio2() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
    for z in [14, 8] {
        m.set_displacement_energy_ev(z, 15.0).unwrap();
        if m.surface_binding_energy_ev(z).is_err() {
            m.set_surface_binding_energy_ev(z, 2.0).unwrap();
        }
    }
    m
}

/// Bare Si (Z = 14) with no displacement or surface binding energy set, as
/// the electron tests use it: they never read either.
pub(crate) fn bare_si() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

/// Energy grid (eV) of the synthetic electron tables.
pub(crate) const GRID: [f64; 5] = [1.0, 10.0, 100.0, 1_000.0, 10_000.0];

/// Cumulative probabilities of the synthetic table rows.
pub(crate) const PROB: [f64; 3] = [0.0, 0.5, 1.0];

/// A synthetic cross-section table on the energy grid `grid`: constant
/// inverse mean free path `rate` (1/m) and `row(E, p)` as the `p`-quantile at
/// energy `E` over [`PROB`]. A zero `rate` gives empty rows (no interaction).
/// The numbers are made up to exercise the code, not material data.
pub(crate) fn synthetic_table_on(
    grid: &[f64],
    axis: SamplingAxis,
    rate: f64,
    row: impl Fn(f64, f64) -> f64,
) -> CrossSectionTable {
    let quantiles = grid
        .iter()
        .map(|&e| {
            if rate == 0.0 {
                vec![]
            } else {
                PROB.iter().map(|&p| row(e, p)).collect()
            }
        })
        .collect();
    CrossSectionTable::new(CrossSectionTableParts {
        model: "synthetic".into(),
        material: "synthetic".into(),
        provenance: "synthetic test table".into(),
        axis,
        energy_ev: grid.to_vec(),
        inverse_mfp_per_m: vec![rate; grid.len()],
        probability: PROB.to_vec(),
        quantiles,
    })
    .unwrap()
}

/// [`synthetic_table_on`] on the standard [`GRID`].
pub(crate) fn synthetic_table(
    axis: SamplingAxis,
    rate: f64,
    row: impl Fn(f64, f64) -> f64,
) -> CrossSectionTable {
    synthetic_table_on(&GRID, axis, rate, row)
}

/// No elastic scattering on the standard [`GRID`]: a zero rate, so the
/// channel never fires.
pub(crate) fn elastic_none() -> CrossSectionTable {
    synthetic_table(SamplingAxis::ElasticPolarAngle, 0.0, |_, _| 0.0)
}

/// Isotropic elastic scattering at mean free path `lambda_m` (m) on the
/// standard [`GRID`]: `cos θ = 1 - 2p`.
pub(crate) fn elastic_isotropic(lambda_m: f64) -> CrossSectionTable {
    synthetic_table(SamplingAxis::ElasticPolarAngle, 1.0 / lambda_m, |_, p| {
        (1.0 - 2.0 * p).acos()
    })
}

/// Inelastic losses uniform on `[0, frac E]` at mean free path `lambda_m`
/// (m), on the standard [`GRID`].
pub(crate) fn inelastic_frac(lambda_m: f64, frac: f64) -> CrossSectionTable {
    synthetic_table(
        SamplingAxis::InelasticEnergyLoss,
        1.0 / lambda_m,
        move |e, p| frac * e * p,
    )
}

/// The two channels of one layer.
pub(crate) fn layer_pair(elastic: CrossSectionTable, inelastic: CrossSectionTable) -> LayerTables {
    LayerTables { elastic, inelastic }
}

/// A layer with isotropic elastic scattering at mean free path `lambda_el`
/// (m) and inelastic losses uniform on `[0, frac E]` at mean free path
/// `lambda_inel` (m).
pub(crate) fn synthetic_scattering(lambda_el: f64, lambda_inel: f64, frac: f64) -> LayerTables {
    layer_pair(
        elastic_isotropic(lambda_el),
        inelastic_frac(lambda_inel, frac),
    )
}

/// A layer with no interaction at all: both channels at a zero rate.
pub(crate) fn empty_layer() -> LayerTables {
    layer_pair(
        elastic_none(),
        synthetic_table(SamplingAxis::InelasticEnergyLoss, 0.0, |_, _| 0.0),
    )
}

/// Synthetic metal band structure: `E_F = 5 eV`, `Φ = 4 eV`, so `U = 9 eV`.
/// Made-up numbers chosen to exercise the code, not material data; the caller
/// passes the provenance string that says so.
pub(crate) fn synthetic_metal(provenance: &str) -> BandStructure {
    BandStructure::new(
        BandModel::Metal {
            fermi_ev: 5.0,
            work_function_ev: 4.0,
        },
        provenance,
    )
    .unwrap()
}

/// Synthetic insulator band structure: `W_v = 6 eV`, `E_g = 3 eV`,
/// `χ = 1 eV`, so `E_F = 7.5 eV` and `U = 10 eV`. Made-up numbers chosen to
/// exercise the code, not material data.
pub(crate) fn synthetic_insulator(provenance: &str) -> BandStructure {
    BandStructure::new(
        BandModel::Insulator {
            valence_band_width_ev: 6.0,
            band_gap_ev: 3.0,
            affinity_ev: 1.0,
        },
        provenance,
    )
    .unwrap()
}

/// Run `f` inside a rayon pool of exactly `n` threads, for the
/// thread-count determinism tests.
pub(crate) fn with_threads<R: Send>(n: usize, f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .unwrap()
        .install(f)
}

/// Relative deviation `|a / b - 1|`.
pub(crate) fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

/// `true` when `a` and `b` agree to relative tolerance `rel`.
pub(crate) fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()).max(f64::MIN_POSITIVE)
}

/// Composite Simpson rule in `ln x` on `[lo, hi]` with `n` (even) intervals.
pub(crate) fn simpson_log(f: &impl Fn(f64) -> f64, lo: f64, hi: f64, n: usize) -> f64 {
    let (a, b) = (lo.ln(), hi.ln());
    let h = (b - a) / n as f64;
    let mut s = 0.0;
    for i in 0..=n {
        let u = a + h * i as f64;
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let x = u.exp().clamp(lo, hi);
        s += w * f(x) * x;
    }
    s * h / 3.0
}

// The crystal helpers below are shared by the `crystal_*` integration tests
// (issue #751). `crystal_si` builds only the plain static lattice: a test that
// needs thermal vibration or an explicit search radius chains onto the result
// rather than growing an option parameter here.

/// Nanometre, m.
const NM: f64 = 1e-9;

/// Beam of `z` at `e` eV with polar angle `tilt_deg` and azimuth
/// `twist_deg`, both in degrees.
pub(crate) fn beam(z: u8, e: f64, tilt_deg: f64, twist_deg: f64, count: u64) -> Beam {
    Beam {
        ion: Ion::new(z).unwrap(),
        energy_ev: e,
        polar_rad: tilt_deg.to_radians(),
        azimuth_rad: twist_deg.to_radians(),
        count,
    }
}

/// A silicon crystal cut on (100), with the in-plane reference direction
/// [010], matching [`beam`]`(.., tilt_deg, twist_deg, ..)`. The static
/// lattice only; chain `.with_thermal(..)` or set `p_max_m` on the result for
/// the variants the individual tests need.
pub(crate) fn crystal_si(tilt_deg: f64, twist_deg: f64) -> CrystalTarget {
    let lat = Lattice::silicon();
    let o = Orientation::new(
        &lat,
        [1, 0, 0],
        [0, 1, 0],
        tilt_deg.to_radians(),
        twist_deg.to_radians(),
        0.0,
    )
    .unwrap();
    CrystalTarget::new(lat, o)
}

/// A [`SummaryTally`] run of beam `b` into `st` at `seed`, with `crystal`
/// attached to its regions if given, `tweak` applied to the config and the
/// run pinned to `threads` rayon threads if given.
pub(crate) fn run_summary(
    st: &Stack,
    b: Beam,
    seed: u64,
    crystal: Option<(CrystalTarget, Vec<usize>)>,
    tweak: impl Fn(&mut BcaConfig),
    threads: Option<usize>,
) -> SummaryTally {
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = seed;
    tweak(&mut cfg);
    let mut bca = Bca::new(b, st, cfg, &ls, zbl_si_table()).unwrap();
    if let Some((c, regions)) = crystal {
        bca = bca.with_crystal(c, &regions).unwrap();
    }
    let go = || bca.run(|| SummaryTally::new(0.5 * NM, 2000)).unwrap();
    match threads {
        None => go(),
        Some(n) => with_threads(n, go),
    }
}

/// Depth below which a fraction `q` of the stopped primaries lie, nm.
pub(crate) fn percentile_nm(t: &SummaryTally, q: f64) -> f64 {
    let tot: u64 = t.depth_hist.iter().sum();
    let mut acc = 0u64;
    for (i, &h) in t.depth_hist.iter().enumerate() {
        acc += h;
        if acc as f64 >= q * tot as f64 {
            return (i as f64 + 0.5) * t.bin_width_m / NM;
        }
    }
    f64::NAN
}

/// Fraction of stopped primaries deeper than `depth_m`.
pub(crate) fn fraction_deeper(t: &SummaryTally, depth_m: f64) -> f64 {
    let tot: u64 = t.depth_hist.iter().sum();
    let first = (depth_m / t.bin_width_m).ceil() as usize;
    t.depth_hist[first.min(t.depth_hist.len())..]
        .iter()
        .sum::<u64>() as f64
        / tot as f64
}
