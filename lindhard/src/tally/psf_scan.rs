//! Position-dependent point-spread function: the energy an electron beam
//! deposits in a thin depth slab as its landing point is scanned across a
//! laterally varying (voxel) target.
//!
//! # What is reported
//!
//! The radial profile of [`super::psf`] is the PSF of a beam on a laterally
//! uniform stack. Over a target with a buried edge or line array the energy
//! left in a film depends on *where* the beam lands, so this module runs the
//! electron transport once per landing point `(y, z)` and reports, for the
//! depth slab `depth_lo <= x < depth_hi`, the energy deposited per incident
//! electron, with its standard error. Each run is a separate
//! [`Transport::run_beam`] (parallel over histories, bit-identical at any
//! thread count); a single run with a wide beam and spatial binning would mix
//! the landing positions and lose the per-position response, so it is not
//! used.
//!
//! The slab energy and its error are those of the
//! [`RadialAccumulator`](super::psf::RadialAccumulator)
//! total (`RadialProfile::total_ev` and `total_std_err_ev`): per-history
//! slab totals `x_h` (the primary and all its secondaries), standard error
//! `(N/(N-1) (Q - S²/N))^(1/2)` of the sum `S = Σ x_h`, `Q = Σ x_h²`, see
//! [`super::psf`]. The energy per electron is `S / N` with error `SE / N`.
//! The radial bins of that accumulator are about the origin of the target
//! frame, not the landing point, and are not reported.
//!
//! # Determinism
//!
//! The landing points are looped sequentially in the order given. The seed of
//! point `i` is derived from the base seed and `i` alone ([`point_seed`]), so
//! it does not depend on the loop order or on the other points, and the report
//! is bit-identical at any thread count.
//!
//! The mixing function is SplitMix64 (G. L. Steele Jr., D. Lea and
//! C. H. Flood, "Fast splittable pseudorandom number generators", OOPSLA 2014,
//! doi:10.1145/2660193.2660195, Fig. 15 / the `mix64` finaliser).

use super::electron::{ElectronTallyConfig, ElectronTallyError, FullElectronTally};
use super::hist::Binning;
use super::psf::{LogRadialBinning, PsfConfig, PsfError};
use crate::electron::beam::BeamSpread;
use crate::electron::transport::{Primary, Transport, TransportError};
use crate::geometry::Geometry;
use serde::{Deserialize, Serialize};

/// Errors of a landing-point scan.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum PsfScanError {
    /// No landing point was given.
    #[error("a PSF scan needs at least one landing point")]
    NoPoints,
    /// A landing coordinate is not finite.
    #[error("landing point {index} is not finite: ({y}, {z}) m")]
    NonFinitePoint {
        /// Index of the point.
        index: usize,
        /// Its `y`, m.
        y: f64,
        /// Its `z`, m.
        z: f64,
    },
    /// Zero histories per point or a zero chunk size.
    #[error("histories per point and chunk size must be positive")]
    NoHistories,
    /// The depth slab is not finite with `lo < hi`.
    #[error(transparent)]
    Slab(#[from] PsfError),
    /// The tally could not be built.
    #[error(transparent)]
    Tally(#[from] ElectronTallyError),
    /// The transport rejected the primary or the beam.
    #[error(transparent)]
    Transport(#[from] TransportError),
}

/// What to scan: the depth slab, the landing points and the statistics.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PsfScanConfig {
    /// Lower depth of the slab, m.
    pub depth_lo_m: f64,
    /// Upper depth of the slab (excluded), m.
    pub depth_hi_m: f64,
    /// Landing points `(y, z)` on the front face, m, in report order. As for
    /// [`BeamSpread::landing_m`] they are offsets from the geometry's entry
    /// point (the centre of the front face of a voxel grid).
    pub landing_m: Vec<[f64; 2]>,
    /// Histories run at each point.
    pub histories_per_point: u64,
    /// Histories per tally chunk (see [`Transport::run_beam`]).
    pub chunk_size: u64,
    /// Base seed; the seed of point `i` is [`point_seed`]`(seed, i)`.
    pub seed: u64,
}

/// The seed of landing point `index` from the base `seed`: SplitMix64 of
/// `seed + (index + 1) · 0x9E3779B97F4A7C15`. A function of the two numbers
/// only (see the module docs for the source).
pub fn point_seed(seed: u64, index: u64) -> u64 {
    let mut z = seed.wrapping_add(index.wrapping_add(1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// The slab response at one landing point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PsfScanPoint {
    /// Landing point `(y, z)`, m, as an offset from the geometry's entry point.
    pub landing_m: [f64; 2],
    /// Seed of this point's run.
    pub seed: u64,
    /// Histories run.
    pub histories: u64,
    /// Energy deposited in the slab per incident electron, eV.
    pub energy_per_electron_ev: f64,
    /// Its standard error, eV (0 for fewer than two histories).
    pub std_err_ev: f64,
}

/// The energy deposited in a depth slab as a function of landing point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PsfScanReport {
    /// Free-text provenance of the run (what target, beam and tables).
    pub source: String,
    /// Primary kinetic energy, eV.
    pub primary_energy_ev: f64,
    /// Lower depth of the slab, m.
    pub depth_lo_m: f64,
    /// Upper depth of the slab (excluded), m.
    pub depth_hi_m: f64,
    /// Base seed.
    pub seed: u64,
    /// One entry per landing point, in the order requested.
    pub points: Vec<PsfScanPoint>,
}

impl PsfScanReport {
    /// The report as CSV, one row per landing point. Columns: `y_m`, `z_m`,
    /// `seed`, `histories`, `energy_per_electron_ev`, `std_err_ev`. Numbers
    /// are written in Rust's shortest round-trip form.
    pub fn to_csv(&self) -> String {
        let mut s = String::from("y_m,z_m,seed,histories,energy_per_electron_ev,std_err_ev\n");
        for p in &self.points {
            s.push_str(&format!(
                "{},{},{},{},{},{}\n",
                p.landing_m[0],
                p.landing_m[1],
                p.seed,
                p.histories,
                p.energy_per_electron_ev,
                p.std_err_ev
            ));
        }
        s
    }
}

/// Scan the landing point of `beam` over `scan.landing_m`: one
/// [`Transport::run_beam`] per point, the point replacing `beam.landing_m`
/// (the beam's profile and divergence are kept). `source` is recorded in the
/// report.
pub fn scan_landing_points<G: Geometry>(
    transport: &Transport<G>,
    primary: &Primary,
    beam: &BeamSpread,
    scan: &PsfScanConfig,
    source: &str,
) -> Result<PsfScanReport, PsfScanError> {
    if scan.landing_m.is_empty() {
        return Err(PsfScanError::NoPoints);
    }
    if scan.histories_per_point == 0 || scan.chunk_size == 0 {
        return Err(PsfScanError::NoHistories);
    }
    if let Some((index, p)) = scan
        .landing_m
        .iter()
        .enumerate()
        .find(|(_, p)| !(p[0].is_finite() && p[1].is_finite()))
    {
        return Err(PsfScanError::NonFinitePoint {
            index,
            y: p[0],
            z: p[1],
        });
    }
    // One wide radial bin: only the slab total (with its per-history error)
    // is used; the bins about the frame origin are not reported.
    let psf = PsfConfig::new(
        LogRadialBinning::new(1.0, 2.0, 1)?,
        scan.depth_lo_m,
        scan.depth_hi_m,
    )?;
    let mut tally_config = ElectronTallyConfig::new(
        Binning::new(0.0, primary.energy_ev, 1).expect("one bin on a positive energy"),
        Binning::new(0.0, std::f64::consts::FRAC_PI_2, 1).expect("one polar bin"),
    );
    tally_config.psf = Some(psf);
    let proto = FullElectronTally::new(transport, tally_config)?;

    let mut points = Vec::with_capacity(scan.landing_m.len());
    for (i, &landing) in scan.landing_m.iter().enumerate() {
        let seed = point_seed(scan.seed, i as u64);
        let spread = BeamSpread {
            landing_m: landing,
            ..*beam
        };
        let run = transport.run_beam(
            seed,
            scan.histories_per_point,
            scan.chunk_size,
            primary,
            &spread,
            || proto.clone(),
        )?;
        let profile = run
            .tally
            .report()
            .deposition
            .psf
            .expect("the PSF accumulator was configured");
        let n = profile.histories as f64;
        points.push(PsfScanPoint {
            landing_m: landing,
            seed,
            histories: profile.histories,
            energy_per_electron_ev: profile.total_ev / n,
            std_err_ev: profile.total_std_err_ev / n,
        });
    }
    Ok(PsfScanReport {
        source: source.to_string(),
        primary_energy_ev: primary.energy_ev,
        depth_lo_m: scan.depth_lo_m,
        depth_hi_m: scan.depth_hi_m,
        seed: scan.seed,
        points,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_seed_depends_on_both_inputs_and_is_stable() {
        assert_ne!(point_seed(1, 0), point_seed(1, 1));
        assert_ne!(point_seed(1, 0), point_seed(2, 0));
        assert_eq!(point_seed(7, 3), point_seed(7, 3));
        // SplitMix64 reference stream from state 0: first output.
        assert_eq!(point_seed(0, 0), 0xE220_A839_7B1D_CDAF);
    }

    #[test]
    fn csv_has_one_row_per_point() {
        let r = PsfScanReport {
            source: "t".into(),
            primary_energy_ev: 1.0,
            depth_lo_m: 0.0,
            depth_hi_m: 1.0,
            seed: 1,
            points: vec![
                PsfScanPoint {
                    landing_m: [1.0, 2.0],
                    seed: 3,
                    histories: 4,
                    energy_per_electron_ev: 5.5,
                    std_err_ev: 0.5,
                };
                2
            ],
        };
        let csv = r.to_csv();
        assert_eq!(csv.lines().count(), 3);
        assert!(csv.starts_with("y_m,z_m,seed,histories,"));
        assert_eq!(csv.lines().nth(1).unwrap(), "1,2,3,4,5.5,0.5");
    }
}
