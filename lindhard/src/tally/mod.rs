//! Tallies: turning particle histories into results.
//!
//! * [`hist`]: fixed-binning histograms with integer counts.
//! * [`moments`]: mergeable one-pass moments (mean, `σ`, `γ`, `β`) with
//!   standard errors.
//! * [`pearson`]: Pearson IV and dual-Pearson representations of a range
//!   profile.
//! * [`ion`]: the [`crate::ion::bca::BcaTally`] that collects range, damage
//!   (NRT model estimates alongside cascade defect counts) and
//!   backscatter/transmission/sputter statistics, and its plain-data
//!   [`IonReport`].
//! * [`electron`]: the [`crate::electron::transport::ElectronTally`] that
//!   collects 3D energy deposition (Cartesian and cylindrical grids),
//!   backscatter and secondary yields split at a configurable energy
//!   (50 eV by default), escape spectra, generation-volume moments and the
//!   energy balance, and its plain-data [`ElectronReport`].
//! * [`psf`]: the radial point-spread function of a pencil beam in a depth
//!   slab (log radial bins, per-history errors) and its double- and
//!   triple-Gaussian fits, with CSV export and serde derives.
//! * [`table_coverage`]: counts of the electron transport's cross-section
//!   rate evaluations inside and beyond each table's energy grid, per layer
//!   and channel (a numerical diagnostic, carried by the electron report).
//!
//! Every accumulator here merges deterministically: integer counts are exact
//! sums, and floating-point accumulators are combined in the order the caller
//! merges them, which [`crate::rng::run_particles`] fixes to chunk order. The
//! report types derive `serde::Serialize`, with units documented per field
//! (SI lengths in metres, energies in eV, angles in radians).

pub mod electron;
pub mod hist;
pub mod ion;
pub mod moments;
pub mod pearson;
pub mod psf;
pub mod table_coverage;

pub use electron::{
    CartesianDeposition, CartesianGrid, CylindricalDeposition, CylindricalGrid, DepositionReport,
    ElectronEnergyBudget, ElectronReport, ElectronTallyConfig, ElectronTallyError,
    ElectronTallyMetadata, EmissionClass, EventCapCounts, FaceEmission, FateCounts,
    FullElectronTally, GenerationVolume, PrimaryStoppingPoints, StoppingPoints, Yields,
    SE_BSE_SPLIT_EV,
};
pub use hist::{Binning, BinningError, Histogram};
pub use ion::{
    CascadeDefects, DamageReport, EscapeReport, FaceEscape, IonReport, IonTally, IonTallyConfig,
    LayerDamage, LayerRange, NrtDamage, RangeReport, SpeciesEscape,
};
pub use moments::{MomentSummary, Moments};
pub use pearson::{type_iv_min_kurtosis, DualPearson, DualPearsonFit, PearsonError, PearsonIv};
pub use psf::{
    fit_psf, GaussianPsf, LogRadialBinning, PsfConfig, PsfError, PsfFit, PsfFitOptions, PsfModel,
    PsfNormalization, PsfReport, PsfResidual, RadialAccumulator, RadialProfile,
};
pub use table_coverage::{LayerTableCoverage, TableCoverageCounts, TableCoverageTally};
