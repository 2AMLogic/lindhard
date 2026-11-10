//! Lateral implanted-species and displacement-damage profiles of a 2-D
//! dynamic run, the film cut-through result and the comparison of a milled
//! profile with published points.
//!
//! # What is tallied
//!
//! [`DynamicRun2d::with_profiles`](super::DynamicRun2d::with_profiles) bins the
//! events of every primary on a [`VoxelTallyGrid`] in the sample frame
//! (`x` the depth from the original front surface, `y` the lateral
//! coordinate), reusing [`VoxelTally`] for the event definitions:
//!
//! * `implanted`: beam particles that came to rest in the voxel;
//! * `displacements`: target atoms displaced in the voxel (recoils created
//!   above the displacement energy), the simulated cascade count. These are
//!   the definitions of the 1-D [`IonReport`](crate::tally::ion::IonReport)
//!   and of the [`VoxelTally`]; the NRT estimate is not tallied. Amorphization
//!   and damage accumulation are not modelled (the count is per ion, summed
//!   over the run).
//!
//! The counts of a step are scaled by the fluence each primary stands for
//! (`F p_j / count` ions/m² for a primary entering column `j`) and by the
//! column width `h`, so a voxel holds **atoms per metre of the invariant
//! axis `z`** (`*_per_m`). Dividing by the voxel cross-section `dx dy` gives
//! the number density (`*_per_m3`). Counts of the substrate are also kept
//! apart (`substrate_displacements_per_m`): they are the events whose region
//! was the substrate layer at the time, which is the damage left in the layer
//! beneath the film. Rejected adaptive attempts are not accumulated, and the
//! sums are made in column order, so profiles are bit-identical at any thread
//! count.
//!
//! # Export
//!
//! [`ProfileSnapshot::write_csv`] follows the conventions of
//! [`VoxelReport::write_csv`](crate::tally::voxel::VoxelReport::write_csv): a
//! `# key = value` header carrying the grid, a column row, then one row per
//! occupied voxel in flat-index order. The snapshot also implements `serde::Serialize`
//! (JSON with any serde format crate).

use std::collections::BTreeMap;
use std::io::{self, Write};

use serde::{Deserialize, Serialize};

use crate::tally::voxel::{VoxelReport, VoxelTallyGrid};

/// Accumulated tallies of one voxel, in atoms per metre of `z`.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct ProfileCell {
    /// Beam atoms at rest here.
    pub implanted_per_m: f64,
    /// Target atoms displaced here (all layers).
    pub displacements_per_m: f64,
    /// Of those, displaced while in the substrate layer.
    pub substrate_displacements_per_m: f64,
}

impl ProfileCell {
    fn add(&mut self, o: &ProfileCell) {
        self.implanted_per_m += o.implanted_per_m;
        self.displacements_per_m += o.displacements_per_m;
        self.substrate_displacements_per_m += o.substrate_displacements_per_m;
    }
}

/// Running profile sums of a run (or of one step before it is accepted).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ProfileAccumulator {
    grid: VoxelTallyGrid,
    cells: BTreeMap<u64, ProfileCell>,
    outside: ProfileCell,
}

impl ProfileAccumulator {
    pub(crate) fn new(grid: VoxelTallyGrid) -> Self {
        Self {
            grid,
            cells: BTreeMap::new(),
            outside: ProfileCell::default(),
        }
    }

    /// Add the counts of all events (`all`) and of the substrate events
    /// (`substrate`) of one column, each count standing for `scale` atoms/m.
    pub(crate) fn add_counts(&mut self, all: &VoxelReport, substrate: &VoxelReport, scale: f64) {
        for v in &all.voxels {
            let c = self.cells.entry(self.grid.flat_index(v.index)).or_default();
            c.implanted_per_m += v.cell.implanted as f64 * scale;
            c.displacements_per_m += v.cell.displacements as f64 * scale;
        }
        self.outside.implanted_per_m += all.outside.implanted as f64 * scale;
        self.outside.displacements_per_m += all.outside.displacements as f64 * scale;
        for v in &substrate.voxels {
            let c = self.cells.entry(self.grid.flat_index(v.index)).or_default();
            c.substrate_displacements_per_m += v.cell.displacements as f64 * scale;
        }
        self.outside.substrate_displacements_per_m +=
            substrate.outside.displacements as f64 * scale;
    }

    pub(crate) fn merge(&mut self, o: &ProfileAccumulator) {
        for (&k, c) in &o.cells {
            self.cells.entry(k).or_default().add(c);
        }
        self.outside.add(&o.outside);
    }

    pub(crate) fn snapshot(&self, beam_z: u8, fluence_m2: f64, ions: u64) -> ProfileSnapshot {
        ProfileSnapshot {
            beam_z,
            peak_fluence_m2: fluence_m2,
            ions_per_column: ions,
            grid: self.grid,
            voxels: self
                .cells
                .iter()
                .map(|(&k, &cell)| (self.grid.unflatten(k), cell))
                .collect(),
            outside: self.outside,
        }
    }
}

/// The profiles at one fluence. See the module docs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProfileSnapshot {
    /// Atomic number of the beam species.
    pub beam_z: u8,
    /// Peak fluence delivered when the snapshot was taken, ions/m².
    pub peak_fluence_m2: f64,
    /// Primaries per illuminated column used.
    pub ions_per_column: u64,
    /// The grid the voxel indices refer to.
    pub grid: VoxelTallyGrid,
    /// Occupied voxels, `[i, j, k]` with their tallies, in flat-index order.
    pub voxels: Vec<([usize; 3], ProfileCell)>,
    /// Events outside the grid.
    pub outside: ProfileCell,
}

impl ProfileSnapshot {
    /// Sum over all voxels and the outside.
    pub fn grand_total(&self) -> ProfileCell {
        let mut t = self.outside;
        for (_, c) in &self.voxels {
            t.add(c);
        }
        t
    }

    /// Voxel cross-section `dx dy` (m²), the divisor from atoms/m to
    /// atoms/m³.
    pub fn voxel_area_m2(&self) -> f64 {
        let s = self.grid.spacing_m();
        s[0] * s[1]
    }

    /// Laterally integrated profile along the depth axis: entry `i` sums the
    /// voxels with `x` index `i`.
    pub fn depth_profile(&self) -> Vec<ProfileCell> {
        let mut out = vec![ProfileCell::default(); self.grid.dims()[0]];
        for (ijk, c) in &self.voxels {
            out[ijk[0]].add(c);
        }
        out
    }

    /// Profile along the lateral axis: entry `j` sums the voxels with `y`
    /// index `j`.
    pub fn lateral_profile(&self) -> Vec<ProfileCell> {
        let mut out = vec![ProfileCell::default(); self.grid.dims()[1]];
        for (ijk, c) in &self.voxels {
            out[ijk[1]].add(c);
        }
        out
    }

    /// Write the CSV export (see the module docs); densities are in
    /// atoms/m³.
    pub fn write_csv<W: Write>(&self, mut w: W) -> io::Result<()> {
        let g = &self.grid;
        let j3 = |a: [f64; 3]| format!("{:e},{:e},{:e}", a[0], a[1], a[2]);
        writeln!(w, "# lindhard 2d profile v1")?;
        writeln!(w, "# beam_z = {}", self.beam_z)?;
        writeln!(w, "# peak_fluence_m2 = {:e}", self.peak_fluence_m2)?;
        writeln!(w, "# ions_per_column = {}", self.ions_per_column)?;
        writeln!(w, "# origin_m = {}", j3(g.origin_m()))?;
        writeln!(w, "# spacing_m = {}", j3(g.spacing_m()))?;
        let d = g.dims();
        writeln!(w, "# dims = {},{},{}", d[0], d[1], d[2])?;
        let o = &self.outside;
        writeln!(
            w,
            "# outside_per_m = {:e},{:e},{:e}",
            o.implanted_per_m, o.displacements_per_m, o.substrate_displacements_per_m
        )?;
        writeln!(
            w,
            "i,j,k,x_m,y_m,z_m,implanted_per_m,displacements_per_m,\
             substrate_displacements_per_m,implanted_per_m3,displacements_per_m3,\
             substrate_displacements_per_m3"
        )?;
        let a = self.voxel_area_m2();
        for (ijk, c) in &self.voxels {
            let m = g.center_m(*ijk);
            writeln!(
                w,
                "{},{},{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e}",
                ijk[0],
                ijk[1],
                ijk[2],
                m[0],
                m[1],
                m[2],
                c.implanted_per_m,
                c.displacements_per_m,
                c.substrate_displacements_per_m,
                c.implanted_per_m / a,
                c.displacements_per_m / a,
                c.substrate_displacements_per_m / a
            )?;
        }
        Ok(())
    }

    /// The CSV export as a string.
    pub fn to_csv(&self) -> String {
        let mut buf = Vec::new();
        self.write_csv(&mut buf).expect("writing to a Vec");
        String::from_utf8(buf).expect("the CSV is ASCII")
    }
}

/// Result of [`DynamicRun2d::run_until_cut_through`](super::DynamicRun2d::run_until_cut_through).
#[derive(Debug, Clone, PartialEq)]
pub struct CutThrough {
    /// Whether an illuminated column became bare substrate within the run.
    pub reached: bool,
    /// The first column that did (the lowest index if several in one step).
    pub column: Option<usize>,
    /// Peak fluence at the end of the last step before the cut-through,
    /// ions/m². The true cut-through fluence lies in
    /// `(fluence_before_m2, fluence_m2]`; the bracket is one step wide, which
    /// the adaptive policy bounds by the allowed composition change.
    pub fluence_before_m2: f64,
    /// Peak fluence at the end of the step in which the cut-through was seen
    /// (the fluence delivered if it was not reached), ions/m².
    pub fluence_m2: f64,
    /// Primaries per column used.
    pub ions_done: u64,
    /// The profiles at the stop, if profiles are on. The damage in the layer
    /// beneath is in `substrate_displacements_per_m`.
    pub profiles: Option<ProfileSnapshot>,
}

/// One published point compared with the simulated profile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointDifference {
    /// Lateral position of the published point (relative to the beam centre),
    /// m.
    pub y_m: f64,
    /// Published depth, m.
    pub published_m: f64,
    /// Simulated depth at that position (linear interpolation), m.
    pub simulated_m: f64,
}

/// A milled profile compared with published points. See [`compare_profile`].
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileComparison {
    /// Per-point values, in the order given.
    pub points: Vec<PointDifference>,
    /// Root-mean-square depth difference, m.
    pub rms_m: f64,
    /// Largest absolute depth difference, m.
    pub max_abs_m: f64,
    /// `rms_m` over the largest published depth.
    pub relative_rms: f64,
    /// The tolerance applied to `relative_rms`.
    pub tolerance: f64,
    /// `relative_rms <= tolerance`.
    pub within_tolerance: bool,
}

/// Compare a simulated profile with published points.
///
/// `sim_y_m` and `sim_depth_m` are the lateral coordinate of the column
/// centres (relative to the beam centre) and the depth of the surface there
/// (increasing `y`); `published` holds `(y, depth)` pairs, m, in the same
/// frame. The simulation is interpolated linearly; a published point outside
/// the simulated range is an error, never extrapolated. The metric is the
/// RMS depth difference over the published depth maximum; the tolerance is
/// the caller's, with nothing fitted: a gap is reported, not removed.
pub fn compare_profile(
    sim_y_m: &[f64],
    sim_depth_m: &[f64],
    published: &[(f64, f64)],
    tolerance: f64,
) -> Result<ProfileComparison, String> {
    if sim_y_m.len() != sim_depth_m.len() || sim_y_m.len() < 2 {
        return Err("the simulated profile needs at least two matching points".into());
    }
    if published.is_empty() {
        return Err("no published points".into());
    }
    if !(tolerance.is_finite() && tolerance > 0.0) {
        return Err("the tolerance must be finite and positive".into());
    }
    if sim_y_m.windows(2).any(|w| w[1] <= w[0]) {
        return Err("the simulated lateral coordinates must increase".into());
    }
    let mut points = Vec::with_capacity(published.len());
    for &(y, d) in published {
        if !(y.is_finite() && d.is_finite()) {
            return Err("published points must be finite".into());
        }
        let (lo, hi) = (sim_y_m[0], sim_y_m[sim_y_m.len() - 1]);
        if y < lo || y > hi {
            return Err(format!(
                "published point y = {y:e} m is outside [{lo:e}, {hi:e}]"
            ));
        }
        let i = sim_y_m
            .partition_point(|&v| v <= y)
            .clamp(1, sim_y_m.len() - 1);
        let t = (y - sim_y_m[i - 1]) / (sim_y_m[i] - sim_y_m[i - 1]);
        let s = sim_depth_m[i - 1] + t * (sim_depth_m[i] - sim_depth_m[i - 1]);
        points.push(PointDifference {
            y_m: y,
            published_m: d,
            simulated_m: s,
        });
    }
    let n = points.len() as f64;
    let rms = (points
        .iter()
        .map(|p| (p.simulated_m - p.published_m).powi(2))
        .sum::<f64>()
        / n)
        .sqrt();
    let max_abs = points
        .iter()
        .map(|p| (p.simulated_m - p.published_m).abs())
        .fold(0.0, f64::max);
    let scale = points
        .iter()
        .map(|p| p.published_m.abs())
        .fold(0.0, f64::max);
    if scale == 0.0 {
        return Err("the largest published depth is zero".into());
    }
    let relative_rms = rms / scale;
    Ok(ProfileComparison {
        points,
        rms_m: rms,
        max_abs_m: max_abs,
        relative_rms,
        tolerance,
        within_tolerance: relative_rms <= tolerance,
    })
}
