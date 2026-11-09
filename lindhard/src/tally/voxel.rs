//! Sparse 3D voxel tallies of implanted atoms, damage and deposited energy,
//! with CSV and VTK export.
//!
//! [`VoxelTally`] implements [`BcaTally`]. It bins events by position on its
//! own regular grid ([`VoxelTallyGrid`]), so it works with any
//! [`crate::geometry::Geometry`] (a [`crate::geometry::VoxelGrid`], a mesh or
//! a CSG solid), not only with a voxel grid; use
//! [`VoxelTallyGrid::from_voxel_grid`] to tally on the voxels of a grid
//! target. Build one prototype with [`VoxelTally::new`], hand clones to
//! [`crate::ion::bca::Bca::run`], and turn the merged tally into a plain-data
//! [`VoxelReport`] with [`VoxelTally::report`].
//!
//! Units: positions and lengths in metres, energies in eV. `x` is the depth.
//!
//! # What is tallied per voxel
//!
//! * `implanted`: beam particles that came to rest in the voxel (the voxel
//!   counterpart of [`crate::tally::ion::RangeReport::stopped`]).
//! * `displacements`: target atoms displaced in the voxel, counted where the
//!   recoil is created (transfer above `E_d`), the simulated cascade count of
//!   [`crate::tally::ion::CascadeDefects::displacements`]. This is the
//!   damage measure of the tally. The NRT / Kinchin-Pease model estimate
//!   ([`crate::ion::damage`]) needs a per-material Lindhard partition and is
//!   not tallied per voxel; the two are never conflated (see
//!   [`crate::ion::damage`] for the convention and its sources, Norgett,
//!   Robinson and Torrens, Nucl. Eng. Des. 33 (1975) 50, and Lindhard et al.,
//!   Mat. Fys. Medd. Dan. Vid. Selsk. 33 (1963) no. 10).
//! * Deposited energy, split by channel so each can be checked against the
//!   engine's [`crate::ion::bca::EnergyBudget`]:
//!   `electronic_ev` (nonlocal and local electronic loss),
//!   `lattice_ev` (nuclear energy left in the lattice: subthreshold and weak
//!   transfers and the binding energy of displaced atoms) and `rest_ev` (the
//!   kinetic energy of a particle when it stops below its cutoff).
//!   [`VoxelCell::deposited_ev`] is their sum. The work done against the
//!   planar surface barrier and the energy carried out of the target are not
//!   deposited in any voxel and are not tallied here.
//!
//! The nonlocal electronic loss of a free flight is spread uniformly along
//! the straight segment (the engine's convention, see
//! [`crate::ion::bca::BcaTally::electronic`]), so a segment that crosses
//! voxel faces is split between the voxels in proportion to the length in
//! each. Everything else is a point event and goes to the voxel containing
//! it. Voxels own the half-open box `[o + i h, o + (i + 1) h)` on each axis;
//! events outside the grid are summed in [`VoxelReport::outside`].
//!
//! # Sparse storage
//!
//! Only voxels that received an event are stored, in a `BTreeMap` keyed on the
//! flat index `i + nx (j + ny k)`, so empty voxels cost no memory and
//! iteration, merging and export follow one fixed order.
//!
//! # Determinism
//!
//! Counts are integers and energies are `f64` sums, merged voxel by voxel in
//! the order [`crate::rng::run_particles`] calls `merge` (chunk order). For a
//! fixed `chunk_size` the result is bit-identical at any thread count
//! (`lindhard/tests/voxel_tally.rs`).
//!
//! # Export formats
//!
//! [`VoxelReport::write_csv`] writes one row per occupied voxel after a
//! `# key = value` header that carries the grid, so the file reads back
//! ([`VoxelReport::read_csv`]) exactly: floats are written in the shortest
//! form that parses back to the same bits. [`VoxelReport::write_vtk`] writes
//! a legacy ASCII VTK 3.0 `UNSTRUCTURED_GRID` of `VTK_VOXEL` cells, one per
//! occupied voxel (so the file is as sparse as the tally), with the tallies
//! as cell scalars; format: Kitware, *The VTK User's Guide*, 11th ed. (2006),
//! "VTK File Formats". The VTK file does not carry the `outside` sums.
//! [`VoxelReport::read_vtk`] reads back what `write_vtk` writes.

use std::collections::BTreeMap;
use std::io::{self, Write};

use serde::{Deserialize, Serialize};

use crate::geometry::VoxelGrid;
use crate::ion::bca::{BcaTally, ElectronicChannel, Face, LatticeDeposit, Particle};

/// Errors from building a [`VoxelTallyGrid`] or reading an exported file.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum VoxelTallyError {
    /// A grid dimension was zero.
    #[error("voxel tally grid dimension on axis {axis} is zero")]
    ZeroDimension {
        /// Axis (0 = x, 1 = y, 2 = z).
        axis: usize,
    },
    /// A spacing was zero, negative or not finite.
    #[error(
        "voxel tally spacing {spacing_m} m on axis {axis} is invalid: must be finite and positive"
    )]
    InvalidSpacing {
        /// Axis (0 = x, 1 = y, 2 = z).
        axis: usize,
        /// Offending value.
        spacing_m: f64,
    },
    /// An origin coordinate was not finite.
    #[error("voxel tally origin {origin_m} m on axis {axis} is not finite")]
    InvalidOrigin {
        /// Axis (0 = x, 1 = y, 2 = z).
        axis: usize,
        /// Offending value.
        origin_m: f64,
    },
    /// The number of voxels does not fit in `u64`.
    #[error("voxel tally dimensions {dims:?} overflow the voxel count")]
    TooManyVoxels {
        /// Requested dimensions.
        dims: [usize; 3],
    },
    /// An exported file could not be read. `line` is 1-based.
    #[error("voxel tally parse error (line {line}): {message}")]
    Parse {
        /// Line number in the file.
        line: usize,
        /// What was wrong.
        message: String,
    },
}

fn parse_err(line: usize, message: impl Into<String>) -> VoxelTallyError {
    VoxelTallyError::Parse {
        line,
        message: message.into(),
    }
}

/// The regular grid a [`VoxelTally`] bins on: `dims` voxels of edge
/// `spacing_m` starting at `origin_m`. Voxel `(i, j, k)` has the flat index
/// `i + nx (j + ny k)` (x fastest, as in [`VoxelGrid`]).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VoxelTallyGrid {
    origin_m: [f64; 3],
    spacing_m: [f64; 3],
    dims: [usize; 3],
}

impl VoxelTallyGrid {
    /// A grid of `dims` voxels of edge lengths `spacing_m` whose lowest corner
    /// is `origin_m`. Every dimension must be nonzero, every spacing finite
    /// and positive, the origin finite, and the voxel count must fit in `u64`.
    pub fn new(
        origin_m: [f64; 3],
        spacing_m: [f64; 3],
        dims: [usize; 3],
    ) -> Result<Self, VoxelTallyError> {
        for axis in 0..3 {
            if dims[axis] == 0 {
                return Err(VoxelTallyError::ZeroDimension { axis });
            }
            let h = spacing_m[axis];
            if !(h.is_finite() && h > 0.0) {
                return Err(VoxelTallyError::InvalidSpacing { axis, spacing_m: h });
            }
            if !origin_m[axis].is_finite() {
                return Err(VoxelTallyError::InvalidOrigin {
                    axis,
                    origin_m: origin_m[axis],
                });
            }
        }
        (dims[0] as u64)
            .checked_mul(dims[1] as u64)
            .and_then(|v| v.checked_mul(dims[2] as u64))
            .ok_or(VoxelTallyError::TooManyVoxels { dims })?;
        Ok(Self {
            origin_m,
            spacing_m,
            dims,
        })
    }

    /// The grid of the voxels of `grid` (origin at the corner of the front
    /// face), so a tally voxel is a target voxel.
    pub fn from_voxel_grid(grid: &VoxelGrid) -> Self {
        Self {
            origin_m: [0.0; 3],
            spacing_m: grid.spacing_m(),
            dims: grid.dims(),
        }
    }

    /// Lowest corner `[x, y, z]`, m.
    pub fn origin_m(&self) -> [f64; 3] {
        self.origin_m
    }

    /// Voxel edge lengths `[hx, hy, hz]`, m.
    pub fn spacing_m(&self) -> [f64; 3] {
        self.spacing_m
    }

    /// Voxels per axis `[nx, ny, nz]`.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Total number of voxels `nx ny nz`.
    pub fn n_voxels(&self) -> u64 {
        self.dims.iter().map(|&n| n as u64).product()
    }

    /// The voxel `[i, j, k]` containing `pos`, or `None` outside the grid
    /// (or for a non-finite position). Voxels are half-open, so a point on an
    /// interior face belongs to the voxel with the higher index and a point
    /// on the far end of the grid is outside.
    pub fn locate(&self, pos: [f64; 3]) -> Option<[usize; 3]> {
        let mut ijk = [0usize; 3];
        for a in 0..3 {
            let s = ((pos[a] - self.origin_m[a]) / self.spacing_m[a]).floor();
            if !(s >= 0.0 && s < self.dims[a] as f64) {
                return None;
            }
            ijk[a] = s as usize;
        }
        Some(ijk)
    }

    /// Flat index `i + nx (j + ny k)` of voxel `ijk`.
    pub fn flat_index(&self, ijk: [usize; 3]) -> u64 {
        ijk[0] as u64 + self.dims[0] as u64 * (ijk[1] as u64 + self.dims[1] as u64 * ijk[2] as u64)
    }

    /// The voxel `[i, j, k]` of flat index `flat` (the inverse of
    /// [`Self::flat_index`]).
    pub fn unflatten(&self, flat: u64) -> [usize; 3] {
        let nx = self.dims[0] as u64;
        let ny = self.dims[1] as u64;
        [
            (flat % nx) as usize,
            ((flat / nx) % ny) as usize,
            (flat / (nx * ny)) as usize,
        ]
    }

    /// Lowest corner of voxel `ijk`, m.
    pub fn corner_m(&self, ijk: [usize; 3]) -> [f64; 3] {
        std::array::from_fn(|a| self.origin_m[a] + ijk[a] as f64 * self.spacing_m[a])
    }

    /// Centre of voxel `ijk`, m.
    pub fn center_m(&self, ijk: [usize; 3]) -> [f64; 3] {
        std::array::from_fn(|a| self.origin_m[a] + (ijk[a] as f64 + 0.5) * self.spacing_m[a])
    }
}

/// The tallies of one voxel (or of the whole grid, for totals).
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct VoxelCell {
    /// Beam particles that came to rest here.
    pub implanted: u64,
    /// Target atoms displaced here (recoils created).
    pub displacements: u64,
    /// Electronic energy loss deposited here, eV.
    pub electronic_ev: f64,
    /// Nuclear energy left in the lattice here, eV.
    pub lattice_ev: f64,
    /// Kinetic energy of particles that stopped here, eV.
    pub rest_ev: f64,
}

impl VoxelCell {
    /// Energy deposited in the target: `electronic_ev + lattice_ev +
    /// rest_ev`, eV.
    pub fn deposited_ev(&self) -> f64 {
        self.electronic_ev + self.lattice_ev + self.rest_ev
    }

    /// Field-wise sum.
    pub fn add(&mut self, o: &VoxelCell) {
        self.implanted += o.implanted;
        self.displacements += o.displacements;
        self.electronic_ev += o.electronic_ev;
        self.lattice_ev += o.lattice_ev;
        self.rest_ev += o.rest_ev;
    }
}

/// A sparse voxel tally; see the [module docs](self).
#[derive(Debug, Clone)]
pub struct VoxelTally {
    grid: VoxelTallyGrid,
    histories: u64,
    cells: BTreeMap<u64, VoxelCell>,
    outside: VoxelCell,
    /// Scratch for splitting a flight segment (not merged).
    cuts: Vec<f64>,
}

impl VoxelTally {
    /// An empty tally on `grid`.
    pub fn new(grid: VoxelTallyGrid) -> Self {
        Self {
            grid,
            histories: 0,
            cells: BTreeMap::new(),
            outside: VoxelCell::default(),
            cuts: Vec::new(),
        }
    }

    /// The grid this tally bins on.
    pub fn grid(&self) -> &VoxelTallyGrid {
        &self.grid
    }

    /// Number of stored (occupied) voxels.
    pub fn occupied(&self) -> usize {
        self.cells.len()
    }

    fn cell_at(&mut self, pos: [f64; 3]) -> &mut VoxelCell {
        match self.grid.locate(pos) {
            Some(ijk) => self.cells.entry(self.grid.flat_index(ijk)).or_default(),
            None => &mut self.outside,
        }
    }

    /// Spread `energy_ev` uniformly along the straight segment `from..to`,
    /// adding each voxel's share with `add`.
    fn spread(
        &mut self,
        from: [f64; 3],
        to: [f64; 3],
        energy_ev: f64,
        add: fn(&mut VoxelCell, f64),
    ) {
        let (a, b) = (self.grid.locate(from), self.grid.locate(to));
        if a.is_some() && a == b {
            add(self.cell_at(from), energy_ev);
            return;
        }
        // Fractions of the segment at which it crosses a voxel face, with
        // faces clamped to the grid (beyond it everything is "outside").
        let mut cuts = std::mem::take(&mut self.cuts);
        cuts.clear();
        cuts.push(0.0);
        for ax in 0..3 {
            let h = self.grid.spacing_m[ax];
            let f0 = (from[ax] - self.grid.origin_m[ax]) / h;
            let f1 = (to[ax] - self.grid.origin_m[ax]) / h;
            if !(f0.is_finite() && f1.is_finite()) || f0 == f1 {
                continue;
            }
            let n = self.grid.dims[ax] as f64;
            let first = (f0.min(f1).floor() + 1.0).max(0.0);
            let last = (f0.max(f1).ceil() - 1.0).min(n);
            let mut k = first;
            while k <= last {
                let t = (k - f0) / (f1 - f0);
                if t > 0.0 && t < 1.0 {
                    cuts.push(t);
                }
                k += 1.0;
            }
        }
        cuts.push(1.0);
        cuts.sort_by(f64::total_cmp);
        for w in cuts.windows(2) {
            let (t0, t1) = (w[0], w[1]);
            if t1 > t0 {
                let tm = 0.5 * (t0 + t1);
                let mid = std::array::from_fn(|i| from[i] + tm * (to[i] - from[i]));
                add(self.cell_at(mid), energy_ev * (t1 - t0));
            }
        }
        self.cuts = cuts;
    }

    /// The merged tally as plain data: occupied voxels in flat-index order.
    pub fn report(&self) -> VoxelReport {
        VoxelReport {
            grid: self.grid,
            histories: self.histories,
            voxels: self
                .cells
                .iter()
                .map(|(&flat, &cell)| VoxelRecord {
                    index: self.grid.unflatten(flat),
                    cell,
                })
                .collect(),
            outside: self.outside,
        }
    }
}

impl BcaTally for VoxelTally {
    fn electronic(
        &mut self,
        p: &Particle,
        from: [f64; 3],
        channel: ElectronicChannel,
        energy_ev: f64,
    ) {
        let add = |c: &mut VoxelCell, e: f64| c.electronic_ev += e;
        match channel {
            ElectronicChannel::NonLocal => self.spread(from, p.pos, energy_ev, add),
            ElectronicChannel::Local => add(self.cell_at(p.pos), energy_ev),
        }
    }

    fn lattice(&mut self, at: [f64; 3], _layer: usize, _kind: LatticeDeposit, energy_ev: f64) {
        self.cell_at(at).lattice_ev += energy_ev;
    }

    fn recoil(&mut self, r: &Particle) {
        self.cell_at(r.pos).displacements += 1;
    }

    fn stopped(&mut self, p: &Particle) {
        let primary = p.is_primary();
        let cell = self.cell_at(p.pos);
        cell.rest_ev += p.energy_ev;
        if primary {
            cell.implanted += 1;
        }
    }

    fn escaped(&mut self, _p: &Particle, _face: Face) {}

    fn end_history(&mut self, _index: u64, _budget: &crate::ion::bca::EnergyBudget) {
        self.histories += 1;
    }

    /// Voxel-by-voxel merge.
    ///
    /// # Panics
    /// If `other` was built on a different grid.
    fn merge(&mut self, o: Self) {
        assert!(
            self.grid == o.grid,
            "merging voxel tallies built on different grids"
        );
        self.histories += o.histories;
        self.outside.add(&o.outside);
        for (flat, cell) in o.cells {
            self.cells.entry(flat).or_default().add(&cell);
        }
    }
}

/// One occupied voxel of a [`VoxelReport`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VoxelRecord {
    /// Voxel indices `[i, j, k]`.
    pub index: [usize; 3],
    /// Its tallies.
    pub cell: VoxelCell,
}

/// A [`VoxelTally`] as plain data. Voxels are in increasing flat-index order
/// (`i + nx (j + ny k)`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VoxelReport {
    /// The grid the voxels refer to.
    pub grid: VoxelTallyGrid,
    /// Primary histories tallied.
    pub histories: u64,
    /// Occupied voxels.
    pub voxels: Vec<VoxelRecord>,
    /// Sum of the events outside the grid.
    pub outside: VoxelCell,
}

const CSV_COLUMNS: [&str; 12] = [
    "i",
    "j",
    "k",
    "x_m",
    "y_m",
    "z_m",
    "implanted",
    "displacements",
    "electronic_ev",
    "lattice_ev",
    "rest_ev",
    "deposited_ev",
];

const VTK_SCALARS: [(&str, &str); 6] = [
    ("implanted", "unsigned_long"),
    ("displacements", "unsigned_long"),
    ("electronic_ev", "double"),
    ("lattice_ev", "double"),
    ("rest_ev", "double"),
    ("deposited_ev", "double"),
];

/// Corner offsets of a `VTK_VOXEL` cell, in VTK point order (x fastest).
const VOXEL_CORNERS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

const VTK_VOXEL: u8 = 11;

fn join3(v: [f64; 3]) -> String {
    format!("{:e},{:e},{:e}", v[0], v[1], v[2])
}

fn parse3<T: std::str::FromStr>(s: &str, line: usize) -> Result<[T; 3], VoxelTallyError> {
    let v: Vec<T> = s
        .split(',')
        .map(|t| t.trim().parse::<T>())
        .collect::<Result<_, _>>()
        .map_err(|_| parse_err(line, format!("cannot parse the triple {s:?}")))?;
    v.try_into()
        .map_err(|_| parse_err(line, format!("expected three values in {s:?}")))
}

fn num<T: std::str::FromStr>(s: &str, line: usize) -> Result<T, VoxelTallyError> {
    s.trim()
        .parse::<T>()
        .map_err(|_| parse_err(line, format!("cannot parse the number {s:?}")))
}

/// Whitespace-separated tokens of a VTK file with their 1-based lines.
struct Tokens<'a> {
    iter: Box<dyn Iterator<Item = (usize, &'a str)> + 'a>,
    last_line: usize,
}

impl<'a> Tokens<'a> {
    fn next(&mut self, what: &str) -> Result<(usize, &'a str), VoxelTallyError> {
        self.iter
            .next()
            .ok_or_else(|| parse_err(self.last_line, format!("expected {what}")))
    }

    fn expect(&mut self, word: &str) -> Result<(), VoxelTallyError> {
        let (n, t) = self.next(word)?;
        if t == word {
            Ok(())
        } else {
            Err(parse_err(n, format!("expected `{word}`, found `{t}`")))
        }
    }
}

impl VoxelReport {
    /// Sum over the occupied voxels (not including [`Self::outside`]), added
    /// in flat-index order.
    pub fn totals(&self) -> VoxelCell {
        let mut t = VoxelCell::default();
        for v in &self.voxels {
            t.add(&v.cell);
        }
        t
    }

    /// Sum over the grid and the outside: everything the tally received.
    pub fn grand_totals(&self) -> VoxelCell {
        let mut t = self.totals();
        t.add(&self.outside);
        t
    }

    /// Per-slab sums along `axis` (0 = x): entry `n` adds every occupied
    /// voxel whose index on that axis is `n`. This is the laterally
    /// integrated profile, comparable to a 1D depth tally of the same bins.
    pub fn profile(&self, axis: usize) -> Vec<VoxelCell> {
        let mut out = vec![VoxelCell::default(); self.grid.dims[axis]];
        for v in &self.voxels {
            out[v.index[axis]].add(&v.cell);
        }
        out
    }

    /// Write the CSV export: a `# key = value` header, a column row, then one
    /// row per occupied voxel (indices, centre in m, tallies; energies in
    /// eV). See the [module docs](self).
    pub fn write_csv<W: Write>(&self, mut w: W) -> io::Result<()> {
        let g = &self.grid;
        let o = &self.outside;
        writeln!(w, "# lindhard voxel tally v1")?;
        writeln!(w, "# origin_m = {}", join3(g.origin_m))?;
        writeln!(w, "# spacing_m = {}", join3(g.spacing_m))?;
        writeln!(w, "# dims = {},{},{}", g.dims[0], g.dims[1], g.dims[2])?;
        writeln!(w, "# histories = {}", self.histories)?;
        writeln!(
            w,
            "# outside = {},{},{:e},{:e},{:e}",
            o.implanted, o.displacements, o.electronic_ev, o.lattice_ev, o.rest_ev
        )?;
        writeln!(w, "{}", CSV_COLUMNS.join(","))?;
        for v in &self.voxels {
            let c = g.center_m(v.index);
            let t = &v.cell;
            writeln!(
                w,
                "{},{},{},{:e},{:e},{:e},{},{},{:e},{:e},{:e},{:e}",
                v.index[0],
                v.index[1],
                v.index[2],
                c[0],
                c[1],
                c[2],
                t.implanted,
                t.displacements,
                t.electronic_ev,
                t.lattice_ev,
                t.rest_ev,
                t.deposited_ev()
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

    /// Read back a file written by [`Self::write_csv`].
    pub fn read_csv(text: &str) -> Result<Self, VoxelTallyError> {
        let mut origin = None;
        let mut spacing = None;
        let mut dims = None;
        let mut histories = None;
        let mut outside = None;
        let mut lines = text.lines().enumerate().map(|(i, l)| (i + 1, l.trim()));
        let mut header_row = None;
        for (n, line) in lines.by_ref() {
            if line.is_empty() {
                continue;
            }
            let Some(meta) = line.strip_prefix('#') else {
                header_row = Some((n, line));
                break;
            };
            let Some((key, value)) = meta.split_once('=') else {
                continue; // the format banner
            };
            match key.trim() {
                "origin_m" => origin = Some(parse3::<f64>(value, n)?),
                "spacing_m" => spacing = Some(parse3::<f64>(value, n)?),
                "dims" => dims = Some(parse3::<usize>(value, n)?),
                "histories" => histories = Some(num::<u64>(value, n)?),
                "outside" => {
                    let f: Vec<&str> = value.split(',').collect();
                    if f.len() != 5 {
                        return Err(parse_err(n, "outside needs five values"));
                    }
                    outside = Some(VoxelCell {
                        implanted: num(f[0], n)?,
                        displacements: num(f[1], n)?,
                        electronic_ev: num(f[2], n)?,
                        lattice_ev: num(f[3], n)?,
                        rest_ev: num(f[4], n)?,
                    });
                }
                other => return Err(parse_err(n, format!("unknown header key {other:?}"))),
            }
        }
        let missing = |what: &str| parse_err(1, format!("header is missing `{what}`"));
        let grid = VoxelTallyGrid::new(
            origin.ok_or_else(|| missing("origin_m"))?,
            spacing.ok_or_else(|| missing("spacing_m"))?,
            dims.ok_or_else(|| missing("dims"))?,
        )?;
        let (hn, header_row) = header_row.ok_or_else(|| parse_err(1, "no column row"))?;
        if header_row != CSV_COLUMNS.join(",") {
            return Err(parse_err(hn, "unexpected column row"));
        }
        let mut voxels: Vec<VoxelRecord> = Vec::new();
        for (n, line) in lines {
            if line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split(',').collect();
            if f.len() != CSV_COLUMNS.len() {
                return Err(parse_err(n, "wrong number of columns"));
            }
            let index = [num(f[0], n)?, num(f[1], n)?, num(f[2], n)?];
            if (0..3).any(|a| index[a] >= grid.dims[a]) {
                return Err(parse_err(n, "voxel index outside the grid"));
            }
            let record = VoxelRecord {
                index,
                cell: VoxelCell {
                    implanted: num(f[6], n)?,
                    displacements: num(f[7], n)?,
                    electronic_ev: num(f[8], n)?,
                    lattice_ev: num(f[9], n)?,
                    rest_ev: num(f[10], n)?,
                },
            };
            if let Some(last) = voxels.last() {
                if grid.flat_index(last.index) >= grid.flat_index(index) {
                    return Err(parse_err(
                        n,
                        "voxels must be in increasing flat-index order",
                    ));
                }
            }
            voxels.push(record);
        }
        Ok(Self {
            grid,
            histories: histories.ok_or_else(|| missing("histories"))?,
            voxels,
            outside: outside.ok_or_else(|| missing("outside"))?,
        })
    }

    /// Write the legacy ASCII VTK export: an `UNSTRUCTURED_GRID` with one
    /// `VTK_VOXEL` cell per occupied voxel and the tallies as cell scalars.
    /// The grid and the history count are in the title line. See the
    /// [module docs](self).
    pub fn write_vtk<W: Write>(&self, mut w: W) -> io::Result<()> {
        let g = &self.grid;
        let title = format!(
            "lindhard voxel tally v1 origin={} spacing={} dims={},{},{} histories={}",
            join3(g.origin_m),
            join3(g.spacing_m),
            g.dims[0],
            g.dims[1],
            g.dims[2],
            self.histories
        );
        if title.len() > 255 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "the VTK title line would exceed 255 characters",
            ));
        }
        let n = self.voxels.len();
        writeln!(w, "# vtk DataFile Version 3.0")?;
        writeln!(w, "{title}")?;
        writeln!(w, "ASCII")?;
        writeln!(w, "DATASET UNSTRUCTURED_GRID")?;
        writeln!(w, "POINTS {} double", 8 * n)?;
        for v in &self.voxels {
            for c in VOXEL_CORNERS {
                let p = g.corner_m(std::array::from_fn(|a| v.index[a] + c[a]));
                writeln!(w, "{:e} {:e} {:e}", p[0], p[1], p[2])?;
            }
        }
        writeln!(w, "CELLS {} {}", n, 9 * n)?;
        for i in 0..n {
            let b = 8 * i;
            writeln!(
                w,
                "8 {} {} {} {} {} {} {} {}",
                b,
                b + 1,
                b + 2,
                b + 3,
                b + 4,
                b + 5,
                b + 6,
                b + 7
            )?;
        }
        writeln!(w, "CELL_TYPES {n}")?;
        for _ in 0..n {
            writeln!(w, "{VTK_VOXEL}")?;
        }
        writeln!(w, "CELL_DATA {n}")?;
        for (name, kind) in VTK_SCALARS {
            writeln!(w, "SCALARS {name} {kind} 1")?;
            writeln!(w, "LOOKUP_TABLE default")?;
            for v in &self.voxels {
                let c = &v.cell;
                match name {
                    "implanted" => writeln!(w, "{}", c.implanted)?,
                    "displacements" => writeln!(w, "{}", c.displacements)?,
                    "electronic_ev" => writeln!(w, "{:e}", c.electronic_ev)?,
                    "lattice_ev" => writeln!(w, "{:e}", c.lattice_ev)?,
                    "rest_ev" => writeln!(w, "{:e}", c.rest_ev)?,
                    _ => writeln!(w, "{:e}", c.deposited_ev())?,
                }
            }
        }
        Ok(())
    }

    /// The VTK export as a string.
    pub fn to_vtk(&self) -> String {
        let mut buf = Vec::new();
        self.write_vtk(&mut buf).expect("writing to a Vec");
        String::from_utf8(buf).expect("the VTK file is ASCII")
    }

    /// Read back a file written by [`Self::write_vtk`]. The result has an
    /// empty [`Self::outside`], which the VTK file does not carry.
    pub fn read_vtk(text: &str) -> Result<Self, VoxelTallyError> {
        let mut lines = text.lines();
        if lines.next().map(str::trim) != Some("# vtk DataFile Version 3.0") {
            return Err(parse_err(1, "not a legacy VTK 3.0 file"));
        }
        let title = lines.next().ok_or_else(|| parse_err(2, "no title line"))?;
        let mut fields = title.split_whitespace();
        if fields.by_ref().take(4).collect::<Vec<_>>() != ["lindhard", "voxel", "tally", "v1"] {
            return Err(parse_err(2, "not a lindhard voxel tally"));
        }
        let (mut origin, mut spacing, mut dims, mut histories) = (None, None, None, None);
        for f in fields {
            match f.split_once('=') {
                Some(("origin", v)) => origin = Some(parse3::<f64>(v, 2)?),
                Some(("spacing", v)) => spacing = Some(parse3::<f64>(v, 2)?),
                Some(("dims", v)) => dims = Some(parse3::<usize>(v, 2)?),
                Some(("histories", v)) => histories = Some(num::<u64>(v, 2)?),
                _ => return Err(parse_err(2, format!("unexpected title field {f:?}"))),
            }
        }
        let missing = |what: &str| parse_err(2, format!("title is missing `{what}`"));
        let grid = VoxelTallyGrid::new(
            origin.ok_or_else(|| missing("origin"))?,
            spacing.ok_or_else(|| missing("spacing"))?,
            dims.ok_or_else(|| missing("dims"))?,
        )?;
        let histories = histories.ok_or_else(|| missing("histories"))?;

        // The rest is a token stream; each token remembers its line.
        let mut tk = Tokens {
            iter: Box::new(
                text.lines()
                    .enumerate()
                    .skip(2)
                    .flat_map(|(i, l)| l.split_whitespace().map(move |w| (i + 1, w))),
            ),
            last_line: text.lines().count().max(1),
        };
        tk.expect("ASCII")?;
        tk.expect("DATASET")?;
        tk.expect("UNSTRUCTURED_GRID")?;
        tk.expect("POINTS")?;
        let (n, t) = tk.next("a point count")?;
        let n_points: usize = num(t, n)?;
        if !n_points.is_multiple_of(8) {
            return Err(parse_err(n, "the point count must be a multiple of 8"));
        }
        tk.expect("double")?;
        let n_cells = n_points / 8;
        let mut corners = Vec::with_capacity(n_cells);
        for _ in 0..n_cells {
            let mut lo = [0.0f64; 3];
            for p in 0..8 {
                for coord in &mut lo {
                    let (n, t) = tk.next("a coordinate")?;
                    let v: f64 = num(t, n)?;
                    if p == 0 {
                        *coord = v;
                    }
                }
            }
            corners.push(lo);
        }
        tk.expect("CELLS")?;
        for want in [n_cells, 9 * n_cells] {
            let (n, t) = tk.next("a cell count")?;
            if num::<usize>(t, n)? != want {
                return Err(parse_err(n, "inconsistent cell counts"));
            }
        }
        for _ in 0..9 * n_cells {
            tk.next("a connectivity entry")?;
        }
        tk.expect("CELL_TYPES")?;
        let (n, t) = tk.next("a cell count")?;
        if num::<usize>(t, n)? != n_cells {
            return Err(parse_err(n, "inconsistent cell counts"));
        }
        for _ in 0..n_cells {
            let (n, t) = tk.next("a cell type")?;
            if num::<u8>(t, n)? != VTK_VOXEL {
                return Err(parse_err(n, "only VTK_VOXEL (11) cells are supported"));
            }
        }
        tk.expect("CELL_DATA")?;
        let (n, t) = tk.next("a cell count")?;
        if num::<usize>(t, n)? != n_cells {
            return Err(parse_err(n, "inconsistent cell counts"));
        }
        let mut cells = vec![VoxelCell::default(); n_cells];
        for (name, kind) in VTK_SCALARS {
            tk.expect("SCALARS")?;
            tk.expect(name)?;
            tk.expect(kind)?;
            tk.expect("1")?;
            tk.expect("LOOKUP_TABLE")?;
            tk.expect("default")?;
            for cell in &mut cells {
                let (n, t) = tk.next("a scalar")?;
                match name {
                    "implanted" => cell.implanted = num(t, n)?,
                    "displacements" => cell.displacements = num(t, n)?,
                    "electronic_ev" => cell.electronic_ev = num(t, n)?,
                    "lattice_ev" => cell.lattice_ev = num(t, n)?,
                    "rest_ev" => cell.rest_ev = num(t, n)?,
                    _ => {} // deposited_ev is derived
                }
            }
        }
        let voxels = corners
            .iter()
            .zip(cells)
            .map(|(lo, cell)| {
                let mut index = [0usize; 3];
                for a in 0..3 {
                    let s = ((lo[a] - grid.origin_m[a]) / grid.spacing_m[a]).round();
                    if !(s >= 0.0 && s < grid.dims[a] as f64) {
                        return Err(parse_err(1, "a cell lies outside the grid"));
                    }
                    index[a] = s as usize;
                }
                Ok(VoxelRecord { index, cell })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            grid,
            histories,
            voxels,
            outside: VoxelCell::default(),
        })
    }
}
