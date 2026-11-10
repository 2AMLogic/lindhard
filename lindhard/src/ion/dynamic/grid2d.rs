//! A two-dimensional (depth–lateral) dynamic target: a laterally finite beam
//! on a surface that recedes non-uniformly with fluence.
//!
//! [`DynamicRun2d`] is the 2-D counterpart of [`DynamicRun`](super::DynamicRun).
//! The target is a row of **columns** of lateral width `h`; each column is a
//! complete 1-D [`CompositionGrid`] (slabs, areal inventories, a
//! [`Relaxation`](super::Relaxation) convention, an optional immutable
//! substrate), and the **surface height profile** is one number per column,
//! the cumulative recession `R_j` of its front (m, positive into the target;
//! negative where redeposition has built the surface up). Nothing is
//! three-dimensional: the target is uniform along the third axis `z`, and a
//! column stands for the target between lateral positions `j h` and
//! `(j+1) h`. Inventories are atoms/m² of column area, so a column behaves
//! exactly as a 1-D run under the beam it receives, and with identical columns
//! under a uniform beam the run reduces to the 1-D result (see "Uniform-beam
//! limit"). Coordinates follow the engine: `x` is the depth, `y` the lateral
//! coordinate (`[0, n h]`), `z` the invariant axis. (The `x–z` plane of the
//! issue that introduced this module is the `x–y` plane here.)
//!
//! Only published descriptions of the physics were used: sputter erosion and
//! the dynamic target follow the 1-D module ([`DynamicRun`](super::DynamicRun),
//! TRIDYN, Moller and Eckstein, Nucl. Instrum. Methods B 2 (1984) 814, and
//! SDTrimSP, Mutzke et al., IPP report 2019-02). The faceted height-field
//! geometry, the redeposition rule and the local-angle rule below are this
//! crate's own statements, not taken from a published code.
//!
//! # Surface and transport geometry
//!
//! The transport geometry of a step is a [`Geometry`] built from the current
//! columns. The **surface** is the piecewise-linear curve through the
//! vertices `(y_j, R_j)` at the column centres `y_j = (j + 1/2) h` (extended
//! periodically, or flat beyond the first and last centre for
//! [`Boundary::Vacuum`]). Material lies at depths `x >= s(y)`. Inside column
//! `j` the layer boundaries are flat at `R_j + (cumulative slab thickness)`;
//! the top slab extends up to the surface. A particle reaching the surface
//! escapes or reflects against the **normal of the surface segment it
//! crosses**, so the angle of incidence of a primary, the surface-barrier
//! refraction and the angular distribution of sputtered atoms all follow
//! from the local surface orientation without any tabulated angular yield.
//! Column faces and layer boundaries are material events as in
//! [`VoxelGrid`](crate::geometry::VoxelGrid) (a face between columns is an
//! event even where the materials coincide; the unused part of the flight is
//! kept, as documented in [`crate::ion::bca`]).
//!
//! **Local-angle rule.** The local surface normal on a segment between two
//! neighbouring columns is that of the line through their two vertices, i.e.
//! the slope is `(R_{j+1} - R_j) / h`; at a column centre the beam's
//! incidence angle is reported ([`Step2d::incidence_rad`]) from the
//! central-difference slope `m_j = (R_{j+1} - R_{j-1}) / (2 h)`, so that
//! `cos(alpha) = (d_x - m_j d_y) / sqrt(1 + m_j^2)` for the beam direction
//! `d` (one-sided at a [`Boundary::Vacuum`] end). The yield of a column then
//! is a result of the transport against the faceted surface: it is the yield
//! at the incidence angle the primaries have on the segments they strike, and
//! the per-column yields are reported ([`Step2d::column_yield`]). This is
//! exact for flat columns and a piecewise-linear approximation of a curved
//! surface otherwise; sub-column curvature is not resolved.
//!
//! # Beam
//!
//! A lateral profile ([`BeamProfile`]) gives the relative fluence `p_j` at
//! each column centre: uniform, or Gaussian of stated full width at half
//! maximum, `p = exp(-d^2 / (2 sigma^2))` with `sigma = FWHM / (2 sqrt(2 ln 2))`
//! and `d` the (nearest-image, for a periodic axis) distance from the beam
//! centre. A column is **illuminated** if `p_j >= PROFILE_CUTOFF` (1e-4), so
//! the far tails are dropped. [`Dynamic2dConfig::fluence_m2`] is the **peak**
//! fluence (at `p = 1`), and the beam's `count` is the number of primaries
//! **per illuminated column** over the run; a primary in column `j` stands for
//! `fluence_m2 p_j / count` ions/m². All primaries of a column enter at the
//! centre of the column, so the beam's lateral resolution is the column width
//! (a wider-than-column spot is represented by many columns, not by sub-column
//! sampling). Primary `k` of column `j` draws from
//! [`stream`](crate::rng::stream)`(seed, j * count + k)`.
//!
//! # Composition changes and relaxation
//!
//! The inventory bookkeeping is that of the 1-D module and is delegated to
//! each column's [`CompositionGrid`], so both relaxation conventions
//! ([`Relaxation::IdealMixing`](super::Relaxation::IdealMixing) and
//! [`Relaxation::FixedNumberDensity`](super::Relaxation::FixedNumberDensity))
//! are supported unchanged: swelling pushes the interior and the back of a
//! column down while its front stays at `R_j`. An atom displaced in one column
//! that comes to rest in another changes both columns (lateral mixing is
//! transport, not diffusion; no surface diffusion, charging or swelling-driven
//! surface motion is modelled).
//!
//! # Erosion
//!
//! As in the 1-D module, a sputtered atom is removed from the **front of the
//! column it was displaced in** (the origin loss is cancelled and `Y_Z` atoms
//! of each element are taken slab by slab from the front, never more than the
//! column holds), and the recession of the column is the thickness the removed
//! inventory occupied under its relaxation convention. A column whose slabs
//! have all been removed is the bare substrate: the substrate is immutable, so
//! that column stops receding there (the hole has reached the substrate).
//!
//! # Redeposition
//!
//! One rule, deterministic and atom-conserving: a sputtered atom travels from
//! its exit point (position and direction after the surface refraction, as
//! reported by the transport) in a **straight line in vacuum** (no gas-phase
//! collisions) and **sticks with probability 1 where it first meets the
//! surface again**. It is then added to the top slab of the column it hit, as
//! one atom of its own element, and the front of that column rises by the
//! thickness it occupies (`R_j` decreases by the thickness the deposited
//! inventory takes under the relaxation convention, so the interfaces below
//! stay where they are). An atom whose line never meets the surface (it
//! leaves above the highest point of the profile, leaves through a
//! [`Boundary::Vacuum`] lateral end, hits a column that is bare substrate, or
//! travels more than `16 n + 16` surface segments) is lost. The rule works from
//! the exit position and direction distribution that the transport already
//! delivers; there is no sticking-coefficient or angular re-emission model.
//! Only atoms displaced from a finite slab are redeposited (an atom of the
//! substrate was never part of any inventory). Per step the run reports the
//! sputtered and redeposited atoms per metre of `z` and their ratio, the
//! redeposited fraction ([`Step2d::redeposited_fraction`]); atoms are
//! conserved: for each element the change of the total inventory of all
//! columns is the implanted atoms minus the sputtered ones plus the
//! redeposited ones.
//!
//! # Determinism
//!
//! Per column the transport is [`Bca::run_range`]; the integer tallies of each
//! column are scaled by the column's weight and summed in column order, and
//! the deltas are applied in `(column, slab, Z)` order, so a run is
//! bit-identical at any thread count. A step that covers primaries
//! `first..first + n` of each column never replays earlier indices; adaptive
//! shrinking follows the 1-D rules ([`StepPolicy`]) with the largest relative
//! change taken over all slabs of all columns. A step that makes a removal
//! exceed what a slab holds is halved down to `min_ions_per_step`, where the
//! removal is capped at what the slab holds and counted in [`Step2d::clamped`]
//! (as in the 1-D run, so inventory is not conserved for that step by up to
//! one step's worth of removal; the last atoms of a column that is being cut
//! through); a fixed policy returns the error instead.
//!
//! # Uniform-beam limit
//!
//! With identical columns, a [`BeamProfile::Uniform`] beam and a periodic
//! lateral axis, every column receives statistically the same history as a
//! 1-D run, and the surface stays flat to within the counting noise: the
//! recession and the composition–depth profile averaged over columns agree
//! with [`DynamicRun`](super::DynamicRun) within the statistical error
//! (committed test `tests/dynamic_2d.rs`). Differences at the noise level of
//! a column grow into the surface roughness a real uniform beam would also
//! develop.

use std::collections::BTreeMap;

use crate::geometry::{Boundary, Exit, ExitOutcome, Face, Flight, Geometry};
use crate::ion::bca::{Bca, BcaConfig, BcaTally, Beam, Particle};
use crate::ion::scattering::ScatteringTable;
use crate::ion::stopping::ElectronicStopping;
use crate::material::Material;

use super::adapter::Yields;
use super::driver::{DynamicRunError, StepPolicy};
use super::{CompositionGrid, DynamicError, InventoryDelta};

/// Columns whose relative fluence is below this are not illuminated.
pub const PROFILE_CUTOFF: f64 = 1.0e-4;

/// Lateral beam profile (relative fluence at a lateral position).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BeamProfile {
    /// The same fluence in every column.
    Uniform,
    /// Gaussian of full width at half maximum `fwhm_m`, centred at `center_m`
    /// (lateral coordinate `y`).
    Gaussian {
        /// Centre, m.
        center_m: f64,
        /// Full width at half maximum, m, finite and positive.
        fwhm_m: f64,
    },
}

/// Settings of a 2-D run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dynamic2dConfig {
    /// Peak fluence (at relative fluence 1), ions/m², represented by `count`
    /// primaries per illuminated column.
    pub fluence_m2: f64,
    /// Step sizes, in primaries per illuminated column.
    pub policy: StepPolicy,
    /// Lateral beam profile.
    pub profile: BeamProfile,
    /// Behaviour at the two lateral ends. [`Boundary::Periodic`] makes the
    /// target an infinite periodic array; [`Boundary::Vacuum`] a finite strip.
    pub lateral: Boundary,
}

/// What one accepted step of a 2-D run did.
#[derive(Debug, Clone, PartialEq)]
pub struct Step2d {
    /// Step number, from 1.
    pub step: u64,
    /// Index of the first primary of the step within each column.
    pub first_index: u64,
    /// Primaries per illuminated column in the step.
    pub ions: u64,
    /// Attempts made at this `first_index`.
    pub attempts: u32,
    /// Largest relative composition change of the step over all slabs.
    pub max_change: f64,
    /// `(column, slab, element)` removals capped at what the slab held (only at
    /// the smallest step of an adaptive policy; 0 otherwise).
    pub clamped: u32,
    /// Event counts of the step, summed over the illuminated columns.
    pub yields: Yields,
    /// Surface position after the step: cumulative recession per column, m
    /// (positive into the target). This is `z_s(y)` of the issue, in depth
    /// convention.
    pub surface_m: Vec<f64>,
    /// Beam incidence angle from the local surface normal at the start of the
    /// step, per column, radians (local-angle rule, module docs).
    pub incidence_rad: Vec<f64>,
    /// Sputtered atoms per primary of the column the beam entered, per column
    /// (0 for columns that are not illuminated).
    pub column_yield: Vec<f64>,
    /// Atoms per metre of `z` sputtered in the step (from finite slabs).
    pub sputtered_atoms_per_m: f64,
    /// Atoms per metre of `z` of those that were redeposited.
    pub redeposited_atoms_per_m: f64,
    /// Beam atoms per metre of `z` that came to rest in a finite slab.
    pub implanted_atoms_per_m: f64,
    /// Atoms that left through a lateral end (a [`Boundary::Vacuum`] axis):
    /// primaries and recoils, event count.
    pub lateral_escapes: u64,
}

impl Step2d {
    /// Redeposited over sputtered atoms in this step (0 if nothing was
    /// sputtered).
    pub fn redeposited_fraction(&self) -> f64 {
        if self.sputtered_atoms_per_m > 0.0 {
            self.redeposited_atoms_per_m / self.sputtered_atoms_per_m
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

/// The transport geometry of one step: faceted height-field surface over
/// columns of flat layers. See the module docs.
#[derive(Debug)]
struct ColumnTarget {
    h: f64,
    n_cols: usize,
    boundary: Boundary,
    /// Surface depth at the column centres.
    surface: Vec<f64>,
    /// Region-indexed materials (region = column offset + layer).
    materials: Vec<Material>,
    /// First region of each column; one extra entry at the end.
    offset: Vec<usize>,
    /// Number of finite slabs per column (the substrate, if any, is the layer
    /// after them).
    n_slabs: Vec<usize>,
    /// `tops[j][k]`: depth of the upper boundary of layer `k` of column `j`
    /// (entry 0 is the nominal front and not used for lookups).
    tops: Vec<Vec<f64>>,
    /// Depth of the back face of each column without a substrate.
    back: Vec<Option<f64>>,
    min_surface: f64,
}

impl ColumnTarget {
    fn new(
        grids: &[CompositionGrid],
        surface: &[f64],
        h: f64,
        boundary: Boundary,
    ) -> Result<Self, DynamicError> {
        let n_cols = grids.len();
        let mut materials = Vec::new();
        let mut offset = Vec::with_capacity(n_cols + 1);
        let mut n_slabs = Vec::with_capacity(n_cols);
        let mut tops = Vec::with_capacity(n_cols);
        let mut back = Vec::with_capacity(n_cols);
        for (j, g) in grids.iter().enumerate() {
            offset.push(materials.len());
            let stack = g.to_stack()?;
            let mut x = surface[j];
            let mut t = Vec::new();
            let mut ns = 0;
            let mut bk = None;
            for (k, layer) in stack.layers().iter().enumerate() {
                t.push(x);
                materials.push(layer.material().clone());
                if layer.thickness_m().is_finite() {
                    x += layer.thickness_m();
                    ns += 1;
                } else {
                    debug_assert_eq!(k + 1, stack.layers().len());
                }
            }
            if ns == stack.layers().len() {
                bk = Some(x);
            }
            n_slabs.push(ns);
            tops.push(t);
            back.push(bk);
        }
        offset.push(materials.len());
        let min_surface = surface.iter().copied().fold(f64::INFINITY, f64::min);
        Ok(Self {
            h,
            n_cols,
            boundary,
            surface: surface.to_vec(),
            materials,
            offset,
            n_slabs,
            tops,
            back,
            min_surface,
        })
    }

    fn width(&self) -> f64 {
        self.n_cols as f64 * self.h
    }

    /// Column and layer of a region index.
    fn split(&self, region: usize) -> (usize, usize) {
        let j = self.offset.partition_point(|&o| o <= region) - 1;
        (j, region - self.offset[j])
    }

    fn col_of(&self, y: f64) -> Option<usize> {
        let l = self.width();
        if !(y >= 0.0 && y <= l) {
            return None;
        }
        if y == l {
            return match self.boundary {
                Boundary::Vacuum => None,
                Boundary::Periodic => Some(0),
            };
        }
        Some(((y / self.h) as usize).min(self.n_cols - 1))
    }

    /// Layer of column `j` at depth `x` (`None` beyond the back face of a
    /// column without substrate). The top layer extends up to the surface.
    fn layer_of(&self, j: usize, x: f64) -> Option<usize> {
        if let Some(b) = self.back[j] {
            if x >= b {
                return None;
            }
        }
        let t = &self.tops[j];
        Some(t[1..].partition_point(|&b| b <= x))
    }

    /// Surface depths at the two neighbouring centres of column `j`.
    fn neighbours(&self, j: usize) -> (f64, f64) {
        let n = self.n_cols;
        let left = if j > 0 {
            self.surface[j - 1]
        } else if self.boundary == Boundary::Periodic {
            self.surface[n - 1]
        } else {
            self.surface[0]
        };
        let right = if j + 1 < n {
            self.surface[j + 1]
        } else if self.boundary == Boundary::Periodic {
            self.surface[0]
        } else {
            self.surface[n - 1]
        };
        (left, right)
    }

    /// Surface depth and slope at lateral position `y` in column `j`.
    fn surface_at(&self, j: usize, y: f64) -> (f64, f64) {
        let yc = (j as f64 + 0.5) * self.h;
        let (l, r) = self.neighbours(j);
        let rj = self.surface[j];
        let m = if y < yc {
            (rj - l) / self.h
        } else {
            (r - rj) / self.h
        };
        (rj + (y - yc) * m, m)
    }

    /// First crossing of the surface from inside (material) to outside within
    /// path length `tmax`, with the outward unit normal there.
    fn surface_crossing(
        &self,
        j: usize,
        pos: [f64; 3],
        dir: [f64; 3],
        tmax: f64,
    ) -> Option<(f64, [f64; 3])> {
        let (x0, y0) = (pos[0], pos[1]);
        let (dx, dy) = (dir[0], dir[1]);
        let yc = (j as f64 + 0.5) * self.h;
        let (left, right) = self.neighbours(j);
        let rj = self.surface[j];
        let t_c = if dy != 0.0 {
            (yc - y0) / dy
        } else {
            f64::INFINITY
        };
        let mut lower = y0 < yc || (y0 == yc && dy < 0.0);
        let mut t0 = 0.0;
        loop {
            let m = if lower {
                (rj - left) / self.h
            } else {
                (right - rj) / self.h
            };
            let switches = (lower && dy > 0.0) || (!lower && dy < 0.0);
            let t1 = if switches && t_c > t0 && t_c < tmax {
                t_c
            } else {
                tmax
            };
            let g0 = x0 + t0 * dx - (rj + (y0 + t0 * dy - yc) * m);
            let slope = dx - m * dy;
            let g1 = g0 + slope * (t1 - t0);
            if slope < 0.0 && g1 < 0.0 {
                let t = if g0 <= 0.0 { t0 } else { t0 - g0 / slope };
                let norm = (1.0 + m * m).sqrt();
                return Some((t.min(tmax), [-1.0 / norm, m / norm, 0.0]));
            }
            if t1 >= tmax {
                return None;
            }
            t0 = t1;
            lower = !lower;
        }
    }

    /// Column that a sputtered atom leaving `pos` along `dir` meets first, if
    /// any (straight line in vacuum; see "Redeposition" in the module docs).
    fn trace_to_surface(&self, pos: [f64; 3], dir: [f64; 3]) -> Option<usize> {
        let (x0, dx, dy) = (pos[0], dir[0], dir[1]);
        if !(x0.is_finite() && dx.is_finite() && dy.is_finite()) {
            return None;
        }
        let mut j = self.col_of(pos[1].clamp(0.0, self.width()))?;
        let y_start = pos[1];
        let mut t0 = 0.0;
        // `y` at the start of the current piece.
        let mut y = y_start;
        let mut first = true;
        for _ in 0..(16 * self.n_cols + 16) {
            let yc = (j as f64 + 0.5) * self.h;
            let lower = y < yc || (y == yc && dy < 0.0);
            let (l, r) = self.neighbours(j);
            let rj = self.surface[j];
            let m = if lower {
                (rj - l) / self.h
            } else {
                (r - rj) / self.h
            };
            let (y_edge, at_centre) = if dy > 0.0 {
                if lower {
                    (yc, true)
                } else {
                    ((j as f64 + 1.0) * self.h, false)
                }
            } else if dy < 0.0 {
                if lower {
                    (j as f64 * self.h, false)
                } else {
                    (yc, true)
                }
            } else {
                (y, true)
            };
            let t1 = if dy == 0.0 {
                f64::INFINITY
            } else {
                t0 + ((y_edge - y) / dy).max(0.0)
            };
            let mut g0 = x0 + t0 * dx - (rj + (y - yc) * m);
            if first {
                g0 = g0.min(0.0);
                first = false;
            }
            let slope = dx - m * dy;
            if slope > 0.0 {
                let t_hit = t0 + (-g0).max(0.0) / slope;
                if t_hit <= t1 {
                    return Some(j);
                }
            }
            if !t1.is_finite() {
                return None;
            }
            // The ray rises above every surface point and keeps rising.
            if dx < 0.0 && x0 + t1 * dx < self.min_surface {
                return None;
            }
            t0 = t1;
            y = y_edge;
            if !at_centre {
                if dy > 0.0 {
                    if j + 1 == self.n_cols {
                        match self.boundary {
                            Boundary::Vacuum => return None,
                            Boundary::Periodic => {
                                j = 0;
                                y = 0.0;
                            }
                        }
                    } else {
                        j += 1;
                    }
                } else if j == 0 {
                    match self.boundary {
                        Boundary::Vacuum => return None,
                        Boundary::Periodic => {
                            j = self.n_cols - 1;
                            y = self.width();
                        }
                    }
                } else {
                    j -= 1;
                }
            }
        }
        None
    }
}

impl Geometry for ColumnTarget {
    fn n_materials(&self) -> usize {
        self.materials.len()
    }

    fn material(&self, index: usize) -> &Material {
        &self.materials[index]
    }

    fn n_regions(&self) -> usize {
        self.materials.len()
    }

    fn material_index(&self, region: usize) -> usize {
        region
    }

    fn locate(&self, pos: [f64; 3]) -> Option<usize> {
        if !pos.iter().all(|c| c.is_finite()) {
            return None;
        }
        let j = self.col_of(pos[1])?;
        let (s, _) = self.surface_at(j, pos[1].min(self.width()));
        if pos[0] < s {
            return None;
        }
        let k = self.layer_of(j, pos[0])?;
        Some(self.offset[j] + k)
    }

    fn entry_point(&self) -> [f64; 3] {
        [self.surface[0] + 1.0e-6 * self.h, 0.5 * self.h, 0.0]
    }

    fn exit(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Option<Exit> {
        match self.flight(region, pos, dir, limit) {
            Flight::Event(e) => Some(e),
            Flight::Clear { .. } => None,
        }
    }

    fn flight(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Flight {
        if !dir.iter().all(|d| d.is_finite()) || dir == [0.0; 3] {
            return Flight::Clear { region };
        }
        #[derive(Clone, Copy)]
        enum Kind {
            Lateral(bool),
            Interface(bool),
            Back,
        }
        let (j, k) = self.split(region);
        let (dx, dy) = (dir[0], dir[1]);
        let mut t_end = f64::INFINITY;
        let mut kind = None;
        let mut consider = |t: f64, kd: Kind| {
            if t < t_end {
                t_end = t;
                kind = Some(kd);
            }
        };
        if dy > 0.0 {
            consider(
                (((j as f64 + 1.0) * self.h - pos[1]) / dy).max(0.0),
                Kind::Lateral(true),
            );
        } else if dy < 0.0 {
            consider(
                ((j as f64 * self.h - pos[1]) / dy).max(0.0),
                Kind::Lateral(false),
            );
        }
        if dx > 0.0 {
            let n_layers = self.offset[j + 1] - self.offset[j];
            if k + 1 < n_layers {
                consider(
                    ((self.tops[j][k + 1] - pos[0]) / dx).max(0.0),
                    Kind::Interface(true),
                );
            } else if let Some(b) = self.back[j] {
                consider(((b - pos[0]) / dx).max(0.0), Kind::Back);
            }
        } else if dx < 0.0 && k > 0 {
            consider(
                ((self.tops[j][k] - pos[0]) / dx).max(0.0),
                Kind::Interface(false),
            );
        }
        let tmax = t_end.min(limit);
        let point = |t: f64| [pos[0] + t * dx, pos[1] + t * dy, pos[2] + t * dir[2]];
        if let Some((t, normal)) = self.surface_crossing(j, pos, dir, tmax) {
            return Flight::Event(Exit {
                distance: t,
                at: point(t),
                outcome: ExitOutcome::Escape {
                    face: Face::Front,
                    normal,
                },
            });
        }
        let Some(kind) = kind.filter(|_| t_end <= limit) else {
            return Flight::Clear { region };
        };
        let mut at = point(t_end);
        let outcome = match kind {
            Kind::Interface(down) => ExitOutcome::Enter {
                region: if down { region + 1 } else { region - 1 },
                pos: at,
            },
            Kind::Back => ExitOutcome::Escape {
                face: Face::Back,
                normal: [1.0, 0.0, 0.0],
            },
            Kind::Lateral(up) => {
                at[1] = if up {
                    (j as f64 + 1.0) * self.h
                } else {
                    j as f64 * self.h
                };
                let mut new_pos = at;
                let next = if up {
                    if j + 1 == self.n_cols {
                        (self.boundary == Boundary::Periodic).then(|| {
                            new_pos[1] = 0.0;
                            0
                        })
                    } else {
                        Some(j + 1)
                    }
                } else if j == 0 {
                    (self.boundary == Boundary::Periodic).then(|| {
                        new_pos[1] = self.width();
                        self.n_cols - 1
                    })
                } else {
                    Some(j - 1)
                };
                let normal = [0.0, if up { 1.0 } else { -1.0 }, 0.0];
                match next {
                    None => ExitOutcome::Escape {
                        face: Face::Side,
                        normal,
                    },
                    Some(nj) => match self.layer_of(nj, at[0]) {
                        Some(nk) => ExitOutcome::Enter {
                            region: self.offset[nj] + nk,
                            pos: new_pos,
                        },
                        // Beyond the back face of the neighbouring column.
                        None => ExitOutcome::Escape {
                            face: Face::Back,
                            normal,
                        },
                    },
                }
            }
        };
        Flight::Event(Exit {
            distance: t_end,
            at,
            outcome,
        })
    }

    fn in_vacuum(&self, pos: [f64; 3]) -> bool {
        match self.col_of(pos[1]) {
            None => true,
            Some(j) => pos[0] < self.surface_at(j, pos[1].min(self.width())).0,
        }
    }
}

// ---------------------------------------------------------------------------
// Tally
// ---------------------------------------------------------------------------

/// Per-column tally of one block of primaries: net atom counts per
/// `(column, slab, Z)`, sputtered atoms per origin, and where they land again.
struct Tally2d<'t> {
    target: &'t ColumnTarget,
    counts: BTreeMap<(usize, usize, u8), i64>,
    sputtered_from: BTreeMap<(usize, usize, u8), u64>,
    redeposited: BTreeMap<(usize, u8), u64>,
    lateral: u64,
    yields: Yields,
}

impl<'t> Tally2d<'t> {
    fn new(target: &'t ColumnTarget) -> Self {
        Self {
            target,
            counts: BTreeMap::new(),
            sputtered_from: BTreeMap::new(),
            redeposited: BTreeMap::new(),
            lateral: 0,
            yields: Yields::default(),
        }
    }

    /// `(column, slab)` of a region if it is a finite slab.
    fn slab(&self, region: usize) -> Option<(usize, usize)> {
        let (j, k) = self.target.split(region);
        (k < self.target.n_slabs[j]).then_some((j, k))
    }
}

impl BcaTally for Tally2d<'_> {
    fn begin_history(&mut self, _index: u64) {
        self.yields.histories += 1;
    }

    fn recoil(&mut self, recoil: &Particle) {
        self.yields.displaced += 1;
        if let Some((j, k)) = self.slab(recoil.layer) {
            *self.counts.entry((j, k, recoil.z)).or_insert(0) -= 1;
        }
    }

    fn stopped(&mut self, p: &Particle) {
        let slab = self.slab(p.layer);
        if let Some((j, k)) = slab {
            *self.counts.entry((j, k, p.z)).or_insert(0) += 1;
        }
        if p.is_primary() {
            if slab.is_some() {
                self.yields.primaries_in_slabs += 1;
            } else {
                self.yields.primaries_in_substrate += 1;
            }
        }
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        match (p.is_primary(), face) {
            (true, Face::Front) => self.yields.backscattered += 1,
            (true, Face::Back) => self.yields.transmitted += 1,
            (false, Face::Front) => {
                *self.yields.sputtered.entry(p.z).or_insert(0) += 1;
                if let Some((j, k)) = self.slab(p.origin_layer) {
                    *self.sputtered_from.entry((j, k, p.z)).or_insert(0) += 1;
                    if let Some(c) = self.target.trace_to_surface(p.pos, p.dir) {
                        if self.target.n_slabs[c] > 0 {
                            *self.redeposited.entry((c, p.z)).or_insert(0) += 1;
                        }
                    }
                }
            }
            (false, Face::Back) => self.yields.recoils_transmitted += 1,
            (_, Face::Side) => self.lateral += 1,
        }
    }

    fn merge(&mut self, other: Self) {
        for (k, n) in other.counts {
            *self.counts.entry(k).or_insert(0) += n;
        }
        for (k, n) in other.sputtered_from {
            *self.sputtered_from.entry(k).or_insert(0) += n;
        }
        for (k, n) in other.redeposited {
            *self.redeposited.entry(k).or_insert(0) += n;
        }
        self.lateral += other.lateral;
        self.yields.add(&other.yields);
    }
}

// ---------------------------------------------------------------------------
// The run
// ---------------------------------------------------------------------------

/// A fluence-stepped 2-D run. See the module documentation.
pub struct DynamicRun2d<'a> {
    grids: Vec<CompositionGrid>,
    surface: Vec<f64>,
    h: f64,
    beam: Beam,
    config: BcaConfig,
    stopping: &'a (dyn ElectronicStopping + Sync),
    table: &'a ScatteringTable,
    cfg: Dynamic2dConfig,
    /// Relative fluence per column (0 where not illuminated).
    weights: Vec<f64>,
    next_index: u64,
    steps: u64,
    current_n: u64,
    cumulative: Yields,
    sputtered_atoms_per_m: f64,
    redeposited_atoms_per_m: f64,
}

impl std::fmt::Debug for DynamicRun2d<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DynamicRun2d")
            .field("columns", &self.grids.len())
            .field("column_width_m", &self.h)
            .field("cfg", &self.cfg)
            .field("next_index", &self.next_index)
            .field("steps", &self.steps)
            .finish_non_exhaustive()
    }
}

impl<'a> DynamicRun2d<'a> {
    /// Set up a run on `columns` (front to back identical bookkeeping to the
    /// 1-D module; all columns must share the lateral width `column_width_m`).
    /// `beam.count` is the number of primaries per illuminated column.
    pub fn new(
        columns: Vec<CompositionGrid>,
        column_width_m: f64,
        beam: Beam,
        config: BcaConfig,
        stopping: &'a (dyn ElectronicStopping + Sync),
        table: &'a ScatteringTable,
        cfg: Dynamic2dConfig,
    ) -> Result<Self, DynamicRunError> {
        let bad = |m: &str| Err(DynamicRunError::Invalid(m.to_string()));
        if columns.is_empty() {
            return bad("at least one column is needed");
        }
        if !(column_width_m.is_finite() && column_width_m > 0.0) {
            return bad("the column width must be finite and positive");
        }
        if beam.count == 0 {
            return bad("the beam needs at least one primary per column");
        }
        if (columns.len() as u64).checked_mul(beam.count).is_none() {
            return bad("columns times primaries overflows the index space");
        }
        if !(cfg.fluence_m2.is_finite() && cfg.fluence_m2 > 0.0) {
            return bad("fluence must be finite and positive");
        }
        let n = columns.len();
        let weights: Vec<f64> = match cfg.profile {
            BeamProfile::Uniform => vec![1.0; n],
            BeamProfile::Gaussian { center_m, fwhm_m } => {
                if !(center_m.is_finite() && fwhm_m.is_finite() && fwhm_m > 0.0) {
                    return bad("the Gaussian profile needs a finite centre and positive FWHM");
                }
                let sigma = fwhm_m / (2.0 * (2.0 * std::f64::consts::LN_2).sqrt());
                let l = n as f64 * column_width_m;
                (0..n)
                    .map(|j| {
                        let mut d = (j as f64 + 0.5) * column_width_m - center_m;
                        if cfg.lateral == Boundary::Periodic {
                            d -= l * (d / l).round();
                        }
                        let p = (-d * d / (2.0 * sigma * sigma)).exp();
                        if p >= PROFILE_CUTOFF {
                            p
                        } else {
                            0.0
                        }
                    })
                    .collect()
            }
        };
        if weights.iter().all(|&w| w == 0.0) {
            return bad("the beam illuminates no column");
        }
        let current_n = match cfg.policy {
            StepPolicy::Fixed { ions_per_step } => {
                if ions_per_step == 0 {
                    return bad("ions_per_step must be at least 1");
                }
                ions_per_step
            }
            StepPolicy::Adaptive {
                max_ions_per_step,
                min_ions_per_step,
                max_change,
            } => {
                if min_ions_per_step == 0 || max_ions_per_step < min_ions_per_step {
                    return bad("need 1 <= min_ions_per_step <= max_ions_per_step");
                }
                if !(max_change.is_finite() && max_change > 0.0) {
                    return bad("max_change must be finite and positive");
                }
                max_ions_per_step
            }
        };
        Ok(Self {
            surface: vec![0.0; n],
            grids: columns,
            h: column_width_m,
            beam,
            config,
            stopping,
            table,
            cfg,
            weights,
            next_index: 0,
            steps: 0,
            current_n,
            cumulative: Yields::default(),
            sputtered_atoms_per_m: 0.0,
            redeposited_atoms_per_m: 0.0,
        })
    }

    /// The current columns.
    pub fn columns(&self) -> &[CompositionGrid] {
        &self.grids
    }

    /// Surface position per column: cumulative recession, m.
    pub fn surface_m(&self) -> &[f64] {
        &self.surface
    }

    /// Relative fluence per column (0 where not illuminated).
    pub fn profile_weights(&self) -> &[f64] {
        &self.weights
    }

    /// Primaries per column used so far.
    pub fn ions_done(&self) -> u64 {
        self.next_index
    }

    /// Peak fluence delivered so far, ions/m².
    pub fn fluence_done_m2(&self) -> f64 {
        self.cfg.fluence_m2 * self.next_index as f64 / self.beam.count as f64
    }

    /// Accepted steps so far.
    pub fn steps_done(&self) -> u64 {
        self.steps
    }

    /// Event counts summed over the accepted steps and columns.
    pub fn cumulative(&self) -> &Yields {
        &self.cumulative
    }

    /// Redeposited over sputtered atoms, over the whole run so far.
    pub fn redeposited_fraction(&self) -> f64 {
        if self.sputtered_atoms_per_m > 0.0 {
            self.redeposited_atoms_per_m / self.sputtered_atoms_per_m
        } else {
            0.0
        }
    }

    /// Beam incidence angle from the local surface normal per column, radians
    /// (the local-angle rule of the module docs) for the current surface.
    pub fn incidence_rad(&self) -> Vec<f64> {
        let n = self.surface.len();
        let d = self.beam.direction();
        (0..n)
            .map(|j| {
                let periodic = self.cfg.lateral == Boundary::Periodic;
                let lo = if j > 0 {
                    Some(j - 1)
                } else if periodic {
                    Some(n - 1)
                } else {
                    None
                };
                let hi = if j + 1 < n {
                    Some(j + 1)
                } else if periodic {
                    Some(0)
                } else {
                    None
                };
                let (a, b, span) = match (lo, hi) {
                    (Some(l), Some(r)) => (self.surface[l], self.surface[r], 2.0),
                    (None, Some(r)) => (self.surface[j], self.surface[r], 1.0),
                    (Some(l), None) => (self.surface[l], self.surface[j], 1.0),
                    (None, None) => (0.0, 0.0, 1.0),
                };
                let m = (b - a) / (span * self.h);
                ((d[0] - m * d[1]) / (1.0 + m * m).sqrt())
                    .clamp(-1.0, 1.0)
                    .acos()
            })
            .collect()
    }

    fn max_change(&self, deltas: &BTreeMap<(usize, usize, u8), f64>) -> f64 {
        let totals: Vec<Vec<f64>> = self
            .grids
            .iter()
            .map(|g| {
                (0..g.n_slabs())
                    .map(|k| {
                        g.inventory(k)
                            .map_or(0.0, |v| v.iter().map(|&(_, a)| a).sum())
                    })
                    .collect()
            })
            .collect();
        deltas
            .iter()
            .map(|(&(j, k, _), d)| d.abs() / totals[j][k])
            .fold(0.0, f64::max)
    }

    /// The columns after applying `delta`; the run state is untouched.
    fn apply_all(
        &self,
        delta: &BTreeMap<(usize, usize, u8), f64>,
    ) -> Result<Vec<CompositionGrid>, DynamicError> {
        let mut grids = self.grids.clone();
        for (j, g) in grids.iter_mut().enumerate() {
            let ds: Vec<InventoryDelta> = delta
                .range((j, 0, 0)..=(j, usize::MAX, u8::MAX))
                .map(|(&(_, slab, z), &d)| InventoryDelta {
                    slab,
                    z,
                    delta_atoms_m2: d,
                })
                .collect();
            g.apply(&ds)?;
        }
        Ok(grids)
    }

    /// Run the next step. `Ok(None)` when all primaries are used. On error the
    /// run state is unchanged.
    pub fn step(&mut self) -> Result<Option<Step2d>, DynamicRunError> {
        let remaining = self.beam.count - self.next_index;
        if remaining == 0 {
            return Ok(None);
        }
        let first = self.next_index;
        let mut n = self.current_n.min(remaining);
        let mut attempts = 0u32;
        let incidence = self.incidence_rad();
        let n_cols = self.grids.len();
        loop {
            attempts += 1;
            let target = ColumnTarget::new(&self.grids, &self.surface, self.h, self.cfg.lateral)?;

            let mut delta: BTreeMap<(usize, usize, u8), f64> = BTreeMap::new();
            let mut sputtered: BTreeMap<(usize, u8), f64> = BTreeMap::new();
            let mut deposited: BTreeMap<(usize, u8), f64> = BTreeMap::new();
            let mut yields = Yields::default();
            let mut column_yield = vec![0.0; n_cols];
            let mut lateral = 0u64;
            let mut implanted = 0.0;
            #[allow(clippy::needless_range_loop)]
            for c in 0..n_cols {
                if self.weights[c] == 0.0 {
                    continue;
                }
                let beam = Beam {
                    count: n,
                    ..self.beam
                };
                let entry = [
                    self.surface[c] + 1.0e-6 * self.h,
                    (c as f64 + 0.5) * self.h,
                    0.0,
                ];
                let bca = Bca::new(beam, &target, self.config, self.stopping, self.table)?
                    .with_entry_point(entry)?;
                let tally: Tally2d<'_> =
                    bca.run_range(c as u64 * self.beam.count + first, n, || {
                        Tally2d::new(&target)
                    })?;
                // Ions/m² one primary of this column stands for.
                let w = self.cfg.fluence_m2 * self.weights[c] / self.beam.count as f64;
                for (&k, &cnt) in &tally.counts {
                    if cnt != 0 {
                        *delta.entry(k).or_insert(0.0) += cnt as f64 * w;
                    }
                }
                for (&(j, k, z), &cnt) in &tally.sputtered_from {
                    // Cancel the origin loss; the atoms come off the front.
                    *delta.entry((j, k, z)).or_insert(0.0) += cnt as f64 * w;
                    *sputtered.entry((j, z)).or_insert(0.0) += cnt as f64 * w;
                }
                for (&k, &cnt) in &tally.redeposited {
                    *deposited.entry(k).or_insert(0.0) += cnt as f64 * w;
                }
                column_yield[c] = tally.yields.sputtered_total() as f64 / n as f64;
                lateral += tally.lateral;
                implanted += tally.yields.primaries_in_slabs as f64 * w * self.h;
                yields.add(&tally.yields);
            }

            // Erosion: take the sputtered atoms from the front of the column
            // they came from, element by element.
            let mut removed: Vec<BTreeMap<u8, f64>> = vec![BTreeMap::new(); n_cols];
            for (&(j, z), &y) in &sputtered {
                let mut rem = y;
                for slab in 0..self.grids[j].n_slabs() {
                    if rem <= 0.0 {
                        break;
                    }
                    let inv = self.grids[j].inventory(slab).unwrap_or_default();
                    let have = inv
                        .iter()
                        .find(|&&(zz, _)| zz == z)
                        .map_or(0.0, |&(_, a)| a);
                    let d = delta.entry((j, slab, z)).or_insert(0.0);
                    let take = (have + *d).max(0.0).min(rem);
                    if take > 0.0 {
                        *d -= take;
                        rem -= take;
                        *removed[j].entry(z).or_insert(0.0) += take;
                    }
                }
            }
            // Redeposition: onto the top slab of the column hit.
            let mut dep_inv: Vec<BTreeMap<u8, f64>> = vec![BTreeMap::new(); n_cols];
            for (&(j, z), &a) in &deposited {
                *delta.entry((j, 0, z)).or_insert(0.0) += a;
                *dep_inv[j].entry(z).or_insert(0.0) += a;
            }
            let mut step_recession = vec![0.0; n_cols];
            for j in 0..n_cols {
                let relax = self.grids[j].relaxation();
                if !removed[j].is_empty() {
                    step_recession[j] += relax.thickness_m(&removed[j])?;
                }
                if !dep_inv[j].is_empty() {
                    step_recession[j] -= relax.thickness_m(&dep_inv[j])?;
                }
            }
            delta.retain(|_, d| *d != 0.0);
            let change = self.max_change(&delta);

            if let StepPolicy::Adaptive {
                min_ions_per_step,
                max_change,
                ..
            } = self.cfg.policy
            {
                if n > min_ions_per_step && change > max_change {
                    let scaled = (n as f64 * max_change / change).floor() as u64;
                    n = scaled.clamp(min_ions_per_step, n - 1);
                    continue;
                }
            }

            let mut clamped = 0u32;
            let new_grids = match self.apply_all(&delta) {
                Ok(g) => g,
                Err(e) => {
                    let StepPolicy::Adaptive {
                        min_ions_per_step, ..
                    } = self.cfg.policy
                    else {
                        return Err(e.into());
                    };
                    if n > min_ions_per_step {
                        n = (n / 2).max(min_ions_per_step);
                        continue;
                    }
                    // Smallest step: a removal beyond what a slab holds is
                    // capped at what it holds (as in the 1-D run).
                    if !matches!(e, DynamicError::NegativeInventory { .. }) {
                        return Err(e.into());
                    }
                    let capped: BTreeMap<(usize, usize, u8), f64> = delta
                        .iter()
                        .map(|(&(j, k, z), &d)| {
                            let have = self.grids[j]
                                .inventory(k)
                                .and_then(|v| v.iter().find(|&&(zz, _)| zz == z).map(|&(_, a)| a))
                                .unwrap_or(0.0);
                            if d < -have {
                                clamped += 1;
                                ((j, k, z), -have)
                            } else {
                                ((j, k, z), d)
                            }
                        })
                        .collect();
                    self.apply_all(&capped)?
                }
            };
            self.grids = new_grids;
            for (r, d) in self.surface.iter_mut().zip(&step_recession) {
                *r += d;
            }
            self.steps += 1;
            self.next_index += n;
            self.cumulative.add(&yields);
            let per_m = |m: &BTreeMap<(usize, u8), f64>| (0.0 + m.values().sum::<f64>()) * self.h;
            let sputtered_apm = per_m(&sputtered);
            let redeposited_apm = per_m(&deposited);
            self.sputtered_atoms_per_m += sputtered_apm;
            self.redeposited_atoms_per_m += redeposited_apm;
            if let StepPolicy::Adaptive {
                max_ions_per_step,
                max_change,
                ..
            } = self.cfg.policy
            {
                self.current_n = if change < 0.5 * max_change {
                    n.saturating_mul(2).min(max_ions_per_step)
                } else {
                    n
                };
            }
            return Ok(Some(Step2d {
                step: self.steps,
                first_index: first,
                ions: n,
                attempts,
                max_change: change,
                clamped,
                yields,
                surface_m: self.surface.clone(),
                incidence_rad: incidence,
                column_yield,
                sputtered_atoms_per_m: sputtered_apm,
                redeposited_atoms_per_m: redeposited_apm,
                implanted_atoms_per_m: implanted,
                lateral_escapes: lateral,
            }));
        }
    }

    /// Run all remaining steps, calling `on_step` after each accepted step.
    pub fn run_all(
        &mut self,
        mut on_step: impl FnMut(&Step2d, &[CompositionGrid]),
    ) -> Result<(), DynamicRunError> {
        while let Some(rec) = self.step()? {
            on_step(&rec, &self.grids);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Stack;
    use crate::ion::dynamic::Relaxation;

    const NM: f64 = 1e-9;

    fn si() -> Material {
        let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
        m.set_displacement_energy_ev(14, 15.0).unwrap();
        m.set_lattice_binding_energy_ev(14, 2.0).unwrap();
        m.set_surface_binding_energy_ev(14, 2.0).unwrap();
        m
    }

    fn target(surface: &[f64], boundary: Boundary) -> ColumnTarget {
        let stack = Stack::new(vec![(si(), 10.0 * NM)], Some(si())).unwrap();
        let g =
            CompositionGrid::from_stack(&stack, Relaxation::ideal_mixing_elemental(&[14]).unwrap())
                .unwrap();
        let grids = vec![g; surface.len()];
        ColumnTarget::new(&grids, surface, 10.0 * NM, boundary).unwrap()
    }

    #[test]
    fn flat_surface_escape_is_vertical_and_exact() {
        let t = target(&[0.0; 3], Boundary::Periodic);
        let r = t.locate([2.0 * NM, 5.0 * NM, 0.0]).unwrap();
        // Straight up, 2 nm to the surface.
        let f = t.flight(r, [2.0 * NM, 5.0 * NM, 0.0], [-1.0, 0.0, 0.0], 1.0);
        let Flight::Event(e) = f else { panic!("{f:?}") };
        assert!((e.distance - 2.0 * NM).abs() < 1e-18);
        match e.outcome {
            ExitOutcome::Escape { face, normal } => {
                assert_eq!(face, Face::Front);
                assert_eq!(normal, [-1.0, 0.0, 0.0]);
            }
            o => panic!("{o:?}"),
        }
        // Down through the layer boundary at 10 nm.
        let f = t.flight(r, [2.0 * NM, 5.0 * NM, 0.0], [1.0, 0.0, 0.0], 1.0);
        let Flight::Event(e) = f else { panic!("{f:?}") };
        assert!((e.distance - 8.0 * NM).abs() < 1e-18);
        assert!(matches!(e.outcome, ExitOutcome::Enter { region, .. } if region == r + 1));
    }

    #[test]
    fn sloped_surface_normal_follows_the_neighbours() {
        // Surface deepens by 10 nm per column: slope 1 (45 degrees).
        let t = target(&[0.0, 10.0 * NM, 20.0 * NM], Boundary::Vacuum);
        // Between centres 0 and 1 (y in [5, 15] nm) the surface is x = y - 5.
        let pos = [12.0 * NM, 10.0 * NM, 0.0];
        let r = t.locate(pos).unwrap();
        let f = t.flight(r, pos, [-1.0, 0.0, 0.0], 1.0);
        let Flight::Event(e) = f else { panic!("{f:?}") };
        assert!((e.distance - 7.0 * NM).abs() < 1e-18, "{}", e.distance);
        let ExitOutcome::Escape { normal, .. } = e.outcome else {
            panic!()
        };
        let s = std::f64::consts::FRAC_1_SQRT_2;
        assert!((normal[0] + s).abs() < 1e-12 && (normal[1] - s).abs() < 1e-12);
        // Above the surface is outside.
        assert_eq!(t.locate([2.0 * NM, 10.0 * NM, 0.0]), None);
        assert!(t.in_vacuum([2.0 * NM, 10.0 * NM, 0.0]));
    }

    #[test]
    fn periodic_wrap_and_vacuum_side() {
        let p = target(&[0.0; 2], Boundary::Periodic);
        let pos = [3.0 * NM, 19.0 * NM, 0.0];
        let r = p.locate(pos).unwrap();
        let Flight::Event(e) = p.flight(r, pos, [0.0, 1.0, 0.0], 1.0) else {
            panic!()
        };
        let ExitOutcome::Enter { region, pos } = e.outcome else {
            panic!()
        };
        assert_eq!(p.split(region).0, 0);
        assert_eq!(pos[1], 0.0);
        let v = target(&[0.0; 2], Boundary::Vacuum);
        let r = v.locate([3.0 * NM, 19.0 * NM, 0.0]).unwrap();
        let Flight::Event(e) = v.flight(r, [3.0 * NM, 19.0 * NM, 0.0], [0.0, 1.0, 0.0], 1.0) else {
            panic!()
        };
        assert!(matches!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Side,
                ..
            }
        ));
    }

    #[test]
    fn redeposition_trace_hits_the_far_wall_of_a_groove() {
        // A groove: depths 0, 40, 0 nm at the centres 5, 15, 25 nm (vacuum
        // ends). The right flank is x = 40 - 4 (y - 15) for y in [15, 25].
        // A ray from the floor, 4 nm left of the flank at 36 nm depth,
        // going up and right meets it at y = 16.2 nm: column 1.
        let t = target(&[0.0, 40.0 * NM, 0.0], Boundary::Vacuum);
        let hit = t.trace_to_surface([36.0 * NM, 15.0 * NM, 0.0], [-0.6, 0.8, 0.0]);
        assert_eq!(hit, Some(1));
        // Straight up: lost.
        assert_eq!(
            t.trace_to_surface([40.0 * NM, 15.0 * NM, 0.0], [-1.0, 0.0, 0.0]),
            None
        );
        // Going up and out over the flat rim (vacuum end): lost.
        assert_eq!(
            t.trace_to_surface([0.0, 25.0 * NM, 0.0], [-0.5, 0.866, 0.0]),
            None
        );
    }
}
