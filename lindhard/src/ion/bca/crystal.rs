//! The crystal flight model: collision partners taken from explicit lattice
//! sites instead of a random impact parameter (step 20b of the M2 crystal
//! plan).
//!
//! # Choosing the model per region
//!
//! [`Geometry`] is unchanged. A [`Bca`] is amorphous everywhere unless
//! [`Bca::with_crystal`] names a [`CrystalTarget`] (lattice and orientation)
//! and the list of **regions** (layers of a stack, voxels of a grid) it
//! fills. A particle whose current region is in that list uses the crystal
//! model; in every other region, and for a [`Bca`] without `with_crystal`, the
//! amorphous model is used unchanged (and bit for bit: the amorphous code path
//! draws the same random numbers and does the same arithmetic as before; the
//! test `tests/crystal_amorphous_identity.rs` compares complete reports).
//! Several calls give several crystals (different lattices or orientations)
//! in different regions.
//!
//! The material of a crystal region still supplies the displacement and
//! binding energies, the stopping and the surface barrier; it must contain
//! every element of the lattice and its atom density must agree with the
//! lattice's to 5 %. The crystal model needs the constant-path convention and
//! no weak collisions (the lattice replaces both). Beam direction and surface
//! stay in the lab frame of the geometry; the orientation's lab-to-crystal
//! rotation ([`Orientation::to_crystal`]) is the only thing read from it (the
//! beam is the [`Beam`](super::Beam)'s, so give the `Beam` the orientation's tilt and
//! twist, as `Orientation::beam_lab` does).
//!
//! # Lattice frame of a history
//!
//! The crystal position of a lab point `r` is `R r + u`, with `R` the
//! orientation's rotation and `u` a translation drawn uniformly from the cube
//! `[0, a)^3` once per primary history (three uniform numbers, before
//! anything else), shared by the primary and all its recoils. A beam of finite
//! width samples all positions of the unit cell, so this is the static
//! equivalent of a random entry point; it also makes the surface termination
//! random. Lattice sites that lie outside the target (in front of the surface)
//! or in a region that is not part of this crystal are not partners.
//!
//! # Collision sequence (after DISPLATH, Tier A)
//!
//! The ion moves in straight segments. For a segment from position `r0` along
//! the unit direction `d` the partners are the lattice sites within `p_max`
//! of the line, ordered by the path distance `s = (r - r0) . d` of their
//! point of closest approach ([`crate::ion::crystal::LatticeSearch`]). The
//! scheme is the one of DISPLATH (<https://github.com/permissionx/DISPLATH>,
//! MIT, commit `7f461141c518305e6a5de8ce7e37dd314cd20ccb`, `src/dynamics.jl`:
//! `GetTargetsFromNeighbor` and `Collision!`, and `SimultaneousCriteria` in
//! `src/geometry.jl`), re-implemented here on [`LatticeSearch`] with
//! the lattice-site list of this crate, and with the project's own scattering
//! tables, kinematics and energy bookkeeping. Attribution:
//! `THIRD_PARTY_LICENSES.md`.
//!
//! * **Free path.** The ion flies straight to the closest-approach point of
//!   the nearest partner (smallest `s`, ties broken by the search's documented
//!   order). There is no free-path length parameter: the path is the crystal's.
//!   The nonlocal electronic loss of the model chosen in the [`BcaConfig`](super::BcaConfig) is
//!   taken over that length at the energy at the start of the segment, as in
//!   the amorphous model; impact-parameter-dependent stopping is a later
//!   step. When no partner lies within `search_length` the ion flies that
//!   distance and searches again (a channel). The geometry is asked about
//!   every segment, so surfaces and interfaces cut it exactly like an
//!   amorphous flight.
//! * **Simultaneous collisions.** Further sites join the nearest one if, for
//!   every partner `t` already accepted, `Δ = s - s_t` (≥ 0) satisfies
//!   `Δ <= q_max`, `p_t² + Δ² <= p_max²` and `p² + Δ² <= p_max²`
//!   (DISPLATH's `SimultaneousCriteria`: sites at similar path distance and
//!   inside a sphere of radius `p_max` around the collision point). DISPLATH
//!   sets `q_max` to the sum of the two atomic radii (about the
//!   nearest-neighbour distance); here the default is a small fraction of it,
//!   see the next section.
//! * **Collision.** Each partner `i` is scattered from the incoming state at
//!   its own impact parameter `p_i` (the distance from the site to the line)
//!   with the scattering table: transfer `T_i = γ E sin²(Θ_i/2)`, the recoil
//!   leaving along `sin(Θ_i/2) d + cos(Θ_i/2) n_i`, with `n_i` the unit vector
//!   from the line to the site (the half-angle form of
//!   [`kinematics::lab_angles`]: the recoil angle is `(π - Θ)/2`). The
//!   projectile takes the momentum that is left, `P = P0 - Σ sqrt(2 m_i T_i)
//!   r_i`, and the energies are scaled by the common factor
//!   `λ = (E - Σ Q_i) / (|P|²/2m + Σ T_i)` so that energy is conserved with
//!   the local electronic losses `Q_i` (DISPLATH's `Collision_!`). For one
//!   partner and no local loss this is the ordinary binary collision
//!   (`λ = 1`). The projectile stays at the closest-approach point of the
//!   nearest partner; the time-delay offsets of Robinson and Torrens (1974)
//!   are not applied (DISPLATH does not apply them either).
//! * **Recoils** start at their lattice site, with the usual displacement
//!   criterion ([`BcaConfig::follow_recoils`](super::BcaConfig::follow_recoils)) and follow the same model. The
//!   sites that were just hit are excluded from the next search (so is the
//!   site the recoil starts on): the ion cannot collide twice in a row with
//!   the same atom.
//!
//! # Choice of the search parameters
//!
//! Measured on this engine (B 5 keV into Si, 4000 ions, ZBL, Lindhard-Scharff,
//! 7 degrees tilt, 22 degrees twist, static lattice; ratios of Rp to the
//! amorphous result), with `nn` the nearest-neighbour distance (2.35 Å):
//!
//! * With `q_max = nn` (DISPLATH's rule) the result depends strongly on
//!   `p_max`: `p_max` = 1.53 Å gives Rp ratio 1.11, but `p_max` = 2.3 Å gives
//!   4.4. The merged collisions of distant sites cancel the kicks of a
//!   channel's walls, so almost every ion channels, at 20 and 30 degrees tilt
//!   as well. That is an artefact of merging over a sphere as large as the
//!   channel, not physics.
//! * With `q_max <= 0.1 nn` the result is independent of `q_max` and, for
//!   `p_max >= 1.5 * 1.53 Å = nn`, of `p_max` within the statistics (Rp ratio
//!   1.08-1.09 for `p_max` = 2.3, 3.1 Å and `q_max` = 0.01-0.1 nn; at
//!   `q_max` = 0.2 nn the artefact starts).
//!
//! So the defaults are `p_max = nn` and `q_max = 0.05 nn`: converged in both,
//! and the cost grows with `p_max`. These are properties of this engine's
//! implementation of the criterion, checked on silicon only; other lattices
//! need the same convergence check.
//!
//! # Measured checks (static lattice, 4000 ions, seed 1, Si, ZBL, Lindhard-Scharff)
//!
//! `tests/crystal_bca.rs` (the statistical ones are `#[ignore]`d; release
//! build). Ratios are crystal over amorphous:
//!
//! | Case | Rp | dRp | Notes |
//! |---|---|---|---|
//! | B 5 keV, 30° tilt, 17° twist | 0.96 | 1.10 | within 10 % / 15 % |
//! | As 30 keV, 30° tilt, 17° twist | 0.91 | 1.12 | within 10 % / 15 % |
//! | B 5 keV, 7° tilt, 22° twist | 1.08 | **1.27** | median 24.2 vs 23.8 nm, 99th percentile 77.8 vs 56.8 nm |
//! | As 30 keV, 7° tilt, 22° twist | 1.04 | **1.95** | median 22.2 vs 25.2 nm, 99th percentile 114.2 vs 50.8 nm |
//! | B 5 keV along <110> | **6.1** (121 vs 20 nm) | 3.1 | 91 % of ions deeper than twice the amorphous Rp (amorphous 5 %) |
//!
//! **Gap: the dRp bound of the 7°/22° check is not met.** Rp agrees within
//! 10 % for both ions, but a static lattice keeps a channeling tail at this
//! orientation (the beam is only 2.6° to 2.7° from a {100} and a {110}
//! plane), which widens dRp by 27 % and 95 %; the bulk of the profile (median,
//! 90th percentile) matches the amorphous one. Thermal vibration (step 21b)
//! is the intended remedy and has not been tried here. At 30°/17°, far from
//! every low-index axis and plane, both bounds hold. The 7°/22° test records
//! the gap instead of asserting the unmet bound.
//!
//! Electronic loss in a channel is the amorphous average (impact-parameter-
//! dependent stopping is Phase 2), which makes the channeled ranges an upper
//! bound on this effect.
//!
//! # What is not here
//!
//! The lattice is **static and perfect** (thermal vibration is step 21b) and
//! stays perfect: a displaced atom leaves no vacancy, so a later particle can
//! still collide with its site (dynamic damage, dechanneling and
//! amorphization are later steps). There is no electronic-stopping dependence
//! on the impact parameter (Phase 2). Boundaries of a periodic voxel axis do
//! not continue the lattice. The equation-level treatment of Robinson and
//! Torrens, Phys. Rev. B 9, 5008 (1974) (MARLOWE, Tier C) was not available
//! to read when this was written; the choices above are DISPLATH's, not
//! quotations of that paper.
//!
//! [`Geometry`]: crate::geometry::Geometry
//! [`LatticeSearch`]: crate::ion::crystal::LatticeSearch

use crate::geometry::Flight;
use crate::ion::crystal::{Lattice, LatticeSearch, Orientation};
use crate::ion::scattering::closest_approach;
use crate::ion::stopping::Ion;
use crate::units::J_PER_EV;

use super::{
    kinematics, Bca, BcaError, BcaTally, ElectronicChannel, EnergyBudget, MeanFreePath, Particle,
    Scratch,
};
pub(super) use crate::ion::crystal::Candidate;

/// A lattice site: its conventional cell and index in the cell.
pub(super) type SiteId = ([i64; 3], u8);

/// A site closer than this fraction of the lattice constant to the particle
/// is the particle's own starting site, not a partner.
const SELF_REL: f64 = 1e-6;

/// Default `q_max`, as a fraction of the nearest-neighbour distance (module
/// docs, "Choice of the search parameters").
pub const DEFAULT_Q_MAX_FRACTION: f64 = 0.05;

/// Largest relative density mismatch between lattice and material.
const DENSITY_TOLERANCE: f64 = 0.05;

/// A crystal for [`Bca::with_crystal`]: the lattice, the orientation of the
/// target in the lab frame and the three search parameters.
#[derive(Debug, Clone)]
pub struct CrystalTarget {
    /// The lattice (cubic; diamond or zincblende).
    pub lattice: Lattice,
    /// Wafer cut and angles. Only its lab-to-crystal rotation is used.
    pub orientation: Orientation,
    /// Largest impact parameter of a partner, m. Default `None`: the
    /// nearest-neighbour distance of the lattice. Steering by the atoms of a
    /// channel wall needs partners out to about that distance; smaller radii
    /// (DISPLATH's examples use `a/3` and the amorphous constant-path radius
    /// `(π N^(2/3))^(-1/2)`, `a/2/π^0.5` for silicon) leave the channeled
    /// fraction too small. See the module docs, "Choice of the search
    /// parameters".
    pub p_max_m: Option<f64>,
    /// Largest path distance between simultaneous partners, m. Default
    /// `None`: [`DEFAULT_Q_MAX_FRACTION`] of the nearest-neighbour distance.
    pub q_max_m: Option<f64>,
    /// Length of one search segment, m. Default `None`: the lattice
    /// constant. Affects speed only (a segment with no partner is followed by
    /// another), never the result.
    pub search_length_m: Option<f64>,
}

impl CrystalTarget {
    /// A crystal with the default search parameters.
    pub fn new(lattice: Lattice, orientation: Orientation) -> Self {
        Self {
            lattice,
            orientation,
            p_max_m: None,
            q_max_m: None,
            search_length_m: None,
        }
    }
}

/// Validated, precomputed crystal data.
#[derive(Debug, Clone)]
pub(super) struct CrystalData {
    search: LatticeSearch,
    orientation: Orientation,
    a: f64,
    p_max: f64,
    q_max: f64,
    search_len: f64,
}

impl CrystalData {
    pub(super) fn lattice_constant(&self) -> f64 {
        self.a
    }

    /// Crystal-frame position of a lab point for the history translation
    /// `shift`.
    fn to_crystal(&self, lab: [f64; 3], shift: [f64; 3]) -> [f64; 3] {
        let c = self.orientation.to_crystal(lab);
        [c[0] + shift[0], c[1] + shift[1], c[2] + shift[2]]
    }

    /// Lab position of a crystal-frame point.
    fn to_lab(&self, c: [f64; 3], shift: [f64; 3]) -> [f64; 3] {
        self.orientation
            .to_lab([c[0] - shift[0], c[1] - shift[1], c[2] - shift[2]])
    }
}

/// One partner of the current collision step.
#[derive(Debug, Clone, Copy)]
pub(super) struct Partner {
    id: SiteId,
    /// Path distance of the closest approach, m.
    s: f64,
    /// Impact parameter, m.
    b: f64,
    /// Unit vector (lab) from the line to the site.
    n: [f64; 3],
    /// Site position (lab), m.
    site: [f64; 3],
    /// Region of the site and element index in its material.
    region: usize,
    j: usize,
    /// Transfer, recoil direction and local loss of the pair collision.
    t: f64,
    rdir: [f64; 3],
    q: f64,
}

/// Nearest-neighbour distance of `lattice`, m.
fn nearest_neighbour(lattice: &Lattice) -> f64 {
    let a = lattice.lattice_constant();
    let sites = lattice.conventional_cell_sites();
    let mut best = f64::INFINITY;
    for (_, f) in &sites {
        for (_, g) in &sites {
            for i in -1..=1 {
                for j in -1..=1 {
                    for k in -1..=1 {
                        let dx = a * (g[0] + f64::from(i) - f[0]);
                        let dy = a * (g[1] + f64::from(j) - f[1]);
                        let dz = a * (g[2] + f64::from(k) - f[2]);
                        let d = (dx * dx + dy * dy + dz * dz).sqrt();
                        if d > 1e-6 * a {
                            best = best.min(d);
                        }
                    }
                }
            }
        }
    }
    best
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl<'a> Bca<'a> {
    /// Use the crystal flight model in `regions` (see the module docs). The
    /// regions are indices into the geometry's regions: layer indices of a
    /// stack, flat voxel indices of a voxel grid.
    ///
    /// # Errors
    /// [`BcaError::InvalidConfig`] if the configuration is not the constant
    /// free path without weak collisions, a region does not exist or already
    /// belongs to a crystal, a search parameter is not finite and positive, a
    /// region's material lacks a lattice element, or its atom density differs
    /// from the lattice's by more than 5 %.
    pub fn with_crystal(
        mut self,
        target: CrystalTarget,
        regions: &[usize],
    ) -> Result<Self, BcaError> {
        let bad = |m: String| BcaError::InvalidConfig(m);
        if !matches!(self.config.mean_free_path, MeanFreePath::Constant)
            || self.config.weak_collisions > 0
        {
            return Err(bad(
                "crystal regions need the constant free path and no weak collisions".into(),
            ));
        }
        if regions.is_empty() {
            return Err(bad("a crystal needs at least one region".into()));
        }
        let lat = &target.lattice;
        let a = lat.lattice_constant();
        let n_lat = lat.atom_number_density();
        let zs: Vec<u8> = {
            let mut v: Vec<u8> = lat.conventional_cell_sites().iter().map(|s| s.0).collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        if self.region_crystal.is_empty() {
            self.region_crystal = vec![None; self.geometry.n_regions()];
        }
        let ci = self.crystals.len();
        for &r in regions {
            if r >= self.geometry.n_regions() {
                return Err(bad(format!("crystal region {r} does not exist")));
            }
            if self.region_crystal[r].is_some() {
                return Err(bad(format!("region {r} already belongs to a crystal")));
            }
            let lay = self.lay(r);
            for &z in &zs {
                if !lay.elems.iter().any(|e| e.z == z) {
                    return Err(bad(format!(
                        "region {r}: the material has no element Z={z} of the lattice"
                    )));
                }
            }
            if (lay.n / n_lat - 1.0).abs() > DENSITY_TOLERANCE {
                return Err(bad(format!(
                    "region {r}: material atom density {:.4e} m^-3 differs from the \
                     lattice's {n_lat:.4e} m^-3 by more than 5 %",
                    lay.n
                )));
            }
        }
        let nn = nearest_neighbour(lat);
        let p_max = target.p_max_m.unwrap_or(nn);
        let q_max = target.q_max_m.unwrap_or(DEFAULT_Q_MAX_FRACTION * nn);
        let search_len = target.search_length_m.unwrap_or(a);
        for (name, v) in [
            ("p_max_m", p_max),
            ("q_max_m", q_max),
            ("search_length_m", search_len),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(bad(format!("{name} = {v} must be finite and positive")));
            }
        }
        for &r in regions {
            self.region_crystal[r] = Some(ci);
        }
        self.crystals.push(CrystalData {
            search: LatticeSearch::new(lat),
            orientation: target.orientation,
            a,
            p_max,
            q_max,
            search_len,
        });
        Ok(self)
    }

    /// The crystal of `region`, if it is a crystal region.
    #[inline]
    pub(super) fn crystal_at(&self, region: usize) -> Option<usize> {
        self.region_crystal.get(region).copied().flatten()
    }

    /// Region and element index of a lattice site that may be a partner:
    /// inside the target, in a region of crystal `ci`.
    fn site_region(&self, ci: usize, site_lab: [f64; 3], z: u8) -> Option<(usize, usize)> {
        let r = self.geometry.locate(site_lab)?;
        if self.crystal_at(r) != Some(ci) {
            return None;
        }
        let j = self.lay(r).elems.iter().position(|e| e.z == z)?;
        Some((r, j))
    }

    /// One step of a particle in crystal `ci`: a flight to the next
    /// collision (or to a surface or interface, or one search segment) and
    /// the collision there. Returns `true` if the particle is finished.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn crystal_step<T: BcaTally>(
        &self,
        ci: usize,
        shift: [f64; 3],
        p: &mut Particle,
        pending: &mut Vec<Particle>,
        budget: &mut EnergyBudget,
        tally: &mut T,
        scratch: &mut Scratch,
    ) -> Result<bool, crate::ion::stopping::StoppingError> {
        let cr = &self.crystals[ci];
        let cutoff = self.species[p.species].cutoff_ev;
        let Scratch {
            cands,
            last,
            targets,
            ..
        } = scratch;
        let o = cr.to_crystal(p.pos, shift);
        let d = cr.orientation.to_crystal(p.dir);
        cr.search
            .search_into(o, d, cr.p_max, cr.search_len + cr.q_max, cands)
            .expect("particle state is finite");

        // The nearest valid partner within one search segment.
        let self_r2 = (SELF_REL * cr.a) * (SELF_REL * cr.a);
        let usable = |c: &Candidate| -> Option<(usize, usize)> {
            if c.s * c.s + c.p * c.p < self_r2 || last.contains(&(c.cell, c.basis)) {
                return None;
            }
            self.site_region(ci, cr.to_lab(c.position, shift), c.z)
        };
        let hit = cands
            .iter()
            .position(|c| c.s <= cr.search_len && usable(c).is_some());
        let limit = hit.map_or(cr.search_len, |i| cands[i].s);

        match self.geometry.flight(p.layer, p.pos, p.dir, limit) {
            Flight::Event(ex) => {
                let from = p.pos;
                p.pos = ex.at;
                self.electronic_nonlocal(p, from, ex.distance, budget, tally)?;
                if p.energy_ev < cutoff {
                    return Ok(false);
                }
                return Ok(self.apply_event(p, ex.outcome, budget, tally));
            }
            Flight::Clear { region } => {
                let from = p.pos;
                Self::advance(p, limit);
                p.layer = region;
                self.electronic_nonlocal(p, from, limit, budget, tally)?;
            }
        }
        let Some(hit) = hit else {
            return Ok(false);
        };
        if p.energy_ev < cutoff || self.crystal_at(p.layer) != Some(ci) {
            return Ok(false);
        }

        // Partners: the nearest, then every later site that is simultaneous
        // with all partners accepted so far.
        targets.clear();
        let s0 = cands[hit].s;
        let sp2 = cr.p_max * cr.p_max;
        for c in &cands[hit..] {
            if c.s - s0 > cr.q_max {
                break;
            }
            let Some((region, j)) = usable(c) else {
                continue;
            };
            let simultaneous = targets.iter().all(|t: &Partner| {
                let dl = c.s - t.s;
                dl <= cr.q_max && t.b * t.b + dl * dl <= sp2 && c.p * c.p + dl * dl <= sp2
            });
            if !simultaneous {
                continue;
            }
            // Unit vector from the line to the site, in the lab frame.
            let q = [
                c.position[0] - o[0] - c.s * d[0],
                c.position[1] - o[1] - c.s * d[1],
                c.position[2] - o[2] - c.s * d[2],
            ];
            let qn = dot(q, q).sqrt();
            let n_c = if qn > 1e-9 * cr.a {
                [q[0] / qn, q[1] / qn, q[2] / qn]
            } else {
                // On the line: any perpendicular does (the same rule as the
                // azimuth of an amorphous collision).
                kinematics::rotate_sc(d, (1.0, 0.0), (1.0, 0.0))
            };
            targets.push(Partner {
                id: (c.cell, c.basis),
                s: c.s,
                b: c.p,
                n: kinematics::normalize(cr.orientation.to_lab(n_c)),
                site: cr.to_lab(c.position, shift),
                region,
                j,
                t: 0.0,
                rdir: [0.0; 3],
                q: 0.0,
            });
        }
        last.clear();
        last.extend(targets.iter().map(|t| t.id));
        self.crystal_collide(p, targets, pending, budget, tally)?;
        Ok(false)
    }

    /// The simultaneous collision of `p` with `targets` (module docs).
    fn crystal_collide<T: BcaTally>(
        &self,
        p: &mut Particle,
        targets: &mut [Partner],
        pending: &mut Vec<Particle>,
        budget: &mut EnergyBudget,
        tally: &mut T,
    ) -> Result<(), crate::ion::stopping::StoppingError> {
        let ns = self.species.len();
        let e0 = p.energy_ev;
        let d0 = p.dir;
        let m_p = p.mass_amu;
        let beta_min = self.table.spec().beta_min;
        let p0 = (2.0 * m_p * e0).sqrt();
        let mut mom = [p0 * d0[0], p0 * d0[1], p0 * d0[2]];
        let (mut sum_t, mut sum_q) = (0.0, 0.0);
        for t in targets.iter_mut() {
            let elem = &self.lay(t.region).elems[t.j];
            let pair = self.pairs[p.species * ns + elem.species];
            let eps = e0 * pair.eps_per_ev;
            let beta = (t.b / pair.a).max(beta_min);
            let tan_half = self.half_angle_tan(eps, beta);
            let c = 1.0 / (1.0 + tan_half * tan_half).sqrt();
            let s = tan_half * c;
            t.t = (pair.gamma * e0 * s * s).min(e0);
            // Recoil leaves at (pi - Theta)/2 from the incoming direction,
            // towards the site: cos = sin(Theta/2), sin = cos(Theta/2).
            t.rdir = [
                s * d0[0] + c * t.n[0],
                s * d0[1] + c * t.n[1],
                s * d0[2] + c * t.n[2],
            ];
            if self.config.electronic == super::ElectronicLoss::EquipartitionLsOr {
                let r_min = closest_approach(self.screening, eps, beta) * pair.a;
                let ion: &Ion = &self.species[p.species].ion;
                let q = self.mix.0.local_loss(ion, elem.z, e0, r_min)? / J_PER_EV;
                t.q = q.min(e0 - sum_q).max(0.0);
                sum_q += t.q;
            }
            let pt = (2.0 * self.species[elem.species].ion.mass_amu() * t.t).sqrt();
            for (m, r) in mom.iter_mut().zip(t.rdir) {
                *m -= pt * r;
            }
            sum_t += t.t;
        }
        // Energy conservation with the local losses: scale the transfers and
        // the projectile's remaining energy by a common factor.
        let pm2 = dot(mom, mom);
        let e_p = pm2 / (2.0 * m_p);
        let need = (e0 - sum_q).max(0.0);
        let lambda = if e_p + sum_t > 0.0 {
            need / (e_p + sum_t)
        } else {
            1.0
        };
        if pm2 > 0.0 {
            p.dir = kinematics::normalize(mom);
        }
        p.energy_ev = lambda * e_p;
        for t in targets.iter() {
            if t.q > 0.0 {
                budget.electronic_local += t.q;
                tally.electronic(p, p.pos, ElectronicChannel::Local, t.q);
            }
        }
        for t in targets.iter() {
            let rdir = kinematics::normalize(t.rdir);
            self.emit_recoil(
                p.generation,
                t.site,
                t.region,
                t.j,
                lambda * t.t,
                || rdir,
                pending,
                budget,
                tally,
            );
        }
        Ok(())
    }
}
