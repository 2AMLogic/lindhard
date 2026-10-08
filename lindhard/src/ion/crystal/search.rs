//! Lattice neighbour search: the sites near a straight path segment, found
//! by walking unit cells (step 20a of the M2 crystal plan).
//!
//! Given an ion position, a direction, a search radius `p_max` and a segment
//! length `L`, [`LatticeSearch::search`] returns every lattice site whose
//! impact parameter against the segment is at most `p_max`, together with
//! its distance along the path, ordered along the path. This is the
//! geometric first step of a binary-collision search on explicit atom
//! positions (Robinson and Torrens, Phys. Rev. B 9, 5008 (1974),
//! doi:10.1103/PhysRevB.9.5008, describe the method; only the paper is used
//! here, no code). Nothing in the binary-collision engine reads this module
//! yet.
//!
//! # Geometry and units
//!
//! Everything is in the **crystal frame, metres** (see the `crystal` module
//! docs, "Frames"). The lattice is the infinite periodic crystal: cell
//! `(i, j, k)` of the conventional cubic cell (edge `a`) holds the atoms
//! `a (cell + f)` for each fractional position `f` in
//! [`Lattice::conventional_cell_sites`]. There is no global atom list: the
//! per-cell site list is built once in [`LatticeSearch::new`], and memory use
//! does not depend on how far from the origin the path lies.
//!
//! For a unit direction `d`, a path start `r0` and a site at `r`:
//!
//! ```text
//! s  = (r - r0) . d             distance along the path
//! p² = | (r - r0) - s d |²      squared impact parameter
//! ```
//!
//! [`path_metrics`] is the single implementation of these two expressions;
//! the cell walk and any reference scan must both go through it so that they
//! agree bit for bit.
//!
//! # Boundary rule
//!
//! A site is a candidate iff `0 <= s <= L` **and** `p² <= p_max²`
//! (evaluated in floating point exactly as written, on the values returned by
//! [`path_metrics`]). Both ends are closed: a site exactly at the start of
//! the segment (`s = 0`, e.g. the atom the ion starts on), exactly at its end
//! (`s = L`) or exactly at the radius (`p² = p_max²`) is included. Sites
//! behind the start (`s < 0`) are never returned. An empty result is not an
//! error.
//!
//! # Ordering
//!
//! Results are sorted by `s` ascending ([`f64::total_cmp`]); equal `s`
//! (many sites at once for rays along low-index directions) is broken by the
//! cell index `(i, j, k)` lexicographically and then the basis index. The
//! order is a pure function of the inputs, so it is independent of thread
//! count; a [`LatticeSearch`] is immutable after construction and holds no
//! shared state.
//!
//! # Cell walk
//!
//! The cell axis most parallel to the path is swept slab by slab. For each
//! slab the path parameter range that can reach it (a cylinder point lies
//! within `p_max` of the axis in every coordinate) is clipped to `[0, L]`,
//! which bounds the other two cell indices. Each cell in that box is skipped
//! if its centre cannot be within `p_max` of the segment (the metrics are
//! 1-Lipschitz in position, so a half-diagonal margin is conservative). All
//! bounds are padded by a relative `1e-6 a` so that rounding never drops a
//! site that the exact comparison would keep.

use super::{CrystalError, Lattice};

/// Relative padding (in cell edges) on every conservative bound.
const PAD_REL: f64 = 1e-6;

/// One lattice site found near the path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Candidate {
    /// Atomic number of the atom on the site.
    pub z: u8,
    /// Integer index `(i, j, k)` of the conventional cell holding the site.
    pub cell: [i64; 3],
    /// Index of the site within the cell, into
    /// [`Lattice::conventional_cell_sites`].
    pub basis: u8,
    /// Site position, crystal frame, metres.
    pub position: [f64; 3],
    /// Distance along the path from its start, metres (`0 <= s <= L`).
    pub s: f64,
    /// Impact parameter, metres (`p <= p_max`).
    pub p: f64,
}

/// Distance along the path and squared impact parameter of `site` for a path
/// starting at `origin` with **unit** direction `dir`:
/// `s = (site - origin) . dir`, `p² = |(site - origin) - s dir|²`.
///
/// This is the one formula used for the boundary comparisons; see the module
/// docs.
pub fn path_metrics(origin: [f64; 3], dir: [f64; 3], site: [f64; 3]) -> (f64, f64) {
    let dr = [
        site[0] - origin[0],
        site[1] - origin[1],
        site[2] - origin[2],
    ];
    let s = dr[0] * dir[0] + dr[1] * dir[1] + dr[2] * dir[2];
    let q = [dr[0] - s * dir[0], dr[1] - s * dir[1], dr[2] - s * dir[2]];
    (s, q[0] * q[0] + q[1] * q[1] + q[2] * q[2])
}

/// Position of basis site `f` (fractions of the cell edge) in cell `cell`.
/// Shared by the search and by reference scans.
pub fn site_position(a: f64, cell: [i64; 3], f: [f64; 3]) -> [f64; 3] {
    [
        a * (cell[0] as f64 + f[0]),
        a * (cell[1] as f64 + f[1]),
        a * (cell[2] as f64 + f[2]),
    ]
}

/// Candidate-site search through an infinite periodic [`Lattice`].
#[derive(Debug, Clone)]
pub struct LatticeSearch {
    a: f64,
    sites: Vec<(u8, [f64; 3])>,
}

fn bad(name: &'static str, why: &'static str) -> CrystalError {
    CrystalError::InvalidSearch { name, why }
}

impl LatticeSearch {
    /// Prepare the per-cell site list of `lattice`.
    pub fn new(lattice: &Lattice) -> Self {
        Self {
            a: lattice.lattice_constant(),
            sites: lattice.conventional_cell_sites(),
        }
    }

    /// Cubic lattice constant, metres.
    pub fn lattice_constant(&self) -> f64 {
        self.a
    }

    /// Candidate sites within `p_max` of the segment of length `length`
    /// starting at `origin` along `direction` (any non-zero length; it is
    /// normalised), ordered along the path. See the module docs for the
    /// boundary and ordering rules.
    ///
    /// # Errors
    /// [`CrystalError::InvalidSearch`] for a non-finite origin, a zero or
    /// non-finite direction, a `p_max` that is not finite and positive, or a
    /// `length` that is not finite and non-negative.
    pub fn search(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        p_max: f64,
        length: f64,
    ) -> Result<Vec<Candidate>, CrystalError> {
        let mut out = Vec::new();
        self.search_into(origin, direction, p_max, length, &mut out)?;
        Ok(out)
    }

    /// As [`search`](Self::search), but clears and fills `out`, reusing its
    /// allocation.
    pub fn search_into(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        p_max: f64,
        length: f64,
        out: &mut Vec<Candidate>,
    ) -> Result<(), CrystalError> {
        if !origin.iter().all(|x| x.is_finite()) {
            return Err(bad("origin", "must be finite"));
        }
        let norm2 = direction.iter().map(|x| x * x).sum::<f64>();
        if !direction.iter().all(|x| x.is_finite()) || !norm2.is_finite() || norm2 == 0.0 {
            return Err(bad("direction", "must be finite and non-zero"));
        }
        if !(p_max.is_finite() && p_max > 0.0) {
            return Err(bad("p_max", "must be finite and positive"));
        }
        if !(length.is_finite() && length >= 0.0) {
            return Err(bad("length", "must be finite and non-negative"));
        }
        out.clear();
        let n = norm2.sqrt();
        let d = [direction[0] / n, direction[1] / n, direction[2] / n];
        let a = self.a;
        let pad = p_max + PAD_REL * a;
        let spad = PAD_REL * a;
        // Half-diagonal of a cell, padded: bound on |site - cell centre|.
        let h = a * 0.75_f64.sqrt() * (1.0 + PAD_REL) + spad;

        // Sweep along the axis most parallel to the path.
        let ax = (0..3)
            .max_by(|&i, &j| d[i].abs().total_cmp(&d[j].abs()).then(j.cmp(&i)))
            .unwrap_or(0);
        let (b, c) = match ax {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        };
        let x0 = origin[ax];
        let x1 = x0 + length * d[ax];
        let (xlo, xhi) = (x0.min(x1) - pad, x0.max(x1) + pad);
        let i_lo = (xlo / a).floor() as i64;
        let i_hi = (xhi / a).floor() as i64;

        for ia in i_lo..=i_hi {
            // Path parameters that can reach this slab (|offset| <= p_max in
            // every coordinate), clipped to the segment.
            let ta = ((ia as f64) * a - pad - x0) / d[ax];
            let tb = (((ia + 1) as f64) * a + pad - x0) / d[ax];
            let (t_lo, t_hi) = if ta <= tb { (ta, tb) } else { (tb, ta) };
            let t_lo = t_lo.max(-spad);
            let t_hi = t_hi.min(length + spad);
            if t_lo > t_hi {
                continue;
            }
            let range = |k: usize| {
                let u = origin[k] + t_lo * d[k];
                let v = origin[k] + t_hi * d[k];
                (
                    ((u.min(v) - pad) / a).floor() as i64,
                    ((u.max(v) + pad) / a).floor() as i64,
                )
            };
            let (jb_lo, jb_hi) = range(b);
            let (jc_lo, jc_hi) = range(c);
            for ib in jb_lo..=jb_hi {
                for ic in jc_lo..=jc_hi {
                    let mut cell = [0i64; 3];
                    cell[ax] = ia;
                    cell[b] = ib;
                    cell[c] = ic;
                    let centre = site_position(a, cell, [0.5; 3]);
                    let (sc, p2c) = path_metrics(origin, d, centre);
                    if sc < -h || sc > length + h {
                        continue;
                    }
                    let lim = p_max + h;
                    if p2c > lim * lim {
                        continue;
                    }
                    for (bi, &(z, f)) in self.sites.iter().enumerate() {
                        let pos = site_position(a, cell, f);
                        let (s, p2) = path_metrics(origin, d, pos);
                        if s >= 0.0 && s <= length && p2 <= p_max * p_max {
                            out.push(Candidate {
                                z,
                                cell,
                                basis: bi as u8,
                                position: pos,
                                s,
                                p: p2.sqrt(),
                            });
                        }
                    }
                }
            }
        }
        out.sort_by(compare);
        Ok(())
    }
}

/// The documented order: `s`, then cell index, then basis index.
pub fn compare(x: &Candidate, y: &Candidate) -> std::cmp::Ordering {
    x.s.total_cmp(&y.s)
        .then(x.cell.cmp(&y.cell))
        .then(x.basis.cmp(&y.basis))
}
