//! Iso-dose contours of a gridded scalar field, with the area and centroid of
//! each closed contour.
//!
//! The contours are found by marching squares (W. E. Lorensen and H. E. Cline,
//! "Marching cubes: a high resolution 3D surface construction algorithm",
//! Comput. Graph. 21, 163 (1987), the 2D case), with linear interpolation
//! along the cell edges. A grid point is *inside* when its value is
//! `>= threshold`.
//!
//! * Every contour is directed with the inside on its left. An outer
//!   boundary of an exposed feature is therefore counter-clockwise (positive
//!   [`Contour::signed_area_m2`]) and the boundary of a hole in a feature is
//!   clockwise (negative). The area of a feature with holes is the sum of the
//!   signed areas of its outer contour and its holes.
//! * The saddle cells (two diagonal inside corners) are resolved by the mean
//!   of the four corner values: inside connects the inside corners when the
//!   mean is `>= threshold`, and separates them otherwise.
//! * **Open contours.** A feature that touches the grid edge has a contour
//!   that leaves the grid. It is returned with `closed = false`, in the
//!   order walked, and has no area and no centroid (they are `None`): the
//!   grid does not say where the feature ends. Enlarge the grid to close it.
//! * The area of a closed contour is the shoelace area of its polygon, and the
//!   centroid is the area centroid of that polygon (Bourke's polygon
//!   formulas; textbook). The polygon is the piecewise-linear contour, so the
//!   area differs from the smooth contour's by `O(h²)` in the grid spacing.
//!
//! Contours are produced in a fixed order (by the cell and edge where each
//! starts), so the output is deterministic.

use std::collections::BTreeMap;

use crate::exposure::{DoseGrid, Grid};

/// One contour.
#[derive(Debug, Clone, PartialEq)]
pub struct Contour {
    /// Vertices `[x, y]`, m. For a closed contour the first vertex is not
    /// repeated at the end.
    pub points: Vec<[f64; 2]>,
    /// True if the contour closes on itself inside the grid.
    pub closed: bool,
    /// Signed shoelace area, m² (positive counter-clockwise, i.e. an outer
    /// boundary); `None` for an open contour.
    pub signed_area_m2: Option<f64>,
    /// Area centroid `[x, y]`, m; `None` for an open contour or one of zero
    /// area.
    pub centroid_m: Option<[f64; 2]>,
}

/// A grid edge: `(i, j, vertical)`; horizontal runs from `(i, j)` to
/// `(i+1, j)`, vertical from `(i, j)` to `(i, j+1)`.
type Edge = (usize, usize, bool);

/// The contours of `dose` at `threshold`.
pub fn contours(dose: &DoseGrid, threshold: f64) -> Vec<Contour> {
    contours_of(&dose.grid, &dose.dose_ev_per_m2, threshold)
}

/// The contours of `values` (row-major on `grid`) at `threshold`.
pub fn contours_of(grid: &Grid, values: &[f64], threshold: f64) -> Vec<Contour> {
    assert_eq!(values.len(), grid.len(), "one value per grid point");
    let (nx, ny) = (grid.nx, grid.ny);
    let v = |i: usize, j: usize| values[j * nx + i];
    // Directed segments, keyed by start edge. Marching-square contours never
    // branch, so each edge starts at most one segment.
    let mut next: BTreeMap<Edge, Edge> = BTreeMap::new();
    for j in 0..ny.saturating_sub(1) {
        for i in 0..nx.saturating_sub(1) {
            let (bl, br, tr, tl) = (v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1));
            let case = usize::from(bl >= threshold)
                | usize::from(br >= threshold) << 1
                | usize::from(tr >= threshold) << 2
                | usize::from(tl >= threshold) << 3;
            let b: Edge = (i, j, false);
            let r: Edge = (i + 1, j, true);
            let t: Edge = (i, j + 1, false);
            let l: Edge = (i, j, true);
            let centre_in = (bl + br + tr + tl) / 4.0 >= threshold;
            let segs: &[(Edge, Edge)] = match case {
                1 => &[(b, l)],
                2 => &[(r, b)],
                3 => &[(r, l)],
                4 => &[(t, r)],
                5 if centre_in => &[(b, r), (t, l)],
                5 => &[(b, l), (t, r)],
                6 => &[(t, b)],
                7 => &[(t, l)],
                8 => &[(l, t)],
                9 => &[(b, t)],
                10 if centre_in => &[(l, b), (r, t)],
                10 => &[(r, b), (l, t)],
                11 => &[(r, t)],
                12 => &[(l, r)],
                13 => &[(b, r)],
                14 => &[(l, b)],
                _ => &[],
            };
            for &(from, to) in segs {
                next.insert(from, to);
            }
        }
    }
    let point = |e: Edge| -> [f64; 2] {
        let (i, j, vertical) = e;
        let (i2, j2) = if vertical { (i, j + 1) } else { (i + 1, j) };
        let (a, b) = (v(i, j), v(i2, j2));
        let t = (threshold - a) / (b - a);
        let [xa, ya] = grid.point(i, j);
        let [xb, yb] = grid.point(i2, j2);
        [xa + t * (xb - xa), ya + t * (yb - ya)]
    };
    let has_incoming: std::collections::BTreeSet<Edge> = next.values().copied().collect();
    let mut visited: std::collections::BTreeSet<Edge> = Default::default();
    let mut out = Vec::new();
    // Open contours start at an edge nothing leads into.
    let starts: Vec<Edge> = next
        .keys()
        .copied()
        .filter(|e| !has_incoming.contains(e))
        .collect();
    for s in starts {
        let mut pts = vec![point(s)];
        visited.insert(s);
        let mut cur = s;
        while let Some(&n) = next.get(&cur) {
            pts.push(point(n));
            visited.insert(n);
            cur = n;
        }
        out.push(open(pts));
    }
    // What remains is closed loops.
    for &s in next.keys() {
        if visited.contains(&s) {
            continue;
        }
        let mut pts = Vec::new();
        let mut cur = s;
        loop {
            visited.insert(cur);
            pts.push(point(cur));
            cur = next[&cur];
            if cur == s {
                break;
            }
        }
        out.push(closed(pts));
    }
    out
}

fn open(points: Vec<[f64; 2]>) -> Contour {
    Contour {
        points,
        closed: false,
        signed_area_m2: None,
        centroid_m: None,
    }
}

fn closed(points: Vec<[f64; 2]>) -> Contour {
    let n = points.len();
    let (mut a2, mut cx, mut cy) = (0.0, 0.0, 0.0);
    for k in 0..n {
        let [x0, y0] = points[k];
        let [x1, y1] = points[(k + 1) % n];
        let cross = x0 * y1 - x1 * y0;
        a2 += cross;
        cx += (x0 + x1) * cross;
        cy += (y0 + y1) * cross;
    }
    let area = a2 / 2.0;
    let centroid = (area != 0.0).then(|| [cx / (6.0 * area), cy / (6.0 * area)]);
    Contour {
        points,
        closed: true,
        signed_area_m2: Some(area),
        centroid_m: centroid,
    }
}

impl Contour {
    /// CSV of a list of contours: `contour,closed,signed_area_m2,centroid_x_m,
    /// centroid_y_m,x_m,y_m`, one row per vertex (area and centroid on every
    /// row of the contour, empty for an open one).
    pub fn csv(contours: &[Contour]) -> String {
        let mut s =
            String::from("contour,closed,signed_area_m2,centroid_x_m,centroid_y_m,x_m,y_m\n");
        let f = |o: Option<f64>| o.map(|v| v.to_string()).unwrap_or_default();
        for (k, c) in contours.iter().enumerate() {
            for p in &c.points {
                s.push_str(&format!(
                    "{k},{},{},{},{},{},{}\n",
                    u8::from(c.closed),
                    f(c.signed_area_m2),
                    f(c.centroid_m.map(|m| m[0])),
                    f(c.centroid_m.map(|m| m[1])),
                    p[0],
                    p[1]
                ));
            }
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(n: usize) -> Grid {
        Grid {
            x0_m: -1.0,
            y0_m: -1.0,
            dx_m: 2.0 / (n - 1) as f64,
            dy_m: 2.0 / (n - 1) as f64,
            nx: n,
            ny: n,
        }
    }

    fn field(g: &Grid, f: impl Fn(f64, f64) -> f64) -> Vec<f64> {
        (0..g.len())
            .map(|k| {
                let [x, y] = g.point(k % g.nx, k / g.nx);
                f(x, y)
            })
            .collect()
    }

    #[test]
    fn circle_area_and_centroid() {
        let g = grid(201);
        let v = field(&g, |x, y| {
            0.3 - ((x - 0.1).powi(2) + (y + 0.2).powi(2)).sqrt()
        });
        let c = contours_of(&g, &v, 0.0);
        assert_eq!(c.len(), 1);
        let a = c[0].signed_area_m2.unwrap();
        assert!(c[0].closed && a > 0.0);
        assert!((a - std::f64::consts::PI * 0.09).abs() / (std::f64::consts::PI * 0.09) < 1e-3);
        let m = c[0].centroid_m.unwrap();
        assert!((m[0] - 0.1).abs() < 1e-3 && (m[1] + 0.2).abs() < 1e-3);
    }

    #[test]
    fn two_features_and_a_hole() {
        let g = grid(161);
        // Two discs, one with a hole (an annulus), via a radial field.
        let v = field(&g, |x, y| {
            let d1 = ((x + 0.5).powi(2) + y * y).sqrt();
            let d2 = ((x - 0.5).powi(2) + y * y).sqrt();
            (0.3 - d1).max(0.1 - (d2 - 0.2).abs())
        });
        let c = contours_of(&g, &v, 0.0);
        let mut areas: Vec<f64> = c.iter().map(|c| c.signed_area_m2.unwrap()).collect();
        areas.sort_by(f64::total_cmp);
        // Annulus: outer radius 0.3 (area π 0.09), inner 0.1 (hole, negative);
        // plus the 0.3 disc.
        assert_eq!(areas.len(), 3, "{areas:?}");
        assert!(areas[0] < 0.0);
        let pi = std::f64::consts::PI;
        assert!((areas[0] + pi * 0.01).abs() < 2e-3);
        assert!((areas[1] - pi * 0.09).abs() < 2e-3 && (areas[2] - pi * 0.09).abs() < 2e-3);
    }

    #[test]
    fn no_contour_when_threshold_above_or_below_all() {
        let g = grid(11);
        let v = field(&g, |x, y| x + y);
        assert!(contours_of(&g, &v, 10.0).is_empty());
        assert!(contours_of(&g, &v, -10.0).is_empty());
    }

    #[test]
    fn feature_touching_the_edge_is_open_without_area() {
        let g = grid(41);
        let v = field(&g, |x, _| 0.2 - x); // inside is x <= 0.2, reaching the edge
        let c = contours_of(&g, &v, 0.0);
        assert_eq!(c.len(), 1);
        assert!(!c[0].closed && c[0].signed_area_m2.is_none() && c[0].centroid_m.is_none());
    }

    #[test]
    fn saddle_is_resolved_by_the_cell_mean() {
        let g = Grid {
            x0_m: 0.0,
            y0_m: 0.0,
            dx_m: 1.0,
            dy_m: 1.0,
            nx: 2,
            ny: 2,
        };
        // BL and TR inside; mean 0.5 >= 0.5 joins them (one open contour per
        // outside corner), mean below separates them (one per inside corner).
        let joined = contours_of(&g, &[1.0, 0.0, 0.0, 1.0], 0.4);
        let split = contours_of(&g, &[1.0, 0.0, 0.0, 1.0], 0.6);
        assert_eq!((joined.len(), split.len()), (2, 2));
        assert_ne!(joined, split);
    }
}
