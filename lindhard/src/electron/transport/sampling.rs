//! The stateless lookups of the electron transport: the rate of a table at an
//! energy, the inverse CDF at a probability, the channel chosen from the
//! rates, and the deflection of a direction by a polar and an azimuth angle.
//! None of them draws from the random stream; the engine
//! ([`Transport`](super::Transport)) draws and passes the uniform in.

use crate::electron::data::CrossSectionTable;
use crate::vec3::{cross, normalize};

/// The channel for `x = u · total`: the first index whose cumulative rate
/// exceeds `x`, skipping zero rates; if rounding leaves `x` at or past the
/// sum, the last channel with a positive rate. With only the first two rates
/// positive this is `x < el` → elastic, else inelastic, as before the
/// insulator and inner-shell channels existed. Only called when the total is
/// positive.
pub(super) fn choose(x: f64, rates: impl IntoIterator<Item = f64>) -> usize {
    let mut cum = 0.0;
    let mut last = 0;
    for (k, r) in rates.into_iter().enumerate() {
        if r > 0.0 {
            cum += r;
            last = k;
            if x < cum {
                return k;
            }
        }
    }
    last
}

/// Where an energy sits on a table's grid: lower row and weight of the upper.
#[derive(Clone, Copy)]
pub(super) struct GridPos {
    i: usize,
    f: f64,
}

/// Inverse mean free path at `e` (linear in energy, constant beyond the ends)
/// and the grid position.
pub(super) fn rate(t: &CrossSectionTable, e: f64) -> (f64, GridPos) {
    let g = t.energy_ev();
    let r = t.inverse_mfp_per_m();
    let n = g.len();
    let pos = if e <= g[0] {
        GridPos { i: 0, f: 0.0 }
    } else if e >= g[n - 1] {
        GridPos { i: n - 2, f: 1.0 }
    } else {
        let i = g.partition_point(|&x| x <= e) - 1;
        GridPos {
            i,
            f: (e - g[i]) / (g[i + 1] - g[i]),
        }
    };
    let v = if pos.f == 0.0 {
        r[pos.i]
    } else if pos.f == 1.0 {
        r[pos.i + 1]
    } else {
        r[pos.i] + pos.f * (r[pos.i + 1] - r[pos.i])
    };
    (v, pos)
}

/// The inverse CDF at probability `u`, interpolated between the bracketing
/// rows. Only called where the rate is positive, so one of the rows exists.
pub(super) fn sample(t: &CrossSectionTable, at: GridPos, u: f64) -> f64 {
    let row = |i: usize| t.inverse_cdf(i, u).ok();
    match (at.f, row(at.i), row(at.i + 1)) {
        (f, Some(a), Some(b)) if f > 0.0 && f < 1.0 => a + f * (b - a),
        (f, Some(a), _) if f < 1.0 => a,
        (_, _, Some(b)) => b,
        (_, Some(a), None) => a,
        (_, None, None) => unreachable!("a positive rate implies a stored row"),
    }
}

/// Rotate the unit vector `d` by polar angle `theta` about itself and azimuth
/// `phi` about `d` (measured from an arbitrary but deterministic reference
/// axis), returning a unit vector.
pub(super) fn deflect(d: [f64; 3], theta: f64, phi: f64) -> [f64; 3] {
    let (st, ct) = theta.sin_cos();
    deflect_cs(d, ct, st, phi)
}

/// [`deflect`] with the polar angle given by its cosine `ct` and sine `st`.
pub(crate) fn deflect_cs(d: [f64; 3], ct: f64, st: f64, phi: f64) -> [f64; 3] {
    let helper = if d[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let e1 = normalize(cross(helper, d));
    let e2 = cross(d, e1);
    let (sp, cp) = phi.sin_cos();
    normalize([
        ct * d[0] + st * (cp * e1[0] + sp * e2[0]),
        ct * d[1] + st * (cp * e1[1] + sp * e2[1]),
        ct * d[2] + st * (cp * e1[2] + sp * e2[2]),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deflect_keeps_unit_length_and_angle() {
        for d in [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            normalize([1.0, 2.0, -3.0]),
        ] {
            for (theta, phi) in [(0.0, 0.0), (0.7, 1.9), (3.0, 5.0)] {
                let n = deflect(d, theta, phi);
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                assert!((len - 1.0).abs() < 1e-12);
                let dot = n[0] * d[0] + n[1] * d[1] + n[2] * d[2];
                assert!((dot - theta.cos()).abs() < 1e-12);
            }
        }
    }
}
