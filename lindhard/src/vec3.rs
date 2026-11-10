//! Shared helpers for 3-vectors stored as `[f64; 3]`.
//!
//! These are plain component-wise expressions; the evaluation order is fixed
//! here so every caller gets bit-identical results (the determinism-across-
//! thread-counts invariant depends on it).

/// A 3-vector in Cartesian components.
pub(crate) type V3 = [f64; 3];

/// `a - b`.
pub(crate) fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Right-handed cross product `a x b`.
pub(crate) fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product `a . b`.
pub(crate) fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean length `|a|`.
pub(crate) fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// `a * s` (component-wise scaling).
pub(crate) fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Normalise a vector to unit length.
pub(crate) fn normalize(v: V3) -> V3 {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}
