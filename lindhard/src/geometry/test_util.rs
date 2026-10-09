//! Helpers shared by the geometry unit tests.

use super::{Exit, ExitOutcome};
use crate::material::Material;

pub(super) fn si() -> Material {
    Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
}

pub(super) fn c() -> Material {
    Material::from_atom_fractions(&[(6, 1.0)], None).unwrap()
}

pub(super) fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

pub(super) fn enter(e: &Exit) -> (usize, [f64; 3]) {
    match e.outcome {
        ExitOutcome::Enter { region, pos } => (region, pos),
        ExitOutcome::Escape { .. } => panic!("expected Enter, got {e:?}"),
    }
}
