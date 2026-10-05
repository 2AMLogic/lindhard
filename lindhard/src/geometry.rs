//! Target geometry. M0 provides the 1D layered target only.
//!
//! # Coordinates
//!
//! `x` is depth, measured from the front surface (`x = 0`) into the target;
//! `y` and `z` are lateral. Lengths are metres (SI, see [`crate::units`]).
//! A [`Stack`] is an ordered list of layers, each homogeneous and laterally
//! infinite (the planar slab model of amorphous-target binary-collision codes,
//! e.g. J. P. Biersack and L. G. Haggmark, Nucl. Instrum. Methods 174 (1980)
//! 257). The last layer may be semi-infinite (a substrate); otherwise the
//! target has a back face at its total thickness.

use crate::material::Material;

/// Errors from building a [`Stack`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum GeometryError {
    /// No layers and no substrate.
    #[error("a stack needs at least one layer or a substrate")]
    Empty,
    /// A layer thickness was zero, negative or not finite.
    #[error("layer {index} has invalid thickness {thickness_m} m: must be finite and positive")]
    InvalidThickness {
        /// Layer index (0 = front).
        index: usize,
        /// Offending value.
        thickness_m: f64,
    },
}

/// One homogeneous layer.
#[derive(Debug, Clone, PartialEq)]
pub struct Layer {
    material: Material,
    /// Depth of the front face, m.
    front_m: f64,
    /// Depth of the back face, m (`INFINITY` for a substrate).
    back_m: f64,
}

impl Layer {
    /// The layer material.
    pub fn material(&self) -> &Material {
        &self.material
    }

    /// Depth of the front face, m.
    pub fn front_m(&self) -> f64 {
        self.front_m
    }

    /// Depth of the back face, m (`f64::INFINITY` for a semi-infinite substrate).
    pub fn back_m(&self) -> f64 {
        self.back_m
    }

    /// Thickness, m (`f64::INFINITY` for a semi-infinite substrate).
    pub fn thickness_m(&self) -> f64 {
        self.back_m - self.front_m
    }
}

/// A 1D layered target: finite layers front to back, optionally ending in a
/// semi-infinite substrate.
#[derive(Debug, Clone, PartialEq)]
pub struct Stack {
    layers: Vec<Layer>,
}

impl Stack {
    /// Build from `(material, thickness in m)` pairs, front first, and an
    /// optional semi-infinite substrate behind them.
    pub fn new(
        layers: Vec<(Material, f64)>,
        substrate: Option<Material>,
    ) -> Result<Self, GeometryError> {
        if layers.is_empty() && substrate.is_none() {
            return Err(GeometryError::Empty);
        }
        let mut out = Vec::with_capacity(layers.len() + 1);
        let mut x = 0.0f64;
        for (index, (material, thickness_m)) in layers.into_iter().enumerate() {
            if !(thickness_m.is_finite() && thickness_m > 0.0) {
                return Err(GeometryError::InvalidThickness { index, thickness_m });
            }
            let back = x + thickness_m;
            out.push(Layer {
                material,
                front_m: x,
                back_m: back,
            });
            x = back;
        }
        if let Some(material) = substrate {
            out.push(Layer {
                material,
                front_m: x,
                back_m: f64::INFINITY,
            });
        }
        Ok(Self { layers: out })
    }

    /// A single semi-infinite material.
    pub fn semi_infinite(material: Material) -> Self {
        Self {
            layers: vec![Layer {
                material,
                front_m: 0.0,
                back_m: f64::INFINITY,
            }],
        }
    }

    /// Layers front to back.
    pub fn layers(&self) -> &[Layer] {
        &self.layers
    }

    /// Depth of the back face, m, or `None` if the last layer is semi-infinite.
    pub fn back_face_m(&self) -> Option<f64> {
        let b = self.layers.last().expect("non-empty").back_m;
        b.is_finite().then_some(b)
    }

    /// Index of the layer containing depth `x`, or `None` outside the target.
    /// A depth exactly on an interface belongs to the deeper layer.
    pub fn layer_at(&self, x: f64) -> Option<usize> {
        if x.is_nan() || x < 0.0 {
            return None;
        }
        self.layers.iter().position(|l| x < l.back_m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn si() -> Material {
        Material::from_atom_fractions(&[(14, 1.0)], None).unwrap()
    }

    #[test]
    fn boundaries_accumulate() {
        let s = Stack::new(vec![(si(), 1e-8), (si(), 2e-8)], Some(si())).unwrap();
        let b: Vec<_> = s
            .layers()
            .iter()
            .map(|l| (l.front_m(), l.back_m()))
            .collect();
        assert_eq!(b[0], (0.0, 1e-8));
        assert_eq!(b[1], (1e-8, 1e-8 + 2e-8));
        assert_eq!(b[2].1, f64::INFINITY);
        assert_eq!(s.back_face_m(), None);
        assert_eq!(s.layer_at(0.0), Some(0));
        assert_eq!(s.layer_at(1e-8), Some(1));
        assert_eq!(s.layer_at(1.0), Some(2));
        assert_eq!(s.layer_at(-1e-12), None);
        let f = Stack::new(vec![(si(), 1e-8)], None).unwrap();
        assert_eq!(f.back_face_m(), Some(1e-8));
        assert_eq!(f.layer_at(1e-8), None);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(Stack::new(vec![], None), Err(GeometryError::Empty));
        for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                Stack::new(vec![(si(), t)], None),
                Err(GeometryError::InvalidThickness { index: 0, .. })
            ));
        }
    }
}
