//! The 1D layered [`Stack`]: ordered homogeneous layers, laterally infinite.
//! See the [module docs](super) for the coordinates and the `Geometry`
//! contract.

use super::{Exit, ExitOutcome, Face, Geometry, GeometryError};
use crate::material::Material;

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

impl Geometry for Stack {
    fn n_materials(&self) -> usize {
        self.layers.len()
    }

    fn material(&self, index: usize) -> &Material {
        &self.layers[index].material
    }

    fn n_regions(&self) -> usize {
        self.layers.len()
    }

    /// Each layer is its own material entry, so layer data is cached per
    /// layer.
    fn material_index(&self, region: usize) -> usize {
        region
    }

    fn locate(&self, pos: [f64; 3]) -> Option<usize> {
        self.layer_at(pos[0])
    }

    fn entry_point(&self) -> [f64; 3] {
        [0.0; 3]
    }

    fn exit(&self, region: usize, pos: [f64; 3], dir: [f64; 3], limit: f64) -> Option<Exit> {
        let layer = &self.layers[region];
        let (d, boundary) = if dir[0] > 0.0 {
            ((layer.back_m - pos[0]) / dir[0], layer.back_m)
        } else if dir[0] < 0.0 {
            ((layer.front_m - pos[0]) / dir[0], layer.front_m)
        } else {
            return None;
        };
        if d.is_nan() || d > limit {
            return None;
        }
        let at = [boundary, pos[1] + d * dir[1], pos[2] + d * dir[2]];
        let outward = dir[0] > 0.0;
        let outcome = if !outward && region == 0 {
            ExitOutcome::Escape {
                face: Face::Front,
                normal: [-1.0, 0.0, 0.0],
            }
        } else if outward && region + 1 == self.layers.len() {
            ExitOutcome::Escape {
                face: Face::Back,
                normal: [1.0, 0.0, 0.0],
            }
        } else {
            ExitOutcome::Enter {
                region: if outward { region + 1 } else { region - 1 },
                pos: at,
            }
        };
        Some(Exit {
            distance: d,
            at,
            outcome,
        })
    }

    fn in_vacuum(&self, pos: [f64; 3]) -> bool {
        pos[0] < 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_util::{c, enter, si, unit};
    use super::*;

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

    #[test]
    fn stack_exit_matches_the_layer_contract() {
        let s = Stack::new(vec![(si(), 1.0), (c(), 2.0)], None).unwrap();
        assert_eq!(s.n_materials(), 2);
        assert_eq!(s.locate([0.0; 3]), Some(0));
        assert_eq!(s.entry_point(), [0.0; 3]);
        assert!(s.in_vacuum([-1e-9, 0.0, 0.0]) && !s.in_vacuum([4.0, 0.0, 0.0]));
        let d = unit([1.0, 1.0, 0.0]);
        let e = s.exit(0, [0.25, 0.0, 0.0], d, 10.0).unwrap();
        assert_eq!(e.at[0], 1.0);
        assert_eq!(enter(&e).0, 1);
        let e = s.exit(1, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 10.0).unwrap();
        assert_eq!(e.distance, 2.0);
        assert!(matches!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Back,
                ..
            }
        ));
        let e = s.exit(0, [0.5, 0.0, 0.0], [-1.0, 0.0, 0.0], 10.0).unwrap();
        assert!(matches!(
            e.outcome,
            ExitOutcome::Escape {
                face: Face::Front,
                normal: [-1.0, 0.0, 0.0]
            }
        ));
        // Parallel to the interfaces, or too short: no event.
        assert!(s.exit(0, [0.5; 3], [0.0, 1.0, 0.0], 1e9).is_none());
        assert!(s.exit(0, [0.5; 3], [1.0, 0.0, 0.0], 0.1).is_none());
        // Substrate: no back face.
        let sub = Stack::semi_infinite(si());
        assert!(sub.exit(0, [0.5; 3], [1.0, 0.0, 0.0], 1e9).is_none());
    }
}
