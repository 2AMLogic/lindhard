//! Dose to named layers, depth-dose through a stack, and interface
//! transmission, for the electron tally.
//!
//! Everything here is bookkeeping on the deposits and interface crossings the
//! transport already reports; nothing changes the transport. The accumulators
//! are opt-in through [`DoseConfig`] on
//! [`ElectronTallyConfig`](super::ElectronTallyConfig), so a run without it is
//! unchanged.
//!
//! # Conventions
//!
//! - **Per-history errors.** A history's contribution to a bin is the sum over
//!   the primary and all its secondaries. The reported mean is the sum over
//!   histories divided by the number of histories, and the standard error is
//!   `sqrt((Σx² − (Σx)²/N) / (N (N − 1)))`. It is `NaN` for fewer than two
//!   histories.
//! - **Dose per fluence.** The target is laterally infinite, so lateral escape
//!   cancels and the energy deposited in a layer per primary, divided by the
//!   layer's areal mass density `ρ t`, is the absorbed dose per unit incident
//!   fluence, in Gy m² (`J kg⁻¹ m²`). Multiply by a fluence in m⁻² for a dose.
//!   A semi-infinite last layer has no finite areal mass, so it is reported as
//!   energy only (its [`LayerDose::dose_gy_m2`] is `None`).
//! - **Depth-dose grid.** Each finite layer gets `bins_per_layer` uniform bins
//!   over its middle and, with [`InterfaceRefine`], `bins` finer uniform bins
//!   in the band of width `width_m` next to each of its faces (at most half
//!   the layer each side). A semi-infinite last layer is binned over
//!   [`DepthDoseConfig::substrate_depth_m`]; energy deeper than that is kept
//!   in [`DepthDose::outside_ev`]. A deposit is binned by the layer it was
//!   made in (the transport's `layer`), clamped to that layer's extent, so the
//!   bins of a layer sum to its layer total to summation error.
//! - **Deposition point.** The same as the rest of the electron tally: an
//!   inelastic loss at the collision point, the sub-threshold remainder at the
//!   stopping point.
//! - **Transmission.** One count per crossing from the transport's
//!   `interface` hook: an electron that crosses an interface several times is
//!   counted each time. Primaries and secondaries are separate entries. The
//!   energy is the kinetic energy on the new side of the face (after
//!   refraction with a step barrier). Forward is towards larger depth.
//! - **Deposits are measured from the band bottom**, like
//!   [`DepositionReport`](super::DepositionReport); with a secondary model in
//!   a layer with a Fermi energy they include the Fermi-sea energy of liberated
//!   conduction electrons.

use serde::{Deserialize, Serialize};

use super::hist::{Binning, Histogram};
use crate::geometry::Stack;
use crate::units::J_PER_EV;

/// Refined depth bins next to each interface.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct InterfaceRefine {
    /// Width of the refined band on each side of an interface, m. It is
    /// limited to half the layer, so the two bands of a thin layer do not
    /// overlap.
    pub width_m: f64,
    /// Uniform bins in each refined band.
    pub bins: usize,
}

/// The 1D depth-dose grid.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DepthDoseConfig {
    /// Uniform bins over the middle of each layer (the whole layer without
    /// refinement), at least 1.
    pub bins_per_layer: usize,
    /// Optional refinement next to each interface and the front face.
    #[serde(default)]
    pub refine: Option<InterfaceRefine>,
    /// Depth binned in a semi-infinite last layer, m; ignored without one.
    #[serde(default)]
    pub substrate_depth_m: f64,
}

/// Settings of the dose, depth-dose and transmission tallies.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DoseConfig {
    /// Depth-dose profile, if wanted. The per-layer totals are always on.
    #[serde(default)]
    pub depth: Option<DepthDoseConfig>,
    /// Count electrons crossing each interface, with spectra on the tally's
    /// escape energy binning.
    #[serde(default)]
    pub transmission: bool,
}

/// An invalid [`DoseConfig`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[error("invalid dose configuration: {0}")]
pub struct DoseError(pub String);

/// A mean with its standard error.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Estimate {
    /// Mean per primary.
    pub mean: f64,
    /// Standard error of the mean (`NaN` for fewer than two histories).
    pub std_err: f64,
}

impl Estimate {
    fn scaled(self, k: f64) -> Self {
        Self {
            mean: self.mean * k,
            std_err: self.std_err * k.abs(),
        }
    }
}

/// Per-history sums for a fixed set of bins. The current history's totals are
/// folded into `Σx` and `Σx²` by [`HistorySums::end_history`].
#[derive(Debug, Clone, PartialEq)]
struct HistorySums {
    histories: u64,
    sum: Vec<f64>,
    sum_sq: Vec<f64>,
    scratch: Vec<f64>,
    touched: Vec<usize>,
}

impl HistorySums {
    fn new(n: usize) -> Self {
        Self {
            histories: 0,
            sum: vec![0.0; n],
            sum_sq: vec![0.0; n],
            scratch: vec![0.0; n],
            touched: Vec::new(),
        }
    }

    fn add(&mut self, i: usize, x: f64) {
        if x <= 0.0 {
            return;
        }
        if self.scratch[i] == 0.0 {
            self.touched.push(i);
        }
        self.scratch[i] += x;
    }

    fn end_history(&mut self) {
        self.histories += 1;
        for &i in &self.touched {
            let s = self.scratch[i];
            self.sum[i] += s;
            self.sum_sq[i] += s * s;
            self.scratch[i] = 0.0;
        }
        self.touched.clear();
    }

    fn merge(&mut self, o: &HistorySums) {
        self.histories += o.histories;
        for (a, b) in self.sum.iter_mut().zip(&o.sum) {
            *a += b;
        }
        for (a, b) in self.sum_sq.iter_mut().zip(&o.sum_sq) {
            *a += b;
        }
    }

    fn estimate(&self, i: usize) -> Estimate {
        let n = self.histories as f64;
        if self.histories == 0 {
            return Estimate {
                mean: 0.0,
                std_err: f64::NAN,
            };
        }
        let mean = self.sum[i] / n;
        let std_err = if self.histories < 2 {
            f64::NAN
        } else {
            let var = (self.sum_sq[i] - self.sum[i] * self.sum[i] / n).max(0.0);
            (var / (n * (n - 1.0))).sqrt()
        };
        Estimate { mean, std_err }
    }
}

/// Static description of one layer, taken from the stack.
#[derive(Debug, Clone, PartialEq)]
struct LayerInfo {
    name: Option<String>,
    front_m: f64,
    back_m: f64,
    density_kg_m3: f64,
    /// Bin edges of this layer in the depth profile (`bins + 1` entries), and
    /// the index of its first bin; empty edges without a depth profile.
    edges: Vec<f64>,
    first_bin: usize,
}

impl LayerInfo {
    fn areal_mass(&self) -> Option<f64> {
        let t = self.back_m - self.front_m;
        t.is_finite().then_some(self.density_kg_m3 * t)
    }
}

fn layer_edges(front: f64, back: f64, finite: bool, cfg: &DepthDoseConfig) -> Vec<f64> {
    let t = back - front;
    let mut edges = Vec::new();
    let band = cfg
        .refine
        .map(|r| r.width_m.min(0.5 * t))
        .unwrap_or(0.0)
        .max(0.0);
    let fine = cfg.refine.map(|r| r.bins).unwrap_or(0);
    let refine_front = band > 0.0 && fine > 0;
    // A semi-infinite layer is refined at its front face only.
    let refine_back = refine_front && finite;
    let mid_lo = if refine_front { front + band } else { front };
    let mid_hi = if refine_back { back - band } else { back };
    if refine_front {
        for i in 0..fine {
            edges.push(front + band * i as f64 / fine as f64);
        }
    }
    for i in 0..cfg.bins_per_layer {
        edges.push(mid_lo + (mid_hi - mid_lo) * i as f64 / cfg.bins_per_layer as f64);
    }
    if refine_back {
        for i in 0..fine {
            edges.push(mid_hi + band * i as f64 / fine as f64);
        }
    }
    edges.push(back);
    // Drop coincident edges that rounding may produce.
    edges.dedup();
    edges
}

/// The accumulators of the dose tallies of one [`FullElectronTally`](super::FullElectronTally).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct DoseAcc {
    config: DoseConfig,
    layers: Vec<LayerInfo>,
    /// Per layer.
    layer_sums: HistorySums,
    /// Per depth bin, then one `outside` entry.
    depth_sums: Option<HistorySums>,
    n_bins: usize,
    /// Per `(interface, direction, class)`: count at `2 k`, energy at `2 k + 1`.
    transmission_sums: Option<HistorySums>,
    spectra: Vec<Histogram>,
}

impl DoseAcc {
    pub(crate) fn new(
        stack: &Stack,
        config: DoseConfig,
        escape_energy: Binning,
    ) -> Result<Self, DoseError> {
        if let Some(d) = &config.depth {
            if d.bins_per_layer == 0 {
                return Err(DoseError("bins_per_layer must be at least 1".into()));
            }
            if let Some(r) = &d.refine {
                if !(r.width_m.is_finite() && r.width_m >= 0.0) {
                    return Err(DoseError(format!(
                        "refine.width_m must be finite and non-negative, got {}",
                        r.width_m
                    )));
                }
            }
            let last_infinite = stack
                .layers()
                .last()
                .is_some_and(|l| !l.back_m().is_finite());
            if last_infinite && !(d.substrate_depth_m.is_finite() && d.substrate_depth_m > 0.0) {
                return Err(DoseError(
                    "substrate_depth_m must be finite and positive for a semi-infinite last layer"
                        .into(),
                ));
            }
        }
        let mut n_bins = 0usize;
        let mut layers = Vec::new();
        for l in stack.layers() {
            let (edges, first_bin) = match &config.depth {
                Some(d) => {
                    let back = if l.back_m().is_finite() {
                        l.back_m()
                    } else {
                        l.front_m() + d.substrate_depth_m
                    };
                    let e = layer_edges(l.front_m(), back, l.back_m().is_finite(), d);
                    let first = n_bins;
                    n_bins += e.len() - 1;
                    (e, first)
                }
                None => (Vec::new(), 0),
            };
            layers.push(LayerInfo {
                name: l.material().name().map(str::to_string),
                front_m: l.front_m(),
                back_m: l.back_m(),
                density_kg_m3: l.material().mass_density(),
                edges,
                first_bin,
            });
        }
        let n_if = layers.len().saturating_sub(1);
        Ok(Self {
            config,
            layer_sums: HistorySums::new(layers.len()),
            depth_sums: config.depth.map(|_| HistorySums::new(n_bins + 1)),
            n_bins,
            transmission_sums: config.transmission.then(|| HistorySums::new(n_if * 4 * 2)),
            spectra: if config.transmission {
                (0..n_if * 4)
                    .map(|_| Histogram::new(escape_energy))
                    .collect()
            } else {
                Vec::new()
            },
            layers,
        })
    }

    pub(crate) fn deposit(&mut self, layer: usize, x: f64, e: f64) {
        self.layer_sums.add(layer, e);
        if let Some(d) = &mut self.depth_sums {
            let info = &self.layers[layer];
            let last = *info.edges.last().expect("depth edges built");
            if x >= last && !info.back_m.is_finite() {
                d.add(self.n_bins, e);
                return;
            }
            let xc = x.clamp(info.edges[0], last);
            let nb = info.edges.len() - 1;
            let i = info.edges[1..nb].partition_point(|&edge| edge <= xc);
            d.add(info.first_bin + i, e);
        }
    }

    /// A crossing of interface `from -> to` (adjacent layers).
    pub(crate) fn crossing(&mut self, from: usize, to: usize, energy_ev: f64, secondary: bool) {
        let Some(t) = &mut self.transmission_sums else {
            return;
        };
        let (interface, backward) = if to > from { (from, 0) } else { (to, 1) };
        let k = (interface * 2 + backward) * 2 + usize::from(secondary);
        t.add(2 * k, 1.0);
        t.add(2 * k + 1, energy_ev);
        self.spectra[k].fill(energy_ev);
    }

    pub(crate) fn end_history(&mut self) {
        self.layer_sums.end_history();
        if let Some(d) = &mut self.depth_sums {
            d.end_history();
        }
        if let Some(t) = &mut self.transmission_sums {
            t.end_history();
        }
    }

    pub(crate) fn merge(&mut self, o: &DoseAcc) {
        self.layer_sums.merge(&o.layer_sums);
        if let (Some(a), Some(b)) = (&mut self.depth_sums, &o.depth_sums) {
            a.merge(b);
        }
        if let (Some(a), Some(b)) = (&mut self.transmission_sums, &o.transmission_sums) {
            a.merge(b);
        }
        for (a, b) in self.spectra.iter_mut().zip(&o.spectra) {
            for (x, y) in a.counts.iter_mut().zip(&b.counts) {
                *x += y;
            }
            a.underflow += b.underflow;
            a.overflow += b.overflow;
        }
    }

    pub(crate) fn report(&self, incident_ev: f64) -> DoseReport {
        let layers = self
            .layers
            .iter()
            .enumerate()
            .map(|(i, l)| {
                let e = self.layer_sums.estimate(i);
                LayerDose {
                    index: i,
                    name: l.name.clone(),
                    thickness_m: (l.back_m - l.front_m)
                        .is_finite()
                        .then_some(l.back_m - l.front_m),
                    areal_mass_kg_m2: l.areal_mass(),
                    energy_ev: e,
                    dose_gy_m2: l.areal_mass().map(|m| e.scaled(J_PER_EV / m)),
                }
            })
            .collect();
        let depth = self.depth_sums.as_ref().map(|d| {
            let mut edges_m = Vec::new();
            let mut layer_of_bin = Vec::new();
            let mut energy_ev = Vec::new();
            let mut ev_per_nm = Vec::new();
            let mut dose = Vec::new();
            for (li, l) in self.layers.iter().enumerate() {
                for b in 0..l.edges.len() - 1 {
                    let (lo, hi) = (l.edges[b], l.edges[b + 1]);
                    let e = d.estimate(l.first_bin + b);
                    let w = hi - lo;
                    edges_m.push(lo);
                    layer_of_bin.push(li);
                    energy_ev.push(e);
                    ev_per_nm.push(e.scaled(1e-9 / w));
                    dose.push(
                        l.back_m
                            .is_finite()
                            .then(|| e.scaled(J_PER_EV / (l.density_kg_m3 * w))),
                    );
                }
            }
            edges_m.push(
                self.layers
                    .last()
                    .and_then(|l| l.edges.last().copied())
                    .unwrap_or(0.0),
            );
            DepthDose {
                edges_m,
                layer_of_bin,
                energy_ev,
                ev_per_nm,
                dose_gy_m2: dose,
                outside_ev: d.estimate(self.n_bins),
            }
        });
        let transmission = self
            .transmission_sums
            .as_ref()
            .map(|t| {
                let n_if = self.layers.len().saturating_sub(1);
                let mut out = Vec::new();
                for interface in 0..n_if {
                    for backward in 0..2 {
                        for secondary in 0..2 {
                            let k = (interface * 2 + backward) * 2 + secondary;
                            let energy = t.estimate(2 * k + 1);
                            let total_energy = t.sum[2 * k + 1];
                            out.push(InterfaceTransmission {
                                interface,
                                from_layer: if backward == 0 {
                                    interface
                                } else {
                                    interface + 1
                                },
                                to_layer: if backward == 0 {
                                    interface + 1
                                } else {
                                    interface
                                },
                                secondary: secondary == 1,
                                count: t.estimate(2 * k),
                                energy_ev: energy,
                                energy_fraction: if incident_ev > 0.0 {
                                    total_energy / incident_ev
                                } else {
                                    0.0
                                },
                                spectrum: self.spectra[k].clone(),
                            });
                        }
                    }
                }
                out
            })
            .unwrap_or_default();
        DoseReport {
            layers,
            depth,
            transmission,
        }
    }
}

/// Energy and dose in one layer, per primary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LayerDose {
    /// Layer index (0 = front).
    pub index: usize,
    /// The material's name, if it has one.
    pub name: Option<String>,
    /// Thickness, m; `None` for a semi-infinite layer.
    pub thickness_m: Option<f64>,
    /// Areal mass density `ρ t`, kg/m²; `None` for a semi-infinite layer.
    pub areal_mass_kg_m2: Option<f64>,
    /// Energy deposited in the layer per primary, eV.
    pub energy_ev: Estimate,
    /// Absorbed dose per unit incident fluence, Gy m²; `None` for a
    /// semi-infinite layer.
    pub dose_gy_m2: Option<Estimate>,
}

/// The depth-dose profile through the stack.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DepthDose {
    /// Bin edges, m: bin `i` is `[edges_m[i], edges_m[i + 1])`, so there is one
    /// more edge than bins. Layer boundaries are shared edges.
    pub edges_m: Vec<f64>,
    /// The layer of each bin.
    pub layer_of_bin: Vec<usize>,
    /// Energy deposited per primary in each bin, eV.
    pub energy_ev: Vec<Estimate>,
    /// The same, per nm of depth, eV per primary per nm.
    pub ev_per_nm: Vec<Estimate>,
    /// Dose per unit incident fluence in each bin, Gy m²; `None` in a
    /// semi-infinite layer.
    pub dose_gy_m2: Vec<Option<Estimate>>,
    /// Energy deposited per primary deeper than the binned depth of a
    /// semi-infinite layer, eV.
    pub outside_ev: Estimate,
}

/// Electrons crossing one interface in one direction.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InterfaceTransmission {
    /// Interface index: between layers `interface` and `interface + 1`.
    pub interface: usize,
    /// Layer left.
    pub from_layer: usize,
    /// Layer entered. Forward is `to_layer > from_layer`.
    pub to_layer: usize,
    /// Secondaries (`true`) or primaries (`false`).
    pub secondary: bool,
    /// Crossings per primary.
    pub count: Estimate,
    /// Kinetic energy carried across per primary, eV.
    pub energy_ev: Estimate,
    /// That energy as a fraction of the incident energy of all primaries.
    pub energy_fraction: f64,
    /// Kinetic energy spectrum of the crossings, on the tally's escape energy
    /// binning (raw counts).
    pub spectrum: Histogram,
}

/// Result of the dose tallies.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DoseReport {
    /// Per layer.
    pub layers: Vec<LayerDose>,
    /// Depth-dose profile, if configured.
    pub depth: Option<DepthDose>,
    /// Interface crossings, if configured: four entries per interface
    /// (forward and back, primary then secondary).
    pub transmission: Vec<InterfaceTransmission>,
}
