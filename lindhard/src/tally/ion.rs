//! The full ion tally: ranges, damage and escapes, wired to the BCA engine.
//!
//! [`IonTally`] implements [`BcaTally`]. Build one prototype with
//! [`IonTally::new`] and hand clones to [`crate::ion::bca::Bca::run`]; the
//! merged tally turns into a plain-data [`IonReport`] with
//! [`IonTally::report`].
//!
//! Units: positions and lengths in metres, energies in eV, angles in radians.
//! Totals are summed over all histories; divide by
//! [`IonReport::histories`] for per-ion values (the escape coefficients are
//! already per ion).
//!
//! # Determinism
//!
//! Counts are integers, moments use the mergeable accumulator of
//! [`super::moments`], and energy sums are plain `f64` sums; all are merged
//! field by field in the order [`crate::rng::run_particles`] calls `merge`
//! (chunk order). For a fixed `chunk_size` the result is bit-identical at any
//! thread count (`lindhard/tests/tally.rs`). Nothing is stored in a hash map.

use serde::{Deserialize, Serialize};

use super::hist::{Binning, Histogram};
use super::moments::{MomentSummary, Moments};
use super::pearson::{DualPearson, DualPearsonFit, PearsonError, PearsonIv};
use crate::geometry::Stack;
use crate::ion::bca::{BcaTally, EnergyBudget, Face, LatticeDeposit, Particle};
use crate::ion::damage::{damage_energy_ev, kinchin_pease, nrt_displacements};
use crate::material::MaterialError;

/// Binnings for an [`IonTally`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct IonTallyConfig {
    /// Depth `x`, m: range histograms (overall and per layer) and the
    /// defect depth histograms.
    pub depth: Binning,
    /// Lateral position `y` and `z`, m, of stopped primaries.
    pub lateral: Binning,
    /// Radial distance `(y^2 + z^2)^(1/2)`, m, of stopped primaries.
    pub radial: Binning,
    /// Kinetic energy outside the target of escaping particles, eV.
    pub escape_energy: Binning,
    /// Polar angle of escaping particles from the outward surface normal,
    /// rad (`0` is along the normal, `π/2` grazing).
    pub escape_polar: Binning,
}

/// Displacement-model estimates (NRT and Kinchin-Pease) from the PKA damage
/// energies. **Not** a count of simulated defects: see [`CascadeDefects`]
/// and [`crate::ion::damage`] for why the two are kept apart.
///
/// In a layer with more than one element, the Lindhard partition uses the
/// PKA's own `Z1`/`A1` and the layer's atom-fraction mean `Z2`/`A2`. NRT is
/// defined for monatomic targets, so this is an approximation and this
/// crate's own convention: see
/// [the compound-target section](crate::ion::damage#compound-and-multi-element-targets-an-approximation)
/// of [`crate::ion::damage`].
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct NrtDamage {
    /// Primary knock-on atoms: target atoms displaced by the beam particle.
    pub pka_count: u64,
    /// Sum of the PKA energies `T` (energy transferred, before the lattice
    /// binding is subtracted), eV.
    pub pka_energy_ev: f64,
    /// Sum of the PKA damage energies `T_dam` (Lindhard partition), eV.
    pub damage_energy_ev: f64,
    /// Sum of the NRT displacements `N_NRT(T_dam)` over PKAs.
    pub nrt_displacements: f64,
    /// Sum of the Kinchin-Pease displacements `N_KP(T_dam)` over PKAs.
    pub kinchin_pease_displacements: f64,
}

impl NrtDamage {
    fn merge(&mut self, o: &NrtDamage) {
        self.pka_count += o.pka_count;
        self.pka_energy_ev += o.pka_energy_ev;
        self.damage_energy_ev += o.damage_energy_ev;
        self.nrt_displacements += o.nrt_displacements;
        self.kinchin_pease_displacements += o.kinchin_pease_displacements;
    }
}

/// Defects counted event by event in the simulated cascades.
///
/// * `displacements`: target atoms given more than `E_d` (recoils created).
/// * `replacements`: a moving atom came to rest, immediately after the
///   collision in which it displaced an atom of its own element, at that
///   atom's site, so it fills the site. In this engine a particle stops when
///   its energy falls below its cutoff, so this count depends on the cutoffs
///   (a recoil cutoff near `E_d` gives the most replacements). When a recoil
///   is not followed (below its cutoff, or with
///   [`BcaConfig::follow_recoils`](crate::ion::bca::BcaConfig::follow_recoils)
///   `= false`), the engine reports `stopped(&r)` straight after
///   `recoil(&r)`, which clears the pending site, so the displacer cannot
///   score a replacement in that collision. With `follow_recoils = false` the
///   replacement count is therefore always 0.
/// * `vacancies = displacements - replacements`.
/// * `interstitials`: recoils that came to rest in the target and did not
///   fill a site. Implanted beam particles are not counted here (see the
///   range tally). Recoils that leave the target create a vacancy and no
///   interstitial, so the two counts differ.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct CascadeDefects {
    /// Atoms displaced.
    pub displacements: u64,
    /// Replacement events.
    pub replacements: u64,
    /// Vacancies left.
    pub vacancies: u64,
    /// Interstitials left.
    pub interstitials: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct DefectCounts {
    displacements: u64,
    replacements: u64,
    interstitials: u64,
}

impl DefectCounts {
    fn merge(&mut self, o: &DefectCounts) {
        self.displacements += o.displacements;
        self.replacements += o.replacements;
        self.interstitials += o.interstitials;
    }

    fn report(&self) -> CascadeDefects {
        CascadeDefects {
            displacements: self.displacements,
            replacements: self.replacements,
            vacancies: self.displacements - self.replacements,
            interstitials: self.interstitials,
        }
    }
}

/// Per-layer constants for the damage models (not merged).
#[derive(Debug, Clone, PartialEq)]
struct LayerConsts {
    /// `(Z, E_d, E_b)` of each element, eV.
    elems: Vec<(u8, f64, f64)>,
    /// Atom-fraction mean atomic number.
    mean_z: f64,
    /// Mean atomic weight, u.
    mean_a: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct LayerAcc {
    stopped: Moments,
    depth_hist: Histogram,
    nrt: NrtDamage,
    defects: DefectCounts,
}

#[derive(Debug, Clone, PartialEq)]
struct FaceAcc {
    count: u64,
    energy_ev: f64,
    energy_hist: Histogram,
    polar_hist: Histogram,
}

impl FaceAcc {
    fn new(c: &IonTallyConfig) -> Self {
        Self {
            count: 0,
            energy_ev: 0.0,
            energy_hist: Histogram::new(c.escape_energy),
            polar_hist: Histogram::new(c.escape_polar),
        }
    }

    fn merge(&mut self, o: &FaceAcc) {
        self.count += o.count;
        self.energy_ev += o.energy_ev;
        self.energy_hist.merge(&o.energy_hist);
        self.polar_hist.merge(&o.polar_hist);
    }

    fn report(&self, histories: u64) -> FaceEscape {
        FaceEscape {
            count: self.count,
            per_ion: per(self.count as f64, histories),
            energy_ev: self.energy_ev,
            mean_energy_ev: if self.count > 0 {
                self.energy_ev / self.count as f64
            } else {
                0.0
            },
            energy_histogram: self.energy_hist.clone(),
            polar_histogram: self.polar_hist.clone(),
        }
    }
}

fn per(x: f64, histories: u64) -> f64 {
    if histories == 0 {
        0.0
    } else {
        x / histories as f64
    }
}

/// The site of the most recent displacement, for replacement detection.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Site {
    pos: [f64; 3],
    z: u8,
    displacer_generation: u32,
}

/// Range, damage and escape tally for the BCA engine. See the module docs.
#[derive(Debug, Clone, PartialEq)]
pub struct IonTally {
    config: IonTallyConfig,
    species_z: Vec<u8>,
    consts: Vec<LayerConsts>,
    histories: u64,
    budget: EnergyBudget,
    depth: Moments,
    lateral_y: Moments,
    lateral_z: Moments,
    radial: Moments,
    depth_hist: Histogram,
    lateral_y_hist: Histogram,
    lateral_z_hist: Histogram,
    radial_hist: Histogram,
    layers: Vec<LayerAcc>,
    displacement_hist: Histogram,
    replacement_hist: Histogram,
    interstitial_hist: Histogram,
    /// `[species] -> (front, back)`.
    escapes: Vec<(FaceAcc, FaceAcc)>,
    /// Per-history scratch: not part of the result, reset every history.
    last_site: Option<Site>,
}

impl IonTally {
    /// Empty tally for `stack`, with the particle species of the engine
    /// (`species_z`, as returned by [`crate::ion::bca::Bca::species_z`]:
    /// index 0 is the beam).
    ///
    /// Every element of every layer must have `E_d` and `E_b` set (the engine
    /// requires the same).
    pub fn new(
        stack: &Stack,
        species_z: &[u8],
        config: IonTallyConfig,
    ) -> Result<Self, MaterialError> {
        let mut consts = Vec::with_capacity(stack.layers().len());
        for l in stack.layers() {
            let m = l.material();
            let mut elems = Vec::new();
            let mut mean_z = 0.0;
            for (z, x) in m.atom_fractions() {
                elems.push((
                    z,
                    m.displacement_energy_ev(z)?,
                    m.lattice_binding_energy_ev(z)?,
                ));
                mean_z += x * f64::from(z);
            }
            consts.push(LayerConsts {
                elems,
                mean_z,
                mean_a: m.mean_atomic_weight(),
            });
        }
        let layers = (0..consts.len())
            .map(|_| LayerAcc {
                stopped: Moments::new(),
                depth_hist: Histogram::new(config.depth),
                nrt: NrtDamage::default(),
                defects: DefectCounts::default(),
            })
            .collect();
        Ok(Self {
            config,
            species_z: species_z.to_vec(),
            consts,
            histories: 0,
            budget: EnergyBudget::default(),
            depth: Moments::new(),
            lateral_y: Moments::new(),
            lateral_z: Moments::new(),
            radial: Moments::new(),
            depth_hist: Histogram::new(config.depth),
            lateral_y_hist: Histogram::new(config.lateral),
            lateral_z_hist: Histogram::new(config.lateral),
            radial_hist: Histogram::new(config.radial),
            layers,
            displacement_hist: Histogram::new(config.depth),
            replacement_hist: Histogram::new(config.depth),
            interstitial_hist: Histogram::new(config.depth),
            escapes: species_z
                .iter()
                .map(|_| (FaceAcc::new(&config), FaceAcc::new(&config)))
                .collect(),
            last_site: None,
        })
    }

    /// Histories tallied so far.
    pub fn histories(&self) -> u64 {
        self.histories
    }

    /// `E_d` and `E_b` of element `z` in layer `layer`.
    fn energies(&self, layer: usize, z: u8) -> (f64, f64) {
        self.consts[layer]
            .elems
            .iter()
            .find(|e| e.0 == z)
            .map(|e| (e.1, e.2))
            .expect("recoils are elements of their layer")
    }

    /// The plain-data result. With `fit_dual_pearson` the depth histogram is
    /// also fitted with [`DualPearson::fit`] (slower: a few thousand
    /// evaluations over the histogram).
    pub fn report(&self, fit_dual_pearson: bool) -> IonReport {
        let h = self.histories;
        let depth = self.depth.summary();
        let (pearson_iv, pearson_iv_error) = match depth {
            None => (None, None),
            Some(s) => match PearsonIv::from_moments(s.mean, s.std_dev, s.skewness, s.kurtosis) {
                Ok(p) => (Some(p), None),
                Err(e) => (None, Some(e)),
            },
        };
        let (dual_pearson, dual_pearson_error) = if fit_dual_pearson {
            match DualPearson::fit(&self.depth_hist) {
                Ok(f) => (Some(f), None),
                Err(e) => (None, Some(e)),
            }
        } else {
            (None, None)
        };
        let range = RangeReport {
            stopped: self.depth.n,
            depth,
            pearson_iv,
            pearson_iv_error,
            dual_pearson,
            dual_pearson_error,
            depth_histogram: self.depth_hist.clone(),
            lateral_y: self.lateral_y.summary(),
            lateral_z: self.lateral_z.summary(),
            radial: self.radial.summary(),
            lateral_y_histogram: self.lateral_y_hist.clone(),
            lateral_z_histogram: self.lateral_z_hist.clone(),
            radial_histogram: self.radial_hist.clone(),
            layers: self
                .layers
                .iter()
                .enumerate()
                .map(|(i, l)| LayerRange {
                    layer: i,
                    stopped: l.stopped.n,
                    depth: l.stopped.summary(),
                    depth_histogram: l.depth_hist.clone(),
                })
                .collect(),
        };
        let mut nrt = NrtDamage::default();
        let mut defects = DefectCounts::default();
        for l in &self.layers {
            nrt.merge(&l.nrt);
            defects.merge(&l.defects);
        }
        let mut vacancy_hist = self.displacement_hist.clone();
        for (v, r) in vacancy_hist
            .counts
            .iter_mut()
            .zip(&self.replacement_hist.counts)
        {
            *v -= r;
        }
        vacancy_hist.underflow -= self.replacement_hist.underflow;
        vacancy_hist.overflow -= self.replacement_hist.overflow;
        let damage = DamageReport {
            nrt,
            cascade: defects.report(),
            layers: self
                .layers
                .iter()
                .enumerate()
                .map(|(i, l)| LayerDamage {
                    layer: i,
                    nrt: l.nrt,
                    cascade: l.defects.report(),
                })
                .collect(),
            vacancy_histogram: vacancy_hist,
            interstitial_histogram: self.interstitial_hist.clone(),
            replacement_histogram: self.replacement_hist.clone(),
        };
        let species: Vec<SpeciesEscape> = self
            .escapes
            .iter()
            .enumerate()
            .map(|(i, (f, b))| SpeciesEscape {
                species: i,
                z: self.species_z[i],
                beam: i == 0,
                front: f.report(h),
                back: b.report(h),
            })
            .collect();
        let sputtered: u64 = self.escapes.iter().skip(1).map(|e| e.0.count).sum();
        let escapes = EscapeReport {
            backscatter_coefficient: per(self.escapes[0].0.count as f64, h),
            transmission_coefficient: per(self.escapes[0].1.count as f64, h),
            energy_reflection_coefficient: if self.budget.incident > 0.0 {
                self.budget.backscattered / self.budget.incident
            } else {
                0.0
            },
            sputter_yield: per(sputtered as f64, h),
            species,
        };
        IonReport {
            histories: h,
            species_z: self.species_z.clone(),
            budget: self.budget,
            range,
            damage,
            escapes,
        }
    }
}

impl BcaTally for IonTally {
    fn begin_history(&mut self, _index: u64) {
        self.last_site = None;
    }

    fn lattice(&mut self, _at: [f64; 3], _layer: usize, kind: LatticeDeposit, _energy_ev: f64) {
        if kind == LatticeDeposit::Subthreshold {
            self.last_site = None;
        }
    }

    fn recoil(&mut self, r: &Particle) {
        let layer = r.layer;
        let (e_d, e_b) = self.energies(layer, r.z);
        let acc = &mut self.layers[layer];
        acc.defects.displacements += 1;
        self.displacement_hist.fill(r.pos[0]);
        if r.generation == 1 {
            // A PKA: the transfer T is the recoil energy plus the binding
            // energy the engine left in the lattice.
            let c = &self.consts[layer];
            let t = r.energy_ev + e_b;
            let t_dam = damage_energy_ev(t, f64::from(r.z), r.mass_amu, c.mean_z, c.mean_a);
            acc.nrt.pka_count += 1;
            acc.nrt.pka_energy_ev += t;
            acc.nrt.damage_energy_ev += t_dam;
            acc.nrt.nrt_displacements += nrt_displacements(t_dam, e_d);
            acc.nrt.kinchin_pease_displacements += kinchin_pease(t_dam, e_d);
        }
        self.last_site = Some(Site {
            pos: r.pos,
            z: r.z,
            displacer_generation: r.generation.saturating_sub(1),
        });
    }

    fn stopped(&mut self, p: &Particle) {
        let replaced = matches!(self.last_site, Some(s)
            if s.pos == p.pos && s.z == p.z && s.displacer_generation == p.generation);
        self.last_site = None;
        let x = p.pos[0];
        if replaced {
            self.layers[p.layer].defects.replacements += 1;
            self.replacement_hist.fill(x);
        } else if !p.is_primary() {
            self.layers[p.layer].defects.interstitials += 1;
            self.interstitial_hist.fill(x);
        }
        if p.is_primary() {
            let (y, z) = (p.pos[1], p.pos[2]);
            let r = y.hypot(z);
            self.depth.push(x);
            self.lateral_y.push(y);
            self.lateral_z.push(z);
            self.radial.push(r);
            self.depth_hist.fill(x);
            self.lateral_y_hist.fill(y);
            self.lateral_z_hist.fill(z);
            self.radial_hist.fill(r);
            let l = &mut self.layers[p.layer];
            l.stopped.push(x);
            l.depth_hist.fill(x);
        }
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        self.last_site = None;
        let polar = p.dir[0].abs().min(1.0).acos();
        let pair = &mut self.escapes[p.species];
        let acc = match face {
            Face::Front => &mut pair.0,
            Face::Back => &mut pair.1,
            // Lateral faces exist only on voxel grids, which this stack
            // tally does not describe.
            Face::Side => return,
        };
        acc.count += 1;
        acc.energy_ev += p.energy_ev;
        acc.energy_hist.fill(p.energy_ev);
        acc.polar_hist.fill(polar);
    }

    fn end_history(&mut self, _index: u64, budget: &EnergyBudget) {
        self.histories += 1;
        self.budget.add(budget);
        self.last_site = None;
    }

    /// Field-wise merge.
    ///
    /// # Panics
    /// If `other` was built with a different stack, species list or
    /// configuration.
    fn merge(&mut self, o: Self) {
        assert!(
            self.config == o.config && self.species_z == o.species_z && self.consts == o.consts,
            "merging ion tallies built for different runs"
        );
        self.histories += o.histories;
        self.budget.add(&o.budget);
        self.depth.merge(&o.depth);
        self.lateral_y.merge(&o.lateral_y);
        self.lateral_z.merge(&o.lateral_z);
        self.radial.merge(&o.radial);
        self.depth_hist.merge(&o.depth_hist);
        self.lateral_y_hist.merge(&o.lateral_y_hist);
        self.lateral_z_hist.merge(&o.lateral_z_hist);
        self.radial_hist.merge(&o.radial_hist);
        for (a, b) in self.layers.iter_mut().zip(&o.layers) {
            a.stopped.merge(&b.stopped);
            a.depth_hist.merge(&b.depth_hist);
            a.nrt.merge(&b.nrt);
            a.defects.merge(&b.defects);
        }
        self.displacement_hist.merge(&o.displacement_hist);
        self.replacement_hist.merge(&o.replacement_hist);
        self.interstitial_hist.merge(&o.interstitial_hist);
        for (a, b) in self.escapes.iter_mut().zip(&o.escapes) {
            a.0.merge(&b.0);
            a.1.merge(&b.1);
        }
    }
}

/// Everything an [`IonTally`] measured, as plain data.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IonReport {
    /// Primary histories.
    pub histories: u64,
    /// Atomic number by species index (0 is the beam).
    pub species_z: Vec<u8>,
    /// Summed energy budget of all histories, eV.
    pub budget: EnergyBudget,
    /// Range distributions of the stopped beam particles.
    pub range: RangeReport,
    /// Damage: model estimates and cascade counts, side by side.
    pub damage: DamageReport,
    /// Backscattering, transmission and sputtering.
    pub escapes: EscapeReport,
}

/// Range distributions of beam particles that came to rest in the target.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RangeReport {
    /// Beam particles at rest in the target.
    pub stopped: u64,
    /// Depth moments (`Rp`, `ΔRp`, `γ`, `β`), m. `None` with fewer than two
    /// stopped particles.
    pub depth: Option<MomentSummary>,
    /// Pearson IV with the depth moments, if they are in the type IV region.
    pub pearson_iv: Option<PearsonIv>,
    /// Why `pearson_iv` is `None`, when there were moments to try.
    pub pearson_iv_error: Option<PearsonError>,
    /// Dual-Pearson fit to the depth histogram, if requested and possible.
    pub dual_pearson: Option<DualPearsonFit>,
    /// Why the requested dual-Pearson fit failed.
    pub dual_pearson_error: Option<PearsonError>,
    /// Depth histogram, m.
    pub depth_histogram: Histogram,
    /// Moments of the lateral coordinate `y`, m.
    pub lateral_y: Option<MomentSummary>,
    /// Moments of the lateral coordinate `z`, m.
    pub lateral_z: Option<MomentSummary>,
    /// Moments of the radial distance from the beam axis through the entry
    /// point, m.
    pub radial: Option<MomentSummary>,
    /// Histogram of `y`, m.
    pub lateral_y_histogram: Histogram,
    /// Histogram of `z`, m.
    pub lateral_z_histogram: Histogram,
    /// Histogram of the radial distance, m.
    pub radial_histogram: Histogram,
    /// The same depth statistics per layer (by where the particle stopped).
    pub layers: Vec<LayerRange>,
}

/// Depth statistics of the beam particles that stopped in one layer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LayerRange {
    /// Layer index (0 is the front layer).
    pub layer: usize,
    /// Beam particles at rest in this layer.
    pub stopped: u64,
    /// Their depth moments, m.
    pub depth: Option<MomentSummary>,
    /// Their depth histogram, m (overall binning).
    pub depth_histogram: Histogram,
}

/// Damage, as two separate kinds of result. See [`crate::ion::damage`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DamageReport {
    /// NRT / Kinchin-Pease model estimates, all layers.
    pub nrt: NrtDamage,
    /// Simulated cascade defect counts, all layers.
    pub cascade: CascadeDefects,
    /// Both, per layer (by the layer of the event).
    pub layers: Vec<LayerDamage>,
    /// Depth histogram of vacancies, m.
    pub vacancy_histogram: Histogram,
    /// Depth histogram of interstitials, m.
    pub interstitial_histogram: Histogram,
    /// Depth histogram of replacements, m.
    pub replacement_histogram: Histogram,
}

/// Damage in one layer.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LayerDamage {
    /// Layer index.
    pub layer: usize,
    /// Model estimates for PKAs created in this layer.
    pub nrt: NrtDamage,
    /// Cascade counts for events in this layer.
    pub cascade: CascadeDefects,
}

/// Particles leaving the target.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EscapeReport {
    /// Beam particles leaving the front face, per incident ion.
    pub backscatter_coefficient: f64,
    /// Beam particles leaving the back face, per incident ion.
    pub transmission_coefficient: f64,
    /// Energy carried out of the front face by beam particles, as a fraction
    /// of the incident energy.
    pub energy_reflection_coefficient: f64,
    /// Target atoms leaving the front face, all species, per incident ion.
    pub sputter_yield: f64,
    /// Per species (index 0 is the beam).
    pub species: Vec<SpeciesEscape>,
}

/// Escapes of one species through each face.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpeciesEscape {
    /// Species index.
    pub species: usize,
    /// Atomic number.
    pub z: u8,
    /// True for the beam species (front escapes are backscattering), false
    /// for a target element (front escapes are sputtering).
    pub beam: bool,
    /// Through the front face (`x = 0`).
    pub front: FaceEscape,
    /// Through the back face of a finite stack.
    pub back: FaceEscape,
}

/// Count and spectra of particles leaving through one face.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FaceEscape {
    /// Particles.
    pub count: u64,
    /// Particles per incident ion (backscatter or transmission coefficient,
    /// or partial sputter yield).
    pub per_ion: f64,
    /// Sum of their energies outside the target, eV.
    pub energy_ev: f64,
    /// Mean energy outside the target, eV (0 if none escaped).
    pub mean_energy_ev: f64,
    /// Energy spectrum, eV.
    pub energy_histogram: Histogram,
    /// Polar angle from the outward normal, rad.
    pub polar_histogram: Histogram,
}
