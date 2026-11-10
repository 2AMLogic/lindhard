//! The photon transport loop.

use crate::geometry::{ExitOutcome, Flight, Geometry};
use crate::material::MaterialError;
use crate::rng::{run_particles, uniform, ParticleRng};
use crate::tally::PhotonTally;
use crate::vec3::{norm, normalize, V3};

use super::data::{AtomPhotonData, PhotonTable, ENERGY_MAX_EV, ENERGY_MIN_EV, Z_MAX};
use super::sample::{rotate, sample_incoherent, CoherentSampler};

/// 1 barn in m².
const BARN_M2: f64 = 1.0e-28;

/// Errors from setting up or running photon transport.
#[derive(Debug, thiserror::Error)]
pub enum PhotonTransportError {
    /// A material contains an element the photon table lacks.
    #[error("the photon table has no data for element Z = {z}")]
    MissingElement {
        /// Atomic number.
        z: u8,
    },
    /// A material is unusable.
    #[error("material error: {0}")]
    Material(#[from] MaterialError),
    /// The source energy is outside the table.
    #[error(
        "source energy {energy_ev} eV is outside the photon table range \
         {ENERGY_MIN_EV}..={ENERGY_MAX_EV} eV"
    )]
    EnergyOutOfRange {
        /// The offending energy, eV.
        energy_ev: f64,
    },
    /// The source direction is not a finite non-zero vector.
    #[error("source direction must be a finite non-zero vector")]
    BadDirection,
    /// The source position is not finite.
    #[error("source position must be finite")]
    BadPosition,
}

/// A monoenergetic, mono-directional point source of photons.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhotonSource {
    /// Photon energy, eV (within the table range, 1 keV to 1 MeV).
    pub energy_ev: f64,
    /// Where the photons start, m; `None` for the geometry's entry point.
    pub position: Option<[f64; 3]>,
    /// Direction of flight (normalised internally); `+x` for [`Self::beam`].
    pub direction: [f64; 3],
}

impl PhotonSource {
    /// A beam of photons of `energy_ev` entering at the geometry's entry
    /// point along `+x`.
    pub fn beam(energy_ev: f64) -> Self {
        Self {
            energy_ev,
            position: None,
            direction: [1.0, 0.0, 0.0],
        }
    }
}

/// The elements of one material: atoms per m³ and their photon data.
#[derive(Debug)]
struct MaterialPhotonData<'a> {
    elements: Vec<(f64, &'a AtomPhotonData)>,
}

impl MaterialPhotonData<'_> {
    /// Macroscopic total cross section, 1/m.
    fn sigma_total(&self, e_ev: f64) -> f64 {
        self.elements
            .iter()
            .map(|(n, a)| n * a.total_barns(e_ev))
            .sum::<f64>()
            * BARN_M2
    }
}

/// Monte Carlo transport of photons through a [`Geometry`].
///
/// See the [module docs](super) for the physics. Built once per run: it holds
/// the per-material element lists and the per-element coherent-angle samplers,
/// all read-only afterwards, so histories run in parallel without locks.
pub struct PhotonTransport<'a> {
    geometry: &'a dyn Geometry,
    materials: Vec<MaterialPhotonData<'a>>,
    coherent: Vec<Option<CoherentSampler>>,
}

impl std::fmt::Debug for PhotonTransport<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PhotonTransport")
            .field("n_materials", &self.materials.len())
            .field("n_regions", &self.geometry.n_regions())
            .finish_non_exhaustive()
    }
}

impl<'a> PhotonTransport<'a> {
    /// Prepare transport through `geometry` with the data in `table`; every
    /// element of every material must be in the table.
    pub fn new(
        geometry: &'a dyn Geometry,
        table: &'a PhotonTable,
    ) -> Result<Self, PhotonTransportError> {
        let mut coherent: Vec<Option<CoherentSampler>> = vec![None; usize::from(Z_MAX) + 1];
        let mut materials = Vec::with_capacity(geometry.n_materials());
        for m in 0..geometry.n_materials() {
            let mat = geometry.material(m);
            let mut elements = Vec::new();
            for c in mat.components() {
                let z = c.z();
                let atom = table
                    .atom(z)
                    .ok_or(PhotonTransportError::MissingElement { z })?;
                coherent[usize::from(z)]
                    .get_or_insert_with(|| CoherentSampler::new(&atom.form_factor));
                elements.push((mat.number_density_of(z)?, atom));
            }
            materials.push(MaterialPhotonData { elements });
        }
        Ok(Self {
            geometry,
            materials,
            coherent,
        })
    }

    /// The geometry.
    pub fn geometry(&self) -> &dyn Geometry {
        self.geometry
    }

    /// Macroscopic cross sections of material `index` at `energy_ev`:
    /// `(photoelectric, incoherent, coherent)`, 1/m.
    pub fn macroscopic_cross_sections(&self, index: usize, energy_ev: f64) -> [f64; 3] {
        let mut s = [0.0; 3];
        for (n, a) in &self.materials[index].elements {
            s[0] += n * a.photoelectric.eval(energy_ev);
            s[1] += n * a.incoherent.eval(energy_ev);
            s[2] += n * a.coherent.eval(energy_ev);
        }
        s.map(|v| v * BARN_M2)
    }

    /// Total linear attenuation coefficient of material `index`, 1/m.
    pub fn attenuation_coefficient(&self, index: usize, energy_ev: f64) -> f64 {
        self.materials[index].sigma_total(energy_ev)
    }

    fn check_source(&self, src: &PhotonSource) -> Result<(), PhotonTransportError> {
        if !(ENERGY_MIN_EV..=ENERGY_MAX_EV).contains(&src.energy_ev) {
            return Err(PhotonTransportError::EnergyOutOfRange {
                energy_ev: src.energy_ev,
            });
        }
        let d = src.direction;
        if !d.iter().all(|v| v.is_finite()) || norm(d) == 0.0 {
            return Err(PhotonTransportError::BadDirection);
        }
        if src
            .position
            .is_some_and(|p| !p.iter().all(|v| v.is_finite()))
        {
            return Err(PhotonTransportError::BadPosition);
        }
        Ok(())
    }

    /// Run `n_photons` histories on the current rayon pool and return the
    /// merged tally. Bit-identical for any thread count (see
    /// [`crate::rng::run_particles`]); `chunk_size` must not depend on it.
    pub fn run(
        &self,
        source: &PhotonSource,
        seed: u64,
        n_photons: u64,
        chunk_size: u64,
    ) -> Result<PhotonTally, PhotonTransportError> {
        self.check_source(source)?;
        let n_regions = self.geometry.n_regions();
        Ok(run_particles(
            seed,
            n_photons,
            chunk_size,
            || PhotonTally::new(n_regions),
            |tally, rng, _index| self.history(source, rng, tally),
            |total, chunk| total.merge(chunk),
        ))
    }

    /// Follow one photon to its end and add it to `tally`. The source must
    /// have passed the checks of [`Self::run`]; use `run` unless a caller
    /// needs to drive the stream itself.
    ///
    /// Flights are exact: the path to the next interaction is drawn from the
    /// current material's total cross section, the geometry says whether a
    /// face comes first, and at a face the photon is moved there and the draw
    /// is repeated in the next region (the exponential law is memoryless).
    pub fn history(&self, source: &PhotonSource, rng: &mut ParticleRng, tally: &mut PhotonTally) {
        tally.n_histories += 1;
        tally.incident_ev += source.energy_ev;
        let mut pos: V3 = source
            .position
            .unwrap_or_else(|| self.geometry.entry_point());
        let mut dir = normalize(source.direction);
        let Some(mut region) = self.geometry.locate(pos) else {
            tally.n_missed += 1;
            tally.missed_ev += source.energy_ev;
            return;
        };
        let mut e = source.energy_ev;
        let mut uncollided = true;
        loop {
            if e < ENERGY_MIN_EV {
                tally.cutoff_ev[region] += e;
                tally.n_cutoff += 1;
                return;
            }
            let mat = &self.materials[self.geometry.material_index(region)];
            let s = -(1.0 - uniform(rng)).ln() / mat.sigma_total(e);
            match self.geometry.flight(region, pos, dir, s) {
                Flight::Event(exit) => match exit.outcome {
                    ExitOutcome::Enter { region: r, pos: p } => {
                        region = r;
                        pos = p;
                    }
                    ExitOutcome::Escape { face, .. } => {
                        tally.record_escape(face, e, uncollided);
                        return;
                    }
                },
                Flight::Clear { region: r } => {
                    region = r;
                    pos = [
                        pos[0] + s * dir[0],
                        pos[1] + s * dir[1],
                        pos[2] + s * dir[2],
                    ];
                    match self.interact(mat, e, rng) {
                        Interaction::Photoelectric => {
                            tally.photoelectric_ev[region] += e;
                            tally.n_photoelectric += 1;
                            return;
                        }
                        Interaction::Incoherent { cos_t, e_out } => {
                            tally.compton_ev[region] += e - e_out;
                            tally.n_incoherent += 1;
                            dir = rotate(dir, cos_t, std::f64::consts::TAU * uniform(rng));
                            e = e_out;
                        }
                        Interaction::Coherent { cos_t } => {
                            tally.n_coherent += 1;
                            dir = rotate(dir, cos_t, std::f64::consts::TAU * uniform(rng));
                        }
                    }
                    uncollided = false;
                }
            }
        }
    }

    /// Choose the element, then the process, of an interaction at `e`.
    fn interact(&self, mat: &MaterialPhotonData<'_>, e: f64, rng: &mut ParticleRng) -> Interaction {
        let total = mat.sigma_total(e) / BARN_M2;
        let target = uniform(rng) * total;
        let mut acc = 0.0;
        let mut chosen = mat.elements.len() - 1;
        for (i, (n, a)) in mat.elements.iter().enumerate() {
            acc += n * a.total_barns(e);
            if target < acc {
                chosen = i;
                break;
            }
        }
        let atom = mat.elements[chosen].1;
        let (pe, inc) = (atom.photoelectric.eval(e), atom.incoherent.eval(e));
        let u = uniform(rng) * atom.total_barns(e);
        if u < pe {
            Interaction::Photoelectric
        } else if u < pe + inc {
            let (cos_t, e_out) = sample_incoherent(rng, e, atom.z, &atom.scattering_function);
            Interaction::Incoherent { cos_t, e_out }
        } else {
            let sampler = self.coherent[usize::from(atom.z)]
                .as_ref()
                .expect("built for every element of every material");
            Interaction::Coherent {
                cos_t: sampler.sample_cos(rng, e),
            }
        }
    }
}

enum Interaction {
    Photoelectric,
    Incoherent { cos_t: f64, e_out: f64 },
    Coherent { cos_t: f64 },
}
