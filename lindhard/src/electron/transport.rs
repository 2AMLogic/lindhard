//! Event-by-event electron transport in a layered stack, on precomputed
//! cross-section tables.
//!
//! This is the low-energy electron Monte Carlo loop of Kieft and Bosch,
//! J. Phys. D: Appl. Phys. 41, 215310 (2008), doi:10.1088/0022-3727/41/21/215310
//! (the scheme on which Nebula is built): an electron at `(r, direction, E)`
//! in a layer flies a free path drawn from the total inverse mean free path
//! (the sum of the elastic and inelastic [`CrossSectionTable`] rates of the
//! layer at `E`), then picks the elastic or the inelastic channel with
//! probability proportional to the two rates, samples the polar angle `θ` or
//! the energy loss `W` from that channel's inverse CDF, and updates its state.
//! The history ends when `E` falls below the cutoff or the electron leaves the
//! target. This implementation is written from the paper; no Nebula source is
//! ported (so there is no `THIRD_PARTY_LICENSES.md` entry for it).
//!
//! # Model
//!
//! - **Free flight.** `s = -ln(1 - u) / (λ_el⁻¹ + λ_inel⁻¹)`, `u` uniform in
//!   `[0, 1)`. Rates are interpolated linearly in energy on each table's own
//!   grid and held constant beyond its ends.
//! - **Layer boundaries are crossed exactly.** If the drawn path is longer
//!   than the distance to the layer face along the flight direction, the
//!   electron moves to the face (its depth is set to the face value exactly),
//!   enters the neighbouring layer and **redraws** the path with that layer's
//!   rates. This is exact because the free path is memoryless.
//! - **Elastic collision.** `θ` from the table, azimuth uniform on `[0, 2π)`;
//!   the energy is unchanged (no recoil energy).
//! - **Inelastic collision.** `W` from the table (clamped to `E`); the energy
//!   becomes `E - W`. The flight direction is unchanged, and no secondary is
//!   made: both belong to later work, see [Hooks](#hooks).
//! - **Inverse CDFs between grid energies** are interpolated linearly in
//!   energy between the two bracketing rows, at the same cumulative
//!   probability (a zero-rate bracketing row is ignored).
//!
//! # Random draw order
//!
//! The order of draws from the per-history stream is part of the contract (the
//! tests replay it). Per flight with a nonzero total rate: one uniform for the
//! path. If a collision happens, one for the channel, then `θ` and the
//! azimuth (elastic) or `W` (inelastic). A path drawn for a flight that ends at
//! a layer face is discarded and a new one is drawn in the next layer. A
//! flight with zero total rate draws nothing.
//!
//! # Reproducibility
//!
//! [`Transport::run`] goes through [`crate::rng::run_particles`]: history `i`
//! uses the stream `(seed, i)`, and chunk tallies are merged in chunk order,
//! so a tally whose [`ElectronTally::merge`] is order-faithful (a sum, or a
//! concatenation) is bit-identical at any thread count.
//!
//! # Hooks
//!
//! The [`ElectronTally`] trait mirrors [`crate::ion::bca::BcaTally`]: every
//! event hook has a no-op default. The extension points for later work are
//! observer hooks that already carry what those models need:
//! [`ElectronTally::interface`] (interface refraction and transmission
//! models), [`ElectronTally::inelastic`] with the energy loss `W` and the
//! post-event state (secondary electron generation), and
//! [`ElectronTally::escaped`] (surface barrier). Phonon and polaron channels
//! would add a third channel next to the two tables.

use serde::Serialize;

use rand_core::Rng;

use crate::electron::data::{CrossSectionTable, SamplingAxis};
use crate::geometry::Stack;
use crate::rng::run_particles;

/// Errors from building a [`Transport`] or starting a run.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum TransportError {
    /// The number of table pairs differs from the number of layers.
    #[error("the stack has {layers} layers but {tables} table pairs were given")]
    LayerCount {
        /// Layers in the stack.
        layers: usize,
        /// Table pairs given.
        tables: usize,
    },
    /// A table has the wrong sampling axis for its slot.
    #[error("layer {layer}: the {slot} table samples {found:?}, expected {expected:?}")]
    WrongAxis {
        /// Layer index.
        layer: usize,
        /// `"elastic"` or `"inelastic"`.
        slot: &'static str,
        /// Axis found.
        found: SamplingAxis,
        /// Axis expected.
        expected: SamplingAxis,
    },
    /// A configuration or primary value is out of range.
    #[error("invalid {what}: {why}")]
    Invalid {
        /// The quantity.
        what: &'static str,
        /// Why it was rejected.
        why: String,
    },
}

fn invalid<T>(what: &'static str, why: impl Into<String>) -> Result<T, TransportError> {
    Err(TransportError::Invalid {
        what,
        why: why.into(),
    })
}

/// What happens when an electron reaches a face of the target. There is no
/// surface barrier or refraction yet: an electron that crosses a face leaves
/// at its current energy and direction and never returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum EscapeRule {
    /// Both the front face (`x = 0`) and the back face of a finite stack
    /// release the electron ([`Fate::Escaped`] with [`Face::Front`] or
    /// [`Face::Back`]).
    BothFaces,
    /// The front face releases the electron; the back face of a finite stack
    /// absorbs it ([`Fate::Absorbed`]).
    FrontOnly,
}

/// Which face of the target an electron left through.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Face {
    /// The front surface, `x = 0`.
    Front,
    /// The back face of a finite stack.
    Back,
}

/// How a history ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fate {
    /// The energy fell below the cutoff inside the target.
    Stopped,
    /// The electron left the target through a face.
    Escaped(Face),
    /// The electron reached the back face under [`EscapeRule::FrontOnly`].
    Absorbed,
    /// The electron has no interaction available and no face to reach (zero
    /// total rate heading parallel to the faces, or away from every face).
    Trapped,
    /// The history hit [`TransportConfig::max_events`] and was cut off.
    EventCap,
}

/// Run configuration: the energy cutoff and the escape condition are recorded
/// in the [`RunMetadata`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TransportConfig {
    /// A history stops when the energy falls below this, eV.
    pub cutoff_ev: f64,
    /// What a face of the target does to an electron reaching it.
    pub escape_rule: EscapeRule,
    /// Safety cap on collisions per history. A history that reaches it ends
    /// with [`Fate::EventCap`].
    pub max_events: u64,
}

impl TransportConfig {
    /// A configuration with the given cutoff, [`EscapeRule::BothFaces`] and a
    /// cap of 10 million collisions per history.
    pub fn new(cutoff_ev: f64) -> Self {
        Self {
            cutoff_ev,
            escape_rule: EscapeRule::BothFaces,
            max_events: 10_000_000,
        }
    }
}

/// The tables of one layer.
#[derive(Debug, Clone, PartialEq)]
pub struct LayerTables {
    /// Elastic table ([`SamplingAxis::ElasticPolarAngle`]).
    pub elastic: CrossSectionTable,
    /// Inelastic table ([`SamplingAxis::InelasticEnergyLoss`]).
    pub inelastic: CrossSectionTable,
}

/// A primary electron: it starts on the front face (`x = 0`, `y = z = 0`) in
/// the first layer.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Primary {
    /// Kinetic energy, eV; must exceed the cutoff.
    pub energy_ev: f64,
    /// Direction of flight (any nonzero finite vector; normalised on use).
    /// Normal incidence is `[1, 0, 0]`; a negative `x` component starts the
    /// electron moving out of the target.
    pub direction: [f64; 3],
}

impl Primary {
    /// Normal incidence at `energy_ev`.
    pub fn normal(energy_ev: f64) -> Self {
        Self {
            energy_ev,
            direction: [1.0, 0.0, 0.0],
        }
    }
}

/// The state of an electron after an event. Positions are `[x, y, z]` in
/// metres, `x` the depth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElectronState {
    /// Position, m.
    pub pos: [f64; 3],
    /// Unit direction of flight.
    pub dir: [f64; 3],
    /// Kinetic energy, eV.
    pub energy_ev: f64,
    /// Index of the layer the electron is in (or was last in, after an
    /// escape).
    pub layer: usize,
}

/// Event hooks called by [`Transport`]. All default to no-ops except `merge`.
/// The `ElectronState` passed is the state **after** the event.
pub trait ElectronTally: Send {
    /// A primary history starts.
    fn begin_history(&mut self, _index: u64, _start: &ElectronState) {}

    /// A flight segment of `length_m` ended at `end.pos` (at a collision, a
    /// layer face or the target face). It started at `from` with energy
    /// `end.energy_ev` (the energy does not change in flight).
    fn step(&mut self, _from: [f64; 3], _end: &ElectronState, _length_m: f64) {}

    /// An elastic collision deflected the electron by `theta` rad.
    fn elastic(&mut self, _after: &ElectronState, _theta: f64) {}

    /// An inelastic collision took `w_ev` from the electron. Secondary
    /// electron generation attaches here.
    fn inelastic(&mut self, _after: &ElectronState, _w_ev: f64) {}

    /// The electron crossed from layer `from_layer` into `to_layer`; `at` is
    /// its state on the face, inside the new layer. Interface refraction and
    /// transmission models attach here.
    fn interface(&mut self, _at: &ElectronState, _from_layer: usize, _to_layer: usize) {}

    /// The electron fell below the cutoff (or was trapped) at `at`.
    fn stopped(&mut self, _at: &ElectronState) {}

    /// The electron left the target through `face`.
    fn escaped(&mut self, _at: &ElectronState, _face: Face) {}

    /// The electron was absorbed at the back face ([`EscapeRule::FrontOnly`]).
    fn absorbed(&mut self, _at: &ElectronState) {}

    /// The history ended with `fate`.
    fn end_history(&mut self, _index: u64, _fate: Fate) {}

    /// Fold `other` (a later chunk) into `self`.
    fn merge(&mut self, other: Self)
    where
        Self: Sized;
}

/// A minimal summary tally: counts, path length and energy bookkeeping.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SummaryTally {
    /// Histories run.
    pub histories: u64,
    /// Histories that stopped below the cutoff.
    pub stopped: u64,
    /// Electrons that left through the front face.
    pub escaped_front: u64,
    /// Electrons that left through the back face.
    pub escaped_back: u64,
    /// Electrons absorbed at the back face.
    pub absorbed: u64,
    /// Trapped electrons.
    pub trapped: u64,
    /// Histories cut off by the event cap.
    pub event_capped: u64,
    /// Elastic collisions.
    pub elastic_events: u64,
    /// Inelastic collisions.
    pub inelastic_events: u64,
    /// Layer-face crossings between layers.
    pub interface_crossings: u64,
    /// Total flight path, m.
    pub path_m: f64,
    /// Total energy lost in inelastic collisions, eV.
    pub inelastic_loss_ev: f64,
    /// Energy carried out of the target by escaped electrons, eV.
    pub escaped_energy_ev: f64,
    /// Energy of stopped electrons when they stopped, eV.
    pub rest_energy_ev: f64,
}

impl ElectronTally for SummaryTally {
    fn step(&mut self, _from: [f64; 3], _end: &ElectronState, length_m: f64) {
        self.path_m += length_m;
    }
    fn elastic(&mut self, _after: &ElectronState, _theta: f64) {
        self.elastic_events += 1;
    }
    fn inelastic(&mut self, _after: &ElectronState, w_ev: f64) {
        self.inelastic_events += 1;
        self.inelastic_loss_ev += w_ev;
    }
    fn interface(&mut self, _at: &ElectronState, _from: usize, _to: usize) {
        self.interface_crossings += 1;
    }
    fn stopped(&mut self, at: &ElectronState) {
        self.rest_energy_ev += at.energy_ev;
    }
    fn escaped(&mut self, at: &ElectronState, _face: Face) {
        self.escaped_energy_ev += at.energy_ev;
    }
    fn end_history(&mut self, _index: u64, fate: Fate) {
        self.histories += 1;
        match fate {
            Fate::Stopped => self.stopped += 1,
            Fate::Escaped(Face::Front) => self.escaped_front += 1,
            Fate::Escaped(Face::Back) => self.escaped_back += 1,
            Fate::Absorbed => self.absorbed += 1,
            Fate::Trapped => self.trapped += 1,
            Fate::EventCap => self.event_capped += 1,
        }
    }
    fn merge(&mut self, o: Self) {
        self.histories += o.histories;
        self.stopped += o.stopped;
        self.escaped_front += o.escaped_front;
        self.escaped_back += o.escaped_back;
        self.absorbed += o.absorbed;
        self.trapped += o.trapped;
        self.event_capped += o.event_capped;
        self.elastic_events += o.elastic_events;
        self.inelastic_events += o.inelastic_events;
        self.interface_crossings += o.interface_crossings;
        self.path_m += o.path_m;
        self.inelastic_loss_ev += o.inelastic_loss_ev;
        self.escaped_energy_ev += o.escaped_energy_ev;
        self.rest_energy_ev += o.rest_energy_ev;
    }
}

/// Identity of the tables used in one layer, for [`RunMetadata`].
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LayerMetadata {
    /// Layer index (0 = front).
    pub index: usize,
    /// Depth of the front face, m.
    pub front_m: f64,
    /// Depth of the back face, m; `None` for a semi-infinite substrate.
    pub back_m: Option<f64>,
    /// Elastic table model identity.
    pub elastic_model: String,
    /// Elastic table provenance.
    pub elastic_provenance: String,
    /// Inelastic table model identity.
    pub inelastic_model: String,
    /// Inelastic table provenance.
    pub inelastic_provenance: String,
}

/// What a run was configured with, for output metadata.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct RunMetadata {
    /// Energy cutoff, eV.
    pub cutoff_ev: f64,
    /// Escape rule at the target faces.
    pub escape_rule: EscapeRule,
    /// Collision cap per history.
    pub max_events: u64,
    /// Run seed.
    pub seed: u64,
    /// Primaries run.
    pub n_histories: u64,
    /// Chunk size of the parallel driver.
    pub chunk_size: u64,
    /// The primary electron.
    pub primary: Primary,
    /// Per-layer table identities.
    pub layers: Vec<LayerMetadata>,
}

/// The result of [`Transport::run`].
#[derive(Debug, Clone, PartialEq)]
pub struct TransportRun<T> {
    /// The merged tally.
    pub tally: T,
    /// The run configuration.
    pub metadata: RunMetadata,
}

/// The electron transport engine: a stack, the tables of each layer and a
/// configuration.
#[derive(Debug, Clone)]
pub struct Transport {
    stack: Stack,
    tables: Vec<LayerTables>,
    config: TransportConfig,
}

impl Transport {
    /// Check and assemble. `tables[i]` belongs to `stack.layers()[i]`.
    pub fn new(
        stack: Stack,
        tables: Vec<LayerTables>,
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        if tables.len() != stack.layers().len() {
            return Err(TransportError::LayerCount {
                layers: stack.layers().len(),
                tables: tables.len(),
            });
        }
        for (layer, t) in tables.iter().enumerate() {
            for (slot, table, expected) in [
                ("elastic", &t.elastic, SamplingAxis::ElasticPolarAngle),
                ("inelastic", &t.inelastic, SamplingAxis::InelasticEnergyLoss),
            ] {
                if table.axis() != expected {
                    return Err(TransportError::WrongAxis {
                        layer,
                        slot,
                        found: table.axis(),
                        expected,
                    });
                }
            }
        }
        if !(config.cutoff_ev.is_finite() && config.cutoff_ev > 0.0) {
            return invalid("cutoff", "must be finite and positive");
        }
        if config.max_events == 0 {
            return invalid("max_events", "must be at least 1");
        }
        Ok(Self {
            stack,
            tables,
            config,
        })
    }

    /// The configuration.
    pub fn config(&self) -> &TransportConfig {
        &self.config
    }

    /// The target.
    pub fn stack(&self) -> &Stack {
        &self.stack
    }

    fn check_primary(&self, p: &Primary) -> Result<[f64; 3], TransportError> {
        if !(p.energy_ev.is_finite() && p.energy_ev > self.config.cutoff_ev) {
            return invalid(
                "primary energy",
                format!("{} eV must exceed the cutoff", p.energy_ev),
            );
        }
        let n = (p.direction[0].powi(2) + p.direction[1].powi(2) + p.direction[2].powi(2)).sqrt();
        if !(n.is_finite() && n > 0.0) {
            return invalid("primary direction", "must be finite and nonzero");
        }
        Ok([p.direction[0] / n, p.direction[1] / n, p.direction[2] / n])
    }

    /// The metadata of a run with these arguments.
    pub fn metadata(
        &self,
        seed: u64,
        n_histories: u64,
        chunk_size: u64,
        primary: &Primary,
    ) -> RunMetadata {
        RunMetadata {
            cutoff_ev: self.config.cutoff_ev,
            escape_rule: self.config.escape_rule,
            max_events: self.config.max_events,
            seed,
            n_histories,
            chunk_size,
            primary: *primary,
            layers: self
                .stack
                .layers()
                .iter()
                .zip(&self.tables)
                .enumerate()
                .map(|(index, (l, t))| LayerMetadata {
                    index,
                    front_m: l.front_m(),
                    back_m: l.back_m().is_finite().then_some(l.back_m()),
                    elastic_model: t.elastic.model().to_string(),
                    elastic_provenance: t.elastic.provenance().to_string(),
                    inelastic_model: t.inelastic.model().to_string(),
                    inelastic_provenance: t.inelastic.provenance().to_string(),
                })
                .collect(),
        }
    }

    /// Run `n_histories` primaries in parallel (on the current rayon pool),
    /// bit-identically at any thread count. `new_tally` makes an empty tally;
    /// one is made per chunk of `chunk_size` histories.
    pub fn run<T, N>(
        &self,
        seed: u64,
        n_histories: u64,
        chunk_size: u64,
        primary: &Primary,
        new_tally: N,
    ) -> Result<TransportRun<T>, TransportError>
    where
        T: ElectronTally,
        N: Fn() -> T + Sync,
    {
        self.check_primary(primary)?;
        let tally = run_particles(
            seed,
            n_histories,
            chunk_size,
            new_tally,
            |t, rng, i| {
                self.history(t, rng, i, primary)
                    .expect("primary validated before the run");
            },
            |total, part| total.merge(part),
        );
        Ok(TransportRun {
            tally,
            metadata: self.metadata(seed, n_histories, chunk_size, primary),
        })
    }

    /// Simulate one primary with the given random stream.
    pub fn history<T, R>(
        &self,
        tally: &mut T,
        rng: &mut R,
        index: u64,
        primary: &Primary,
    ) -> Result<Fate, TransportError>
    where
        T: ElectronTally,
        R: Rng,
    {
        let dir = self.check_primary(primary)?;
        let mut st = ElectronState {
            pos: [0.0; 3],
            dir,
            energy_ev: primary.energy_ev,
            layer: 0,
        };
        tally.begin_history(index, &st);
        let fate = self.follow(tally, rng, &mut st);
        tally.end_history(index, fate);
        Ok(fate)
    }

    fn follow<T: ElectronTally, R: Rng>(
        &self,
        tally: &mut T,
        rng: &mut R,
        st: &mut ElectronState,
    ) -> Fate {
        let layers = self.stack.layers();
        let mut events = 0u64;
        loop {
            if events >= self.config.max_events {
                return Fate::EventCap;
            }
            let tabs = &self.tables[st.layer];
            let (el, el_at) = rate(&tabs.elastic, st.energy_ev);
            let (inel, inel_at) = rate(&tabs.inelastic, st.energy_ev);
            let total = el + inel;
            let layer = &layers[st.layer];
            let mu = st.dir[0];
            let to_face = if mu > 0.0 {
                ((layer.back_m() - st.pos[0]) / mu).max(0.0)
            } else if mu < 0.0 {
                ((layer.front_m() - st.pos[0]) / mu).max(0.0)
            } else {
                f64::INFINITY
            };
            let free = if total > 0.0 {
                -(1.0 - uniform(rng)).ln() / total
            } else {
                f64::INFINITY
            };
            if free.is_infinite() && to_face.is_infinite() {
                tally.stopped(st);
                return Fate::Trapped;
            }
            let from = st.pos;
            if free < to_face {
                for k in 0..3 {
                    st.pos[k] += free * st.dir[k];
                }
                tally.step(from, st, free);
                events += 1;
                if uniform(rng) * total < el {
                    let u = uniform(rng);
                    let theta = sample(&tabs.elastic, el_at, u).clamp(0.0, std::f64::consts::PI);
                    let phi = std::f64::consts::TAU * uniform(rng);
                    st.dir = deflect(st.dir, theta, phi);
                    tally.elastic(st, theta);
                } else {
                    let u = uniform(rng);
                    let w = sample(&tabs.inelastic, inel_at, u).clamp(0.0, st.energy_ev);
                    st.energy_ev -= w;
                    tally.inelastic(st, w);
                    if st.energy_ev < self.config.cutoff_ev {
                        tally.stopped(st);
                        return Fate::Stopped;
                    }
                }
                continue;
            }
            // Reach the face exactly; the path is redrawn in the next layer.
            let face_x = if mu > 0.0 {
                layer.back_m()
            } else {
                layer.front_m()
            };
            for k in 0..3 {
                st.pos[k] += to_face * st.dir[k];
            }
            st.pos[0] = face_x;
            tally.step(from, st, to_face);
            let last = layers.len() - 1;
            if mu > 0.0 {
                if st.layer < last {
                    st.layer += 1;
                    tally.interface(st, st.layer - 1, st.layer);
                } else if self.config.escape_rule == EscapeRule::BothFaces {
                    tally.escaped(st, Face::Back);
                    return Fate::Escaped(Face::Back);
                } else {
                    tally.absorbed(st);
                    return Fate::Absorbed;
                }
            } else if st.layer > 0 {
                st.layer -= 1;
                tally.interface(st, st.layer + 1, st.layer);
            } else {
                tally.escaped(st, Face::Front);
                return Fate::Escaped(Face::Front);
            }
        }
    }
}

/// Uniform on `[0, 1)` with 53 random bits.
fn uniform<R: Rng>(rng: &mut R) -> f64 {
    (rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// Where an energy sits on a table's grid: lower row and weight of the upper.
#[derive(Clone, Copy)]
struct GridPos {
    i: usize,
    f: f64,
}

/// Inverse mean free path at `e` (linear in energy, constant beyond the ends)
/// and the grid position.
fn rate(t: &CrossSectionTable, e: f64) -> (f64, GridPos) {
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
fn sample(t: &CrossSectionTable, at: GridPos, u: f64) -> f64 {
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
fn deflect(d: [f64; 3], theta: f64, phi: f64) -> [f64; 3] {
    let helper = if d[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let e1 = normalize(cross(helper, d));
    let e2 = cross(d, e1);
    let (st, ct) = theta.sin_cos();
    let (sp, cp) = phi.sin_cos();
    normalize([
        ct * d[0] + st * (cp * e1[0] + sp * e2[0]),
        ct * d[1] + st * (cp * e1[1] + sp * e2[1]),
        ct * d[2] + st * (cp * e1[2] + sp * e2[2]),
    ])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
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
