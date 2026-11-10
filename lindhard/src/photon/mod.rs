//! Photon transport: photons through a [`Geometry`](crate::geometry::Geometry)
//! with photoelectric absorption, incoherent (Compton) and coherent (Rayleigh)
//! scattering, from EPDL2017 cross sections.
//!
//! This is a stand-alone channel: a photon source is configured explicitly
//! ([`PhotonSource`]) and nothing in the electron or ion engines changes when
//! none is, so every existing output stays bit-identical.
//!
//! # Physics
//!
//! * Cross sections per atom: photoelectric (total photoionization),
//!   incoherent and coherent, from the committed EPDL2017 table
//!   ([`PhotonTable::epdl2017`], `Z = 1..=92`, 1 keV to 1 MeV; see
//!   [`data`] and `docs/data-provenance.md`). A material's macroscopic cross
//!   section is the sum over its elements of number density times the atomic
//!   value (independent atoms; no molecular or crystal effects).
//! * Free flight: exact, as in the electron transport's layer crossings: draw
//!   the path from the current material's total cross section, truncate at the
//!   first face (or surface), and redraw in the next region.
//! * Incoherent scattering: angle from `KN(E, θ) · S(x, Z)` (Klein-Nishina
//!   times the EPDL scattering function, so the same binding correction as the
//!   tabulated cross section; Kahn's rejection method for the Klein-Nishina factor),
//!   scattered energy from Compton kinematics `E' = E / (1 + α(1 - cos θ))`. No
//!   Doppler broadening: `E'` is the free-electron value.
//! * Coherent scattering: angle from `(1 + cos²θ)/2 · F(x, Z)²`, no energy
//!   change. The anomalous scattering factors are not used.
//! * Pair production is out of scope; tables stop at 1 MeV, below its
//!   threshold, and a source above that is refused.
//!
//! # Photoelectron rule
//!
//! **A photoelectric absorption deposits the whole photon energy at the
//! interaction point, and the electron of a Compton scattering deposits `E - E'`
//! there.** No electron (photo-, Compton recoil) and no fluorescence photon or
//! Auger electron is transported; the binding energy is not separated from the
//! kinetic energy. The consequences for the energy balance:
//!
//! * `incident = deposited + escaped` holds exactly (to rounding): every
//!   photon-borne joule is either deposited locally or leaves as a photon.
//! * The deposit is photon-borne dose only; it is the electron-borne dose that
//!   a later child (electron hand-back, relaxation) will split off it. Near a
//!   layer boundary this over-counts the layer that holds the interaction when
//!   the electron range reaches the neighbour, and it never counts energy that
//!   a fluorescence photon would carry away (the K-shell fluorescence of a
//!   heavy element above its edge can be a large fraction of the absorbed
//!   energy; that energy is counted as deposited here).
//! * A photon whose energy falls below the table's lowest energy (1 keV) after
//!   scattering is absorbed on the spot (tallied separately as `cutoff`).
//!
//! # Determinism
//!
//! Histories draw from the counter-based streams of [`crate::rng`]; the tally
//! merges in chunk order, so a run is bit-identical on any thread count
//! (`lindhard/tests/photon_determinism.rs`).

pub mod data;
mod sample;
mod transport;

pub use data::{AtomPhotonData, Curve, PhotonDataError, PhotonTable};
pub use transport::{PhotonSource, PhotonTransport, PhotonTransportError};
