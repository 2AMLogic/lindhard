//! The photon-borne tally: energy deposited by photon interactions per region
//! (split from electron-borne dose, which this tally never contains), photon
//! escape by face, and the energy balance of a photon run. Collected by
//! [`crate::photon::PhotonTransport`].
//!
//! Every accumulator merges deterministically: counts are exact sums and
//! energies are added in the order [`crate::rng::run_particles`] merges chunks
//! (chunk order), so a run is bit-identical on any thread count.
//!
//! Energies are in eV. A photon-borne deposit is the energy handed to
//! electrons that are *not* transported (the photoelectron rule in the
//! [`crate::photon`] module docs): the photoelectric deposit is the whole
//! photon energy, the Compton deposit is `E - E'`. The laterally infinite
//! [`crate::geometry::Stack`] has no mass per layer, so the tally keeps energy
//! (eV, total over the run) rather than Gy; divide by `n_histories` for the
//! deposit per incident photon, and by the layer mass for dose.

use serde::Serialize;

use crate::geometry::Face;

/// Photons that left through one face.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize)]
pub struct PhotonFaceEscape {
    /// Number of photons.
    pub n: u64,
    /// Of those, photons that never interacted (the primary, unscattered
    /// beam).
    pub n_uncollided: u64,
    /// Total energy carried out, eV.
    pub energy_ev: f64,
}

impl PhotonFaceEscape {
    fn merge(&mut self, o: &Self) {
        self.n += o.n;
        self.n_uncollided += o.n_uncollided;
        self.energy_ev += o.energy_ev;
    }
}

/// Mergeable photon tally. See the [module docs](self).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PhotonTally {
    /// Histories run.
    pub n_histories: u64,
    /// Sum of the source photon energies, eV.
    pub incident_ev: f64,
    /// Per region: energy of photoelectric absorptions (the whole photon
    /// energy, deposited locally), eV.
    pub photoelectric_ev: Vec<f64>,
    /// Per region: energy given to Compton recoil electrons (`E - E'`,
    /// deposited locally), eV.
    pub compton_ev: Vec<f64>,
    /// Per region: energy of photons that fell below the lowest table energy
    /// and were absorbed on the spot, eV.
    pub cutoff_ev: Vec<f64>,
    /// Escape through the front face.
    pub front: PhotonFaceEscape,
    /// Escape through the back face.
    pub back: PhotonFaceEscape,
    /// Escape through a lateral face of a voxel grid.
    pub side: PhotonFaceEscape,
    /// Number of source photons that began outside the target.
    pub n_missed: u64,
    /// Energy of those, eV (they never enter the target).
    pub missed_ev: f64,
    /// Photoelectric absorptions.
    pub n_photoelectric: u64,
    /// Incoherent (Compton) scatterings.
    pub n_incoherent: u64,
    /// Coherent (Rayleigh) scatterings.
    pub n_coherent: u64,
    /// Photons absorbed below the table energy.
    pub n_cutoff: u64,
}

/// `incident = deposited + escaped + missed`, with the residual.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct PhotonEnergyBalance {
    /// Source energy, eV.
    pub incident_ev: f64,
    /// Photon-borne deposits over all regions, eV.
    pub deposited_ev: f64,
    /// Energy carried out through all faces, eV.
    pub escaped_ev: f64,
    /// Energy of photons that never entered, eV.
    pub missed_ev: f64,
    /// `incident - deposited - escaped - missed`, eV (rounding only).
    pub residual_ev: f64,
}

impl PhotonEnergyBalance {
    /// `|residual| / incident` (0 for an empty run).
    pub fn relative_residual(&self) -> f64 {
        if self.incident_ev > 0.0 {
            self.residual_ev.abs() / self.incident_ev
        } else {
            0.0
        }
    }
}

impl PhotonTally {
    /// An empty tally for a geometry with `n_regions` regions.
    pub fn new(n_regions: usize) -> Self {
        Self {
            n_histories: 0,
            incident_ev: 0.0,
            photoelectric_ev: vec![0.0; n_regions],
            compton_ev: vec![0.0; n_regions],
            cutoff_ev: vec![0.0; n_regions],
            front: PhotonFaceEscape::default(),
            back: PhotonFaceEscape::default(),
            side: PhotonFaceEscape::default(),
            n_missed: 0,
            missed_ev: 0.0,
            n_photoelectric: 0,
            n_incoherent: 0,
            n_coherent: 0,
            n_cutoff: 0,
        }
    }

    /// Number of regions.
    pub fn n_regions(&self) -> usize {
        self.photoelectric_ev.len()
    }

    /// Photon-borne energy deposited in `region`, eV (all three channels).
    pub fn deposited_ev(&self, region: usize) -> f64 {
        self.photoelectric_ev[region] + self.compton_ev[region] + self.cutoff_ev[region]
    }

    /// Photon-borne energy deposited over all regions, eV.
    pub fn total_deposited_ev(&self) -> f64 {
        (0..self.n_regions()).map(|r| self.deposited_ev(r)).sum()
    }

    /// Record a photon leaving through `face` with `energy_ev`.
    pub(crate) fn record_escape(&mut self, face: Face, energy_ev: f64, uncollided: bool) {
        let f = match face {
            Face::Front => &mut self.front,
            Face::Back => &mut self.back,
            Face::Side => &mut self.side,
        };
        f.n += 1;
        f.n_uncollided += u64::from(uncollided);
        f.energy_ev += energy_ev;
    }

    /// The energy balance of the run.
    pub fn energy_balance(&self) -> PhotonEnergyBalance {
        let deposited_ev = self.total_deposited_ev();
        let escaped_ev = self.front.energy_ev + self.back.energy_ev + self.side.energy_ev;
        PhotonEnergyBalance {
            incident_ev: self.incident_ev,
            deposited_ev,
            escaped_ev,
            missed_ev: self.missed_ev,
            residual_ev: self.incident_ev - deposited_ev - escaped_ev - self.missed_ev,
        }
    }

    /// Fold `other` (same number of regions) into `self`.
    ///
    /// # Panics
    ///
    /// If the region counts differ.
    pub fn merge(&mut self, other: Self) {
        assert_eq!(self.n_regions(), other.n_regions(), "region count mismatch");
        self.n_histories += other.n_histories;
        self.incident_ev += other.incident_ev;
        for (a, b) in self
            .photoelectric_ev
            .iter_mut()
            .zip(&other.photoelectric_ev)
        {
            *a += b;
        }
        for (a, b) in self.compton_ev.iter_mut().zip(&other.compton_ev) {
            *a += b;
        }
        for (a, b) in self.cutoff_ev.iter_mut().zip(&other.cutoff_ev) {
            *a += b;
        }
        self.front.merge(&other.front);
        self.back.merge(&other.back);
        self.side.merge(&other.side);
        self.n_missed += other.n_missed;
        self.missed_ev += other.missed_ev;
        self.n_photoelectric += other.n_photoelectric;
        self.n_incoherent += other.n_incoherent;
        self.n_coherent += other.n_coherent;
        self.n_cutoff += other.n_cutoff;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_adds_everything_and_balance_closes() {
        let mut a = PhotonTally::new(2);
        a.n_histories = 2;
        a.incident_ev = 200.0;
        a.photoelectric_ev[0] = 60.0;
        a.compton_ev[1] = 20.0;
        a.record_escape(Face::Back, 120.0, false);
        let mut b = PhotonTally::new(2);
        b.n_histories = 1;
        b.incident_ev = 100.0;
        b.cutoff_ev[1] = 10.0;
        b.record_escape(Face::Front, 90.0, true);
        a.merge(b);
        assert_eq!(a.n_histories, 3);
        assert_eq!(a.deposited_ev(1), 30.0);
        assert_eq!(a.front.n_uncollided, 1);
        let bal = a.energy_balance();
        assert_eq!(bal.deposited_ev, 90.0);
        assert_eq!(bal.escaped_ev, 210.0);
        assert_eq!(bal.residual_ev, 0.0);
        assert_eq!(bal.relative_residual(), 0.0);
    }

    #[test]
    #[should_panic(expected = "region count")]
    fn merge_refuses_a_different_geometry() {
        PhotonTally::new(1).merge(PhotonTally::new(2));
    }
}
