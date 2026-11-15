//! Opt-in electronic-loss straggling in BCA transport
//! ([`BcaConfig::straggling`], input `physics.straggling`).
//!
//! The Bohr variance is sampled per free flight and clamped to `[0, E]`
//! (`ion::stopping::straggling` module docs). The clamp depresses the sampled
//! variance whenever the per-flight standard deviation is comparable to the
//! per-flight mean loss, which for the constant free-path convention
//! (`lambda = N^(-1/3)`, here about 0.21 nm for carbon) it always is at MeV
//! energies: a probe of this setup with the constant convention measured
//! about 0.42 of the Bohr variance. The transmission test below therefore
//! uses the energy-dependent free path at `min_cm_angle_deg = 0.05`, whose
//! flights in this foil are tens of nanometres long, so a flight's Gaussian
//! sits far from the clamp and the sampled variance reproduces the Bohr
//! value to about two percent. That dependence is the documented truncation
//! bias of the model, not a defect of the sampling: the variance per unit
//! path length is the same either way.

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, EnergyBudget, MeanFreePath, SummaryTally, Thermal,
};
use lindhard::ion::crystal::debye::THETA_D_SI;
use lindhard::ion::stopping::bethe::BetheBloch;
use lindhard::ion::stopping::bethe::EffectiveCharge;
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::straggling::{variance_per_atom, StragglingModel};
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use lindhard::tally::{Binning, IonReport, IonTally, IonTallyConfig};
use lindhard::units::J_PER_EV;

mod common;
use common::{beam, crystal_si, si, with_threads};

const NM: f64 = 1e-9;

/// Carbon (tabulated density) with test energies, like the other tests:
/// the engine only needs some value, not data.
fn carbon() -> Material {
    let mut m = Material::from_atom_fractions(&[(6, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(6, 15.0).unwrap();
    m.set_surface_binding_energy_ev(6, 2.0).unwrap();
    m
}

/// A foil of `thickness_nm` of carbon and nothing else, so primaries leave
/// through the back face and are transmitted.
fn foil(thickness_nm: f64) -> Stack {
    Stack::new(vec![(carbon(), thickness_nm * NM)], None).unwrap()
}

/// Per-history transmitted energy (zero for a primary that never leaves
/// through the back face), plus the largest relative energy-budget residual.
#[derive(Default, Debug, PartialEq)]
struct Transmitted {
    energies: Vec<f64>,
    max_rel_residual: f64,
}

impl Transmitted {
    fn mean(&self) -> f64 {
        self.energies.iter().sum::<f64>() / self.energies.len() as f64
    }

    /// Sample variance of the transmitted energies.
    fn variance(&self) -> f64 {
        let mean = self.mean();
        let n = self.energies.len() as f64;
        self.energies
            .iter()
            .map(|e| (e - mean).powi(2))
            .sum::<f64>()
            / (n - 1.0)
    }
}

impl BcaTally for Transmitted {
    fn begin_history(&mut self, _index: u64) {
        self.energies.push(0.0);
    }

    fn end_history(&mut self, _index: u64, budget: &EnergyBudget) {
        if let Some(e) = self.energies.last_mut() {
            *e = budget.transmitted;
        }
        let r = (budget.residual() / budget.incident).abs();
        if r > self.max_rel_residual {
            self.max_rel_residual = r;
        }
    }

    fn merge(&mut self, other: Self) {
        self.energies.extend_from_slice(&other.energies);
        self.max_rel_residual = self.max_rel_residual.max(other.max_rel_residual);
    }
}

/// Bohr variance of the energy lost in a foil of `thickness_nm`, in eV²:
/// the per-atom variance times the areal atom density of the foil.
fn bohr_variance_ev2(
    ion: &Ion,
    target_z: u8,
    energy_ev: f64,
    n_density: f64,
    thickness_nm: f64,
) -> f64 {
    let omega2 = variance_per_atom(
        StragglingModel::Bohr,
        EffectiveCharge::Bare,
        ion,
        target_z,
        energy_ev,
    )
    .unwrap();
    omega2 * n_density * thickness_nm * NM / (J_PER_EV * J_PER_EV)
}

#[test]
fn transmitted_variance_matches_bohr_in_a_thin_foil() {
    let bb = BetheBloch::new();
    let thickness_nm = 100.0;
    let stack = foil(thickness_nm);
    let n_density = carbon().atom_number_density();
    let mut cfg = BcaConfig::new(1.0e5, 1.0);
    cfg.seed = 2026;
    cfg.straggling = Some(StragglingModel::Bohr);
    cfg.follow_recoils = false;
    // Flights of tens of nanometres (see the module docs above); the small
    // angle also keeps every sampled nuclear transfer subthreshold, so the
    // transmitted width is the electronic (Bohr) width.
    cfg.mean_free_path = MeanFreePath::EnergyDependent {
        min_cm_angle_rad: 0.05f64.to_radians(),
    };
    // 1 MeV proton on carbon: in Bohr's validity range by two orders of
    // magnitude (v = 141 v0 >> v0 Z1^(2/3) = v0) and in Bethe-Bloch's
    // (which covers the whole path: the mean loss is 6 keV of 1 MeV).
    // Recoils are not followed so the stopping model is never asked about a
    // low-energy carbon ion, which it does not cover.
    let beam = Beam::normal(Ion::proton(), 1.0e6, 20_000);
    let bca = Bca::new(beam, &stack, cfg, &bb, common::zbl_si_table()).unwrap();
    let t = bca.run(Transmitted::default).unwrap();

    // Every primary is transmitted (the C range at 1 MeV is tens of um).
    assert_eq!(t.energies.len(), 20_000);
    assert!(
        t.energies.iter().all(|&e| e > 0.0),
        "some primary not transmitted"
    );
    assert!(t.max_rel_residual < 1e-12, "{}", t.max_rel_residual);

    let measured = t.variance();
    let bohr = bohr_variance_ev2(&Ion::proton(), 6, 1.0e6, n_density, thickness_nm);
    // Statistical: the standard error of a variance estimated from 20_000
    // samples is sqrt(2/(N-1)) = 1% of the value. Systematic: the clamp
    // truncation and the nuclear path effects of the setup (the probe in the
    // module docs measured -1.7% net). The bound is stated at 10%.
    let rel = measured / bohr;
    assert!(
        (rel - 1.0).abs() < 0.10,
        "measured variance {measured:.3e} eV^2, Bohr {bohr:.3e} eV^2, ratio {rel}"
    );

    // The clamp must not inflate the mean loss either: with the truncation
    // inactive at this flight length the sampled mean loss matches the
    // deterministic mean loss to a few percent (probe measured +0.6%).
    let mut cfg0 = BcaConfig::new(1.0e5, 1.0);
    cfg0.seed = 2026;
    cfg0.follow_recoils = false;
    cfg0.mean_free_path = MeanFreePath::EnergyDependent {
        min_cm_angle_rad: 0.05f64.to_radians(),
    };
    let beam0 = Beam::normal(Ion::proton(), 1.0e6, 5_000);
    let bca0 = Bca::new(beam0, &stack, cfg0, &bb, common::zbl_si_table()).unwrap();
    let s0 = bca0.run(|| SummaryTally::new(NM, 10)).unwrap();
    let mean_loss_off = s0.budget.electronic_nonlocal / s0.histories as f64;
    let mean_loss_on = 1.0e6 - t.mean();
    let rel_loss = mean_loss_on / mean_loss_off - 1.0;
    assert!(
        rel_loss.abs() < 0.05,
        "mean loss with straggling {mean_loss_on:.1} eV, without {mean_loss_off:.1} eV \
         (ratio {rel_loss:+.4})"
    );
}

#[test]
fn straggling_run_is_bit_identical_across_thread_counts() {
    let bb = BetheBloch::new();
    let stack = foil(100.0);
    let mut cfg = BcaConfig::new(1.0e5, 1.0);
    cfg.seed = 2026;
    cfg.straggling = Some(StragglingModel::Bohr);
    cfg.follow_recoils = false;
    cfg.mean_free_path = MeanFreePath::EnergyDependent {
        min_cm_angle_rad: 0.05f64.to_radians(),
    };
    let beam = Beam::normal(Ion::proton(), 1.0e6, 2_000);
    let bca = Bca::new(beam, &stack, cfg, &bb, common::zbl_si_table()).unwrap();
    let r1 = with_threads(1, || bca.run(Transmitted::default).unwrap());
    for n in [2, 4] {
        let r = with_threads(n, || bca.run(Transmitted::default).unwrap());
        assert_eq!(r1, r, "straggling run differed at {n} threads");
    }
}

/// The complete serialised `IonReport` of a 10 keV As beam (40 ions, tilt
/// 7 deg, twist 22 deg) into a semi-infinite Si crystal, with `straggling`
/// and `thermal` switched and the run pinned to `threads` rayon threads —
/// the same whole-report comparison `crystal_thermal.rs` uses. The stopping
/// model is Lindhard-Scharff, valid here for the beam and every recoil
/// (Bethe-Bloch is not, at these energies).
fn crystal_report_json(
    straggling: Option<StragglingModel>,
    thermal: Option<Thermal>,
    threads: usize,
) -> String {
    let st = Stack::semi_infinite(si());
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = 11;
    cfg.straggling = straggling;
    let crystal = match thermal {
        Some(t) => crystal_si(7.0, 22.0).with_thermal(t),
        None => crystal_si(7.0, 22.0),
    };
    let ls = LindhardScharff::new();
    let bca = Bca::new(
        beam(33, 1.0e4, 7.0, 22.0, 40),
        &st,
        cfg,
        &ls,
        common::zbl_si_table(),
    )
    .unwrap()
    .with_crystal(crystal, &[0])
    .unwrap();
    let species = bca.species_z();
    let tc = IonTallyConfig {
        depth: Binning::new(0.0, 400.0 * NM, 400).unwrap(),
        lateral: Binning::new(-100.0 * NM, 100.0 * NM, 100).unwrap(),
        radial: Binning::new(0.0, 100.0 * NM, 50).unwrap(),
        escape_energy: Binning::new(0.0, 5.0e4, 50).unwrap(),
        escape_polar: Binning::new(0.0, std::f64::consts::FRAC_PI_2, 18).unwrap(),
    };
    let tally = with_threads(threads, || {
        bca.run(|| IonTally::new(&st, &species, tc).unwrap())
            .unwrap()
    });
    serde_json::to_string(&tally.report(false)).unwrap()
}

/// The lattice flight path samples straggling from the history's own
/// transport stream (the `electronic_nonlocal` call sites in `bca::crystal`,
/// which receive the transport `rng`, not the thermal one). Regression
/// guard for that plumbing: the complete report of a static and a thermal
/// Si lattice is bit-identical at 1, 2 and 4 threads, and switching
/// straggling on changes the report. A refactor of `crystal_step` that
/// silently stopped sampling fails the last assertions; one that reached
/// for the thermal stream on the static lattice fails loudly, because with
/// no thermal crystal there is no thermal stream to draw from. (A
/// wrong-stream swap on a *thermal* lattice is not caught here — the
/// thermal stream is as deterministic per history as the transport one —
/// that guard was deliberately rejected as a golden-tally baseline;
/// issue #830.)
#[test]
fn crystal_straggling_run_is_bit_identical_across_thread_counts() {
    let thermal = Thermal::new(300.0, THETA_D_SI);
    for (label, thermal) in [("static lattice", None), ("thermal lattice", Some(thermal))] {
        let on1 = crystal_report_json(Some(StragglingModel::Bohr), thermal, 1);
        // Histories ran and lost energy electronically: the comparisons are
        // not short-circuiting on an empty tally.
        let rep: IonReport = serde_json::from_str(&on1).unwrap();
        assert_eq!(rep.histories, 40, "{label}: histories ran");
        assert!(rep.budget.electronic_nonlocal > 0.0, "{label}: lost energy");
        for n in [2, 4] {
            assert_eq!(
                on1,
                crystal_report_json(Some(StragglingModel::Bohr), thermal, n),
                "{label}: straggling run differed at {n} threads"
            );
        }
        assert_ne!(
            on1,
            crystal_report_json(None, thermal, 1),
            "{label}: straggling changed nothing"
        );
    }
}

/// Where the per-flight standard deviation is comparable to the per-flight
/// mean loss (short flights, low energy) the clamp binds on most flights.
/// The truncated loss still never leaves the ion's energy: every history
/// conserves energy exactly, and the sampled losses stay non-negative.
#[test]
fn clamp_binds_on_short_flights_and_energy_is_conserved() {
    let ls = LindhardScharff::new();
    let stack = foil(2.0);
    let e0 = 5.0e3;
    let mut cfg = BcaConfig::new(50.0, 1.0);
    cfg.seed = 7;
    cfg.straggling = Some(StragglingModel::Bohr);
    // Constant free path (lambda = N^(-1/3), about 0.2 nm): at 5 keV the
    // mean loss per flight is under 1 eV while the Bohr standard deviation
    // per flight is tens of eV, so the lower clamp binds almost always.
    // Lindhard-Scharff is valid here for the proton and every recoil energy.
    cfg.mean_free_path = MeanFreePath::Constant;
    let beam = Beam::normal(Ion::proton(), e0, 2_000);
    let bca = Bca::new(beam, &stack, cfg, &ls, common::zbl_si_table()).unwrap();
    let t = bca.run(Transmitted::default).unwrap();
    assert!(
        t.max_rel_residual < 1e-12,
        "clamped straggling broke conservation: {}",
        t.max_rel_residual
    );
    // Transmitted accumulates only the exit energies of the primary and the
    // recoils, so `>= 0` holds by construction and does not test the clamp;
    // the bound that bites is the upper one. A negative sampled loss would
    // be an energy gain, and a history that gained energy would leave with
    // more than the beam put in (the relative slack covers float
    // accumulation over the additions).
    assert!(
        t.energies
            .iter()
            .all(|&e| e.is_finite() && e >= 0.0 && e <= e0 * (1.0 + 1e-9)),
        "a history left with more than the beam energy: {:?}",
        t.energies.iter().filter(|&&e| e > e0).collect::<Vec<_>>()
    );
}
