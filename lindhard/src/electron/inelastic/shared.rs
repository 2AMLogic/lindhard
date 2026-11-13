//! Argument checks and unit conversions shared by the inelastic models.
//!
//! The single-pole Penn model ([`super::penn`]), the full Penn algorithm
//! ([`super::full_penn`]) and the Mermin-ELF model ([`super::mermin`]) take
//! the same arguments — an electron energy, a Fermi energy, a relative
//! integration tolerance — and all return the inverse IMFP and the stopping
//! power as an [`InelasticPoint`] built from atomic units. One definition of
//! each check and each conversion lives here, so the three cannot drift apart
//! in their validity rules or in the wording of the errors they return (#745).
//!
//! Nothing here is public: the models keep their own public methods and doc
//! comments, and the errors are the same `ElectronDataError::Invalid` values
//! they have always returned.

use super::penn::InelasticPoint;
use crate::constants::BOHR_RADIUS;
use crate::electron::data::ElectronDataError;
use crate::electron::elastic::HARTREE_EV;

type Result<T> = std::result::Result<T, ElectronDataError>;

/// The loosest relative integration tolerance any of the models accepts.
/// The tightest one is per model: the single pole evaluates a closed-form
/// kernel and goes to [`super::penn::DEFAULT_RELATIVE_TOLERANCE`] and below,
/// whereas the full Penn and Mermin models nest quadratures and cannot
/// resolve better than `1e-9`.
pub(super) const MAX_RELATIVE_TOLERANCE: f64 = 1e-2;

/// `v`, an energy in eV described by `what`, must be finite and positive.
pub(super) fn check_energy(what: &'static str, v: f64) -> Result<()> {
    if v.is_finite() && v > 0.0 {
        Ok(())
    } else {
        Err(ElectronDataError::Invalid {
            what,
            reason: format!("must be finite and positive, got {v} eV"),
        })
    }
}

/// `fermi_ev` must be finite and non-negative; returns it in Hartree, the
/// unit the models hold it in.
pub(super) fn fermi_energy_hartree(fermi_ev: f64) -> Result<f64> {
    if !(fermi_ev.is_finite() && fermi_ev >= 0.0) {
        return Err(ElectronDataError::Invalid {
            what: "Fermi energy",
            reason: format!("must be finite and non-negative, got {fermi_ev} eV"),
        });
    }
    Ok(fermi_ev / HARTREE_EV)
}

/// `rel_tol` must lie in `[min_tol, MAX_RELATIVE_TOLERANCE]`. The bounds are
/// written into the message in exponent form, as the models' doc comments
/// write them.
pub(super) fn check_relative_tolerance(rel_tol: f64, min_tol: f64) -> Result<()> {
    if !(min_tol..=MAX_RELATIVE_TOLERANCE).contains(&rel_tol) {
        return Err(ElectronDataError::Invalid {
            what: "integration tolerance",
            reason: format!("must be in [{min_tol:e}, {MAX_RELATIVE_TOLERANCE:e}], got {rel_tol}"),
        });
    }
    Ok(())
}

/// The point at kinetic energy `energy_ev` (eV) from an inverse IMFP
/// `inverse_imfp_au` (per bohr) and a stopping power `stopping_au` (hartree
/// per bohr).
pub(super) fn inelastic_point_from_au(
    energy_ev: f64,
    inverse_imfp_au: f64,
    stopping_au: f64,
) -> InelasticPoint {
    InelasticPoint {
        energy_ev,
        inverse_imfp_per_m: inverse_imfp_au / BOHR_RADIUS,
        stopping_ev_per_m: stopping_au * HARTREE_EV / BOHR_RADIUS,
    }
}

#[cfg(test)]
mod tests {
    use crate::electron::data::{ElectronDataError, OpticalElf};
    use crate::electron::inelastic::{
        DrudeLorentz, DrudeLorentzOscillator, FullPenn, MerminFitOptions, MerminPenn,
        SinglePolePenn,
    };

    fn elf() -> OpticalElf {
        DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])
            .unwrap()
            .to_optical_elf("synthetic Drude plasmon", 1.0, 1e3, 40)
            .unwrap()
    }

    fn reason(e: ElectronDataError) -> (&'static str, String) {
        match e {
            ElectronDataError::Invalid { what, reason } => (what, reason),
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    fn mermin() -> MerminPenn {
        let options = MerminFitOptions {
            n_oscillators: 1,
            max_iterations: 5,
            ..Default::default()
        };
        MerminPenn::fit(elf(), &options).unwrap()
    }

    /// The three models reject a non-positive electron energy with the same
    /// `what` and the same wording (#745).
    #[test]
    fn the_three_models_reject_a_bad_energy_alike() {
        let single = SinglePolePenn::new(elf());
        let full = FullPenn::new(elf());
        let mermin = mermin();
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let want = (
                "electron energy",
                format!("must be finite and positive, got {bad} eV"),
            );
            for (name, got) in [
                ("single-pole", single.imfp_and_stopping(bad).unwrap_err()),
                ("full", full.imfp_and_stopping(bad).unwrap_err()),
                ("mermin", mermin.imfp_and_stopping(bad).unwrap_err()),
            ] {
                assert_eq!(reason(got), want, "{name} at {bad}");
            }
        }
    }

    /// Same for a negative or non-finite Fermi energy (#745).
    #[test]
    fn the_three_models_reject_a_bad_fermi_energy_alike() {
        for bad in [-1.0, f64::NAN, f64::NEG_INFINITY] {
            let want = (
                "Fermi energy",
                format!("must be finite and non-negative, got {bad} eV"),
            );
            let single = SinglePolePenn::new(elf()).with_fermi_energy_ev(bad);
            let full = FullPenn::new(elf()).with_fermi_energy_ev(bad);
            let mermin = mermin().with_fermi_energy_ev(bad);
            for (name, got) in [
                ("single-pole", single.err().unwrap()),
                ("full", full.err().unwrap()),
                ("mermin", mermin.err().unwrap()),
            ] {
                assert_eq!(reason(got), want, "{name} at {bad}");
            }
        }
    }

    /// The tolerance range differs by its lower bound only, and the bounds
    /// are spelled in the message exactly as the doc comments spell them
    /// (#745): `1e-12` for the single pole, `1e-9` for the other two.
    #[test]
    fn the_tolerance_message_names_each_models_lower_bound() {
        let bad = 1.0;
        for (name, low, got) in [
            (
                "single-pole",
                "1e-12",
                SinglePolePenn::new(elf())
                    .with_relative_tolerance(bad)
                    .err()
                    .unwrap(),
            ),
            (
                "full",
                "1e-9",
                FullPenn::new(elf())
                    .with_relative_tolerance(bad)
                    .err()
                    .unwrap(),
            ),
            (
                "mermin",
                "1e-9",
                mermin().with_relative_tolerance(bad).err().unwrap(),
            ),
        ] {
            assert_eq!(
                reason(got),
                (
                    "integration tolerance",
                    format!("must be in [{low}, 1e-2], got {bad}"),
                ),
                "{name}"
            );
        }
        // The lower bound is enforced, not merely quoted.
        assert!(SinglePolePenn::new(elf())
            .with_relative_tolerance(1e-12)
            .is_ok());
        assert!(FullPenn::new(elf()).with_relative_tolerance(1e-12).is_err());
        assert!(mermin().with_relative_tolerance(1e-9).is_ok());
    }
}
