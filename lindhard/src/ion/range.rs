//! Nuclear stopping in physical units and the continuous-slowing-down (CSDA)
//! range, for stopping tables.
//!
//! # Conventions
//!
//! * Stopping cross sections are per target atom in **J m²**, as in
//!   [`crate::ion::stopping`]; ranges are in metres.
//! * The nuclear stopping cross section is the reduced `s_n(eps)` of
//!   [`nuclear_stopping_reduced`] in physical units (Lindhard, Scharff and
//!   Schiott 1963; Ziegler, Biersack and Littmark 1985, ch. 2):
//!
//!   ```text
//!   S_n = 4 pi a Z1 Z2 e^2/(4 pi eps0) * M1/(M1 + M2) * s_n(eps)
//!   ```
//!
//!   which follows from `S_n = Int T dsigma` with `T = gamma E sin^2(theta/2)`
//!   and `eps = a E M2 / (Z1 Z2 e^2 (M1 + M2))`.
//! * The CSDA range is `R(E) = Int_{E0}^{E} dE' / (N S(E'))` with `S` the
//!   total (electronic plus nuclear) stopping cross section per atom and `N`
//!   the atom number density. It is the mean path length of an ion that loses
//!   energy continuously; it ignores fluctuations of the energy loss and is
//!   not the projected range.
//!
//! Every function is a pure function of its arguments, so results are
//! bitwise reproducible across threads.

use super::potential::Potential;
use super::scattering::nuclear_stopping_reduced;
use crate::constants::COULOMB_E2;
use crate::units::J_PER_EV;
use std::f64::consts::PI;

/// Largest reduced impact parameter of the nuclear stopping integral.
const BETA_MAX: f64 = 1.0e3;
/// Simpson intervals per decade of impact parameter.
const BETA_PER_DECADE: usize = 40;

/// Nuclear stopping cross section per target atom, J m², of an ion of mass
/// `m1_amu` and the target atom of mass `m2_amu` in the pair potential
/// `potential`, at laboratory energy `energy_ev`.
///
/// Returns `None` when an input is not finite and positive.
pub fn nuclear_stopping_cross_section(
    potential: &Potential,
    m1_amu: f64,
    m2_amu: f64,
    energy_ev: f64,
) -> Option<f64> {
    let ok = |x: f64| x.is_finite() && x > 0.0;
    if !(ok(m1_amu) && ok(m2_amu) && ok(energy_ev)) {
        return None;
    }
    let e_cm = Potential::cm_energy(energy_ev * J_PER_EV, m1_amu, m2_amu);
    let eps = potential.reduced_energy(e_cm);
    // The head-on part neglected below beta_min is ~ beta_min^2, and the
    // Rutherford-like regime reaches out to ~ 1/eps (convergence checked in
    // `ion::scattering` tests).
    let beta_min = 1e-4f64.min(1e-3 / eps);
    let sn = nuclear_stopping_reduced(
        potential.screening,
        eps,
        beta_min,
        BETA_MAX,
        BETA_PER_DECADE,
    );
    Some(
        4.0 * PI * potential.screening_length() * potential.z1 * potential.z2 * COULOMB_E2 * m1_amu
            / (m1_amu + m2_amu)
            * sn,
    )
}

/// Cells per decade of energy of the CSDA integral.
const CELLS_PER_DECADE: f64 = 8.0;

/// Nodes and weights of 5-point Gauss-Legendre quadrature on `[-1, 1]`.
const GL5: [(f64, f64); 5] = [
    (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
    (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
    (0.0, 0.568_888_888_888_888_9),
    (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
    (0.906_179_845_938_664, 0.236_926_885_056_189_1),
];

/// CSDA range, metres, at each of `energies_ev` (ascending, each at least
/// `from_ev`), integrated from `from_ev`.
///
/// Energies are in eV. `total_stopping(E)` is the total stopping cross section
/// per atom at `E` (J m², positive) and `number_density` the atom density (m⁻³). The integral
/// runs over `ln E` in cells of fixed size measured from `from_ev`, by
/// 5-point Gauss-Legendre in each, so the range at an energy depends on
/// `from_ev` and the stopping, never on the other energies asked for. The
/// error of `total_stopping` is returned as is.
///
/// # Panics
///
/// If `energies_ev` is not ascending or has a value below `from_ev`, or
/// `from_ev` or `number_density` is not finite and positive.
pub fn csda_ranges<E>(
    mut total_stopping: impl FnMut(f64) -> Result<f64, E>,
    number_density: f64,
    from_ev: f64,
    energies_ev: &[f64],
) -> Result<Vec<f64>, E> {
    assert!(from_ev.is_finite() && from_ev > 0.0, "from_ev");
    assert!(
        number_density.is_finite() && number_density > 0.0,
        "number_density"
    );
    assert!(
        energies_ev.windows(2).all(|w| w[0] <= w[1]) && energies_ev.iter().all(|&e| e >= from_ev),
        "energies must be ascending and not below from_ev"
    );
    let du = std::f64::consts::LN_10 / CELLS_PER_DECADE;
    // Integral over u = ln E from u_a to u_b of E / (N S(E)) du.
    let mut cell = |ua: f64, ub: f64| -> Result<f64, E> {
        let (mid, half) = (0.5 * (ua + ub), 0.5 * (ub - ua));
        let mut acc = 0.0;
        for (x, w) in GL5 {
            let e = (mid + half * x).exp();
            acc += w * e * J_PER_EV / (number_density * total_stopping(e)?);
        }
        Ok(acc * half)
    };
    let u0 = from_ev.ln();
    // cumulative[k]: range from `from_ev` to the node k cells above it.
    let mut cumulative = vec![0.0];
    let mut out = Vec::with_capacity(energies_ev.len());
    for &e in energies_ev {
        let u = e.ln();
        let k = (((u - u0) / du).floor().max(0.0)) as usize;
        while cumulative.len() <= k {
            let j = cumulative.len() - 1;
            let step = cell(u0 + j as f64 * du, u0 + (j + 1) as f64 * du)?;
            cumulative.push(cumulative[j] + step);
        }
        let node = u0 + k as f64 * du;
        let rest = if u > node { cell(node, u)? } else { 0.0 };
        out.push(cumulative[k] + rest);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ion::potential::Screening;

    /// For `S = c sqrt(E)`, `R = 2 (sqrt(E) - sqrt(E0)) / (N c)` (E in joules).
    #[test]
    fn csda_range_matches_closed_form_for_sqrt_energy_stopping() {
        let (c, n, e0) = (3.0e-30, 5.0e28, 10.0);
        let energies = [10.0, 11.5, 100.0, 1.0e3, 7.7e4, 1.0e6];
        let r = csda_ranges::<()>(|e| Ok(c * e.sqrt()), n, e0, &energies).unwrap();
        for (&e, &got) in energies.iter().zip(&r) {
            let want = 2.0 * (e.sqrt() - e0.sqrt()) * J_PER_EV / (n * c);
            assert!((got / want - 1.0).abs() < 1e-10 || (got == 0.0 && want == 0.0));
        }
    }

    #[test]
    fn csda_range_does_not_depend_on_the_other_energies() {
        let s = |e: f64| Ok::<_, ()>(2.0e-30 * e.powf(0.4) + 1e-31);
        let all = csda_ranges(s, 6e28, 1.0, &[1.0, 3.0, 50.0, 4000.0]).unwrap();
        let one = csda_ranges(s, 6e28, 1.0, &[4000.0]).unwrap();
        assert_eq!(all[3].to_bits(), one[0].to_bits());
    }

    #[test]
    fn csda_range_returns_the_stopping_error() {
        let r = csda_ranges(
            |e| if e > 50.0 { Err("out") } else { Ok(1e-30) },
            1e28,
            1.0,
            &[100.0],
        );
        assert_eq!(r, Err("out"));
    }

    /// Reduced and physical forms agree: `S_n = 4 pi a Z1 Z2 e^2 M1/(M1+M2) s_n`.
    #[test]
    fn nuclear_stopping_is_the_reduced_value_in_physical_units() {
        let pot = Potential::new(Screening::ZblUniversal, 5.0, 14.0);
        let (m1, m2, e) = (11.009, 28.086, 2.0e3);
        let eps = pot.reduced_energy(Potential::cm_energy(e * J_PER_EV, m1, m2));
        let sn = nuclear_stopping_reduced(
            Screening::ZblUniversal,
            eps,
            1e-4f64.min(1e-3 / eps),
            BETA_MAX,
            BETA_PER_DECADE,
        );
        let want =
            4.0 * PI * pot.screening_length() * 5.0 * 14.0 * COULOMB_E2 * m1 / (m1 + m2) * sn;
        let got = nuclear_stopping_cross_section(&pot, m1, m2, e).unwrap();
        assert!((got / want - 1.0).abs() < 1e-14);
        // Order of magnitude: B in Si at 2 keV is a few 1e-15 eV cm^2 scale.
        let ev_1e15 = crate::ion::stopping::to_ev_1e15_cm2(got);
        assert!(ev_1e15 > 10.0 && ev_1e15 < 1.0e3, "{ev_1e15}");
    }

    #[test]
    fn nuclear_stopping_rejects_bad_input() {
        let pot = Potential::new(Screening::ZblUniversal, 5.0, 14.0);
        for (m1, m2, e) in [(0.0, 28.0, 1e3), (11.0, f64::NAN, 1e3), (11.0, 28.0, -1.0)] {
            assert!(nuclear_stopping_cross_section(&pot, m1, m2, e).is_none());
        }
    }
}
