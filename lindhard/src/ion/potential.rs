//! Screened-Coulomb interatomic potentials.
//!
//! The interaction of two atoms with atomic numbers `Z1`, `Z2` at separation
//! `r` is written
//!
//! ```text
//! V(r) = Z1 Z2 e^2 / (4 pi eps0 r) * phi(r / a)
//! ```
//!
//! where `phi` is a universal screening function (`phi(0) = 1`, decreasing to 0)
//! and `a` a screening length that depends on `Z1` and `Z2`. Everything in the
//! scattering code is expressed in the reduced variables `x = r / a` and
//! `eps = a E_cm / (Z1 Z2 e^2 / (4 pi eps0))`, in which `phi` is the only
//! input; the screening length only converts to and from SI.
//!
//! Coefficients follow the cited papers. The primary sources have not been
//! re-read digit by digit for this module. The ZBL, Kr-C and Moliere sets have
//! been cross-checked against an independent MIT-licensed implementation, and
//! the Lenz-Jensen set has not been cross-checked at all; see the
//! `docs/data-provenance.md` rows for the status of each set.

use crate::constants::{BOHR_RADIUS, COULOMB_E2};

/// Screening function `phi(x)`, `x = r / a`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Screening {
    /// ZBL "universal" screening: four exponentials fitted to Hartree-Fock-Slater
    /// solid-state pair potentials. J. F. Ziegler, J. P. Biersack, U. Littmark,
    /// *The Stopping and Range of Ions in Solids* (Pergamon, 1985), Ch. 2.
    ZblUniversal,
    /// Kr-C potential: three exponentials fitted to Hartree-Fock Kr-Kr pair
    /// interaction. W. D. Wilson, L. G. Haggmark, J. P. Biersack,
    /// Phys. Rev. B 15, 2458 (1977).
    KrC,
    /// Moliere's approximation to the Thomas-Fermi function: three exponentials.
    /// G. Moliere, Z. Naturforsch. A 2, 133 (1947).
    Moliere,
    /// Lenz-Jensen screening, a polynomial-times-exponential approximation of
    /// the Thomas-Fermi-Jensen statistical model. W. Lenz, Z. Phys. 77, 713
    /// (1932); H. Jensen, Z. Phys. 77, 722 (1932). Polynomial form as tabulated
    /// in e.g. Nastasi, Mayer, Hirvonen, *Ion-Solid Interactions* (CUP, 1996).
    LenzJensen,
}

/// Prefactors and exponents `(c_i, b_i)` of the sum-of-exponentials screening
/// functions, `phi(x) = Sum_i c_i exp(-b_i x)`.
///
/// Cross-checked digit for digit (Kr-C, Moliere) or to the printed rounding
/// (ZBL, whose 4-digit values here are the commonly printed rounding of
/// 0.18175, 0.50986, 0.28022, 0.028171 / 3.19980, 0.94229, 0.40290, 0.20162)
/// against `include/screened_coulomb.h` of the MIT-licensed `ir2-lab/screened_coulomb`
/// (commit f84c3c8), a component of OpenTRIM. This is a secondary source, not
/// the cited papers.
const ZBL_C: [(f64, f64); 4] = [
    (0.181_8, 3.2),
    (0.509_9, 0.942_3),
    (0.280_2, 0.402_9),
    (0.028_17, 0.201_6),
];
const KRC_C: [(f64, f64); 3] = [
    (0.190_945, 0.278_544),
    (0.473_674, 0.637_174),
    (0.335_381, 1.919_249),
];
const MOLIERE_C: [(f64, f64); 3] = [(0.35, 0.3), (0.55, 1.2), (0.10, 6.0)];

/// Lenz-Jensen: `phi = P(q) exp(-q)`, `q = sqrt(LJ_Q2 x)`,
/// `P(q) = 1 + q + 0.3344 q^2 + 0.0485 q^3 + 0.002647 q^4`.
///
/// Not verified against any source: no permissively licensed implementation
/// carries this set, and the papers were not available when it was entered.
/// The argument scale is the least certain value (`9.67` here; a scale of
/// `3.11126^2 = 9.680` is also plausible). See `docs/data-provenance.md`.
const LJ_P: [f64; 5] = [1.0, 1.0, 0.3344, 0.0485, 0.002_647];
const LJ_Q2: f64 = 9.67;

fn sum_exp(c: &[(f64, f64)], x: f64) -> f64 {
    c.iter().map(|&(a, b)| a * (-b * x).exp()).sum()
}
fn sum_exp_d(c: &[(f64, f64)], x: f64) -> f64 {
    c.iter().map(|&(a, b)| -a * b * (-b * x).exp()).sum()
}

impl Screening {
    /// All implemented screening functions.
    pub const ALL: [Screening; 4] = [
        Screening::ZblUniversal,
        Screening::KrC,
        Screening::Moliere,
        Screening::LenzJensen,
    ];

    /// Screening function value `phi(x)`.
    pub fn phi(self, x: f64) -> f64 {
        match self {
            Screening::ZblUniversal => sum_exp(&ZBL_C, x),
            Screening::KrC => sum_exp(&KRC_C, x),
            Screening::Moliere => sum_exp(&MOLIERE_C, x),
            Screening::LenzJensen => {
                let q = (LJ_Q2 * x).sqrt();
                horner(&LJ_P, q) * (-q).exp()
            }
        }
    }

    /// Derivative `dphi/dx`.
    pub fn dphi(self, x: f64) -> f64 {
        match self {
            Screening::ZblUniversal => sum_exp_d(&ZBL_C, x),
            Screening::KrC => sum_exp_d(&KRC_C, x),
            Screening::Moliere => sum_exp_d(&MOLIERE_C, x),
            Screening::LenzJensen => {
                if x <= 0.0 {
                    return f64::NEG_INFINITY;
                }
                let q = (LJ_Q2 * x).sqrt();
                // d/dq [P e^-q] = (P' - P) e^-q ;  dq/dx = q / (2x).
                let p = horner(&LJ_P, q);
                let dp = LJ_P
                    .iter()
                    .enumerate()
                    .skip(1)
                    .rev()
                    .fold(0.0, |acc, (i, &c)| acc * q + c * i as f64);
                (dp - p) * (-q).exp() * q / (2.0 * x)
            }
        }
    }

    /// The `(c_i, b_i)` terms of a sum-of-exponentials screening function
    /// `phi(x) = Sum_i c_i exp(-b_i x)`, or `None` for Lenz-Jensen. Used by the
    /// analytic stopping asymptote in the scattering tests.
    #[cfg(test)]
    pub(crate) fn exponential_terms(self) -> Option<&'static [(f64, f64)]> {
        match self {
            Screening::ZblUniversal => Some(&ZBL_C),
            Screening::KrC => Some(&KRC_C),
            Screening::Moliere => Some(&MOLIERE_C),
            Screening::LenzJensen => None,
        }
    }

    /// The screening length conventionally paired with this function (the
    /// pairing used where the function was introduced or fitted):
    /// ZBL universal with the universal length, Kr-C and Moliere with the
    /// Firsov length, Lenz-Jensen with the Lindhard (Thomas-Fermi) length.
    pub fn default_length(self) -> ScreeningLength {
        match self {
            Screening::ZblUniversal => ScreeningLength::Universal,
            Screening::KrC | Screening::Moliere => ScreeningLength::Firsov,
            Screening::LenzJensen => ScreeningLength::Lindhard,
        }
    }
}

fn horner(c: &[f64], q: f64) -> f64 {
    c.iter().rev().fold(0.0, |acc, &k| acc * q + k)
}

/// Screening length `a(Z1, Z2)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ScreeningLength {
    /// Universal length, `a_U = 0.8854 a0 / (Z1^0.23 + Z2^0.23)`.
    /// Ziegler, Biersack, Littmark (1985), Ch. 2.
    Universal,
    /// Firsov length, `a_F = 0.8853 a0 / (sqrt(Z1) + sqrt(Z2))^(2/3)`.
    /// O. B. Firsov, Sov. Phys. JETP 6, 534 (1958).
    Firsov,
    /// Lindhard (Thomas-Fermi) length,
    /// `a_L = 0.8853 a0 / (Z1^(2/3) + Z2^(2/3))^(1/2)`.
    /// J. Lindhard, M. Scharff, H. E. Schiott, Mat. Fys. Medd. Dan. Vid. Selsk.
    /// 33(14) (1963).
    Lindhard,
}

/// The Thomas-Fermi constant `(9 pi^2 / 128)^(1/3) = 0.8853`, which multiplies
/// the Bohr radius in the Firsov and Lindhard lengths.
pub fn thomas_fermi_constant() -> f64 {
    (9.0 * std::f64::consts::PI * std::f64::consts::PI / 128.0).cbrt()
}

impl ScreeningLength {
    /// Screening length in metres for atomic numbers `z1`, `z2`.
    pub fn metres(self, z1: f64, z2: f64) -> f64 {
        let a0 = BOHR_RADIUS;
        match self {
            ScreeningLength::Universal => 0.8854 * a0 / (z1.powf(0.23) + z2.powf(0.23)),
            ScreeningLength::Firsov => {
                thomas_fermi_constant() * a0 / (z1.sqrt() + z2.sqrt()).powf(2.0 / 3.0)
            }
            ScreeningLength::Lindhard => {
                thomas_fermi_constant() * a0 / (z1.powf(2.0 / 3.0) + z2.powf(2.0 / 3.0)).sqrt()
            }
        }
    }
}

/// A screened-Coulomb pair potential for a specific `(Z1, Z2)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Potential {
    /// Screening function.
    pub screening: Screening,
    /// Screening length convention.
    pub length: ScreeningLength,
    /// Projectile atomic number.
    pub z1: f64,
    /// Target atomic number.
    pub z2: f64,
}

impl Potential {
    /// Potential with the conventional length for the chosen screening function.
    pub fn new(screening: Screening, z1: f64, z2: f64) -> Self {
        Self {
            screening,
            length: screening.default_length(),
            z1,
            z2,
        }
    }

    /// Override the screening length.
    pub fn with_length(mut self, length: ScreeningLength) -> Self {
        self.length = length;
        self
    }

    /// Screening length `a`, metres.
    pub fn screening_length(&self) -> f64 {
        self.length.metres(self.z1, self.z2)
    }

    /// Potential energy at separation `r` (metres), joules.
    pub fn energy(&self, r: f64) -> f64 {
        self.z1 * self.z2 * COULOMB_E2 / r * self.screening.phi(r / self.screening_length())
    }

    /// Reduced centre-of-mass energy `eps = a E_cm / (Z1 Z2 e^2/(4 pi eps0))`,
    /// from `E_cm` in joules. (Lindhard-Scharff-Schiott reduced energy.)
    pub fn reduced_energy(&self, e_cm_j: f64) -> f64 {
        self.screening_length() * e_cm_j / (self.z1 * self.z2 * COULOMB_E2)
    }

    /// Centre-of-mass energy in joules from lab energy `e_lab_j` for projectile
    /// mass `m1` and target mass `m2` (any common unit): `E_cm = E M2/(M1+M2)`.
    pub fn cm_energy(e_lab_j: f64, m1: f64, m2: f64) -> f64 {
        e_lab_j * m2 / (m1 + m2)
    }

    /// Impact parameter in metres from reduced `beta = b / a`.
    pub fn impact_parameter(&self, beta: f64) -> f64 {
        beta * self.screening_length()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phi_is_one_at_origin_and_decreasing() {
        for s in Screening::ALL {
            let p0 = s.phi(1e-12);
            assert!((p0 - 1.0).abs() < 5e-3, "{s:?} phi(0) = {p0}");
            let mut prev = p0;
            for k in 1..200 {
                let v = s.phi(k as f64 * 0.1);
                assert!(v < prev && v > 0.0, "{s:?} not decreasing at {k}");
                prev = v;
            }
        }
    }

    #[test]
    fn derivative_matches_finite_difference() {
        for s in Screening::ALL {
            for &x in &[0.05, 0.3, 1.0, 3.0, 10.0] {
                let h = 1e-6 * x;
                let fd = (s.phi(x + h) - s.phi(x - h)) / (2.0 * h);
                let an = s.dphi(x);
                assert!((fd - an).abs() < 1e-6 * an.abs().max(1e-3), "{s:?} x={x}");
            }
        }
    }

    #[test]
    fn functions_agree_with_thomas_fermi_scale() {
        // All four approximate the Thomas-Fermi function, phi_TF(1) ~ 0.40.
        for s in Screening::ALL {
            let v = s.phi(1.0);
            assert!((0.3..0.5).contains(&v), "{s:?} phi(1) = {v}");
        }
    }

    #[test]
    fn screening_lengths_match_closed_forms() {
        // Explicit closed forms for Z1 = Z2 = 14.
        let a0 = BOHR_RADIUS;
        let tf = thomas_fermi_constant();
        assert!((tf - 0.8853).abs() < 1e-4);
        let a = ScreeningLength::Lindhard.metres(14.0, 14.0);
        let expect = tf * a0 / (2.0 * 14f64.powf(2.0 / 3.0)).sqrt();
        assert!((a / expect - 1.0).abs() < 1e-12);
        let a = ScreeningLength::Firsov.metres(14.0, 14.0);
        let expect = tf * a0 / (2.0 * 14f64.sqrt()).powf(2.0 / 3.0);
        assert!((a / expect - 1.0).abs() < 1e-12);
    }

    #[test]
    fn exponential_sets_are_normalised() {
        // phi(0) = 1 requires the prefactors to sum to 1.
        let sum = |c: &[(f64, f64)]| c.iter().map(|&(a, _)| a).sum::<f64>();
        assert!((sum(&ZBL_C) - 1.0).abs() < 1e-4);
        assert!((sum(&KRC_C) - 1.0).abs() < 1e-5);
        assert!((sum(&MOLIERE_C) - 1.0).abs() < 1e-12);
        // The ZBL exponents are 3.2, 0.9423, 0.4029, 0.2016.
        assert_eq!(ZBL_C[2].1, 0.4029);
    }

    #[test]
    fn reduced_energy_round_trip() {
        let p = Potential::new(Screening::ZblUniversal, 5.0, 14.0);
        let e = 1.0e3 * crate::units::J_PER_EV;
        let eps = p.reduced_energy(e);
        assert!(eps > 0.0 && eps < 1.0);
        // V(r) at r = a: Z1 Z2 e^2/a * phi(1), in units of E_cm = eps^-1 ...
        let a = p.screening_length();
        let v = p.energy(a);
        let expect = p.z1 * p.z2 * COULOMB_E2 / a * p.screening.phi(1.0);
        assert!((v / expect - 1.0).abs() < 1e-12);
    }
}
