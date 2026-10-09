//! Table I of Salvat, Martinez, Mayol & Parellada, "Analytical
//! Dirac-Hartree-Fock-Slater screening function for atoms (Z = 1-92)",
//! Phys. Rev. A 36, 467 (1987), doi:10.1103/PhysRevA.36.467, pp. 470-471:
//! the parameters `A_1`, `A_2`, `alpha_1`, `alpha_2`, `alpha_3` of the
//! screening function `phi(r) = sum_i A_i exp(-alpha_i r)` (Eq. (11)), `r` in
//! bohr, `alpha_i` in 1/bohr. `A_3` is not printed: `A_3 = 1 - A_1 - A_2`.
//!
//! Elements marked with an asterisk in Table I ("DHFS radial expected values
//! inconsistent with conditions (15)") were fitted with `A_3 = 0` (p. 471:
//! "parameters for those elements have been determined by setting
//! `A_3 = 0`"); their `alpha_3` column is blank and is stored here as `None`.
//! They are H through P (Z = 1..15), Ca, Sc, Se, Br, Kr and Xe.
//!
//! Source copy: the open copy of the paper in the University of Barcelona
//! repository, <https://diposit.ub.edu/items/fc6162de-89b1-445e-bcd9-11b5c5b3eaf5>
//! (PDF bitstream fbde82c8-6d92-4de0-8e96-bdf4ee698a47, 1 194 975 bytes, MD5
//! 29d2f37598cc70edd957f946866b3887; an image-only scan). Transcribed from
//! the page images rendered with poppler `pdftoppm` (journal pp. 470-471 =
//! PDF pages 4-5), then read a second time row by row from 400 dpi crops;
//! the two reads agree on every value (issue #130). Every value below is
//! copied as printed; none was entered from memory. See
//! `docs/data-provenance.md`.

/// One row of Table I.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Row {
    /// Atomic number.
    pub z: u32,
    /// `A_1` (dimensionless).
    pub a1: f64,
    /// `A_2` (dimensionless).
    pub a2: f64,
    /// `alpha_1`, 1/bohr.
    pub alpha1: f64,
    /// `alpha_2`, 1/bohr.
    pub alpha2: f64,
    /// `alpha_3`, 1/bohr; `None` for the asterisked two-term rows
    /// (`A_3 = 0`).
    pub alpha3: Option<f64>,
}

/// Table I, Z = 1..92 in order (`TABLE_I[z - 1]` is element `z`). Row
/// comments give the element symbol, with `*` for the asterisked rows.
// Rh `alpha_2 = 2.7183` is the printed fit value, not Euler's number.
#[allow(clippy::approx_constant)]
#[rustfmt::skip]
pub(super) const TABLE_I: [Row; 92] = [
    // H*
    Row { z: 1, a1: -184.39, a2: 185.39, alpha1: 2.0027, alpha2: 1.9973, alpha3: None },
    // He*
    Row { z: 2, a1: -0.2259, a2: 1.2259, alpha1: 5.5272, alpha2: 2.3992, alpha3: None },
    // Li*
    Row { z: 3, a1: 0.6045, a2: 0.3955, alpha1: 2.8174, alpha2: 0.6625, alpha3: None },
    // Be*
    Row { z: 4, a1: 0.3278, a2: 0.6722, alpha1: 4.5430, alpha2: 0.9852, alpha3: None },
    // B*
    Row { z: 5, a1: 0.2327, a2: 0.7673, alpha1: 5.9900, alpha2: 1.2135, alpha3: None },
    // C*
    Row { z: 6, a1: 0.1537, a2: 0.8463, alpha1: 8.0404, alpha2: 1.4913, alpha3: None },
    // N*
    Row { z: 7, a1: 0.0996, a2: 0.9004, alpha1: 10.812, alpha2: 1.7687, alpha3: None },
    // O*
    Row { z: 8, a1: 0.0625, a2: 0.9375, alpha1: 14.823, alpha2: 2.0403, alpha3: None },
    // F*
    Row { z: 9, a1: 0.0368, a2: 0.9632, alpha1: 21.400, alpha2: 2.3060, alpha3: None },
    // Ne*
    Row { z: 10, a1: 0.0188, a2: 0.9812, alpha1: 34.999, alpha2: 2.5662, alpha3: None },
    // Na*
    Row { z: 11, a1: 0.7444, a2: 0.2556, alpha1: 4.1205, alpha2: 0.8718, alpha3: None },
    // Mg*
    Row { z: 12, a1: 0.6423, a2: 0.3577, alpha1: 4.7266, alpha2: 1.0025, alpha3: None },
    // Al*
    Row { z: 13, a1: 0.6002, a2: 0.3998, alpha1: 5.1405, alpha2: 1.0153, alpha3: None },
    // Si*
    Row { z: 14, a1: 0.5160, a2: 0.4840, alpha1: 5.8492, alpha2: 1.1732, alpha3: None },
    // P*
    Row { z: 15, a1: 0.4387, a2: 0.5613, alpha1: 6.6707, alpha2: 1.3410, alpha3: None },
    // S
    Row { z: 16, a1: 0.5459, a2: -0.5333, alpha1: 6.3703, alpha2: 2.5517, alpha3: Some(1.6753) },
    // Cl
    Row { z: 17, a1: 0.7249, a2: -0.7548, alpha1: 6.2118, alpha2: 3.3883, alpha3: Some(1.8596) },
    // Ar
    Row { z: 18, a1: 2.1912, a2: -2.2852, alpha1: 5.5470, alpha2: 4.5687, alpha3: Some(2.0446) },
    // K
    Row { z: 19, a1: 0.0486, a2: 0.7759, alpha1: 30.260, alpha2: 3.1243, alpha3: Some(0.7326) },
    // Ca*
    Row { z: 20, a1: 0.5800, a2: 0.4200, alpha1: 6.3218, alpha2: 1.0094, alpha3: None },
    // Sc*
    Row { z: 21, a1: 0.5543, a2: 0.4457, alpha1: 6.6328, alpha2: 1.1023, alpha3: None },
    // Ti
    Row { z: 22, a1: 0.0112, a2: 0.6832, alpha1: 99.757, alpha2: 4.1286, alpha3: Some(1.0090) },
    // V
    Row { z: 23, a1: 0.0318, a2: 0.6753, alpha1: 42.533, alpha2: 3.9404, alpha3: Some(1.0533) },
    // Cr
    Row { z: 24, a1: 0.1075, a2: 0.7162, alpha1: 18.959, alpha2: 3.0638, alpha3: Some(1.0014) },
    // Mn
    Row { z: 25, a1: 0.0498, a2: 0.6866, alpha1: 31.864, alpha2: 3.7811, alpha3: Some(1.1279) },
    // Fe
    Row { z: 26, a1: 0.0512, a2: 0.6995, alpha1: 31.825, alpha2: 3.7716, alpha3: Some(1.1606) },
    // Co
    Row { z: 27, a1: 0.0500, a2: 0.7142, alpha1: 32.915, alpha2: 3.7908, alpha3: Some(1.1915) },
    // Ni
    Row { z: 28, a1: 0.0474, a2: 0.7294, alpha1: 34.758, alpha2: 3.8299, alpha3: Some(1.2209) },
    // Cu
    Row { z: 29, a1: 0.0771, a2: 0.7951, alpha1: 25.326, alpha2: 3.3928, alpha3: Some(1.1426) },
    // Zn
    Row { z: 30, a1: 0.0400, a2: 0.7590, alpha1: 40.343, alpha2: 3.9465, alpha3: Some(1.2759) },
    // Ga
    Row { z: 31, a1: 0.1083, a2: 0.7489, alpha1: 20.192, alpha2: 3.4733, alpha3: Some(1.0064) },
    // Ge
    Row { z: 32, a1: 0.0610, a2: 0.7157, alpha1: 29.200, alpha2: 4.1252, alpha3: Some(1.1845) },
    // As
    Row { z: 33, a1: 0.0212, a2: 0.6709, alpha1: 62.487, alpha2: 4.9502, alpha3: Some(1.3582) },
    // Se*
    Row { z: 34, a1: 0.4836, a2: 0.5164, alpha1: 8.7824, alpha2: 1.6967, alpha3: None },
    // Br*
    Row { z: 35, a1: 0.4504, a2: 0.5496, alpha1: 9.3348, alpha2: 1.7900, alpha3: None },
    // Kr*
    Row { z: 36, a1: 0.4190, a2: 0.5810, alpha1: 9.9142, alpha2: 1.8835, alpha3: None },
    // Rb
    Row { z: 37, a1: 0.1734, a2: 0.7253, alpha1: 17.166, alpha2: 3.1103, alpha3: Some(0.7177) },
    // Sr
    Row { z: 38, a1: 0.0336, a2: 0.7816, alpha1: 55.208, alpha2: 4.2842, alpha3: Some(0.8578) },
    // Y
    Row { z: 39, a1: 0.0689, a2: 0.7202, alpha1: 31.366, alpha2: 4.2412, alpha3: Some(0.9472) },
    // Zr
    Row { z: 40, a1: 0.1176, a2: 0.6581, alpha1: 22.054, alpha2: 4.0325, alpha3: Some(1.0181) },
    // Nb
    Row { z: 41, a1: 0.2257, a2: 0.5821, alpha1: 14.240, alpha2: 2.9702, alpha3: Some(1.0170) },
    // Mo
    Row { z: 42, a1: 0.2693, a2: 0.5763, alpha1: 14.044, alpha2: 2.8611, alpha3: Some(1.0591) },
    // Tc
    Row { z: 43, a1: 0.2201, a2: 0.5618, alpha1: 15.918, alpha2: 3.3672, alpha3: Some(1.1548) },
    // Ru
    Row { z: 44, a1: 0.2751, a2: 0.5943, alpha1: 14.314, alpha2: 2.7370, alpha3: Some(1.1092) },
    // Rh
    Row { z: 45, a1: 0.2711, a2: 0.6119, alpha1: 14.654, alpha2: 2.7183, alpha3: Some(1.1234) },
    // Pd
    Row { z: 46, a1: 0.2784, a2: 0.6067, alpha1: 14.645, alpha2: 2.6155, alpha3: Some(1.4318) },
    // Ag
    Row { z: 47, a1: 0.2562, a2: 0.6505, alpha1: 15.588, alpha2: 2.7412, alpha3: Some(1.1408) },
    // Cd
    Row { z: 48, a1: 0.2271, a2: 0.6155, alpha1: 16.914, alpha2: 3.0841, alpha3: Some(1.2619) },
    // In
    Row { z: 49, a1: 0.2492, a2: 0.6440, alpha1: 16.155, alpha2: 2.8819, alpha3: Some(0.9942) },
    // Sn
    Row { z: 50, a1: 0.2153, a2: 0.6115, alpha1: 17.793, alpha2: 3.2937, alpha3: Some(1.1478) },
    // Sb
    Row { z: 51, a1: 0.1806, a2: 0.5767, alpha1: 19.875, alpha2: 3.8092, alpha3: Some(1.2829) },
    // Te
    Row { z: 52, a1: 0.1308, a2: 0.5504, alpha1: 24.154, alpha2: 4.6119, alpha3: Some(1.4195) },
    // I
    Row { z: 53, a1: 0.0588, a2: 0.5482, alpha1: 39.996, alpha2: 5.9132, alpha3: Some(1.5471) },
    // Xe*
    Row { z: 54, a1: 0.4451, a2: 0.5549, alpha1: 11.805, alpha2: 1.7967, alpha3: None },
    // Cs
    Row { z: 55, a1: 0.2708, a2: 0.6524, alpha1: 16.591, alpha2: 2.6964, alpha3: Some(0.6814) },
    // Ba
    Row { z: 56, a1: 0.1728, a2: 0.6845, alpha1: 22.397, alpha2: 3.4595, alpha3: Some(0.8073) },
    // La
    Row { z: 57, a1: 0.1947, a2: 0.6384, alpha1: 20.764, alpha2: 3.4657, alpha3: Some(0.8911) },
    // Ce
    Row { z: 58, a1: 0.1913, a2: 0.6467, alpha1: 21.235, alpha2: 3.4819, alpha3: Some(0.9011) },
    // Pr
    Row { z: 59, a1: 0.1868, a2: 0.6558, alpha1: 21.803, alpha2: 3.5098, alpha3: Some(0.9106) },
    // Nd
    Row { z: 60, a1: 0.1665, a2: 0.7057, alpha1: 23.949, alpha2: 3.5199, alpha3: Some(0.8486) },
    // Pm
    Row { z: 61, a1: 0.1624, a2: 0.7133, alpha1: 24.598, alpha2: 3.5560, alpha3: Some(0.8569) },
    // Sm
    Row { z: 62, a1: 0.1580, a2: 0.7210, alpha1: 25.297, alpha2: 3.5963, alpha3: Some(0.8650) },
    // Eu
    Row { z: 63, a1: 0.1538, a2: 0.7284, alpha1: 26.017, alpha2: 3.6383, alpha3: Some(0.8731) },
    // Gd
    Row { z: 64, a1: 0.1587, a2: 0.7024, alpha1: 25.497, alpha2: 3.7364, alpha3: Some(0.9550) },
    // Tb
    Row { z: 65, a1: 0.1453, a2: 0.7426, alpha1: 27.547, alpha2: 3.7288, alpha3: Some(0.8890) },
    // Dy
    Row { z: 66, a1: 0.1413, a2: 0.7494, alpha1: 28.346, alpha2: 3.7763, alpha3: Some(0.8969) },
    // Ho
    Row { z: 67, a1: 0.1374, a2: 0.7558, alpha1: 29.160, alpha2: 3.8244, alpha3: Some(0.9048) },
    // Er
    Row { z: 68, a1: 0.1336, a2: 0.7619, alpha1: 29.990, alpha2: 3.8734, alpha3: Some(0.9128) },
    // Tm
    Row { z: 69, a1: 0.1299, a2: 0.7680, alpha1: 30.835, alpha2: 3.9233, alpha3: Some(0.9203) },
    // Yb
    Row { z: 70, a1: 0.1267, a2: 0.7734, alpha1: 31.681, alpha2: 3.9727, alpha3: Some(0.9288) },
    // Lu
    Row { z: 71, a1: 0.1288, a2: 0.7528, alpha1: 31.353, alpha2: 4.0904, alpha3: Some(1.0072) },
    // Hf
    Row { z: 72, a1: 0.1303, a2: 0.7324, alpha1: 31.217, alpha2: 4.2049, alpha3: Some(1.0946) },
    // Ta
    Row { z: 73, a1: 0.1384, a2: 0.7096, alpha1: 30.077, alpha2: 4.2492, alpha3: Some(1.1697) },
    // W
    Row { z: 74, a1: 0.1500, a2: 0.6871, alpha1: 28.630, alpha2: 4.2426, alpha3: Some(1.2340) },
    // Re
    Row { z: 75, a1: 0.1608, a2: 0.6659, alpha1: 27.568, alpha2: 4.2341, alpha3: Some(1.2970) },
    // Os
    Row { z: 76, a1: 0.1722, a2: 0.6468, alpha1: 26.586, alpha2: 4.1999, alpha3: Some(1.3535) },
    // Ir
    Row { z: 77, a1: 0.1834, a2: 0.6306, alpha1: 25.734, alpha2: 4.1462, alpha3: Some(1.4037) },
    // Pt
    Row { z: 78, a1: 0.2230, a2: 0.6176, alpha1: 22.994, alpha2: 3.7346, alpha3: Some(1.4428) },
    // Au
    Row { z: 79, a1: 0.2289, a2: 0.6114, alpha1: 22.864, alpha2: 3.6914, alpha3: Some(1.4886) },
    // Hg
    Row { z: 80, a1: 0.2098, a2: 0.6004, alpha1: 24.408, alpha2: 3.9643, alpha3: Some(1.5343) },
    // Tl
    Row { z: 81, a1: 0.2708, a2: 0.6428, alpha1: 20.941, alpha2: 3.2456, alpha3: Some(1.1121) },
    // Pb
    Row { z: 82, a1: 0.2380, a2: 0.6308, alpha1: 22.987, alpha2: 3.6217, alpha3: Some(1.2373) },
    // Bi
    Row { z: 83, a1: 0.2288, a2: 0.6220, alpha1: 23.792, alpha2: 3.7796, alpha3: Some(1.2534) },
    // Po
    Row { z: 84, a1: 0.1941, a2: 0.6105, alpha1: 26.695, alpha2: 4.2582, alpha3: Some(1.3577) },
    // At
    Row { z: 85, a1: 0.1500, a2: 0.6031, alpha1: 31.840, alpha2: 4.9285, alpha3: Some(1.4683) },
    // Rn
    Row { z: 86, a1: 0.0955, a2: 0.6060, alpha1: 43.489, alpha2: 5.8520, alpha3: Some(1.5736) },
    // Fr
    Row { z: 87, a1: 0.3192, a2: 0.6233, alpha1: 20.015, alpha2: 2.9091, alpha3: Some(0.7207) },
    // Ra
    Row { z: 88, a1: 0.2404, a2: 0.6567, alpha1: 24.501, alpha2: 3.5524, alpha3: Some(0.8376) },
    // Ac
    Row { z: 89, a1: 0.2266, a2: 0.6422, alpha1: 25.684, alpha2: 3.7922, alpha3: Some(0.9335) },
    // Th
    Row { z: 90, a1: 0.2176, a2: 0.6240, alpha1: 26.554, alpha2: 4.0044, alpha3: Some(1.0238) },
    // Pa
    Row { z: 91, a1: 0.2413, a2: 0.6304, alpha1: 25.193, alpha2: 3.6780, alpha3: Some(0.9699) },
    // U
    Row { z: 92, a1: 0.2448, a2: 0.6298, alpha1: 25.252, alpha2: 3.6397, alpha3: Some(0.9825) },
];
