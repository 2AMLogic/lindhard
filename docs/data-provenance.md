# Data provenance

Every numeric dataset in this tree gets a row here: element masses, fit
coefficients, optical data, experimental points. A dataset without a row is a
CI-review blocker.

The allowed origins are:

1. **Computed from a published formula.** Cite the paper; the code is ours.
2. **Experimental measurement.** Cite the paper or database record.
   Measurements are facts; we cite them.
3. **A permissive source** (Tier A in [`../CONTRIBUTING.md`](../CONTRIBUTING.md)),
   with its notice kept in `THIRD_PARTY_LICENSES.md`.

The following are **never allowed**:

- SRIM tables, including any that other projects pass along.
- ICRU report tables.
- NIST SRD tables, until their terms have been reviewed and recorded here.
- Anything whose origin is not known.

| Dataset | Path | Origin | Citation | Licence / terms | Added |
|---|---|---|---|---|---|
| Standard atomic weights, Z=1..92 | `lindhard/src/elements.rs` | Measurement compilation (IUPAC) | CIAAW standard atomic weights 2021 (<https://www.ciaaw.org>); Meija et al., Pure Appl. Chem. 88, 265 (2016). Conventional values for interval elements; mass number of longest-lived isotope for Tc, Pm, Po, At, Rn, Fr, Ra, Ac (no standard weight; flagged `weight_is_mass_number`). Not yet verified against the cited table; spot-check before relying on it | Facts, cited | 2026-10-04 |
| Elemental solid densities | `lindhard/src/elements.rs` | Experimental / handbook compilation | Haynes (ed.), CRC Handbook of Chemistry and Physics, "Physical Constants of the Elements", about 20 C; `None` for gases, liquids and elements without data. Values entered from memory of the handbook by the contributor and not yet verified: every value needs a spot-check against the cited edition before it is relied on. Reviewer-flagged as possibly wrong: Sc (2.989 here; CRC may list 2.985) and Ra (5.0 here; CRC may list 5.5) | Facts, cited | 2026-10-04 |
| Default surface binding energy `E_s` (cohesive energy) | `lindhard/src/elements.rs` | Experimental (thermochemistry) | Kittel, Introduction to Solid State Physics, 8th ed., ch. 3 cohesive energies; use as `E_s` is the Sigmund 1969 convention (Phys. Rev. 184, 383). Convention, not a measured surface barrier. Missing elements are unset by design. Not yet verified against the cited edition; spot-check before relying on it | Facts, cited | 2026-10-04 |
| Default displacement energy `E_d` | `lindhard/src/elements.rs` | Convention | ASTM E521 recommended values (the NRT convention; Norgett, Robinson, Torrens, Nucl. Eng. Des. 33, 50 (1975)). A convention for damage accounting; not yet verified against the standard; spot-check each value before relying on it. Missing elements are unset by design | Convention, cited | 2026-10-04 |
| Default lattice binding energy `E_b` | `lindhard/src/material.rs` | Convention | 0 eV, the common amorphous-target BCA convention; no data | None | 2026-10-04 |
| Physical constants | `lindhard/src/constants.rs` | Measurement compilation (CODATA) | CODATA 2022, Tiesinga et al., Rev. Mod. Phys. 97, 025002 (2025). Not yet verified digit by digit against the published table | Facts, cited | 2026-10-04 |
| Lindhard-Scharff stopping constants: 0.8853 (screening length), 0.0793 (`k_L`), `ξ_e = Z1^(1/6)` | `lindhard/src/ion/stopping/lindhard_scharff.rs` | Computed from a published formula | Lindhard & Scharff, Phys. Rev. 124, 128 (1961); Lindhard, Scharff, Schiott, Mat. Fys. Medd. 33 (14) (1963). Checked against the dimensional LS form in a unit test (1 %); not yet checked digit by digit against the papers | Formula, cited | 2026-10-04 |
| Oen-Robinson local-loss constants: `c = 0.3`, Firsov-type length `0.8853 a0 (Z1^(1/2)+Z2^(1/2))^(-2/3)` | `lindhard/src/ion/stopping/oen_robinson.rs` | Computed from a published formula | Oen & Robinson, NIM 132, 647 (1976). Entered from the contributor's memory of the paper and **not verified**; spot-check against the paper before relying on it | Formula, cited | 2026-10-04 |
| Bethe-Bloch: Bloch rule `I = 10 eV · Z`, Bloch series, Barkas empirical effective charge `125 β / Z1^(2/3)` | `lindhard/src/ion/stopping/bethe.rs` | Computed from a published formula | Bloch, Ann. Phys. 408, 285 (1933) (rule and series); Barkas, Nuclear Research Emulsions I (1963) and Northcliffe, Ann. Rev. Nucl. Sci. 13, 67 (1963) (effective charge). Not verified against the originals. No shell, density-effect or Barkas-term coefficients are used | Formula, cited | 2026-10-04 |
| Bohr straggling, relativistic factor `(1 - β²/2)/(1 - β²)` | `lindhard/src/ion/stopping/straggling.rs` | Computed from a published formula | Bohr, Mat. Fys. Medd. 18 (8) (1948); Bethe & Livingston, Rev. Mod. Phys. 9, 245 (1937). Not verified against the originals | Formula, cited | 2026-10-04 |
| Proton-in-Si Bethe hand-calculation values (1, 10, 100 MeV) | `lindhard/src/ion/stopping/bethe.rs` (tests) | Computed from the Bethe formula with `I = 140 eV` and CODATA constants | Our own arithmetic; not taken from PSTAR/ICRU/SRIM | Computed | 2026-10-04 |
| Screening function coefficients (ZBL universal, Kr-C, Moliere, Lenz-Jensen) | `lindhard/src/ion/potential.rs` | Computed from published formulas | ZBL: Ziegler, Biersack, Littmark (1985) Ch. 2; Kr-C: Wilson, Haggmark, Biersack, Phys. Rev. B 15, 2458 (1977); Moliere: Z. Naturforsch. A 2, 133 (1947); Lenz-Jensen: Lenz, Z. Phys. 77, 713 (1932), Jensen, Z. Phys. 77, 722 (1932), polynomial form as in Nastasi, Mayer, Hirvonen (1996). Entered from the contributor's memory of the papers, not checked digit by digit against them: verify all coefficients. The Lenz-Jensen argument scale `q = sqrt(9.67 x)` is the least certain (small-x slope is 1.60 vs the Thomas-Fermi 1.588); check it first | Formulas, cited | 2026-10-04 |
| Screening lengths (universal, Firsov, Lindhard) | `lindhard/src/ion/potential.rs` | Computed from published formulas | Universal: ZBL (1985) Ch. 2 (0.8854 a0 / (Z1^0.23 + Z2^0.23)); Firsov: Sov. Phys. JETP 6, 534 (1958); Lindhard: Lindhard, Scharff, Schiott, Mat. Fys. Medd. Dan. Vid. Selsk. 33(14) (1963). 0.8853 is computed as (9 pi^2/128)^(1/3). Not yet verified against the papers | Formulas, cited | 2026-10-04 |
| Magic-formula constants C1..C5 (ZBL and Moliere sets) | `lindhard/src/ion/scattering.rs` | Computed from published formula | Biersack and Haggmark, Nucl. Instrum. Methods 174 (1980) 257; ZBL (1985) Ch. 2. Entered from memory, not verified against the paper. Sanity check only: the ZBL set reproduces the quadrature `cos(theta/2)` to 0.025 over the tested grid (`magic_formula_tracks_quadrature_for_zbl`); the Moliere set is untested | Formulas, cited | 2026-10-04 |
| ZBL universal reduced nuclear stopping fit (test reference only) | `lindhard/src/ion/scattering.rs` (tests) | Published fit formula | ZBL (1985) Ch. 2: `ln(1+1.1383 eps) / (2 (eps + 0.01321 eps^0.21226 + 0.19593 eps^0.5))`, `eps <= 30`; `ln(eps)/(2 eps)` above. Entered from memory, not verified. Our quadrature agrees to 0.3 % for eps <= 0.3 and deviates up to 5.5 % at eps around 30 to 100: either the recalled fit or the published fit's high-energy accuracy is at fault; verify the coefficients against the book | Formula, cited | 2026-10-04 |

## Open questions

- IAEA stopping database reuse terms: unread (bot wall on 2026-10-04).
  Resolve before ingesting points.
- NIST ESTAR/PSTAR terms under 15 USC 290e: unresolved.
- Optical data for dielectric-function models (Si, SiO₂, resists, metals):
  identify sources whose terms allow redistribution.
