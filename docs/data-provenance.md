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
| Standard atomic weights, Z=1..92 | `lindhard/src/elements.rs` | Measurement compilation (IUPAC) | CIAAW standard atomic weights 2021 (<https://www.ciaaw.org>); Meija et al., Pure Appl. Chem. 88, 265 (2016). Conventional values for interval elements; mass number of longest-lived isotope for Tc, Pm, Po, At, Rn, Fr, Ra, Ac (no standard weight; flagged `weight_is_mass_number`) | Facts, cited | 2026-10-04 |
| Elemental solid densities | `lindhard/src/elements.rs` | Experimental / handbook compilation | Haynes (ed.), CRC Handbook of Chemistry and Physics, "Physical Constants of the Elements", about 20 C; `None` for gases, liquids and elements without data. Values entered from the handbook by the contributor; a reviewer should spot-check against the current edition | Facts, cited | 2026-10-04 |
| Default surface binding energy `E_s` (cohesive energy) | `lindhard/src/elements.rs` | Experimental (thermochemistry) | Kittel, Introduction to Solid State Physics, 8th ed., ch. 3 cohesive energies; use as `E_s` is the Sigmund 1969 convention (Phys. Rev. 184, 383). Convention, not a measured surface barrier. Missing elements are unset by design | Facts, cited | 2026-10-04 |
| Default displacement energy `E_d` | `lindhard/src/elements.rs` | Convention | ASTM E521 recommended values (the NRT convention; Norgett, Robinson, Torrens, Nucl. Eng. Des. 33, 50 (1975)). A convention for damage accounting; reviewer to verify each value against the standard. Missing elements are unset by design | Convention, cited | 2026-10-04 |
| Default lattice binding energy `E_b` | `lindhard/src/material.rs` | Convention | 0 eV, the common amorphous-target BCA convention; no data | None | 2026-10-04 |
| Physical constants | `lindhard/src/constants.rs` | Measurement compilation (CODATA) | CODATA 2022, Tiesinga et al., Rev. Mod. Phys. 97, 025002 (2025) | Facts, cited | 2026-10-04 |

## Open questions

- IAEA stopping database reuse terms: unread (bot wall on 2026-10-04).
  Resolve before ingesting points.
- NIST ESTAR/PSTAR terms under 15 USC 290e: unresolved.
- Optical data for dielectric-function models (Si, SiO₂, resists, metals):
  identify sources whose terms allow redistribution.
