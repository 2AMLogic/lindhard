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
| *(none yet)* | | | | | |

## Open questions

- IAEA stopping database reuse terms: unread (bot wall on 2026-10-04).
  Resolve before ingesting points.
- NIST ESTAR/PSTAR terms under 15 USC 290e: unresolved.
- Optical data for dielectric-function models (Si, SiO₂, resists, metals):
  identify sources whose terms allow redistribution.
