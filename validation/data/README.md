# Experimental reference data (level 3)

Published measurements that `lindhard` is compared against
([`docs/validation.md`](../../docs/validation.md) section 3). Everything here
is a measured fact cited to its paper. **Never** SRIM/TRIM output, tables
derived from it, or any other code's output, however well known
([`CONTRIBUTING.md`](../../CONTRIBUTING.md), clean-room tiers).

## Status

**No datasets yet.** The first range datasets (B, P, As in amorphous Si) are
tracked by #51. A dataset is added only when its
numbers can be read from the cited paper itself (a table, or a figure
digitized as described below). A dataset recalled from memory or copied from
a secondary compilation whose own source is unclear does not go in.

## Layout

```
validation/data/ranges/<id>.json      one dataset: metadata + moments
validation/data/ranges/<id>.csv       optional depth profile (depth_nm,value)
```

`validation/experiments/run.py` reads every `ranges/*.json`, runs `lindhard`
on the same ion, energy, tilt and target, and writes
`validation/experiments/results.json`; `validation/run.sh` turns that into the
level-3 table in `docs/validation.md`.

## Dataset schema (`ranges/<id>.json`)

```json
{
  "id": "b_5keV_asi_<firstauthor><year>",
  "ion": "B",
  "mass_amu": 11.009,
  "energy_ev": 5000.0,
  "tilt_deg": 7.0,
  "target": "Si",
  "target_state": "amorphous: pre-amorphized by Si+ or Ge+ implantation (state the dose)",
  "measured": {
    "method": "SIMS",
    "rp_nm": 0.0,
    "rp_unc_nm": 0.0,
    "drp_nm": 0.0,
    "drp_unc_nm": 0.0,
    "dose_cm2": 0.0
  },
  "profile_csv": null,
  "citation": "Authors, Journal Volume, Page (Year)",
  "doi": "10.xxxx/...",
  "source_location": "Table 2, row 3 (or Fig. 4, open squares)",
  "extraction": "transcribed from a table | digitized from a figure: tool, axis calibration, estimated digitizing uncertainty",
  "terms": "Facts, cited",
  "added": "YYYY-MM-DD"
}
```

Rules:

- `rp_unc_nm` and `drp_unc_nm` include the paper's stated uncertainty (for
  example the SIMS depth-scale calibration) and, for digitized values, the
  digitizing uncertainty. If the paper states none, say so in `extraction`
  and estimate it there; do not leave it zero.
- Moments must be the paper's own (or computed by us from the paper's
  profile, saying so in `extraction`). Profiles above a dose where the
  target state changes, or with a surface/oxide artefact, are flagged in
  `target_state`.
- Each dataset gets a row in [`docs/data-provenance.md`](../../docs/data-provenance.md)
  in the same PR.
