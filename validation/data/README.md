# Experimental reference data (level 3)

Published measurements that `lindhard` is compared against
([`docs/validation.md`](../../docs/validation.md) section 3). Everything here
is a measured fact cited to its paper. **Never** SRIM/TRIM output, tables
derived from it, or any other code's output, however well known
([`CONTRIBUTING.md`](../../CONTRIBUTING.md), clean-room tiers).

## Status

- **B in amorphous Si:** `ranges/b_{1,2,3,5,10,20}keV_asi_wach1982.json`,
  measured ranges of Wach and Wittmaack (1982) digitized from the compilation
  of Wittmaack and Mutzke (2017), Fig. 8 (#51). The script is in
  `digitize/`.
- **P and As in amorphous Si:** open gaps (#51). What was searched, and why
  the one source found could not be used, is in
  [`docs/data-provenance.md`](../../docs/data-provenance.md).

A dataset is added only when its numbers can be read from an accessible
source (a table, or a figure digitized as described below) whose terms allow
storing them. A dataset recalled from memory, or copied from a secondary
compilation whose own source is unclear, does not go in. A secondary source
is acceptable when it names the measurement and holds measured values only
(no SRIM/TRIM or other code output mixed in); cite both, and say which was
seen.

## Layout

```
validation/data/ranges/<id>.json      one dataset: metadata + moments
validation/data/ranges/<id>.csv       optional depth profile (depth_nm,value)
validation/data/digitize/             scripts that digitized a figure (not run
                                      by the harness; may need numpy/Pillow)
```

`validation/experiments/run.py` reads every `ranges/*.json`, checks it (see
"Enforced" below), runs `lindhard` on the same ion, energy, tilt and target,
and writes
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
  "short_citation": "Author (Year)",
  "doi": "10.xxxx/...",
  "primary_doi": "10.xxxx/... (if read from a secondary source)",
  "url": "https://... (where the source was read)",
  "source_location": "Table 2, row 3 (or Fig. 4, open squares)",
  "extraction": "transcribed from a table | digitized from a figure: tool, axis calibration, estimated digitizing uncertainty",
  "terms": "Facts, cited",
  "added": "YYYY-MM-DD"
}
```

Rules:

- `drp_nm` and `drp_unc_nm` are both `null` when the source gives no ΔRp;
  the level-3 table then shows only lindhard's value.
- `rp_unc_nm` and `drp_unc_nm` include the paper's stated uncertainty (for
  example the SIMS depth-scale calibration) and, for digitized values, the
  digitizing uncertainty. If the paper states none, say so in `extraction`
  and estimate it there; do not leave it zero.
- Moments must be the paper's own (or computed by us from the paper's
  profile, saying so in `extraction`). Profiles above a dose where the
  target state changes, or with a surface/oxide artefact, are flagged in
  `target_state`.
- Each dataset gets a row in [`docs/data-provenance.md`](../../docs/data-provenance.md)
  in the same PR, naming its `id`, with the source's licence or terms. If the
  terms do not allow storing the numbers, store none: record the citation and
  a note on how to fetch and derive them instead.

## Enforced

`validation/experiments/run.py` (and `--check`, which `validation/run.sh`
always runs) exits with an error naming the file and field if a dataset:

- lacks a non-empty `id` (equal to the file name), `ion`, `target`,
  `target_state`, `citation`, `source_location`, `extraction`, `terms` or
  `added`, or both of `doi` and `url`;
- has no positive `energy_ev`, `measured.rp_nm` or `measured.rp_unc_nm`, no
  `measured.method`, or a ΔRp without a positive uncertainty;
- has no row in `docs/data-provenance.md` that names its `id`.
