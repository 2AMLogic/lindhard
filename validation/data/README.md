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
- **Sputter yields, Ar on Cu:** `sputtering/ar_cu_sputter_*.json`, 15
  measured sets (47 points, 0.2 to 10 keV, normal incidence) digitized from
  the compilation of Matsunami et al., IPPJ-AM-32 (1983), and cross-checked
  against Yamamura and Tawara, NIFS-DATA-23 (1995) (#69). Stored under an
  operator decision (#69) despite the compilations' cover notes; no figure
  image or fitted curve is stored.
- **Sputter yields, Ar on Si, Ag and Au:** `sputtering/ar_si_sputter_*.json`
  (4 sets, 19 points), `ar_ag_sputter_*.json` (9 sets, 34 points) and
  `ar_au_sputter_*.json` (12 sets, 44 points), same compilation, range and
  decision (#70). Every point has a manual second read
  (`digitize/secondread_ar_*_am32.json`), and the symbols that NIFS-DATA-23
  separates a second compilation read (`digitize/crosscheck_ar_*_nifs23.json`).
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
validation/data/sputtering/<id>.json  one measured sputter-yield set: metadata + points
validation/data/digitize/             scripts that digitized a figure (not run
                                      by the harness; may need numpy/Pillow/scipy),
                                      and their cross-check records
```

`validation/experiments/run.py` reads every `ranges/*.json` and
`sputtering/*.json`, checks it (see
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

## Sputter-yield schema (`sputtering/<id>.json`)

One file per *original measurement set* (one reference, one ion-target
pair), holding all its points. `id` is the file name, e.g.
`ar_cu_sputter_<firstauthor><year>`.

```json
{
  "kind": "sputter_yield",
  "id": "ar_cu_sputter_<firstauthor><year>",
  "ion": "Ar", "mass_amu": null,
  "target": "Cu",
  "target_state": "as stated by the source, or 'not stated' and why",
  "incidence_deg": 0.0,
  "points": [
    { "energy_ev": 1000.0, "energy_unc_rel": 0.02, "yield": 2.0, "yield_unc_rel": 0.02,
      "note": "optional, e.g. why a symbol identity is inferred" }
  ],
  "original_reference": "Authors, Journal Volume, Page (Year)",
  "original_doi": "10.xxxx/... or null",
  "compilation": "the secondary source that was read",
  "compilation_figure": "Ar -> Cu",
  "compilation_pdf_page": 118,
  "compilation_symbol": "letter used for this reference in the figure",
  "compilation_table_ref": "where the compilation lists the reference",
  "url": "where the compilation was read",
  "crosscheck": { "compilation": "...", "figure": "...", "pdf_page": 0, "symbol": "...", "agreement": "...",
                  "second_read": "optional: a second symbol-location method in the same figure" },
  "extraction": "tool, axis calibration, symbol-location method, uncertainty",
  "reliability_note": "why the compilation treats the set as reliable",
  "terms": "Facts, cited; ...",
  "added": "YYYY-MM-DD"
}
```

Rules:

- Measured points only: never a compilation's fitted curve, and never a
  point the compilation marks as calculated (e.g. ACAT in NIFS-DATA-23).
- `energy_unc_rel` and `yield_unc_rel` are the digitizing uncertainties,
  never zero; say in `extraction` that the source gives no measurement
  uncertainty when it gives none.
- `incidence_deg` is required (normal incidence = 0); a later angular
  extension must not be mistaken for normal incidence.
- A symbol whose identity cannot be decided is not stored; the digitizing
  script lists it and why.
- Every set digitized from IPPJ-AM-32 is re-read independently from
  NIFS-DATA-23 where that figure separates it; the record is
  `digitize/crosscheck_ar_<target>_nifs23.json`, summarized in each
  `crosscheck` (`agreement`). A point the second compilation cannot separate
  is located a second time in the same figure by another method (a glyph box
  read by eye against the template centre); the record is
  `digitize/secondread_ar_<target>_am32.json`, summarized in
  `crosscheck.second_read` (Si, Ag and Au, #70).
- A second read beyond the combined uncertainty is either explained by a rule
  written in the digitizing script and replayable on the committed record, or
  reported as unexplained; beyond twice the combined uncertainty the point
  gets `"flag": "disagrees_between_compilations"` and is left out of the
  statistics (`energy_groups` in `validation/experiments/run.py`). No stored
  point is flagged at present.

## Enforced

`validation/experiments/run.py` (and `--check`, which `validation/run.sh`
always runs) exits with an error naming the file and field if a dataset:

- lacks a non-empty `id` (equal to the file name), `ion`, `target`,
  `target_state`, `citation`, `source_location`, `extraction`, `terms` or
  `added`, or both of `doi` and `url`;
- has no positive `energy_ev`, `measured.rp_nm` or `measured.rp_unc_nm`, no
  `measured.method`, or a ΔRp without a positive uncertainty;
- has no row in `docs/data-provenance.md` that names its `id`.

and, for a sputter-yield dataset, if it:

- has a `kind` other than `sputter_yield`, or lacks a non-empty `id` (equal
  to the file name), `ion`, `target`, `target_state`, `original_reference`,
  `compilation`, `compilation_figure`, `extraction`, `terms` or `added`, or
  both of `original_doi` and `url`;
- has no positive `compilation_pdf_page`, no numeric `incidence_deg`, no
  points, or a point without a positive `energy_ev`, `yield`,
  `energy_unc_rel` or `yield_unc_rel`;
- has no row in `docs/data-provenance.md` that names its `id`.

`validation/experiments/test_run.py` (run by `validation/run.sh`) checks
that each of these failures is caught.
