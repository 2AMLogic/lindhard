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
- **Optical ELF, Al, Cu, glassy C and Au:**
  `optical/{al,cu,c,au}_elf_hagemann1975.toml`, ELF converted by
  `optical/ingest_hagemann.py` from Hagemann, Gudat and Kunz (1975) as
  transcribed in the CC0 refractiveindex.info database (Al and Cu #98, C and
  Au #148), with a second read against the scanned tables of DESY report
  SR-74/7 (`optical/secondread_hagemann1975.json`). Sum-rule test:
  `lindhard/tests/optical_sumrule.rs` (gated on the linear interpolation `elf()` serves: Al N_eff and P_eff,
  Cu N_eff, C N_eff and Au N_eff miss 5 % and are pinned with their
  explanation; power-law segments are kept as a labelled comparison, where C
  N_eff also misses). Si: open gap.
- **Electron backscatter coefficients, C, Al, Si, Cu and Au:**
  `backscatter/eta_<el>_<author><year>.json`, 57 measured sets (one per
  original measurement and element, 380 points, 0.1 to 102 keV), transcribed
  from the tables of D. C. Joy's database (revision 01-01, 2001, read in its
  Internet Archive copy) by `backscatter/ingest_joy2001.py`, and re-read from
  the independent Akbari (2022) database
  (`backscatter/crosscheck_akbari2022.json`) (#148).
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
validation/data/backscatter/<id>.json one measured electron backscatter set: metadata + points
validation/data/optical/<id>.toml     measured optical energy-loss functions (ELF)
                                      and the scripts that ingested them
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

## Backscatter-coefficient schema (`backscatter/<id>.json`)

One file per *original measurement set* and target, holding all its points;
sets of the same (E, Z) are never merged or averaged in the files.
`validation/experiments/backscatter.py` takes the median at run time.

```json
{
  "kind": "backscatter_coefficient",
  "id": "eta_<el>_<firstauthor><year>",
  "target": "Al", "z": 13,
  "target_state": "as stated by the source, or 'not stated'",
  "incidence_deg": 0.0,
  "incidence_note": "what the source says about the geometry",
  "eta_threshold_note": "the SE/BSE energy split, if the source states it",
  "points": [ { "energy_ev": 10000.0, "eta": 0.15, "eta_unc_abs": 0.005 } ],
  "original_reference": "Authors, Journal Volume, Page (Year)",
  "original_doi": "10.xxxx/... or null",
  "original_reference_note": "Crossref check, discrepancies, 'not opened'",
  "compilation": "the secondary source that was read",
  "compilation_set": "where in it (heading, data set number)",
  "compilation_reference_number": 2,
  "url": "where the compilation was read",
  "crosscheck": { "compilation": "...", "agreement": "...", "record": "..." },
  "extraction": "how the numbers were taken, and what the uncertainty is",
  "terms": "Facts, cited; ...",
  "added": "YYYY-MM-DD"
}
```

Rules:

- `eta_unc_abs` is never zero. When the source gives no measurement
  uncertainty, it is half a unit of the last printed digit and `extraction`
  says so.
- `incidence_deg` is required, as for sputter yields; only normal-incidence
  sets enter the comparison.
- A set whose original measurement cannot be named is not stored (the
  Akbari-only sets, `docs/data-provenance.md`).

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

and, for a backscatter dataset (`validation/experiments/backscatter.py
--check`), if it has a `kind` other than `backscatter_coefficient`, lacks a
non-empty `id` (equal to the file name), `target`, `target_state`,
`original_reference`, `compilation`, `compilation_set`, `extraction`, `terms`
or `added`, or both of `original_doi` and `url`, has no positive integer `z`,
no numeric `incidence_deg`, no points, a point without a positive
`energy_ev` or `eta_unc_abs` or with `eta` outside (0, 1), or no row in
`docs/data-provenance.md` that names its `id`.

`validation/experiments/test_run.py` (run by `validation/run.sh`) checks
that each of these failures is caught.

What the checker does **not** enforce (not schema rules at present): the
nominal plotting range of a figure (the stored points span 196.3 to
10 020 eV, and seven lie just outside 200-10 000 eV), any
particular wording of free-text fields such as `terms`, `target_state` or
`original_reference` (only that the required ones are non-empty), the
presence of `crosscheck` or `reliability_note`, and the absence of
additional top-level keys. Tightening any of these needs its own schema
decision, not an incidental check.

## Sensitivity scenarios (sputter yields, #78)

The level-3 sputter summary in `docs/validation.md` compares lindhard with
the median of **every** stored set; that is the baseline and it does not
change. Below it, a generated table repeats the comparison with sets whose
attribution or target state is doubtful left out. The scenarios are the
named table `SENSITIVITY_SCENARIOS` in `validation/update_docs.py`: each
row names existing dataset `id`s of one target and the place where the
caveat is stored (the dataset's `target_state` or `original_reference`, or
the "Titles that name another system" paragraph of `docs/validation.md`).
Nothing is inferred from prose, and the datasets carry no exclusion marker.
`update_docs.py` stops with a message naming the scenario if an `id` is
unknown, listed twice, or belongs to another target, if a caveat source is
not one of those three, or if two scenarios are identical.

For each scenario the remaining sets of the target are filtered **first**,
then regrouped by the unchanged `energy_groups` rule (2 % grouping, normal
incidence, flagged points left out, median). Leaving a set out can
therefore change a group's median and its representative energy, or remove
the group. The code yields are the committed ones (`results.json`, and the
RustBCA summaries where present): exact at a run energy, otherwise
interpolated linearly in log E - log Y between the two bracketing run
energies, for positive yields only. That interpolation approximates the
committed code curve; it is not a new simulation and not an uncertainty
estimate. Precomputed ratios are never interpolated and nothing is
extrapolated: an energy outside the runs is listed as unsupported and the
row is marked partial; a scenario with no remaining energies, or a code
that was not run, is reported as unavailable. The no-exclusion row is
recomputed by the same path and must reproduce the committed baseline to
its printed precision, or `update_docs.py` stops.
