# Experimental stopping and range papers: bibliography

`bibliography.csv` lists the experimental stopping-power and range papers
that srim.org's citation files record, one row per paper (1230 rows). Titles,
authors, journal references and years are facts about published papers; the
compiler's comment says what each paper measured. No measured value is in this
directory.

Released under CC0 1.0 (`LICENSE`), separate from the MIT code of the rest of
the repository, as `CONTRIBUTING.md` ("The srim.org experimental-data
compilation") allows. The citation PDFs and page text are not copied into the
tree.

## Where it came from

srim.org (unmaintained; its author has died) links one citation file per
ion-target pair from `PlotCit.htm` (`Citations/Cit <ZZ>-<ZZ>.pdf`, 202 files,
fetched 2026-10-29) and one per compound from `SRIM/Compounds.htm`
(`SRIM/Compounds/Citations/<compound>.pdf`, 148 files). `fetch_bibliography.py`
fetches them (at most one request per second, into a scratch directory outside
the tree), runs `pdftotext -layout` and parses the entries. The `Cit <Z>-All`
and `Cit All-<Z>` aggregates repeat per-pair entries; they are not a source of
rows and serve only as a cross-check.

## Columns

| column | meaning |
|---|---|
| `citation_key` | the compiler's key, exactly as printed (for example `1962-Ewin`); not unique, since 73 keys are used for more than one paper |
| `year` | publication year as printed |
| `authors` | author list as printed, joined with single spaces |
| `title` | title; text order comes from the PDF layout and can be scrambled in a few entries |
| `journal` | journal reference as printed |
| `comment` | the compiler's note: quantity (S stopping, R range, d prefix for straggling), energy range, ions and targets |
| `old_ref` | the old reference number; empty in the compound files, which do not print one |
| `pairs` | ion-target pairs the paper is listed under, as the file codes `<ion Z>-<target Z>` (`01-14` is H in Si), joined with `;` |
| `target` | compound names the paper is listed under, exactly as the compound file names them (`Al2O3`, `B4C`, `Air-H`), joined with `;`; empty for papers listed only under elemental pairs |
| `doi` | Crossref DOI, filled only when exactly one candidate matches the normalized title, the year, the journal and, where the citation prints one, the first page; empty otherwise (275 of 1230 filled). Never guessed |

A paper is one row (same key and normalized title), however many pairs and
compounds list it, so `pairs` and `target` can both be filled.

## Checks

- Per file, the number of `Comment :` lines (one per entry), the number of
  entry-start lines and the number of parsed plus unparsed entries must agree;
  the script exits non-zero otherwise. All 350 files pass.
- Entries that do not parse are written verbatim to `bibliography-unparsed.txt`
  instead of being dropped: two, the same paper (`1975-Tsai`)
  listed in two compound files, whose title the PDF layout interleaves with other text.
- No PDF is in the tree. `test_fetch_bibliography.py` tests the
  parser on a synthetic fixture (`python3 -B -m unittest test_fetch_bibliography`).

## Rebuilding

    python3 fetch_bibliography.py            # fetch, parse, Crossref lookup (about 40 min)
    python3 fetch_bibliography.py --no-fetch # reuse the scratch directory

Requires poppler's `pdftotext`. The scratch directory defaults to
`$TMPDIR/lindhard-srim-bibliography`; the Crossref answers are cached there, so
an interrupted run resumes.
