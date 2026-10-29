#!/usr/bin/env python3
"""Convert the 18-column digitizer output into a `lindhard-stopping-data/1` file.

    python3 -I convert_to_loader.py <srim_org_stopping.csv> <srim_org_bibliography.csv> <out.csv>

Issue #499. The loader (`docs/stopping-data-format.md`) reads a small fixed
layout; this script maps the digitizer's columns onto it and leaves the
digitizer's own output untouched:

* `ion` (H, D, He) becomes `ion_z` and `ion_mass_amu` (H 1.007825, D 2.014102,
  He 4.002603 u: atomic masses, as the loader's ion mass is).
* `energy_kev_per_amu` becomes `energy_ev` = keV/u x mass x 1000.
* `stopping_ev_cm2_per_1e15_atoms` is already in the loader's unit;
  `u_stopping` (an absolute half-pixel sigma in that unit) becomes the sigma.
  The relative energy uncertainty `u_energy_rel` has no loader column and is
  dropped.
* `citation_key` is the key of `citation_keys` when attribution is `unique`.
  For `shared-colour` and `partial-legend` rows it keeps every candidate
  key, joined with `|` as in the source, in the source order; no candidate is
  picked, because the point belongs to one of them (or, for `partial-legend`,
  possibly to a paper not listed). `citation` spells the candidates out.
* `origin` is `digitized-compilation`.
* Rows of any other origin (the panel-recovered `residual-reconstructed`
  points) and ions without a mass above (Li) are skipped, not converted: the
  loader's origin set is closed, and no Li mass is chosen here.

Row order is the input order; the output is a pure function of the two inputs
(no clock, no environment), so a rerun is byte-identical.
"""
import csv
import sys

MASS_AMU = {"H": 1.007825, "D": 2.014102, "He": 4.002603}
ION_Z = {"H": 1, "D": 1, "He": 2}
COLUMNS = (
    "ion_z,ion_mass_amu,target,energy_ev,s_e_ev_1e15_cm2,sigma_ev_1e15_cm2,"
    "method,citation,citation_key,origin"
).split(",")
HEADER = """\
# format: lindhard-stopping-data/1
# name: srim-org-h-he-digitized
# description: H, D and He stopping points read off the plots of the srim.org measured-data compilation (not copied from a table); each point has a half-pixel uncertainty and the candidate papers named in citation_key
"""


def reference(key, bibliography):
    row = bibliography.get(key)
    if row is None:
        return f"{key} (paper not resolved)"
    return f'{row["authors"]}, {row["journal"]}'


def convert(source, bibliography_path, out):
    with open(bibliography_path, newline="", encoding="utf-8") as f:
        bibliography = {r["key"]: r for r in csv.DictReader(f)}
    with open(source, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out.write(HEADER)
    w = csv.writer(out, lineterminator="\n")
    w.writerow(COLUMNS)
    n = 0
    for r in rows:
        if r["origin"] != "digitized-compilation" or r["ion"] not in MASS_AMU:
            continue
        mass = MASS_AMU[r["ion"]]
        keys = r["citation_keys"].split("|")
        attribution = r["attribution"]
        refs = "; ".join(reference(k, bibliography) for k in keys)
        if attribution == "unique":
            citation = f"Digitized from the srim.org compilation: {refs}"
        else:
            citation = (
                f"Digitized from the srim.org compilation ({attribution}); "
                f"one of: {refs}"
            )
        w.writerow(
            [
                ION_Z[r["ion"]],
                f"{mass:.6f}",
                r["target"],
                f"{float(r['energy_kev_per_amu']) * mass * 1000.0:.7g}",
                r["stopping_ev_cm2_per_1e15_atoms"],
                r["u_stopping"],
                f"digitized plot, {r['method']}",
                citation,
                "|".join(keys),
                r["origin"],
            ]
        )
        n += 1
    return n


def main(argv):
    if len(argv) != 4:
        raise SystemExit(__doc__)
    with open(argv[3], "w", newline="", encoding="utf-8") as out:
        n = convert(argv[1], argv[2], out)
    print(f"{n} points written to {argv[3]}")


if __name__ == "__main__":
    main(sys.argv)
