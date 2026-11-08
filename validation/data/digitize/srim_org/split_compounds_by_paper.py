#!/usr/bin/env python3
"""Write one `lindhard-stopping-data/1` file per (ion, compound, phase, paper).

    python3 -I split_compounds_by_paper.py <srim_org_compounds.csv> <bibliography.csv> <out_dir>
    python3 -I split_compounds_by_paper.py --self-test

Issue #525. Companion of `split_by_paper.py` for the measured compound points
of `srim_org_compounds.csv` (#498). Only the figures in `FIGURES` are used
(H in water, H and He in polystyrene). A row is kept when its `key_status` is
`exact` and its `citation_key` is a single key (no `|`): a `|`-joined,
`unconfirmed` or `corrected` key is never assigned to one paper. Groups with
fewer than `MIN_POINTS` rows, or whose key is not in the bibliography, are not
written.

Units: the water figure is per molecule (`eV/(1e15 molecules/cm2)`), so value
and sigma are divided by 3 to give the per-average-atom loader unit. The
polystyrene figures are already per atom. Targets are `H2O` and `C8H8`. The
energy is keV/amu times the ion mass (`convert_to_loader.MASS_AMU`) times 1000.
The loader has no phase column, so the phase is part of the file name and of
`name`/`description`. The output is a pure function of the inputs (rows keep
input order), so a rerun is byte-identical. File names:
`srim_org_<ion>_<target>[-<phase>]_<key>.csv`, target in lower case.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import convert_to_loader as c  # noqa: E402

MIN_POINTS = 5
# figure -> (formula, atoms per unit of the figure's value, expected value unit)
FIGURES = {
    "stop_h_water": ("H2O", 3, "eV/(1e15 molecules/cm2)"),
    "stop_h_polystyrene": ("C8H8", 1, "eV/(1e15 atoms/cm2)"),
    "stop_he_polystyrene": ("C8H8", 1, "eV/(1e15 atoms/cm2)"),
}


def per_atom(x, atoms):
    """Value or sigma per average atom, formatted as the files print it."""
    return f"{float(x) / atoms:.7g}"


def select(rows):
    """Group the usable rows by (figure, phase, key), in first-seen order."""
    groups = {}
    for r in rows:
        if (
            r["figure"] in FIGURES
            and r["origin"] == "digitized-compilation"
            and r["key_status"] == "exact"
            and r["citation_key"]
            and "|" not in r["citation_key"]
        ):
            if r["value_unit"] != FIGURES[r["figure"]][2]:
                raise SystemExit(f"unexpected unit {r['value_unit']!r} in {r['figure']}")
            groups.setdefault((r["figure"], r["phase"], r["citation_key"]), []).append(r)
    return groups


def system_tag(ion, formula, phase):
    return f"{ion.lower()}-{formula.lower()}" + (f"-{phase}" if phase else "")


def file_name(ion, formula, phase, key):
    tail = f"-{phase}" if phase else ""
    return f"srim_org_{ion.lower()}_{formula.lower()}{tail}_{key}.csv"


def header(ion, formula, phase, key, citation, n):
    where = f"{formula} ({phase})" if phase else formula
    return (
        "# format: lindhard-stopping-data/1\n"
        f"# name: srim-org-{system_tag(ion, formula, phase)}-{key.lower()}\n"
        f"# description: {ion} in {where}, {n} points of one paper ({citation}) "
        "read off the srim.org compilation plot; each has a half-pixel stopping "
        "uncertainty; the stopping is per average atom of the formula; the energy "
        "uncertainty (relative, in the compilation file) is not carried\n"
    )


def row_of(r, formula, atoms, key, ref):
    ion = r["ion"]
    mass = c.MASS_AMU[ion]
    return [
        c.ION_Z[ion],
        f"{mass:.6f}",
        formula,
        f"{float(r['energy']) * mass * 1000.0:.7g}",
        per_atom(r["value"], atoms),
        per_atom(r["u_value"], atoms),
        f"digitized plot, {r['method']}",
        f"Digitized from the srim.org compilation: {ref}",
        key,
        r["origin"],
    ]


def reference(key, bib):
    row = bib.get(key)
    if row is None:
        return None
    return f'{row["authors"]}, {row["journal"]}'


def write_all(rows, bib, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for (figure, phase, key), pts in select(rows).items():
        ref = reference(key, bib)
        if len(pts) < MIN_POINTS or ref is None:
            continue
        formula, atoms, _ = FIGURES[figure]
        ion = pts[0]["ion"]
        path = os.path.join(out_dir, file_name(ion, formula, phase, key))
        with open(path, "w", newline="", encoding="utf-8") as out:
            out.write(header(ion, formula, phase, key, ref, len(pts)))
            w = csv.writer(out, lineterminator="\n")
            w.writerow(c.COLUMNS)
            for r in pts:
                w.writerow(row_of(r, formula, atoms, key, ref))
        print(f"{len(pts)} points -> {path}")


def self_test():
    """The water conversion on one known row: 22.18 eV/(1e15 molecules/cm2) at
    28.43 keV/amu (1990-Mitt, H in water vapour) is 7.393333 per atom at
    28.65 keV."""
    r = {
        "figure": "stop_h_water", "ion": "H", "energy": "28.43", "value": "22.18",
        "u_value": "0.055", "method": "m", "origin": "digitized-compilation",
    }
    row = row_of(r, "H2O", 3, "1990-Mitt", "ref")
    assert row[2] == "H2O"
    assert row[4] == "7.393333", row[4]
    assert row[5] == "0.01833333", row[5]
    assert row[3] == f"{28.43 * c.MASS_AMU['H'] * 1000.0:.7g}"
    p = dict(r, figure="stop_h_polystyrene")
    assert row_of(p, "C8H8", 1, "k", "ref")[4] == "22.18"
    print("self-test ok")


def main(argv):
    if argv[1:] == ["--self-test"]:
        return self_test()
    if len(argv) != 4:
        raise SystemExit(__doc__)
    with open(argv[2], newline="", encoding="utf-8") as f:
        bib = {}
        for r in csv.DictReader(f):
            bib.setdefault(r["citation_key"], r)
    with open(argv[1], newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    write_all(rows, bib, argv[3])


if __name__ == "__main__":
    main(sys.argv)
