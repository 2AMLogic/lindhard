#!/usr/bin/env python3
"""Write one `lindhard-stopping-data/1` file per (ion, target, paper).

    python3 -I split_by_paper.py <srim_org_stopping.csv> <srim_org_bibliography.csv> <out_dir>

Issue #488. A file holds the points of one paper in one ion-target system, so
its residuals can be read per paper. Only rows that `convert_to_loader.py`
would convert and whose attribution is `unique` are used: a point whose colour
is shared by several candidate papers is never assigned to one of them. Groups
with fewer than `MIN_POINTS` such points, or whose key is not in the
bibliography (`ref<N>`), are not written.

Only the systems in `SYSTEMS` (H and He in Si, Cu and Au) are considered.
Column mapping, masses and the citation text come from `convert_to_loader.py`.
The output is a pure function of the two inputs (rows keep input order), so a
rerun is byte-identical. File names: `srim_org_<ion>_<target>_<key>.csv`.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import convert_to_loader as c  # noqa: E402

SYSTEMS = [(i, t) for i in ("H", "He") for t in ("Si", "Cu", "Au")]
MIN_POINTS = 5
KEYS = {
    ("H", "Si"): ["1977-Cemb", "1976-Grah", "1988-Mert"],
    ("He", "Si"): ["1973-Lin", "1973-Zieg", "1981-Pear"],
    ("He", "Cu"): ["1977-Ande3", "1975-Nomu"],
}


def header(ion, target, key, citation, n):
    return (
        "# format: lindhard-stopping-data/1\n"
        f"# name: srim-org-{ion.lower()}-{target.lower()}-{key.lower()}\n"
        f"# description: {ion} in {target}, {n} points of one paper ({citation}) "
        "read off the srim.org compilation plot; each has a half-pixel stopping "
        "uncertainty; the energy uncertainty (relative, in the compilation file) "
        "is not carried\n"
    )


def main(argv):
    if len(argv) != 4:
        raise SystemExit(__doc__)
    with open(argv[2], newline="", encoding="utf-8") as f:
        bib = {r["key"]: r for r in csv.DictReader(f)}
    with open(argv[1], newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    os.makedirs(argv[3], exist_ok=True)
    for ion, target in SYSTEMS:
        groups = {}
        for r in rows:
            if (
                r["origin"] == "digitized-compilation"
                and r["attribution"] == "unique"
                and r["ion"] == ion
                and r["target"] == target
            ):
                groups.setdefault(r["citation_keys"], []).append(r)
        for key in KEYS.get((ion, target), []):
            pts = groups.get(key, [])
            if len(pts) < MIN_POINTS or key not in bib:
                raise SystemExit(f"{ion} in {target}, {key}: not enough points")
            mass = c.MASS_AMU[ion]
            path = os.path.join(
                argv[3], f"srim_org_{ion.lower()}_{target.lower()}_{key}.csv"
            )
            ref = c.reference(key, bib)
            with open(path, "w", newline="", encoding="utf-8") as out:
                out.write(header(ion, target, key, ref, len(pts)))
                w = csv.writer(out, lineterminator="\n")
                w.writerow(c.COLUMNS)
                for r in pts:
                    w.writerow(
                        [
                            c.ION_Z[ion],
                            f"{mass:.6f}",
                            target,
                            f"{float(r['energy_kev_per_amu']) * mass * 1000.0:.7g}",
                            r["stopping_ev_cm2_per_1e15_atoms"],
                            r["u_stopping"],
                            f"digitized plot, {r['method']}",
                            f"Digitized from the srim.org compilation: {ref}",
                            key,
                            r["origin"],
                        ]
                    )
            print(f"{len(pts)} points -> {path}")


if __name__ == "__main__":
    main(sys.argv)
