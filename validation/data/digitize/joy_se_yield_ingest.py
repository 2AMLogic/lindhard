#!/usr/bin/env python3
"""Transcribe the measured secondary-electron yields of Al, Cu, Si and Au from
D. C. Joy, "A database on electron-solid interactions" (Revision 01-01, 2001;
cf. Scanning 17, 270 (1995), doi:10.1002/sca.4950170501), Section 1(a), into
one JSON file per *original measurement set* (issue #149).

Stdlib only. The compilation is a Word 97 document that Joy distributed at
<http://web.utk.edu/~srcutk/database.doc>. That host no longer answers; the
copy read here is the Internet Archive capture of 2010:

    curl -L -o joy_database.doc \\
        "https://web.archive.org/web/2010id_/http://web.utk.edu/~srcutk/database.doc"
    validation/data/digitize/joy_se_yield_ingest.py joy_database.doc \\
        --out validation/data/se_yield [--check]

The document is not stored in this tree. This is a *transcription of tables*,
not a digitization: the text runs of the file are scanned for the printable
lines Word stores (no Word library), the "Aluminum Z=13", "Silicon Z=14",
"Silicon (Z=14) amorphous", "Copper Z=29" and "Gold (Z=79)" blocks of
Section 1(a) are cut out (each ends at the next element's heading), every
"Reference (n)" header opens a data set, and every "E(keV) <tab> SE yield"
line pair is read as written. Reference numbers are resolved through the
document's own Reference List, parsed by the same script. Values are
kept exactly as printed (including duplicate energies and the printed
precision); nothing is converted except keV to eV (x 1000, exact in decimal).

`--check` recounts: it re-reads the document a second time with a different
line-splitting rule (the raw bytes split on every non-printable byte, rather
than on the paragraph mark) and fails if the two reads differ, and fails if a
block has no data set or a data set has no points. It does not replace the
by-eye comparison of the committed JSON with the document that the PR records.

Compilation facts that matter for use (Joy, Introduction, "Secondary Yields"):
the energy range counted as secondary is "usually now taken to be 0 to 50 eV,
although in some early work 0 to 70 or even 0 to 100 eV was used"; Reference 20
(Dawson 1966) is flagged by the compilation as taken for cutoffs from 30 to
100 eV. The compilation does not state the angle of incidence of any set, and
this script does not assume one.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

# (heading regex, next-heading regex, element, short id, target_state text)
BLOCKS = [
    ("Aluminum Z=13", "Silicon Z=14", "Al", "al", "as compiled; surface state not stated by the compilation"),
    ("Silicon Z=14", "Silicon (Z=14)  amorphous", "Si", "si", "as compiled; crystalline or not, and surface state, not stated by the compilation"),
    ("Silicon (Z=14)  amorphous", "Phosphorus Z=15", "Si", "si_amorphous", "amorphous silicon (the compilation's heading 'Silicon (Z=14) amorphous'); surface state not stated"),
    ("Copper Z=29", "Zinc (Z=30)", "Cu", "cu", "as compiled; surface state not stated by the compilation"),
    ("Gold (Z=79)", "Mercury (Z=80)", "Au", "au", "as compiled; surface state not stated by the compilation"),
]

NUM = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s+([0-9]*\.?[0-9]+)\s*$")
REF = re.compile(r"Reference\s*\(?\s*(\d+)\s*\)?")


def lines_by_paragraph(raw: bytes) -> list[str]:
    """Printable runs, split on the paragraph mark (read 1)."""
    out = []
    for chunk in re.split(rb"[\r\x07]", raw):
        for run in re.findall(rb"[\x09\x20-\x7e]+", chunk):
            out.append(run.decode("ascii"))
    return out


def lines_by_nonprintable(raw: bytes) -> list[str]:
    """Printable runs split on every other byte (read 2)."""
    return [m.decode("ascii") for m in re.findall(rb"[\x09\x20-\x7e]{1,}", raw)]


def section1_start(lines: list[str]) -> int:
    for i, ln in enumerate(lines):
        if ln.strip() == "Secondary Electron Yield Data":
            return i
    raise SystemExit("error: 'Secondary Electron Yield Data' heading not found")


def parse_references(lines: list[str]) -> dict[int, str]:
    refs: dict[int, str] = {}
    start = next(i for i, ln in enumerate(lines) if ln.strip() == "Reference List")
    end = section1_start(lines)
    # The document prints entry 42 (Glupe and Melhorne) without its number, so
    # that text is appended to entry 41; no data set used here cites 41 or 42.
    cur = None
    for ln in lines[start + 1 : end]:
        m = re.match(r"^\s*(\d+)\.?\s+(\S.*)$", ln)
        if m and 1 <= int(m.group(1)) <= 125 and (cur is None or cur < int(m.group(1)) <= cur + 2):
            cur = int(m.group(1))
            refs[cur] = m.group(2).strip()
        elif cur is not None and ln.strip() and not ln.startswith("Microanalysis"):
            refs[cur] += " " + ln.strip()
    return refs


def cut(lines: list[str], head: str, nxt: str, begin: int) -> list[str]:
    def norm(s: str) -> str:
        return re.sub(r"\s+", " ", s.strip())

    i = next(k for k in range(begin, len(lines)) if norm(lines[k]) == norm(head))
    j = next(k for k in range(i + 1, len(lines)) if norm(lines[k]).startswith(norm(nxt)))
    return lines[i + 1 : j]


def parse_block(block: list[str]) -> list[dict]:
    sets = []
    cur = None
    for ln in block:
        m = REF.search(ln)
        if m and "E(keV)" not in ln:
            cur = {"reference": int(m.group(1)), "points": []}
            sets.append(cur)
            continue
        m = REF.search(ln)
        if m and cur is None:
            cur = {"reference": int(m.group(1)), "points": []}
            sets.append(cur)
            continue
        n = NUM.match(ln)
        if n and cur is not None:
            cur["points"].append((n.group(1), n.group(2)))
    return sets


def ingest(raw: bytes, lines_fn) -> tuple[dict[int, str], list[tuple]]:
    lines = lines_fn(raw)
    refs = parse_references(lines)
    begin = section1_start(lines)
    out = []
    for head, nxt, el, tag, state in BLOCKS:
        sets = parse_block(cut(lines, head, nxt, begin))
        out.append((el, tag, state, sets))
    return refs, out


def short_ref(text: str) -> tuple[str, str]:
    """(first author's surname lowercased ascii, year) from a reference-list entry."""
    first = re.match(r"([A-Za-z\- ]+?)(?: [A-Z](?: [A-Z])?[,.]| [A-Z]{1,3}\b)", text)
    surname = (first.group(1) if first else text.split()[0]).split()[-1].strip(",.").lower()
    y = re.search(r"\((\d{4})\)|,\s*(\d{4})\b|\b(19\d\d)\b", text)
    year = (y.group(1) or y.group(2) or y.group(3)) if y else "nd"
    return surname, year


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("doc", type=Path)
    ap.add_argument("--out", type=Path, default=None, help="directory to write the JSON files to")
    ap.add_argument("--check", action="store_true", help="re-read with a second line rule and compare")
    args = ap.parse_args()
    raw = args.doc.read_bytes()
    refs, blocks = ingest(raw, lines_by_paragraph)
    if args.check:
        refs2, blocks2 = ingest(raw, lines_by_nonprintable)
        # read 2 may split a line in two; compare the numeric content only.
        if [(e, t, [(s["reference"], s["points"]) for s in sets]) for e, t, _, sets in blocks] != [
            (e, t, [(s["reference"], s["points"]) for s in sets]) for e, t, _, sets in blocks2
        ]:
            print("error: the two reads of the document differ", file=sys.stderr)
            return 1
    n_files = 0
    for el, tag, state, sets in blocks:
        if not sets:
            print(f"error: no data set for {el}", file=sys.stderr)
            return 1
        seen: dict[str, int] = {}
        for k, s in enumerate(sets, start=1):
            if not s["points"]:
                print(f"error: {el} data set #{k} (ref {s['reference']}) has no points", file=sys.stderr)
                return 1
            r = refs.get(s["reference"])
            if r is None:
                print(f"error: reference {s['reference']} not in the reference list", file=sys.stderr)
                return 1
            surname, year = short_ref(r)
            base = f"se_{tag}_{surname}{year}"
            seen[base] = seen.get(base, 0) + 1
            ident = base if seen[base] == 1 else f"{base}_{seen[base]}"
            rec = {
                "kind": "se_yield",
                "id": ident,
                "target": el,
                "target_state": state,
                "incidence_deg": None,
                "incidence_note": "not stated by the compilation for any set; normal incidence is NOT confirmed by the source",
                "points": [
                    {"energy_ev": round(float(e) * 1000.0, 6), "energy_kev_printed": e, "yield": float(y), "yield_printed": y}
                    for e, y in s["points"]
                ],
                "original_reference": r,
                "compilation_reference_number": s["reference"],
                "compilation": "D. C. Joy, A database on electron-solid interactions, Revision 01-01 (2001); cf. Scanning 17, 270 (1995), doi:10.1002/sca.4950170501",
                "compilation_location": f"Section 1(a), '{el}' heading, data set #{k} of {len(sets)} (document headings as printed)",
                "url": "https://web.archive.org/web/2010id_/http://web.utk.edu/~srcutk/database.doc",
                "extraction": "transcribed from the compilation's two-column table (E in keV, SE yield) by validation/data/digitize/joy_se_yield_ingest.py; values as printed, keV converted to eV; no curve fitted, nothing interpolated",
                "se_energy_range": "not stated per set; compilation: usually 0 to 50 eV, early work 0 to 70 or 100 eV; reference 20 flagged 30 to 100 eV",
                "terms": "Facts, cited (operator ruling #34). The compilation asks users to acknowledge it with its URL (http://web.utk.edu/~srcutk/); the compilation itself is not stored here.",
                "added": "2026-10-07",
            }
            if args.out:
                args.out.mkdir(parents=True, exist_ok=True)
                (args.out / f"{ident}.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n")
            n_files += 1
            print(f"{ident}: {len(s['points'])} points, ref {s['reference']}")
    print(f"{n_files} data sets", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
