#!/usr/bin/env python3
"""Transcribe the measured electron backscatter coefficients of C, Al, Si, Cu
and Au from D. C. Joy's database into one JSON file per original measurement
set, and cross-check every set against a second, independent compilation (#148).

Primary source (read, not recalled): D. C. Joy, "A Database of
Electron-Solid Interactions", revision 01-01 (2001), Microsoft Word file
`database.doc` (488 960 bytes, SHA-256 9445f8dd...afd25fd7), archived by the
Internet Archive on 2004-03-20 at
http://web.archive.org/web/20040320063526/http://web.utk.edu:80/~srcutk/database.doc
Section 2(a), "Backscattered Electron Yield Data ... for elements": per
element, numbered data sets ("Data set #n Reference (k)", the first set
unnumbered) of (E in keV, BS yield), and the numbered "Reference List" at the
front of the file. The published description of the database is D. C. Joy,
Scanning 17, 270 (1995), doi:10.1002/sca.4950170501 (not opened).

Cross-check (second read): F. Akbari, "A comprehensive open-access database
of electron backscattering coefficients for energies ranging from 0.1 keV to
15 MeV", Zenodo, doi:10.5281/zenodo.7810951 (CC-BY 4.0), archive
"Energy dependence (csv).rar" (17 143 bytes, SHA-256 e35d1155...1d7846987),
one CSV per target (eta in %). Its sets are numbered per table by its own
reference list, which could not be opened, so it supplies no attribution; it
is used only to re-read Joy's values (each Joy set is matched to the Akbari set
sharing most points) and its values are not stored as data.

Usage (fetch both files outside the tree, unpack the RAR with any tool):
    ingest_joy2001.py database.doc "<dir with C.csv Al.csv Si.csv Si-a.csv Cu.csv Au.csv>" OUTDIR
Writes OUTDIR/eta_<el>_<author><year>.json and OUTDIR/crosscheck_akbari2022.json
and prints the provenance rows (markdown) on stdout. Standard library only.
The inputs are untrusted data: parsed as text, never executed.
"""
import csv
import hashlib
import json
import re
import sys
from decimal import Decimal
from pathlib import Path

ADDED = "2026-10-07"
DOC_SHA256 = "9445f8ddbebb72892a30e2e3bb5c77dac3b4d5b36790287e7c0ce658afd25fd7"
DOC_URL = (
    "http://web.archive.org/web/20040320063526/http://web.utk.edu:80/~srcutk/database.doc"
)
COMPILATION = (
    "D. C. Joy, 'A Database of Electron-Solid Interactions', revision 01-01 (2001), "
    "database.doc (SHA-256 " + DOC_SHA256 + "), Section 2(a) 'Backscattered Electron "
    "Yield Data, for elements'; published description: D. C. Joy, Scanning 17, 270 "
    "(1995), doi:10.1002/sca.4950170501 (not opened)"
)
AKBARI = (
    "F. Akbari, Zenodo doi:10.5281/zenodo.7810951 (CC-BY 4.0), 'Energy dependence "
    "(csv).rar' (SHA-256 e35d11553720b6d122e129928993448e7f6456f189a21ad640ec15b1d7846987)"
)

# Element headings of Section 2(a), as printed, and the Akbari CSV of each.
ELEMENTS = {
    "Carbon (Z=6)": ("C", 6, "C.csv"),
    "Aluminum (Z=13)": ("Al", 13, "Al.csv"),
    "Silicon (Z=14)": ("Si", 14, "Si.csv"),
    "Silicon (Z=14) Amorphous": ("Si", 14, "Si-a.csv"),
    "Copper (Z=29)": ("Cu", 29, "Cu.csv"),
    "Gold (Z=79)": ("Au", 79, "Au.csv"),
}

# Joy's reference numbers used by these elements: (id stem, reference as
# printed in Joy's Reference List, DOI checked against Crossref on 2026-10-07
# or None, note). A DOI is given only where Crossref's volume, first page,
# year and authors agree with Joy's entry.
REFS = {
    1: ("hunger1979", "H.-J. Hunger and L. Kuchler, phys. stat. sol. (a) 56, K45 (1979)",
        "10.1002/pssa.2210560157",
        "Crossref title: 'Measurements of the electron backscattering coefficient for quantitative EPMA in the energy range of 4 to 40 keV'"),
    2: ("reimer1980", "L. Reimer and C. Tolkamp, Scanning 3, 35 (1980)",
        "10.1002/sca.4950030105",
        "Crossref spells the second author Tollkamp; title 'Measuring the backscattering coefficient and secondary electron yield inside a scanning electron microscope'"),
    3: ("moncrieff1976", "D. A. Moncrieff and P. R. Barker, Scanning 1, 195 (1976)",
        "10.1002/sca.4950010307",
        "Crossref gives the year as 1978 (Scanning 1, 195-197, 'Secondary electron emission in the scanning electron microscope'); volume and page agree, Joy prints 1976"),
    4: ("bongeler1993", "R. Bongeler, U. Golla, M. Kussens, L. Reimer, B. Schendler, R. Senkel and M. Spranck, Scanning 15, 1 (1993)",
        "10.1002/sca.4950150102",
        "Crossref: Boengeler, Golla, Kaessens, Reimer, Schindler, Senkel, Spranck, 'Electron-specimen interactions in low-voltage scanning electron microscopy'"),
    5: ("shimizu1974", "R. Shimizu, J. Appl. Phys. 45, 2107 (1974)",
        "10.1063/1.1663552",
        "Crossref title: 'Secondary electron yield with primary electron beam of kilo-electron-volts'"),
    6: ("bishop1963", "H. E. Bishop, PhD thesis, University of Cambridge (1963)", None,
        "thesis, no DOI; Joy's introduction cites 'Bishop H, (1966), Ph.D Thesis University of Cambridge'"),
    7: ("philibert1963", "J. Philibert and E. Weinryb, Proc. 3rd Conf. on X-ray Optics and Microanalysis, p. 163 (1963); also E. Weinryb and J. Philibert, C. R. Acad. Sci. 256, 4535 (1964)", None,
        "conference proceedings, no DOI found"),
    8: ("heinrich1966", "K. F. J. Heinrich, Proc. 4th Conf. on X-ray Optics and Microanalysis, ed. R. Castaing et al. (Hermann, Paris), p. 159 (1966)", None,
        "conference proceedings, no DOI found"),
    13: ("neubert1980", "G. Neubert and S. Rogaschewski, phys. stat. sol. (a) 59, 35 (1980)",
         "10.1002/pssa.2210590104",
         "Crossref title: 'Backscattering coefficient measurements of 15 to 60 keV electrons for solids at various angles of incidence'"),
    14: ("kanter1961", "M. Kanter, Phys. Rev. 121, 1677 (1961)", None,
         "page as printed by Joy; Crossref has no Phys. Rev. 121 article at p. 1677. Kanter's two 1961 articles in that volume are pp. 677-681 (doi:10.1103/physrev.121.677) and 681-684 (doi:10.1103/physrev.121.681); which one holds these values is not settled, so no DOI is given"),
    15: ("drescher1970", "H. Drescher, L. Reimer and M. Seidel, Z. angew. Physik 29, 331 (1970)", None,
         "no Crossref record found"),
    22: ("cosslett1965", "V. E. Cosslett and R. N. Thomas, Brit. J. Appl. Phys. 16, 774 (1965)", None,
         "page as printed by Joy; Crossref has Cosslett and Thomas, Brit. J. Appl. Phys. 16, 779-796 (1965), 'Multiple scattering of 5-30 keV electrons in evaporated metal films III: Backscattering and absorption', doi:10.1088/0508-3443/16/6/303, probably this paper (first page differs), so no DOI is asserted"),
    23: ("koshikawa1973", "T. Koshikawa and R. Shimizu, J. Phys. D: Appl. Phys. 6, 1369 (1973)",
         "10.1088/0022-3727/6/11/312",
         "Crossref title: 'Secondary electron and backscattering measurements for polycrystalline copper with a spherical retarding-field analyser'"),
    34: ("sternglass1954", "E. J. Sternglass, Phys. Rev. 95, 345 (1954)",
         "10.1103/physrev.95.345",
         "Crossref title: 'Backscattering of Kilovolt Electrons from Solids'"),
    35: ("palluel1947", "P. Palluel, Compt. Rend. 224, 1492 (1947)", None,
         "no Crossref record found"),
    68: ("wittry1966", "D. B. Wittry, Proc. 4th Conf. on X-ray Optics and Microanalysis, ed. R. Castaing et al. (Hermann, Paris), p. 168 (1966)", None,
         "conference proceedings, no DOI found"),
    86: ("joy1996", "D. C. Joy and C. S. Joy, SEMATECH Report TT# 96063130A-TR (August 1996)", None,
         "technical report, no DOI"),
    106: ("bronstein1969", "I. M. Bronstein and B. S. Fraiman, Vtorichnaya Elektronnaya Emissiya (Nauka, Moskva, 1969), p. 340", None,
          "book, no DOI"),
    107: ("elgomati1997", "M. M. El Gomati and A. M. D. Assad, Proc. 5th European Workshop on Microbeam Analysis (1997), in press when compiled", None,
          "proceedings 'in press' in the compilation, no DOI found"),
}


# Short attribution for the provenance rows.
SHORT = {
    1: "Hunger and Kuchler (1979)", 2: "Reimer and Tolkamp (1980)", 3: "Moncrieff and Barker (1976)",
    4: "Bongeler et al. (1993)", 5: "Shimizu (1974)", 6: "Bishop (1963)", 7: "Philibert and Weinryb (1963)",
    8: "Heinrich (1966)", 13: "Neubert and Rogaschewski (1980)", 14: "Kanter (1961)",
    15: "Drescher, Reimer and Seidel (1970)", 22: "Cosslett and Thomas (1965)",
    23: "Koshikawa and Shimizu (1973)", 34: "Sternglass (1954)", 35: "Palluel (1947)",
    68: "Wittry (1966)", 86: "Joy and Joy (1996)", 106: "Bronstein and Fraiman (1969)",
    107: "El Gomati and Assad (1997)",
}


def text_runs(doc: bytes) -> list[str]:
    """Runs of printable ASCII and tab, i.e. what `strings -n 1` prints: the
    Word file stores this text uncompressed, one paragraph per run."""
    return [m.group().decode("ascii") for m in re.finditer(rb"[\x20-\x7e\t]+", doc)]


def joy_sets(lines: list[str]) -> list[dict]:
    start = lines.index("Backscattered Electron Yield Data")
    end = next(i for i in range(start, len(lines)) if re.match(r"\s*SECTION THREE", lines[i], re.I))
    heading = re.compile(r"^\s*[A-Z][a-z]+.*\(Z\s*=\s*\d+\)")
    out, cur, ds = [], None, None
    for i in range(start, end):
        s = lines[i].strip()
        if heading.match(lines[i]):
            cur = s if s in ELEMENTS else None
            ds = None
            continue
        if cur is None:
            continue
        m = re.match(r"^(?:Data set\s*#\s*(\d+)\s*)?Reference\s*\((\d+)\)$", s)
        if m:
            ds = {"heading": cur, "set": int(m.group(1) or 1), "ref": int(m.group(2)),
                  "header_text": s, "points": []}
            out.append(ds)
            continue
        if not s or s.startswith("E(keV)"):
            continue
        m = re.match(r"^([0-9.]+)\s+([0-9.]+)$", s)
        if m and ds is not None:
            ds["points"].append((m.group(1), m.group(2)))
            continue
        sys.exit(f"unparsed line in {cur}: {lines[i]!r}")
    return out


def akbari_sets(path: Path) -> list[dict]:
    rows = list(csv.reader(path.read_text(encoding="utf-8-sig").splitlines()))
    out, cur = [], None
    for r in rows[3:]:
        r = [c.strip() for c in r] + [""] * 5
        if not r[0]:
            continue
        if r[2]:
            cur = {"ref": r[2], "comment": " ".join(c for c in r[3:5] if c), "rows": []}
            out.append(cur)
        cur["rows"].append((Decimal(r[0]), Decimal(r[1])))
    return out


def half_ulp(printed: str) -> float:
    """Half a unit of the last printed digit: the rounding of a table entry."""
    decimals = len(printed.split(".")[1]) if "." in printed else 0
    return float(Decimal(5) / Decimal(10) ** (decimals + 1))


def main(doc_path: str, akbari_dir: str, out_dir: str) -> None:
    doc = Path(doc_path).read_bytes()
    if hashlib.sha256(doc).hexdigest() != DOC_SHA256:
        sys.exit("database.doc is not the archived revision 01-01 file (SHA-256 differs)")
    sets = joy_sets(text_runs(doc))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cross, prov = [], []
    for s in sets:
        el, z, akbari_csv = ELEMENTS[s["heading"]]
        amorphous = "Amorphous" in s["heading"]
        stem, ref_text, doi, ref_note = REFS[s["ref"]]
        sid = f"eta_{el.lower()}_{stem}" + ("_amorphous" if amorphous else "")
        # Second read: the Akbari set sharing most (E, eta) points.
        ak = akbari_sets(Path(akbari_dir) / akbari_csv)
        best = max(
            range(len(ak)),
            key=lambda j: sum(
                1 for e, v in s["points"]
                if any(ae == Decimal(e) and abs(av / 100 - Decimal(v)) < Decimal("0.0006")
                       for ae, av in ak[j]["rows"])
            ),
        )
        amap = {ae: av / 100 for ae, av in ak[best]["rows"]}
        compared = []
        for e, v in s["points"]:
            av = amap.get(Decimal(e))
            compared.append({
                "energy_kev": e, "eta_joy": v,
                "eta_akbari": None if av is None else str(av.normalize()),
                "agree": av is not None and abs(av - Decimal(v)) <= Decimal("0.0005"),
            })
        n_agree = sum(c["agree"] for c in compared)
        disagree = [c for c in compared if not c["agree"]]
        agreement = (
            f"{n_agree} of {len(compared)} points re-read identically (to the 0.1 % rounding of "
            f"the Akbari table) in Akbari set {ak[best]['ref']}"
            + (f" ('{ak[best]['comment']}')" if ak[best]["comment"] else "")
            + f" of {akbari_csv}"
        )
        if disagree:
            agreement += "; differing: " + "; ".join(
                f"{c['energy_kev']} keV Joy {c['eta_joy']} vs Akbari "
                + (c["eta_akbari"] if c["eta_akbari"] is not None else "no point at this energy")
                for c in disagree
            ) + " (Joy's value is stored)"
        cross.append({"id": sid, "joy_set": s["header_text"], "akbari_csv": akbari_csv,
                      "akbari_set": ak[best]["ref"], "akbari_comment": ak[best]["comment"],
                      "points": compared})
        points = []
        for e, v in s["points"]:
            points.append({
                "energy_ev": float(Decimal(e) * 1000),
                "eta": float(Decimal(v)),
                "eta_unc_abs": half_ulp(v),
            })
        label = f"{s['heading']}, " + ("first set, 'Reference (%d)'" % s["ref"] if s["set"] == 1
                                       else f"'Data set #{s['set']} Reference ({s['ref']})'")
        d = {
            "kind": "backscatter_coefficient",
            "id": sid,
            "target": el,
            "z": z,
            "target_state": (
                "amorphous Si, as headed by the compilation ('Silicon (Z=14) Amorphous')"
                if amorphous else
                "not stated by the compilation (no sample preparation, surface condition or thickness is given per set)"
            ),
            "incidence_deg": 0.0,
            "incidence_note": "not stated by the compilation, which tabulates the yield against energy only; taken as normal incidence on a bulk target",
            "eta_threshold_note": "the compilation does not state the energy threshold that separates backscattered from secondary electrons per set; the conventional 50 eV is assumed",
            "points": points,
            "original_reference": ref_text,
            "original_doi": doi,
            "original_reference_note": ref_note + "; the original was not opened",
            "compilation": COMPILATION,
            "compilation_set": label,
            "compilation_reference_number": s["ref"],
            "url": DOC_URL,
            "crosscheck": {"compilation": AKBARI, "agreement": agreement,
                           "record": "validation/data/backscatter/crosscheck_akbari2022.json"},
            "extraction": (
                "transcribed from the compilation's table (text of the Word file, extracted by "
                "validation/data/backscatter/ingest_joy2001.py); E in keV times 1000, eta as printed. "
                "The compilation gives no measurement uncertainty; eta_unc_abs is half a unit of the "
                "last printed digit (the table's rounding), not an experimental error"
            ),
            "terms": "Facts, cited (operator ruling #34: published measurements may be committed as cited reference facts)",
            "added": ADDED,
        }
        (out / f"{sid}.json").write_text(json.dumps(d, indent=2, ensure_ascii=False) + "\n")
        es = [float(Decimal(e)) for e, _ in s["points"]]
        prov.append(
            f"| Measured backscatter coefficient of {el}"
            + (" (amorphous)" if amorphous else "")
            + f", set of {SHORT[s['ref']]} (`{sid}`) | "
            f"`validation/data/backscatter/{sid}.json` | Experimental measurement, transcribed "
            f"from a table of a secondary compilation | Measurement: {ref_text}"
            + (f", doi:{doi}" if doi else "") + f"; **not seen** ({ref_note}). Read from: Joy "
            f"database rev. 01-01 (2001), Section 2(a), {label}, {len(points)} points, "
            f"{min(es):g} to {max(es):g} keV (Wayback copy, SHA-256 in the file). Second read: "
            f"Akbari (2022) Zenodo table, {n_agree} of {len(points)} points identical "
            f"(`crosscheck_akbari2022.json`) | Facts, cited (#34) | {ADDED} |"
        )
    (out / "crosscheck_akbari2022.json").write_text(json.dumps({
        "format": "lindhard-backscatter-crosscheck/1",
        "primary": COMPILATION + "; " + DOC_URL,
        "second_read": AKBARI,
        "method": (
            "Each Joy set is matched to the Akbari set (same element table) that shares most "
            "(E, eta) points; every Joy point is then looked up at the same printed energy. "
            "agree: |eta_akbari - eta_joy| <= 0.0005 (Akbari prints eta in % to 0.1). "
            "Joy's values are the stored ones; Akbari's are recorded here only as the re-read"
        ),
        "sets": cross,
    }, indent=1) + "\n")
    print("\n".join(prov))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
