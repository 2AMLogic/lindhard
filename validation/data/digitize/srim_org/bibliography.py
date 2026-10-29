"""Write the bibliography for the keys used in the compilation CSV.

Usage: python3 -I bibliography.py <mirror-site-dir> <compilation.csv> <out.csv>

Bibliographic facts only (key, citation number, authors, year, title,
journal reference), parsed from the per-system citation files
Citations/Cit zz-ZZ.pdf. The files' comment lines are not copied.
"""
import csv, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import citpdf


def split(entry):
    """(year, title, journal) from an entry's text lines."""
    body = [l for l in entry["lines"] if not l.startswith("Comment")]
    text = " ".join(body)
    m = re.match(r"(\d{4})\s+(.*)", text)
    year, rest = (m.group(1), m.group(2)) if m else ("", text)
    t = re.match(r"\s*'(.*?)'\s*(.*)", rest)
    title, journal = (t.group(1), t.group(2)) if t else ("", rest)
    clean = lambda s: re.sub(r"\s+", " ", s.replace("þ", "\"")).strip()
    return year, clean(title), clean(journal)


def main(site, comp, out):
    want = {}  # key -> set of (ion, target) pairs it appears on
    for row in csv.DictReader(open(comp)):
        for k in row["citation_keys"].split("|"):
            want.setdefault(k, set()).add(row["plot"][4:8])
    found = {}
    for e in citpdf.index(site).values():
        if e["key"] in want and e["key"] not in found:
            found[e["key"]] = e
    w = csv.writer(open(out, "w", newline=""), lineterminator="\n")
    w.writerow(["key", "old_ref", "authors", "year", "title", "journal"])
    for k in sorted(found):
        e = found[k]
        authors = "; ".join(a.strip() for a in re.split(r"\s{3,}", e["authors"]) if a.strip())
        w.writerow([k, e["ref"] or "", authors, *split(e)])
    missing = sorted(k for k in want if k not in found)
    print(f"{len(found)} entries; {len(missing)} unresolved old citation numbers (kept as ref<N>): {' '.join(missing)}", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:4])
