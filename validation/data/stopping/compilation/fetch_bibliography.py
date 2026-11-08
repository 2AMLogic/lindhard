#!/usr/bin/env python3
"""Build `bibliography.csv` from the srim.org experimental-paper citation files.

srim.org lists, for every ion-target pair (`Citations/Cit <ZZ>-<ZZ>.pdf`,
linked from `PlotCit.htm`) and for 148 compounds
(`SRIM/Compounds/Citations/<compound>.pdf`, linked from
`SRIM/Compounds.htm`), the experimental stopping and range papers its compiler
knew. Each entry gives a citation key, year, authors, title, journal reference,
an old reference number and a comment on what was measured. Those are facts
about published papers; CONTRIBUTING.md ("The srim.org experimental-data
compilation") allows them in the tree as CC0. The PDFs themselves are not
committed: they are fetched into a scratch directory outside the tree.

Usage: fetch_bibliography.py [--scratch DIR] [--out-dir DIR] [--no-fetch]
                             [--no-doi] [--max-doi N]

  fetch (default)  download the index pages and citation PDFs into SCRATCH at
                   most once per second, with a browser User-Agent and a
                   Referer (srim.org answers 403 to bare scripted requests);
                   files already in SCRATCH are not fetched again
  --no-fetch       parse what is already in SCRATCH
  --no-doi         skip the Crossref lookup (the `doi` column stays empty
                   except for DOIs already cached)
  --max-doi N      look up at most N uncached DOIs this run (the cache in
                   SCRATCH makes the lookup resumable)

Outputs, in OUT_DIR: `bibliography.csv` (one row per paper) and
`bibliography-unparsed.txt` (every entry the parser could not read, verbatim,
so that nothing is dropped silently). The script exits non-zero when a
per-file entry count does not match the parsed plus unparsed entries.

The `Cit <Z>-All` and `Cit All-<Z>` aggregates repeat per-pair entries, so
they are never a source of rows; they are only a cross-check of totals.
`pdftotext -layout` (poppler) must be on PATH. The PDFs are untrusted data:
they are only ever passed to `pdftotext`, whose text output is parsed as plain
text.
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

SITE = "http://www.srim.org/"
REFERER = SITE + "PlotCit.htm"
USER_AGENT = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"
)
CROSSREF = "https://api.crossref.org/works"
MAILTO = "lindhard-bibliography@users.noreply.github.com"
FETCH_DELAY_S = 1.0

HERE = os.path.dirname(os.path.abspath(__file__))

FIELDS = [
    "citation_key", "year", "authors", "title", "journal", "comment",
    "old_ref", "pairs", "target", "doi",
]

# Page furniture repeated on every PDF page.
FURNITURE = re.compile(
    r"^\s*(\(C\) \d{4} SRIM\.org|Experimental Stopping and Range Papers|Page \d+ of \d+|"
    r"Stopping for Ion\b|Citations for Target\b|Pub\.\s+Citation|Year\s+Authors, Title)"
)
# First line of an entry: authors, then (two or more spaces) the citation key.
# Keys are printed as `YYYY-Name`, sometimes with a space or apostrophe.
ENTRY_START = re.compile(r"^\s+(\S.*?)\s+(\d{4}-[A-Za-z]\S*(?: \S+)?)\s*$")
YEAR_LINE = re.compile(r"^\s*(\d{4})\s+(['\"‘’“”].*)$")
TRAILING_REF = re.compile(r"^(.*?)\s{3,}(\d{4})\s*$")
LONE_REF = re.compile(r"^\s+(\d{4})\s*$")
COMMENT = re.compile(r"^\s*Comment\s*:\s*(.*)$")


# --- fetching ---------------------------------------------------------------

_last_request = [0.0]


def http_get(url, referer=REFERER, accept=None):
    wait = FETCH_DELAY_S - (time.monotonic() - _last_request[0])
    if wait > 0:
        time.sleep(wait)
    headers = {"User-Agent": USER_AGENT, "Referer": referer}
    if accept:
        headers["Accept"] = accept
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.read()
    finally:
        _last_request[0] = time.monotonic()


def fetch_cached(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = http_get(url)
    with open(path + ".part", "wb") as f:
        f.write(data)
    os.replace(path + ".part", path)


def index_links(html, pattern):
    """Distinct (href, name) pairs for hrefs matching PATTERN, in page order."""
    seen, out = set(), []
    for m in re.finditer(r"""href\s*=\s*["']?([^"'\s>]+)""", html, re.I):
        href = m.group(1)
        mm = re.search(pattern, href)
        if mm and href not in seen:
            seen.add(href)
            out.append((href, mm.group(1)))
    return out


def discover(scratch, fetch):
    """([(pair code, url, path)], [(aggregate code, url, path)], [(compound, url, path)])."""
    plotcit = os.path.join(scratch, "PlotCit.htm")
    compounds = os.path.join(scratch, "Compounds.htm")
    if fetch:
        fetch_cached(SITE + "PlotCit.htm", plotcit)
        fetch_cached(SITE + "SRIM/Compounds.htm", compounds)
    pairs, aggs, comps = [], [], []
    with open(plotcit, encoding="latin-1") as f:
        html = f.read()
    for href, name in index_links(html, r"^Citations/Cit%20([0-9A-Za-z]+-[0-9A-Za-z]+)\.pdf$"):
        url = SITE + href
        path = os.path.join(scratch, "pairs", f"Cit {name}.pdf")
        if re.fullmatch(r"\d\d-\d\d", name):
            pairs.append((name, url, path))
        else:
            aggs.append((name, url, path))
    with open(compounds, encoding="latin-1") as f:
        html = f.read()
    for href, name in index_links(html, r"Compounds/Citations/([^/]+)\.pdf$"):
        url = SITE + "SRIM/" + urllib.parse.quote(href, safe="%/")
        name = urllib.parse.unquote(name)
        comps.append((name, url, os.path.join(scratch, "compounds", f"{name}.pdf")))
    return pairs, aggs, comps


def pdf_text(path):
    txt = path[:-4] + ".txt"
    if not os.path.exists(txt) or os.path.getmtime(txt) < os.path.getmtime(path):
        subprocess.run(["pdftotext", "-layout", path, txt], check=True)
    with open(txt, encoding="utf-8", errors="replace") as f:
        return f.read()


# --- parsing ----------------------------------------------------------------

def clean_lines(text):
    return [l.rstrip() for l in text.replace("\f", "\n").split("\n") if not FURNITURE.match(l)]


def split_entries(text):
    """Entry blocks (lists of lines); the first line of each is its ENTRY_START line."""
    blocks, cur = [], None
    for line in clean_lines(text):
        if ENTRY_START.match(line) and not COMMENT.match(line):
            cur = [line]
            blocks.append(cur)
        elif cur is not None:
            cur.append(line)
    return blocks


def parse_block(block):
    """dict of fields, or None if the entry does not have the expected shape."""
    m = ENTRY_START.match(block[0])
    authors = " ".join(m.group(1).split())
    key = m.group(2)
    body = block[1:]
    # Authors may continue on lines before the year line.
    i = 0
    while i < len(body) and not YEAR_LINE.match(body[i]):
        if COMMENT.match(body[i]) or not body[i].strip():
            return None
        authors += " " + " ".join(body[i].split())
        i += 1
    if i == len(body):
        return None
    ym = YEAR_LINE.match(body[i])
    year, title_parts, ref = ym.group(1), [], ""
    line, i = ym.group(2), i + 1
    glued = None  # an empty title (`''`) followed on the same line by the journal
    gm = re.match(r"^(['\"\u2019\u201d]{2})(\S.*)$", line)
    if gm:
        line, glued = gm.group(1), gm.group(2)
    while not glued:
        tm = TRAILING_REF.match(line)
        if tm:
            line, ref = tm.group(1), tm.group(2)
        title_parts.append(" ".join(line.split()))
        if re.search(r"['\"’”]\s*$", line) or i >= len(body):
            break
        if COMMENT.match(body[i]):
            return None
        line, i = body[i], i + 1
    title = " ".join(title_parts).strip()
    if len(title) >= 2 and title[0] in "'\"‘“":
        title = title[1:]
    if title and title[-1] in "'\"’”":
        title = title[:-1]
    journal, comment, in_comment = [], [], False
    if glued:
        tm = TRAILING_REF.match(glued)
        if tm:
            glued, ref = tm.group(1), tm.group(2)
        journal.append(" ".join(glued.split()))
    for line in body[i:]:
        if not line.strip():
            if in_comment:
                break
            continue
        cm = COMMENT.match(line)
        if cm:
            in_comment = True
            comment.append(cm.group(1))
            continue
        lm = LONE_REF.match(line)
        if lm and not in_comment:
            ref = ref or lm.group(1)
            continue
        tm = TRAILING_REF.match(line)
        if tm and not in_comment:
            line, ref = tm.group(1), ref or tm.group(2)
        (comment if in_comment else journal).append(" ".join(line.split()))
    journal = " ".join(journal).strip()
    if not (key and year and journal):
        return None
    return {
        "citation_key": key,
        "year": year,
        "authors": authors,
        "title": title,
        "journal": journal,
        "comment": " ".join(comment).strip(),
        "old_ref": ref,
    }


def count_entries(text):
    """Entry count made independently of the field parser: `Comment :` lines
    (every entry has exactly one) must equal the number of entry-start lines."""
    lines = clean_lines(text)
    return (sum(1 for l in lines if COMMENT.match(l)),
            sum(1 for l in lines if ENTRY_START.match(l) and not COMMENT.match(l)))


def norm(s):
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


# --- Crossref ---------------------------------------------------------------

STOP = {"of", "the", "and", "a", "an", "for", "in"}


def journal_name(ref):
    """Journal part of a reference such as `Phys. Rev. B, 10, 3781-88 (1974)`."""
    return re.split(r",\s*(?:[A-Z]{1,3}-?\d|\d)", ref, maxsplit=1)[0].strip()


def first_page(ref):
    """First page of a reference such as `Phys. Rev., 84, 489-92 (1951)` ('' if none)."""
    m = re.findall(r",\s*[A-Za-z]?(\d+)\s*(?:-\s*\d+)?\s*\(\d{4}\)", ref)
    return m[-1] if m else ""


def journal_matches(cited, full_titles):
    """Every word of the printed journal abbreviation is the start of the
    corresponding word of the Crossref journal title (Phys. Rev. Letters ~
    Physical Review Letters)."""
    a = [w for w in norm(cited).split() if w not in STOP]
    if not a:
        return False
    for ft in full_titles:
        b = [w for w in norm(ft).split() if w not in STOP]
        if a == b or (len(a) == len(b) and all(y.startswith(x) for x, y in zip(a, b))):
            return True
    return False


def crossref_doi(row, cache, budget):
    """DOI if exactly one Crossref candidate matches title, year, journal and
    (where the citation prints one) first page; else ''."""
    ck = row["citation_key"] + "|" + norm(row["title"])
    if ck not in cache:
        if budget[0] <= 0:
            return ""
        budget[0] -= 1
        q = urllib.parse.urlencode({
            "query.bibliographic": row["title"],
            "query.author": row["authors"].split(",")[0],
            "rows": 5,
            "select": "DOI,title,issued,container-title,short-container-title,page",
            "mailto": MAILTO,
        })
        try:
            data = json.loads(http_get(f"{CROSSREF}?{q}", referer="https://www.crossref.org/"))
            cache[ck] = data["message"]["items"]
        except (urllib.error.URLError, ValueError, KeyError, TimeoutError) as e:
            print(f"crossref: {row['citation_key']}: {e}", file=sys.stderr)
            return ""
    want_title, want_year = norm(row["title"]), int(row["year"])
    jname, page = journal_name(row["journal"]), first_page(row["journal"])
    hits = set()
    for it in cache[ck]:
        titles = it.get("title") or []
        issued = (it.get("issued") or {}).get("date-parts") or [[None]]
        containers = (it.get("container-title") or []) + (it.get("short-container-title") or [])
        if (any(norm(t) == want_title for t in titles) and issued[0][0] == want_year
                and journal_matches(jname, containers)
                and (not page or re.sub(r"\D.*", "", re.sub(r"^\D+", "", it.get("page", ""))) == page)):
            hits.add(it["DOI"].lower())
    return hits.pop() if len(hits) == 1 else ""


# --- main -------------------------------------------------------------------

def merge(rows, parsed, pair=None, target=None):
    k = (parsed["citation_key"], norm(parsed["title"]))
    r = rows.setdefault(k, dict(parsed, pairs=[], target=[], doi=""))
    if pair and pair not in r["pairs"]:
        r["pairs"].append(pair)
    if target and target not in r["target"]:
        r["target"].append(target)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scratch", default=os.path.join(tempfile.gettempdir(), "lindhard-srim-bibliography"))
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--no-fetch", action="store_true")
    ap.add_argument("--no-doi", action="store_true")
    ap.add_argument("--max-doi", type=int, default=10**9)
    args = ap.parse_args()
    fetch = not args.no_fetch

    pairs, aggs, comps = discover(args.scratch, fetch)
    print(f"{len(pairs)} pair files, {len(aggs)} aggregates, {len(comps)} compound files")
    if fetch:
        for n, (_, url, path) in enumerate(pairs + comps + aggs, 1):
            try:
                fetch_cached(url, path)
            except (urllib.error.URLError, TimeoutError) as e:
                if (_, url, path) in aggs:  # only a cross-check; a dead link is not fatal
                    print(f"warning: aggregate not fetched, {url}: {e}", file=sys.stderr)
                    continue
                sys.exit(f"fetch failed for {url}: {e}")
            if n % 50 == 0:
                print(f"fetched {n} files", flush=True)

    rows, unparsed, bad, per_ion = {}, [], 0, {}
    sources = [("pair", c, p) for c, _, p in pairs] + [("target", c, p) for c, _, p in comps]
    for kind, code, path in sources:
        text = pdf_text(path)
        n_comment, n_start = count_entries(text)
        ok = 0
        for block in split_entries(text):
            parsed = parse_block(block)
            if parsed is None:
                unparsed.append((kind, code, block))
                continue
            ok += 1
            merge(rows, parsed, pair=code if kind == "pair" else None,
                  target=code if kind == "target" else None)
        if not (n_comment == n_start == ok + sum(1 for k, c, _ in unparsed if (k, c) == (kind, code))):
            bad += 1
            print(f"COUNT MISMATCH {kind} {code}: Comment lines {n_comment}, "
                  f"entry starts {n_start}, parsed {ok}", file=sys.stderr)

    # Cross-check of the Cit <Z>-All aggregates (one line per distinct paper
    # for ion Z) against the distinct papers found in the per-pair files.
    for k, r in rows.items():
        for pr in r["pairs"]:
            per_ion.setdefault(pr.split("-")[0], set()).add(k)
    for code, _, path in aggs:
        if not os.path.exists(path):
            continue
        z, other = code.split("-")
        if z != "All" and z in per_ion:
            _, n = count_entries(pdf_text(path))
            if n < len(per_ion[z]):  # compound targets add entries beyond the pairs
                print(f"aggregate {code}: {n} entries, {len(per_ion[z])} distinct papers in the pair files", file=sys.stderr)

    cache_path = os.path.join(args.scratch, "crossref-cache.json")
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    budget = [0 if args.no_doi else args.max_doi]
    try:
        for r in rows.values():
            r["doi"] = crossref_doi(r, cache, budget)
    finally:
        os.makedirs(args.scratch, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(cache, f)

    out = sorted(rows.values(), key=lambda r: (r["year"], r["citation_key"], norm(r["title"])))
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "bibliography.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, FIELDS, lineterminator="\n")
        w.writeheader()
        for r in out:
            w.writerow(dict(r, pairs=";".join(sorted(r["pairs"])), target=";".join(sorted(r["target"]))))
    with open(os.path.join(args.out_dir, "bibliography-unparsed.txt"), "w", encoding="utf-8") as f:
        f.write("# Entries the parser could not read, verbatim from `pdftotext -layout`.\n")
        for kind, code, block in unparsed:
            f.write(f"\n=== {kind} {code} ===\n" + "\n".join(block).rstrip() + "\n")
    print(f"{len(out)} papers written, {len(unparsed)} entries unparsed, "
          f"{sum(1 for r in out if r['doi'])} with a DOI, {bad} files with a count mismatch")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
