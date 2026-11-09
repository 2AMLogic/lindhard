"""Parse a srim.org Citations/Cit zz-ZZ.pdf into entries: citation key,
authors, old citation number ("ref", the number the plot legends print) and
the entry's text lines."""
import re, subprocess, sys, json

KEY = re.compile(r"\s(\d{4}-[A-Za-z][A-Za-z0-9]*)\s*$")
REF = re.compile(r"^\s+(\d{3,5})\s*$")
TRAILING_REF = re.compile(r"\s{8,}(\d{4})\s*$")


def index(site):
    """ref -> entry over every per-system citation file. Ref numbers are global
    to the compilation, so a legend row whose ref is missing from its own
    system's file is usually listed in another."""
    import glob
    out = {}
    for pdf in sorted(glob.glob(f"{site}/Citations/Cit *.pdf")):
        for e in parse(pdf):
            if e["ref"] and e["ref"] not in out:
                out[e["ref"]] = e
    return out


def parse(pdf):
    text = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, check=True).stdout.decode("latin-1")
    entries, cur = [], None
    for line in text.splitlines():
        m = KEY.search(line)
        if m and not line.lstrip().startswith(("Comment", "Pub.")):
            cur = {"key": m.group(1), "authors": line[:m.start()].strip(), "ref": None, "lines": []}
            entries.append(cur)
            continue
        if cur is None:
            continue
        m = REF.match(line)
        if m and cur["ref"] is None:
            cur["ref"] = int(m.group(1))
            continue
        # The number is sometimes printed far right on a text line instead.
        t = TRAILING_REF.search(line)
        if t and cur["ref"] is None and not line.lstrip().startswith("Comment"):
            cur["ref"] = int(t.group(1))
            line = line[:t.start()]
        if line.strip():
            cur["lines"].append(line.strip())
    return entries


if __name__ == "__main__":
    e = parse(sys.argv[1])
    print(len(e), "entries;", sum(1 for x in e if x["ref"] is None), "without ref")
    for x in e:
        if x["key"] in sys.argv[2:]:
            print(json.dumps(x))
