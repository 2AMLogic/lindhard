"""The one-off srim.org figures: measured ranges (SRIM/SRIMPICS/RANGES) and
stopping in compounds (SRIM/Compounds). Each has its own layout, so each is
described once in figures.json (axis labels read off the figure, the legend
rows as printed, clean examples of each marker); everything else is measured
from the image and checked:

- axes: the given label values are snapped to the figure's grid lines, and
  every other grid line must then land on a round value;
- legend: either rows found in the image on a fixed pitch and coloured from
  the plotting program's palette (which must agree with the colours seen),
  or rows given explicitly (symbol, colour, target, phase) where the figure
  marks data sets by letter or by target;
- citations: a printed key is attributed only if a citation list on the site
  has it and that entry names this ion and this target;
- markers: templates are the examples' own ink plus every isolated mark of
  about their size (for letters, assigned to the example letter it overlaps
  best), or drawn by hand where every marker touches a curve; they are
  matched in each legend colour, tolerating cover by other ink, and the best
  match wins where several fit one place; counts per colour are checked
  against the legend;
- where markers overlap too much for a blind match, a figure may list
  seeds (markers read by eye at 10x, then placed by the best template match
  within 4 px, or placed by eye where nothing fits) and rejects (matches
  read by eye as part of a neighbour). Every such point says so.

Usage: python3 -I figures.py <mirror-site-dir> <work-dir> [--ranges=<csv>]
       [--compounds=<csv>] [figure-id ...]
"""
import csv, json, math, os, re, sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import citpdf
import ppm

WHITE, BLACK, GRID = (255, 255, 255), (0, 0, 0), (85, 255, 255)
# The plotting program colours legend rows (and their markers) from this
# cycle, row by row; a column may start part way through it.
PALETTE = ["#0000FF", "#FF0000", "#FF00FF", "#009200", "#AA4900", "#FF9255", "#D524AA", "#FF0000",
           "#5500FF", "#55FF55", "#55FFFF", "#8092AA", "#5592FF", "#5500AA", "#AA4955", "#00FF00",
           "#55B600", "#FF00FF", "#00FFFF", "#000000"]
# Digits this legend font draws alike at its size.
CONFUSABLE = {"6": "8", "8": "69", "9": "8", "1": "7", "7": "1", "3": "5", "5": "36"}
METHOD = "srim_org-figures/1"
HERE = os.path.dirname(os.path.abspath(__file__))


def hexc(c):
    return "#%02X%02X%02X" % c


def rgb(h):
    return tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))


def frame(r, w, h):
    """The plot box: the innermost long black lines around the middle."""
    rows = [y for y in range(h) if sum(1 for x in range(w) if r[y][x] == BLACK) > w * 0.4]
    cols = [x for x in range(w) if sum(1 for y in range(h) if r[y][x] == BLACK) > h * 0.4]
    y0 = max(y for y in rows if y < h / 2)
    y1 = min(y for y in rows if y > h / 2)
    x0 = max(x for x in cols if x < w / 2)
    x1 = min(x for x in cols if x > w / 2)
    return x0, x1, y0, y1


def gridlines(r, box, colour, frac=0.3):
    x0, x1, y0, y1 = box
    vx = [x for x in range(x0 + 1, x1) if sum(1 for y in range(y0, y1) if r[y][x] == colour) > (y1 - y0) * frac]
    hy = [y for y in range(y0 + 1, y1) if sum(1 for x in range(x0, x1) if r[y][x] == colour) > (x1 - x0) * frac]
    return vx, hy


def snap(px, lines, edges, tol=6):
    """The grid line (or frame edge) nearest a label's position."""
    cands = list(lines) + list(edges)
    best = min(cands, key=lambda L: abs(L - px))
    if abs(best - px) > tol:
        raise ValueError("label at %s has no grid line within %d px" % (px, tol))
    return best


def axis(spec, lines, edges):
    """px -> value for one axis from two labelled grid lines. Returns
    (to_value, report): report lists every grid line's value and how far it
    is from a round value, in px."""
    (v1, p1), (v2, p2) = spec["refs"]
    a, b = snap(p1, lines, edges), snap(p2, lines, edges)
    log = spec["scale"] == "log"
    f = math.log10 if log else (lambda v: v)
    slope = (f(v2) - f(v1)) / (b - a)

    def to_value(px):
        u = f(v1) + (px - a) * slope
        return 10 ** u if log else u

    worst = 0.0
    for L in lines:
        v = to_value(L)
        if log:
            k = math.floor(math.log10(v) + 1e-9)
            m = v / 10 ** k
            nice = min(range(1, 11), key=lambda n: abs(m - n)) * 10 ** k
        else:
            step = spec["step"]
            nice = round(v / step) * step
        worst = max(worst, abs((f(nice) - f(v)) / slope))
    return to_value, {"refs_px": [a, b], "per_px": slope, "grid_worst_px": round(worst, 2), "lines": len(lines)}


def legend_rows(r, col, vx, hy, grid):
    """Colour of each text row in one legend column. Rows sit on a fixed
    pitch: col = [x_from, x_to, first_top, pitch, n]. Grid pixels are
    ignored; a row whose only ink is the grid colour off the grid lines is in
    that colour. Raises if a row is empty or the row after the last is not."""
    xa, xb, top, pitch, n = col
    skip_x, skip_y = set(vx), set(hy)

    def band(t):
        return Counter(r[y][x] for y in range(t, t + 7) for x in range(xa, xb)
                       if x not in skip_x and y not in skip_y and r[y][x] != WHITE)

    out = []
    for k in range(n + 1):
        t = round(top + k * pitch)
        b = band(t)
        if k == n:
            if sum(b.values()) > 2:
                raise ValueError("legend column %s: ink below its last row (at %d)" % (col, t))
            break
        if not b:
            raise ValueError("legend column %s: row %d (at %d) is empty" % (col, k, t))
        other = [c for c, _ in b.most_common() if c != grid]
        out.append((t, other[0] if other else grid))
    return out


def template(r, ex):
    """The example marker's ink: offsets from its box, and its colour."""
    x, y, w, h = ex
    c = Counter(r[yy][xx] for yy in range(y, y + h) for xx in range(x, x + w)
                if r[yy][xx] not in (WHITE,)).most_common(1)[0][0]
    return [(i, j) for j in range(h) for i in range(w) if r[y + j][x + i] == c], c


def components(r, box, colours, join):
    """Connected ink of one colour (pixels up to join apart are connected):
    (colour, x, y, w, h, frozenset of offsets)."""
    x0, x1, y0, y1 = box
    seen, out = set(), []
    for y in range(y0 + 1, y1):
        for x in range(x0 + 1, x1):
            c = r[y][x]
            if c not in colours or (x, y) in seen:
                continue
            stack, px = [(x, y)], []
            seen.add((x, y))
            while stack:
                a, b = stack.pop()
                px.append((a, b))
                for da in range(-join, join + 1):
                    for db in range(-join, join + 1):
                        q = (a + da, b + db)
                        if q not in seen and x0 < q[0] < x1 and y0 < q[1] < y1 and r[q[1]][q[0]] == c:
                            seen.add(q)
                            stack.append(q)
            xs, ys = [p[0] for p in px], [p[1] for p in px]
            mx, my = min(xs), min(ys)
            out.append((c, mx, my, max(xs) - mx + 1, max(ys) - my + 1,
                        frozenset((a - mx, b - my) for a, b in px)))
    return out


def learn_templates(r, box, colours, example, join, exclude):
    """The raster variants of one marker: the example's ink, plus the pixel
    mask of every isolated component of about its size and ink (isolated:
    nothing but white or grid around its box). Each such component is a
    marker drawn in the clear, so its own mask is a variant."""
    tmpl, c = template(r, example)
    _, _, w, h = example
    masks = Counter()
    for col, x, y, bw, bh, m in components(r, box, colours, join):
        if not (abs(bw - w) <= 2 and abs(bh - h) <= 1 and 0.75 * len(tmpl) <= len(m) <= 1.25 * len(tmpl)):
            continue
        if any(ex[0] <= x + bw / 2 <= ex[2] and ex[1] <= y + bh / 2 <= ex[3] for ex in exclude):
            continue
        ring = [r[yy][xx] for yy in range(y - 2, y + bh + 2) for xx in range(x - 2, x + bw + 2)
                if not (y <= yy < y + bh and x <= xx < x + bw)]
        if all(k in (WHITE, GRID) for k in ring):
            masks[m] += 1
    out = [sorted(tmpl)]
    for m, n in masks.items():
        if sorted(m) not in out:
            out.append(sorted(m))
    return out


def iou(a, b):
    """Best overlap of two pixel masks over shifts of up to one pixel."""
    a = set(a)
    best = 0.0
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            m = {(i + dx, j + dy) for i, j in b}
            best = max(best, len(a & m) / len(a | m))
    return best


def letter_variants(r, box, examples, join, exclude, min_iou=0.6):
    """Variants for letter markers, which share one size so cannot be told
    apart by it: every isolated component of a letter's colour and about a
    letter's size is a variant of the example letter it overlaps best, if
    that overlap is at least min_iou. examples: {label: [(mask, colour)]}."""
    out = {k: [sorted(m) for m, _ in v] for k, v in examples.items()}
    sizes = [(max(i for i, _ in m) + 1, max(j for _, j in m) + 1) for v in examples.values() for m, _ in v]
    colours = {c for v in examples.values() for _, c in v}
    for col, x, y, bw, bh, m in components(r, box, colours, join):
        if not any(abs(bw - w) <= 1 and abs(bh - h) <= 1 for w, h in sizes):
            continue
        if any(ex[0] <= x + bw / 2 <= ex[2] and ex[1] <= y + bh / 2 <= ex[3] for ex in exclude):
            continue
        ring = [r[yy][xx] for yy in range(y - 2, y + bh + 2) for xx in range(x - 2, x + bw + 2)
                if not (y <= yy < y + bh and x <= xx < x + bw)]
        if not all(k in (WHITE, GRID) for k in ring):
            continue
        score, label = max((iou(t, m), k) for k, v in examples.items() for t, c in v if c == col)
        if score >= min_iou and sorted(m) not in out[label]:
            out[label].append(sorted(m))
    return out


def find(r, box, tmpl, w, h, colours, exclude, need=0.8, max_white=0.12):
    x0, x1, y0, y1 = box
    n = len(tmpl)
    # A cheap first look: every third template pixel decides the colour.
    probe = tmpl[::3]
    hits = []
    for y in range(y0 + 1, y1 - h):
        for x in range(x0 + 1, x1 - w):
            ink = Counter(r[y + j][x + i] for i, j in probe)
            c = max((k for k in ink if k in colours), key=lambda k: ink[k], default=None)
            if c is None or ink[c] < need * len(probe) - 1:
                continue
            on = sum(1 for i, j in tmpl if r[y + j][x + i] == c)
            if on < need * n:
                continue
            white = sum(1 for i, j in tmpl if r[y + j][x + i] == WHITE)
            if white > max_white * n:
                continue
            if any(ex[0] <= x + w / 2 <= ex[2] and ex[1] <= y + h / 2 <= ex[3] for ex in exclude):
                continue
            # Rank by the fraction of the template matched, less the ink of
            # that colour in and just around its box that it does not explain:
            # an F inside an E matches all of the F but leaves the E's bottom
            # bar over, and an F placed low on an E leaves its top bar over.
            extra = sum(1 for j in range(-2, h + 2) for i in range(-2, w + 2)
                        if r[y + j][x + i] == c) - on
            hits.append(((on - extra) / n, x, y, c))
    hits.sort(reverse=True)
    kept = []
    for s, x, y, c in hits:
        if all(c != k[3] or abs(x - k[1]) > 2 or abs(y - k[2]) > 2 for k in kept):
            kept.append((s, x, y, c))
    return [(sc, x, y, c) for sc, x, y, c in kept]


def names_system(entry, words):
    text = " ".join(entry["lines"])
    return all(re.search(wd, text) for wd in words)


def all_citations(site, cache_path):
    """key -> entry over every citation list on the site, cached as JSON."""
    if cache_path and os.path.exists(cache_path):
        return json.load(open(cache_path))
    out = {}
    for d in ("Citations", "SRIM/Compounds/Citations"):
        for f in sorted(os.listdir(os.path.join(site, d))):
            if f.endswith(".pdf"):
                try:
                    for e in citpdf.parse(os.path.join(site, d, f)):
                        out.setdefault(e["key"], dict(e, file=d + "/" + f))
                except Exception as err:  # one unreadable list must not stop the rest
                    print("%s: %s" % (f, err), file=sys.stderr)
    if cache_path:
        json.dump(out, open(cache_path, "w"))
    return out


def variants(key):
    """Keys one confusable digit, or one author letter, away from key."""
    m = re.match(r"(\d{4})-(.*)$", key)
    if not m:
        return set()
    year, rest = m.groups()
    out = set()
    for i, d in enumerate(year):
        for alt in CONFUSABLE.get(d, ""):
            out.add(year[:i] + alt + year[i + 1:] + "-" + rest)
    return out


def edit1(a, b):
    if abs(len(a) - len(b)) > 1:
        return False
    if len(a) == len(b):
        return sum(x != y for x, y in zip(a, b)) <= 1
    if len(a) > len(b):
        a, b = b, a
    return any(a == b[:i] + b[i + 1:] for i in range(len(b)))


def resolve(printed, local, everywhere, words, local_is_system=False):
    """(key, status, entry). status: exact | corrected (a confusable digit or
    one author letter) | unconfirmed (the key is in a citation list kept for
    this one system, but its entry does not name the system) | key-mismatch
    (the key exists, for another system) | not-found."""
    ok = {k: e for k, e in local.items() if names_system(e, words)}
    if printed in ok:
        return printed, "exact", ok[printed]
    cands = [k for k in variants(printed) if k in ok]
    if not cands:
        year, _, author = printed.partition("-")
        cands = [k for k in ok if k.startswith(year + "-") and edit1(k[len(year) + 1:], author)]
    if len(cands) == 1:
        return cands[0], "corrected", ok[cands[0]]
    if local_is_system and printed in local:
        return printed, "unconfirmed", local[printed]
    if printed in everywhere:
        return printed, "key-mismatch", everywhere.get(printed) or local.get(printed)
    return printed, "not-found", None


def run_figure(site, fig, cit_cache):
    w, h, r = ppm.load(os.path.join(site, fig["file"]))
    box = frame(r, w, h)
    grid = rgb(fig.get("grid", "#55FFFF"))
    vx, hy = gridlines(r, box, grid)
    fx, rx = axis(fig["x"], vx, (box[0], box[1]))
    fy, ry = axis(fig["y"], hy, (box[2], box[3]))

    # Legend. Palette mode: rows found on a fixed pitch, coloured from the
    # palette cycle, zipped with the transcription; every row has the
    # figure's one marker symbol. Explicit mode: each row names its symbol,
    # colour and anything else that row decides (target, phase).
    rows = []
    leg = fig["legend"]
    if "columns" in leg:
        sym = fig["markers"][0]["label"]
        for col, printed in zip(leg["columns"], leg["rows"]):
            xa, xb, top, pitch, start = col
            found = legend_rows(r, [xa, xb, top, pitch, len(printed)], vx, hy, grid)
            for k, ((t, c), (key, n)) in enumerate(zip(found, printed)):
                rows.append({"symbol": sym, "key_printed": key, "count": n,
                             "colour": PALETTE[(start + k) % len(PALETTE)], "colour_seen": hexc(c), "top": t})
        agree = sum(row["colour"] == row["colour_seen"] for row in rows)
        if agree < 0.9 * len(rows):
            raise ValueError("%s: legend colours disagree with the palette (%d of %d)" % (fig["id"], agree, len(rows)))
    else:
        for row in leg["rows"]:
            rows.append(dict(row, key_printed=row.get("key")))

    # Citations: the printed key, checked against the site's lists.
    local = {}
    for p in fig.get("citations", []):
        if p not in cit_cache:
            cit_cache[p] = citpdf.parse(os.path.join(site, p))
        for e in cit_cache[p]:
            local.setdefault(e["key"], e)
    for row in rows:
        if not row["key_printed"]:
            # No key: the figure prints its one reference instead.
            row.update(key=None, status="not-printed", reference=row.get("reference"))
            continue
        words = row.get("system_words", fig["system_words"])
        key, status, e = resolve(row["key_printed"], local, cit_cache["*"], words,
                                 fig.get("citations_are_system", False))
        row.update(key=key, status=status, reference=" ".join(e["lines"]) if e else None)

    # Markers.
    colours = {rgb(row["colour"]) for row in rows} | {rgb(c) for c in fig.get("extra_colours", [])}
    # Every symbol's templates (its examples, plus learned variants), matched
    # in the colours its legend rows use. A marker matched by several
    # templates, of one symbol or of several, is one point: the best match
    # wins (more ink matched, less white), so an 'E' is not also an 'F'.
    hits, variants, by_label = [], {}, {}
    letters = {}
    if fig.get("letters"):
        letters = letter_variants(r, box, {sym["label"]: [template(r, ex) for ex in sym["examples"]]
                                           for sym in fig["markers"]},
                                  fig.get("join", 1), fig.get("exclude", []))
    for sym in fig["markers"]:
        cols = {rgb(row["colour"]) for row in rows if row["symbol"] == sym["label"]}
        tmpls = letters.get(sym["label"], [])
        if "ascii" in sym:
            # Templates drawn by hand, for markers no curve leaves clear: one
            # per raster variant of the glyph.
            tmpls = [sorted((i, j) for j, line in enumerate(g) for i, ch in enumerate(line) if ch == "#")
                     for g in sym["ascii"]]
        for ex in ([] if letters or "ascii" in sym else sym.get("examples", [sym.get("example")])):
            for t in learn_templates(r, box, cols, ex, sym.get("join", 1), fig.get("exclude", [])):
                if t not in tmpls:
                    tmpls.append(t)
        variants[sym["label"]] = len(tmpls)
        by_label[sym["label"]] = (tmpls, cols, sym.get("anchor"))
        for k, t in enumerate(tmpls):
            tw = max(i for i, _ in t) + 1
            th = max(j for _, j in t) + 1
            # The plotted point: the template's box centre, unless the template
            # is only part of the label (then "anchor" is the label's centre).
            ax, ay = sym.get("anchor", [tw / 2, th / 2])
            for sc, x, y, c in find(r, box, t, tw, th, cols, fig.get("exclude", []),
                                    fig.get("need", 0.8), fig.get("max_white", 0.12)):
                hits.append((sc, x + ax, y + ay, c, sym["label"], k, x, y, "auto"))
    hits.sort(key=lambda t: -t[0])
    kept = []
    ddx, ddy = fig.get("dedupe_px", [3, 3])
    for h_ in hits:
        if all(h_[3] != k_[3] or abs(h_[1] - k_[1]) > ddx or abs(h_[2] - k_[2]) > ddy for k_ in kept):
            kept.append(h_)
    # Seeds: markers read by eye where letters overlap each other or a curve
    # too much for a blind match. Each is placed at the best match of its
    # letter's templates within 4 px of where it was read, so the position
    # still comes from the raster, not from the reading. A seed may name
    # several letters ("B|D") where the visible part of a marker cannot tell
    # them apart: the best match picks one and the point records the choice.
    # A seed marked "eye" has no template to snap to (a letter with no clean
    # example, or one so overlapped that no template fits), so its box centre
    # is the reading, and an eye seed naming several letters keeps them all
    # as its symbol. Any other fourth element is a note kept with the point.
    #
    # Rejects: blind matches read by eye as no marker (part of a neighbour),
    # each [x, y, why] at the match's top-left.
    kept = [h_ + (None,) for h_ in kept
            if not any(abs(h_[6] - rx_) <= 2 and abs(h_[7] - ry_) <= 2 for rx_, ry_, _ in fig.get("reject", []))]
    for seed in fig.get("seeds", []):
        label, sx, sy = seed[:3]
        note = next((e for e in seed[3:] if e != "eye"), None)
        labels = label.split("|")
        if "eye" in seed[3:]:
            c = rgb(next(row["colour"] for row in rows if row["symbol"] == labels[0]))
            kept.append((None, sx, sy, c, label, None, sx - 4.5, sy - 6, "eye", note))
            continue
        best = None
        for lab in labels:
            tmpls, cols, anchor = by_label[lab]
            for k, t in enumerate(tmpls):
                tw = max(i for i, _ in t) + 1
                th = max(j for _, j in t) + 1
                ax, ay = anchor or [tw / 2, th / 2]
                for c in cols:
                    for y in range(sy - 4, sy + 5):
                        for x in range(sx - 4, sx + 5):
                            on = sum(1 for i, j in t if r[y + j][x + i] == c)
                            extra = sum(1 for j in range(-2, th + 2) for i in range(-2, tw + 2)
                                        if r[y + j][x + i] == c) - on
                            sc = (on - extra) / len(t)
                            if best is None or sc > best[0]:
                                best = (sc, x + ax, y + ay, c, lab, k, x, y, "seeded",
                                        note or (("one of " + label) if len(labels) > 1 else None))
        # A seed also corrects a blind match it lands on (an E crossed by a
        # curve can match a B better than an E).
        kept = [k_ for k_ in kept if not (best[3] == k_[3] and abs(best[1] - k_[1]) <= ddx
                                          and abs(best[2] - k_[2]) <= ddy)]
        kept.append(best)
    pts = []
    for sc, cx_, cy_, c, label, k, x, y, how, note in kept:
        mine = [row for row in rows if row["symbol"] in label.split("|") and row["colour"] == hexc(c)]
        pts.append({"symbol": label, "colour": hexc(c), "px": [x, y], "variant": k, "how": how,
                    "score": None if sc is None else round(sc, 3), "note": note,
                    "x": fx(cx_), "y": fy(cy_), "rows": [rows.index(row) for row in mine]})
    count = Counter((p["symbol"], p["colour"]) for p in pts)
    # A legend that prints no counts leaves nothing to check against (None).
    want = defaultdict(int)
    for row in rows:
        if row["count"] is None:
            want[(row["symbol"], row["colour"])] = None
        elif want[(row["symbol"], row["colour"])] is not None:
            want[(row["symbol"], row["colour"])] += row["count"]
    check = [{"symbol": s_, "colour": c, "legend": want.get((s_, c), 0), "found": count.get((s_, c), 0)}
             for s_, c in sorted(set(want) | set(count))]
    return {"id": fig["id"], "file": fig["file"], "frame": box, "x_axis": rx, "y_axis": ry,
            "legend": rows, "count_check": check, "points": pts, "variants": variants,
            "legend_points": sum(v or 0 for v in want.values()), "found": len(pts)}


COLUMNS = ["figure", "quantity", "ion_z", "ion", "target", "phase", "energy", "energy_unit", "value",
           "value_unit", "u_energy_rel", "u_value", "symbol", "colour", "key_printed", "citation_key",
           "key_status", "reference", "candidate_key", "how", "note", "pixel_x", "pixel_y", "origin", "method"]


def sig(v, n=4):
    return round(v, n - 1 - int(math.floor(math.log10(abs(v))))) if v else v


def export(figs, outs, path):
    """One CSV row per point. Uncertainty is read-off only: half a pixel,
    the anchor's (u_anchor_px; 0.7 px by default, from the range figures'
    round energies), and 1.5 px more for a point placed by eye. A point whose
    colour (and symbol) several legend rows share lists every such row: the
    key, status and reference columns are "|"-joined in the same row order."""
    w = csv.writer(open(path, "w", newline=""), lineterminator="\n")
    w.writerow(COLUMNS)
    n = 0
    for fig in figs:
        o = outs[fig["id"]]
        for p in sorted(o["points"], key=lambda p: (p["colour"], p["x"], p["y"])):
            rows = [o["legend"][i] for i in p["rows"]]
            join = lambda k: "|".join(sorted({str(row.get(k) or "") for row in rows} - {""}))
            each = lambda k: "|".join(str(row.get(k) or "") for row in rows)
            u_px = math.hypot(0.5, fig.get("u_anchor_px", 0.7))
            if p["how"] == "eye":
                u_px = math.hypot(u_px, 1.5)
            ux, uy = abs(o["x_axis"]["per_px"]) * u_px, abs(o["y_axis"]["per_px"]) * u_px
            u_e = math.log(10) * ux if fig["x"]["scale"] == "log" else ux / p["x"]
            u_v = p["y"] * math.log(10) * uy if fig["y"]["scale"] == "log" else uy
            w.writerow([fig["id"], fig["quantity"], fig["ion_z"], fig["ion"], join("target") or fig["target"],
                        join("phase"), sig(p["x"]), fig["x"]["unit"], sig(p["y"]), fig["y"]["unit"],
                        sig(u_e, 2), sig(u_v, 2), p["symbol"], p["colour"], each("key_printed"),
                        each("key"), each("status"), each("reference"), join("candidate_key"), p["how"],
                        p["note"] or "", p["px"][0], p["px"][1], "digitized-compilation", METHOD])
            n += 1
    print("%d points -> %s" % (n, path), file=sys.stderr)


def main():
    site, out = sys.argv[1], sys.argv[2]
    opts = dict(a[2:].split("=", 1) for a in sys.argv[3:] if a.startswith("--"))
    only = {a for a in sys.argv[3:] if not a.startswith("--")}
    figs = [f for f in json.load(open(os.path.join(HERE, "figures.json"))) if not only or f["id"] in only]
    os.makedirs(out, exist_ok=True)
    cache = {"*": all_citations(site, os.path.join(out, "citations.json"))}
    outs = {}
    for fig in figs:
        o = outs[fig["id"]] = run_figure(site, fig, cache)
        json.dump(o, open(os.path.join(out, fig["id"] + ".json"), "w"), indent=1)
        st = Counter(r["status"] for r in o["legend"])
        print("%s: %d markers, legend %d; grid worst %.2f/%.2f px; keys %s"
              % (fig["id"], o["found"], o["legend_points"], o["x_axis"]["grid_worst_px"],
                 o["y_axis"]["grid_worst_px"], dict(st)), file=sys.stderr)
    # Ranges and stopping in compounds go to separate files.
    for opt, quantity in (("ranges", "range"), ("compounds", "stopping_cross_section")):
        if opt in opts:
            export([f for f in figs if f["quantity"] == quantity], outs, opts[opt])


if __name__ == "__main__":
    main()
