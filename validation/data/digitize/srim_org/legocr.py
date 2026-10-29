"""Read srim.org legend rows with templates of the plot's own bitmap font.

A legend row is: marker, ion label, citation key, total points, points below
and at/above 1 MeV/u, and the old citation number. Glyphs are 8 px tall and
separated by blank (all-white) columns; fields by wider gaps."""
WHITE = (255, 255, 255)
import json as _json
import os as _os
DIGITS = {k: v for k, v in _json.load(open(_os.path.join(_os.path.dirname(__file__), "digits.json"))).items()}


def field_mask(r, a, b, top, c):
    """1 = ink of this row's colour, 0 = white, None = under other ink."""
    return [[(1 if r[top + j][x] == c else (0 if r[top + j][x] == WHITE else None)) for x in range(a, b)]
            for j in range(8)]


def score(m, s):
    """Lowest mismatch count for drawing string s anywhere in field mask m,
    with 0-2 blank columns between glyphs; ink outside the glyphs counts too."""
    W = len(m[0])
    ink_total = sum(1 for j in range(8) for i in range(W) if m[j][i] == 1)
    best = None
    # dp over glyph index: dict end_x -> (cost, ink_covered)
    states = {}
    for bm in DIGITS[s[0]]:
        w0 = len(bm[0])
        for x in range(0, W - w0 + 1):
            v = glyph_cost(m, bm, x)
            if x + w0 not in states or v[0] - v[1] < states[x + w0][0] - states[x + w0][1]:
                states[x + w0] = v
    for ch in s[1:]:
        nxt = {}
        for bm in DIGITS[ch]:
            tw = len(bm[0])
            for end, (cost, cov) in states.items():
                for g in range(0, 6):
                    x = end + g
                    if x + tw > W:
                        continue
                    gc, gcov = glyph_cost(m, bm, x)
                    v = (cost + gc, cov + gcov)
                    if x + tw not in nxt or v[0] - v[1] < nxt[x + tw][0] - nxt[x + tw][1]:
                        nxt[x + tw] = v
        states = nxt
    for cost, cov in states.values():
        tot = cost + (ink_total - cov)
        if best is None or tot < best:
            best = tot
    return best if best is not None else 999


def glyph_cost(m, bm, x):
    cost = cov = 0
    for j in range(8):
        for i in range(len(bm[0])):
            v = m[j][x + i]
            if v is None:
                continue
            on = bm[j][i] == "O"
            if v == 1:
                cov += 1
                if not on:
                    cost += 1
            elif on:
                cost += 1
    return cost, cov


def viterbi(m, maxlen=5):
    """Best digit string for field mask m: glyphs left to right, 0-5 blank
    columns between them, leftover ink penalised. Returns (string, cost)."""
    W = len(m[0])
    ink_total = sum(1 for j in range(8) for i in range(W) if m[j][i] == 1)
    # Free leading blank: a glyph may start anywhere. state: end_x -> (cost-cov, cost, cov, s)
    states = {}
    variants = [(ch, bm) for ch, bms in DIGITS.items() for bm in bms]
    for ch, bm in variants:
        tw = len(bm[0])
        for x in range(0, W - tw + 1):
            c, v = glyph_cost(m, bm, x)
            k = x + tw
            if k not in states or c - v < states[k][0]:
                states[k] = (c - v, c, v, ch)
    finals = list(states.values())
    for _ in range(maxlen - 1):
        nxt = {}
        for end, (_, c0, v0, s0) in states.items():
            for ch, bm in variants:
                tw = len(bm[0])
                for g in range(0, 6):
                    x = end + g
                    if x + tw > W:
                        break
                    c, v = glyph_cost(m, bm, x)
                    k = x + tw
                    val = (c0 + c) - (v0 + v)
                    if k not in nxt or val < nxt[k][0]:
                        nxt[k] = (val, c0 + c, v0 + v, s0 + ch)
        states = nxt
        finals += list(states.values())
    if not finals:
        return "", 999
    best = min(finals, key=lambda t: t[1] + (ink_total - t[2]))
    return best[3], best[1] + (ink_total - best[2])


def alternatives(s):
    """s and every string one digit substitution, insertion or deletion away."""
    out = {s}
    D = "0123456789"
    for i in range(len(s)):
        out.add(s[:i] + s[i + 1:])
        for ch in D:
            out.add(s[:i] + ch + s[i + 1:])
    for i in range(len(s) + 1):
        for ch in D:
            out.add(s[:i] + ch + s[i:])
    return {t for t in out if t and len(t) <= 4 and (t == "0" or t[0] != "0")}


# The legend layout is fixed relative to the frame: field spans (total, < 1,
# >= 1, old ref) as offsets from the frame's right edge, and the depth of the
# heading rows below its top edge.
FIELD_SPANS = [(-108, -90), (-88, -61), (-56, -38), (-37, -1)]
HEADER_PX = 44
ROW0_PX, ROW_PITCH = 48, 10
GRIDC = (85, 255, 255)


def grid_rows(r, x0, x1, y0, y1, vx, hy):
    """(top, colour) for each legend row on the fixed grid. A row's colour is
    the commonest ink in its number columns; three empty rows end the legend."""
    from collections import Counter
    vxs, hys = set(vx), set(hy)
    rows, empty = [], 0
    t = y0 + ROW0_PX
    while t + 8 < y1 and empty < 3:
        cnt = Counter()
        for a, b in FIELD_SPANS:
            for x in range(x1 + a, x1 + b):
                if x in vxs:
                    continue
                for j in range(8):
                    if t + j in hys:
                        continue
                    p = r[t + j][x]
                    if p != WHITE:
                        cnt[p] += 1
        if sum(cnt.values()) >= 6:
            rows.append((t, cnt.most_common(1)[0][0]))
            empty = 0
        else:
            empty += 1
        t += ROW_PITCH
    return rows


def read_legend(r, x0, x1, y0, y1, vx, hy, refs=(), width=200):
    out = []
    rows = grid_rows(r, x0, x1, y0, y1, vx, hy)
    for top, c in rows:
        masks = []
        for a, b in FIELD_SPANS:
            a, b = x1 + a, x1 + b
            m = field_mask(r, a, b, top, c)
            if c == GRIDC:
                for x in vx:
                    if a <= x < b:
                        for j in range(8):
                            m[j][x - a] = None
            masks.append(m)
        res = decode_fields(masks, refs)
        if res is None:
            continue
        res.update(top=top, colour="#%02X%02X%02X" % c)
        out.append(res)
    return out


def decode_fields(masks, refs):
    free = [viterbi(m) for m in masks]
    if not all(s for s, _ in free):
        return None
    # A long citation key can spill into the total column and add a leading
    # digit; the < 1 and >= 1 columns are clear of it, so when the total read
    # ends in their sum, the sum is the total.
    lt, ge = free[1][0], free[2][0]
    if lt.isdigit() and ge.isdigit():
        sm = str(int(lt) + int(ge))
        if free[0][0] != sm and free[0][0].endswith(sm):
            free[0] = (sm, free[0][1])
    cand = []
    for m, (s, cost) in zip(masks[:3], free[:3]):
        opts = {t: score(m, t) for t in alternatives(s)}
        cand.append(sorted(opts.items(), key=lambda kv: kv[1])[:12])
    best = None
    for t, ct in cand[0]:
        for l, cl in cand[1]:
            for g, cg in cand[2]:
                if int(t) == int(l) + int(g):
                    tot = ct + cl + cg
                    if best is None or tot < best[0]:
                        best = (tot, t, l, g)
    free_cost = sum(c for _, c in free[:3])
    out = {"free": [s for s, _ in free], "free_cost": free_cost}
    if best is not None:
        out.update(total=int(best[1]), lt1=int(best[2]), ge1=int(best[3]), count_cost=best[0],
                   count_corrected=[best[1], best[2], best[3]] != [s for s, _ in free[:3]])
    else:
        out.update(total=None, lt1=None, ge1=None, count_cost=None, count_corrected=False)
    s, cost = free[3]
    ref, ref_cost = int(s), cost
    if refs and ref not in refs:
        opts = [(score(masks[3], t), t) for t in alternatives(s) if int(t) in refs]
        if opts:
            oc, ot = min(opts)
            if oc <= cost + 6:
                ref, ref_cost = int(ot), oc
    out.update(ref=ref, ref_cost=ref_cost, ref_in_pdf=(ref in refs) if refs else None)
    return out
