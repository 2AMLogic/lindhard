"""Markers drawn as the ion's name only ("He", "Li"): match the legend's own
symbol bitmap in each colour. Pixels under other ink are unknown; white where
the glyph has ink is evidence against."""
WHITE, BLACK = (255, 255, 255), (0, 0, 0)


def stroke(r, x, y, c, n=8, need=6, white=1):
    col = [r[y + j][x] for j in range(n)]
    return col.count(c) >= need and col.count(WHITE) <= white


def find_H(r, x0, x1, y0, y1, colours, tail=True, max_inner=10, need=4, max_white=2):
    """'H' glyphs (8 tall, strokes 5 apart, crossbar) in each colour; with
    tail, also some ink of the same colour right after (the 'e' of 'He')."""
    hits = []
    for y in range(y0 + 1, y1 - 9):
        for x in range(x0 + 1, x1 - 14):
            c = r[y][x]
            if c not in colours:
                continue
            if not (stroke(r, x, y, c, need=need, white=max_white) and stroke(r, x + 5, y, c, need=need, white=max_white)):
                continue
            if r[y][x - 1] == c and r[y + 4][x - 1] == c:
                continue  # part of a thicker stroke
            bars = [j for j in range(2, 6) if sum(1 for i in range(1, 5) if r[y + j][x + i] == c) >= 3]
            if not bars:
                continue
            j = max(bars, key=lambda b: sum(1 for i in range(1, 5) if r[y + b][x + i] == c))
            inner = sum(1 for jj in range(8) if jj != j for i in range(1, 5) if r[y + jj][x + i] == c)
            if inner > max_inner:
                continue
            if tail:
                ink = sum(1 for jj in range(2, 8) for i in range(7, 13) if r[y + jj][x + i] == c)
                if ink < 5:
                    continue
            hits.append((x, y, c, j))
    kept = []
    for h in hits:
        if all(h[2] != k[2] or abs(h[0] - k[0]) > 3 or abs(h[1] - k[1]) > 3 for k in kept):
            kept.append(h)
    return kept
