"""Minimal PPM (P6) reader and helpers for the srim.org plot digitizer."""
import subprocess


def load(gif_path):
    """Return (width, height, rows) where rows[y][x] is an (r, g, b) tuple."""
    data = subprocess.run(["magick", gif_path, "-depth", "8", "ppm:-"], capture_output=True, check=True).stdout
    # Header: P6\n<w> <h>\n<max>\n then raw bytes.
    parts, i = [], 0
    while len(parts) < 4:
        while data[i:i + 1].isspace():
            i += 1
        if data[i:i + 1] == b"#":
            while data[i:i + 1] != b"\n":
                i += 1
            continue
        j = i
        while not data[j:j + 1].isspace():
            j += 1
        parts.append(data[i:j])
        i = j
    i += 1
    w, h = int(parts[1]), int(parts[2])
    px = data[i:]
    rows = [[tuple(px[(y * w + x) * 3:(y * w + x) * 3 + 3]) for x in range(w)] for y in range(h)]
    return w, h, rows


def ascii_art(rows, x0, y0, x1, y1, legend=None):
    """Print a region, one character per colour (legend maps colour -> char)."""
    legend = dict(legend or {})
    legend.setdefault((255, 255, 255), ".")
    chars = iter("abcdefghijklmnopqrstuvwxyz")
    for y in range(y0, y1):
        line = ""
        for x in range(x0, x1):
            c = rows[y][x]
            if c not in legend:
                legend[c] = next(chars)
            line += legend[c]
        print(f"{y:4d} {line}")
    return legend
