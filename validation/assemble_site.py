#!/usr/bin/env python3
"""Put the landing page (site/) over the root of the rendered book.

mdBook renders book/src/introduction.md to index.html. The landing page takes
that URL; the introduction stays reachable at introduction.html, so no
published book URL breaks. The trajectory simulator (SIMULATOR_PAGE) ships in
every build, so a missing simulator page is an error, not a fallback. The API
docs link (API_PAGE, published from target/doc by book.yml) is optional: it is
rendered only when its target exists, so the page never carries a dead link.

Usage: validation/assemble_site.py [BOOK_DIR]     (default: book/book)
Python standard library only.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"

# The one place the simulator page's path lives.
SIMULATOR_PAGE = "calculator/simulator.html"
API_PAGE = "api/lindhard/index.html"


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "book" / "book"
    if not (out / "introduction.html").is_file():
        print(f"error: {out} is not a rendered book (no introduction.html)", file=sys.stderr)
        return 1
    html = (SITE / "index.html").read_text(encoding="utf-8")
    if not (out / SIMULATOR_PAGE).is_file():
        print(f"error: {out} lacks the simulator page {SIMULATOR_PAGE}", file=sys.stderr)
        return 1
    sim = f'<a href="{SIMULATOR_PAGE}">Trajectory simulator</a>: ion trajectories in your browser.'
    if (out / API_PAGE).is_file():
        api = f'<a href="{API_PAGE}">API documentation</a> of the <code>lindhard</code> library (rustdoc)'
    else:
        api = "API documentation: <code>cargo doc -p lindhard --open</code>"
    for marker, text in (("<!--SIMULATOR_SLOT-->", sim), ("<!--API_LINK-->", api)):
        if marker not in html:
            print(f"error: site/index.html lacks {marker}", file=sys.stderr)
            return 1
        html = html.replace(marker, text)
    (out / "index.html").write_text(html, encoding="utf-8")
    shutil.copyfile(SITE / "style.css", out / "landing.css")
    print(f"landing page written to {out / 'index.html'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
