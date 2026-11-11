#!/usr/bin/env python3
"""Put the landing page (site/) over the root of the rendered book.

mdBook renders book/src/introduction.md to index.html. The landing page takes
that URL; the introduction stays reachable at introduction.html, so no
published book URL breaks. The trajectory simulator (SIMULATOR_PAGE) ships in
every build, so a missing simulator page is an error, not a fallback. The API
docs link (API_PAGE, published from target/doc by book.yml) is optional: it is
rendered only when its target exists, so the page never carries a dead link.

It also puts the shared site chrome on every surface (site/tokens.css,
site/chrome.css, site/header.html, site/footer.html): the files are copied to
the root of the output and to calculator/, the landing page and the calculator
pages get the header and footer at their markers, and the book pages (which
carry the header from book/theme/header.hbs) get the footer, and every page its
aria-current link. The tool pages still open standalone from a checkout: their
markers are plain comments there.

Usage: validation/assemble_site.py [BOOK_DIR]     (default: book/book)
Python standard library only.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"

# The one place the simulator page's path lives.
SIMULATOR_PAGE = "calculator/simulator.html"
API_PAGE = "api/lindhard/index.html"

CALCULATOR_DIR = "calculator"
# Pages (under CALCULATOR_DIR) that carry the markers below.
TOOL_PAGES = ("index.html", "simulator.html")
ROOT_TOKEN = "{{path_to_root}}"  # mdBook's own variable, so header.hbs is a byte copy
HEAD_LINKS = '<link rel="stylesheet" href="tokens.css">\n<link rel="stylesheet" href="chrome.css">'


def section_of(rel: str) -> str | None:
    """The header link that is the current page, by the page's path in the site."""
    if rel == "index.html":
        return "home"
    if rel == f"{CALCULATOR_DIR}/index.html":
        return "calculator"
    if rel == f"{CALCULATOR_DIR}/simulator.html" or rel == SIMULATOR_PAGE:
        return "simulator"
    if rel in ("validation.html", "status.html"):
        return rel.removesuffix(".html")
    if rel.startswith(("api/", f"{CALCULATOR_DIR}/")):
        return None
    return "manual"


def root_of(rel: str) -> str:
    return "../" * rel.count("/")


def mark_current(html: str, rel: str) -> str:
    """Turn data-cur="key" into aria-current="page" on the current section's link."""
    current = section_of(rel)
    return re.sub(
        r' data-cur="(\w+)"', lambda m: ' aria-current="page"' if m.group(1) == current else "", html
    )


def chrome(html: str, rel: str) -> str:
    """Fill the SITE_* markers of a page that has them."""
    header = (SITE / "header.html").read_text(encoding="utf-8")
    footer = (SITE / "footer.html").read_text(encoding="utf-8")
    for marker, text in (
        ("<!--SITE_HEAD-->", HEAD_LINKS),
        ("<!--SITE_HEADER-->", header),
        ("<!--SITE_FOOTER-->", footer),
    ):
        html = html.replace(marker, text)
    return mark_current(html.replace(ROOT_TOKEN, root_of(rel)), rel)


def finish_book_pages(out: Path) -> int:
    """Footer and aria-current on every book page that carries the site bar."""
    footer = (SITE / "footer.html").read_text(encoding="utf-8")
    count = 0
    for page in sorted(out.rglob("*.html")):
        rel = page.relative_to(out).as_posix()
        if rel.startswith(("api/", f"{CALCULATOR_DIR}/")) or rel == "index.html":
            continue
        html = page.read_text(encoding="utf-8")
        if 'class="lh-bar"' not in html:
            continue  # not rendered by mdBook with the theme (for example redirects)
        if "</main>" not in html:
            print(f"error: {rel} has the site bar but no </main> to put the footer after", file=sys.stderr)
            return -1
        html = html.replace("</main>", "</main>\n" + footer.replace(ROOT_TOKEN, root_of(rel)), 1)
        page.write_text(mark_current(html, rel), encoding="utf-8")
        count += 1
    return count


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
    (out / "index.html").write_text(chrome(html, "index.html"), encoding="utf-8")
    shutil.copyfile(SITE / "style.css", out / "landing.css")
    for name in ("tokens.css", "chrome.css"):
        shutil.copyfile(SITE / name, out / name)
        if (out / CALCULATOR_DIR).is_dir():
            shutil.copyfile(SITE / name, out / CALCULATOR_DIR / name)
    for name in TOOL_PAGES:
        page = out / CALCULATOR_DIR / name
        if page.is_file():
            page.write_text(
                chrome(page.read_text(encoding="utf-8"), f"{CALCULATOR_DIR}/{name}"), encoding="utf-8"
            )
    pages = finish_book_pages(out)
    if pages < 0:
        return 1
    print(f"site bar and footer on {pages} book pages")
    print(f"landing page written to {out / 'index.html'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
