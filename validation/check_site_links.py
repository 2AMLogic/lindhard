#!/usr/bin/env python3
"""Check the rendered site: landing page, a book page and the two tools.

Every relative href/src in each of those pages must resolve to a file in the built
output (a directory link needs an index.html), and the book URLs that existed
before the landing page took over the root must still resolve. External links
are not fetched (flaky in CI). Anchors are checked against the ids of the
target page when it is local.

Usage: validation/check_site_links.py [BOOK_DIR]     (default: book/book)
Exit status 1 lists every failure. Python standard library only.
"""

from __future__ import annotations

import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]

# URLs published before the landing page; they must keep resolving.
EXISTING = ["introduction.html", "models/overview.html", "guide/install.html", "guide/first-run.html"]
# Pages the landing page links to that must exist in the built book.
EXISTING += ["validation.html", "status.html", "about.html"]
# The landing page's hero figure, gallery thumbnails and examples.
EXISTING += ["gallery/index.html", "gallery/generated/quickstart.svg"]
# Ids the landing page's sections must keep (anchors other pages may link to).
LANDING_IDS = ["try", "what", "validation", "physics", "examples", "start", "status"]


class Collector(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.refs: list[str] = []
        self.ids: set[str] = set()

    def handle_starttag(self, tag, attrs):
        for k, v in attrs:
            if k in ("href", "src") and v:
                self.refs.append(v)
            if k == "id" and v:
                self.ids.add(v)


def parse(path: Path) -> Collector:
    c = Collector()
    c.feed(path.read_text(encoding="utf-8"))
    return c


# One page per surface: landing, book, range calculator, simulator.
SURFACES = ["index.html", "introduction.html", "calculator/index.html", "calculator/simulator.html"]
# Tool pages whose own CSS and script must hold no colour literal (tokens only).
TOOL_SOURCES = ["lindhard-wasm/web/simulator.html"]
HEX_COLOUR = re.compile(r"(?<![\w&\"'=/-])#(?:[0-9a-fA-F]{8}|[0-9a-fA-F]{6}|[0-9a-fA-F]{3,4})\b(?![-\w])")


def check_page(out: Path, rel: str, bad: list[str]) -> int:
    """Every relative href/src of the page resolves, relative to the page."""
    path = out / rel
    if not path.is_file():
        bad.append(f"{rel} does not exist")
        return 0
    page = parse(path)
    for ref in page.refs:
        u = urlsplit(ref)
        if u.scheme or u.netloc:
            continue
        if not u.path:  # same-page anchor
            if u.fragment and u.fragment not in page.ids:
                bad.append(f"{rel}: {ref}: no such id on the page")
            continue
        if u.path.startswith("/"):
            bad.append(f"{rel}: {ref}: root-relative link breaks under /lindhard/")
            continue
        target = (path.parent / unquote(u.path)).resolve()
        if target.is_dir():
            target = target / "index.html"
        if not target.is_file():
            bad.append(f"{rel}: {ref}: not found")
        elif u.fragment and target.suffix == ".html" and u.fragment not in parse(target).ids:
            bad.append(f"{rel}: {ref}: no id {u.fragment!r} in {target.name}")
    return len(page.refs)


def check_landing(out: Path, bad: list[str]) -> None:
    """The landing page keeps its sections, copy buttons and no unfilled build slot."""
    path = out / "index.html"
    ids = parse(path).ids
    for i in LANDING_IDS:
        if i not in ids:
            bad.append(f"index.html: section id {i!r} missing")
    text = path.read_text(encoding="utf-8")
    if "<!--" in text and re.search(r"<!--[A-Z_]+-->", text):
        bad.append("index.html: an unfilled build slot is left in the page")
    if text.count('class="codeblock"') < 2:
        bad.append("index.html: the build and run blocks are missing")


def check_chrome_sources(bad: list[str]) -> None:
    """The book's header is a byte copy of the shared fragment; tool pages use tokens only."""
    header = (ROOT / "site" / "header.html").read_text(encoding="utf-8")
    if (ROOT / "book" / "theme" / "header.hbs").read_text(encoding="utf-8") != header:
        bad.append("book/theme/header.hbs differs from site/header.html (copy it)")
    for rel in TOOL_SOURCES:
        for n, line in enumerate((ROOT / rel).read_text(encoding="utf-8").splitlines(), 1):
            if HEX_COLOUR.search(line):
                bad.append(f"{rel}:{n}: colour literal; use a --lh- token from site/tokens.css")


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "book" / "book"
    if not (out / "index.html").is_file():
        print(f"error: {out / 'index.html'} does not exist", file=sys.stderr)
        return 1
    bad: list[str] = []
    checked = sum(check_page(out, rel, bad) for rel in SURFACES)
    for url in EXISTING:
        if not (out / url).is_file():
            bad.append(f"existing book URL {url} no longer resolves")
    check_chrome_sources(bad)
    check_landing(out, bad)
    if bad:
        print("site link check failed:", file=sys.stderr)
        for b in bad:
            print(f"  {b}", file=sys.stderr)
        return 1
    print(f"ok: {checked} references checked on {len(SURFACES)} pages in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
