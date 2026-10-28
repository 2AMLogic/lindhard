#!/usr/bin/env python3
"""Check the landing page of the rendered book.

Every relative href/src in index.html must resolve to a file in the built
output (a directory link needs an index.html), and the book URLs that existed
before the landing page took over the root must still resolve. External links
are not fetched (flaky in CI). Anchors are checked against the ids of the
target page when it is local.

Usage: validation/check_site_links.py [BOOK_DIR]     (default: book/book)
Exit status 1 lists every failure. Python standard library only.
"""

from __future__ import annotations

import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]

# URLs published before the landing page; they must keep resolving.
EXISTING = ["introduction.html", "models/overview.html", "guide/install.html", "guide/first-run.html"]


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


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "book" / "book"
    index = out / "index.html"
    if not index.is_file():
        print(f"error: {index} does not exist", file=sys.stderr)
        return 1
    page = parse(index)
    bad: list[str] = []
    for ref in page.refs:
        u = urlsplit(ref)
        if u.scheme or u.netloc:
            continue
        if not u.path:  # same-page anchor
            if u.fragment and u.fragment not in page.ids:
                bad.append(f"{ref}: no such id on index.html")
            continue
        if u.path.startswith("/"):
            bad.append(f"{ref}: root-relative link breaks under /lindhard/")
            continue
        target = out / unquote(u.path)
        if target.is_dir():
            target = target / "index.html"
        if not target.is_file():
            bad.append(f"{ref}: {target.relative_to(out)} not found")
        elif u.fragment and target.suffix == ".html" and u.fragment not in parse(target).ids:
            bad.append(f"{ref}: no id {u.fragment!r} in {target.relative_to(out)}")
    for url in EXISTING:
        if not (out / url).is_file():
            bad.append(f"existing book URL {url} no longer resolves")
    if bad:
        print("landing page link check failed:", file=sys.stderr)
        for b in bad:
            print(f"  {b}", file=sys.stderr)
        return 1
    print(f"ok: {len(page.refs)} references checked in {index}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
