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

import html as htmlmod
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"

# The one place the simulator page's path lives.
SIMULATOR_PAGE = "calculator/simulator.html"
API_PAGE = "api/lindhard/index.html"

# Landing-page sources that are generated or published elsewhere.
VALIDATION_MD = ROOT / "book" / "src" / "validation.md"
CARGO_TOML = ROOT / "Cargo.toml"
QUICKSTART_SVG = "gallery/generated/quickstart.svg"
QUICKSTART_MD = ROOT / "book" / "src" / "gallery" / "generated" / "quickstart.md"
# Rows of the level 1 table and bullets of "Known deviations" the landing page
# quotes, chosen by the start of their text. A label that is not found fails
# the build, so rewording a row in validation.md cannot silently drop it.
TRUST_ROWS = (
    "Quadrature scattering angle",
    "Projected range, 50 keV P into Si",
    "Lindhard partition solver",
    "Per-history energy budget",
)
DEVIATIONS = (
    "ZBL universal fit",
    "Range against measurement",
    "Dynamic composition",
)

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


def inline(text: str) -> str:
    """Escape text and turn the **bold** and `code` of the book's Markdown into HTML."""
    text = htmlmod.escape(text, quote=False)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", text)


def workspace_version() -> str:
    m = re.search(r"^\[workspace\.package\]\s*$.*?^version\s*=\s*\"([^\"]+)\"", CARGO_TOML.read_text(encoding="utf-8"), re.S | re.M)
    if not m:
        raise ValueError("no [workspace.package] version in Cargo.toml")
    return m.group(1)


def section(md: str, heading: str) -> str:
    m = re.search(rf"^## {re.escape(heading)}[^\n]*\n(.*?)(?=^## |\Z)", md, re.S | re.M)
    if not m:
        raise ValueError(f"validation.md has no section '{heading}'")
    return m.group(1)


def trust_rows(md: str) -> str:
    """HTML rows (compared / against / deviation) picked by label from the level 1 table."""
    table = [
        [c.strip() for c in line.strip().strip("|").split("|")]
        for line in section(md, "Level 1").splitlines()
        if line.startswith("|")
    ][2:]  # drop the header and the rule
    out = []
    for label in TRUST_ROWS:
        row = next((r for r in table if r[0].startswith(label)), None)
        if row is None or len(row) < 3:
            raise ValueError(f"validation.md level 1 table has no row starting '{label}'")
        out.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in row[:3]) + "</tr>")
    return "\n".join(out)


def known_deviations(md: str) -> str:
    """One list item per chosen deviation: its label and first sentence."""
    items = re.findall(r"^- (.*?)(?=^- |\Z)", section(md, "Known deviations"), re.S | re.M)
    out = []
    for label in DEVIATIONS:
        raw = next((" ".join(i.split()) for i in items if i.startswith(f"**{label}")), None)
        if raw is None:
            raise ValueError(f"validation.md Known deviations has no item starting '**{label}'")
        first = re.split(r"(?<=[.])\s", raw, maxsplit=1)[0]
        out.append(f"<li>{inline(first)}</li>")
    return "\n".join(out)


def hero_figure(out: Path) -> str:
    """The quick-start figure and its generated caption, unmodified."""
    if not (out / QUICKSTART_SVG).is_file():
        raise ValueError(f"{out} lacks {QUICKSTART_SVG} (run validation/gallery.py)")
    caption = QUICKSTART_MD.read_text(encoding="utf-8").strip()
    caption = re.sub(r"`([^`]+)`", r"<code>\1</code>", caption)  # the caption's only Markdown
    return (
        '<figure class="hero-figure">'
        f'<a href="{SIMULATOR_PAGE}"><img src="{QUICKSTART_SVG}" width="640" height="530" '
        'alt="Paths of the first ions above the depth profile of 5 keV boron in silicon"></a>'
        f"<figcaption>{caption}</figcaption></figure>"
    )


def favicon(out: Path) -> str:
    """Links to the favicons mdBook wrote for the book, for those that exist."""
    links = []
    for name, kind in (("favicon.svg", "image/svg+xml"), ("favicon.png", "image/png")):
        if (out / name).is_file():
            links.append(f'<link rel="icon" href="{name}" type="{kind}">')
    return "\n".join(links)


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "book" / "book"
    if not (out / "introduction.html").is_file():
        print(f"error: {out} is not a rendered book (no introduction.html)", file=sys.stderr)
        return 1
    html = (SITE / "index.html").read_text(encoding="utf-8")
    if not (out / SIMULATOR_PAGE).is_file():
        print(f"error: {out} lacks the simulator page {SIMULATOR_PAGE}", file=sys.stderr)
        return 1
    sim = f'<a href="{SIMULATOR_PAGE}">Trajectory simulator</a>'
    if (out / API_PAGE).is_file():
        api = f'<a href="{API_PAGE}">API documentation</a> (rustdoc)'
    else:
        api = "API documentation: <code>cargo doc -p lindhard --open</code>"
    md = VALIDATION_MD.read_text(encoding="utf-8")
    try:
        generated = (
            ("<!--VERSION-->", f"<code>{htmlmod.escape(workspace_version())}</code>"),
            ("<!--TRUST_ROWS-->", trust_rows(md)),
            ("<!--KNOWN_DEVIATIONS-->", known_deviations(md)),
            ("<!--HERO_FIGURE-->", hero_figure(out)),
            ("<!--FAVICON-->", favicon(out)),
        )
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    for marker, text in generated + (("<!--SIMULATOR_SLOT-->", sim), ("<!--API_LINK-->", api)):
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
