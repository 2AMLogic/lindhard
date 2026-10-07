#!/usr/bin/env python3
"""Check that every variant of the model-selecting enums is named in the
physics manual (book/src/models/*.md).

Each enum below is read from its Rust source: the variants are extracted from
the `pub enum <Name> { ... }` block, so a variant that is added or renamed in
the code and not documented makes this check fail. For every variant the
manual must contain the Rust spelling `<Enum>::<Variant>` in backticks, and for
enums that the TOML input exposes (serde `rename_all = "kebab-case"`) also the
TOML spelling `<key> = "<kebab-case>"` in backticks. Both spellings are
required so a reader can go from either one to the page.

The page template (book/src/models/_template.md) does not count.

Usage: validation/check_manual_coverage.py [--list]

--list prints the spellings that are looked for, one per line, and exits 0.
Exit status 1 lists every missing spelling (CI).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUAL = ROOT / "book" / "src" / "models"
TEMPLATE = "_template.md"

# (source file, enum name, TOML key under [physics] or None if the enum is
# not part of the input schema). Add a row when a new model-selecting enum
# lands; the manual page for it then becomes a CI requirement.
ENUMS: list[tuple[str, str, str | None]] = [
    ("lindhard/src/input.rs", "PotentialChoice", "potential"),
    ("lindhard/src/ion/potential.rs", "Screening", None),
    ("lindhard/src/input.rs", "LengthChoice", "screening_length"),
    ("lindhard/src/ion/potential.rs", "ScreeningLength", None),
    ("lindhard/src/input.rs", "StoppingChoice", "stopping"),
    ("lindhard/src/ion/stopping/bethe.rs", "EffectiveCharge", None),
    ("lindhard/src/ion/stopping/straggling.rs", "StragglingModel", None),
    ("lindhard/src/input.rs", "FreePathChoice", "free_path"),
    ("lindhard/src/ion/bca/mod.rs", "MeanFreePath", None),
    ("lindhard/src/ion/bca/mod.rs", "ElectronicLoss", None),
]

VARIANT_RE = re.compile(r"^([A-Z][A-Za-z0-9]*)\s*(?:[,({=]|$)")
RENAME_RE = re.compile(r'#\[serde\([^)]*\brename\s*=\s*"([^"]+)"')


class CoverageError(Exception):
    pass


def kebab(name: str) -> str:
    """serde's kebab-case of a CamelCase variant name (`LenzJensen` ->
    `lenz-jensen`, `EquipartitionLsOr` -> `equipartition-ls-or`)."""
    out = []
    for i, ch in enumerate(name):
        if ch.isupper() and i > 0:
            out.append("-")
        out.append(ch.lower())
    return "".join(out)


def strip_comment(line: str) -> str:
    # Enough for enum bodies: no string literals contain `//` there.
    return line.split("//", 1)[0]


def enum_variants(source: str, enum: str) -> list[tuple[str, str | None]]:
    """Variants of `pub enum <enum>` as (name, explicit serde rename)."""
    m = re.search(rf"^\s*pub enum {enum}\s*\{{", source, re.MULTILINE)
    if m is None:
        raise CoverageError(f"pub enum {enum} not found")
    depth = 1
    variants: list[tuple[str, str | None]] = []
    pending_rename: str | None = None
    for raw in source[m.end() :].splitlines():
        line = strip_comment(raw).strip()
        if depth == 1 and line.startswith("#["):
            r = RENAME_RE.search(line)
            if r:
                pending_rename = r.group(1)
        elif depth == 1 and line:
            v = VARIANT_RE.match(line)
            if v:
                variants.append((v.group(1), pending_rename))
                pending_rename = None
        depth += line.count("{") - line.count("}")
        if depth <= 0:
            break
    else:
        raise CoverageError(f"pub enum {enum}: no closing brace")
    if not variants:
        raise CoverageError(f"pub enum {enum}: no variants found")
    return variants


def kebab_case_enum(source: str, enum: str) -> bool:
    """True if the enum carries `#[serde(rename_all = "kebab-case")]`."""
    m = re.search(
        rf'#\[serde\(rename_all\s*=\s*"kebab-case"\)\]\s*\n(?:\s*#\[[^\n]*\n)*\s*pub enum {enum}\b',
        source,
    )
    return m is not None


def required_spellings() -> list[tuple[str, str]]:
    """(spelling, where it comes from) pairs the manual must contain."""
    req: list[tuple[str, str]] = []
    for path, enum, key in ENUMS:
        src_path = ROOT / path
        try:
            source = src_path.read_text(encoding="utf-8")
        except OSError as e:
            raise CoverageError(f"{path}: {e}") from e
        variants = enum_variants(source, enum)
        if key is not None and not kebab_case_enum(source, enum):
            raise CoverageError(
                f'{path}: {enum} has a TOML key but no rename_all = "kebab-case"; '
                "update ENUMS in validation/check_manual_coverage.py"
            )
        for name, rename in variants:
            where = f"{path}: {enum}::{name}"
            req.append((f"`{enum}::{name}`", where))
            if key is not None:
                toml = rename if rename is not None else kebab(name)
                req.append((f'`{key} = "{toml}"`', where))
    return req


def manual_text() -> str:
    pages = sorted(p for p in MANUAL.glob("**/*.md") if p.name != TEMPLATE)
    if not pages:
        raise CoverageError(f"no manual pages under {MANUAL.relative_to(ROOT)}")
    return "\n".join(p.read_text(encoding="utf-8") for p in pages)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--list", action="store_true", help="print the required spellings")
    args = ap.parse_args(argv)
    try:
        req = required_spellings()
        if args.list:
            for spelling, where in req:
                print(f"{spelling}\t{where}")
            return 0
        text = manual_text()
    except CoverageError as e:
        print(f"check_manual_coverage: {e}", file=sys.stderr)
        return 1
    missing = [(s, w) for s, w in req if s not in text]
    if missing:
        print(
            "check_manual_coverage: the physics manual (book/src/models/) does not "
            "mention these model variants:",
            file=sys.stderr,
        )
        for spelling, where in missing:
            print(f"  {spelling}  ({where})", file=sys.stderr)
        print(
            "Add them to the model's page (see book/src/models/_template.md).",
            file=sys.stderr,
        )
        return 1
    print(f"manual coverage: {len(req)} spellings of {len(ENUMS)} enums, all present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
