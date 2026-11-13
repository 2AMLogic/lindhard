"""The checks every level-3 dataset obeys, and the loop that applies them.

The level-3 experiment scripts (`validation/experiments/run.py`,
`backscatter.py`, `se_yield.py`) each read a directory of JSON datasets, run a
kind-specific checker over every file, print each problem and exit 2 if any
file failed. The parts that do not depend on the kind live here: the three
rules of `validation/data/README.md` that hold for every dataset, and the
loading loop itself.

The third rule is the project rule that **every dataset has a row in
`docs/data-provenance.md`** (`CLAUDE.md`); keeping it in one place means a
change to it is made once rather than once per experiment.

Standard library only, like the rest of the harness.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Callable, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lindhard_cli  # noqa: E402

REPO = lindhard_cli.REPO
PROVENANCE = REPO / "docs" / "data-provenance.md"


def is_number(x) -> bool:
    """Whether `x` is a finite number. `True` and `False` are not numbers here."""
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def check_common(path: Path, d: dict, required_text: list[str], provenance: str) -> list[str]:
    """The problems with `d` that do not depend on its kind, as messages.

    Three rules, in the wording the experiment scripts have always used:
    every key of `required_text` is a non-empty string, `id` is the file name
    without its suffix, and `provenance` (the text of
    `docs/data-provenance.md`) names the `id`.
    """
    errs = []
    for k in required_text:
        if not isinstance(d.get(k), str) or not d[k].strip():
            errs.append(f"`{k}` missing or empty")
    if d.get("id") != path.stem:
        errs.append(f"`id` {d.get('id')!r} differs from the file name {path.stem!r}")
    if isinstance(d.get("id"), str) and d["id"] not in provenance:
        errs.append(f"no row in docs/data-provenance.md names `{d['id']}`")
    return errs


def load_checked(
    files: Iterable[Path],
    check: Callable[[Path, dict, str], list[str]],
    label: str,
) -> list[tuple[Path, dict]]:
    """Every dataset of `files`, read and checked, as `(path, dict)` pairs.

    `check` is the kind's checker, called as `check(path, d, provenance)` with
    the text of `docs/data-provenance.md`; it returns one message per problem.
    Every problem is printed as `error: <path>: <message>`, and a file that is
    not valid JSON is a problem of its own. If any file failed, prints
    `error: <label> check failed (rules: validation/data/README.md)` and exits
    2; otherwise returns the pairs in the order `files` gave them.
    """
    provenance = PROVENANCE.read_text()
    out, failed = [], False
    for f in files:
        rel = f.relative_to(REPO) if f.is_relative_to(REPO) else f
        try:
            d = json.loads(f.read_text())
        except json.JSONDecodeError as e:
            print(f"error: {rel}: not valid JSON: {e}", file=sys.stderr)
            failed = True
            continue
        errs = check(f, d, provenance)
        for e in errs:
            print(f"error: {rel}: {e}", file=sys.stderr)
        failed |= bool(errs)
        out.append((f, d))
    if failed:
        print(f"error: {label} check failed (rules: validation/data/README.md)", file=sys.stderr)
        sys.exit(2)
    return out
