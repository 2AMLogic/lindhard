"""The oracle-harness scaffolding the two level-2 runners share
(`validation/oracles/run.py` and `validation/oracles/run_electron.py`).

An *oracle* is a third-party program the user installs outside this tree and
that the harness runs unmodified; it is never vendored, and nothing in CI
needs one (CONTRIBUTING.md, "Oracles"). Each is configured by environment
variables named after the oracle:

    <PREFIX>_BIN        the executable (the `env` attribute below)
    <PREFIX>_SRC        optional: the git clone it was built from, for the
                        recorded version and commit
    <PREFIX>_VERSION, <PREFIX>_COMMIT, <PREFIX>_BUILD
                        optional overrides of the recorded version, commit
                        and build description

`<PREFIX>` is `env` without its trailing `_BIN`; an oracle whose variables do
not follow that convention sets `prefix` explicitly.

`version()` is `$<PREFIX>_VERSION`, else `git describe` of `$<PREFIX>_SRC`,
else the first line of `<binary> --version`, else `"unknown"`. The last step
assumes the program answers `--version`: an oracle that does not (because it
has no such option, or answers a different one) is expected to override
`version()` — `Nebula` and `OpenTRIM` both do.

Standard library only (Python 3.9 or newer), like the rest of the harness.
File hashing lives next door, in `lindhard_cli.sha256_file`.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


class NotConfigured(Exception):
    """An oracle the user has not installed or configured. The runner prints
    it and carries on: a missing oracle is a skip, never a failure."""


def git_output(src: str | None, *args: str) -> str | None:
    """`git -C src <args>` stripped, or None if `src` is unset, git is
    missing, slow or unhappy, or the output is empty."""
    if not src:
        return None
    try:
        out = subprocess.run(["git", "-C", src, *args], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def round_floats(x):
    """Every float in a JSON-ready structure to six significant digits, so a
    rerun that only changes the last bits leaves the committed summary
    byte-identical."""
    if isinstance(x, float):
        return float(f"{x:.6g}")
    if isinstance(x, dict):
        return {k: round_floats(v) for k, v in x.items()}
    if isinstance(x, list):
        return [round_floats(v) for v in x]
    return x


class Oracle:
    """One external program: where its executable is, and how it identifies
    itself in a summary. Each runner adds its own `run` contract."""

    name = "?"
    env = "?"
    # The environment-variable prefix; None = `env` without its `_BIN`.
    prefix: str | None = None
    source_url = "?"
    license_note = ""

    def _var(self, suffix: str) -> str | None:
        prefix = self.prefix if self.prefix is not None else self.env.removesuffix("_BIN")
        return os.environ.get(prefix + suffix)

    def binary(self) -> Path:
        b = os.environ.get(self.env)
        if not b:
            raise NotConfigured(f"{self.name}: not configured (set {self.env}); skipped")
        p = Path(b)
        if not p.exists():
            raise NotConfigured(f"{self.name}: {self.env}={b} does not exist; skipped")
        self.check(p)
        return p

    def check(self, binary: Path) -> None:
        """Extra prerequisites; raise NotConfigured if one is missing."""

    def version(self, binary: Path) -> str:
        """`$<PREFIX>_VERSION`, else `git describe` of `$<PREFIX>_SRC`, else
        the first line of `--version`. An oracle without `--version` overrides
        this (see the module docstring)."""
        v = self._var("_VERSION") or git_output(self._var("_SRC"), "describe", "--tags", "--always")
        if v:
            return v
        try:
            out = subprocess.run([str(binary), "--version"], capture_output=True, text=True, timeout=30)
            line = (out.stdout or out.stderr).strip().splitlines()
            return line[0] if line else "unknown"
        except (OSError, subprocess.TimeoutExpired):
            return "unknown"

    def commit(self) -> str:
        return self._var("_COMMIT") or git_output(self._var("_SRC"), "rev-parse", "HEAD") or "unknown"

    def build(self) -> str:
        return self._var("_BUILD") or "unknown"

    def applies(self, problem: dict) -> str | None:
        """None if the oracle can run the problem, else the reason it cannot.
        Consulted by the electron runner (`run_electron.py`); the ion runner
        matches every oracle to every problem and never calls it."""
        return None
