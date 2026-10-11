"""Shared fixtures: the example inputs and the `lindhard` command to compare with."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
EXAMPLES = REPO / "examples"


def _target_dir() -> Path:
    env = os.environ.get("CARGO_TARGET_DIR")
    return Path(env) if env else REPO / "target"


@pytest.fixture(scope="session")
def cli() -> Path:
    """Path of the `lindhard` binary: $LINDHARD_BIN, an existing build, or a fresh
    `cargo build -p lindhard-cli`."""
    env = os.environ.get("LINDHARD_BIN")
    if env:
        return Path(env)
    for profile in ("debug", "release"):
        p = _target_dir() / profile / "lindhard"
        if p.exists():
            return p
    if shutil.which("cargo") is None:
        pytest.skip("no lindhard binary and no cargo to build one")
    subprocess.run(
        ["cargo", "build", "-p", "lindhard-cli", "-j", "2"], cwd=REPO, check=True
    )
    return _target_dir() / "debug" / "lindhard"


@pytest.fixture(params=sorted(EXAMPLES.glob("*.toml")), ids=lambda p: p.stem)
def example(request) -> Path:
    return request.param
