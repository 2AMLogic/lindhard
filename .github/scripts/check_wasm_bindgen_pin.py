#!/usr/bin/env python3
"""Assert the workflows' wasm-bindgen CLI version equals the crate's `=` pin.

The wasm-bindgen command-line tool refuses a module built by a different
version, and Dependabot bumps only the `=` pin in lindhard-wasm/Cargo.toml —
the WASM_BINDGEN_VERSION env pins in .github/workflows/wasm.yml and
.github/workflows/book.yml must be bumped by hand or every WASM build fails
40 s in with an opaque "rust Wasm file schema version" error (issue #868).
This check turns that drift into an immediate failure naming every value
and file, in either direction.

Python standard library only; regex parse, no TOML/YAML dependency.
Run: python3 .github/scripts/check_wasm_bindgen_pin.py
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")

CRATE_TOML = os.path.join(ROOT, "lindhard-wasm", "Cargo.toml")
WORKFLOWS = ("wasm.yml", "book.yml")

# A `wasm-bindgen = "=..."` dependency line. Anchored so that the key must
# start the line: `wasm-bindgen-futures = ...` is a different crate and the
# prose comments never carry `= "..."` right after the name.
PIN_RE = re.compile(r'^[ \t]*wasm-bindgen[ \t]*=[ \t]*"([^"]*)"[ \t]*$', re.MULTILINE)
# The env definition, not the `${{ env.WASM_BINDGEN_VERSION }}` /
# `$WASM_BINDGEN_VERSION` uses in run steps (no colon after the name there).
VERSION_RE = re.compile(r"^[ \t]*WASM_BINDGEN_VERSION:[ \t]*(.*?)[ \t]*$", re.MULTILINE)
SEMVER_RE = re.compile(r"\d+\.\d+\.\d+$")


def parse_crate_pin(text):
    """The pinned release ("0.2.129") from Cargo.toml text, or None.

    Only the exact-pin form `=0.2.<n>` counts: anything else (a range, no
    dependency line, two contradicting pins) is unparseable for this gate's
    purpose and must fail loudly, not silently pass.
    """
    specs = PIN_RE.findall(text)
    if len(set(specs)) != 1:
        return None
    match = re.fullmatch(r"=(\d+\.\d+\.\d+)", specs[0])
    return match.group(1) if match else None


def parse_workflow_version(text):
    """The WASM_BINDGEN_VERSION release from workflow text, or None."""
    raws = VERSION_RE.findall(text)
    if len(set(raws)) != 1:
        return None
    raw = raws[0].strip()
    if raw and raw[0] in "\"'":
        end = raw.find(raw[0], 1)
        if end == -1:
            return None
        version = raw[1:end]
    else:
        version = raw.split("#", 1)[0].strip()
    if not SEMVER_RE.fullmatch(version):
        return None
    return version


def check(crate_text, workflow_texts):
    """Compare the pins; None on agreement, else an error naming every file."""
    if not workflow_texts:
        return "check_wasm_bindgen_pin: no workflow files given"
    crate = parse_crate_pin(crate_text)
    if crate is None:
        return ("lindhard-wasm/Cargo.toml: cannot parse the wasm-bindgen `=` pin"
                "\n  expected exactly one `wasm-bindgen = \"=<major>.<minor>.<patch>\"`"
                " dependency line")
    versions = {}
    for name, text in workflow_texts.items():
        version = parse_workflow_version(text)
        if version is None:
            return (f"{name}: cannot parse WASM_BINDGEN_VERSION"
                    "\n  expected exactly one"
                    " `WASM_BINDGEN_VERSION: \"<major>.<minor>.<patch>\"` env value")
        versions[name] = version
    if all(version == crate for version in versions.values()):
        return None
    lines = ["wasm-bindgen pin mismatch: the workflows' CLI version must"
             " equal the crate pin"]
    for name in sorted(versions):
        lines.append(f"  {name}: {versions[name]}")
    lines.append(f"  lindhard-wasm/Cargo.toml: ={crate}")
    return "\n".join(lines)


def collect(root=ROOT):
    """Read the pinned files from a checkout; (crate text, {name: text})."""
    with open(os.path.join(root, "lindhard-wasm", "Cargo.toml")) as f:
        crate_text = f.read()
    workflows = {}
    for name in WORKFLOWS:
        with open(os.path.join(root, ".github", "workflows", name)) as f:
            workflows[name] = f.read()
    return crate_text, workflows


def main():
    crate_text, workflows = collect()
    error = check(crate_text, workflows)
    if error is None:
        print("wasm-bindgen pin consistent across wasm.yml, book.yml and"
              f" lindhard-wasm/Cargo.toml: {parse_crate_pin(crate_text)}")
        return 0
    print("::error::" + error.splitlines()[0], file=sys.stderr)
    print(error, file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
