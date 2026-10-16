#!/usr/bin/env bash
# Build the WebAssembly package of lindhard-wasm.
#
#   lindhard-wasm/build.sh <target> <out-dir>
#
# <target> is a wasm-bindgen target: `web` (ES module for the browser page) or
# `nodejs` (CommonJS for the Node tests). Needs the `wasm32-unknown-unknown`
# Rust target (`rustup target add wasm32-unknown-unknown`) and the
# `wasm-bindgen` command-line tool at the version pinned in
# lindhard-wasm/Cargo.toml (`=0.2.x`); the tool refuses to run on a module
# built with a different version. The README lists the pinned release.
set -euo pipefail

target="${1:?usage: build.sh <web|nodejs> <out-dir>}"
out="${2:?usage: build.sh <web|nodejs> <out-dir>}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"

cd "$root"
cargo build --release --locked --target wasm32-unknown-unknown -p lindhard-wasm
wasm="$(cargo metadata --format-version 1 --no-deps |
  python3 -c 'import json,sys; print(json.load(sys.stdin)["target_directory"])')"
wasm-bindgen --target "$target" --out-dir "$out" \
  "$wasm/wasm32-unknown-unknown/release/lindhard_wasm.wasm"
echo "built $target package in $out"
