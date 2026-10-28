# lindhard-wasm

WebAssembly front end for [lindhard](../README.md): a single-thread ion range
calculator (`calculate(ion, energy_ev, target, ions, seed)` returning `Rp`,
`ΔRp` and a depth histogram) and the static page in [`web/`](web/index.html)
that drives it. The library crate stays front-end agnostic; this crate is a
thin layer over it (`wasm-bindgen` is its only extra dependency) and is never
published to crates.io. The page ships with the book (see
[the user-guide chapter](../book/src/guide/range-calculator.md)).

## Build

Needs the `wasm32-unknown-unknown` target and the `wasm-bindgen` command-line
tool at the exact version pinned in `Cargo.toml` (release `0.2.129`;
CI downloads it with a checksum, see `.github/workflows/wasm.yml`).

```sh
rustup target add wasm32-unknown-unknown
lindhard-wasm/build.sh web    out/pkg    # ES module for the page
lindhard-wasm/build.sh nodejs out/node   # CommonJS for the Node tests
```

## Test

```sh
cargo test -p lindhard-wasm                                    # native build vs recorded result
LINDHARD_WASM_PKG=out/node node --test lindhard-wasm/tests/node/parity.test.mjs   # WASM build vs the same record
```

Both read `tests/b_5keV_si.golden.json`: 10^4 ions of 5 keV B into Si, seed 1.
The tolerance (one standard error) and its reason are in `tests/parity.rs`.

## Notes

- `lindhard-cli`'s run driver is not used: it times itself with
  `std::time::Instant`, which panics on `wasm32-unknown-unknown`.
- rayon runs on the calling thread there, so results do not depend on the
  thread count (a tested invariant of the library).
