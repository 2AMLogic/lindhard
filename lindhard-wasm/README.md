# lindhard-wasm

WebAssembly front end for [lindhard](../README.md): a single-thread ion range
calculator (`calculate(ion, energy_ev, target, ions, seed)` returning `Rp`,
`ΔRp` and a depth histogram) and the static page in [`web/`](web/index.html)
that drives it. The library crate stays front-end agnostic; this crate is a
thin layer over it (`wasm-bindgen` and `serde_json` are its only extra dependencies) and is never
published to crates.io. The page ships with the book (see
[the user-guide chapter](../book/src/guide/range-calculator.md)).

## Simulation session

`Simulation` is the full engine behind a handle that a Web Worker drives in
small steps (the in-browser simulator). It is plain Rust in
[`src/session.rs`](src/session.rs) with a thin `#[wasm_bindgen]` wrapper in
`src/lib.rs`:

```js
const sim = new Simulation(inputToml, recordHistories, maxPoints); // throws Error on bad input
while (!sim.done()) {
  sim.step(200);                       // returns the histories done so far
  const report = sim.report_json();    // IonReport of everything run so far
  const tracks = sim.take_trajectories_json(); // recorded histories since the last call
}
sim.warnings();                        // the resolver's warnings (string array)
elements_json();                       // the element table, from lindhard::elements
```

- `inputToml` is a run description as in [`docs/cli.md`](../docs/cli.md):
  compound and multi-layer targets, `tilt_deg`/`azimuth_deg`, beam divergence,
  `physics.follow_recoils`, per-element energy overrides, the `[tally]` bins.
  `run.ions` is the target count, at most 1,000,000 per session.
- `[[crystal]]` targets and `[dynamic]` composition are rejected with an error:
  the session transports an amorphous, static target.
- `step(count)` runs at least `count` more histories, **rounded up to a whole
  engine chunk** (`chunk_size`, 64 by default; the last chunk may be shorter)
  and returns the total done. This is what makes the result independent of the
  step sizes: the tally's floating-point moments are merged chunk by chunk in
  index order, so `report_json()` and the trajectories are byte-identical for
  any sequence of step sizes (tested natively and under Node).
- Trajectories are recorded for history indices below `recordHistories`, at
  most `maxPoints` points per history. `take_trajectories_json()` returns
  `{"histories": [{"index", "truncated", "tracks": [{"species", "generation",
  "parent", "end", "points": [[x, y, z, energy_ev], ...]}]}]}` (positions in
  m) and forgets what it returned.
- `report_json()` is `IonTally::report` serialized with serde; the dual
  Pearson fit is made only if the input sets `tally.dual_pearson = true`.
- `elements_json()` lists `symbol`, `z`, `atomic_weight`,
  `weight_is_mass_number`, `density_g_cm3`, `default_displacement_ev` and
  `default_surface_binding_ev` (`null` where the table records none). The
  table has no element names, isotope masses or default binding energies.

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
LINDHARD_WASM_PKG=out/node node --test lindhard-wasm/tests/node/   # WASM build vs the same records
```

Both read `tests/b_5keV_si.golden.json`: 10^4 ions of 5 keV B into Si, seed 1.
The tolerance (one standard error) and its reason are in `tests/parity.rs`.
The session has the same pair: `tests/session.rs` (which also holds the
chunk-invariance and error tests) and `tests/node/session.test.mjs` check
`Simulation` against `tests/as_50keV_si_sio2.golden.json`, 1000 ions of the
`examples/as_50keV_si_sio2.toml` compound two-layer target with recoils
followed, with the same tolerance. CI runs both Node files and prints the size
of the browser `.wasm` (raw and gzipped) in the log of `.github/workflows/wasm.yml`.

## Notes

- `lindhard-cli`'s run driver is not used: it times itself with
  `std::time::Instant`, which panics on `wasm32-unknown-unknown`.
- rayon runs on the calling thread there, so results do not depend on the
  thread count (a tested invariant of the library).
