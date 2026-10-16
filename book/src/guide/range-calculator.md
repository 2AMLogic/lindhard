# Range calculator in the browser

The [range calculator](../calculator/index.html) runs the ion transport of
this book in your browser, as WebAssembly on one thread. Give it an ion, an
energy and an elemental target; it returns the projected range
$R_p$, the straggle $\Delta R_p$ and a depth profile. It is a static page:
there is no server, and nothing you enter leaves your machine.

A run of $10^4$ 5 keV B ions into Si takes a few seconds in a current browser.
The calculator uses the choices of `examples/b_5keV_si.toml` (amorphous
target at its tabulated density, ZBL screening, Lindhard–Scharff electronic
stopping, 7° tilt), and a 15 eV displacement energy for a target element the
element table gives none. It is not a replacement for the command line: it has
no layers, compounds, crystals or electron mode.

## Build it yourself

```sh
rustup target add wasm32-unknown-unknown
# wasm-bindgen command-line tool, release 0.2.129 (the version pinned in
# lindhard-wasm/Cargo.toml)
lindhard-wasm/build.sh web out/pkg
cp lindhard-wasm/web/index.html out/index.html
python3 -m http.server -d out        # then open http://localhost:8000
```

## Agreement with the native build

The same input and seed give the same $R_p$ and $\Delta R_p$ as the native
build, to within one standard error (in practice to about $10^{-15}$
relative): `cargo test -p lindhard-wasm` checks the native build and
`lindhard-wasm/tests/node/parity.test.mjs` the WebAssembly build against one
recorded result. The tolerance and its reason are in
`lindhard-wasm/tests/parity.rs`: `wasm32-unknown-unknown` uses Rust's bundled
libm where a native build uses the platform's, the two differ in the last bit,
and a collision cascade amplifies that until individual histories diverge into
statistically equivalent ones.
