# Range calculator in the browser

The [range calculator](../calculator/index.html) runs the ion transport of
this book in your browser, as WebAssembly in a Web Worker. Give it an ion, an
energy and an elemental target; it returns the projected range
$R_p$, the straggle $\Delta R_p$ and a depth profile. It is a static page:
there is no server, and nothing you enter leaves your machine.

The calculator is the [simulator](simulator.md) in its Quick range state: the
same page, engine and plot, with only the ion, energy and target shown and the
other settings under "More options". Like the simulator it opens with an
example run of 1000 ions with seed 1 (10 keV B into Si), labelled as such.
A run of $10^4$ ions takes a few seconds in a current browser; the page runs at
most $10^5$ ions. The depth profile is plotted as counts per depth bin.
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
cp lindhard-wasm/web/*.js lindhard-wasm/web/simulator.html out/
# The Quick range page is the simulator page with data-view="quick" on the body.
sed 's|^<body>$|<body data-view="quick">|' lindhard-wasm/web/simulator.html > out/index.html
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
