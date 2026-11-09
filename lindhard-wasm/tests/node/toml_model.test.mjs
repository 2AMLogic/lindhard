// Tests of the simulator page's form model (web/toml_model.js): the setup to
// TOML to setup round trip, rejection of keys the form cannot represent, and
// the check that the checked-in file tests/page_generated.toml (which
// tests/cli_equivalence.rs runs natively) is what the page generates.
//
//   node --test lindhard-wasm/tests/node/toml_model.test.mjs

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import {
  ModelError,
  defaultSetup,
  fromToml,
  parseToml,
  toToml,
  validate,
} from "../../web/toml_model.js";

const read = (rel) => readFileSync(new URL(rel, import.meta.url), "utf8");

// The rows of elements_json() the tests need.
const row = (symbol, z, density, eD, eS) => ({
  symbol,
  z,
  density_g_cm3: density,
  default_displacement_ev: eD,
  default_surface_binding_ev: eS,
});
const ELEMENTS = [
  row("B", 5, 2.34, null, null),
  row("O", 8, null, null, null),
  row("Si", 14, 2.33, 15, 4.63),
  row("Ar", 18, null, null, null),
  row("Cu", 29, 8.96, 19, 3.49),
  row("As", 33, 5.73, null, null),
];

// The setup of tests/page_generated.toml: As into 10 nm SiO2 on Si, the case
// of examples/as_50keV_si_sio2.toml, with the cascade followed.
export function pageSetup() {
  const s = defaultSetup();
  s.beam = { ion: "As", massAmu: null, energyEv: 50000, tiltDeg: 7 };
  s.layers = [
    {
      name: "SiO2",
      thicknessNm: 10,
      densityGcm3: 2.2,
      elements: [
        { symbol: "Si", fraction: 1, eDEv: 15, eBEv: null, eSEv: null },
        { symbol: "O", fraction: 2, eDEv: 20, eBEv: null, eSEv: 2 },
      ],
    },
    {
      name: "Si",
      thicknessNm: null,
      densityGcm3: null,
      elements: [{ symbol: "Si", fraction: 1, eDEv: 15, eBEv: null, eSEv: null }],
    },
  ];
  s.physics.followRecoils = true;
  s.run = { ions: 200, seed: 1 };
  s.tally = { depthBinNm: 1, depthBins: 150 };
  return s;
}

test("the page-generated TOML file is what the model writes", () => {
  assert.equal(toToml(pageSetup()), read("../page_generated.toml"));
});

test("setup to TOML to setup is the identity", () => {
  for (const s of [defaultSetup(), pageSetup()]) {
    assert.deepEqual(fromToml(toToml(s)), s);
  }
  const odd = pageSetup();
  odd.beam.massAmu = 10.0129;
  odd.beam.tiltDeg = 60;
  odd.layers[0].elements[0].eBEv = 0.5;
  odd.physics = { potential: "kr-c", stopping: "bethe-bloch", followRecoils: false, primaryCutoffEv: 2.5, recoilCutoffEv: 1 };
  odd.tally = { depthBinNm: 0.25, depthBins: 40 };
  assert.deepEqual(fromToml(toToml(odd)), odd);
});

test("a hand-written input file loads (the shipped example)", () => {
  const s = fromToml(read("../../../examples/as_50keV_si_sio2.toml"));
  assert.equal(s.layers.length, 2);
  assert.equal(s.layers[0].thicknessNm, 10);
  assert.equal(s.layers[1].thicknessNm, null);
  // [physics.energies.Si] is applied to the Si of both layers.
  assert.equal(s.layers[0].elements[0].eDEv, 15);
  assert.equal(s.layers[1].elements[0].eDEv, 15);
  assert.equal(s.layers[0].elements[1].eSEv, 2);
  // An element symbol as a material is the pure element.
  assert.equal(s.layers[1].densityGcm3, null);
  // The Cu example sets tally keys the form does not expose: rejected by name.
  const cuText = read("../../../examples/ar_1keV_cu.toml");
  assert.throws(() => fromToml(cuText), /escape_energy_max_ev/);
  const cu = fromToml(cuText.slice(0, cuText.indexOf("# Escape spectra")));
  assert.equal(cu.layers[0].elements[0].symbol, "Cu");
  assert.equal(cu.physics.potential, "kr-c");
});

test("keys the form cannot represent are rejected, not dropped", () => {
  const base = toToml(defaultSetup());
  const cases = [
    [base.replace("[beam]\n", "[beam]\nazimuth_deg = 10.0\n"), /azimuth_deg/],
    [base.replace("[run]\n", "[run]\nthreads = 2\n"), /threads/],
    [base + "\n[dynamic]\nfluence_cm2 = 1e16\nions_per_step = 10\n", /dynamic/],
    [base + "\n[[crystal]]\nlayers = [0]\n", /crystal/],
    [base.replace("[tally]\n", "[tally]\ndual_pearson = true\n"), /dual_pearson/],
    [base.replace('free_path = "constant"', 'free_path = "energy-dependent"'), /free_path/],
    [base.replace("atom_fraction", "mass_fraction"), /mass_fraction/],
  ];
  for (const [text, want] of cases) {
    assert.throws(() => fromToml(text), (e) => e instanceof ModelError && want.test(e.message), String(want));
  }
});

test("a physics energy for an element not in the target is an error", () => {
  const text = toToml(defaultSetup()) + "\n[physics.energies.Cu]\ne_d_ev = 19.0\n";
  assert.throws(() => fromToml(text), /applies to no element/);
});

test("the TOML reader handles comments, dotted keys, multi-line arrays and quotes", () => {
  const t = parseToml(`
# comment
a.b = 1.5e3 # trailing
[x."y z"]
list = [
  1, 2, # two
  3,
]
inline = { p = "q\\n", r = [true, false] }
[[arr]]
k = 1
[[arr]]
k = 2
`);
  assert.equal(t.a.b, 1500);
  assert.deepEqual(t.x["y z"].list, [1, 2, 3]);
  assert.deepEqual(t.x["y z"].inline, { p: "q\n", r: [true, false] });
  assert.deepEqual(t.arr, [{ k: 1 }, { k: 2 }]);
  assert.throws(() => parseToml("a = "), ModelError);
  assert.throws(() => parseToml("a = 1\na = 2"), /duplicate/);
  assert.throws(() => parseToml('a = "open'), /unterminated/);
});

test("validate flags empty values and invalid layers", () => {
  assert.deepEqual(validate(pageSetup(), ELEMENTS), []);
  assert.deepEqual(validate(defaultSetup(), ELEMENTS), []);
  const paths = (s) => validate(s, ELEMENTS).map((p) => p.path);

  const noEd = defaultSetup();
  noEd.beam.ion = "As";
  noEd.layers[0].elements = [{ symbol: "O", fraction: 1, eDEv: null, eBEv: null, eSEv: null }];
  assert.ok(paths(noEd).includes("layers[0].elements[0].eD"));
  assert.ok(paths(noEd).includes("layers[0].density"));

  const compound = pageSetup();
  compound.layers[0].densityGcm3 = null;
  assert.deepEqual(paths(compound), ["layers[0].density"]);

  const zero = pageSetup();
  zero.layers[0].thicknessNm = 0;
  assert.deepEqual(paths(zero), ["layers[0].thickness"]);

  const emptyLayer = defaultSetup();
  emptyLayer.layers[0].elements = [];
  assert.ok(paths(emptyLayer).includes("layers[0].elements"));

  const middle = pageSetup();
  middle.layers.reverse();
  assert.ok(paths(middle).includes("layers[0].thickness"));

  const dup = pageSetup();
  dup.layers[1].name = "SiO2";
  assert.ok(paths(dup).includes("layers[1].name"));

  const tilt = defaultSetup();
  tilt.beam.tiltDeg = 90;
  tilt.beam.energyEv = -1;
  tilt.run.ions = 0;
  assert.deepEqual(paths(tilt).sort(), ["beam.energy", "beam.tilt", "run.ions"]);

  const unknown = defaultSetup();
  unknown.beam.ion = "Xx";
  assert.deepEqual(paths(unknown), ["beam.ion"]);
});
