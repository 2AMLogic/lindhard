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
  PRESETS,
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
  s.beam = { ion: "As", massAmu: null, energyEv: 50000, tiltDeg: 7, azimuthDeg: 0 };
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

// The setup of tests/page_crystal.toml: pageSetup with the substrate on the
// Si lattice, 7° tilt and 22° azimuth, 200 ions (what
// tests/cli_equivalence.rs runs natively).
export function pageCrystalSetup() {
  const s = pageSetup();
  s.beam.azimuthDeg = 22;
  s.crystal = { layers: [1], preset: "Si", normal: [0, 0, 1], reference: [0, 1, 0], waferRotationDeg: 0, temperatureK: 300 };
  s.run = { ions: 200, seed: 1 };
  return s;
}

test("the page-generated TOML file is what the model writes", () => {
  assert.equal(toToml(pageSetup()), read("../page_generated.toml"));
});

test("the page-generated crystal TOML file is what the model writes", () => {
  assert.equal(toToml(pageCrystalSetup()), read("../page_crystal.toml"));
});

test("setup to TOML to setup is the identity", () => {
  for (const s of [defaultSetup(), pageSetup(), pageCrystalSetup()]) {
    assert.deepEqual(fromToml(toToml(s)), s);
  }
  const odd = pageSetup();
  odd.beam.massAmu = 10.0129;
  odd.beam.tiltDeg = 60;
  odd.beam.azimuthDeg = 35;
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
  const crystal = (body) => base + `\n[[crystal]]\n${body}\n`;
  const cases = [
    [base.replace("[run]\n", "[run]\nthreads = 2\n"), /threads/],
    [base + "\n[dynamic]\nfluence_cm2 = 1e16\nions_per_step = 10\n", /dynamic/],
    [
      crystal('layers = [0]\npreset = "Si"\nnormal = [0, 0, 1]\nreference = [0, 1, 0]\np_max_nm = 0.2\n'),
      /p_max_nm/,
    ],
    [crystal('layers = [0]\npreset = "Si"\n[crystal.thermal]\ntemperature_k = 300.0\ndebye_temperature_k = 640.0\n'), /debye_temperature_k/],
    [crystal('layers = [0]\npreset = "Si"\n\n[[crystal]]\nlayers = [0]\npreset = "Ge"\n'), /2 \[\[crystal\]\]/],
    [base.replace("[tally]\n", "[tally]\ndual_pearson = true\n"), /dual_pearson/],
    [base.replace('free_path = "constant"', 'free_path = "energy-dependent"'), /free_path/],
    [base.replace("atom_fraction", "mass_fraction"), /mass_fraction/],
  ];
  for (const [text, want] of cases) {
    assert.throws(() => fromToml(text), (e) => e instanceof ModelError && want.test(e.message), String(want));
  }
});

// The body of a [[crystal]] table the form accepts, on defaultSetup's one-layer
// Si substrate (layer index 0).
const crystalBody = 'layers = [0]\npreset = "Si"\nnormal = [0, 0, 1]\nreference = [0, 1, 0]';
const crystalFile = (body) => toToml(defaultSetup()) + `\n[[crystal]]\n${body}\n`;

test("the shipped crystal example loads and round-trips", () => {
  const s = fromToml(read("../../../examples/b_5keV_si_crystal.toml"));
  // The beam azimuth is the crystal's twist: the form now holds it.
  assert.equal(s.beam.azimuthDeg, 22);
  assert.deepEqual(s.crystal.layers, [0]);
  assert.equal(s.crystal.preset, "Si");
  assert.deepEqual(s.crystal.normal, [0, 0, 1]);
  assert.deepEqual(s.crystal.reference, [0, 1, 0]);
  assert.equal(s.crystal.waferRotationDeg, 0);
  assert.equal(s.crystal.temperatureK, 300);
  assert.deepEqual(fromToml(toToml(s)), s);
  // A crystal-free file still gives crystal: null.
  assert.equal(fromToml(toToml(pageSetup())).crystal, null);
});

test("a crystal setup round-trips through the model", () => {
  const s = defaultSetup();
  s.beam.azimuthDeg = 22;
  s.crystal = { layers: [0], preset: "3C-SiC", normal: [0, 0, 1], reference: [0, 1, 0], waferRotationDeg: 15, temperatureK: null };
  assert.deepEqual(fromToml(toToml(s)), s);
  s.crystal.temperatureK = 300;
  assert.deepEqual(fromToml(toToml(s)), s);
  assert.ok(toToml(s).includes("[crystal.thermal]"));
  s.crystal.temperatureK = null;
  assert.ok(!toToml(s).includes("[crystal.thermal]"));
  // A zero azimuth is the engine default (the key is absent, as in every
  // amorphous input): not emitted.
  assert.ok(!toToml(pageSetup()).includes("azimuth_deg"));
});

test("an orientation the engine would reject is refused at load, with the engine's meaning", () => {
  // Zone law: [1,1,1] is not in the (001) surface plane (0 + 0 + 1 = 1).
  assert.throws(
    () => fromToml(crystalFile(crystalBody.replace("reference = [0, 1, 0]", "reference = [1, 1, 1]"))),
    (e) => e instanceof ModelError && /reference/.test(e.message) && /h u \+ k v \+ l w = 1, must be 0/.test(e.message),
  );
  assert.throws(
    () => fromToml(crystalFile(crystalBody.replace("normal = [0, 0, 1]", "normal = [0, 0, 0]"))),
    (e) => e instanceof ModelError && /normal/.test(e.message) && /all zero/.test(e.message),
  );
  assert.throws(
    () => fromToml(crystalFile(crystalBody.replace("reference = [0, 1, 0]", "reference = [0, 0, 0]"))),
    (e) => e instanceof ModelError && /reference/.test(e.message) && /all zero/.test(e.message),
  );
  // Layer indices: empty, out of range, assigned twice.
  assert.throws(() => fromToml(crystalFile(crystalBody.replace("layers = [0]", "layers = []"))), (e) => e instanceof ModelError && /layers.*at least one/.test(e.message));
  assert.throws(
    () => fromToml(crystalFile(crystalBody.replace("layers = [0]", "layers = [1]"))),
    (e) => e instanceof ModelError && /does not exist/.test(e.message),
  );
  assert.throws(
    () => fromToml(crystalFile(crystalBody.replace("layers = [0]", "layers = [0, 0]"))),
    (e) => e instanceof ModelError && /already assigned/.test(e.message),
  );
  // Non-integer index triples are refused too.
  assert.throws(
    () => fromToml(crystalFile(crystalBody.replace("normal = [0, 0, 1]", "normal = [0, 0, 1.5]"))),
    (e) => e instanceof ModelError && /normal.*three integers/.test(e.message),
  );
});

test("validate flags crystal problems", () => {
  const s = defaultSetup();
  s.beam.azimuthDeg = 22;
  s.crystal = { layers: [0], preset: "Si", normal: [0, 0, 1], reference: [0, 1, 0], waferRotationDeg: 0, temperatureK: 300 };
  assert.deepEqual(validate(s, ELEMENTS), []);
  const paths = (x) => validate(x, ELEMENTS).map((p) => p.path);

  const zone = { ...s, crystal: { ...s.crystal, reference: [1, 1, 1] } };
  assert.deepEqual(paths(zone), ["crystal.reference"]);

  const empty = { ...s, crystal: { ...s.crystal, layers: [] } };
  assert.deepEqual(paths(empty), ["crystal.layers"]);

  const outOfRange = { ...s, layers: [{ ...s.layers[0], name: "Si-1", thicknessNm: 100 }, { ...s.layers[0], thicknessNm: null }], crystal: { ...s.crystal, layers: [2] } };
  assert.deepEqual(paths(outOfRange), ["crystal.layers"]);

  const preset = { ...s, crystal: { ...s.crystal, preset: "NaCl" } };
  assert.deepEqual(paths(preset), ["crystal.preset"]);

  const temp = { ...s, crystal: { ...s.crystal, temperatureK: -5 } };
  assert.deepEqual(paths(temp), ["crystal.temperature"]);

  const dup = { ...s, crystal: { ...s.crystal, layers: [0, 0] } };
  assert.deepEqual(paths(dup), ["crystal.layers"]);

  const az = { ...s, beam: { ...s.beam, azimuthDeg: NaN } };
  assert.deepEqual(paths(az), ["beam.azimuth"]);
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

test("Si and SiO2 presets validate without a displacement-energy error; GaAs asks for Ga and As", () => {
  // Element rows as the page gets them: no tabulated E_d for Si, O, Ga, As.
  const table = ["Si", "O", "Ga", "As", "B", "Cu", "Au"].map((symbol) => ({
    symbol, default_displacement_ev: null, default_surface_binding_ev: 3, density_gcm3: 2,
  }));
  const eDErrors = (key) => {
    const s = defaultSetup();
    s.layers = [{ name: key, thicknessNm: null, densityGcm3: PRESETS[key].density ?? 2,
      elements: PRESETS[key].elements.map((e) => ({ ...e })) }];
    return validate(s, table).filter((b) => b.path.endsWith(".eD"));
  };
  for (const key of ["Si", "SiO2"]) assert.deepEqual(eDErrors(key), [], key);
  assert.equal(PRESETS.Si.elements[0].eDEv, 15);
  assert.equal(PRESETS.SiO2.elements[0].eDEv, 15);
  assert.equal(eDErrors("GaAs").length, 2);
  assert.match(PRESETS.GaAs.note, /enter one/);
});
