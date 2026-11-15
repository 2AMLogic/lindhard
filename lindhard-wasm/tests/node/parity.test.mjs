// The WASM half of the native-versus-WASM parity test, and the smoke test
// that the built package loads and runs the golden fixture.
//
// Run by `lindhard-wasm/build.sh nodejs <dir>` followed by
//   LINDHARD_WASM_PKG=<dir> node --test lindhard-wasm/tests/node/parity.test.mjs
// (the CI workflow .github/workflows/wasm.yml does exactly that).
//
// The tolerance and its reason are in lindhard-wasm/tests/parity.rs, which
// checks the native build against the same golden file: Rp and dRp within one
// standard error of the golden values, the stopped count within the binomial
// scatter.

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { resolve } from "node:path";
import test from "node:test";

const pkg = process.env.LINDHARD_WASM_PKG;
assert.ok(pkg, "set LINDHARD_WASM_PKG to the directory build.sh wrote");
const wasm = createRequire(import.meta.url)(resolve(pkg, "lindhard_wasm.js"));
const read = (rel) => readFileSync(new URL(rel, import.meta.url), "utf8");
const inputToml = read("../b_5keV_si.toml");
const golden = JSON.parse(read("../b_5keV_si.golden.json"));

// Run `toml` to completion and return its final report.
function runAll(toml) {
  const sim = new wasm.Simulation(toml, 0, 10);
  try {
    while (!sim.done()) sim.step(1000);
    return JSON.parse(sim.report_json());
  } finally {
    sim.free();
  }
}

test("the package loads and reports the library version", () => {
  assert.match(wasm.version(), /^\d+\.\d+\.\d+/);
});

test("B into Si matches the native golden record", () => {
  const report = runAll(inputToml);
  const depth = report.range.depth;
  assert.ok(
    Math.abs(depth.mean * 1e9 - golden.rp_nm) <= golden.rp_std_err_nm,
    `Rp ${depth.mean * 1e9} vs ${golden.rp_nm}`,
  );
  assert.ok(
    Math.abs(depth.std_dev * 1e9 - golden.drp_nm) <= golden.drp_std_err_nm,
    `dRp ${depth.std_dev * 1e9} vs ${golden.drp_nm}`,
  );
  const sigma = Math.sqrt(golden.stopped * (1 - golden.stopped / golden.ions));
  assert.ok(Math.abs(report.range.stopped - golden.stopped) <= sigma);
  assert.equal(report.histories, golden.ions);
  // The depth histogram accounts for every stopped primary.
  const h = report.range.depth_histogram;
  const binned = h.counts.reduce((a, b) => a + b, 0) + h.underflow + h.overflow;
  assert.equal(binned, report.range.stopped);
});

test("invalid input throws an Error with the library's message", () => {
  assert.throws(
    () => new wasm.Simulation(inputToml.replace('ion = "B"', 'ion = "Xx"'), 0, 10),
    /unknown element/,
  );
  assert.throws(
    () => new wasm.Simulation(inputToml.replace("ions = 10000", "ions = 2000000"), 0, 10),
    /ions must be/,
  );
  // The module is still usable after the errors.
  assert.equal(runAll(inputToml).histories, golden.ions);
});
