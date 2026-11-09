// The WASM half of the native-versus-WASM parity test, and the smoke test
// that the built package loads and returns a profile.
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
const golden = JSON.parse(
  readFileSync(new URL("../b_5keV_si.golden.json", import.meta.url), "utf8"),
);

test("the package loads and reports the library version", () => {
  assert.match(wasm.version(), /^\d+\.\d+\.\d+/);
});

test("B into Si matches the native golden record", () => {
  const r = wasm.calculate(
    golden.ion,
    golden.energy_ev,
    golden.target,
    golden.ions,
    golden.seed,
  );
  assert.ok(
    Math.abs(r.rp_nm - golden.rp_nm) <= golden.rp_std_err_nm,
    `Rp ${r.rp_nm} vs ${golden.rp_nm}`,
  );
  assert.ok(
    Math.abs(r.drp_nm - golden.drp_nm) <= golden.drp_std_err_nm,
    `dRp ${r.drp_nm} vs ${golden.drp_nm}`,
  );
  const sigma = Math.sqrt(golden.stopped * (1 - golden.stopped / golden.ions));
  assert.ok(Math.abs(r.stopped - golden.stopped) <= sigma);
  assert.equal(r.histories, golden.ions);
});

test("the profile is a histogram of the stopped ions", () => {
  const r = wasm.calculate("B", 5000, "Si", 2000, 7);
  const counts = Array.from(r.counts);
  assert.equal(counts.length, 100);
  assert.equal(
    counts.reduce((a, b) => a + b, 0),
    r.stopped,
  );
  assert.ok(r.bin_nm > 0);
});

test("invalid input throws an Error with the library's message", () => {
  assert.throws(() => wasm.calculate("Xx", 5000, "Si", 100, 1), /unknown element/);
  assert.throws(() => wasm.calculate("B", 5000, "Si", 0, 1), /ions must be/);
});
