// The WASM half of the native-versus-WASM session parity test, and the Node
// checks of the Simulation API: run the same input as the native test
// (tests/session.rs, native_session_matches_the_golden_record) through the
// WASM build and compare report_json() with the native golden record, with
// the same tolerance and for the same reason (stated in tests/session.rs and
// tests/parity.rs).
//
//   LINDHARD_WASM_PKG=<dir> node --test lindhard-wasm/tests/node/session.test.mjs

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { resolve } from "node:path";
import test from "node:test";

const pkg = process.env.LINDHARD_WASM_PKG;
assert.ok(pkg, "set LINDHARD_WASM_PKG to the directory build.sh wrote");
const wasm = createRequire(import.meta.url)(resolve(pkg, "lindhard_wasm.js"));
const read = (rel) => readFileSync(new URL(rel, import.meta.url), "utf8");
const golden = JSON.parse(read("../as_50keV_si_sio2.golden.json"));

// The input of tests/session.rs: the example with the golden ion count and
// the recoil cascade followed.
function input(ions) {
  const t = read("../../../examples/as_50keV_si_sio2.toml").replace(
    "ions = 2000",
    `ions = ${ions}`,
  );
  assert.ok(t.includes(`ions = ${ions}`), "example changed");
  return t.includes("follow_recoils")
    ? t
    : t.replace("[physics]\n", "[physics]\nfollow_recoils = true\n");
}

function runAll(toml, record, maxPoints) {
  const sim = new wasm.Simulation(toml, record, maxPoints);
  try {
    while (!sim.done()) sim.step(1000);
    return { report: JSON.parse(sim.report_json()), sim };
  } finally {
    sim.free();
  }
}

test("the compound two-layer report matches the native golden record", () => {
  const { report } = runAll(input(golden.ions), 0, 10);
  const d = report.range.depth;
  assert.ok(
    Math.abs(d.mean - golden.depth_mean_m) <= golden.depth_mean_std_err_m,
    `mean ${d.mean} vs ${golden.depth_mean_m}`,
  );
  assert.ok(
    Math.abs(d.std_dev - golden.depth_std_dev_m) <=
      golden.depth_std_dev_std_err_m,
    `std_dev ${d.std_dev} vs ${golden.depth_std_dev_m}`,
  );
  const sigma = Math.max(
    1,
    Math.sqrt(golden.stopped * (1 - golden.stopped / golden.ions)),
  );
  assert.ok(Math.abs(report.range.stopped - golden.stopped) <= sigma);
  assert.equal(report.histories, golden.ions);
});

test("chunking does not change the report or the trajectories", () => {
  const run = (steps) => {
    const sim = new wasm.Simulation(input(300), 10, 200);
    const traj = [];
    for (let i = 0; !sim.done(); i++) {
      sim.step(steps[Math.min(i, steps.length - 1)]);
      traj.push(...JSON.parse(sim.take_trajectories_json()).histories);
    }
    const out = [sim.report_json(), JSON.stringify(traj)];
    sim.free();
    return out;
  };
  const [reportA, trajA] = run([300]);
  const [reportB, trajB] = run([1, 7, 100, 300]);
  assert.equal(reportA, reportB);
  assert.equal(trajA, trajB);
  assert.equal(JSON.parse(trajA).length, 10);
});

test("input errors throw an Error, they do not abort the module", () => {
  assert.throws(() => new wasm.Simulation("this is [not toml", 0, 10), Error);
  assert.throws(
    () => new wasm.Simulation(input(10).replace('ion = "As"', 'ion = "Xx"'), 0, 10),
    Error,
  );
  assert.throws(
    () =>
      new wasm.Simulation(
        input(10).replace("energy_ev = 50000.0", "energy_ev = -5.0"),
        0,
        10,
      ),
    Error,
  );
  // The module is still usable after the errors.
  assert.ok(runAll(input(64), 0, 10).report.histories === 64);
});

test("elements_json comes from the library table", () => {
  const rows = JSON.parse(wasm.elements_json());
  assert.equal(rows.length, 92);
  const si = rows.find((r) => r.symbol === "Si");
  assert.equal(si.z, 14);
  assert.ok(si.atomic_weight > 28 && si.atomic_weight < 28.2);
});
