// The worker controller (web/worker.js) over the real built Simulation, not a
// stand-in: a small run to the end, then a stop part-way. Guards the seam the
// stand-in tests in worker.test.mjs cannot see, such as a method the
// controller calls that the WASM class does not export.
//
//   LINDHARD_WASM_PKG=<dir> node --test lindhard-wasm/tests/node/worker_real.test.mjs

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { resolve } from "node:path";
import test from "node:test";
import { createController } from "../../web/worker.js";

const pkg = process.env.LINDHARD_WASM_PKG;
assert.ok(pkg, "set LINDHARD_WASM_PKG to the directory build.sh wrote");
const wasm = createRequire(import.meta.url)(resolve(pkg, "lindhard_wasm.js"));

const IONS = 200;
const toml = readFileSync(
  new URL("../../../examples/as_50keV_si_sio2.toml", import.meta.url),
  "utf8",
).replace("ions = 2000", `ions = ${IONS}`);
assert.ok(toml.includes(`ions = ${IONS}`), "example changed");

const rig = () => {
  const out = [];
  return { c: createController(wasm, (m) => out.push(m)), out };
};

test("a real session runs to the end through the controller", async () => {
  const { c, out } = rig();
  c.handle({ type: "init", toml, ions: IONS });
  assert.equal(out[0].type, "ready", JSON.stringify(out[0]));
  c.handle({ type: "run" });
  await c.idle;
  const last = out[out.length - 1];
  assert.equal(last.type, "finished", JSON.stringify(last));
  assert.equal(last.done, IONS);
  assert.equal(last.total, IONS);
  assert.ok(out.some((m) => m.type === "progress"));
});

test("stop on a real session reports stopped, not an error", async () => {
  const { c, out } = rig();
  c.handle({ type: "init", toml, ions: IONS });
  c.handle({ type: "run" });
  c.handle({ type: "pause" });
  await c.idle;
  c.handle({ type: "stop" });
  const last = out[out.length - 1];
  assert.equal(last.type, "stopped", JSON.stringify(last));
  assert.equal(last.total, IONS);
});
