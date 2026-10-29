// The worker message protocol (web/worker.js) against a stand-in session:
// progress after every step, Pause within one step, Resume, Stop, errors.
//
//   node --test lindhard-wasm/tests/node/worker.test.mjs

import assert from "node:assert/strict";
import test from "node:test";
import { createController } from "../../web/worker.js";

class FakeSimulation {
  static throwOnNew = false;
  constructor(toml) {
    if (FakeSimulation.throwOnNew) throw new Error("bad input");
    this.toml = toml;
    this.n = 0;
    this.total = 1000;
  }
  step(count) {
    this.n = Math.min(this.total, this.n + Math.max(count, 64));
    return this.n;
  }
  done() { return this.n >= this.total; }
  histories_done() { return BigInt(this.n); }
  warnings() { return ["w"]; }
  report_json() { return JSON.stringify({ histories: this.n }); }
  take_trajectories_json() { return '{"histories":[]}'; }
  free() { this.freed = true; }
}

const rig = () => {
  const out = [];
  const c = createController({ Simulation: FakeSimulation }, (m) => out.push(m));
  return { c, out };
};
const tick = () => new Promise((r) => setTimeout(r, 5));

test("init, run to the end: ready, progress per step, finished", async () => {
  const { c, out } = rig();
  c.handle({ type: "init", toml: "x", ions: 1000 });
  assert.deepEqual(out[0], { type: "ready", warnings: ["w"] });
  c.handle({ type: "run" });
  await c.idle;
  const kinds = out.map((m) => m.type);
  assert.equal(kinds[kinds.length - 1], "finished");
  assert.ok(kinds.filter((k) => k === "progress").length >= 1);
  const last = out[out.length - 1];
  assert.equal(last.done, 1000);
  assert.equal(last.total, 1000);
  assert.equal(last.report.histories, 1000);
  assert.equal(c.state, "idle");
});

test("pause stops further steps; resume continues to the end", async () => {
  const { c, out } = rig();
  c.handle({ type: "init", toml: "x", ions: 1000 });
  c.handle({ type: "run" });
  c.handle({ type: "pause" }); // handled before the loop's first yield ends
  await c.idle;
  const paused = out.filter((m) => m.type === "progress").length;
  assert.equal(c.state, "paused");
  await tick();
  assert.equal(out.filter((m) => m.type === "progress").length, paused);
  assert.ok(!out.some((m) => m.type === "finished"));
  c.handle({ type: "resume" });
  await c.idle;
  assert.equal(out[out.length - 1].type, "finished");
});

test("stop reports the histories done and releases the session", async () => {
  const { c, out } = rig();
  c.handle({ type: "init", toml: "x", ions: 1000 });
  c.handle({ type: "run" });
  c.handle({ type: "pause" });
  await c.idle;
  c.handle({ type: "stop" });
  const stopped = out[out.length - 1];
  assert.equal(stopped.type, "stopped");
  assert.equal(stopped.total, 1000);
  assert.ok(stopped.done > 0 && stopped.done < 1000);
  c.handle({ type: "run" }); // no session: ignored
  assert.equal(c.state, "idle");
});

test("an input error comes back as an error message", () => {
  const { c, out } = rig();
  FakeSimulation.throwOnNew = true;
  c.handle({ type: "init", toml: "x", ions: 1 });
  FakeSimulation.throwOnNew = false;
  assert.deepEqual(out, [{ type: "error", message: "bad input" }]);
});
