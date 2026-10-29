// The Worker that owns the WASM `Simulation`: the main thread only draws.
//
// Messages to the worker
//   { type: "init", toml, ions, recordHistories, maxPoints }  build a session (ions: run.ions)
//   { type: "run" } | { type: "pause" } | { type: "resume" } | { type: "stop" }
// Messages from the worker
//   { type: "ready", warnings }
//   { type: "progress", done, total, report }   after every step
//   { type: "finished", done, total, report } | { type: "stopped", done, total }
//   { type: "error", message }
//
// `createController` holds the protocol and takes the WASM module and a
// `post` function, so the Node tests drive it without a browser.

const TARGET_STEP_MS = 50;

export function createController(wasm, post) {
  let sim = null;
  let state = "idle"; // idle | running | paused
  let loop = null;
  let chunk = 64; // histories per step, adapted to the measured rate

  const free = () => {
    sim?.free?.();
    sim = null;
  };
  const fail = (e) => {
    state = "idle";
    free();
    post({ type: "error", message: String(e?.message ?? e) });
  };

  // A resume before the previous loop has noticed its pause starts a new
  // loop; the generation retires the old one.
  let generation = 0;
  async function drive() {
    const mine = ++generation;
    while (state === "running" && mine === generation) {
      const t0 = performance.now();
      const before = sim.histories_done();
      let done;
      try {
        done = sim.step(chunk);
        post({ type: "progress", done, total: total(), report: JSON.parse(sim.report_json()), trajectories: JSON.parse(sim.take_trajectories_json()).histories });
      } catch (e) {
        return fail(e);
      }
      const ms = Math.max(performance.now() - t0, 0.5);
      const ran = Number(done) - Number(before);
      chunk = Math.max(1, Math.min(100000, Math.round((ran * TARGET_STEP_MS) / ms)));
      if (sim.done()) {
        state = "idle";
        post({ type: "finished", done, total: total(), report: JSON.parse(sim.report_json()) });
        return;
      }
      // Yield to the event loop so a pause or stop message is handled
      // between steps: Pause takes effect within one step.
      await new Promise((r) => setTimeout(r, 0));
    }
  }
  let totalIons = 0;
  const total = () => totalIons;

  return {
    handle(msg) {
      try {
        switch (msg.type) {
          case "init":
            state = "idle";
            free();
            sim = new wasm.Simulation(msg.toml, msg.recordHistories ?? 0, msg.maxPoints ?? 1000);
            totalIons = msg.ions;
            chunk = 64;
            post({ type: "ready", warnings: Array.from(sim.warnings()) });
            break;
          case "run":
          case "resume":
            if (sim && state !== "running") {
              state = "running";
              loop = drive();
            }
            break;
          case "pause":
            if (state === "running") state = "paused";
            break;
          case "stop": {
            const done = sim ? Number(sim.histories_done()) : 0;
            state = "idle";
            free();
            post({ type: "stopped", done, total: totalIons });
            break;
          }
        }
      } catch (e) {
        fail(e);
      }
    },
    get state() {
      return state;
    },
    get idle() {
      return loop;
    },
  };
}

if (typeof WorkerGlobalScope !== "undefined" && self instanceof WorkerGlobalScope) {
  const wasm = await import("./pkg/lindhard_wasm.js");
  await wasm.default();
  const c = createController(wasm, (m) => self.postMessage(m));
  self.onmessage = (e) => c.handle(e.data);
}
