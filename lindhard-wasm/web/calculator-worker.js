// Runs the range calculator off the main thread so the page stays responsive
// during a long run. One message in, one message out.
//
//   in:  { ion, energyEv, target, ions, seed }
//   out: { ok: true, result } | { ok: false, message, field? }
//        field is "ion" or "target" when the element symbol is unknown.
//
// `calculate` is the same single-threaded call as before, so a given seed
// gives the same numbers.

import init, { calculate, elements_json } from "./pkg/lindhard_wasm.js";

const ready = init().then(() => new Set(JSON.parse(elements_json()).map((e) => e.symbol.toLowerCase())));

self.onmessage = async (event) => {
  const { ion, energyEv, target, ions, seed } = event.data;
  try {
    const symbols = await ready;
    for (const [field, sym] of [["ion", ion], ["target", target]]) {
      if (!symbols.has(sym.toLowerCase())) {
        self.postMessage({ ok: false, field, message: `Unknown element symbol “${sym}”.` });
        return;
      }
    }
    const r = calculate(ion, energyEv, target, ions, seed);
    const result = {
      rp_nm: r.rp_nm, drp_nm: r.drp_nm,
      rp_std_err_nm: r.rp_std_err_nm, drp_std_err_nm: r.drp_std_err_nm,
      stopped: r.stopped, escaped: r.escaped, histories: r.histories,
      bin_nm: r.bin_nm, counts: Array.from(r.counts), warnings: r.warnings,
    };
    r.free?.();
    self.postMessage({ ok: true, result });
  } catch (err) {
    self.postMessage({ ok: false, message: String(err?.message ?? err) });
  }
};
