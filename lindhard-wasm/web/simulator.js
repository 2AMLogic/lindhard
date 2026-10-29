// The simulator page: Setup (form <-> lindhard TOML) and a minimal Run view.
// The Worker (worker.js) owns the WASM session; this file only draws.

import init, { elements_json } from "./pkg/lindhard_wasm.js";
import { POTENTIALS, STOPPINGS, MAX_IONS, ModelError, defaultSetup, fromToml, toToml, validate } from "./toml_model.js";

const $ = (id) => document.getElementById(id);
const UNIT_EV = { eV: 1, keV: 1e3, MeV: 1e6 };
const UNIT_NM = { "Å": 0.1, nm: 1, "µm": 1000 };

// Citations of the choices (docs: book/src/models).
const CITE = {
  zbl: "ZBL universal screening: Ziegler, Biersack and Littmark (1985), ch. 2",
  "kr-c": "Kr-C: Wilson, Haggmark and Biersack (1977)",
  moliere: "Molière: Molière (1947)",
  "lenz-jensen": "Lenz-Jensen: Lenz (1932), Jensen (1932)",
  "lindhard-scharff": "Lindhard-Scharff: Lindhard and Scharff (1961); Lindhard, Scharff and Schiøtt (1963)",
  "bethe-bloch": "Bethe-Bloch: Bethe (1930, 1932), Fano (1963); a high-velocity model, the engine warns outside its range",
  "equipartition-ls-or": "Lindhard-Scharff and Oen-Robinson equipartition mix: Oen and Robinson (1976); see the stopping-model notes",
};

// Presets. SiO2 density: the usual handbook density of vitreous silica (CRC
// Handbook; docs/data-provenance.md, example inputs). GaAs: 5.321 g/cm^3,
// NBS Monograph 25 sec. 3 (docs/data-provenance.md, lattice constants). Si,
// Cu, Au: the element table's densities. The O energies are the illustrative
// model parameters of examples/as_50keV_si_sio2.toml, not recommendations.
const el = (symbol, fraction = 1, eDEv = null, eSEv = null) => ({ symbol, fraction, eDEv, eBEv: null, eSEv });
const PRESETS = {
  Si: { density: null, elements: [el("Si")] },
  Cu: { density: null, elements: [el("Cu")] },
  Au: { density: null, elements: [el("Au")] },
  SiO2: { density: 2.2, elements: [el("Si"), el("O", 2, 20, 2)], note: "O energies are illustrative model parameters" },
  GaAs: { density: 5.321, elements: [el("Ga"), el("As")] },
};

let elements = [];
let table = new Map();
let setup = defaultSetup();
let unitE = "keV";
const unitT = []; // per-layer thickness unit
let worker = null;

function h(tag, props = {}, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (k === "on") for (const [ev, fn] of Object.entries(v)) e.addEventListener(ev, fn);
    else if (k in e && k !== "list") e[k] = v;
    else e.setAttribute(k, v);
  }
  e.append(...kids);
  return e;
}
const numOrNull = (s) => (s.trim() === "" ? null : Number(s));
const show = (x) => (x == null ? "" : String(x));

function status(text, error = false) {
  $("status").textContent = text;
  $("status").className = error ? "error" : "";
}

// ---- form <-> setup --------------------------------------------------------

function renderLayers() {
  const root = $("layers");
  root.replaceChildren(...setup.layers.map(layerView));
  if (root.firstChild === null) root.append(h("p", { className: "flag" }, "No layers: add one."));
}

function layerView(l, i) {
  const last = i === setup.layers.length - 1;
  const u = (unitT[i] ??= "nm");
  const upd = (fn) => (e) => { fn(e.target); refresh(); };
  const thick = h("input", {
    type: "number", min: 0, step: "any", id: `thick-${i}`, disabled: l.thicknessNm == null,
    value: l.thicknessNm == null ? "" : String(l.thicknessNm / UNIT_NM[u]),
    on: { input: upd((t) => (l.thicknessNm = Number(t.value) * UNIT_NM[unitT[i]])) },
  });
  const unit = h("select", {
    id: `tunit-${i}`,
    on: { change: (e) => { unitT[i] = e.target.value; if (l.thicknessNm != null) thick.value = String(l.thicknessNm / UNIT_NM[unitT[i]]); } },
  }, ...Object.keys(UNIT_NM).map((k) => h("option", { selected: k === u }, k)));
  const semi = h("input", {
    type: "checkbox", id: `semi-${i}`, checked: l.thicknessNm == null, disabled: !last,
    on: { change: upd((c) => { l.thicknessNm = c.checked ? null : 100; renderLayers(); }) },
  });
  const pure = l.elements.length === 1;
  const tab = table.get(l.elements[0]?.symbol);
  const density = h("input", {
    type: "number", min: 0, step: "any", id: `dens-${i}`, value: show(l.densityGcm3),
    placeholder: pure && tab?.density_g_cm3 != null ? String(tab.density_g_cm3) : "required",
    on: { input: upd((t) => (l.densityGcm3 = numOrNull(t.value))) },
  });
  const densSource = pure && l.densityGcm3 == null && tab?.density_g_cm3 != null
    ? "tabulated density of the element (X-Ray Data Booklet, via the element table)"
    : l.densityGcm3 != null ? "entered" : "none: enter one";

  const rows = l.elements.map((e, j) => {
    const row = table.get(e.symbol);
    const field = (key, tabKey, label) => h("td", {}, h("input", {
      type: "number", min: 0, step: "any", "aria-label": `${label} of ${e.symbol} in ${l.name}`,
      value: show(e[key]), placeholder: tabKey === null ? "0" : row?.[tabKey] != null ? String(row[tabKey]) : "required",
      on: { input: upd((t) => (e[key] = numOrNull(t.value))) },
    }));
    return h("tr", {},
      h("td", {}, h("input", { value: e.symbol, size: 3, list: "symbols", "aria-label": `element ${j + 1} of ${l.name}`, on: { input: upd((t) => { e.symbol = t.value.trim(); renderLayers(); }) } })),
      h("td", {}, h("input", { type: "number", min: 0, step: "any", value: String(e.fraction), "aria-label": `amount of ${e.symbol} in ${l.name}`, on: { input: upd((t) => (e.fraction = Number(t.value))) } })),
      field("eDEv", "default_displacement_ev", "displacement energy"),
      field("eBEv", null, "lattice binding energy"),
      field("eSEv", "default_surface_binding_ev", "surface binding energy"),
      h("td", {}, h("button", { type: "button", "aria-label": `remove ${e.symbol} from ${l.name}`, on: { click: () => { l.elements.splice(j, 1); renderLayers(); refresh(); } } }, "×")),
    );
  });

  return h("div", { className: "layer" },
    h("h3", {}, `Layer ${i + 1}`),
    h("div", { className: "row" },
      h("label", {}, "Name", h("input", { value: l.name, "aria-label": `name of layer ${i + 1}`, on: { input: upd((t) => (l.name = t.value.trim())) } })),
      h("label", {}, "Thickness", thick),
      h("label", {}, "Unit", unit),
      h("label", { className: "inline" }, semi, "Semi-infinite (substrate)"),
      h("label", {}, "Density (g/cm³)", density),
    ),
    h("p", { className: "note" }, `Density source: ${densSource}.`),
    h("table", { className: "elements" },
      h("tr", {}, ...["Element", "Atom fraction", "E_d (eV)", "E_b (eV)", "E_s (eV)", ""].map((t) => h("th", { scope: "col" }, t))),
      ...rows),
    h("div", { className: "tools" },
      h("button", { type: "button", on: { click: () => { l.elements.push(el("Si")); renderLayers(); refresh(); } } }, "Add element"),
      h("button", { type: "button", disabled: i === 0, on: { click: () => move(i, -1) } }, "Move up"),
      h("button", { type: "button", disabled: last, on: { click: () => move(i, 1) } }, "Move down"),
      h("button", { type: "button", on: { click: () => { setup.layers.splice(i, 1); unitT.splice(i, 1); renderLayers(); refresh(); } } }, "Delete layer")),
  );
}

function move(i, d) {
  const L = setup.layers;
  [L[i], L[i + d]] = [L[i + d], L[i]];
  [unitT[i], unitT[i + d]] = [unitT[i + d], unitT[i]];
  // Only the last layer can be a substrate.
  for (const [k, l] of L.entries()) if (l.thicknessNm == null && k !== L.length - 1) l.thicknessNm = 100;
  renderLayers();
  refresh();
}

function uniqueName(base) {
  let n = base, k = 2;
  while (setup.layers.some((l) => l.name === n)) n = `${base}-${k++}`;
  return n;
}

function addPreset(key) {
  const p = PRESETS[key];
  const last = setup.layers[setup.layers.length - 1];
  if (last && last.thicknessNm == null) last.thicknessNm = 100; // a new layer goes below
  setup.layers.push({
    name: uniqueName(key), thicknessNm: null, densityGcm3: p.density,
    elements: p.elements.map((e) => ({ ...e })),
  });
  renderLayers();
  refresh();
}

function readScalars() {
  const b = setup.beam;
  b.ion = $("ion").value.trim();
  b.massAmu = numOrNull($("mass").value);
  unitE = $("energy-unit").value;
  b.energyEv = Number($("energy").value) * UNIT_EV[unitE];
  b.tiltDeg = Number($("tilt").value);
  const p = setup.physics;
  p.followRecoils = $("mode").value === "full";
  p.potential = $("potential").value;
  p.stopping = $("stopping").value;
  p.primaryCutoffEv = Number($("pcut").value);
  p.recoilCutoffEv = Number($("rcut").value);
  setup.run.ions = Number($("ions").value);
  setup.run.seed = Number($("seed").value);
  setup.tally.depthBinNm = Number($("bin").value);
  setup.tally.depthBins = Number($("bins").value);
}

function writeScalars() {
  const b = setup.beam, p = setup.physics;
  $("ion").value = b.ion;
  $("mass").value = show(b.massAmu);
  $("energy-unit").value = unitE;
  $("energy").value = String(b.energyEv / UNIT_EV[unitE]);
  $("tilt").value = String(b.tiltDeg);
  $("mode").value = p.followRecoils ? "full" : "quick";
  $("potential").value = p.potential;
  $("stopping").value = p.stopping;
  $("pcut").value = String(p.primaryCutoffEv);
  $("rcut").value = String(p.recoilCutoffEv);
  $("ions").value = String(setup.run.ions);
  $("seed").value = String(setup.run.seed);
  $("bin").value = String(setup.tally.depthBinNm);
  $("bins").value = String(setup.tally.depthBins);
}

let problems = [];
function refresh() {
  readScalars();
  const row = table.get(setup.beam.ion);
  $("mass").placeholder = row ? `${row.atomic_weight} (standard atomic weight)` : "";
  $("mass-note").textContent = row?.weight_is_mass_number
    ? `${row.symbol} has no standard atomic weight; the mass number of its longest-lived isotope is used.` : "";
  $("model-cite").textContent = `${CITE[setup.physics.potential]}. ${CITE[setup.physics.stopping]}.`;
  problems = validate(setup, elements);
  $("problems").replaceChildren(...problems.map((p) => h("li", {}, `${p.path}: ${p.message}`)));
  $("problems-box").hidden = problems.length === 0;
  $("toml").value = toToml(setup);
  $("btn-run").disabled = problems.length > 0 || running();
}

// ---- input file ------------------------------------------------------------

$("download").addEventListener("click", () => {
  const a = h("a", { href: URL.createObjectURL(new Blob([toToml(setup)], { type: "text/plain" })), download: "lindhard-input.toml" });
  a.click();
  URL.revokeObjectURL(a.href);
});

$("load").addEventListener("change", async (e) => {
  const file = e.target.files[0];
  if (!file) return;
  try {
    setup = fromToml(await file.text());
    unitT.length = 0;
    writeScalars();
    renderLayers();
    refresh();
    status(`Loaded ${file.name}.`);
  } catch (err) {
    if (!(err instanceof ModelError)) throw err;
    status(`${file.name} was not loaded: ${err.message}`, true);
  }
  e.target.value = "";
});

// ---- run -------------------------------------------------------------------

const running = () => worker !== null;
const fmt = (x) => (x >= 100 ? x.toFixed(1) : x >= 10 ? x.toFixed(2) : x.toPrecision(3));
const nm = (m) => m * 1e9;

function setButtons(state) {
  $("btn-run").disabled = state !== "idle" || problems.length > 0;
  $("btn-pause").disabled = state !== "running";
  $("btn-resume").disabled = state !== "paused";
  $("btn-stop").disabled = state === "idle";
}

function counters(done, total, report) {
  $("progress").value = total ? done / total : 0;
  $("progress-text").textContent = `${done} of ${total} ions`;
  $("c-done").textContent = `${done} / ${total}`;
  const d = report?.range?.depth;
  $("c-rp").textContent = d ? `${fmt(nm(d.mean))} nm (± ${fmt(nm(d.mean_std_err))})` : "–";
  $("c-drp").textContent = d ? `${fmt(nm(d.std_dev))} nm (± ${fmt(nm(d.std_dev_std_err))})` : "–";
  $("c-stopped").textContent = report ? String(report.range.stopped) : "–";
}

function startRun() {
  refresh();
  if (problems.length) return;
  worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
  setButtons("running");
  status("Starting…");
  worker.onmessage = ({ data: m }) => {
    switch (m.type) {
      case "ready":
        $("warnings").textContent = m.warnings.join(" ");
        worker.postMessage({ type: "run" });
        status("Running…");
        break;
      case "progress":
        counters(m.done, m.total, m.report);
        break;
      case "finished":
        counters(m.done, m.total, m.report);
        status("Done (± is one standard error).");
        endRun();
        break;
      case "stopped":
        status(`Stopped after ${m.done} of ${m.total} ions.`);
        endRun();
        break;
      case "error":
        status(m.message, true);
        endRun();
        break;
    }
  };
  worker.onerror = (e) => { status(`The worker failed: ${e.message}`, true); endRun(); };
  worker.postMessage({ type: "init", toml: toToml(setup), ions: setup.run.ions, recordHistories: 0, maxPoints: 1000 });
  counters(0, setup.run.ions, null);
}

function endRun() {
  worker?.terminate();
  worker = null;
  setButtons("idle");
}

$("btn-run").addEventListener("click", startRun);
$("btn-pause").addEventListener("click", () => { worker.postMessage({ type: "pause" }); setButtons("paused"); status("Paused."); });
$("btn-resume").addEventListener("click", () => { worker.postMessage({ type: "resume" }); setButtons("running"); status("Running…"); });
$("btn-stop").addEventListener("click", () => worker.postMessage({ type: "stop" }));

// ---- tabs ------------------------------------------------------------------

function tab(which) {
  for (const t of ["setup", "run"]) {
    $(t).hidden = t !== which;
    $(`tab-${t}`).setAttribute("aria-selected", String(t === which));
  }
}
$("tab-setup").addEventListener("click", () => tab("setup"));
$("tab-run").addEventListener("click", () => tab("run"));
for (const [id, other] of [["tab-setup", "tab-run"], ["tab-run", "tab-setup"]]) {
  $(id).addEventListener("keydown", (e) => {
    if (e.key === "ArrowRight" || e.key === "ArrowLeft") { $(other).focus(); $(other).click(); }
  });
}

// ---- start -----------------------------------------------------------------

try {
  await init();
  elements = JSON.parse(elements_json());
  table = new Map(elements.map((e) => [e.symbol, e]));
  $("symbols").replaceChildren(...elements.map((e) => h("option", { value: e.symbol })));
  $("potential").replaceChildren(...POTENTIALS.map((p) => h("option", {}, p)));
  $("stopping").replaceChildren(...STOPPINGS.map((p) => h("option", {}, p)));
  $("preset").replaceChildren(...Object.keys(PRESETS).map((p) => h("option", {}, p)));
  $("ions").max = String(MAX_IONS);
  $("add-layer").addEventListener("click", () => addPreset("Si"));
  $("add-preset").addEventListener("click", () => addPreset($("preset").value));
  for (const id of ["ion", "mass", "energy", "energy-unit", "tilt", "mode", "ions", "seed", "bin", "bins", "potential", "stopping", "pcut", "rcut"]) {
    $(id).addEventListener("input", refresh);
    $(id).addEventListener("change", refresh);
  }
  writeScalars();
  renderLayers();
  refresh();
  setButtons("idle");
  status("Ready: B at 10 keV into Si.");
} catch (err) {
  status(`The engine could not be loaded: ${err.message ?? err}`, true);
}
