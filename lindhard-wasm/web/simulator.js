// The simulator page: Setup (form <-> lindhard TOML) and a minimal Run view.
// The Worker (worker.js) owns the WASM session; this file only draws.

import init, { elements_json } from "./pkg/lindhard_wasm.js";
import { POTENTIALS, STOPPINGS, MAX_IONS, ModelError, defaultSetup, fromToml, toToml, validate } from "./toml_model.js";
import { buildPlots, caveatLines, dashForSeries, niceTicks, plotExtent, reportJson, toCsv } from "./plots.js";
import {
  INTERACTIVE_MAX_IONS, RECORDED_HISTORIES, MAX_POINTS_PER_TRACK, anyTruncated, dashFor, deriveCounters,
  layerSpans, segmentBounds, speciesLabels, trajectorySegments,
} from "./run_view.js";

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
  if (setup.run.ions > INTERACTIVE_MAX_IONS) {
    problems.push({ path: "run.ions", message: `the page runs at most ${INTERACTIVE_MAX_IONS} ions; use lindhard-cli for more` });
  }
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

const pm = (x, e) => `${fmt(nm(x))} nm (± ${fmt(nm(e))})`;
const dash = (m) => (m ? pm(m.mean, m.err) : "–");
const pct = (f) => (f == null ? "–" : `${(100 * f).toFixed(1)} %`);
const perIon = (x) => (x == null ? "–" : x.toPrecision(3));

function counters(done, total, report) {
  const c = deriveCounters(report, { done, total, followRecoils: setup.physics.followRecoils });
  $("progress").value = c.fraction;
  $("progress-text").textContent = `${done} of ${total} ions`;
  $("c-done").textContent = `${done} / ${total}`;
  $("c-rp").textContent = dash(c.depth);
  $("c-drp").textContent = c.depth ? pm(c.depth.sd, c.depth.sdErr) : "–";
  $("c-ly").textContent = c.lateralY ? `${fmt(nm(c.lateralY.sd))} nm (± ${fmt(nm(c.lateralY.sdErr))})` : "–";
  $("c-lz").textContent = c.lateralZ ? `${fmt(nm(c.lateralZ.sd))} nm (± ${fmt(nm(c.lateralZ.sdErr))})` : "–";
  $("c-rad").textContent = dash(c.radial);
  $("c-stopped").textContent = c.stopped == null ? "–" : String(c.stopped);
  $("c-back").textContent = perIon(c.backscatter);
  $("c-trans").textContent = perIon(c.transmission);
  $("c-sput").textContent = perIon(c.sputterYield);
  $("c-vac").textContent = perIon(c.vacancies);
  $("c-vac-mode").textContent = report ? `(${c.modeLabel})` : "";
  $("c-eion").textContent = pct(c.partition?.ionization);
  $("c-enuc").textContent = pct(c.partition?.nuclear);
  $("c-eesc").textContent = pct(c.partition?.escaped);
  $("c-eoth").textContent = pct(c.partition?.other);
}

// ---- trajectory view -------------------------------------------------------

// Okabe-Ito colours, safe for the common colour-vision deficiencies; the
// second cue is the dash pattern (generation) and an end marker (species).
const SPECIES_COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"];
const MARKERS = ["circle", "square", "triangle", "diamond", "cross", "plus"];
let histories = []; // recorded so far in this run
let symbolOfZ = () => undefined;

function marker(ctx, kind, x, y, r) {
  ctx.beginPath();
  if (kind === "circle") ctx.arc(x, y, r, 0, 2 * Math.PI);
  else if (kind === "square") ctx.rect(x - r, y - r, 2 * r, 2 * r);
  else if (kind === "triangle") { ctx.moveTo(x, y - r); ctx.lineTo(x + r, y + r); ctx.lineTo(x - r, y + r); ctx.closePath(); }
  else if (kind === "diamond") { ctx.moveTo(x, y - r); ctx.lineTo(x + r, y); ctx.lineTo(x, y + r); ctx.lineTo(x - r, y); ctx.closePath(); }
  else if (kind === "cross") { ctx.moveTo(x - r, y - r); ctx.lineTo(x + r, y + r); ctx.moveTo(x + r, y - r); ctx.lineTo(x - r, y + r); }
  else { ctx.moveTo(x - r, y); ctx.lineTo(x + r, y); ctx.moveTo(x, y - r); ctx.lineTo(x, y + r); }
  ctx.fill();
  ctx.stroke();
}

function drawTrajectories() {
  const cv = $("traj"), ctx = cv.getContext("2d");
  const css = getComputedStyle(document.documentElement);
  const fg = css.getPropertyValue("--fg").trim() || "#000", bg = css.getPropertyValue("--bg").trim() || "#fff";
  const panel = css.getPropertyValue("--panel").trim() || "#eee", muted = css.getPropertyValue("--muted").trim() || "#666";
  const view = $("traj-view").value;
  const segs = trajectorySegments(histories, { view, cascades: $("traj-cascades").checked });
  const bd = segmentBounds(segs, view);
  const W = cv.width, H = cv.height, pad = 28;
  const sx = (u) => pad + ((u - bd.u0) / (bd.u1 - bd.u0)) * (W - 2 * pad);
  const sy = (v) => pad + ((v - bd.v0) / (bd.v1 - bd.v0)) * (H - 2 * pad);
  ctx.setLineDash([]);
  ctx.fillStyle = bg;
  ctx.fillRect(0, 0, W, H);
  ctx.font = "12px system-ui, sans-serif";
  ctx.lineWidth = 1;
  if (view === "side") {
    // Layers shaded alternately, boundaries as lines, names at the top.
    layerSpans(setup.layers).forEach((l, i) => {
      const x0 = sx(Math.max(l.top, bd.u0)), x1 = sx(Math.min(l.bottom, bd.u1));
      if (x1 <= x0) return;
      ctx.fillStyle = i % 2 ? bg : panel;
      ctx.fillRect(x0, 0, x1 - x0, H);
      ctx.strokeStyle = muted;
      ctx.beginPath(); ctx.moveTo(x0, 0); ctx.lineTo(x0, H); ctx.stroke();
      ctx.fillStyle = muted;
      ctx.fillText(l.name, x0 + 4, 14);
    });
  }
  const axes = view === "side" ? ["depth", "y"] : ["y", "z"];
  ctx.fillStyle = muted;
  ctx.fillText(`${axes[0]} → ${fmt(nm(bd.u1))} nm`, W - 160, H - 8);
  ctx.fillText(`${axes[1]} ±${fmt(nm(bd.v1))} nm`, 4, H - 8);
  ctx.lineWidth = 1.2;
  for (const s of segs) {
    ctx.strokeStyle = SPECIES_COLOURS[s.species % SPECIES_COLOURS.length];
    ctx.setLineDash(dashFor(s.generation));
    ctx.beginPath(); ctx.moveTo(sx(s.a[0]), sy(s.a[1])); ctx.lineTo(sx(s.b[0]), sy(s.b[1])); ctx.stroke();
  }
  ctx.setLineDash([]);
  // End markers: one per track, shape by species.
  for (const h of histories) for (const t of h.tracks) {
    if (t.generation > 0 && !$("traj-cascades").checked) continue;
    const q = t.points[t.points.length - 1];
    if (!q) continue;
    const [u, v] = view === "side" ? [q[0], q[1]] : [q[1], q[2]];
    ctx.fillStyle = SPECIES_COLOURS[t.species % SPECIES_COLOURS.length];
    ctx.strokeStyle = fg;
    marker(ctx, MARKERS[t.species % MARKERS.length], sx(u), sy(v), t.generation ? 2.5 : 4);
  }
  const n = histories.length;
  $("traj-note").textContent = n === 0
    ? `No trajectories yet. The first ${RECORDED_HISTORIES} ion histories of a run are recorded.`
    : `${n} recorded histor${n === 1 ? "y" : "ies"} (the first ${RECORDED_HISTORIES} of the run), ${segs.length} segments.` +
      (anyTruncated(histories) ? " Some histories were cut at the point limit of each track, so their tails are missing." : "");
}

function drawLegend(report) {
  const names = speciesLabels(report, symbolOfZ);
  const item = (text, swatch) => h("li", {}, swatch, ` ${text}`);
  const line = (i, g) => {
    const c = document.createElement("canvas");
    c.width = 34; c.height = 12;
    c.style.cssText = "width:34px;height:12px;border:0";
    const x = c.getContext("2d");
    x.strokeStyle = SPECIES_COLOURS[i % SPECIES_COLOURS.length]; x.fillStyle = x.strokeStyle; x.lineWidth = 2;
    x.setLineDash(dashFor(g));
    x.beginPath(); x.moveTo(0, 6); x.lineTo(26, 6); x.stroke();
    x.setLineDash([]);
    x.strokeStyle = "#888";
    marker(x, MARKERS[i % MARKERS.length], 29, 6, 3);
    return c;
  };
  $("traj-legend").replaceChildren(
    ...names.map((n, i) => item(`${n}: ${MARKERS[i % MARKERS.length]} at the track's end`, line(i, 0))),
    item("solid: primary ion", line(0, 0)),
    item("dashed: first recoils", line(0, 1)),
    item("dotted: later generations", line(0, 2)),
  );
}

// ---- distribution plots and exports ------------------------------------------

let lastReport = null;
let plotSpecs = [];
const plotViews = new Map(); // id -> {canvas, legend, note}

function download(name, blob) {
  const a = h("a", { href: URL.createObjectURL(blob), download: name });
  a.click();
  URL.revokeObjectURL(a.href);
}

function plotView(spec) {
  let v = plotViews.get(spec.id);
  if (v) return v;
  const canvas = h("canvas", { width: 900, height: 300, role: "img", "aria-label": spec.title });
  const legend = h("ul", { className: "legend" });
  const note = h("p", { className: "note" });
  const csv = h("button", {
    type: "button",
    on: { click: () => {
      const s = plotSpecs.find((p) => p.id === spec.id);
      if (s) download(`lindhard-${s.id}.csv`, new Blob([toCsv(s.table)], { type: "text/csv" }));
    } },
  }, "Download CSV");
  const box = h("div", { className: "plot" }, h("h3", {}, spec.title), canvas, legend, note, h("div", { className: "tools" }, csv));
  $("plots").append(box);
  v = { canvas, legend, note };
  plotViews.set(spec.id, v);
  return v;
}

function legendSwatch(i) {
  const c = document.createElement("canvas");
  c.width = 34; c.height = 12;
  c.style.cssText = "width:34px;height:12px;border:0";
  const x = c.getContext("2d");
  x.strokeStyle = SPECIES_COLOURS[i % SPECIES_COLOURS.length]; x.lineWidth = 2;
  x.setLineDash(dashForSeries(i));
  x.beginPath(); x.moveTo(0, 6); x.lineTo(32, 6); x.stroke();
  return c;
}

function drawPlot(spec) {
  const { canvas, legend, note } = plotView(spec);
  const ctx = canvas.getContext("2d");
  const css = getComputedStyle(document.documentElement);
  const fg = css.getPropertyValue("--fg").trim() || "#000", bg = css.getPropertyValue("--bg").trim() || "#fff";
  const muted = css.getPropertyValue("--muted").trim() || "#666";
  const W = canvas.width, H = canvas.height, L = 64, R = 12, T = 12, B = 40;
  const ext = plotExtent(spec);
  const yTicks = niceTicks(0, ext.y1, 5);
  const y1 = Math.max(ext.y1, yTicks[yTicks.length - 1]);
  const sx = (x) => L + ((x - ext.x0) / (ext.x1 - ext.x0)) * (W - L - R);
  const sy = (y) => H - B - (y / y1) * (H - B - T);
  ctx.setLineDash([]);
  ctx.fillStyle = bg; ctx.fillRect(0, 0, W, H);
  ctx.font = "12px system-ui, sans-serif";
  ctx.lineWidth = 1; ctx.strokeStyle = muted; ctx.fillStyle = muted;
  for (const t of yTicks) {
    ctx.beginPath(); ctx.moveTo(L - 4, sy(t)); ctx.lineTo(L, sy(t)); ctx.stroke();
    ctx.textAlign = "right"; ctx.fillText(String(t), L - 7, sy(t) + 4);
  }
  ctx.beginPath(); ctx.moveTo(L, T); ctx.lineTo(L, H - B); ctx.lineTo(W - R, H - B); ctx.stroke();
  ctx.textAlign = "center";
  if (spec.kind === "bars") {
    const n = spec.categories.length, bw = (W - L - R) / Math.max(n, 1);
    spec.categories.forEach((c, i) => {
      ctx.fillStyle = muted;
      ctx.fillText(c, L + (i + 0.5) * bw, H - B + 16, bw - 4);
      ctx.fillStyle = SPECIES_COLOURS[i % SPECIES_COLOURS.length];
      const top = sy(spec.series[0].y[i]);
      ctx.fillRect(L + i * bw + bw * 0.15, top, bw * 0.7, H - B - top);
      ctx.strokeStyle = fg; ctx.strokeRect(L + i * bw + bw * 0.15, top, bw * 0.7, H - B - top);
    });
  } else {
    for (const t of niceTicks(ext.x0, ext.x1, 6)) {
      ctx.strokeStyle = muted; ctx.fillStyle = muted;
      ctx.beginPath(); ctx.moveTo(sx(t), H - B); ctx.lineTo(sx(t), H - B + 4); ctx.stroke();
      ctx.fillText(String(t), sx(t), H - B + 16);
    }
    ctx.lineWidth = 1.6;
    for (const s of spec.series) {
      ctx.strokeStyle = SPECIES_COLOURS[s.style % SPECIES_COLOURS.length];
      ctx.setLineDash(dashForSeries(s.style));
      ctx.beginPath();
      if (s.line) {
        s.y.forEach((y, i) => {
          const x = sx((s.bins.lo[i] + s.bins.hi[i]) / 2);
          if (i === 0) ctx.moveTo(x, sy(y)); else ctx.lineTo(x, sy(y));
        });
      } else {
        s.y.forEach((y, i) => {
          if (i === 0) ctx.moveTo(sx(s.bins.lo[i]), sy(y)); else ctx.lineTo(sx(s.bins.lo[i]), sy(y));
          ctx.lineTo(sx(s.bins.hi[i]), sy(y));
        });
      }
      ctx.stroke();
    }
    ctx.setLineDash([]);
  }
  ctx.fillStyle = fg; ctx.textAlign = "center";
  ctx.fillText(spec.xLabel, (L + W - R) / 2, H - 6);
  ctx.save(); ctx.translate(14, (T + H - B) / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText(spec.yLabel, 0, 0); ctx.restore();
  if (spec.empty) { ctx.fillStyle = muted; ctx.fillText("No data yet", (L + W - R) / 2, (T + H - B) / 2); }
  legend.replaceChildren(...(spec.kind === "hist"
    ? spec.series.map((s) => h("li", {}, legendSwatch(s.style), ` ${s.label}`)) : []));
  const flows = spec.series.filter((s) => s.bins && !s.overlay && (s.bins.underflow || s.bins.overflow))
    .map((s) => `${s.label}: ${s.bins.underflow} below, ${s.bins.overflow} above the plotted range`);
  note.textContent = [spec.note, ...flows].filter(Boolean).join(" ");
}

function drawPlots(report) {
  lastReport = report;
  if (report === null) { $("plots").replaceChildren(); plotViews.clear(); }
  plotSpecs = buildPlots(report, { layerNames: setup.layers.map((l) => l.name), symbolOfZ });
  for (const s of plotSpecs) drawPlot(s);
  $("caveats").replaceChildren(...caveatLines(report, setup.physics).map((t) => h("li", {}, t)));
  $("dl-json").disabled = report === null;
}

$("dl-json").addEventListener("click", () => {
  if (lastReport) download("lindhard-report.json", new Blob([reportJson(lastReport)], { type: "application/json" }));
});
$("dl-png").addEventListener("click", () => {
  $("traj").toBlob((b) => b && download("lindhard-trajectories.png", b), "image/png");
});

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
        drawPlots(m.report);
        histories.push(...m.trajectories);
        if (histories.length === m.trajectories.length) drawLegend(m.report);
        if (m.trajectories.length) drawTrajectories();
        break;
      case "finished":
        counters(m.done, m.total, m.report);
        drawPlots(m.report);
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
  worker.postMessage({ type: "init", toml: toToml(setup), ions: setup.run.ions, recordHistories: RECORDED_HISTORIES, maxPoints: MAX_POINTS_PER_TRACK });
  histories = [];
  counters(0, setup.run.ions, null);
  drawTrajectories();
  drawPlots(null);
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
  $("ions").max = String(Math.min(MAX_IONS, INTERACTIVE_MAX_IONS));
  $("cap-note").textContent = `The page runs up to ${INTERACTIVE_MAX_IONS} ions and records the paths of the first ${RECORDED_HISTORIES}, so a run stays responsive; the command-line tool has no such limit.`;
  symbolOfZ = (z) => elements.find((e) => e.z === z)?.symbol;
  $("traj-view").addEventListener("change", drawTrajectories);
  $("traj-cascades").addEventListener("change", drawTrajectories);
  drawTrajectories();
  drawPlots(null);
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
