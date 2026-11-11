// The simulator page: Setup (form <-> lindhard TOML) and a minimal Run view.
// The Worker (worker.js) owns the WASM session; this file only draws.

import init, { elements_json } from "./pkg/lindhard_wasm.js";
import { POTENTIALS, STOPPINGS, MAX_IONS, ModelError, PRESETS, el, defaultSetup, fromToml, toToml, validate } from "./toml_model.js";
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
    id: `tunit-${i}`, "aria-label": `thickness unit of layer ${i + 1}`,
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

  const cell = (cap, input) => h("div", { className: "cell" }, h("span", { className: "cap", "aria-hidden": "true" }, ...cap), input);
  const rows = l.elements.map((e, j) => {
    const row = table.get(e.symbol);
    const field = (key, tabKey, label, cap) => cell(cap, h("input", {
      type: "number", min: 0, step: "any", "aria-label": `${label} of ${e.symbol} in ${l.name}`,
      value: show(e[key]), placeholder: tabKey === null ? "0" : row?.[tabKey] != null ? String(row[tabKey]) : "required",
      on: { input: upd((t) => (e[key] = numOrNull(t.value))) },
    }));
    return h("div", { className: "el-row" },
      cell(["Element"], h("input", { value: e.symbol, size: 3, list: "symbols", "aria-label": `element ${j + 1} of ${l.name}`, on: { input: upd((t) => { e.symbol = t.value.trim(); renderLayers(); }) } })),
      cell(["Atom fraction"], h("input", { type: "number", min: 0, step: "any", value: String(e.fraction), "aria-label": `amount of ${e.symbol} in ${l.name}`, on: { input: upd((t) => (e.fraction = Number(t.value))) } })),
      field("eDEv", "default_displacement_ev", "displacement energy", edCap("d")),
      field("eBEv", null, "lattice binding energy", edCap("b")),
      field("eSEv", "default_surface_binding_ev", "surface binding energy", edCap("s")),
      h("button", { type: "button", className: "icon rm", title: `Remove ${e.symbol}`, "aria-label": `remove ${e.symbol} from ${l.name}`, on: { click: () => { l.elements.splice(j, 1); renderLayers(); refresh(); } } }, "×"),
    );
  });
  const head = h("div", { className: "el-head", "aria-hidden": "true" },
    h("span", {}, "Element"), h("span", {}, "Fraction"),
    h("span", {}, ...edCap("d")), h("span", {}, ...edCap("b")), h("span", {}, ...edCap("s")), h("span"));

  const iconBtn = (text, name, disabled, fn) => h("button", { type: "button", className: "icon", title: name, "aria-label": name, disabled, on: { click: fn } }, text);
  return h("div", { className: "layer" },
    h("div", { className: "layer-title" },
      h("h3", {}, `Layer ${i + 1}`),
      h("button", { type: "button", on: { click: () => { l.elements.push(el("Si")); renderLayers(); refresh(); } } }, "Add element"),
      iconBtn("↑", "Move up", i === 0, () => move(i, -1)),
      iconBtn("↓", "Move down", last, () => move(i, 1)),
      iconBtn("×", "Delete layer", false, () => { setup.layers.splice(i, 1); unitT.splice(i, 1); renderLayers(); refresh(); })),
    h("div", { className: "layer-head" },
      h("label", {}, "Name", h("input", { value: l.name, "aria-label": `name of layer ${i + 1}`, on: { input: upd((t) => (l.name = t.value.trim())) } })),
      h("label", {}, "Thickness", h("span", { className: "group" }, thick, unit)),
      h("label", { className: "inline" }, semi, "Semi-infinite (substrate)"),
      h("label", {}, "Density (g/cm³)", density),
    ),
    head,
    ...rows,
  );
}

// E_d / E_b / E_s with subscripts; the legend under the layers names them.
function edCap(sub) {
  return ["E", h("sub", {}, sub), " (eV)"];
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

function showPresetNote(key) {
  const note = PRESETS[key]?.note;
  $("preset-note").textContent = note ? `${key}: ${note}.` : "";
}

const CUSTOM = "Custom";

function addPreset(key) {
  const p = PRESETS[key === CUSTOM ? "Si" : key];
  const last = setup.layers[setup.layers.length - 1];
  if (last && last.thicknessNm == null) last.thicknessNm = 100; // a new layer goes below
  setup.layers.push({
    name: uniqueName(key === CUSTOM ? "layer" : key), thicknessNm: null, densityGcm3: p.density,
    elements: p.elements.map((e) => ({ ...e })),
  });
  showPresetNote(key);
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
  p.followRecoils = $("mode-full").checked;
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
  $("mode-full").checked = p.followRecoils;
  $("mode-quick").checked = !p.followRecoils;
  $("potential").value = p.potential;
  $("stopping").value = p.stopping;
  $("pcut").value = String(p.primaryCutoffEv);
  $("rcut").value = String(p.recoilCutoffEv);
  $("ions").value = String(setup.run.ions);
  $("seed").value = String(setup.run.seed);
  $("bin").value = String(setup.tally.depthBinNm);
  $("bins").value = String(setup.tally.depthBins);
}

// ---- layer stack drawing ---------------------------------------------------

const SVG_NS = "http://www.w3.org/2000/svg";
function svg(tag, attrs = {}, ...kids) {
  const e = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  e.append(...kids);
  return e;
}

const STACK = { x0: 64, top: 20, height: 80, minW: 62, finiteW: 210, substrateW: 74, beamLen: 50 };
const thicknessText = (nmValue) => (nmValue >= 1000 ? `${+(nmValue / 1000).toPrecision(3)} µm` : `${+nmValue.toPrecision(3)} nm`);

// Layers left to right from the surface at x0, beam arriving at the tilt angle
// from the normal (horizontal here). Widths follow thickness, with a minimum so
// a thin layer stays readable; then the drawing says it is not to scale.
function layerStackSvg() {
  const { x0, top, height, minW, finiteW, substrateW, beamLen } = STACK;
  const layers = setup.layers;
  const finite = layers.filter((l) => l.thicknessNm != null);
  const total = finite.reduce((a, l) => a + (Number.isFinite(l.thicknessNm) ? Math.max(l.thicknessNm, 0) : 0), 0);
  let notToScale = false;
  const widths = layers.map((l) => {
    if (l.thicknessNm == null) return substrateW;
    const t = Number.isFinite(l.thicknessNm) ? Math.max(l.thicknessNm, 0) : 0;
    const w = total > 0 ? (t / total) * finiteW : minW;
    if (w < minW) { notToScale = true; return minW; }
    return w;
  });
  const sumW = widths.reduce((a, b) => a + b, 0);
  const tilt = Number.isFinite(setup.beam.tiltDeg) ? Math.min(Math.max(setup.beam.tiltDeg, 0), 89.9) : 0;
  const rad = (tilt * Math.PI) / 180;
  const yMid = top + height / 2;
  const vbW = x0 + Math.max(sumW, 40) + 12;
  const text = layers.length === 0
    ? "No layers."
    : layers.map((l, i) => `${i + 1}: ${l.name || "unnamed"}, ${l.thicknessNm == null ? "substrate (semi-infinite)" : thicknessText(l.thicknessNm)}`).join("; ") + ".";
  const desc = `Beam at ${tilt}° from the surface normal, entering the surface first. Layers from the surface: ${text}${notToScale ? " Not to scale." : ""}`;

  const root = svg("svg", { viewBox: `0 0 ${vbW} ${top + height + 24}`, role: "img", "aria-labelledby": "stack-title", "aria-describedby": "stack-desc" },
    svg("title", { id: "stack-title" }, "Layer stack"),
    svg("desc", { id: "stack-desc" }, desc),
    svg("defs", {}, svg("marker", { id: "stack-arrow", viewBox: "0 0 8 8", refX: 7, refY: 4, markerWidth: 7, markerHeight: 7, orient: "auto" },
      svg("path", { d: "M0 0 L8 4 L0 8 z", class: "stack-beam" }))));
  root.append(svg("line", {
    x1: x0 - beamLen * Math.cos(rad), y1: yMid - beamLen * Math.sin(rad), x2: x0 - 1, y2: yMid,
    class: "stack-beam", "marker-end": "url(#stack-arrow)",
  }));
  let x = x0;
  layers.forEach((l, i) => {
    const w = widths[i], sub = l.thicknessNm == null;
    const edge = sub
      ? `L${x + w} ${top} ` + Array.from({ length: 6 }, (_, k) => `L${x + w - (k % 2 ? 0 : 7)} ${top + ((k + 1) * height) / 6}`).join(" ") + ` L${x} ${top + height} Z`
      : `L${x + w} ${top} L${x + w} ${top + height} L${x} ${top + height} Z`;
    root.append(svg("path", { d: `M${x} ${top} ${edge}`, class: `stack-block${i % 2 ? " alt" : ""}` }));
    const fit = Math.max(3, Math.floor(w / 6));
    const name = l.name || "unnamed";
    const short = name.length > fit ? `${name.slice(0, fit - 1)}…` : name;
    const label = svg("text", { x: x + w / 2 - (sub ? 3 : 0), y: yMid - 2, class: "stack-text" }, short);
    label.append(svg("title", {}, name));
    root.append(label);
    root.append(svg("text", { x: x + w / 2 - (sub ? 3 : 0), y: yMid + 13, class: "stack-text" }, sub ? "substrate" : thicknessText(l.thicknessNm)));
    x += w;
  });
  if (layers.length === 0) root.append(svg("text", { x: x0 + 20, y: yMid + 4, class: "stack-text", "text-anchor": "start" }, "No layers"));
  if (notToScale) root.append(svg("text", { x: vbW - 4, y: top + height + 16, class: "stack-text", style: "text-anchor:end" }, "not to scale"));
  return root;
}

let problems = [];
function refresh() {
  readScalars();
  const row = table.get(setup.beam.ion);
  $("mass").placeholder = row ? String(row.atomic_weight) : "";
  $("mass-note").textContent = "Leave the mass blank to use the element's standard atomic weight (the placeholder); enter a mass to simulate one isotope." +
    (row?.weight_is_mass_number
      ? ` ${row.symbol} has no standard atomic weight; the mass number of its longest-lived isotope is used.` : "");
  $("model-cite").textContent = `${CITE[setup.physics.potential]}. ${CITE[setup.physics.stopping]}.`;
  $("physics-summary").textContent = `${setup.physics.potential} potential, ${setup.physics.stopping} stopping. ${$("model-cite").textContent}`;
  $("stack").replaceChildren(layerStackSvg());
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
  $("r-erec").hidden = c.partition?.recoilsNotFollowed == null;
  $("c-erec").textContent = pct(c.partition?.recoilsNotFollowed);
  $("c-eoth").textContent = pct(c.partition?.other);
}

// ---- trajectory view -------------------------------------------------------

// Series colours are CSS custom properties (--s0..--s5, per colour scheme in
// simulator.html): Okabe-Ito, safe for the common colour-vision deficiencies;
// the second cue is the dash pattern (generation) and an end marker (species).
// `SPECIES_COLOURS` is re-read from the page at the start of every draw.
const SERIES_COUNT = 6;
let SPECIES_COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"];
function readPalette() {
  const css = getComputedStyle(document.documentElement);
  SPECIES_COLOURS = SPECIES_COLOURS.map((old, i) => css.getPropertyValue(`--s${i}`).trim() || old);
}

/**
 * Size a canvas bitmap from its CSS box and devicePixelRatio and scale the
 * context so callers draw in CSS pixels. Returns {w, h} in CSS pixels, or null
 * when the canvas has no size (hidden or collapsed); nothing is drawn then.
 */
function fitCanvas(canvas) {
  const w = canvas.clientWidth, h = canvas.clientHeight; // content box, border excluded
  if (w < 1 || h < 1) return null;
  const dpr = window.devicePixelRatio || 1;
  const bw = Math.round(w * dpr), bh = Math.round(h * dpr);
  if (canvas.width !== bw) canvas.width = bw;
  if (canvas.height !== bh) canvas.height = bh;
  canvas.getContext("2d").setTransform(dpr, 0, 0, dpr, 0, 0);
  return { w, h };
}

const EXPORT_W = 900, EXPORT_H = 360, EXPORT_SCALE = 2; // PNG export size in CSS px and pixel ratio
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
  readPalette();
  const cv = $("traj");
  const size = fitCanvas(cv);
  if (size) paintTrajectories(cv.getContext("2d"), size.w, size.h);
  updateTrajectoryNote();
}

const hiddenSpecies = new Set(); // species whose tracks the legend toggles have hidden

// Paint the trajectory view into a W x H box of CSS pixels. Draw order: layer
// bands (to the canvas edges), recoils (thin, translucent), primaries, then
// markers on the primaries' end points only.
function paintTrajectories(ctx, W, H) {
  const css = getComputedStyle(document.documentElement);
  const fg = css.getPropertyValue("--fg").trim() || "#000", bg = css.getPropertyValue("--bg").trim() || "#fff";
  const panel = css.getPropertyValue("--panel").trim() || "#eee", muted = css.getPropertyValue("--muted").trim() || "#666";
  const view = $("traj-view").value;
  const cascades = $("traj-cascades").checked;
  const segs = trajectorySegments(histories, { view, cascades, hidden: hiddenSpecies });
  const L = 50, R = 12, T = 20, B = 30; // margins hold the tick labels and axis titles
  const pw = Math.max(W - L - R, 1), ph = Math.max(H - T - B, 1);
  const bd = segmentBounds(segs, view, pw / ph);
  const sx = (u) => L + ((u - bd.u0) / (bd.u1 - bd.u0)) * pw;
  const sy = (v) => T + ((v - bd.v0) / (bd.v1 - bd.v0)) * ph;
  const clampX = (x) => Math.min(Math.max(x, 0), W);
  ctx.setLineDash([]);
  ctx.fillStyle = bg;
  ctx.fillRect(0, 0, W, H);
  ctx.font = "12px system-ui, sans-serif";
  ctx.lineWidth = 1;
  ctx.textAlign = "left";
  if (view === "side") {
    // Bands run the full height and out to the canvas edges: vacuum in front of
    // the surface, then the layers shaded alternately, names at the top.
    const spans = layerSpans(setup.layers);
    ctx.fillStyle = muted;
    if (spans.length && sx(0) >= 44) ctx.fillText("vacuum", 4, 14);
    spans.forEach((l, i) => {
      const x0 = clampX(sx(l.top)), x1 = clampX(sx(l.bottom));
      if (x1 <= x0) return;
      ctx.fillStyle = panel;
      ctx.fillRect(x0, 0, x1 - x0, H);
      if (i % 2) { ctx.globalAlpha = 0.5; ctx.fillStyle = bg; ctx.fillRect(x0, 0, x1 - x0, H); ctx.globalAlpha = 1; }
      ctx.strokeStyle = muted;
      ctx.beginPath(); ctx.moveTo(x0, 0); ctx.lineTo(x0, H); ctx.stroke();
      ctx.fillStyle = muted;
      ctx.fillText(l.name, x0 + 4, 14, Math.max(x1 - x0 - 8, 0));
    });
    const last = spans[spans.length - 1];
    if (last && Number.isFinite(last.bottom) && sx(last.bottom) < W - 44) ctx.fillText("vacuum", sx(last.bottom) + 4, 14);
  }
  // Axes with ticks, in nm.
  const [uName, vName] = view === "side" ? ["depth", "y"] : ["y", "z"];
  ctx.strokeStyle = muted; ctx.fillStyle = muted;
  ctx.beginPath(); ctx.moveTo(L, T); ctx.lineTo(L, T + ph); ctx.lineTo(L + pw, T + ph); ctx.stroke();
  ctx.textAlign = "center";
  for (const t of niceTicks(nm(bd.u0), nm(bd.u1), Math.max(Math.floor(pw / 80), 2))) {
    const x = sx(t * 1e-9);
    ctx.beginPath(); ctx.moveTo(x, T + ph); ctx.lineTo(x, T + ph + 4); ctx.stroke();
    ctx.fillText(String(t), x, T + ph + 16);
  }
  ctx.textAlign = "right";
  for (const t of niceTicks(nm(bd.v0), nm(bd.v1), Math.max(Math.floor(ph / 40), 2))) {
    const y = sy(t * 1e-9);
    ctx.beginPath(); ctx.moveTo(L - 4, y); ctx.lineTo(L, y); ctx.stroke();
    ctx.fillText(String(t), L - 7, y + 4);
  }
  ctx.textAlign = "center";
  ctx.fillText(`${uName} (nm)`, L + pw / 2, H - 4);
  ctx.save(); ctx.translate(11, T + ph / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(`${vName} (nm)`, 0, 0); ctx.restore();
  ctx.textAlign = "right";
  ctx.fillText(bd.lateralScale >= 0.999 ? "scale 1:1" : `lateral scale ×${bd.lateralScale.toPrecision(2)} (not 1:1)`, W - R, 14);
  ctx.textAlign = "left";
  // Tracks: one path per species and generation, clipped to the plot area.
  ctx.save();
  ctx.beginPath(); ctx.rect(L, T, pw, ph); ctx.clip();
  const groups = new Map();
  for (const s of segs) {
    const key = s.generation * 1000 + s.species;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(s);
  }
  // Later generations first, so the primaries are on top.
  for (const key of [...groups.keys()].sort((a, b) => b - a)) {
    const list = groups.get(key), generation = Math.floor(key / 1000), species = key % 1000;
    ctx.strokeStyle = SPECIES_COLOURS[species % SPECIES_COLOURS.length];
    ctx.globalAlpha = generation ? 0.45 : 1;
    ctx.lineWidth = generation ? 0.8 : 2;
    ctx.setLineDash(dashFor(generation));
    ctx.beginPath();
    for (const s of list) { ctx.moveTo(sx(s.a[0]), sy(s.a[1])); ctx.lineTo(sx(s.b[0]), sy(s.b[1])); }
    ctx.stroke();
  }
  ctx.globalAlpha = 1;
  ctx.setLineDash([]);
  ctx.lineWidth = 1;
  // Markers on the end points of primary tracks only.
  for (const hist of histories) for (const t of hist.tracks) {
    if (t.generation > 0 || hiddenSpecies.has(t.species)) continue;
    const q = t.points[t.points.length - 1];
    if (!q) continue;
    const [u, v] = view === "side" ? [q[0], q[1]] : [q[1], q[2]];
    ctx.fillStyle = SPECIES_COLOURS[t.species % SPECIES_COLOURS.length];
    ctx.strokeStyle = fg;
    marker(ctx, MARKERS[t.species % MARKERS.length], sx(u), sy(v), 4);
  }
  ctx.restore();
}

function updateTrajectoryNote() {
  const segs = trajectorySegments(histories, { view: $("traj-view").value, cascades: $("traj-cascades").checked, hidden: hiddenSpecies });
  const n = histories.length;
  $("traj-note").textContent = n === 0
    ? `No trajectories yet. The first ${RECORDED_HISTORIES} ion histories of a run are recorded.`
    : `${n} recorded histor${n === 1 ? "y" : "ies"} (the first ${RECORDED_HISTORIES} of the run), ${segs.length} segments.` +
      (anyTruncated(histories) ? " Some histories were cut at the point limit of each track, so their tails are missing." : "");
}

let legendReport = null;
function swatchCanvas() {
  const c = document.createElement("canvas");
  const dpr = window.devicePixelRatio || 1;
  c.width = Math.round(34 * dpr); c.height = Math.round(12 * dpr);
  c.style.cssText = "width:34px;height:12px;border:0";
  c.setAttribute("aria-hidden", "true");
  const x = c.getContext("2d");
  x.setTransform(dpr, 0, 0, dpr, 0, 0);
  return [c, x];
}

function drawLegend(report) {
  readPalette();
  legendReport = report;
  const names = speciesLabels(report, symbolOfZ);
  const item = (text, swatch) => h("li", {}, swatch, ` ${text}`);
  const line = (i, g) => {
    const [c, x] = swatchCanvas();
    x.strokeStyle = SPECIES_COLOURS[i % SPECIES_COLOURS.length]; x.lineWidth = 2;
    x.setLineDash(dashFor(g));
    x.beginPath(); x.moveTo(0, 6); x.lineTo(32, 6); x.stroke();
    return c;
  };
  // Each species has a checkbox: unchecking hides its tracks (the non-colour cue
  // between species); the dash patterns below tell the generations apart.
  const toggle = (n, i) => {
    const box = h("input", { type: "checkbox", checked: !hiddenSpecies.has(i), on: { change: (e) => {
      if (e.target.checked) hiddenSpecies.delete(i); else hiddenSpecies.add(i);
      drawTrajectories();
    } } });
    return h("li", {}, h("label", { className: "inline" }, box, " ", line(i, 0), ` ${n}`));
  };
  $("traj-legend").replaceChildren(
    ...names.map(toggle),
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
  const canvas = h("canvas", { role: "img", "aria-label": spec.title });
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
  resizeObserver?.observe(canvas);
  plotViews.set(spec.id, v);
  return v;
}

function legendSwatch(i, dash = i) {
  const [c, x] = swatchCanvas();
  x.strokeStyle = SPECIES_COLOURS[i % SPECIES_COLOURS.length]; x.lineWidth = 2;
  x.setLineDash(dashForSeries(dash));
  x.beginPath(); x.moveTo(0, 6); x.lineTo(32, 6); x.stroke();
  return c;
}

function drawPlot(spec) {
  const { canvas, legend, note } = plotView(spec);
  const size = fitCanvas(canvas);
  const ctx = canvas.getContext("2d");
  const css = getComputedStyle(document.documentElement);
  const fg = css.getPropertyValue("--fg").trim() || "#000", bg = css.getPropertyValue("--bg").trim() || "#fff";
  const muted = css.getPropertyValue("--muted").trim() || "#666";
  const W = size ? size.w : 900, H = size ? size.h : 300, L = 64, R = 28, T = 12, B = 40;
  const ext = plotExtent(spec);
  const yTicks = niceTicks(0, ext.y1, 5, { integer: true }); // counts are whole numbers
  const y1 = Math.max(ext.y1, yTicks[yTicks.length - 1]);
  const sx = (x) => L + ((x - ext.x0) / (ext.x1 - ext.x0)) * (W - L - R);
  const sy = (y) => H - B - (y / y1) * (H - B - T);
  ctx.setLineDash([]);
  ctx.fillStyle = bg; ctx.fillRect(0, 0, W, H);
  ctx.font = "12px system-ui, sans-serif";
  ctx.lineWidth = 1; ctx.strokeStyle = muted; ctx.fillStyle = muted;
  const gridLine = (x0, y0, x1, y1b) => {
    ctx.save(); ctx.globalAlpha = 0.2; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1b); ctx.stroke(); ctx.restore();
  };
  for (const t of yTicks) {
    gridLine(L, sy(t), W - R, sy(t));
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
      const sp = spec.speciesOf?.[i] ?? i;
      ctx.fillStyle = SPECIES_COLOURS[sp % SPECIES_COLOURS.length];
      const top = sy(spec.series[0].y[i]);
      ctx.fillRect(L + i * bw + bw * 0.15, top, bw * 0.7, H - B - top);
      ctx.strokeStyle = fg; ctx.strokeRect(L + i * bw + bw * 0.15, top, bw * 0.7, H - B - top);
      ctx.fillStyle = fg;
      ctx.fillText(String(spec.series[0].y[i]), L + (i + 0.5) * bw, top - 4, bw - 4);
    });
  } else {
    for (const t of niceTicks(ext.x0, ext.x1, 6)) {
      ctx.strokeStyle = muted; ctx.fillStyle = muted;
      gridLine(sx(t), T, sx(t), H - B);
      ctx.beginPath(); ctx.moveTo(sx(t), H - B); ctx.lineTo(sx(t), H - B + 4); ctx.stroke();
      ctx.fillText(String(t), sx(t), H - B + 16);
    }
    ctx.lineWidth = 1.6;
    // Clip series to the plot area: the x range is trimmed below the data's full range.
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, W - L - R, H - T - B); ctx.clip();
    for (const s of spec.series) {
      ctx.strokeStyle = SPECIES_COLOURS[s.style % SPECIES_COLOURS.length];
      ctx.setLineDash(dashForSeries(s.dashIndex ?? s.style));
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
    ctx.restore();
    ctx.setLineDash([]);
    for (const m of spec.markers ?? []) {
      if (!(m.x >= ext.x0 && m.x <= ext.x1)) continue;
      ctx.strokeStyle = fg; ctx.fillStyle = fg; ctx.lineWidth = 1; ctx.setLineDash([4, 3]);
      ctx.beginPath(); ctx.moveTo(sx(m.x), T); ctx.lineTo(sx(m.x), H - B); ctx.stroke();
      ctx.setLineDash([]); ctx.textAlign = "left";
      ctx.fillText(`${m.label} = ${Number(m.x.toPrecision(3))} nm`, sx(m.x) + 4, T + 10);
      ctx.textAlign = "center";
    }
  }
  ctx.fillStyle = fg; ctx.textAlign = "center";
  ctx.fillText(spec.xLabel, (L + W - R) / 2, H - 6);
  ctx.save(); ctx.translate(14, (T + H - B) / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText(spec.yLabel, 0, 0); ctx.restore();
  if (spec.empty) { ctx.fillStyle = muted; ctx.fillText("No data yet", (L + W - R) / 2, (T + H - B) / 2); }
  legend.replaceChildren(...(spec.kind === "hist"
    ? spec.series.map((s) => h("li", {}, legendSwatch(s.style, s.dashIndex ?? s.style), ` ${s.label}`)) : []));
  const flows = spec.series.filter((s) => s.bins && !s.overlay && (s.bins.underflow || s.bins.overflow))
    .map((s) => `${s.label}: ${s.bins.underflow} below, ${s.bins.overflow} above the plotted range`);
  const trim = ext.trimmed
    ? `The x axis is trimmed to ${ext.x1}; the tally extends to ${ext.full} and the CSV holds all of it.` : "";
  note.textContent = [spec.note, trim, ...flows].filter(Boolean).join(" ");
}

function drawPlots(report) {
  lastReport = report;
  if (report === null) { $("plots").replaceChildren(); plotViews.clear(); }
  readPalette();
  plotSpecs = buildPlots(report, { layerNames: setup.layers.map((l) => l.name), symbolOfZ });
  for (const s of plotSpecs) drawPlot(s);
  $("caveats").replaceChildren(...caveatLines(report, setup.physics).map((t) => h("li", {}, t)));
  $("dl-json").disabled = report === null;
}

$("dl-json").addEventListener("click", () => {
  if (lastReport) download("lindhard-report.json", new Blob([reportJson(lastReport)], { type: "application/json" }));
});
$("dl-png").addEventListener("click", () => {
  // Fixed size (1800 x 720 px), independent of the window and the pixel ratio.
  readPalette();
  const out = document.createElement("canvas");
  out.width = EXPORT_W * EXPORT_SCALE; out.height = EXPORT_H * EXPORT_SCALE;
  const ctx = out.getContext("2d");
  ctx.setTransform(EXPORT_SCALE, 0, 0, EXPORT_SCALE, 0, 0);
  paintTrajectories(ctx, EXPORT_W, EXPORT_H);
  out.toBlob((b) => b && download("lindhard-trajectories.png", b), "image/png");
});

// ---- redraw on resize, colour-scheme change and when a hidden view is shown ----

function redrawAll() {
  drawTrajectories();
  if (legendReport) drawLegend(legendReport);
  readPalette();
  for (const s of plotSpecs) drawPlot(s);
}

// Latest unpainted progress message; painted by the next animation frame.
let pendingProgress = null, progressDrawQueued = false;
// Paint the unpainted progress report, if any, and clear it. Called by the
// animation frame and, synchronously, by `finished` and `stopped` so the last
// batch (legend, trajectories, plots) is never dropped.
function flushProgress() {
  const p = pendingProgress;
  if (!p) return;
  pendingProgress = null;
  drawPlots(p.report);
  if (p.firstTrajectories) drawLegend(p.report);
  if (p.newTrajectories) drawTrajectories();
}
function queueProgressDraw() {
  if (progressDrawQueued) return;
  progressDrawQueued = true;
  requestAnimationFrame(() => {
    progressDrawQueued = false;
    flushProgress();
  });
}

let redrawQueued = false;
function queueRedraw() {
  if (redrawQueued) return;
  redrawQueued = true;
  requestAnimationFrame(() => { redrawQueued = false; redrawAll(); });
}
const resizeObserver = typeof ResizeObserver === "function" ? new ResizeObserver(queueRedraw) : null;
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", queueRedraw);
document.addEventListener("visibilitychange", () => { if (!document.hidden) queueRedraw(); });

if ("serviceWorker" in navigator) {
  // Offline support: see sw.js. Relative URL so the scope is this directory on
  // both hosts (site root and the GitHub Pages path prefix).
  navigator.serviceWorker.register("./sw.js").catch(() => {});
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
        // Counting and recording are cheap and happen for every message; the
        // drawing is coalesced to one per animation frame, so a slow draw
        // cannot let progress messages queue up behind it (which would keep
        // the counters moving after Pause).
        counters(m.done, m.total, m.report);
        histories.push(...m.trajectories);
        pendingProgress = {
          report: m.report,
          firstTrajectories: pendingProgress?.firstTrajectories || histories.length === m.trajectories.length,
          newTrajectories: (pendingProgress?.newTrajectories || 0) + m.trajectories.length,
        };
        queueProgressDraw();
        break;
      case "finished":
        flushProgress();
        counters(m.done, m.total, m.report);
        drawPlots(m.report);
        status("Done (± is one standard error).");
        endRun();
        break;
      case "stopped":
        flushProgress();
        status(`Stopped after ${m.done} of ${m.total} ions.`);
        endRun();
        break;
      case "error":
        status(m.message, true);
        endRun();
        break;
    }
  };
  worker.onerror = (e) => {
    status(e.message
      ? `The worker failed: ${e.message}`
      : "The worker script could not be loaded; if you are offline, reload the page once online.", true);
    endRun();
  };
  worker.postMessage({ type: "init", toml: toToml(setup), ions: setup.run.ions, recordHistories: RECORDED_HISTORIES, maxPoints: MAX_POINTS_PER_TRACK });
  histories = []; pendingProgress = null;
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
  $("preset").replaceChildren(...[CUSTOM, ...Object.keys(PRESETS)].map((p) => h("option", {}, p)));
  $("ions").max = String(Math.min(MAX_IONS, INTERACTIVE_MAX_IONS));
  $("cap-note").textContent = `The page runs up to ${INTERACTIVE_MAX_IONS} ions and records the paths of the first ${RECORDED_HISTORIES}, so a run stays responsive; the command-line tool has no such limit.`;
  symbolOfZ = (z) => elements.find((e) => e.z === z)?.symbol;
  resizeObserver?.observe($("traj"));
  $("traj-view").addEventListener("change", drawTrajectories);
  $("traj-cascades").addEventListener("change", drawTrajectories);
  drawTrajectories();
  drawPlots(null);
  $("preset").addEventListener("change", () => showPresetNote($("preset").value));
  showPresetNote($("preset").value);
  $("add-layer").addEventListener("click", () => addPreset($("preset").value));
  for (const id of ["ion", "mass", "energy", "energy-unit", "tilt", "mode-quick", "mode-full", "ions", "seed", "bin", "bins", "potential", "stopping", "pcut", "rcut"]) {
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
