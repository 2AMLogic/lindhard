import test from "node:test";
import assert from "node:assert/strict";
import {
  buildPlots, caveatLines, csvField, dashForSeries, niceTicks, pearsonExpected, pearsonPdf, plotExtent, reportJson, toCsv,
} from "../../web/plots.js";

const hist = (lo, hi, counts, underflow = 0, overflow = 0) => ({
  binning: { lo, hi, bins: counts.length }, counts, underflow, overflow,
});
const face = (counts) => ({
  count: counts.reduce((a, b) => a + b, 0), per_ion: 0.1, energy_ev: 5, mean_energy_ev: 2.5,
  energy_histogram: hist(0, 100, counts), polar_histogram: hist(0, Math.PI / 2, counts),
});
const empty = () => hist(0, 1, [0, 0]);

// A Gaussian-like Pearson IV (nu = 0, m large) so the overlay is checkable.
const pearson = { mean: 2e-8, std_dev: 1e-8, skewness: 0, kurtosis: 3, m: 50, nu: 0, a: 7e-8, lambda: 2e-8, log_norm: 0 };

function report(extra = {}) {
  return {
    histories: 100,
    species_z: [5, 14],
    range: {
      stopped: 90,
      depth_histogram: hist(0, 4e-8, [10, 40, 30, 10], 2, 3),
      lateral_y_histogram: hist(-2e-8, 2e-8, [1, 2, 2, 1]),
      lateral_z_histogram: hist(-2e-8, 2e-8, [2, 1, 1, 2]),
      radial_histogram: hist(0, 2e-8, [3, 3]),
      layers: [
        { layer: 0, stopped: 60, depth_histogram: hist(0, 4e-8, [10, 30, 0, 0]) },
        { layer: 1, stopped: 30, depth_histogram: hist(0, 4e-8, [0, 10, 30, 10]) },
      ],
      pearson_iv: null,
      ...extra.range,
    },
    damage: { vacancy_histogram: hist(0, 4e-8, [1, 2, 3, 4]), interstitial_histogram: empty(), replacement_histogram: empty() },
    escapes: {
      species: [
        { species: 0, z: 5, beam: true, front: face([2, 1]), back: face([0, 0]) },
        { species: 1, z: 14, beam: false, front: face([1, 0]), back: face([0, 0]) },
      ],
    },
    ...extra.top,
  };
}
const names = { 5: "B", 14: "Si" };
const opts = { layerNames: ["SiO2", "Si"], symbolOfZ: (z) => names[z] };
const get = (plots, id) => plots.find((p) => p.id === id);

test("six plots, in a stable order, from one report", () => {
  assert.deepEqual(buildPlots(report(), opts).map((p) => p.id), ["range", "lateral", "damage", "escape-species", "escape-energy", "escape-polar"]);
  assert.deepEqual(buildPlots(null), []);
});

test("range plot: total plus one series per layer, in nm, with the layer names", () => {
  const p = get(buildPlots(report(), opts), "range");
  assert.deepEqual(p.series.map((s) => s.label), ["all stopped ions", "SiO2", "Si"]);
  assert.deepEqual(p.series[0].bins.lo, [0, 10, 20, 30]);
  assert.deepEqual(p.series[0].bins.hi, [10, 20, 30, 40]);
  assert.equal(p.series[0].bins.total, 95);
  assert.match(p.note, /No Pearson IV overlay/);
});

test("a single layer draws only the total", () => {
  const r = report();
  r.range.layers = [r.range.layers[0]];
  assert.equal(get(buildPlots(r, opts), "range").series.length, 1);
});

test("Pearson IV overlay is added when present, scaled to the histogram total", () => {
  const p = get(buildPlots(report({ range: { pearson_iv: pearson } }), opts), "range");
  const fit = p.series.at(-1);
  assert.equal(fit.overlay, true);
  assert.equal(fit.y.length, 4);
  assert.ok(fit.y.every(Number.isFinite));
  // symmetric about lambda = 20 nm: bins 1 and 2 are mirror images
  assert.ok(Math.abs(fit.y[1] - fit.y[2]) < 1e-9 * fit.y[1]);
  assert.match(p.note, /Pearson IV density/);
});

test("pearson pdf peaks at lambda for nu = 0 and a non-finite fit gives no overlay", () => {
  assert.ok(pearsonPdf(pearson, 2e-8) > pearsonPdf(pearson, 3e-8));
  assert.equal(pearsonExpected({ ...pearson, a: 0 }, hist(0, 1, [1])), null);
  assert.equal(pearsonExpected(null, hist(0, 1, [1])), null);
});

test("no pearson_iv but an error explains the missing overlay", () => {
  const p = get(buildPlots(report({ range: { pearson_iv_error: { OutsideTypeIv: {} } } }), opts), "range");
  assert.match(p.note, /outside the type IV region/);
  assert.equal(p.series.some((s) => s.overlay), false);
});

test("lateral and radial series keep their own binning", () => {
  const p = get(buildPlots(report(), opts), "lateral");
  assert.deepEqual(p.series.map((s) => s.label), ["lateral y", "lateral z", "radial"]);
  assert.equal(p.series[0].bins.lo[0], -20);
  assert.equal(p.series[2].bins.lo[0], 0);
});

test("CSV is long format with under- and overflow rows, and quotes correctly", () => {
  const p = get(buildPlots(report(), opts), "range");
  const lines = toCsv(p.table).split("\r\n");
  assert.equal(lines[0], "series,depth_nm_lo,depth_nm_hi,count");
  assert.equal(lines[1], "all stopped ions,0,10,10");
  assert.ok(lines.includes("all stopped ions,underflow,,2"));
  assert.ok(lines.includes("all stopped ions,overflow,,3"));
  assert.equal(lines.at(-1), "");
  assert.equal(csvField('a,"b"'), '"a,""b"""');
  assert.equal(csvField(null), "");
  assert.equal(csvField(NaN), "");
  assert.equal(csvField(0.1), "0.1");
});

test("the overlay rows carry no under/overflow rows", () => {
  const p = get(buildPlots(report({ range: { pearson_iv: pearson } }), opts), "range");
  const csv = toCsv(p.table);
  assert.match(csv, /Pearson IV fit \(expected counts\),0,10,/);
  assert.doesNotMatch(csv, /Pearson IV fit \(expected counts\),underflow/);
});

test("escape species: one bar per species and face, labelled by element", () => {
  const p = get(buildPlots(report(), opts), "escape-species");
  assert.equal(p.kind, "bars");
  assert.deepEqual(p.categories, ["B (beam), front", "B (beam), back", "Si, front", "Si, back"]);
  assert.deepEqual(p.series[0].y, [3, 0, 1, 0]);
  assert.equal(toCsv(p.table).split("\r\n")[1], '"B (beam), front",3,0.1,2.5');
});

test("escape spectra: only faces with escapes, angle in degrees", () => {
  const plots = buildPlots(report(), opts);
  const en = get(plots, "escape-energy");
  assert.deepEqual(en.series.map((s) => s.label), ["B (beam), front", "Si, front"]);
  assert.deepEqual(en.series.map((s) => s.style), [0, 1]);
  const pol = get(plots, "escape-polar");
  assert.ok(Math.abs(pol.series[0].bins.hi.at(-1) - 90) < 1e-9);
});

test("empty cases draw without error: zero stopped ions, thin film, no escapes", () => {
  const r = report();
  r.range.depth_histogram = hist(0, 4e-8, [0, 0, 0, 0]);
  r.range.layers = [];
  r.escapes.species = [];
  const plots = buildPlots(r, opts);
  assert.equal(get(plots, "range").empty, true);
  assert.equal(get(plots, "escape-species").empty, true);
  assert.equal(get(plots, "escape-energy").empty, true);
  for (const p of plots) {
    const ext = plotExtent(p);
    assert.ok(ext.x1 > ext.x0 && ext.y1 > 0, p.id);
    toCsv(p.table);
  }
});

test("all ions transmitted: back-face escapes are drawn, the range plot is empty", () => {
  const r = report();
  r.range.depth_histogram = hist(0, 4e-8, [0, 0, 0, 0]);
  r.range.layers = [];
  r.escapes.species[0].back = face([5, 5]);
  const plots = buildPlots(r, opts);
  assert.equal(get(plots, "range").empty, true);
  assert.ok(get(plots, "escape-energy").series.some((s) => s.label === "B (beam), back"));
});

test("missing sections (before the first report fields) do not throw", () => {
  const plots = buildPlots({ histories: 0 }, {});
  for (const p of plots) { plotExtent(p); toCsv(p.table); }
});

test("series styles give distinct dash patterns", () => {
  assert.notDeepEqual(dashForSeries(0), dashForSeries(1));
  assert.deepEqual(dashForSeries(4), dashForSeries(0));
});

test("nice ticks cover the range with 1, 2, 5 steps", () => {
  assert.deepEqual(niceTicks(0, 100, 5), [0, 20, 40, 60, 80, 100]);
  assert.deepEqual(niceTicks(0, 1, 4), [0, 0.2, 0.4, 0.6, 0.8, 1]);
  assert.deepEqual(niceTicks(3, 3, 5), [3]);
});

test("report JSON round-trips and caveats name the models and the validation doc", () => {
  const r = report();
  assert.deepEqual(JSON.parse(reportJson(r)), r);
  const c = caveatLines(r, { potential: "zbl", stopping: "lindhard-scharff", followRecoils: false }).join("\n");
  assert.match(c, /100 ion histories/);
  assert.match(c, /zbl potential/);
  assert.match(c, /quick mode/);
  assert.match(c, /docs\/validation\.md/);
  assert.match(caveatLines(null, null)[0], /no ions run yet/);
});
