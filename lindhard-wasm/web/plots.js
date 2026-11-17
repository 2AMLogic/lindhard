// Pure data-shaping for the simulator's distribution plots and exports: no DOM,
// so the Node tests cover it. `simulator.js` draws what these functions return.
//
// Input is the `IonReport` JSON of `Simulation.report_json()` (field names as
// in lindhard/src/tally/ion.rs and hist.rs). Depths and lateral distances are
// in metres and are shown in nm; energies are in eV; polar angles are in
// radians and are shown in degrees.

const NM = 1e9;
const DEG = 180 / Math.PI;

// Non-colour cue: series index picks a dash pattern as well as a colour.
export const PLOT_DASHES = [[], [6, 3], [2, 3], [8, 3, 2, 3]];
export const dashForSeries = (i) => PLOT_DASHES[i % PLOT_DASHES.length];

// One histogram as plain arrays. `scale` converts the binned unit to the
// displayed unit. Bins are [lo, hi) as in lindhard's `Binning`.
// Drop floating-point noise from an edge (12 significant digits).
const tidy = (x) => Number(x.toPrecision(12));

export function histBins(hist, scale = 1) {
  const b = hist?.binning;
  const counts = hist?.counts ?? [];
  if (!b || !counts.length) return null;
  const w = (b.hi - b.lo) / b.bins;
  const lo = [], hi = [];
  for (let i = 0; i < counts.length; i++) {
    lo.push(tidy((b.lo + i * w) * scale));
    hi.push(tidy((b.lo + (i + 1) * w) * scale));
  }
  const inRange = counts.reduce((s, c) => s + c, 0);
  const underflow = hist.underflow ?? 0, overflow = hist.overflow ?? 0;
  return { lo, hi, counts, underflow, overflow, inRange, total: inRange + underflow + overflow };
}

// Expected counts per bin of a Pearson IV density (Heinrich 2004; see
// lindhard/src/tally/pearson.rs) scaled to the histogram total, by the value of
// the density at each bin centre times the bin width. `p` carries the fit
// parameters of `RangeReport.pearson_iv`.
export function pearsonPdf(p, x) {
  const t = (x - p.lambda) / p.a;
  return Math.exp(p.log_norm - p.m * Math.log1p(t * t) - p.nu * Math.atan(t));
}

export function pearsonExpected(p, hist) {
  const b = hist?.binning;
  if (!p || !(p.a > 0) || !b || !hist.counts?.length) return null;
  const w = (b.hi - b.lo) / b.bins;
  const total = histBins(hist).total;
  const y = hist.counts.map((_, i) => total * w * pearsonPdf(p, b.lo + (i + 0.5) * w));
  return y.every(Number.isFinite) ? y : null;
}

const histSeries = (label, hist, scale, style, extra = {}) => {
  const bins = histBins(hist, scale);
  return bins ? { label, kind: "hist", style, bins, y: bins.counts, ...extra } : null;
};

const layerLabel = (layerNames, l) => layerNames?.[l.layer] || `layer ${l.layer + 1}`;

// The long-format table of a histogram plot: one row per bin per series, then
// the under- and overflow rows of each series (so nothing is dropped).
function histTable(series, xName) {
  const rows = [];
  for (const s of series) {
    s.bins.lo.forEach((lo, i) => rows.push([s.label, lo, s.bins.hi[i], s.y[i]]));
    if (s.overlay) continue;
    rows.push([s.label, "underflow", "", s.bins.underflow]);
    rows.push([s.label, "overflow", "", s.bins.overflow]);
  }
  return { columns: ["series", `${xName}_lo`, `${xName}_hi`, "count"], rows };
}

function histPlot(id, title, xLabel, xName, series, note, markers = []) {
  const s = series.filter(Boolean);
  const data = s.filter((x) => !x.overlay);
  return {
    id, title, xLabel, yLabel: "counts per bin", kind: "hist", series: s,
    empty: data.every((x) => x.bins.total === 0), note, markers,
    table: histTable(s, xName),
  };
}

// The plots, each with its series, labels and CSV table. `layerNames` names
// the layers by index; `symbolOfZ` maps an atomic number to a symbol.
export function buildPlots(report, { layerNames = [], symbolOfZ = () => undefined } = {}) {
  if (!report) return [];
  const r = report.range ?? {};
  const d = report.damage ?? {};
  const e = report.escapes ?? {};

  // 1. Ion range distribution.
  const range = [histSeries("all stopped ions", r.depth_histogram, NM, 0)];
  if ((r.layers?.length ?? 0) > 1) {
    r.layers.forEach((l, i) => range.push(histSeries(layerLabel(layerNames, l), l.depth_histogram, NM, i + 1)));
  }
  const fit = r.pearson_iv ? pearsonExpected(r.pearson_iv, r.depth_histogram) : null;
  if (fit && range[0]) {
    range.push({
      label: "Pearson IV fit (expected counts)", kind: "hist", overlay: true, line: true,
      style: range.length, bins: range[0].bins, y: fit,
    });
  }
  const rangeNote = r.pearson_iv
    ? "The line is the Pearson IV density of the four moments of the depth distribution (Hofker 1975; Heinrich 2004)."
    : r.pearson_iv_error
      ? "No Pearson IV overlay: the moments of this distribution are outside the type IV region or too few ions have stopped."
      : "No Pearson IV overlay yet.";

  // Rp, the mean projected range, marked on the range plot.
  const rp = r.depth?.mean;
  const rangeMarkers = Number.isFinite(rp) ? [{ x: rp * NM, label: "Rp" }] : [];

  // 2. Lateral (signed) and radial (non-negative) distributions, as two plots.
  const lateral = [
    histSeries("lateral y", r.lateral_y_histogram, NM, 0),
    histSeries("lateral z", r.lateral_z_histogram, NM, 1),
  ];
  const radial = [histSeries("radial", r.radial_histogram, NM, 0)];

  // 3. Damage versus depth: displacement events the engine counted. In quick
  // mode only the first recoil of each ion is created, so these are first-generation counts.
  const damage = [
    histSeries("vacancies", d.vacancy_histogram, NM, 0),
    histSeries("interstitials", d.interstitial_histogram, NM, 1),
    histSeries("replacements", d.replacement_histogram, NM, 2),
  ];

  // 4. Escapes: counts per species and face, then the energy and polar spectra.
  const names = (report.species_z ?? []).map((z, i) => `${symbolOfZ(z) ?? `Z=${z}`}${i === 0 ? " (beam)" : ""}`);
  const faces = [["front", "front face (backscattered, sputtered)"], ["back", "back face (transmitted)"]];
  const sp = e.species ?? [];
  const cats = [], vals = [], perIon = [], meanE = [], speciesOf = [];
  for (const s of sp) {
    for (const [face] of faces) {
      const f = s[face];
      if (!f) continue;
      cats.push(`${names[s.species] ?? `species ${s.species}`}, ${face}`);
      speciesOf.push(s.species);
      vals.push(f.count);
      perIon.push(f.per_ion);
      meanE.push(f.mean_energy_ev);
    }
  }
  const species = {
    id: "escape-species", title: "Escapes by species", xLabel: "species and face", yLabel: "particles escaped",
    kind: "bars", categories: cats, speciesOf, series: [{ label: "escaped particles", kind: "bars", style: 0, y: vals }],
    empty: vals.every((v) => !v), note: "Front: leaves through the entry face. Back: leaves through the far face of a finite stack.",
    table: {
      columns: ["species_face", "count", "per_ion", "mean_energy_ev"],
      rows: cats.map((c, i) => [c, vals[i], perIon[i], meanE[i]]),
    },
  };
  const spectrum = (key, scale) => {
    const out = [];
    let i = 0;
    for (const s of sp) {
      faces.forEach(([face, long], fi) => {
        const h = s[face]?.[key];
        const label = `${names[s.species] ?? `species ${s.species}`}, ${face}`;
        // Colour is the species (as in every other plot), dash is the face.
        const ser = histSeries(label, h, scale, s.species, { dashIndex: fi, title: long });
        if (ser && ser.bins.total > 0) { out.push(ser); i++; }
      });
    }
    return out;
  };
  const energy = histPlot("escape-energy", "Energy of escaping particles", "energy (eV)", "energy_ev",
    spectrum("energy_histogram", 1), "Only species and faces with escapes are drawn.");
  const polar = histPlot("escape-polar", "Polar angle of escaping particles", "polar angle (degrees)", "polar_deg",
    spectrum("polar_histogram", DEG), "Angle to the surface normal.");
  if (!energy.series.length) energy.empty = true;
  if (!polar.series.length) polar.empty = true;

  return [
    histPlot("range", "Ion range distribution", "depth (nm)", "depth_nm", range, rangeNote, rangeMarkers),
    histPlot("lateral", "Lateral distribution", "position along the axis (nm)", "distance_nm", lateral,
      "Lateral y and z are signed distances from the beam axis."),
    histPlot("radial", "Radial distribution", "distance from the beam axis (nm)", "distance_nm", radial,
      "Radial is the non-negative distance from the beam axis."),
    histPlot("damage", "Damage versus depth", "depth (nm)", "depth_nm", damage,
      "Displacement events counted by the engine. In quick mode only the first recoil of each ion is created (no cascade), so this is primary-collision damage; see the NRT estimate in the counters."),
    species, energy, polar,
  ];
}

// CSV text of a table: RFC 4180 quoting, "\r\n" line ends, numbers as written
// by JavaScript (round-trip precision).
export function csvField(v) {
  if (v === null || v === undefined) return "";
  const s = typeof v === "number" ? (Number.isFinite(v) ? String(v) : "") : String(v);
  return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

export function toCsv(table) {
  return [table.columns, ...table.rows].map((r) => r.map(csvField).join(",")).join("\r\n") + "\r\n";
}

// The full report as JSON text for download.
export const reportJson = (report) => JSON.stringify(report, null, 2) + "\n";

// "Nice" axis ticks: about `n` values at 1, 2 or 5 times a power of ten that
// cover [lo, hi]. Returns the tick values.
// With `integer`, the step is a whole number of at least 1 (count axes).
export function niceTicks(lo, hi, n = 5, { integer = false } = {}) {
  if (!(Number.isFinite(lo) && Number.isFinite(hi)) || hi <= lo) return [lo];
  const raw = (hi - lo) / Math.max(n, 1);
  const mag = 10 ** Math.floor(Math.log10(raw));
  const f = raw / mag;
  let step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
  if (integer) step = Math.max(1, Math.round(step));
  const ticks = [];
  for (let t = Math.ceil(lo / step - 1e-9) * step; t <= hi + step * 1e-9; t += step) {
    ticks.push(Number(t.toPrecision(12)));
  }
  return ticks;
}

// Extent of a plot's data: {x0, x1, y1}, y from zero. Falls back to a unit
// box for an empty plot so the axes still draw. A histogram's upper x limit is
// trimmed to the upper edge of its last non-empty bin plus a margin, rounded
// up to a tick; then `trimmed` is true and `full` is the untrimmed upper limit.
// Nothing is trimmed when a series has overflow counts (the data run past the
// plotted range) or when the plot is empty. The CSV always holds the full tally.
export function plotExtent(spec) {
  if (spec.kind === "bars") {
    return { x0: 0, x1: Math.max(spec.categories.length, 1), y1: Math.max(1, ...spec.series[0].y) };
  }
  let x0 = Infinity, x1 = -Infinity, y1 = 0, dataEnd = -Infinity, overflow = false;
  for (const s of spec.series) {
    x0 = Math.min(x0, s.bins.lo[0]);
    x1 = Math.max(x1, s.bins.hi[s.bins.hi.length - 1]);
    for (const v of s.y) y1 = Math.max(y1, v);
    if (s.overlay) continue;
    if (s.bins.overflow > 0) overflow = true;
    for (let i = s.y.length - 1; i >= 0; i--) {
      if (s.counts?.[i] > 0 || s.y[i] > 0) { dataEnd = Math.max(dataEnd, s.bins.hi[i]); break; }
    }
  }
  if (!(x1 > x0)) return { x0: 0, x1: 1, y1: 1, trimmed: false, full: 1 };
  const out = { x0, x1, y1: y1 > 0 ? y1 : 1, trimmed: false, full: x1 };
  if (overflow || !(dataEnd > x0) || dataEnd >= x1) return out;
  const target = Math.min(x1, dataEnd + 0.05 * (dataEnd - x0));
  const ticks = niceTicks(x0, target, 6);
  const step = ticks.length > 1 ? ticks[1] - ticks[0] : 0;
  const last = ticks[ticks.length - 1];
  const rounded = Number((last < target - step * 1e-9 && step > 0 ? last + step : last).toPrecision(12));
  if (rounded < x1) { out.x1 = rounded; out.trimmed = true; }
  return out;
}

// The caveat lines shown with the results: statistical error, models in use,
// validation status. `physics` is the setup's physics block.
export function caveatLines(report, physics) {
  const n = report?.histories ?? 0;
  const lines = [
    n > 0
      ? `Statistical error: ${n} ion histories so far. Errors shown are one standard error; histogram bins carry a counting error of about the square root of the count, so sparse bins are noisy.`
      : "Statistical error: no ions run yet; results carry the error of the finite sample.",
  ];
  if (physics) {
    lines.push(`Models in use: ${physics.potential} potential, ${physics.stopping} electronic stopping, ${physics.followRecoils ? "full recoil cascades" : "quick mode (NRT damage estimate, no cascades)"}.`);
  }
  lines.push("Validation status: see docs/validation.md for what has been compared with measurements and where the models are known to be weak.");
  return lines;
}
