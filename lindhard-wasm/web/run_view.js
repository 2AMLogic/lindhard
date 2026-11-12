// Pure data-shaping for the simulator's Run view: no DOM, so the Node tests
// cover it. `simulator.js` draws what these functions return.
//
// Input is the `IonReport` JSON of `Simulation.report_json()` (field names as
// in lindhard/src/tally/ion.rs) and the trajectory histories of
// `take_trajectories_json()` (lindhard-wasm/src/session.rs). Positions are in
// metres, x the depth.

// Interactive limits. The session itself accepts up to 1e6 ions; a browser run
// that must stay responsive and keep a legible trajectory view does not.
export const INTERACTIVE_MAX_IONS = 100000;
export const RECORDED_HISTORIES = 50; // the first histories whose tracks are kept
export const MAX_POINTS_PER_TRACK = 1000;

const finite = (x) => (typeof x === "number" && Number.isFinite(x) ? x : null);

// A moment summary as {mean, err, sd} in metres, or null before any sample.
function moment(m) {
  if (!m) return null;
  return { mean: m.mean, err: m.mean_std_err, sd: m.std_dev, sdErr: m.std_dev_std_err, n: m.n };
}

// Energy partition as fractions of the incident energy. Ionization is the
// electronic loss (nonlocal plus local); nuclear is the energy left in the
// lattice (`budget.lattice`); escaped is what leaves through the faces; other
// is the surface-barrier work and the kinetic energy of particles at rest.
// In quick mode (`followRecoils` false) `budget.rest` is almost entirely the
// kinetic energy of displaced recoils that were handed on but not followed, so
// it is reported apart as `recoilsNotFollowed` (null in full-cascade mode,
// where `rest` is cascade particles coming to rest and stays in `other`).
// The five parts sum to 1.
export function energyPartition(budget, followRecoils = true) {
  const inc = finite(budget?.incident);
  if (!inc || inc <= 0) return null;
  const ionization = (budget.electronic_nonlocal ?? 0) + (budget.electronic_local ?? 0);
  const nuclear = budget.lattice ?? 0;
  const escaped = (budget.backscattered ?? 0) + (budget.sputtered ?? 0) + (budget.transmitted ?? 0) + (budget.lateral ?? 0);
  const rest = budget.rest ?? 0;
  const other = (budget.surface_barrier ?? 0) + (followRecoils ? rest : 0);
  return {
    ionization: ionization / inc,
    nuclear: nuclear / inc,
    escaped: escaped / inc,
    other: other / inc,
    recoilsNotFollowed: followRecoils ? null : rest / inc,
  };
}

// Vacancies per ion. In full-cascade mode the followed cascade's count; in
// quick mode the NRT displacement estimate (the report carries no cascade).
export function vacanciesPerIon(report, followRecoils) {
  const n = report?.histories;
  if (!n) return null;
  const d = report.damage;
  const total = followRecoils ? d?.cascade?.vacancies : d?.nrt?.nrt_displacements;
  return finite(total) === null ? null : total / n;
}

// Everything the counters table shows, from one report. Absent values are null.
export function deriveCounters(report, { done, total, followRecoils }) {
  const r = report?.range;
  const e = report?.escapes;
  return {
    done,
    total,
    fraction: total ? done / total : 0,
    stopped: r ? r.stopped : null,
    depth: moment(r?.depth),
    lateralY: moment(r?.lateral_y),
    lateralZ: moment(r?.lateral_z),
    radial: moment(r?.radial),
    backscatter: e ? e.backscatter_coefficient : null,
    transmission: e ? e.transmission_coefficient : null,
    sputterYield: e ? e.sputter_yield : null,
    vacancies: vacanciesPerIon(report, followRecoils),
    partition: energyPartition(report?.budget, followRecoils),
    modeLabel: followRecoils ? "full cascades" : "quick (NRT estimate)",
  };
}

// Depths (m) of the interior layer boundaries, from thicknesses in nm with
// null for the substrate. Returns [{name, top, bottom}] with bottom = Infinity
// for the substrate.
export function layerSpans(layers) {
  const spans = [];
  let top = 0;
  for (const l of layers) {
    const bottom = l.thicknessNm == null ? Infinity : top + l.thicknessNm * 1e-9;
    spans.push({ name: l.name, top, bottom });
    top = bottom;
  }
  return spans;
}

// Non-colour cue: species picks a colour index, generation a dash pattern.
export const DASHES = [[], [6, 3], [2, 3]];
export const dashFor = (generation) => DASHES[Math.min(generation, DASHES.length - 1)];

// Segments to draw, in metres, as {a: [u, v], b: [u, v], species, generation}
// for a view: "side" is depth across, lateral y down; "end" is y across, z down;
// "3d" keeps all three coordinates, a and b being [x, y, z].
// Cascade recoils (generation > 0) are left out unless `cascades`; species in
// the `hidden` set are left out whatever their generation.
export function trajectorySegments(histories, { view = "side", cascades = true, hidden = new Set() } = {}) {
  const [iu, iv] = view === "end" ? [1, 2] : [0, 1];
  const pick = view === "3d" ? (p) => [p[0], p[1], p[2]] : (p) => [p[iu], p[iv]];
  const segs = [];
  for (const h of histories) {
    for (const t of h.tracks) {
      if (t.generation > 0 && !cascades) continue;
      if (hidden.has(t.species)) continue;
      for (let k = 1; k < t.points.length; k++) {
        const p = t.points[k - 1], q = t.points[k];
        segs.push({ a: pick(p), b: pick(q), species: t.species, generation: t.generation });
      }
    }
  }
  return segs;
}

// Extent {u0, u1, v0, v1, lateralScale} of a segment list, padded; the depth
// axis starts at the surface. A symmetric lateral range keeps the beam axis
// centred. `aspect` is the width/height of the plot area in pixels; with it the
// bounds are drawn at equal scale on both axes (1:1) and lateralScale is 1. The
// side view cannot always be 1:1: when the lateral spread does not fit, the
// lateral range is widened to hold it and lateralScale (< 1) is the lateral
// pixels per metre over the depth pixels per metre, to be printed on the plot.
export function segmentBounds(segs, view = "side", aspect = undefined) {
  const sized = Number.isFinite(aspect) && aspect > 0;
  let u0 = view === "side" ? 0 : Infinity, u1 = -Infinity, vmax = 0;
  for (const s of segs) {
    for (const p of [s.a, s.b]) {
      u0 = Math.min(u0, p[0]);
      u1 = Math.max(u1, p[0]);
      vmax = Math.max(vmax, Math.abs(p[1]));
    }
  }
  if (!segs.length) {
    // 1e-8 deep by 2w wide is 1:1 on a plot of this aspect when w = 5e-9 / aspect.
    const w = sized ? 5e-9 / aspect : 5e-9;
    return { u0: 0, u1: 1e-8, v0: -w, v1: w, lateralScale: 1 };
  }
  if (view === "end") {
    let m = 0;
    for (const s of segs) for (const p of [s.a, s.b]) m = Math.max(m, Math.abs(p[0]));
    const w = Math.max(m, vmax, 1e-10);
    const wu = sized ? w * Math.max(aspect, 1) : w, wv = sized ? w * Math.max(1 / aspect, 1) : w;
    return { u0: -wu, u1: wu, v0: -wv, v1: wv, lateralScale: 1 };
  }
  const span = Math.max(u1 - u0, 1e-10);
  const uEnd = u1 + span * 0.02;
  if (sized) {
    const fit = (uEnd - u0) / aspect / 2; // lateral half-range that makes the axes 1:1
    const w = Math.max(vmax, fit, 1e-10);
    return { u0, u1: uEnd, v0: -w, v1: w, lateralScale: fit / w };
  }
  const w = Math.max(vmax, span / 4, 1e-10);
  return { u0, u1: uEnd, v0: -w, v1: w, lateralScale: 1 };
}

// ---- 3D view ---------------------------------------------------------------

// The box {x0, x1, w} that holds a 3D segment list: depth from x0 (the surface
// at 0, or less if something left through it) to x1, and a lateral half-width
// w about the beam axis that is the same for y and z.
export function boxBounds(segs) {
  if (!segs.length) return { x0: 0, x1: 1e-8, w: 5e-9 };
  let x0 = 0, x1 = -Infinity, w = 0;
  for (const s of segs) {
    for (const p of [s.a, s.b]) {
      x0 = Math.min(x0, p[0]);
      x1 = Math.max(x1, p[0]);
      w = Math.max(w, Math.abs(p[1]), Math.abs(p[2]));
    }
  }
  const span = Math.max(x1 - x0, 1e-10);
  return { x0, x1: x1 + span * 0.02, w: Math.max(w * 1.05, span / 8, 1e-10) };
}

export const ORBIT_HOME = { yaw: -0.55, pitch: 0.3 };
const HALF_PI = Math.PI / 2;
export const clampPitch = (pitch) => Math.min(Math.max(pitch, -HALF_PI), HALF_PI);

// An orthographic projector of the box onto a pw x ph plot area: the same
// pixels per metre along all three axes, so the view is 1:1:1 at any angle.
// `yaw` turns the scene about the z axis (which points up the screen) and
// `pitch` tips the top towards the viewer (a view from above); at yaw 0 and pitch 0 depth runs
// to the right, z up and y away from the viewer. The scale depends on the
// pitch but not on the yaw, so the picture does not pulse while it turns.
// Returns {scale (px/m), project(p) -> [px, py, away]}: px and py from the
// top left of the plot area, `away` the distance behind the box centre in m.
export function orbitProjector({ yaw, pitch }, box, pw, ph) {
  const dx = box.x1 - box.x0, d = 2 * box.w;
  const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const across = Math.hypot(dx, d); // the widest the box can be, at any yaw
  const scale = Math.min(pw / across, ph / (d * Math.abs(cp) + across * Math.abs(sp)));
  const mx = (box.x0 + box.x1) / 2;
  return {
    scale,
    project([x, y, z]) {
      const X = x - mx;
      const right = X * cy - y * sy, away = X * sy + y * cy;
      const up = z * cp + away * sp;
      return [pw / 2 + scale * right, ph / 2 - scale * up, away * cp - z * sp];
    },
  };
}

// The 12 edges of the box as {axis, a, b}: axis 0 runs along the depth, 1
// along y, 2 along z; a and b are the end points in metres.
export function boxEdges(box) {
  const lo = [box.x0, -box.w, -box.w], hi = [box.x1, box.w, box.w];
  const edges = [];
  for (let axis = 0; axis < 3; axis++) {
    const [j, k] = [0, 1, 2].filter((i) => i !== axis);
    for (const vj of [lo[j], hi[j]]) for (const vk of [lo[k], hi[k]]) {
      const a = [0, 0, 0], b = [0, 0, 0];
      a[axis] = lo[axis]; b[axis] = hi[axis];
      a[j] = b[j] = vj; a[k] = b[k] = vk;
      edges.push({ axis, a, b });
    }
  }
  return edges;
}

// True if any recorded history was cut by the point cap.
export const anyTruncated = (histories) => histories.some((h) => h.truncated);

// The species legend: index -> label. Species 0 is the beam; the others are
// target elements, which the report names by atomic number (`species_z`).
export function speciesLabels(report, symbolOfZ) {
  const z = report?.species_z ?? [];
  return z.map((zi, i) => `${symbolOfZ(zi) ?? `Z=${zi}`}${i === 0 ? " (beam)" : ""}`);
}

// A value and its standard error as strings sharing one decimal place: the
// error to two significant figures, the value rounded to the same place
// (46.38 +- 0.999 -> "46.4", "1.0"). With no usable error the value gets
// three significant figures. Both are in whatever unit the caller passes.
export function roundToError(value, err) {
  if (finite(value) === null) return { value: "–", err: null };
  if (finite(err) === null || err < 0) return { value: Number(value.toPrecision(3)).toString(), err: null };
  if (err === 0) return { value: Number(value.toPrecision(3)).toString(), err: "0" };
  const exp = (x) => Math.floor(Math.log10(x));
  let step = 10 ** (exp(err) - 1);
  let rounded = Math.round(err / step) * step;
  step = 10 ** (exp(rounded) - 1); // 0.999 rounds up to 1.00: two figures are now 1.0
  rounded = Math.round(err / step) * step;
  const places = Math.max(0, -Math.round(Math.log10(step)));
  const v = Math.round(value / step) * step;
  return { value: v.toFixed(places), err: rounded.toFixed(places) };
}
