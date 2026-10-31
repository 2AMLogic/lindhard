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
export function energyPartition(budget) {
  const inc = finite(budget?.incident);
  if (!inc || inc <= 0) return null;
  const ionization = (budget.electronic_nonlocal ?? 0) + (budget.electronic_local ?? 0);
  const nuclear = budget.lattice ?? 0;
  const escaped = (budget.backscattered ?? 0) + (budget.sputtered ?? 0) + (budget.transmitted ?? 0) + (budget.lateral ?? 0);
  const other = (budget.surface_barrier ?? 0) + (budget.rest ?? 0);
  return { ionization: ionization / inc, nuclear: nuclear / inc, escaped: escaped / inc, other: other / inc };
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
    partition: energyPartition(report?.budget),
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
// for a view: "side" is depth across, lateral y down; "end" is y across, z down.
// Cascade recoils (generation > 0) are left out unless `cascades`.
export function trajectorySegments(histories, { view = "side", cascades = true } = {}) {
  const [iu, iv] = view === "end" ? [1, 2] : [0, 1];
  const segs = [];
  for (const h of histories) {
    for (const t of h.tracks) {
      if (t.generation > 0 && !cascades) continue;
      for (let k = 1; k < t.points.length; k++) {
        const p = t.points[k - 1], q = t.points[k];
        segs.push({ a: [p[iu], p[iv]], b: [q[iu], q[iv]], species: t.species, generation: t.generation });
      }
    }
  }
  return segs;
}

// Extent {u0, u1, v0, v1} of a segment list, padded; the depth axis starts at
// the surface. A symmetric lateral range keeps the beam axis centred.
export function segmentBounds(segs, view = "side") {
  let u0 = view === "side" ? 0 : Infinity, u1 = -Infinity, vmax = 0;
  for (const s of segs) {
    for (const p of [s.a, s.b]) {
      u0 = Math.min(u0, p[0]);
      u1 = Math.max(u1, p[0]);
      vmax = Math.max(vmax, Math.abs(p[1]));
    }
  }
  if (!segs.length) return { u0: 0, u1: 1e-8, v0: -5e-9, v1: 5e-9 };
  if (view === "end") {
    let m = 0;
    for (const s of segs) for (const p of [s.a, s.b]) m = Math.max(m, Math.abs(p[0]));
    const w = Math.max(m, vmax, 1e-10);
    return { u0: -w, u1: w, v0: -w, v1: w };
  }
  const span = Math.max(u1 - u0, 1e-10);
  const w = Math.max(vmax, span / 4, 1e-10);
  return { u0, u1: u1 + span * 0.02, v0: -w, v1: w };
}

// True if any recorded history was cut by the point cap.
export const anyTruncated = (histories) => histories.some((h) => h.truncated);

// The species legend: index -> label. Species 0 is the beam; the others are
// target elements, which the report names by atomic number (`species_z`).
export function speciesLabels(report, symbolOfZ) {
  const z = report?.species_z ?? [];
  return z.map((zi, i) => `${symbolOfZ(zi) ?? `Z=${zi}`}${i === 0 ? " (beam)" : ""}`);
}
