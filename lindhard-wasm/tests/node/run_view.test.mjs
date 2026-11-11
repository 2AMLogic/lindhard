// The Run view's pure data-shaping (web/run_view.js) against hand-written
// report fixtures with the field names of lindhard::tally::ion::IonReport.
//
//   node --test lindhard-wasm/tests/node/run_view.test.mjs

import assert from "node:assert/strict";
import test from "node:test";
import {
  INTERACTIVE_MAX_IONS, RECORDED_HISTORIES, anyTruncated, dashFor, deriveCounters,
  energyPartition, layerSpans, roundToError, segmentBounds, speciesLabels, trajectorySegments, vacanciesPerIon,
} from "../../web/run_view.js";

const m = (mean, sd) => ({ n: 100, mean, std_dev: sd, skewness: 0, kurtosis: 3, mean_std_err: sd / 10, std_dev_std_err: sd / 14 });
const report = {
  histories: 100,
  species_z: [5, 14],
  budget: { incident: 10000, electronic_nonlocal: 6000, electronic_local: 500, lattice: 2500, surface_barrier: 10, backscattered: 400, sputtered: 20, transmitted: 0, rest: 570 },
  range: { stopped: 99, depth: m(3e-8, 1e-8), lateral_y: m(0, 5e-9), lateral_z: m(1e-10, 5e-9), radial: m(6e-9, 3e-9) },
  damage: { nrt: { nrt_displacements: 500 }, cascade: { vacancies: 300 } },
  escapes: { backscatter_coefficient: 0.04, transmission_coefficient: 0, energy_reflection_coefficient: 0.04, sputter_yield: 0.5 },
};

test("the interactive limits are below the session's 1e6 and a sane history count", () => {
  assert.ok(INTERACTIVE_MAX_IONS < 1_000_000);
  assert.ok(RECORDED_HISTORIES > 0 && RECORDED_HISTORIES <= INTERACTIVE_MAX_IONS);
});

test("energy partition adds up to the incident energy", () => {
  const p = energyPartition(report.budget);
  assert.ok(Math.abs(p.ionization - 0.65) < 1e-12);
  assert.ok(Math.abs(p.nuclear - 0.25) < 1e-12);
  assert.ok(Math.abs(p.escaped - 0.042) < 1e-12);
  assert.ok(Math.abs(p.ionization + p.nuclear + p.escaped + p.other - 1) < 1e-12);
  assert.equal(energyPartition({ incident: 0 }), null);
  assert.equal(energyPartition(undefined), null);
});

test("quick mode splits unfollowed-recoil energy out and both modes sum to 100 %", () => {
  const q = energyPartition(report.budget, false);
  assert.ok(Math.abs(q.recoilsNotFollowed - 0.057) < 1e-12);
  assert.ok(Math.abs(q.other - 0.001) < 1e-12);
  const sum = (p) => p.ionization + p.nuclear + p.escaped + p.other + (p.recoilsNotFollowed ?? 0);
  assert.ok(Math.abs(sum(q) - 1) < 1e-9);
  const f = energyPartition(report.budget, true);
  assert.equal(f.recoilsNotFollowed, null);
  assert.ok(Math.abs(f.other - 0.058) < 1e-12);
  assert.ok(Math.abs(sum(f) - 1) < 1e-9);
  const z = energyPartition({ ...report.budget, rest: 0, surface_barrier: 0, lattice: 3080 }, true);
  assert.ok(Math.abs(sum(z) - 1) < 1e-9);
  assert.equal(deriveCounters(report, { done: 1, total: 1, followRecoils: false }).partition.recoilsNotFollowed, 0.057);
  assert.equal(deriveCounters({ ...report, budget: null }, { done: 1, total: 1, followRecoils: false }).partition, null);
});

test("vacancies per ion use the cascade in full mode and NRT in quick mode", () => {
  assert.equal(vacanciesPerIon(report, true), 3);
  assert.equal(vacanciesPerIon(report, false), 5);
  assert.equal(vacanciesPerIon({ histories: 0 }, true), null);
});

test("counters from a report", () => {
  const c = deriveCounters(report, { done: 100, total: 400, followRecoils: true });
  assert.equal(c.fraction, 0.25);
  assert.equal(c.stopped, 99);
  assert.equal(c.depth.mean, 3e-8);
  assert.equal(c.depth.err, 1e-9);
  assert.equal(c.lateralZ.mean, 1e-10);
  assert.equal(c.radial.sd, 3e-9);
  assert.equal(c.backscatter, 0.04);
  assert.equal(c.sputterYield, 0.5);
  assert.equal(c.vacancies, 3);
  assert.equal(c.modeLabel, "full cascades");
});

test("counters before the first report are empty, not an error", () => {
  const c = deriveCounters(null, { done: 0, total: 10, followRecoils: false });
  assert.equal(c.depth, null);
  assert.equal(c.partition, null);
  assert.equal(c.backscatter, null);
  assert.equal(c.vacancies, null);
  assert.equal(deriveCounters(null, { done: 0, total: 0, followRecoils: false }).fraction, 0);
});

test("a run with no stopped ions has null moments", () => {
  const r = { ...report, range: { stopped: 0, depth: null, lateral_y: null, lateral_z: null, radial: null } };
  const c = deriveCounters(r, { done: 100, total: 100, followRecoils: false });
  assert.equal(c.depth, null);
  assert.equal(c.stopped, 0);
});

test("layer spans accumulate thickness and end in a substrate", () => {
  const s = layerSpans([{ name: "oxide", thicknessNm: 10 }, { name: "Si", thicknessNm: null }]);
  assert.deepEqual(s.map((x) => x.name), ["oxide", "Si"]);
  assert.equal(s[0].top, 0);
  assert.ok(Math.abs(s[0].bottom - 1e-8) < 1e-20);
  assert.ok(Math.abs(s[1].top - 1e-8) < 1e-20);
  assert.equal(s[1].bottom, Infinity);
});

const hist = (truncated = false) => ({
  index: 0, truncated,
  tracks: [
    { species: 0, generation: 0, parent: null, end: "stopped", points: [[0, 0, 0, 1e4], [1e-9, 2e-10, -1e-10, 5e3], [3e-9, 5e-10, 2e-10, 0]] },
    { species: 1, generation: 1, parent: 0, end: "stopped", points: [[1e-9, 2e-10, -1e-10, 40], [1.5e-9, 3e-10, 0, 0]] },
  ],
});

test("trajectory segments: side view uses depth and y, end view uses y and z", () => {
  const side = trajectorySegments([hist()], { view: "side" });
  assert.equal(side.length, 3);
  assert.deepEqual(side[0].a, [0, 0]);
  assert.deepEqual(side[0].b, [1e-9, 2e-10]);
  const end = trajectorySegments([hist()], { view: "end" });
  assert.deepEqual(end[0].b, [2e-10, -1e-10]);
});

test("recoil cascades are toggleable and primaries always shown", () => {
  const withC = trajectorySegments([hist()], { cascades: true });
  const without = trajectorySegments([hist()], { cascades: false });
  assert.equal(withC.length, 3);
  assert.equal(without.length, 2);
  assert.ok(without.every((s) => s.generation === 0));
});

test("a track with fewer than two points and no histories draw nothing", () => {
  assert.deepEqual(trajectorySegments([], {}), []);
  const h = { index: 0, truncated: false, tracks: [{ species: 0, generation: 0, parent: null, end: null, points: [[0, 0, 0, 1]] }] };
  assert.deepEqual(trajectorySegments([h], {}), []);
});

test("bounds: depth starts at the surface and lateral is symmetric", () => {
  const b = segmentBounds(trajectorySegments([hist()], { view: "side" }), "side");
  assert.equal(b.u0, 0);
  assert.ok(b.u1 >= 3e-9);
  assert.equal(b.v0, -b.v1);
  assert.ok(b.v1 >= 5e-10);
  const e = segmentBounds(trajectorySegments([hist()], { view: "end" }), "end");
  assert.equal(e.u0, -e.u1);
  assert.equal(e.v0, -e.v1);
  const empty = segmentBounds([], "side");
  assert.ok(empty.u1 > empty.u0 && empty.v1 > empty.v0);
});

test("bounds with a plot aspect are 1:1, or report the lateral scale", () => {
  const aspect = 3.6;
  const e = segmentBounds(trajectorySegments([hist()], { view: "end" }), "end", aspect);
  assert.ok(Math.abs((e.u1 - e.u0) / (e.v1 - e.v0) - aspect) < 1e-9);
  assert.equal(e.lateralScale, 1);
  // A narrow cascade fits: equal metres per pixel on both axes.
  const narrow = [{ a: [0, 0], b: [1e-7, 1e-9], species: 0, generation: 0 }];
  const s = segmentBounds(narrow, "side", aspect);
  assert.ok(Math.abs((s.u1 - s.u0) / (s.v1 - s.v0) - aspect) < 1e-9);
  assert.equal(s.lateralScale, 1);
  // A wide one does not: the lateral range holds it and the factor is reported.
  const wide = [{ a: [0, 0], b: [1e-8, 4e-8], species: 0, generation: 0 }];
  const w = segmentBounds(wide, "side", aspect);
  assert.ok(w.v1 >= 4e-8 && w.lateralScale < 1);
  assert.ok(Math.abs(w.lateralScale - ((w.u1 - w.u0) / aspect / 2) / w.v1) < 1e-12);
});

test("hidden species leave no segments", () => {
  assert.equal(trajectorySegments([hist()], { hidden: new Set([0]) }).filter((s) => s.species === 0).length, 0);
  assert.equal(trajectorySegments([hist()], { hidden: new Set() }).length, 3);
});

test("generations get different dash patterns (a non-colour cue)", () => {
  assert.notDeepEqual(dashFor(0), dashFor(1));
  assert.notDeepEqual(dashFor(1), dashFor(2));
  assert.deepEqual(dashFor(7), dashFor(2));
});

test("truncation flag and species labels", () => {
  assert.equal(anyTruncated([hist(), hist(true)]), true);
  assert.equal(anyTruncated([hist()]), false);
  const names = { 5: "B", 14: "Si" };
  assert.deepEqual(speciesLabels(report, (z) => names[z]), ["B (beam)", "Si"]);
  assert.deepEqual(speciesLabels(report, () => undefined), ["Z=5 (beam)", "Z=14"]);
});

test("a value and its error share a decimal place, the error to two figures", () => {
  assert.deepEqual(roundToError(46.38, 0.999), { value: "46.4", err: "1.0" });
  assert.deepEqual(roundToError(46.38, 0.131), { value: "46.38", err: "0.13" });
  assert.deepEqual(roundToError(56789.4, 1234), { value: "56800", err: "1200" });
  assert.deepEqual(roundToError(100, 0.1234), { value: "100.00", err: "0.12" });
  assert.deepEqual(roundToError(5, 0), { value: "5", err: "0" });
  assert.deepEqual(roundToError(3, null), { value: "3", err: null });
  assert.deepEqual(roundToError(NaN, 1), { value: "–", err: null });
});
