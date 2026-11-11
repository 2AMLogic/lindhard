// The simulator page, driven as a reader would: set up a beam and a target,
// run it in the Web Worker, and read the counters, plots and exports.
//
// The cases follow the live-verification list of issue #469. Ion counts are
// small so the suite stays quick; the physics is tested natively and under
// Node (lindhard-wasm/tests), and here the question is whether the page works.
import { readFileSync } from "node:fs";
import { test, expect, canvasColours, cellNumber } from "./fixtures.mjs";

const DONE = "Done (± is one standard error).";

async function open(page) {
  await page.goto("calculator/simulator.html");
  await expect(page.locator("#status")).toHaveText("Ready: B at 10 keV into Si.");
}

/** Fill the scalar fields of the Setup tab. */
async function configure(page, { ion, energyKeV, tilt, ions = 300, mode } = {}) {
  if (ion) await page.locator("#ion").fill(ion);
  if (energyKeV != null) await page.locator("#energy").fill(String(energyKeV));
  if (tilt != null) await page.locator("#tilt").fill(String(tilt));
  if (mode) await page.locator("#mode").selectOption(mode);
  await page.locator("#ions").fill(String(ions));
}

/**
 * Add a preset layer below the others. The element table has no displacement
 * energy for some elements (Si among them), and the page asks for one: give
 * the 15 eV the default layer uses.
 */
async function addPreset(page, preset) {
  await page.locator("#preset").selectOption(preset);
  await page.locator("#add-preset").click();
  for (const field of await page.locator('#layers input[aria-label^="displacement energy"][placeholder="required"]').all()) {
    if ((await field.inputValue()) === "") await field.fill("15");
  }
}

/** Replace the default Si substrate with a preset substrate. */
async function substrate(page, preset) {
  await page.getByRole("button", { name: "Delete layer" }).click();
  await addPreset(page, preset);
}

async function run(page, { timeout = 60_000 } = {}) {
  await expect(page.locator("#problems-box")).toBeHidden();
  await page.locator("#tab-run").click();
  await page.locator("#btn-run").click();
  await expect(page.locator("#status")).toHaveText(DONE, { timeout });
}

async function expectResults(page, ions) {
  await expect(page.locator("#c-done")).toHaveText(`${ions} / ${ions}`);
  await expect(page.locator("#progress")).toHaveJSProperty("value", 1);
  expect(await canvasColours(page.locator("#traj"))).toBeGreaterThan(3);
  await expect(page.locator("#traj-note")).toContainText("recorded histor");
  expect(await page.locator("#traj-legend li").count()).toBeGreaterThanOrEqual(4);
  const plots = page.locator("#plots .plot");
  await expect(plots).toHaveCount(6);
  for (const canvas of await plots.locator("canvas[role=img]").all()) {
    expect(await canvasColours(canvas)).toBeGreaterThan(2);
  }
  expect(await page.locator("#caveats li").count()).toBeGreaterThan(0);
  await expect(page.locator("#dl-json")).toBeEnabled();
  await expect(page.locator("#btn-run")).toBeEnabled();
  await expect(page.locator("#btn-stop")).toBeDisabled();
}

test("the page opens ready to run B at 10 keV into Si", async ({ page }) => {
  await open(page);
  await expect(page.locator("#ion")).toHaveValue("B");
  await expect(page.locator("#energy")).toHaveValue("10");
  await expect(page.locator("#energy-unit")).toHaveValue("keV");
  await expect(page.locator("#layers .layer")).toHaveCount(1);
  await expect(page.locator("#toml")).toHaveValue(/ion = "B"/);
  await expect(page.locator("#problems-box")).toBeHidden();
  await expect(page.locator("#cap-note")).toContainText("ions");
});

test("a short run paints the legend and trajectories as soon as it is done", async ({ page }) => {
  // The last progress report and `finished` can arrive in one animation
  // frame; the final draw must flush the pending report, not discard it.
  await open(page);
  await configure(page, { ion: "B", energyKeV: 10 });
  await run(page);
  expect(await page.locator("#traj-legend li").count()).toBeGreaterThanOrEqual(1);
  expect(await canvasColours(page.locator("#traj"))).toBeGreaterThan(3);
});

for (const [ion, energyKeV] of [["B", 10], ["P", 30], ["As", 50]]) {
  test(`${ion} into Si runs to the end`, async ({ page }) => {
    await open(page);
    await configure(page, { ion, energyKeV });
    await run(page);
    await expectResults(page, 300);
    expect(await cellNumber(page, "c-rp")).toBeGreaterThan(1);
    expect(await cellNumber(page, "c-drp")).toBeGreaterThan(0);
    expect(await cellNumber(page, "c-stopped")).toBeGreaterThan(200);
    await expect(page.locator("#c-vac-mode")).toContainText("(");
  });
}

test("As into SiO2 on Si runs with two layers built in the form", async ({ page }) => {
  await open(page);
  await configure(page, { ion: "As", energyKeV: 50 });
  await addPreset(page, "SiO2");
  await page.locator("#layers .layer").nth(1).getByRole("button", { name: "Move up" }).click();
  await page.locator("#thick-0").fill("10");
  await page.locator("#semi-1").check();
  await expect(page.locator("#layers .layer")).toHaveCount(2);
  await expect(page.locator("#toml")).toHaveValue(/substrate = "Si"[\s\S]*material = "SiO2"\s+thickness_nm = 10\.0/);
  await run(page);
  await expectResults(page, 300);
  expect(await cellNumber(page, "c-rp")).toBeGreaterThan(10);
});

test("the repository's two-layer example loads from its file and runs", async ({ page }) => {
  await open(page);
  const example = new URL("../../../examples/as_50keV_si_sio2.toml", import.meta.url);
  await page.locator("#load").setInputFiles({
    name: "as_50keV_si_sio2.toml", mimeType: "text/plain", buffer: readFileSync(example),
  });
  await expect(page.locator("#status")).toHaveText("Loaded as_50keV_si_sio2.toml.");
  await expect(page.locator("#ion")).toHaveValue("As");
  await expect(page.locator("#layers .layer")).toHaveCount(2);
  await page.locator("#ions").fill("300");
  await run(page);
  await expectResults(page, 300);
});

test("a file that is not a lindhard input is refused with a message", async ({ page }) => {
  await open(page);
  await page.locator("#load").setInputFiles({ name: "notes.toml", mimeType: "text/plain", buffer: Buffer.from("this is = = not toml") });
  await expect(page.locator("#status")).toHaveClass("error");
  await expect(page.locator("#status")).toContainText("notes.toml was not loaded");
  await expect(page.locator("#ion")).toHaveValue("B");
});

test("Ar into Cu with full cascades sputters", async ({ page }) => {
  await open(page);
  await configure(page, { ion: "Ar", energyKeV: 1, tilt: 0, mode: "full" });
  await substrate(page, "Cu");
  await run(page);
  await expectResults(page, 300);
  expect(await cellNumber(page, "c-sput")).toBeGreaterThan(0);
  expect(await cellNumber(page, "c-vac")).toBeGreaterThan(0);
  // Recoils are drawn: hiding the cascades leaves fewer segments.
  const segments = async () => Number((await page.locator("#traj-note").innerText()).match(/(\d+) segments/)[1]);
  const withCascades = await segments();
  await page.locator("#traj-cascades").uncheck();
  expect(await segments()).toBeLessThan(withCascades);
  await page.locator("#traj-view").selectOption("end");
  expect(await canvasColours(page.locator("#traj"))).toBeGreaterThan(3);
});

test("a thin film that transmits every ion draws without error", async ({ page }) => {
  await open(page);
  await configure(page, { ion: "H", energyKeV: 500, tilt: 0 });
  await page.locator("#semi-0").uncheck();
  await page.locator("#thick-0").fill("10");
  await run(page);
  await expect(page.locator("#c-done")).toHaveText("300 / 300");
  expect(await cellNumber(page, "c-trans")).toBeGreaterThan(0.99);
  await expect(page.locator("#c-stopped")).toHaveText("0");
  await expect(page.locator("#plots .plot")).toHaveCount(6);
});

test("oblique incidence at 60 degrees runs", async ({ page }) => {
  await open(page);
  await configure(page, { tilt: 60 });
  await run(page);
  await expectResults(page, 300);
});

test("the same seed gives the same counters, and another seed does not", async ({ page }) => {
  await open(page);
  await configure(page, {});
  const once = async () => {
    await run(page);
    const text = await page.locator("#counters").innerText();
    await page.locator("#tab-setup").click();
    return text;
  };
  const first = await once();
  expect(await once()).toBe(first);
  await page.locator("#seed").fill("7");
  expect(await once()).not.toBe(first);
});

test("pause holds the run, resume continues it, stop ends it", async ({ page }) => {
  await open(page);
  await configure(page, { ion: "Ar", energyKeV: 5, tilt: 0, mode: "full", ions: 100000 });
  await substrate(page, "Cu");
  await page.locator("#tab-run").click();
  await page.locator("#btn-run").click();
  await expect(page.locator("#status")).toHaveText("Running…");
  const done = async () => Number((await page.locator("#c-done").innerText()).split("/")[0]);
  await expect.poll(done).toBeGreaterThan(0);

  await page.locator("#btn-pause").click();
  await expect(page.locator("#status")).toHaveText("Paused.");
  // Pause takes effect within one step of the worker (about 50 ms).
  await page.waitForTimeout(500);
  const held = await done();
  await page.waitForTimeout(700);
  expect(await done()).toBe(held);
  expect(held).toBeLessThan(100000);
  await expect(page.locator("#btn-resume")).toBeEnabled();

  await page.locator("#btn-resume").click();
  await expect.poll(done).toBeGreaterThan(held);

  await page.locator("#btn-stop").click();
  await expect(page.locator("#status")).toContainText(/^Stopped after \d+ of 100000 ions\.$/);
  await expect(page.locator("#btn-run")).toBeEnabled();
});

test("the input file, the report, the plots and the picture download", async ({ page }) => {
  await open(page);
  await configure(page, {});
  const save = async (click) => {
    const [download] = await Promise.all([page.waitForEvent("download"), click()]);
    return { name: download.suggestedFilename(), body: readFileSync(await download.path()) };
  };

  const toml = await save(() => page.locator("#download").click());
  expect(toml.name).toBe("lindhard-input.toml");
  expect(toml.body.toString()).toBe(await page.locator("#toml").inputValue());

  await run(page);

  const report = await save(() => page.locator("#dl-json").click());
  expect(report.name).toBe("lindhard-report.json");
  expect(Object.keys(JSON.parse(report.body.toString())).length).toBeGreaterThan(0);

  const png = await save(() => page.locator("#dl-png").click());
  expect(png.name).toBe("lindhard-trajectories.png");
  expect(png.body.subarray(0, 8).toString("latin1")).toBe("\x89PNG\r\n\x1a\n");

  const buttons = page.locator("#plots .plot").getByRole("button", { name: "Download CSV" });
  await expect(buttons).toHaveCount(6);
  for (const button of await buttons.all()) {
    const csv = await save(() => button.click());
    expect(csv.name).toMatch(/^lindhard-[a-z-]+\.csv$/);
    const lines = csv.body.toString().trim().split("\n");
    expect(lines.length, csv.name).toBeGreaterThan(1);
    expect(lines[0], csv.name).toContain(",");
  }
});

test("a setup the page cannot run is explained and Run is disabled", async ({ page }) => {
  await open(page);
  await page.locator("#ions").fill("200000");
  await expect(page.locator("#problems-box")).toBeVisible();
  await expect(page.locator("#problems")).toContainText("run.ions");
  await page.locator("#tab-run").click();
  await expect(page.locator("#btn-run")).toBeDisabled();
  await page.locator("#tab-setup").click();
  await page.locator("#ions").fill("300");
  await page.locator("#ion").fill("Xx");
  await expect(page.locator("#problems-box")).toBeVisible();
  await page.locator("#ion").fill("B");
  await expect(page.locator("#problems-box")).toBeHidden();
});

test("the tabs and the run work from the keyboard", async ({ page }) => {
  await open(page);
  await configure(page, {});
  await page.locator("#tab-setup").focus();
  await page.keyboard.press("ArrowRight");
  await expect(page.locator("#tab-run")).toBeFocused();
  await expect(page.locator("#tab-run")).toHaveAttribute("aria-selected", "true");
  await expect(page.locator("#run")).toBeVisible();
  await expect(page.locator("#setup")).toBeHidden();
  await page.locator("#btn-run").focus();
  await page.keyboard.press("Enter");
  await expect(page.locator("#status")).toHaveText(DONE, { timeout: 60_000 });
  await page.locator("#tab-run").focus();
  await page.keyboard.press("ArrowLeft");
  await expect(page.locator("#setup")).toBeVisible();
});

for (const scheme of ["light", "dark"]) {
  test(`the ${scheme} theme reaches the page and the canvases`, async ({ page }) => {
    await page.emulateMedia({ colorScheme: scheme });
    await open(page);
    await configure(page, {});
    await run(page);
    await expectResults(page, 300);
    const luminance = (rgb) => rgb.reduce((a, b) => a + b, 0) / 3;
    const body = await page.evaluate(() => getComputedStyle(document.body).backgroundColor.match(/\d+/g).map(Number));
    // The corner of a plot is its background.
    const corner = await page.locator("#plots canvas[role=img]").first().evaluate((c) => Array.from(c.getContext("2d").getImageData(c.width - 2, 2, 1, 1).data.slice(0, 3)));
    if (scheme === "dark") {
      expect(luminance(body)).toBeLessThan(80);
      expect(luminance(corner)).toBeLessThan(80);
    } else {
      expect(luminance(body)).toBeGreaterThan(200);
      expect(luminance(corner)).toBeGreaterThan(200);
    }
  });
}

test("species are told apart without colour: markers and dash patterns are named", async ({ page }) => {
  await open(page);
  await configure(page, { ion: "As", energyKeV: 50, mode: "full" });
  await addPreset(page, "SiO2");
  await run(page);
  const legend = await page.locator("#traj-legend").innerText();
  for (const word of ["circle", "square", "solid", "dashed", "dotted"]) expect(legend).toContain(word);
});

test("the page is usable at 1280 px", async ({ page }) => {
  await page.setViewportSize({ width: 1280, height: 800 });
  await open(page);
  await configure(page, {});
  await run(page);
  await expectResults(page, 300);
  expect(await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth)).toBeLessThanOrEqual(0);
});

test("Setup and Run fit a phone screen", async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await open(page);
  const overflow = () => page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
  expect(await overflow()).toBeLessThanOrEqual(0);
  await configure(page, {});
  await run(page);
  expect(await overflow()).toBeLessThanOrEqual(0);
  await expect(page.locator("#traj")).toBeVisible();
});

// Issue #469 asks that the page work offline after its first load. The
// service worker (issue #653) precaches the files each run's Worker fetches.
test("a second run works after the network is gone", async ({ page, context }) => {
  await open(page);
  // The worker answers fetches only once it is active and controls the page.
  await page.evaluate(() => navigator.serviceWorker.ready);
  await page.waitForFunction(() => navigator.serviceWorker.controller !== null);
  await configure(page, {});
  await run(page);
  await context.setOffline(true);
  await page.locator("#btn-run").click();
  await expect(page.locator("#status")).toHaveText(DONE, { timeout: 60_000 });
  await expectResults(page, 300);
});

// Issue #674: canvases are sized from their CSS box and the pixel ratio, and
// redrawn when the box or the colour scheme changes.
const CANVASES = "#traj, #plots canvas[role=img]";

/** For each canvas: bitmap size, CSS content-box size and the ratio in force. */
const canvasSizes = (page) =>
  page.evaluate((sel) => Array.from(document.querySelectorAll(sel), (c) => ({
    w: c.width, h: c.height, cw: c.clientWidth, ch: c.clientHeight, dpr: window.devicePixelRatio,
    right: c.getBoundingClientRect().right,
  })), CANVASES);

test.describe("canvases at device pixel ratio 2", () => {
  test.use({ deviceScaleFactor: 2, viewport: { width: 390, height: 844 } });

  test("the bitmap is the CSS size times the pixel ratio and nothing overflows at 390 px", async ({ page }) => {
    await open(page);
    await configure(page, {});
    await run(page);
    await expectResults(page, 300);
    const sizes = await canvasSizes(page);
    expect(sizes.length).toBeGreaterThan(1);
    for (const c of sizes) {
      expect(c.dpr).toBe(2);
      expect(c.w).toBe(Math.round(c.cw * 2));
      expect(c.h).toBe(Math.round(c.ch * 2));
      expect(c.right).toBeLessThanOrEqual(390);
    }
  });
});

test("resizing the viewport redraws every canvas at the new size", async ({ page }) => {
  await page.setViewportSize({ width: 1280, height: 800 });
  await open(page);
  await configure(page, {});
  await run(page);
  await expectResults(page, 300);
  const wide = await canvasSizes(page);
  await page.setViewportSize({ width: 390, height: 844 });
  await expect.poll(async () => (await canvasSizes(page)).every((c) => c.w === c.cw * c.dpr)).toBe(true);
  const narrow = await canvasSizes(page);
  expect(narrow[0].w).toBeLessThan(wide[0].w);
  for (const c of narrow) expect(c.right).toBeLessThanOrEqual(390);
  expect(await canvasColours(page.locator("#traj"))).toBeGreaterThan(3);
});

test("switching the colour scheme redraws the canvases without a reload", async ({ page }) => {
  await page.emulateMedia({ colorScheme: "light" });
  await open(page);
  await configure(page, {});
  await run(page);
  await expectResults(page, 300);
  // The corner of a canvas is its background.
  const corner = () => page.locator(CANVASES).evaluateAll((cs) =>
    cs.map((c) => Array.from(c.getContext("2d").getImageData(c.width - 2, 2, 1, 1).data.slice(0, 3)).join(",")));
  const light = await corner();
  await page.emulateMedia({ colorScheme: "dark" });
  await expect.poll(async () => (await corner()).every((p, i) => p !== light[i])).toBe(true);
  await page.emulateMedia({ colorScheme: "light" });
  await expect.poll(async () => (await corner()).every((p, i) => p === light[i])).toBe(true);
});

test("the PNG export has a fixed size, whatever the window", async ({ page }) => {
  const size = async (w) => {
    await page.setViewportSize({ width: w, height: 800 });
    const [download] = await Promise.all([page.waitForEvent("download"), page.locator("#dl-png").click()]);
    const body = readFileSync(await download.path());
    return [body.readUInt32BE(16), body.readUInt32BE(20)]; // IHDR width, height
  };
  await page.setViewportSize({ width: 1280, height: 800 });
  await open(page);
  await configure(page, {});
  await run(page);
  const big = await size(1280);
  expect(big).toEqual([1800, 720]);
  expect(await size(390)).toEqual(big);
});
