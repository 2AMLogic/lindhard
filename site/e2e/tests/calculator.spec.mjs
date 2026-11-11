// The range calculator is the simulator page in its Quick range state
// (calculator/index.html): the engine loads in the browser and reproduces the
// native result recorded in lindhard-wasm/tests/b_5keV_si.golden.json.
import { readFileSync } from "node:fs";
import { test, expect, canvasColours } from "./fixtures.mjs";

const golden = JSON.parse(
  readFileSync(new URL("../../../lindhard-wasm/tests/b_5keV_si.golden.json", import.meta.url), "utf8"),
);
const DONE = "Done (± is one standard error).";

test.beforeEach(async ({ page }) => {
  await page.goto("calculator/index.html");
  // The example run on load finishes first.
  await expect(page.locator("#status")).toHaveText(DONE, { timeout: 60_000 });
});

/** Leading number of a tile such as "24.4 nm ± 0.13 nm". */
async function tileNumber(page, id) {
  const text = await page.locator(`#${id}`).innerText();
  const m = text.match(/-?\d+(\.\d+)?(e[+-]?\d+)?/i);
  expect(m, `#${id} holds a number, got "${text}"`).not.toBeNull();
  return Number(m[0]);
}

async function run(page) {
  await page.locator("#btn-run").click();
  await expect(page.locator("#status")).toHaveText(DONE, { timeout: 60_000 });
}

test("the Quick range state shows ion, energy, target and Run; the rest is under More options", async ({ page }) => {
  await expect(page.locator("#q-ion")).toHaveValue(golden.ion);
  // The default is 10 keV, as in the full simulator; the golden record is 5 keV.
  await expect(page.locator("#q-energy")).toHaveValue("10");
  await expect(page.locator("#q-target")).toHaveValue(golden.target);
  await expect(page.locator("#btn-run")).toBeVisible();
  await expect(page.locator("#more-options")).not.toHaveAttribute("open", "");
  await expect(page.locator("#tilt")).toBeHidden();
  await page.locator("#more-options > summary").click();
  await expect(page.locator("#ions")).toBeVisible();
  await expect(page.locator("#tilt")).toBeVisible();
  // Results: the tiles and the range plot only.
  await expect(page.locator("#counters")).toBeVisible();
  await expect(page.locator("#plots .plot:visible")).toHaveCount(1);
  await expect(page.locator("#traj")).toBeHidden();
  await expect(page.locator("#partition-bar")).toBeHidden();
});

test("the page opens with the labelled example run", async ({ page }) => {
  await expect(page.locator("#example-label")).toHaveText("Example run: 1000 ions, seed 1.");
  await expect(page.locator("#ions")).toHaveValue("1000");
  await expect(page.locator("#seed")).toHaveValue(String(golden.seed));
  await expect(page.locator("#c-done")).toHaveText("1000 / 1000");
  expect(await canvasColours(page.locator('#plots .plot[data-plot="range"] canvas'))).toBeGreaterThan(3);
});

test("B into Si matches the native golden record", async ({ page }) => {
  await page.locator("#q-energy").fill(String(golden.energy_ev / 1000));
  await page.locator("#more-options > summary").click();
  await page.locator("#ions").fill(String(golden.ions));
  await page.locator("#seed").fill(String(golden.seed));
  await run(page);
  await expect(page.locator("#c-done")).toHaveText(`${golden.ions} / ${golden.ions}`);
  // The tolerance of lindhard-wasm/tests/parity.rs, widened by the digits the tiles print.
  expect(Math.abs((await tileNumber(page, "c-rp")) - golden.rp_nm)).toBeLessThanOrEqual(golden.rp_std_err_nm + 0.005);
  expect(Math.abs((await tileNumber(page, "c-drp")) - golden.drp_nm)).toBeLessThanOrEqual(golden.drp_std_err_nm + 0.005);
  const sigma = Math.sqrt(golden.stopped * (1 - golden.stopped / golden.ions));
  expect(Math.abs(Number(await page.locator("#c-stopped").innerText()) - golden.stopped)).toBeLessThanOrEqual(sigma);
  expect(await canvasColours(page.locator('#plots .plot[data-plot="range"] canvas'))).toBeGreaterThan(3);
});

test("the same seed gives the same result twice", async ({ page }) => {
  const tiles = () => page.locator("#counters").innerText();
  const first = await tiles();
  await run(page);
  expect(await tiles()).toBe(first); // the example is 1000 ions, seed 1
  await page.locator("#more-options > summary").click();
  await page.locator("#seed").fill("2");
  await run(page);
  expect(await tiles()).not.toBe(first);
});

test("another ion and target run", async ({ page }) => {
  await page.locator("#q-ion").fill("As");
  await page.locator("#q-energy").fill("50");
  await page.locator("#q-target").fill("Ge");
  await run(page);
  expect(await tileNumber(page, "c-rp")).toBeGreaterThan(0);
});

test("an unknown element is reported, not thrown", async ({ page }) => {
  await page.locator("#q-ion").fill("Xx");
  await expect(page.locator("#problems-box")).toBeVisible();
  await expect(page.locator("#problems")).toContainText('unknown element "Xx"');
  await expect(page.locator("#q-ion")).toHaveAttribute("aria-invalid", "true");
  await expect(page.locator("#btn-run")).toBeDisabled();
  await page.locator("#q-ion").fill("B");
  await expect(page.locator("#q-ion")).not.toHaveAttribute("aria-invalid", "true");
  await expect(page.locator("#btn-run")).toBeEnabled();
});

test("the plot is readable in dark mode", async ({ page }) => {
  await page.emulateMedia({ colorScheme: "dark" });
  await run(page);
  // WCAG contrast of the lightest canvas pixel (axis text) against the page background.
  const ratio = await page.evaluate(() => {
    const lum = ([r, g, b]) => {
      const f = (c) => ((c /= 255) <= 0.03928 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
      return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b);
    };
    const bg = getComputedStyle(document.body).backgroundColor.match(/\d+/g).slice(0, 3).map(Number);
    const c = document.querySelector('#plots .plot[data-plot="range"] canvas');
    const d = c.getContext("2d").getImageData(0, 0, c.width, c.height).data;
    let best = 0;
    for (let i = 0; i < d.length; i += 4) if (d[i + 3] === 255) best = Math.max(best, lum([d[i], d[i + 1], d[i + 2]]));
    const l = lum(bg);
    return (Math.max(best, l) + 0.05) / (Math.min(best, l) + 0.05);
  });
  expect(ratio).toBeGreaterThanOrEqual(4.5);
});
