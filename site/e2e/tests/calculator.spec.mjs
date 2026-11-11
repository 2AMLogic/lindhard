// The range calculator: the engine loads in the browser and reproduces the
// native result recorded in lindhard-wasm/tests/b_5keV_si.golden.json.
import { readFileSync } from "node:fs";
import { test, expect, canvasColours, cellNumber } from "./fixtures.mjs";

const golden = JSON.parse(
  readFileSync(new URL("../../../lindhard-wasm/tests/b_5keV_si.golden.json", import.meta.url), "utf8"),
);

test.beforeEach(async ({ page }) => {
  await page.goto("calculator/index.html");
  await expect(page.locator("#status")).toHaveText("Ready.");
});

test("the defaults are the golden case except the energy", async ({ page }) => {
  await expect(page.locator("#ion")).toHaveValue(golden.ion);
  // The page defaults to 10 keV, as the simulator does; the golden record is 5 keV.
  await expect(page.locator("#energy")).toHaveValue("10");
  await expect(page.locator("#target")).toHaveValue(golden.target);
  await expect(page.locator("#ions-n")).toHaveValue(String(golden.ions));
  await expect(page.locator("#seed")).toHaveValue(String(golden.seed));
});

test("B into Si matches the native golden record", async ({ page }) => {
  await page.locator("#energy").fill(String(golden.energy_ev / 1000));
  await page.locator("#run").click();
  await expect(page.locator("#status")).toContainText("done in", { timeout: 60_000 });
  // The tolerance of lindhard-wasm/tests/parity.rs, widened by the two or
  // three digits the page prints.
  expect(Math.abs((await cellNumber(page, "rp")) - golden.rp_nm)).toBeLessThanOrEqual(golden.rp_std_err_nm + 0.005);
  expect(Math.abs((await cellNumber(page, "drp")) - golden.drp_nm)).toBeLessThanOrEqual(golden.drp_std_err_nm + 0.005);
  const sigma = Math.sqrt(golden.stopped * (1 - golden.stopped / golden.ions));
  expect(Math.abs((await cellNumber(page, "stopped")) - golden.stopped)).toBeLessThanOrEqual(sigma);
  await expect(page.locator("#stopped")).toContainText(`of ${golden.ions} ions`);
  await expect(page.locator("#plot")).toBeVisible();
  expect(await canvasColours(page.locator("#plot"))).toBeGreaterThan(3);
});

test("the same seed gives the same result twice", async ({ page }) => {
  await page.locator("#ions-n").fill("1000");
  const run = async () => {
    await page.locator("#run").click();
    await expect(page.locator("#status")).toContainText("done in", { timeout: 60_000 });
    return page.locator("#result").innerText();
  };
  const first = await run();
  expect(await run()).toBe(first);
  await page.locator("#seed").fill("2");
  expect(await run()).not.toBe(first);
});

test("another ion and target run", async ({ page }) => {
  await page.locator("#ion").fill("As");
  await page.locator("#energy").fill("50");
  await page.locator("#target").fill("Ge");
  await page.locator("#ions-n").fill("1000");
  await page.locator("#run").click();
  await expect(page.locator("#status")).toContainText("done in", { timeout: 60_000 });
  expect(await cellNumber(page, "rp")).toBeGreaterThan(0);
});

test("an unknown element is reported, not thrown", async ({ page }) => {
  await page.locator("#ion").fill("Xx");
  await page.locator("#run").click();
  await expect(page.locator("#status")).toHaveClass("error");
  await expect(page.locator("#status")).toHaveText("Unknown element symbol \u201cXx\u201d.");
  await expect(page.locator("#ion")).toHaveAttribute("aria-invalid", "true");
  await page.locator("#ion").fill("B");
  await expect(page.locator("#ion")).not.toHaveAttribute("aria-invalid", "true");
  await expect(page.locator("#run")).toBeEnabled();
});

test("the plot is readable in dark mode", async ({ page }) => {
  await page.emulateMedia({ colorScheme: "dark" });
  await page.locator("#ions-n").fill("1000");
  await page.locator("#run").click();
  await expect(page.locator("#status")).toContainText("done in", { timeout: 60_000 });
  // WCAG contrast of the lightest canvas pixel (axis text) against the page background.
  const ratio = await page.evaluate(() => {
    const lum = ([r, g, b]) => {
      const f = (c) => ((c /= 255) <= 0.03928 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
      return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b);
    };
    const bg = getComputedStyle(document.body).backgroundColor.match(/\d+/g).slice(0, 3).map(Number);
    const c = document.getElementById("plot");
    const d = c.getContext("2d").getImageData(0, 0, c.width, c.height).data;
    let best = 0;
    for (let i = 0; i < d.length; i += 4) if (d[i + 3] === 255) best = Math.max(best, lum([d[i], d[i + 1], d[i + 2]]));
    const l = lum(bg);
    return (Math.max(best, l) + 0.05) / (Math.min(best, l) + 0.05);
  });
  expect(ratio).toBeGreaterThanOrEqual(4.5);
});
