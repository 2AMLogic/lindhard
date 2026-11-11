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

test("the defaults are the golden case", async ({ page }) => {
  await expect(page.locator("#ion")).toHaveValue(golden.ion);
  await expect(page.locator("#energy")).toHaveValue(String(golden.energy_ev / 1000));
  await expect(page.locator("#target")).toHaveValue(golden.target);
  await expect(page.locator("#ions-n")).toHaveValue(String(golden.ions));
  await expect(page.locator("#seed")).toHaveValue(String(golden.seed));
});

test("B into Si matches the native golden record", async ({ page }) => {
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
  await expect(page.locator("#status")).toContainText(/unknown element/i);
  await expect(page.locator("#run")).toBeEnabled();
});
