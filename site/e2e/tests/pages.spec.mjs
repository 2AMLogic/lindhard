// Every published page loads cleanly and the landing page's links resolve.
import { test, expect } from "./fixtures.mjs";

const PAGES = [
  ["index.html", /lindhard/],
  ["introduction.html", /./],
  ["guide/install.html", /./],
  ["guide/first-run.html", /./],
  ["models/bca.html", /./],
  ["models/stopping.html", /./],
  ["validation.html", /./],
  ["status.html", /./],
  ["about.html", /./],
  ["calculator/index.html", /range calculator/i],
  ["calculator/simulator.html", /simulator/i],
  ["api/lindhard/index.html", /lindhard/],
];

for (const [path, title] of PAGES) {
  test(`${path} loads`, async ({ page }) => {
    const response = await page.goto(path);
    expect(response.status()).toBe(200);
    await expect(page).toHaveTitle(title);
    await expect(page.locator("h1").first()).toBeVisible();
  });
}

test("the landing page says what the project is and links to the tools", async ({ page, request }) => {
  await page.goto("index.html");
  await expect(page.locator("h1")).toHaveText("lindhard");
  await expect(page.locator('a[href="calculator/index.html"]').first()).toBeVisible();
  await expect(page.locator('a[href="validation.html"]').first()).toBeVisible();
  // Nav entry and the "Try it" entry both link the simulator, and it resolves.
  const sim = page.locator('a[href="calculator/simulator.html"]');
  await expect(sim).toHaveCount(2);
  await expect(sim.first()).toBeVisible();
  await expect(sim.last()).toBeVisible();
  const r = await request.get("calculator/simulator.html");
  expect(r.status()).toBe(200);
  await expect(page.getByText(/planned/i)).toHaveCount(0);
});

test("every same-site link on the landing page resolves", async ({ page, request }) => {
  await page.goto("index.html");
  const here = new URL(page.url());
  const hrefs = await page.locator("a[href]").evaluateAll((as) => as.map((a) => a.href));
  const local = [...new Set(hrefs.map((h) => h.split("#")[0]))].filter((h) => new URL(h).origin === here.origin);
  expect(local.length).toBeGreaterThan(5);
  for (const url of local) {
    const r = await request.get(url);
    expect(r.status(), url).toBe(200);
  }
});

test("the landing page follows the reader's colour scheme", async ({ page }) => {
  const background = () => page.evaluate(() => getComputedStyle(document.body).backgroundColor);
  await page.emulateMedia({ colorScheme: "light" });
  await page.goto("index.html");
  const light = await background();
  await page.emulateMedia({ colorScheme: "dark" });
  expect(await background()).not.toBe(light);
});

test("the landing page fits a phone screen", async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto("index.html");
  const overflow = await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
  expect(overflow).toBeLessThanOrEqual(0);
});

test("the WebAssembly module is served with its media type", async ({ request }) => {
  const r = await request.get("calculator/pkg/lindhard_wasm_bg.wasm");
  expect(r.status()).toBe(200);
  expect(r.headers()["content-type"]).toContain("application/wasm");
});
