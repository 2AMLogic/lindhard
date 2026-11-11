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

test("the landing page says what the project is and links to the tools", async ({ page }) => {
  await page.goto("index.html");
  await expect(page.locator("h1")).toHaveText("lindhard");
  await expect(page.locator('a[href="calculator/index.html"]').first()).toBeVisible();
  await expect(page.locator('a[href="validation.html"]').first()).toBeVisible();
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

// ---- the shared site bar and footer ----------------------------------------

// One page per surface.
const SURFACES = [
  ["landing page", "index.html"],
  ["book", "guide/install.html"],
  ["range calculator", "calculator/index.html"],
  ["simulator", "calculator/simulator.html"],
];
const BAR_LINKS = ["Manual", "Simulator", "Range calculator", "Validation", "Status", "GitHub"];

for (const [name, path] of SURFACES) {
  test(`the ${name} has the site bar and footer`, async ({ page }) => {
    await page.goto(path);
    const bar = page.locator("header.lh-bar");
    await expect(bar).toBeVisible();
    await expect(bar.getByRole("navigation", { name: "Site" })).toBeVisible();
    await expect(bar.locator("a")).toHaveText(["lindhard", ...BAR_LINKS]);
    await expect(page.locator("footer.lh-footer")).toHaveCount(1);
    await expect(page.locator("footer.lh-footer a").first()).toBeVisible();
  });

  test(`the wordmark leads home from the ${name}`, async ({ page }) => {
    await page.goto(path);
    await page.locator("header.lh-bar .lh-wordmark").click();
    await expect(page).toHaveURL(/\/index\.html$|\/$/);
    await expect(page.locator("h1")).toHaveText("lindhard");
  });

  test(`the site bar of the ${name} fits a phone screen`, async ({ page }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await page.goto(path);
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
    expect(overflow).toBeLessThanOrEqual(0);
  });
}

test("the site bar marks the current section", async ({ page }) => {
  const current = async () => page.locator('header.lh-bar a[aria-current="page"]').allInnerTexts();
  for (const [path, label] of [
    ["guide/install.html", "Manual"],
    ["validation.html", "Validation"],
    ["status.html", "Status"],
    ["calculator/index.html", "Range calculator"],
    ["calculator/simulator.html", "Simulator"],
    ["index.html", "lindhard"],
  ]) {
    await page.goto(path);
    expect(await current(), path).toEqual([label]);
  }
});

test("with a dark system setting every surface is dark on first load, with the same background", async ({ browser }) => {
  const backgrounds = [];
  for (const [, path] of SURFACES) {
    const context = await browser.newContext({ colorScheme: "dark", baseURL: test.info().project.use.baseURL });
    const page = await context.newPage();
    await page.goto(path);
    backgrounds.push(
      await page.evaluate(() => {
        for (const el of [document.body, document.documentElement]) {
          const c = getComputedStyle(el).backgroundColor;
          if (c !== "rgba(0, 0, 0, 0)" && c !== "transparent") return c;
        }
        return "none";
      }),
    );
    await context.close();
  }
  expect(new Set(backgrounds).size, backgrounds.join(" | ")).toBe(1);
  const [r, g, b] = backgrounds[0].match(/\d+/g).map(Number);
  expect((r + g + b) / 3, backgrounds[0]).toBeLessThan(80);
});
