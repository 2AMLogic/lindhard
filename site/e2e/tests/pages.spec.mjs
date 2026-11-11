// Every published page loads cleanly and the landing page's links resolve.
import { readFileSync } from "node:fs";
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
  // Nav entry, hero button, hero figure and the "Try it" card all link the simulator.
  const sim = page.locator('a[href="calculator/simulator.html"]');
  expect(await sim.count()).toBeGreaterThanOrEqual(4);
  await expect(sim.first()).toBeVisible();
  await expect(sim.last()).toBeVisible();
  const r = await request.get("calculator/simulator.html");
  expect(r.status()).toBe(200);
  await expect(page.getByText(/planned/i)).toHaveCount(0);
});

test("at 1440 px the first screen shows the description, both buttons and the figure", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto("index.html");
  const inView = async (loc) => {
    await expect(loc).toBeVisible();
    const box = await loc.boundingBox();
    expect(box.y + box.height, "bottom edge within the first screen").toBeLessThanOrEqual(900);
  };
  await inView(page.locator(".tagline"));
  await inView(page.getByRole("link", { name: "Open the simulator" }));
  await inView(page.getByRole("link", { name: "Read the manual" }));
  await inView(page.locator(".hero-figure img"));
});

test("the hero figure loads, links to the simulator and has a generated caption", async ({ page, request }) => {
  await page.goto("index.html");
  const img = page.locator(".hero-figure img");
  await expect.poll(() => img.evaluate((i) => i.complete && i.naturalWidth > 0)).toBe(true);
  await expect(page.locator(".hero-figure a")).toHaveAttribute("href", "calculator/simulator.html");
  const caption = await page.locator(".hero-figure figcaption").innerText();
  for (const word of [/5 keV B into Si/, /\d+ ions/, /seed \d+/, /models:/]) expect(caption).toMatch(word);
  const generated = await request.get("gallery/generated/quickstart.md");
  expect(generated.status()).toBe(200);
});

test("the metadata line carries the workspace version", async ({ page }) => {
  const cargo = readFileSync(new URL("../../../Cargo.toml", import.meta.url), "utf8");
  const version = cargo.match(/\[workspace\.package\][^\[]*?version\s*=\s*"([^"]+)"/)[1];
  await page.goto("index.html");
  await expect(page.locator(".meta")).toContainText(version);
});

test("the trust section shows checked rows, known deviations and links to the full record", async ({ page }) => {
  await page.goto("index.html");
  const callout = page.locator("section[aria-labelledby=validation] .callout");
  expect(await callout.locator("table.trust tbody tr").count()).toBeGreaterThanOrEqual(3);
  expect(await callout.locator(".deviations li").count()).toBeGreaterThanOrEqual(2);
  await expect(callout.locator('a[href="validation.html"]')).toBeVisible();
  await expect(callout.locator('a[href$="data-provenance.md"]')).toBeVisible();
});

test("the examples show four thumbnails that link to the gallery pages", async ({ page }) => {
  await page.goto("index.html");
  const thumbs = page.locator(".cards.thumbs a.card");
  await expect(thumbs).toHaveCount(4);
  for (const img of await thumbs.locator("img").all()) {
    // The thumbnails are loading="lazy"; scroll each into view so it is requested.
    await img.scrollIntoViewIfNeeded();
    await expect.poll(() => img.evaluate((i) => i.complete && i.naturalWidth > 0)).toBe(true);
  }
  await expect(page.locator('a[href="gallery/index.html"]')).toBeVisible();
});

test("the build and run blocks scroll on a phone and copy to the clipboard", async ({ page, context, browserName }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto("index.html");
  const block = page.locator(".codeblock").first();
  const pre = block.locator("pre");
  expect(await pre.evaluate((p) => getComputedStyle(p).overflowX)).toBe("auto");
  await expect(page.locator(".codeblock")).toHaveCount(2);
  await expect(block.locator("button.copy")).toBeVisible();
  if (browserName === "chromium") {
    await context.grantPermissions(["clipboard-read", "clipboard-write"]);
    await block.locator("button.copy").click();
    expect(await page.evaluate(() => navigator.clipboard.readText())).toContain("cargo build --release");
  }
});

test("without a clipboard API the copy buttons stay hidden", async ({ page }) => {
  await page.addInitScript(() => Object.defineProperty(navigator, "clipboard", { value: undefined }));
  await page.goto("index.html");
  expect(await page.locator("button.copy:visible").count()).toBe(0);
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
