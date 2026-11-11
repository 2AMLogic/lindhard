// Shared fixtures and helpers.
//
// Every test fails if the page logs a console error, throws, or makes a
// request that fails or returns 4xx/5xx, so "the page loaded" also means
// "nothing was missing". Requests to other origins are not counted: the pages
// link to the repository but load nothing from it.
import { test as base, expect } from "@playwright/test";

// mdBook's book.js sets `window.onunload` (to keep the page out of the
// back/forward cache); Chromium reports that as a deprecation, which is not a
// defect of the site.
const KNOWN_CONSOLE = [/^Permissions policy violation: unload is not allowed/];

export const test = base.extend({
  problems: [
    async ({ page, baseURL }, use) => {
      const problems = [];
      const ours = (url) => url.startsWith(new URL(baseURL).origin);
      page.on("console", (m) => {
        if (m.type() === "error" && !KNOWN_CONSOLE.some((k) => k.test(m.text()))) problems.push(`console error: ${m.text()}`);
      });
      page.on("pageerror", (e) => problems.push(`uncaught: ${e.message}`));
      page.on("requestfailed", (r) => {
        if (ours(r.url())) problems.push(`request failed: ${r.url()} (${r.failure()?.errorText})`);
      });
      page.on("response", (r) => {
        if (ours(r.url()) && r.status() >= 400) problems.push(`HTTP ${r.status()}: ${r.url()}`);
      });
      await use(problems);
      expect(problems, "the page reported problems").toEqual([]);
    },
    { auto: true },
  ],
});
export { expect };

/** Number of distinct colours on a canvas: 1 means blank. */
export async function canvasColours(locator) {
  return locator.evaluate((c) => {
    const d = c.getContext("2d").getImageData(0, 0, c.width, c.height).data;
    const seen = new Set();
    for (let i = 0; i < d.length; i += 4) seen.add((d[i] << 16) | (d[i + 1] << 8) | d[i + 2]);
    return seen.size;
  });
}

/** Leading number of a cell such as "24.39 nm (± 0.131)" or "0.412". */
export async function cellNumber(page, id) {
  const text = await page.locator(`#${id}`).innerText();
  const m = text.match(/-?\d+(\.\d+)?(e[+-]?\d+)?/i);
  expect(m, `#${id} holds a number, got "${text}"`).not.toBeNull();
  return Number(m[0]);
}
