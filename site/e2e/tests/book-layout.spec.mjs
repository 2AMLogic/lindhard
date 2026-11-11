// The book's tables and display equations fit their page (book/theme/lindhard.css).
//
// At desktop width a table may be wider than its column only if it is one of
// the known content cases below; on a phone no page's content pane may scroll
// sideways (a table or an equation scrolls inside its own box instead).
import { test, expect } from "./fixtures.mjs";

// Tables that are wider than their wrapper at 1440 px because of their
// content, not their styling: "page|first header cell". None remain (#663).
const KNOWN_WIDE = new Set();

// Every chapter, read from the book's own table of contents.
let pages;
test.beforeAll(async ({ playwright, baseURL }) => {
  const ctx = await playwright.request.newContext({ baseURL });
  const toc = await (await ctx.get("toc.html")).text();
  await ctx.dispose();
  pages = [...new Set([...toc.matchAll(/href="([^"#:]+\.html)"/g)].map((m) => m[1]))].filter((p) => p !== "print.html" && p !== "toc.html");
});

/** Wait until MathJax (when it loaded) has typeset the page. */
async function settle(page) {
  await page.evaluate(() => new Promise((done) => (window.MathJax?.Hub ? window.MathJax.Hub.Queue(done) : done())));
}

test("the book has pages to scan", () => {
  expect(pages.length).toBeGreaterThan(20);
});

test("no book table is wider than its column at 1440 px, except the known cases", async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  const wide = [];
  for (const path of pages) {
    await page.goto(path);
    await settle(page);
    const tables = await page.locator(".content .table-wrapper").evaluateAll((ws) =>
      ws.map((w) => ({
        head: w.querySelector("th")?.textContent.trim() ?? "",
        over: w.scrollWidth - w.clientWidth,
      })),
    );
    for (const t of tables) if (t.over > 1 && !KNOWN_WIDE.has(`${path}|${t.head}`)) wide.push(`${path}|${t.head} (+${t.over}px)`);
  }
  expect(wide, "tables that scroll sideways at 1440 px").toEqual([]);
});

test("no book page scrolls sideways at 390 px", async ({ page }) => {
  test.setTimeout(300_000);
  await page.setViewportSize({ width: 390, height: 844 });
  const over = [];
  for (const path of pages) {
    await page.goto(path);
    await settle(page);
    const px = await page.evaluate(() => {
      const pane = document.querySelector(".content") ?? document.body;
      return Math.max(pane.scrollWidth - pane.clientWidth, document.documentElement.scrollWidth - document.documentElement.clientWidth);
    });
    if (px > 0) over.push(`${path} (+${px}px)`);
  }
  expect(over, "pages whose content pane scrolls sideways").toEqual([]);
});
