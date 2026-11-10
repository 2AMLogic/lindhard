// Browser tests of the public site.
//
//   SITE_URL  test a deployed site (for example https://lindhard.org). No
//             local server is started.
//   SITE_DIR  otherwise, serve this built tree (default ../dist, which
//             scripts/deploy-site.sh --no-deploy writes) on a local port.
import { defineConfig, devices } from "@playwright/test";
import { fileURLToPath } from "node:url";

const live = process.env.SITE_URL;
const dir = process.env.SITE_DIR ?? fileURLToPath(new URL("../dist", import.meta.url));
const port = Number(process.env.SITE_PORT ?? 4173);

export default defineConfig({
  testDir: "tests",
  fullyParallel: true,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 1 : 0,
  timeout: 90_000,
  expect: { timeout: 15_000 },
  reporter: process.env.CI ? [["list"], ["html", { open: "never" }]] : "list",
  use: {
    // A trailing slash, so the tests' relative paths resolve under a base
    // path as well (the GitHub Pages copy is served from /lindhard/).
    baseURL: live ? live.replace(/\/*$/, "/") : `http://127.0.0.1:${port}/`,
    trace: "retain-on-failure",
    screenshot: "only-on-failure",
  },
  webServer: live
    ? undefined
    : {
        command: `python3 -m http.server ${port} --bind 127.0.0.1 --directory "${dir}"`,
        url: `http://127.0.0.1:${port}/index.html`,
        reuseExistingServer: !process.env.CI,
      },
  projects: [
    { name: "chromium", use: { ...devices["Desktop Chrome"] } },
    { name: "firefox", use: { ...devices["Desktop Firefox"] } },
    { name: "webkit", use: { ...devices["Desktop Safari"] } },
  ],
});
