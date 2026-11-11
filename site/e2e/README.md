# Browser tests of the site

[Playwright](https://playwright.dev) tests that open the built site in
Chromium, Firefox and WebKit and use it as a reader would.

| File | What it checks |
|---|---|
| `tests/pages.spec.mjs` | Each published page loads with no console error and no missing file; the landing page's links resolve; it follows the colour scheme and fits a phone; the same site bar and footer on the landing page, a book page and both tools, the wordmark leads home, `aria-current` marks the section, and all four surfaces are dark with matching backgrounds under a dark system setting. |
| `tests/calculator.spec.mjs` | The range calculator loads the engine and reproduces the native result in `lindhard-wasm/tests/b_5keV_si.golden.json`, within the tolerance of `lindhard-wasm/tests/parity.rs`. |
| `tests/simulator.spec.mjs` | The simulator: the cases of issue #469 (B, P and As into Si; As into SiO₂ on Si; Ar into Cu with cascades; a transmitting film; 60° incidence), pause, resume (one toggle) and stop, the 3D trajectory view (drag, keys, reset), the desktop two-column layout (Setup's panels stacked in its column) and the phone tabs, number rounding, loading an input file, every download, validation messages, keyboard use, both themes and a phone screen. |
| `tests/book-layout.spec.mjs` | The book's tables and equations fit (`book/theme/lindhard.css`): at 1440 px no table is wider than its column except the listed content cases; at 390 px no page's content pane scrolls sideways; at 768 px the stopping-residuals tables do not scroll. |

Every test also fails if the page logs a console error, throws, or requests
something that fails (`tests/fixtures.mjs`).

## Running

```sh
scripts/deploy-site.sh --no-deploy        # builds site/dist
cd site/e2e
npm ci
npx playwright install chromium firefox webkit
npx playwright test                       # serves ../dist on a local port
npx playwright test --project=chromium simulator
```

To test a deployed copy, name it: `SITE_URL=https://lindhard.org npx playwright test`.
`SITE_DIR` points the local server at another built tree.

CI runs the suite against the built tree on every pull request and push
(`.github/workflows/book.yml`, job `e2e`); neither host is published unless it
passes. After a publish, the same suite runs against the public site
(job `e2e-live`).

## Known gap

`a second run works after the network is gone` is marked `fixme`: the page
does not yet work offline (#653). The test is there so that fixing it is a one-word
change.

The physics is not tested here. Native tests and the Node tests in
`lindhard-wasm/tests` do that; these tests ask whether the pages work.
