// The offline cache's service worker (web/sw.js), run in a stub of the
// service-worker globals.
//
//   node --test lindhard-wasm/tests/node/sw.test.mjs

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import vm from "node:vm";

const source = readFileSync(new URL("../../web/sw.js", import.meta.url), "utf8");
const ORIGIN = "https://example.test";
const SCOPE = `${ORIGIN}/calculator/`;

/** A response that followed a redirect, as fetch() reports one. */
function redirected(body) {
  const res = new Response(body, { status: 200, headers: { "content-type": "text/html" } });
  Object.defineProperty(res, "redirected", { value: true });
  return res;
}

/** Load sw.js with `fetch` answered by `answer(url)`; returns its listeners and its cache. */
function load(answer) {
  const store = new Map();
  const key = (r) => new URL(typeof r === "string" ? r : r.url, SCOPE).href;
  const cache = {
    put: async (r, res) => { store.set(key(r), res); },
    match: async (r) => store.get(key(r)),
    addAll: async () => { throw new Error("cache.addAll stores redirected responses"); },
  };
  const listeners = {};
  const context = vm.createContext({
    URL, Response, Promise, Error,
    Request: class { constructor(url) { this.url = new URL(url, SCOPE).href; this.method = "GET"; } },
    caches: { open: async () => cache, keys: async () => [], delete: async () => true },
    fetch: async (r) => answer(key(r)),
    self: {
      location: new URL("sw.js", SCOPE),
      addEventListener: (name, fn) => { listeners[name] = fn; },
      skipWaiting: async () => {},
      clients: { claim: async () => {} },
    },
  });
  vm.runInContext(source, context);
  return { listeners, store };
}

/** Dispatch an event and wait for what it was told to wait for. */
async function dispatch(listener, event = {}) {
  let done;
  listener({ ...event, waitUntil: (p) => { done = p; }, respondWith: (p) => { done = p; } });
  return done;
}

test("install stores plain responses where the host redirected", async () => {
  const { listeners, store } = load((url) => (url.endsWith(".html") ? redirected(`page ${url}`) : new Response(`file ${url}`)));
  await dispatch(listeners.install);
  const page = store.get(`${SCOPE}simulator.html`);
  assert.ok(page, "simulator.html is precached under its own URL");
  assert.equal(page.redirected, false);
  assert.equal(page.status, 200);
  assert.equal(page.headers.get("content-type"), "text/html");
  assert.equal(await page.text(), `page ${SCOPE}simulator.html`);
  for (const [url, res] of store) assert.equal(res.redirected, false, url);
  assert.ok(store.size >= 10, "every precache entry is stored");
});

test("install fails when a precached file is missing", async () => {
  const { listeners } = load((url) => (url.endsWith("worker.js") ? new Response("", { status: 404 }) : new Response("ok")));
  await assert.rejects(dispatch(listeners.install), /worker\.js: HTTP 404/);
});

test("a refreshed copy that followed a redirect is stored plain", async () => {
  const { listeners, store } = load(() => {
    const res = redirected("fresh");
    Object.defineProperty(res, "type", { value: "basic" });
    res.clone = () => { const c = redirected("fresh"); Object.defineProperty(c, "type", { value: "basic" }); return c; };
    return res;
  });
  const request = { method: "GET", url: `${SCOPE}plots.js` };
  await dispatch(listeners.fetch, { request });
  const stored = store.get(`${SCOPE}plots.js`);
  assert.equal(stored.redirected, false);
  assert.equal(await stored.text(), "fresh");
});
