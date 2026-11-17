// Offline cache for the calculator pages (issue #653, request #469).
//
// Precaches the calculator's static files on install and answers same-origin
// GETs in this directory's scope from the cache, refreshing the cached copy
// from the network when it can. The Worker's own module fetches (worker.js,
// the wasm-bindgen glue, the .wasm) pass through here, which is why this is a
// service worker and not a Cache API call from the page.
//
// BUILD_ID is replaced by scripts/deploy-site.sh with a hash of the files, so
// every changed deploy gets a new cache name and `activate` deletes the old
// ones. Served unreplaced (local dev), it is a fixed name: reload twice.
const BUILD_ID = "__BUILD_ID__";
const CACHE = `lindhard-calculator-${BUILD_ID}`;
const PRECACHE = [
  "./",
  "index.html",
  "simulator.html",
  "simulator.js",
  "worker.js",
  "toml_model.js",
  "run_view.js",
  "plots.js",
  "tokens.css",
  "chrome.css",
  "pkg/lindhard_wasm.js",
  "pkg/lindhard_wasm_bg.wasm",
];

// A host may answer "simulator.html" with a redirect to "simulator" (Cloudflare
// Pages does). A response that followed a redirect must not answer a
// navigation: the browser fails the load ("Response served by service worker
// has redirections"). So a redirected response is rebuilt as a plain one, with
// the same body, status and headers, before it is stored.
function unredirected(res) {
  if (!res.redirected) return res;
  return new Response(res.body, { status: res.status, statusText: res.statusText, headers: res.headers });
}

self.addEventListener("install", (event) => {
  event.waitUntil(
    caches.open(CACHE)
      // Not cache.addAll: it stores the responses as fetched, redirects included.
      .then((cache) => Promise.all(PRECACHE.map(async (u) => {
        const res = await fetch(new Request(u, { cache: "reload" }));
        if (!res.ok) throw new Error(`precache ${u}: HTTP ${res.status}`);
        await cache.put(u, unredirected(res));
      })))
      .then(() => self.skipWaiting()),
  );
});

self.addEventListener("activate", (event) => {
  event.waitUntil(
    caches.keys()
      .then((keys) => Promise.all(
        keys.filter((k) => k.startsWith("lindhard-calculator-") && k !== CACHE)
          .map((k) => caches.delete(k)),
      ))
      .then(() => self.clients.claim()),
  );
});

self.addEventListener("fetch", (event) => {
  const req = event.request;
  if (req.method !== "GET") return;
  const url = new URL(req.url);
  if (url.origin !== self.location.origin) return;
  event.respondWith((async () => {
    const cache = await caches.open(CACHE);
    const cached = await cache.match(req, { ignoreSearch: true });
    const refresh = fetch(req).then((res) => {
      if (res.ok && res.type === "basic") cache.put(req, unredirected(res.clone()));
      return res;
    });
    if (cached) {
      refresh.catch(() => {}); // offline: the cached copy is the answer
      return cached;
    }
    return refresh;
  })());
});
