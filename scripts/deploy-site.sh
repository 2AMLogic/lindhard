#!/usr/bin/env bash
# Build the combined public site (landing page, book, range calculator,
# simulator, gallery, validation, API docs) into site/dist, and optionally
# publish it to Cloudflare Pages.
#
#   scripts/deploy-site.sh --no-deploy     build only (what CI runs on a PR)
#   scripts/deploy-site.sh                 build, then publish
#   scripts/deploy-site.sh --publish-only  publish an existing site/dist
#
# The build is the one .github/workflows/book.yml uses for GitHub Pages, so
# the two hosts serve the same tree and published URLs stay valid.
#
# Needs: a stable Rust toolchain with the wasm32-unknown-unknown target,
# wasm-bindgen at the version pinned in lindhard-wasm/Cargo.toml, mdbook, and
# python3. Publishing also needs node (npx) for the Cloudflare Pages CLI.
#
# Publishing credentials come from the environment (CLOUDFLARE_API_TOKEN and
# CLOUDFLARE_ACCOUNT_ID). If they are unset, the script sources an owner-only
# env file named by LINDHARD_CF_ENV_FILE (default ~/.cloudflare/lindhard.env),
# kept outside the repository. The Pages project name is LINDHARD_PAGES_PROJECT
# (default "lindhard"). Missing credentials are an error for a publish.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
dist="$root/site/dist"
mode=all

for arg in "$@"; do
  case "$arg" in
    --no-deploy) mode=build ;;
    --publish-only) mode=publish ;;
    -h|--help) sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

build() {
  cd "$root"
  local web=lindhard-wasm/web calc=book/src/calculator
  "$web/../build.sh" web "$calc/pkg"
  cp "$web"/index.html "$web"/simulator.html "$web"/simulator.js "$web"/worker.js \
     "$web"/toml_model.js "$web"/run_view.js "$web"/plots.js "$calc/"
  # Service worker (offline cache): stamp a hash of the shipped files into the
  # cache name so a changed deploy replaces the old cache.
  local build_id
  rm -f "$calc/sw.js"
  build_id="$(cat "$calc"/*.js "$calc"/*.html "$calc"/pkg/lindhard_wasm.js "$calc"/pkg/lindhard_wasm_bg.wasm | sha256sum | cut -c1-16)"
  sed "s/__BUILD_ID__/$build_id/" "$web/sw.js" > "$calc/sw.js"
  mdbook build book
  cargo doc -p lindhard --no-deps
  mkdir -p book/book/api
  cp -r target/doc/. book/book/api/
  python3 validation/assemble_site.py book/book
  python3 validation/check_site_links.py book/book
  rm -rf "$dist"
  mkdir -p "$dist"
  cp -r book/book/. "$dist/"
  echo "site built in $dist"
}

publish() {
  [ -f "$dist/index.html" ] || { echo "error: $dist has no index.html; build first" >&2; exit 1; }
  if [ -z "${CLOUDFLARE_API_TOKEN:-}" ]; then
    local env_file="${LINDHARD_CF_ENV_FILE:-$HOME/.cloudflare/lindhard.env}"
    [ -f "$env_file" ] || { echo "error: no CLOUDFLARE_API_TOKEN and no env file $env_file" >&2; exit 1; }
    set -a; . "$env_file"; set +a
  fi
  : "${CLOUDFLARE_API_TOKEN:?not set}" "${CLOUDFLARE_ACCOUNT_ID:?not set}"
  cd "$root"
  npx --yes wrangler@4 pages deploy "$dist" \
    --project-name "${LINDHARD_PAGES_PROJECT:-lindhard}" --branch main
}

case "$mode" in
  build) build ;;
  publish) publish ;;
  all) build; publish ;;
esac
