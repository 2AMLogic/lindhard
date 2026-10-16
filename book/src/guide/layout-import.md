# Layout cross-section import (finding)

**Status: blocked upstream.** No import code exists. This page records what
the pinned layout tool offers, so the next step is a precise upstream request
rather than a second layout reader in this tree.

## The pin

Layout reading is not reimplemented in lindhard; it is meant to come from the
[`klayout-tools`](https://pypi.org/project/klayout-tools/) package at one
pinned version. The pin lives in one file, `tools/klayout-tools.pin`
(currently `0.7.0`, the latest PyPI release when this was written). Bump it
deliberately, in its own change, and re-run the check below.

Run a pinned release without touching any installed tool:

```sh
uvx --from "klayout-tools==$(cat tools/klayout-tools.pin)" klt --help
```

## What 0.7.0 offers

Checked by running each verb on a two-rectangle layout written with `klt draw`.

| Verb | Gives | Polygons? |
|---|---|---|
| `klt layers --flattened --format json` | layer/datatype list, shape counts, per-layer bounding box | no |
| `klt stats --per-layer --format json` | per-layer drawn area (`area_um2`, a raw sum with overlaps double-counted; `merged_area_um2` is the union), polygon and vertex counts | no |
| `klt clip --region '[l,b,r,t]' -o out.gds` | a GDSII/OASIS stream of the region | only as a stream |
| `klt draw` | writes a stream from JSON polygons | write side only |

The Python modules of the package (`klayout_tools._layout` and friends) are
underscore-private and not documented as an interface, and `klayout.db` is a
transitive dependency, not part of the package's contract. None of them is a
supported polygon export.

## Finding

The pinned release has **no machine-readable path from a layout region to
polygons**. The only way to get geometry out of a region is `klt clip`, which
writes a stream, and reading that stream back would be a second GDS/OASIS
reader in this tree. That is ruled out, so the cross-section import (process
stack, 2D/3D target build, area cross-check against the layout tool) is not
started.

## What is needed upstream

A read-side verb or documented Python function in `klayout-tools` that, for a
layout, an optional `--region` and a layer list, returns per-layer merged
polygons as JSON:

- units of micrometres and the database unit, as `klt layers` already reports;
- one entry per `(layer, datatype)` with a list of polygons, each a hull plus
  holes as vertex lists, hierarchy flattened and merged (same rule as
  `merged_area_um2`);
- a `schema_version`, so the seam can fail on an unsupported version.

With that, the import needs only a process-stack input format (layer to
material, thickness, mask/etch rule), its validation, the target builders on
the Phase 1 geometry kernels, and a test comparing imported areas with
`klt stats` for a fixture layout.
