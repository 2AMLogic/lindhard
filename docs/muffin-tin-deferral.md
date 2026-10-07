# Muffin-tin potential for condensed targets: deferred

Issue #92 asked for either a muffin-tin option on the elastic solver
(`electron::elastic`, #17) or a sourced note arguing for deferral. This is the
note. **Decision: deferred.** The decision is recorded in
[`architecture.md`](architecture.md); the gap is a row in
[`data-provenance.md`](data-provenance.md).

## What is deferred

In a solid the atomic potentials of neighbours overlap. The muffin-tin
construction replaces each free-atom potential by one that is cut off at a
radius set by the number density, and (in the usual form) shifted so that it is
continuous there. The elastic solver
([`lindhard/src/electron/elastic.rs`](../lindhard/src/electron/elastic.rs))
currently solves the radial Dirac equation for the free-atom screened
potential only.

## Why: no readable definition

A muffin-tin option needs two conventions that a "sigma_el goes down at low
energy" test cannot check: the truncation radius and the offset applied to the
potential inside it. The project's source for the construction is Salvat,
Jablonski and Powell, Comput. Phys. Commun. 165, 157 (2005),
doi:10.1016/j.cpc.2004.09.006. It is closed access. Checked 2026-10-07:
Unpaywall reports `is_oa: false`, `best_oa_location: null`,
`has_repository_copy: false`; the publisher page answers HTTP 403. A Crossref
search and an arXiv search for muffin-tin treatments of elastic electron
scattering found no open text that defines the construction as used for
electron Monte Carlo cross sections (the Crossref hits were band-structure and
multiple-scattering papers, which are a different construction and were not
used). No source defining the offset and the truncation was opened, so none is
reconstructed from memory.

## Size of the effect versus energy: not quantified here

No published comparison of free-atom and muffin-tin elastic cross sections
could be opened, so this note gives **no number**. What can be said without a
source, from the structure of the problem alone:

- The truncation only removes the part of the potential beyond the
  truncation radius. A projectile whose wavelength is short compared with that
  radius and whose energy is large compared with the offset is little
  affected, so the effect should fade with energy; the size and the sign of
  the change at a given low energy are not asserted here.
- The issue's premise (free-atom potentials overestimate low-energy elastic
  cross sections) is unverified by us.

## What the missing option leaves

Elastic cross sections of lindhard are free-atom values. At low energy they
are expected to differ from solid-state values by an amount we cannot give.
Consumers needing low-energy accuracy should treat the free-atom
approximation as an uncontrolled error there until the option is built.

## Pieces that are plain geometry (ready for whoever builds it)

The Wigner-Seitz radius from the number density `n = 1 / V_atom` is
`R_ws = (3 / (4 pi n))^(1/3)`, the radius of a sphere holding one atom's
volume. `crate::elements::element(z).density_kg_m3()` and
`ion::dynamic::atomic_volume_from_density(z)` already supply `n`. Whether the
muffin-tin radius is `R_ws` or a touching-sphere radius is part of the
undefined convention above. A wrapper type implementing `ScreenedPotential`
around any inner potential would need no solver change: its `breakpoints()`
must include the truncation radius and its `matching_radius` must not exceed
it. Compounds and per-layer densities are out of scope for a first version.

## To lift the deferral

A reader with access to Salvat, Jablonski and Powell (2005) (or another source
that defines the construction) records the exact equations and page, and a
published free-atom versus muffin-tin comparison for estimating the size.
Then implement `MuffinTin<P: ScreenedPotential>`, test it with the Yukawa
stand-in (the DHFS coefficients, #17, are also a gap), and measure the
direction and size of the effect honestly. Follow-ups that depend on this:
#90 (table provenance must record whether a muffin-tin was applied), #91
(exchange and polarization: decide the order of application) and #93 (say
whether reference values include a muffin-tin).
