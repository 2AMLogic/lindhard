# Simulator in the browser

The [simulator](../calculator/simulator.html) is the layered-target companion
of the [range calculator](range-calculator.md). It runs the same engine as
WebAssembly in a Web Worker, so the page stays responsive during a run, and it
is static: there is no server and nothing you enter leaves your machine.

The page has a **Setup** view and a **Run** view. The Run view shows live
counters, a trajectory view, the energy partition, live distribution plots and
downloads.

## Run

- **Counters**: ions done, R<sub>p</sub> and ΔR<sub>p</sub>, lateral (y, z) and
  radial spread, ions stopped, backscattered and transmitted per ion, sputtered
  atoms per ion, vacancies per ion (the followed cascade in full-cascade mode,
  the NRT estimate in quick mode), and the share of the incident energy that went
  to electrons, to the lattice, out of the target, and elsewhere. Brackets give one
  standard error.
- **Trajectories**: the paths of the first 50 ion histories, in a side view
  (depth across, lateral y down, layers shaded and named) or an end-on y × z view.
  Primary paths are always drawn; recoil cascades can be switched off. Species are
  told apart by colour and by the shape of the marker at the track's end, and the
  generation by the dash pattern (solid primary, dashed first recoils, dotted later
  ones). A history cut at the per-track point limit is flagged.
- **Distributions**, updated as the run proceeds: the ion range distribution
  (total and, with more than one layer, per layer, with the Pearson IV density of
  the four depth moments drawn over it when it exists; Hofker 1975, Heinrich 2004),
  the lateral and radial distributions, damage versus depth (vacancies,
  interstitials, replacements; zero in quick mode), and escapes (counts per species
  and face, with energy and polar-angle spectra). Series are told apart by colour
  and by dash pattern, with a legend; counts that fall outside the plotted range
  are stated under each plot.
- **Downloads**: each plot's data as CSV (one row per bin per series, with the
  under- and overflow counts as rows), the full report as JSON, and the trajectory
  view as a PNG. Results are shown with their statistical error, the models in use,
  and a pointer to `docs/validation.md`.
- **Limit**: the page runs up to 100 000 ions so that it stays responsive; use
  `lindhard run` for larger runs.

## Setup

- **Ion**: element, mass (the standard atomic weight unless set), energy in
  eV, keV or MeV, and the tilt from the surface normal (7° by default).
- **Target layers**, surface first: a name, a thickness (Å, nm or µm) or
  "semi-infinite" for the last layer, a composition by atom fraction, a density
  and the displacement, lattice-binding and surface-binding energies of each
  element. A blank energy uses the element table's value, shown in the field;
  a pure element takes its tabulated density; a compound needs a density.
  Layers can be added, deleted and reordered, and a preset fills one layer
  (Si, Cu, Au, SiO₂, GaAs). The preset densities are the ones recorded in
  [`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md);
  the O energies of the SiO₂ preset are the illustrative model parameters of
  `examples/as_50keV_si_sio2.toml`, not recommendations.
- **Calculation**: "ion distribution, quick damage" (`follow_recoils = false`)
  or "full cascades" (`follow_recoils = true`), the ion count and the seed.
- **Physics**: the potential and the electronic-stopping model, each with its
  citation, and the two energy cutoffs.

## The input file

The page shows the [TOML input](input.md) its setup produces, with Download and
Load buttons. The file runs unchanged under `lindhard run`, and gives the same
report for the same seed and ion count (checked natively in
`lindhard-wasm/tests/cli_equivalence.rs`). Loading a file fills the form. The
form exposes a subset of the input keys; a key it cannot represent (a crystal,
a dynamic target, a beam azimuth, `mass_fraction`, and so on) is reported by
name and the file is not loaded, rather than silently dropped.

## Limits

Amorphous, static targets only. Results carry statistical error (one standard
error is shown) and the models listed on the page; the validation status is in
[`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md).
