# Simulator in the browser

The [simulator](../calculator/simulator.html) is the layered-target companion
of the [range calculator](range-calculator.md). It runs the same engine as
WebAssembly in a Web Worker, so the page stays responsive during a run, and it
is static: there is no server and nothing you enter leaves your machine.

This first release has the **Setup** view and a minimal **Run** view. The
trajectory view, distribution plots and exports are still to come.

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
