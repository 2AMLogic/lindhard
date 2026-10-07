# Introduction

lindhard is a clean-room, MIT-licensed Rust engine for Monte Carlo transport
of ions in matter, in the binary-collision approximation (BCA). For a beam of
ions entering a layered amorphous target it computes where the ions stop,
how much damage their cascades make, and what is sputtered and backscattered.
An electron engine is planned and not yet started.

This book has two parts.

- **Physics manual.** One page per model: the equations, the assumptions,
  the range over which the model is meant to be used, and the papers it comes
  from. Each page names the input keys and the library types that select the
  model, so you can go from a line in an input file to its physics and back.
- **User guide.** Installing, a first run from start to finish (5 keV boron
  into silicon), the reference for the TOML input, and the output formats.

The source of truth for every equation is the code and its doc comments
(`cargo doc -p lindhard --open`); this manual is written from them and
checked against them in CI where that can be automated (every model variant
must be named in the manual). The design notes, the validation results and
the data provenance table stay in the repository's
[`docs/`](https://github.com/2AMLogic/lindhard/tree/main/docs) directory:

- [`docs/validation.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/validation.md):
  analytic checks, code-to-code comparisons and comparisons with measurement,
  with the known deviations.
- [`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md):
  where every number in the tree comes from, and how far it has been
  verified.
- [`docs/architecture.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/architecture.md):
  the crate layout and the milestone plan.

lindhard is early, pre-release software. The project's rules (the clean-room
license tiers in particular) are in
[`CONTRIBUTING.md`](https://github.com/2AMLogic/lindhard/blob/main/CONTRIBUTING.md).
