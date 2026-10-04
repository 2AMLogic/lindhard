# lindhard

**A clean-room, MIT-licensed Rust engine for Monte Carlo transport of ions and
electrons in matter.** It computes implantation range profiles, damage,
sputtering and backscatter for ions, and energy deposition for low-energy
electrons. SRIM/TRIM covers the same problem class; lindhard uses none of its
code or data.

> **Status: scaffold.** Nothing is implemented yet. See [`WORK_PLAN.md`](WORK_PLAN.md)
> for milestones and [`docs/architecture.md`](docs/architecture.md) for the design.

## Why another code

Some open codes already handle parts of this problem well. The gaps are in
license and scope:

- **SRIM/TRIM** is closed source, Windows-only and non-commercial. It treats
  every target as amorphous and as fixed, so the target never changes with
  dose.
- **Crystalline implant simulation** is not available in any open, maintained
  code. That covers channeling under tilt and twist, screen oxides, and damage
  that builds up to amorphization. Today it means IMSIL, Crystal-TRIM, MARLOWE
  or commercial TCAD.
- **Ion codes that are open** are either copyleft (RustBCA, iradina, IM3D) or
  use stopping tables taken from SRIM (OpenTRIM).
- **Low-energy electron transport** (eV–50 keV, for e-beam lithography
  proximity functions, SEM and EBIC) has one permissive open code, Nebula, and
  it has stalled. The reference physics lives inside Geant4 or AGPL codes.

lindhard aims to be **permissive, embeddable, fast, reproducible and
validated**:

- MIT license, with no copyleft dependency anywhere in the build.
- A library first; the CLI and later Python/WASM bindings sit on top of it.
- Every ion runs on its own seeded random stream, so results do not depend on
  the thread count.
- Every model is checked against analytic limits, independent codes and
  published measurements ([`docs/validation.md`](docs/validation.md)).

[`docs/prior-art.md`](docs/prior-art.md) is the survey this project starts
from.

The name honours Jens Lindhard. His LSS theory is the standard theory of ion
range, he did the foundational work on channeling, and his dielectric function
is the starting point for inelastic electron scattering. Between them those
cover both halves of this project.

## Layout

| Path | What |
|---|---|
| `lindhard/` | The library: materials, potentials, stopping, transport, tallies |
| `lindhard-cli/` | `lindhard` binary: TOML in, JSON/CSV out |
| `docs/` | Architecture, prior art, validation plan, data provenance |

## Contributing

Read [`CONTRIBUTING.md`](CONTRIBUTING.md) first. It sets out the **clean-room
rule**: what may be learned from other codes, what may be ported, and what may
never enter this tree.

## License

MIT. See [`LICENSE`](LICENSE).
