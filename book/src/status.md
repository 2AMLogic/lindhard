# Status and roadmap

lindhard is early, pre-release software (version 0.0.1; no release is
published yet, so it is built from source). State as of 2026-11-01:

| Milestone | Scope | State |
|---|---|---|
| M0 | Amorphous ion core: materials, potentials, scattering integral, electronic stopping, layered BCA with cascades, damage and range tallies, CLI, validation harness | Most mature part; validated as the [validation page](validation.md) reports |
| M1 | Low-energy electron engine: Mott elastic scattering, dielectric inelastic scattering, secondary electrons, layered stacks, energy-deposition maps | Under active development; the [resist example](gallery/resist-on-si.md) runs it with synthetic electron data |
| M2 | Crystalline targets: lattice, thermal vibration, tilt and twist, dose-dependent damage | Under active development |
| M3 | Dynamic composition, 2D/3D geometry, Python wheels, WebAssembly | Dynamic composition and the browser [calculator](guide/range-calculator.md) and [simulator](guide/simulator.md) exist; the rest is in progress |

The plan of record is
[`WORK_PLAN.md`](https://github.com/2AMLogic/lindhard/blob/main/WORK_PLAN.md)
and the design is
[`docs/architecture.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/architecture.md).
