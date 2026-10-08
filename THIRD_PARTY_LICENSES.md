# Third-party licenses

Code ported from Tier A (permissive) sources, as described in
`CONTRIBUTING.md`, is listed here with its upstream notice. Dependency
licenses are gated separately by `cargo deny` (`deny.toml`).

## Nebula and cstool (Nebula-simulator, TU Delft)

- Nebula, <https://github.com/Nebula-simulator/nebula>, commit
  `a50a8e83207980ed221ac83cfe874751d39b86a9`:
  - `source/physics/boundary_intersect.h` (the potential step: transmission,
    refraction, reflection), ported to `lindhard/src/electron/boundary.rs`
    (`cross_step`);
  - `source/physics/kieft/inelastic.h` (secondary-electron energy, direction
    and primary deflection, and the clamp of the loss to `E - E_F`), ported
    to `lindhard/src/electron/secondary.rs` (`kieft_bosch`) and
    `lindhard/src/electron/transport.rs`;
  - `source/drivers/cpu/cpu_driver.inl` (the stopping threshold measured from
    the vacuum level), followed in `lindhard/src/electron/transport.rs`
    (`CutoffReference::VacuumLevel`).
- cstool, <https://github.com/Nebula-simulator/cstool>, commit
  `0c739eb3fcc3fe5297e74c601ac4a9546db596cf`:
  `cstool/input_data/band_structure.py` (the metal and insulator band models:
  Fermi energy and inner potential), ported to
  `lindhard/src/electron/boundary.rs` (`BandModel`, `BandStructure`). No
  cstool data file is copied.

Both repositories carry the same licence:

```text
BSD 3-Clause License

Copyright (c) 2020, Nebula-simulator
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## DISPLATH (permissionx)

- DISPLATH, <https://github.com/permissionx/DISPLATH>, commit
  `7f461141c518305e6a5de8ce7e37dd314cd20ccb`:
  - `src/dynamics.jl` (`GetTargetsFromNeighbor`: nearest target by path
    distance to the closest-approach point, then the simultaneous partners;
    `Collision_!`: the momentum balance of simultaneous partners and the
    common energy-scaling factor) and `src/geometry.jl`
    (`SimultaneousCriteria`), re-implemented in Rust (not translated line by
    line) in `lindhard/src/ion/bca/crystal.rs` (`Bca::crystal_step`,
    `Bca::crystal_collide`) on this crate's own lattice search, scattering
    tables, kinematics and energy bookkeeping. `examples/archive/Dynamic_load/Si/main.jl`
    and `examples/3D_implantation_B-in-Si/main.jl` were read for the sizes of
    its `pMax` (cited in the doc comments; the defaults here differ). No
    DISPLATH data file is copied.

```text
MIT License

Copyright (c) 2024 裴幂許Permission

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
