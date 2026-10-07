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
