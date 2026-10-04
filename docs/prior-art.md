# Prior art

This survey of open and closed codes was taken on 2026-10-04 and is the
starting point for this project. Activity figures come from the GitHub API on
that date. **[unverified]** marks claims that were not checked against a
primary source. For what each license tier allows us to do with a code, see
[`../CONTRIBUTING.md`](../CONTRIBUTING.md).

## Ions: stopping, range, implantation

| Code | License | Lang | Activity (2026-10) | Strengths | Gaps |
|---|---|---|---|---|---|
| SRIM/TRIM (Ziegler) | Closed, non-commercial | VB6/Win | SRIM-2013, unmaintained | De-facto standard; very large stopping set | Amorphous and static targets; 1D layers; damage-count controversy (Stoller 2013); overestimates S_e for slow heavy ions in light targets |
| RustBCA (lcpp-org) | GPL-3.0 | Rust | v3.0.0 2026-08, active, one maintainer | 0D–3D mesh geometry; many potentials; seeded RNG; Python/C/Fortran bindings | No crystals or channeling; no dynamic composition; not on PyPI; aimed at fusion sputtering |
| OpenTRIM (ir2-lab) | MIT | C++/Qt | v1.2.0-rc3 2026-09, very active | GUI, CLI, library, Python; 3D grid; damage focus | Electronic stopping is SRIM-2013 tables; appears amorphous/static **[unverified]** |
| iradina | GPL-3.0 | C | 2024-12 | Static 3D, nanostructures | Copyleft; quiet |
| IM3D | GPL-3.0 | C | Dormant since 2021 | Fast 3D CSG | Copyleft; dormant |
| DISPLATH | MIT | Julia | 2026-06, 3 stars | BCA on explicit atoms: crystals, 2D materials | No paper; unproven |
| Corteo | GPL | C | Low **[unverified]** | Ion-beam-analysis spectra | Not for implants |
| SDTrimSP | Paid licence | Fortran | Licensed by IPP | Dynamic 1D/2D/3D | Not open |
| TRIDYN / TRI3DYN | Non-commercial licence | Fortran | HZDR | Dynamic composition | Not open |
| IMSIL (Hobler) | Free on request, closed | Fortran | TU Wien | Crystalline, damage accumulation, 2D | Not open |
| Crystal-TRIM | Executable under licence | — | HZDR | Crystalline Si/Ge/C | Not open |
| MARLOWE | Export-controlled (RSICC) | Fortran | ORNL | Crystalline BCA | Not open |
| Sentaurus MC / Victory | Commercial | — | — | Calibrated crystalline implant | Commercial |

## Stopping-power libraries and data

| Source | License | Notes |
|---|---|---|
| libdEdx (APTG) | GPL-3.0 | PSTAR/ASTAR/MSTAR/ICRU49/73/Bethe |
| CATIMA | AGPL-3.0 | ATIMA-style high-energy ions; pycatima on PyPI |
| JIBAL | GPL-2.0 | Ion-beam-analysis stopping and straggling |
| ESPNN | MIT | Neural network trained on the IAEA database |
| IAEA stopping database | Open access; reuse terms **[unverified]** | Experimental; about 1,500 systems, under 30% well known |
| NIST PSTAR/ASTAR/ESTAR | NIST SRD: may be copyrighted (15 USC 290e) | Treat as Tier C until reviewed |
| ICRU 49/73/90 | Sold reports | Tier C |
| ZBL, Lindhard-Scharff, Oen-Robinson, Bethe-Bloch, Biersack-Varelas | Published formulas | Freely implementable |

## Electrons: low-energy transport

| Code | License | Notes |
|---|---|---|
| Nebula + cstool (TU Delft) | BSD-3 | 0–50 keV; dielectric-function inelastic, Mott elastic; CUDA. Stalled since 2021–24 |
| Geant4 MicroElec | Geant4 licence | Reference low-energy Si/SiO₂/metal models, locked inside Geant4 |
| PENELOPE / penEasy | Permissive notice, registration | Fortran; best-validated general e⁻ physics, ~50 eV–1 GeV |
| PenRed | AGPL-3.0 | Modern C++ PENELOPE |
| EGSnrc | AGPL-3.0 | Hard 1 keV cutoff, so wrong for this regime |
| CASINO | Closed | SEM/EBL standard on Windows |
| PHITS | Closed, free | ETS/ETSART track structure for semiconductors |
| GenISys BEAMER/TRACER | Commercial | Industry e-beam lithography proximity correction (PEC) |
| PyEPICS | BSD-3 | EEDL/EADL/EPDL to HDF5 |

**No Rust electron transport code exists**, and crates.io has nothing in this
space.

## What open source does not provide (lindhard's targets)

1. A permissive, embeddable, fast ion BCA engine with independently sourced
   stopping data.
2. An open, validated **crystalline** implant simulator: channeling under tilt
   and twist, screen oxides, and damage accumulation leading to amorphization.
3. A permissive low-energy **electron** engine that produces e-beam
   lithography point-spread functions, SE/BSE yields and generation volumes.
4. Targets whose composition changes with fluence (dynamic composition).
