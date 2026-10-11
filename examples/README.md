# Examples

Input files for the `lindhard` command. The schema is in
[`../docs/cli.md`](../docs/cli.md); energies are in eV, lengths in nm and
angles in degrees.

| File | What it shows |
|---|---|
| [`b_5keV_si.toml`](b_5keV_si.toml) | 5 keV B into amorphous Si, 7° off normal: a pure-element substrate, a chosen `E_d`, and a dual-Pearson fit of the depth profile |
| [`as_50keV_si_sio2.toml`](as_50keV_si_sio2.toml) | 50 keV As into Si through a 10 nm SiO₂ screen: a compound material, a finite layer over a substrate, and per-element energy overrides |
| [`ar_1keV_cu.toml`](ar_1keV_cu.toml) | 1 keV Ar onto Cu at normal incidence: backscattering and sputtering, with escape spectra and lateral profiles |
| [`b_5keV_si_crystal.toml`](b_5keV_si_crystal.toml) | 5 keV B into crystalline Si on the explicit (001) lattice with an in-plane `[010]` reference, 7° tilt, 22° azimuth and 300 K Debye vibration (`[[crystal]]`); not a validated channeling prediction |
| [`b_5keV_si_crystal_divergence.toml`](b_5keV_si_crystal_divergence.toml) | The crystalline-Si example with a `[beam.divergence]` table: Gaussian spread of 0.1° per plane about the nominal tilt; an input capability, not a validated channeling prediction |
| [`as_50keV_sio2_on_si_crystal.toml`](as_50keV_sio2_on_si_crystal.toml) | 50 keV As through a 10 nm amorphous SiO₂ layer into a crystalline Si substrate: `[[crystal]]` names only layer 1, the oxide stays amorphous; shows layer assignment, not interface physics |
| [`dynamic/as_1keV_si_film.toml`](dynamic/as_1keV_si_film.toml) | 1 keV As into a 10 nm Si film on Si with a `[dynamic]` table: the target composition and thickness change with fluence (adaptive steps, 2 nm slabs); writes `dynamic_*.csv` |
| [`electron/e_10keV_si.toml`](electron/e_10keV_si.toml) | 10 keV electrons into Si with an `[electron]` table: Mott elastic (stand-in potential), single-pole Penn inelastic, Kieft-Bosch secondaries, the surface barrier and a cylindrical deposition grid; writes `electron_summary.json` and `electron_*.csv`. Its optical ELF ([`electron/synthetic_plasmon_elf.toml`](electron/synthetic_plasmon_elf.toml)) and band parameters are **synthetic**, not Si data |

```sh
cargo run -p lindhard-cli -- check examples/b_5keV_si.toml
cargo run --release -p lindhard-cli -- run examples/b_5keV_si.toml --out out/b_5keV_si
```

The `E_d` values, and the O energies in the SiO₂ example, are illustrative
model parameters, not recommendations.
