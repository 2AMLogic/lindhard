# Examples

Input files for the `lindhard` command. The schema is in
[`../docs/cli.md`](../docs/cli.md); energies are in eV, lengths in nm and
angles in degrees.

| File | What it shows |
|---|---|
| [`b_5keV_si.toml`](b_5keV_si.toml) | 5 keV B into amorphous Si, 7° off normal: a pure-element substrate, a chosen `E_d`, and a dual-Pearson fit of the depth profile |
| [`as_50keV_si_sio2.toml`](as_50keV_si_sio2.toml) | 50 keV As into Si through a 10 nm SiO₂ screen: a compound material, a finite layer over a substrate, and per-element energy overrides |
| [`ar_1keV_cu.toml`](ar_1keV_cu.toml) | 1 keV Ar onto Cu at normal incidence: backscattering and sputtering, with escape spectra and lateral profiles |

```sh
cargo run -p lindhard-cli -- check examples/b_5keV_si.toml
cargo run --release -p lindhard-cli -- run examples/b_5keV_si.toml --out out/b_5keV_si
```

The `E_d` values, and the O energies in the SiO₂ example, are illustrative
model parameters, not recommendations.
