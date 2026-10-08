# Crystal targets: lattice, wafer orientation and beam divergence

This page states the conventions of `lindhard::ion::crystal` (steps 19a and
19b of the M2 crystal plan): how a cubic or hexagonal lattice is described, how the beam direction in
the lab becomes a direction in the crystal, and how beam divergence is
sampled. It is a data model only; the transport engine does not read it yet.
The code docs (`lattice`, `orientation` and `divergence` modules) carry the
same statements next to the code.

![Lab frame, wafer cut and beam angles](img/crystal-orientation.svg)

## Lattice

`Lattice` holds the face-centred cubic primitive (Bravais) vectors, a basis of
sites with their atomic numbers, and the cubic lattice constant together with
the temperature that constant refers to.

| Structure | Space group | Basis (lattice coordinates of the primitive vectors) | Source |
|---|---|---|---|
| Diamond (A4): Si, Ge | Fd-3m (227) | `(1/8, 1/8, 1/8)`, `(7/8, 7/8, 7/8)` | Mehl et al. (2017), p. 616 |
| Zincblende (B3): GaAs, 3C-SiC | F-43m (216) | first species `(0, 0, 0)` (Zn, 4a); second `(1/4, 1/4, 1/4)` (S, 4c) | Mehl et al. (2017), p. 543 |

Primitive vectors: `a1 = (0, a/2, a/2)`, `a2 = (a/2, 0, a/2)`,
`a3 = (a/2, a/2, 0)`, so the primitive cell holds 2 atoms in `a³/4` and the
atom number density is `8/a³`. Reference: M. J. Mehl et al., "The AFLOW
Library of Crystallographic Prototypes: Part 1", Comput. Mater. Sci. 136, S1
(2017), doi:10.1016/j.commatsci.2017.01.017 (arXiv:1607.02532v1).

The presets and their lattice constants (every value has a row in
[`data-provenance.md`](data-provenance.md)):

| Preset | `a` | Temperature | Source |
|---|---|---|---|
| `Lattice::silicon()` | 5.431 020 511(89) Å | 22.5 °C, in vacuum | CODATA 2022, Table XXXIV |
| `Lattice::germanium()` | 5.6576 Å | 25 °C | NBS Circular 539, Vol. I (1953), pp. 18-19 |
| `Lattice::gallium_arsenide()` | 5.652 Å | 25 °C | NBS Monograph 25, Section 3 (1964), p. 33 |
| `Lattice::silicon_carbide_3c()` | 4.3596 Å | 297 K | Ioffe NSM archive, citing Taylor and Jones (1960); primary not opened |

The lattice is rigid: no thermal expansion is applied, and the temperature
is metadata of the cited value. Thermal vibration is the separate Debye
model in `ion::crystal::debye`.

Miller indices `(hkl)` and directions `[uvw]` refer to the conventional cubic
cell. The plane normal is the reciprocal lattice vector
`g = h b1 + k b2 + l b3` with `b_i = 2π (c_j × c_k) / V` (Mehl et al.,
eq. (10)); a direction is `u c1 + v c2 + w c3`. Since `c_i · b_j = 2π δ_ij`
(eq. (9)), `[uvw]` lies in `(hkl)` exactly when `hu + kv + lw = 0` (the zone
law).

## Frames and angles

* **Crystal frame**: Cartesian axes along the cube edges; `[100]` is
  `(1, 0, 0)`.
* **Lab frame**: the frame of the geometry module. `+x` is depth, into the
  target; `y` and `z` are lateral.
* **Wafer cut** `(hkl)`: its normal `n` is the crystal direction of lab `+x`,
  so `(hkl)` names the **inward** normal and a beam at zero tilt travels along
  `n`. Pass `(-h, -k, -l)` to name the outward normal instead; this matters
  only for polar faces.
* **Reference direction** `[uvw]`: a lattice direction in the surface
  (`hu + kv + lw = 0`, checked). Its unit vector `r` is the crystal direction
  of lab `+y` at zero wafer rotation. The third axis is `t = n × r`, lab `+z`.
* **Tilt** `θ`: polar angle of the beam from `+x`, in `[0°, 90°)`.
* **Twist** `φ`: azimuth of the beam about `+x`, in the surface, from `+y`
  towards `+z`. The lab beam is `(cos θ, sin θ cos φ, sin θ sin φ)`, the same
  formula as `ion::bca::Beam::direction` with `polar_rad = θ`,
  `azimuth_rad = φ`.
* **Wafer rotation** `ω`: the wafer turns about `+x` by `ω`, right-handed, so
  `r` moves from `+y` towards `+z`.

The crystal-to-lab rotation is `M = R_x(ω) B`, with `B` the matrix of rows
`n`, `r`, `t`; the lab-to-crystal rotation is `Mᵀ`. The beam in the crystal
frame is

```text
d = cos θ n + sin θ cos(φ − ω) r + sin θ sin(φ − ω) t.
```

**The reference direction is explicit.** Tools differ on what twist is
measured from. Here it is always the `[uvw]` the caller passes; there is no
hidden default. For example, Nordlund, Djurabekova and Hobler, Phys. Rev. B
94, 214109 (2016), doi:10.1103/PhysRevB.94.214109, Sec. II B, describe the
incidence direction on a `[001]` surface by tilt `θ` and twist `ϕ`, with
`[011]` at `θ = 45°, ϕ = 0°` and `[111]` at `θ = 54.73°, ϕ = 45°`, so their
twist zero is an in-plane `<100>` axis. With `n = [001]`, `r = [010]` this
module gives `[011]` and `[-111]` (a `<111>` axis; the paper does not fix the
sense of `ϕ`). A wafer flat or notch direction can be passed as `r` instead.

### Worked example: tilt 7°, twist 22° on (100) Si

With `n = [100]`, `r = [010]`, `t = [001]` and `ω = 0`:
`d = (cos 7°, sin 7° cos 22°, sin 7° sin 22°)
= (0.992546, 0.112995, 0.045653)`. With the reference along a `<110>`,
`r = [011]/√2`, `t = [0, −1, 1]/√2`:
`d = (cos 7°, sin 7° (cos 22° − sin 22°)/√2, sin 7° (cos 22° + sin 22°)/√2)
= (0.992546, 0.047618, 0.112181)`. Both are tested to 1e-12
(`lindhard/tests/crystal.rs`).

## Beam divergence

`Divergence` deflects the central beam direction by a polar angle `ϑ` and an
azimuth `ψ`; `ψ` is uniform on `[0, 2π)`.

| Model | Meaning | CDF of `ϑ` | Draw |
|---|---|---|---|
| `Gaussian { sigma_rad: σ }` | independent normal deviations of s.d. `σ` in two orthogonal planes | `1 − exp(−ϑ²/2σ²)` (Rayleigh) | `ϑ = σ √(−2 ln(1 − U))` |
| `UniformCone { half_angle_rad: α }` | uniform in solid angle for `ϑ ≤ α` | `(1 − cos ϑ)/(1 − cos α)` | `ϑ = 2 asin(√U sin(α/2))` |

The Rayleigh law follows from writing the bivariate normal density in polar
coordinates; the cone law from the solid-angle element `sin ϑ dϑ dψ`. The
Gaussian treats the deviations as plane angles (the small-angle reading).

Draws come from the particle's own random stream, `rng::stream(seed,
index)`: exactly two per particle (`U`, then `V` for `ψ = 2πV`), none for
`Divergence::None`. The sampled directions are therefore bit-identical at
any thread count (tested). A one-sample Kolmogorov-Smirnov test at the 1 %
level checks both models against their CDFs (statistic as in the
NIST/SEMATECH e-Handbook, section 1.3.5.16; p-value from Kolmogorov's
limiting distribution, Marsaglia, Tsang and Wang, J. Stat. Softw. 8(18)
(2003), section 3).

## Hexagonal lattices

Wurtzite and the SiC polytypes 4H and 6H (space group P6₃mc, No. 186) use
the hexagonal primitive vectors of Mehl et al. (2017), pp. 423, 425 and 427:

```text
a1 = (a/2, −√3 a/2, 0),   a2 = (a/2, √3 a/2, 0),   a3 = (0, 0, c).
```

So in the crystal frame `z` is the `c` axis `[0001]` and `x` is `[11-20]`
(`a1 + a2`). The hexagonal cell is primitive and is the conventional cell for
indices; its volume is `(√3/2) a² c`.

| Structure | Prototype | Sites (lattice coordinates of `a1, a2, a3`) | Source |
|---|---|---|---|
| Wurtzite (B4), `Lattice::wurtzite` | `AB_hP4_186_b_b` | first species `(1/3, 2/3, 0)`, `(2/3, 1/3, 1/2)`; second `(1/3, 2/3, u)`, `(2/3, 1/3, 1/2 + u)` | Mehl et al. (2017), p. 425 |
| 4H (B5) | `AB_hP8_186_ab_ab` | Si, C on 2a `(0, 0, z)`, `(0, 0, 1/2 + z)` and 2b `(1/3, 2/3, z)`, `(2/3, 1/3, 1/2 + z)`; four `z` | Mehl et al. (2017), p. 423 |
| 6H (B6) | `AB_hP12_186_ab_a2b` | Si, C on 2a and two sets of 2b; six `z` | Mehl et al. (2017), p. 427 |

`u` is the cation-anion bond length along `[0001]` in units of `c`
(Bernardini, Fiorentini and Vanderbilt, Phys. Rev. B 56, R10024 (1997)).
`Lattice::polytype(z_a, z_b, "ABCB", a, c, u, T)` builds the ideal crystal of
any close-packed stacking: with `n` letters, bilayer `k` has its cation on
column A = `(0, 0)`, B = `(1/3, 2/3)` or C = `(2/3, 1/3)` at height `k c/n`
and its anion straight above at `k c/n + 2u c/n` (our generalisation: `u` in
units of the two-bilayer height, so `3/8` is ideal for every `n`, and `"BC"`
is the wurtzite above). The Ramsdell number is the count of bilayers per `c`.

| Preset | `a` | `c` | Internal parameters | Temperature | Source |
|---|---|---|---|---|---|
| `Lattice::gallium_nitride()` | 3.189 Å | 5.178 Å | `u = 1/4 + a²/(3c²)` (**placeholder**, equal bond lengths) | 300 K | Ioffe NSM archive, citing Qian et al. (1996); primary not opened |
| `Lattice::silicon_carbide_4h()` | 3.08051 Å | 10.0848 Å | `z = 0, 0.18784, 0.24982, 0.43671` | not stated | Bauer et al., Acta Cryst. A 57, 60 (2001), via Mehl et al. (2017), pp. 423, 733 |
| `Lattice::silicon_carbide_6h()` | 3.08129 Å | 15.11976 Å | `z = 0, 0.1254, 0.16675, 0.29215, 0.8335, −0.0415` | not stated | Bauer et al. (2001), via Mehl et al. (2017), pp. 427, 734 |

No measured GaN `u` could be opened, so the GaN preset uses the geometric
value at which all four bonds of an atom have one length (derived in the
`lattice` docs; `3/8` at `c/a = √(8/3)`). It is a placeholder, recorded as a
gap in [`data-provenance.md`](data-provenance.md); pass a measured `u` to
`Lattice::wurtzite` to replace it.

### Four-index (Miller-Bravais) indices

* The basal axes are `a1`, `a2` and `a3' = −(a1 + a2)`; the fourth axis is
  `c`.
* A plane `(hkil)` needs `h + k + i = 0`; `(hkl)` are its Miller indices and
  its normal is `h b1 + k b2 + l b3` (Wikipedia "Miller index", revision
  1378193236, section "Hexagonal and rhombohedral structures").
* A direction `[uvtw]` is the vector `u a1 + v a2 + t a3' + w c` with
  `u + v + t = 0`; substituting `a3'` gives the three-index direction
  `[u − t, v − t, w]` (derived). The four-index zone law
  `hu + kv + it + lw = 0` is then the three-index one.

So `[0001]` is `(0, 0, 1)` and `[11-20]` is `(1, 0, 0)` in the crystal frame
(tested to 1e-12 for all three presets). `Orientation::new_miller_bravais`
takes the wafer cut `(hkil)` and the reference `[uvtw]`, for example
`(0001)` with `[11-20]` (c-plane) or `(11-20)` with `[0001]` (a-plane), and
otherwise follows the conventions above.

The lattice search walks rectangular cells: the cube for cubic lattices, and
for hexagonal ones the orthohexagonal cell spanned by `a1 + a2`, `a2 − a1`
and `c` (edges `a`, `√3 a`, `c`), which holds each basis site twice.

## Not here yet

The use of the collision search in the transport engine (including thermal
displacements) and any input-file schema for crystals are later steps.
D. S. Gemmell, Rev. Mod. Phys. 46, 129 (1974), is the background reference
of #19 and #179 for channeling; it was not opened for either step and
nothing here rests on it.
