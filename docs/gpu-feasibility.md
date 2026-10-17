# GPU feasibility study

**Decision: defer.** No GPU backend now. The study could not measure a
hardware GPU (the host that ran it has none), so the throughput half of the
case is open, and the determinism half says a GPU path cannot be bit-identical
to the CPU engine. Revisit when the trigger at the end of this page is met.

`architecture.md` calls a GPU backend "a later, measured decision, not an
assumption". This page is that decision's record. It is a study, not a
commitment.

## The kernel

The prototype is the engine's hot inner kernel: one bilinear lookup in the
log-`tan(theta/2)` scattering table (`ScatteringTable::half_angle_tan`) plus
one Lindhard-Scharff stopping evaluation (`k_L sqrt(eps)`, Lindhard and
Scharff, Phys. Rev. 124, 128 (1961)), per `(eps, beta)` point. Inputs are the
same deterministic log-uniform sample as `benches/scattering.rs`, ZBL
universal screening, Si on Si, 32 points per decade (a 161 x 161 table).

It is a separate crate, [`lindhard-gpu-proto/`](../lindhard-gpu-proto/). It is
not a workspace member (its manifest carries an empty `[workspace]` table),
so default builds, `cargo test --workspace`, CI, `Cargo.lock` and `deny.toml`
never see it or wgpu. Run it explicitly:

```sh
cargo run --release --manifest-path lindhard-gpu-proto/Cargo.toml
```

## Measurement

Method: median of at least 3 repeats (CPU) or 7 repeats (wgpu) after one
warm-up dispatch. "cpu 1T" is the `lindhard` library, f64, one thread. "wgpu
kernel" is dispatch to completion with data already resident; "wgpu e2e"
adds the input upload and the output read-back. `max rel dev` is the largest
relative difference of the f32 shader output from the f64 CPU output.

Hardware: Intel Xeon Platinum 8488C, 14 logical CPUs visible, a shared AWS
dispatch worker at load average about 3.6, rustc stable, wgpu 30.0.1,
2026-10-16. **There is no GPU on this host.** wgpu selected the only adapter
it found, `llvmpipe` (Mesa's software rasteriser, `DeviceType::Cpu`, GL
backend), which runs the WGSL shader on the host CPU cores.

| points | cpu 1T (Mpt/s) | wgpu kernel on llvmpipe (Mpt/s) | wgpu e2e on llvmpipe (Mpt/s) | max rel dev |
|---:|---:|---:|---:|---:|
| 16 384 | 28.3 | 194 | 112 | 9.8e-6 |
| 262 144 | 28.0 | 585 | 241 | 1.1e-5 |
| 2 097 152 | 23.9 | 765 | 126 | 1.2e-5 |

What these numbers do and do not show:

* The wgpu columns are **not GPU throughput**. They are a multi-threaded
  software renderer on 14 CPU threads against one scalar CPU thread, so the
  ratio (up to about 30x) is a thread-count and SIMD ratio, not evidence
  about any GPU. The program prints a warning when the adapter is not a GPU.
* They do show that the kernel is expressible in WGSL, runs correctly (the
  f32 result agrees with f64 to about 1e-5 relative, which is the f32
  rounding of `ln`/`exp` and the table, not a logic difference) and that the
  per-dispatch transfer cost is not negligible: end-to-end throughput is a
  fraction of kernel throughput at every size.
* **No measured GPU speedup exists in this study.** Any number quoted for a
  discrete GPU from this work would be an extrapolation, and none is given.
  Measuring on real hardware is the open item below.
* The kernel is a per-collision lookup. In the engine it sits inside a
  history loop with branching (free flight, recoil bookkeeping, tally
  updates, a cascade that spawns variable numbers of recoils). A kernel-only
  rate is an upper bound on what a ported engine could reach, not a forecast
  of it.

## Comparison

| | wgpu (WGSL) | CubeCL | CUDA (cudarc or similar) |
|---|---|---|---|
| Hardware | Vulkan, Metal, DX12, GL; one code path | CUDA, ROCm/HIP, wgpu; Rust kernels | NVIDIA only |
| Float64 | Not in WGSL; f32 only | Backend dependent; f64 on CUDA | Yes |
| Licence (crate) | MIT OR Apache-2.0 | MIT OR Apache-2.0 (`cubecl`, `cubecl-wgpu`, `cubecl-cuda`) | `cudarc` MIT OR Apache-2.0 |
| Licence (runtime) | Platform graphics drivers | CUDA backend needs NVIDIA's runtime | NVIDIA CUDA toolkit and driver: proprietary EULA, not redistributed by us but required to run |
| Determinism | f32; transcendental functions (`log`, `exp`, `sqrt`) are not required to be correctly rounded and differ between vendors and drivers | As wgpu on wgpu; as CUDA on CUDA | Same binary on same GPU model and driver is repeatable; across architectures `__expf`-class results can differ unless flags pin them |
| Measured throughput | Not measured on a GPU (see above) | Not measured | Not measured (no CUDA device, and no code written) |

### Determinism

The engine's invariant is that results are identical across thread counts. It
holds because each history draws from its own ChaCha8 stream keyed on
`(seed, particle index)` (`lindhard/src/rng.rs`) and because tallies are
merged in chunk-index order, not completion order.

* **RNG streams survive a port.** A counter-based generator needs no shared
  state: a GPU thread can own history `i` and seed ChaCha8 from
  `(seed, i)`, so the random numbers each history sees are the same as on the
  CPU. Philox or ChaCha both work in WGSL/CUDA with integer-only arithmetic.
* **Bit-identical results do not.** The CPU engine is f64 and its math library
  is the platform's. WGSL has no f64, GPU `log`/`exp`/`sqrt` are not
  correctly rounded, and a history is a chaotic sequence of branches (a
  one-ulp change in a scattering angle can change which recoil is created).
  A GPU history would therefore diverge from the CPU history of the same
  index. The measured f32 deviation above is about 1e-5 per lookup, which
  compounds over a cascade. A GPU path could be made repeatable against
  itself (same device, same driver, fixed reduction order) but **cannot be
  made equal to the CPU path**, and it could not promise equality across GPU
  vendors without a software-implemented math library.
* **Reduction order.** Tallies accumulated with atomics are order-dependent in
  floating point. A deterministic GPU path needs per-workgroup partial tallies
  merged in index order, or integer (fixed-point) accumulation.

So a GPU backend would have to be a separate, opt-in engine with its own
statistical validation against the CPU engine (agreement within Monte Carlo
error, not bit-for-bit), and the thread-count determinism test would apply to
the CPU engine only.

### Licences against the clean-room tiers

`CONTRIBUTING.md` tiers apply to *sources we learn from*; this study read no
code from any tier B or C project and ported nothing. Tier A's Nebula (BSD-3)
has a CUDA layout that could be read and ported with attribution; it was not
consulted here.

Dependencies of the prototype were checked on 2026-10-16 with `cargo
metadata` over `lindhard-gpu-proto/Cargo.toml` (150 packages including
`lindhard` itself and its own dependencies). SPDX expressions found:

| count | licence expression |
|---:|---|
| 97 | MIT OR Apache-2.0 |
| 16 | MIT |
| 9 | Apache-2.0 OR MIT |
| 8 | Zlib OR Apache-2.0 OR MIT |
| 5 | Apache-2.0 (`codespan-reporting`, `gl_generator`, `glutin_wgl_sys`, `khronos_api`, `spirv`) |
| 2 | Zlib (`foldhash`, `slotmap`) |
| 2 | MIT OR Apache-2.0 OR Zlib |
| 2 + 2 | MIT/Apache-2.0 and Apache-2.0/MIT (older spelling of the same dual licence) |
| 2 | Unlicense OR MIT |
| 2 | BSD-2-Clause OR Apache-2.0 OR MIT |
| 1 | `(MIT OR Apache-2.0) AND Unicode-3.0` (`unicode-ident`) |
| 1 | **ISC** (`libloading`) |

No copyleft licence appears. One item needs action before wgpu could join the
main tree: **`ISC` is not in `deny.toml`'s allow list** (the others are).
ISC is permissive and MIT-compatible, so it is an allow-list addition, not a
blocker, but the addition belongs in the PR that introduces the dependency.
`cargo deny` was not installed on the host, so the `deny.toml` check was done
by comparing this list to the allow list by hand, not by running the tool.
The CUDA toolkit is the one non-open piece in any CUDA route; it is a runtime
requirement, not tree content.

## Go / no-go / defer

**Defer.**

* *Go* would need a measured speedup on real hardware. This study has none.
* *No-go* is also unearned: the kernel ports cleanly, the dependencies are
  MIT-compatible, and counter-based RNG streams carry over.
* The costs are known without a GPU: an f32-only (wgpu) or vendor-locked
  (CUDA) second engine, no bit-for-bit agreement with the CPU engine, and a
  divergent-branch workload (cascades) that is the hard case for GPUs.

**Trigger to revisit:** a profile of a real workload showing the CPU engine is
the bottleneck at the scale users run, *and* the missing measurement below.

**Open item (needs a human with hardware):** run `lindhard-gpu-proto` on a
discrete GPU (and, for a CUDA comparison, write the equivalent `cudarc`
kernel) and add a row set to the table above. Only the wgpu prototype's
adapter line needs to read as a real GPU for the numbers to count.
