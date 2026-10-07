//! Run a resolved input with a `[dynamic]` section: the fluence stepping loop
//! of `lindhard::ion::dynamic` on a rayon pool of the requested size. Shared
//! by the CLI and anything else that calls it, so the numbers are the same.

use std::collections::BTreeSet;
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use lindhard::elements::element;
use lindhard::geometry::Stack;
use lindhard::input::{DynamicSpec, RelaxationChoice, Resolved};
use lindhard::ion::bca::Beam;
use lindhard::ion::dynamic::{
    atomic_volume_from_density, CompositionGrid, DynamicConfig, DynamicRun, Relaxation, StepPolicy,
    Yields,
};
use lindhard::ion::potential::Potential;
use lindhard::ion::scattering::ScatteringTable;

use crate::output::RunInfo;

/// One slab of the target at one moment.
#[derive(Debug, Clone, PartialEq)]
pub struct SlabRow {
    /// Depth of the front of the slab from the fixed front surface, m.
    pub front_m: f64,
    /// Depth of the back of the slab, m.
    pub back_m: f64,
    /// Areal inventory `(Z, atoms/m^2)`, sorted by `Z`.
    pub inventory: Vec<(u8, f64)>,
}

/// The state after one step (step 0: before any ions).
#[derive(Debug, Clone, PartialEq)]
pub struct StepRow {
    /// Step number; 0 is the initial target.
    pub step: u64,
    /// Global index of the first primary of the step (0 for step 0).
    pub first_index: u64,
    /// Primaries in the step.
    pub ions: u64,
    /// Primaries used up to and including this step.
    pub ions_done: u64,
    /// Fluence delivered up to and including this step, ions/m^2.
    pub fluence_m2: f64,
    /// Attempts at this step (more than 1: the adaptive bound rejected some).
    pub attempts: u32,
    /// Largest relative composition change of the step.
    pub max_change: f64,
    /// Removals capped at what a slab held.
    pub clamped: u32,
    /// Slabs that emptied in this step.
    pub removed_slabs: usize,
    /// Event counts summed over steps 1..=this one.
    pub cumulative: Yields,
    /// The slabs after the step.
    pub slabs: Vec<SlabRow>,
}

/// Everything a finished dynamic run produced.
pub struct DynamicSimulation {
    /// Atomic numbers of every species the run can contain, sorted: the beam
    /// species and every element of the target.
    pub species: Vec<u8>,
    /// Step 0 (the initial target) and every accepted step.
    pub steps: Vec<StepRow>,
    /// Attempts that were rejected and retried, over the whole run.
    pub rejected_attempts: u64,
    /// The scattering table the run used.
    pub table: ScatteringTable,
    /// Thread count and timings (the only thread-dependent part of the output).
    pub info: RunInfo,
}

fn slab_rows(grid: &CompositionGrid) -> Vec<SlabRow> {
    let mut front = 0.0;
    grid.thicknesses_m()
        .iter()
        .enumerate()
        .map(|(i, &t)| {
            let row = SlabRow {
                front_m: front,
                back_m: front + t,
                inventory: grid.inventory(i).expect("slab exists"),
            };
            front += t;
            row
        })
        .collect()
}

/// Split every finite layer into slabs of at most `slab_m`.
fn split_stack(stack: &Stack, slab_m: Option<f64>) -> Result<Stack> {
    let mut layers = Vec::new();
    let mut substrate = None;
    for l in stack.layers() {
        if !l.thickness_m().is_finite() {
            substrate = Some(l.material().clone());
            continue;
        }
        let n = match slab_m {
            Some(s) => (l.thickness_m() / s - 1e-9).ceil().max(1.0) as usize,
            None => 1,
        };
        if n > 100_000 {
            bail!("dynamic.slab_nm would make {n} slabs from one layer; use a larger value");
        }
        let t = l.thickness_m() / n as f64;
        for _ in 0..n {
            layers.push((l.material().clone(), t));
        }
    }
    Stack::new(layers, substrate).context("building the slabs")
}

fn relaxation(d: &DynamicSpec, species: &[u8]) -> Result<Relaxation> {
    match d.relaxation {
        RelaxationChoice::FixedNumberDensity => {
            Relaxation::fixed_number_density(d.number_density_cm3.expect("validated") * 1e6)
                .map_err(Into::into)
        }
        RelaxationChoice::IdealMixing => {
            let mut v = Vec::new();
            for &z in species {
                let sym = element(z).expect("valid").symbol;
                let vol = match d.atomic_volume_nm3.get(sym) {
                    Some(nm3) => nm3 * 1e-27,
                    None => atomic_volume_from_density(z).map_err(|e| {
                        anyhow!("{e}; set dynamic.atomic_volume_nm3.{sym} (nm^3/atom)")
                    })?,
                };
                v.push((z, vol));
            }
            Relaxation::ideal_mixing(&v).map_err(Into::into)
        }
    }
}

/// Run `r` (which must have a `[dynamic]` section) on `threads` workers
/// (`None`: all available cores). The thread count never changes the results.
pub fn simulate_dynamic(r: &Resolved, threads: Option<usize>) -> Result<DynamicSimulation> {
    let d = r
        .input
        .dynamic
        .as_ref()
        .context("the input has no [dynamic] section")?;
    if threads == Some(0) {
        bail!("threads must be at least 1");
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.unwrap_or(0))
        .build()
        .context("building the thread pool")?;

    let t0 = Instant::now();
    let z2 = r.layers[0].material.components()[0].z();
    let pot = Potential::new(r.screening, f64::from(r.beam.ion.z()), f64::from(z2))
        .with_length(r.screening_length);
    let table = ScatteringTable::build(&pot, &r.table_spec);
    let table_build_s = t0.elapsed().as_secs_f64();

    let species: Vec<u8> = {
        let mut s = BTreeSet::new();
        s.insert(r.beam.ion.z());
        for l in &r.layers {
            s.extend(l.material.components().iter().map(|c| c.z()));
        }
        s.into_iter().collect()
    };

    let stack = split_stack(&r.stack, d.slab_nm.map(|nm| nm * 1e-9))?;
    let mut grid = CompositionGrid::from_stack(&stack, relaxation(d, &species)?)
        .context("building the grid")?;
    // An element that enters a slab by implantation or recoil mixing needs the
    // energies the engine requires of every element in a layer.
    for &z in &species {
        let sym = element(z).expect("valid").symbol;
        let o = d
            .energies
            .get(sym)
            .or_else(|| r.input.physics.energies.get(sym));
        grid.seed_energies(
            z,
            o.and_then(|o| o.e_d_ev),
            o.and_then(|o| o.e_b_ev),
            o.and_then(|o| o.e_s_ev),
        );
    }

    let policy = match d.max_change {
        None => StepPolicy::Fixed {
            ions_per_step: d.ions_per_step,
        },
        Some(max_change) => StepPolicy::Adaptive {
            max_ions_per_step: d.ions_per_step,
            min_ions_per_step: d.min_ions_per_step,
            max_change,
        },
    };
    let stopping = r.stopping_model();
    let beam = Beam {
        count: r.input.run.ions,
        ..r.beam
    };
    let mut run = DynamicRun::new(
        grid,
        beam,
        r.config,
        &*stopping,
        &table,
        DynamicConfig {
            fluence_m2: d.fluence_cm2 * 1e4,
            policy,
        },
    )
    .context("setting up the fluence loop")?;

    let mut steps = vec![StepRow {
        step: 0,
        first_index: 0,
        ions: 0,
        ions_done: 0,
        fluence_m2: 0.0,
        attempts: 0,
        max_change: 0.0,
        clamped: 0,
        removed_slabs: 0,
        cumulative: Yields::default(),
        slabs: slab_rows(run.grid()),
    }];
    let mut rejected = 0u64;
    let t1 = Instant::now();
    pool.install(|| {
        while let Some(rec) = run.step().with_context(|| {
            format!(
                "fluence step {} (primaries from {})",
                run.steps_done() + 1,
                run.ions_done()
            )
        })? {
            rejected += u64::from(rec.attempts - 1);
            steps.push(StepRow {
                step: rec.step,
                first_index: rec.first_index,
                ions: rec.ions,
                ions_done: run.ions_done(),
                fluence_m2: run.fluence_done_m2(),
                attempts: rec.attempts,
                max_change: rec.max_change,
                clamped: rec.clamped,
                removed_slabs: rec.removed_slabs.len(),
                cumulative: run.cumulative().clone(),
                slabs: slab_rows(run.grid()),
            });
        }
        anyhow::Ok(())
    })?;
    let transport_s = t1.elapsed().as_secs_f64();
    let info = RunInfo {
        threads: pool.current_num_threads(),
        table_build_s,
        transport_s,
        ions_per_s: run.ions_done() as f64 / transport_s.max(f64::MIN_POSITIVE),
    };
    Ok(DynamicSimulation {
        species,
        steps,
        rejected_attempts: rejected,
        table,
        info,
    })
}
