//! Wall time and accuracy checks of a full-Penn inelastic table build on a
//! measured optical ELF (issue #256).
//!
//! ```text
//! cargo run --release -p lindhard --example penn_full_build_time -- \
//!     <elf.toml> <E_min eV> <E_max eV> <points per decade> [build|grid|imfp] [stride]
//! ```
//!
//! for example `validation/data/optical/al_elf_hagemann1975.toml 51 30000 20
//! build` (the grid of the backscatter comparison). Threads come from
//! `RAYON_NUM_THREADS`. `stride` keeps every n-th grid energy (the ends are
//! kept), for quick probes.
//!
//! * `build` times one `build_inelastic_table_for_model` with the full Penn
//!   model (the DIIMFP grid, the IMFP pass, the loss densities and the
//!   inverse CDFs), then prints, per row, the inverse IMFP and the stopping
//!   power of the table (`λ⁻¹ ⟨W⟩`) against the model's own
//!   `imfp_and_stopping` (computed after the timed build).
//! * `grid` times `FullPenn::diimfp_grid` for the grid energies and
//!   compares it with the direct DIIMFP (`FullPenn::diimfp_evaluation` at a
//!   tolerance of 1e-7) at seven losses spread in `ln W` in every row. The
//!   grid's counts of unconverged momentum panels and loss evaluations
//!   (#298) are printed with those of the direct DIIMFP (the last column of
//!   each line: loss evaluations whose `ω_p` integral did not converge,
//!   #368).
//! * `imfp` times `imfp_and_stopping_evaluation` at every grid energy (one
//!   thread); the last two columns are the unconverged cells of the outer
//!   `ω_p` integral and their error relative to `λ⁻¹` (#298, #368).
//!
//! Two more modes look at single evaluations of the expanded loss function
//! `Im[-1/ε(q, ω)]` (#298), with the work of its `ω_p` integral
//! (`FullPenn::loss_function_evaluation`):
//!
//! ```text
//! ... penn_full_build_time -- <elf.toml> point <q a.u.> <loss eV> [tolerance]
//! ... penn_full_build_time -- <elf.toml> scan <n_q> <tolerance> <loss eV>...
//! ```
//!
//! `point` times one evaluation (for example
//! `validation/data/optical/si_elf_yang2019.toml point 1.480290203 160`);
//! `scan` evaluates at `n_q` log-spaced momenta from 0.005 to 40 a.u. at
//! each loss and prints the largest work and the number of evaluations that
//! did not converge.

use lindhard::electron::data::OpticalElf;
use lindhard::electron::elastic::table::log_energy_grid;
use lindhard::electron::inelastic::table::{
    build_inelastic_table_for_model, stopping_power_ev_per_m, InelasticTableOptions,
};
use lindhard::electron::inelastic::{FullPenn, PennAlgorithm, PennInelastic};
use lindhard::material::Material;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 5 {
        return Err(
            "usage: <elf.toml> <E_min> <E_max> <points per decade> [build|grid|imfp] [stride], \
             <elf.toml> point <q a.u.> <loss eV> [tolerance], \
             or <elf.toml> scan <n_q> <tolerance> <loss eV>..."
                .into(),
        );
    }
    let elf = OpticalElf::from_toml_file(&a[1])?;
    match a.get(2).map(String::as_str) {
        Some("point") => return point(elf, &a[3..]),
        Some("scan") => return scan(elf, &a[3..]),
        _ => {}
    }
    let (lo, hi, ppd): (f64, f64, f64) = (a[2].parse()?, a[3].parse()?, a[4].parse()?);
    let mode = a.get(5).map_or("build", String::as_str);
    let stride: usize = a.get(6).map_or(Ok(1), |s| s.parse())?;
    let full_grid = log_energy_grid(lo, hi, ppd)?;
    let n = full_grid.len();
    let grid: Vec<f64> = full_grid
        .iter()
        .enumerate()
        .filter(|(i, _)| i % stride == 0 || *i == n - 1)
        .map(|(_, &e)| e)
        .collect();
    eprintln!(
        "{}: {} energies {lo}..{hi} eV, {} threads",
        elf.material(),
        grid.len(),
        rayon::current_num_threads()
    );
    let t0 = Instant::now();
    match mode {
        "imfp" => {
            let fp = FullPenn::new(elf);
            let (mut cells, mut points, mut worst_err) = (0u64, 0usize, 0.0f64);
            for &e in &grid {
                let t = Instant::now();
                let ev = fp.imfp_and_stopping_evaluation(e)?;
                let p = ev.point;
                cells += ev.unconverged_cells;
                points += usize::from(ev.unconverged_cells > 0);
                worst_err = worst_err.max(ev.unconverged_error);
                println!(
                    "{e:.6} {:.12e} {:.12e} {:.3} {} {:.1e}",
                    p.inverse_imfp_per_m,
                    p.stopping_ev_per_m,
                    t.elapsed().as_secs_f64(),
                    ev.unconverged_cells,
                    ev.unconverged_error
                );
            }
            eprintln!(
                "total {:.2} s; outer omega_p integral: {cells} unconverged cells at {points} of \
                 {} energies, largest relative error {worst_err:.1e}",
                t0.elapsed().as_secs_f64(),
                grid.len()
            );
        }
        "grid" => {
            let fp = FullPenn::new(elf.clone());
            let g = fp.diimfp_grid(&grid)?;
            eprintln!(
                "grid: {} loss nodes, {} momentum panels ({} unconverged, largest \
                 relative error {:.1e}; {} loss evaluations unconverged), {} unresolved \
                 cells, {:.2} s",
                g.len(),
                g.panel_count(),
                g.unconverged_profile_panels(),
                g.largest_unconverged_profile_error(),
                g.unconverged_profile_loss_evaluations(),
                g.unresolved_cells(),
                t0.elapsed().as_secs_f64()
            );
            if let Some((lo, hi)) = g.unresolved_loss_range_ev() {
                eprintln!("unresolved cells between {lo:.6e} and {hi:.6e} eV loss");
            }
            let reference = fp.with_relative_tolerance(1e-7)?;
            let w_lo = elf.energy_ev()[0];
            let mut worst: f64 = 0.0;
            let (mut unconverged, mut points, mut total) = (0u64, 0usize, 0usize);
            for &e in &grid {
                for j in 0..7 {
                    let f = (j as f64 + 0.37) / 7.0;
                    let w = (w_lo.ln() + f * (e.ln() - w_lo.ln())).exp();
                    let d = reference.diimfp_evaluation(e, w)?;
                    let direct = d.value;
                    unconverged += d.unconverged_loss_evaluations;
                    points += usize::from(d.unconverged_loss_evaluations > 0);
                    total += 1;
                    let from_grid = g.diimfp_per_m_ev(e, w).ok_or("not covered")?;
                    let rel = from_grid / direct - 1.0;
                    worst = worst.max(rel.abs());
                    println!(
                        "{e:.6} {w:.6e} {direct:.9e} {from_grid:.9e} {rel:+.3e} {}",
                        d.unconverged_loss_evaluations
                    );
                }
            }
            eprintln!(
                "max |grid/direct - 1| {worst:.2e}; direct DIIMFP: {unconverged} unconverged \
                 loss evaluations at {points} of {total} points"
            );
        }
        "build" => {
            let model = PennInelastic::new(PennAlgorithm::Full, elf);
            // the material only names the table; the ELF decides the rows
            let mat = Material::from_atom_fractions(&[(13, 1.0)], None)?;
            let t = build_inelastic_table_for_model(
                &model,
                &mat,
                &InelasticTableOptions::new(grid.clone()),
            )?;
            let build = t0.elapsed().as_secs_f64();
            eprintln!(
                "table: {} rows, {} probability points, built in {build:.2} s",
                t.energy_ev().len(),
                t.probability().len()
            );
            let mut worst: f64 = 0.0;
            for (i, &e) in t.energy_ev().iter().enumerate() {
                let s_table = stopping_power_ev_per_m(&t, i).unwrap_or(0.0);
                let p = model.imfp_and_stopping(e)?;
                let rel = s_table / p.stopping_ev_per_m - 1.0;
                worst = worst.max(rel.abs());
                println!(
                    "{e:.6} {:.12e} {s_table:.9e} {:.9e} {rel:+.3e}",
                    t.inverse_mfp_per_m()[i],
                    p.stopping_ev_per_m
                );
            }
            eprintln!("max |S_table/S_model - 1| {worst:.2e}");
        }
        other => return Err(format!("unknown mode {other:?}").into()),
    }
    Ok(())
}

/// `point <q a.u.> <loss eV> [tolerance]`: one evaluation of the expanded
/// loss function, timed (one thread).
fn point(elf: OpticalElf, a: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let (q_au, w): (f64, f64) = (
        a.first().ok_or("missing q")?.parse()?,
        a.get(1).ok_or("missing loss")?.parse()?,
    );
    let tol: f64 = a.get(2).map_or(Ok(1e-4), |s| s.parse())?;
    let fp = FullPenn::new(elf).with_relative_tolerance(tol)?;
    let q = q_au / lindhard::constants::BOHR_RADIUS;
    let t = Instant::now();
    let e = fp.loss_function_evaluation(q, w);
    let secs = t.elapsed().as_secs_f64();
    println!(
        "q {q_au} a.u., loss {w} eV, tol {tol:e}: Im[-1/eps] = {:.12e}, {secs:.4} s, \
         {} Lindhard evaluations, largest interval {} rules, converged {} \
         ({} unconverged cells, error {:.2e})",
        e.value,
        e.lindhard_evaluations,
        e.largest_interval_rules,
        e.converged(),
        e.unconverged_cells,
        e.unconverged_error
    );
    Ok(())
}

/// `scan <n_q> <tolerance> <loss eV>...`: the loss function at `n_q`
/// log-spaced momenta from 0.005 to 40 a.u. at each loss, with the largest
/// work and the number of unconverged evaluations (one thread).
fn scan(elf: OpticalElf, a: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let n: usize = a.first().ok_or("missing n_q")?.parse()?;
    let tol: f64 = a.get(1).ok_or("missing tolerance")?.parse()?;
    let fp = FullPenn::new(elf).with_relative_tolerance(tol)?;
    for w in &a[2..] {
        let w: f64 = w.parse()?;
        let t = Instant::now();
        let (mut evals, mut worst_evals, mut worst_rules, mut unconverged) = (0u64, 0u64, 0u64, 0);
        let mut worst_q = 0.0;
        for i in 0..n {
            let q_au = 0.005 * (40.0f64 / 0.005).powf(i as f64 / (n - 1).max(1) as f64);
            let e = fp.loss_function_evaluation(q_au / lindhard::constants::BOHR_RADIUS, w);
            evals += e.lindhard_evaluations;
            if e.lindhard_evaluations > worst_evals {
                worst_evals = e.lindhard_evaluations;
                worst_q = q_au;
            }
            worst_rules = worst_rules.max(e.largest_interval_rules);
            unconverged += usize::from(!e.converged());
        }
        println!(
            "loss {w} eV, tol {tol:e}: {n} momenta, {:.3} s, mean {:.0} evaluations, \
             most {worst_evals} (q {worst_q:.6} a.u.), largest interval {worst_rules} rules, \
             {unconverged} unconverged",
            t.elapsed().as_secs_f64(),
            evals as f64 / n as f64
        );
    }
    Ok(())
}
