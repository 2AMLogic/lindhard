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
//!   compares it with the direct DIIMFP (`FullPenn::diimfp_per_m_ev` at a
//!   tolerance of 1e-7) at seven losses spread in `ln W` in every row.
//! * `imfp` times `imfp_and_stopping` at every grid energy (one thread).

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
            "usage: <elf.toml> <E_min> <E_max> <points per decade> [build|grid|imfp] [stride]"
                .into(),
        );
    }
    let elf = OpticalElf::from_toml_file(&a[1])?;
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
            for &e in &grid {
                let t = Instant::now();
                let p = fp.imfp_and_stopping(e)?;
                println!(
                    "{e:.6} {:.12e} {:.12e} {:.3}",
                    p.inverse_imfp_per_m,
                    p.stopping_ev_per_m,
                    t.elapsed().as_secs_f64()
                );
            }
            eprintln!("total {:.2} s", t0.elapsed().as_secs_f64());
        }
        "grid" => {
            let fp = FullPenn::new(elf.clone());
            let g = fp.diimfp_grid(&grid)?;
            eprintln!(
                "grid: {} loss nodes, {} momentum panels, {} unresolved cells, {:.2} s",
                g.len(),
                g.panel_count(),
                g.unresolved_cells(),
                t0.elapsed().as_secs_f64()
            );
            let reference = fp.with_relative_tolerance(1e-7)?;
            let w_lo = elf.energy_ev()[0];
            let mut worst: f64 = 0.0;
            for &e in &grid {
                for j in 0..7 {
                    let f = (j as f64 + 0.37) / 7.0;
                    let w = (w_lo.ln() + f * (e.ln() - w_lo.ln())).exp();
                    let direct = reference.diimfp_per_m_ev(e, w)?;
                    let from_grid = g.diimfp_per_m_ev(e, w).ok_or("not covered")?;
                    let rel = from_grid / direct - 1.0;
                    worst = worst.max(rel.abs());
                    println!("{e:.6} {w:.6e} {direct:.9e} {from_grid:.9e} {rel:+.3e}");
                }
            }
            eprintln!("max |grid/direct - 1| {worst:.2e}");
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
