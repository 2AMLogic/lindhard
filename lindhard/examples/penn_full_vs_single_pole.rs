//! Side-by-side stopping power and inverse IMFP of the single-pole and the
//! full Penn algorithms for the synthetic Drude plasmon (E_p = 20 eV,
//! γ = 5 eV; not physical data).
//!
//! ```text
//! cargo run --release -p lindhard --example penn_full_vs_single_pole
//! ```
//!
//! The full model costs about a second per energy in a release build.

use lindhard::electron::inelastic::{
    DrudeLorentz, DrudeLorentzOscillator, PennAlgorithm, PennInelastic,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let elf = DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])?.to_optical_elf(
        "synthetic Drude plasmon",
        1e-2,
        1e5,
        1000,
    )?;
    let spa = PennInelastic::new(PennAlgorithm::SinglePole, elf.clone());
    let fpa = PennInelastic::new(PennAlgorithm::Full, elf);
    println!("single pole: {}", spa.model_identity());
    println!("full:        {}\n", fpa.model_identity());
    println!(
        "{:>9} {:>14} {:>14} {:>8} {:>14} {:>14} {:>8}",
        "E / eV", "S_spa eV/nm", "S_full eV/nm", "ratio", "lam_spa nm", "lam_full nm", "ratio"
    );
    for e in [
        20.0, 30.0, 50.0, 100.0, 200.0, 500.0, 1e3, 3e3, 1e4, 2e4, 5e4,
    ] {
        let (a, b) = (spa.imfp_and_stopping(e)?, fpa.imfp_and_stopping(e)?);
        println!(
            "{e:9.1} {:14.5} {:14.5} {:8.4} {:14.5} {:14.5} {:8.4}",
            a.stopping_ev_per_m * 1e-9,
            b.stopping_ev_per_m * 1e-9,
            b.stopping_ev_per_m / a.stopping_ev_per_m,
            a.imfp_m() * 1e9,
            b.imfp_m() * 1e9,
            b.imfp_m() / a.imfp_m(),
        );
    }
    Ok(())
}
