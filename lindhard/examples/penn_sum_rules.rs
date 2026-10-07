//! Print the sum-rule report of an optical ELF, and the single-pole Penn IMFP
//! and stopping power built from it.
//!
//! ```text
//! cargo run --release -p lindhard --example penn_sum_rules
//! cargo run --release -p lindhard --example penn_sum_rules -- my_elf.toml [target_density_per_m3]
//! ```
//!
//! With no argument the ELF is a synthetic Drude plasmon (E_p = 20 eV,
//! γ = 5 eV, not physical data; no optical data is shipped with lindhard). A
//! path reads an ELF in the `OpticalElf` TOML form (`material`, `provenance`,
//! `energy_ev`, `elf`); a second argument gives the density of atoms or
//! molecules (m⁻³) so that the report includes N_eff per target unit.

use lindhard::electron::data::OpticalElf;
use lindhard::electron::inelastic::{
    DrudeLorentz, DrudeLorentzOscillator, SinglePolePenn, SumRuleReport,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (elf, report) = match args.first() {
        Some(path) => {
            let elf = OpticalElf::from_toml_file(path)?;
            let report = match args.get(1) {
                Some(n) => SumRuleReport::with_target_density(&elf, n.parse()?)?,
                None => SumRuleReport::new(&elf),
            };
            (elf, report)
        }
        None => {
            let model = DrudeLorentz::new(vec![DrudeLorentzOscillator::plasmon(20.0, 5.0)])?;
            let elf = model.to_optical_elf("synthetic Drude plasmon", 1e-3, 1e5, 4000)?;
            // Two electrons per (fictitious) target unit, so N_eff is shown.
            let report = SumRuleReport::with_target_density(
                &elf,
                SumRuleReport::new(&elf).electron_density_per_m3 / 2.0,
            )?;
            println!(
                "Closed forms of the fixture: f-sum {:.6e} eV^2, P_eff {:.6}, E_p {:.6} eV\n",
                model.f_sum_ev2(),
                model.p_eff(),
                model.plasma_energy_ev()
            );
            (elf, report)
        }
    };
    println!("{report}");
    if let Some(n) = report.target_density_per_m3 {
        println!("N_eff(W), electrons per target unit (density {n:.4e} m^-3):");
        for w in [10.0, 30.0, 100.0, 300.0, 1e3, 1e4] {
            if let Some(x) = report.effective_electrons_up_to(w) {
                println!("  W = {w:>8} eV  N_eff = {x:.5}");
            }
        }
        println!();
    }

    let penn = SinglePolePenn::new(elf);
    println!("Single-pole Penn model (nonrelativistic kinematics, E_F = 0):");
    println!(
        "  {:>10}  {:>12}  {:>14}  {:>14}",
        "E (eV)", "IMFP (nm)", "S (eV/nm)", "Bethe S (eV/nm)"
    );
    for e in [20.0, 50.0, 100.0, 300.0, 1e3, 1e4, 5e4] {
        let p = penn.imfp_and_stopping(e)?;
        let bethe = penn.bethe_stopping_ev_per_m(e)?;
        println!(
            "  {:>10}  {:>12.4}  {:>14.4}  {:>14.4}",
            e,
            p.imfp_m() * 1e9,
            p.stopping_ev_per_m * 1e-9,
            bethe * 1e-9
        );
    }
    Ok(())
}
