//! The inelastic table builders (parent module docs): the options, the
//! single-pole, any-model and shell-channel builders, and the rows they
//! share.

#[cfg(doc)]
use super::super::full_penn::{DiimfpGrid, FullPenn};
use super::super::inner_shell::ShellResolvedChannels;
use super::super::model::PennInelastic;
use super::super::penn::SinglePolePenn;
use super::density::{loss_density, loss_quantile, refine_probability, LinearDensity, LossModel};
use super::refine::{refinement_provenance, row_rate, with_refined_grid, RateRefinement};
use super::{
    check_no_exchange, invalid, EnergyAxis, InelasticTableError, DEFAULT_DENSITY_TOLERANCE,
    DEFAULT_LOGISTIC_STEP, DEFAULT_PROBABILITY_TAIL, DEFAULT_REFINE_TOLERANCE,
};
use crate::electron::data::{
    CrossSectionTable, CrossSectionTableParts, SamplingAxis, ShellChannelTable,
};
use crate::electron::elastic::table::{logistic_probability_grid, material_identity};
use crate::material::Material;
use rayon::prelude::*;

/// Grids and tolerances of a table build.
#[derive(Debug, Clone, PartialEq)]
pub struct InelasticTableOptions {
    /// Incident kinetic energies on the axis of [`Self::axis`], eV (strictly
    /// increasing, >= 2 points): above the model's Fermi level by default,
    /// above the band bottom with [`EnergyAxis::BandBottom`].
    pub energy_ev: Vec<f64>,
    /// The axis `energy_ev` is on ([`EnergyAxis::ModelFermiLevel`] by
    /// default).
    pub axis: EnergyAxis,
    /// Starting cumulative-probability grid of the stored inverse CDFs (from
    /// exactly 0 to exactly 1); the final grid when `refine_tolerance` is
    /// `None`.
    pub probability: Vec<f64>,
    /// Tolerance of the loss-density model ([`DEFAULT_DENSITY_TOLERANCE`]).
    pub density_tolerance: f64,
    /// Target of the probability-grid refinement, relative to the mean loss
    /// of each row (`None`: no refinement).
    pub refine_tolerance: Option<f64>,
    /// Extra rows that resolve the inelastic rate just above the Fermi level
    /// (module docs, "Rate refinement"; #339). `None` (the default): the
    /// table has exactly the rows of `energy_ev`. Only valid with
    /// [`EnergyAxis::BandBottom`]; refused on the model's own axis.
    pub rate_refinement: Option<RateRefinement>,
}

impl InelasticTableOptions {
    /// The defaults on the given energy grid.
    pub fn new(energy_ev: Vec<f64>) -> Self {
        Self {
            energy_ev,
            axis: EnergyAxis::ModelFermiLevel,
            probability: default_probability_grid(),
            density_tolerance: DEFAULT_DENSITY_TOLERANCE,
            refine_tolerance: Some(DEFAULT_REFINE_TOLERANCE),
            rate_refinement: None,
        }
    }

    /// The same options with the table stored on `axis`.
    pub fn with_axis(mut self, axis: EnergyAxis) -> Self {
        self.axis = axis;
        self
    }

    /// The same options with the rate refinement `refinement` (module docs,
    /// "Rate refinement"); needs [`EnergyAxis::BandBottom`].
    pub fn with_rate_refinement(mut self, refinement: RateRefinement) -> Self {
        self.rate_refinement = Some(refinement);
        self
    }
}

/// The default starting probability grid: [`logistic_probability_grid`] with
/// [`DEFAULT_PROBABILITY_TAIL`] and [`DEFAULT_LOGISTIC_STEP`].
pub fn default_probability_grid() -> Vec<f64> {
    logistic_probability_grid(DEFAULT_PROBABILITY_TAIL, DEFAULT_LOGISTIC_STEP)
        .expect("valid default grid")
}

fn check_options(o: &InelasticTableOptions) -> Result<(), InelasticTableError> {
    let p = &o.probability;
    let ok = p.len() >= 2
        && p[0] == 0.0
        && p[p.len() - 1] == 1.0
        && p.iter().all(|x| x.is_finite())
        && p.windows(2).all(|w| w[1] > w[0]);
    if !ok {
        return Err(invalid(
            "probability grid",
            "need >= 2 finite, strictly increasing points from exactly 0 to exactly 1",
        ));
    }
    if !(o.density_tolerance.is_finite() && (1.0e-9..=1.0e-2).contains(&o.density_tolerance)) {
        return Err(invalid(
            "density tolerance",
            format!("{} (need 1e-9 ..= 1e-2)", o.density_tolerance),
        ));
    }
    if let Some(t) = o.refine_tolerance {
        if !(t.is_finite() && t > 0.0) {
            return Err(invalid(
                "refine tolerance",
                format!("{t} (need finite, > 0)"),
            ));
        }
    }
    if let Some(r) = &o.rate_refinement {
        if o.axis != EnergyAxis::BandBottom {
            return Err(invalid(
                "rate refinement",
                "only the band-bottom energy axis is refined (EnergyAxis::BandBottom)",
            ));
        }
        let (lo, hi) = r.window_ev;
        if !(lo.is_finite() && hi.is_finite() && lo > 0.0 && hi > lo) {
            return Err(invalid(
                "rate refinement window",
                format!("{lo}..{hi} eV above E_F (need finite, 0 < low < high)"),
            ));
        }
        if !(r.tolerance.is_finite() && r.tolerance > 0.0) {
            return Err(invalid(
                "rate refinement tolerance",
                format!("{} (need finite, > 0)", r.tolerance),
            ));
        }
    }
    Ok(())
}

/// Build the inelastic energy-loss table of `penn` on the grid of `options`.
///
/// `material` only labels the table (its identity string); the density that
/// sets the absolute rate is inside the optical ELF of `penn`
/// ([`SinglePolePenn::imfp_and_stopping`] returns `λ⁻¹` in m⁻¹ with no
/// further input). Energies are on the axis of `options.axis`: kinetic
/// energies above the model's Fermi level by default, above the band bottom
/// with [`EnergyAxis::BandBottom`] (module docs, "Energy axis").
pub fn build_inelastic_table(
    penn: &SinglePolePenn,
    material: &Material,
    options: &InelasticTableOptions,
) -> Result<CrossSectionTable, InelasticTableError> {
    check_no_exchange(penn)?;
    check_options(options)?;
    let (options, added) = with_refined_grid(penn, options)?;
    let options = options.as_ref();
    let (inverse_mfp_per_m, probability, quantiles) = loss_rows(penn, options)?;
    let energy = &options.energy_ev;
    let elf = penn.optical_elf();
    let model = format!(
        "lindhard {} electron::inelastic::table: single-pole Penn DIIMFP (nonrelativistic), \
         W inverse CDF from an adaptive piecewise-linear density in ln W",
        env!("CARGO_PKG_VERSION")
    );
    let provenance = format!(
        "computed by lindhard {} from published formulas (see electron::inelastic docs); \
         optical ELF: {} ({}); Fermi energy {} eV; integration tolerance {}; \
         density tolerance {}; {} energies {}..{} eV; {} probability points",
        env!("CARGO_PKG_VERSION"),
        elf.material(),
        elf.provenance(),
        penn.fermi_energy_ev(),
        penn.relative_tolerance(),
        options.density_tolerance,
        energy.len(),
        energy[0],
        energy[energy.len() - 1],
        probability.len()
    ) + &axis_provenance(penn, options)
        + &refinement_provenance(options, added);
    Ok(CrossSectionTable::new(CrossSectionTableParts {
        model,
        material: material_identity(material),
        provenance,
        axis: SamplingAxis::InelasticEnergyLoss,
        energy_ev: energy.clone(),
        inverse_mfp_per_m,
        probability,
        quantiles,
    })?)
}

/// Build the inelastic energy-loss table of any of the three ELF-extension
/// models ([`PennInelastic`]: single-pole Penn, full Penn or Mermin-ELF) on
/// the grid of `options`, by the same procedure as [`build_inelastic_table`]
/// (module docs), with that model's DIIMFP and inverse IMFP in place of the
/// single pole's. For [`PennInelastic::SinglePole`] the rows are those of
/// [`build_inelastic_table`]; only the `model` and `provenance` strings
/// differ, and they carry [`PennInelastic::model_identity`], so a table
/// names the algorithm it was built with.
///
/// The full Penn and Mermin models evaluate their DIIMFP by numerical
/// integration, so a table of either costs far more than a single-pole table
/// on the same grid. For the full model the rows read the DIIMFP from a
/// [`DiimfpGrid`] built once for the grid energies
/// ([`FullPenn::diimfp_grid`]: interpolation error held to the model
/// tolerance), which is what makes such a table affordable (#256); the
/// inverse mean free paths are still the model's own.
pub fn build_inelastic_table_for_model(
    model: &PennInelastic,
    material: &Material,
    options: &InelasticTableOptions,
) -> Result<CrossSectionTable, InelasticTableError> {
    if let PennInelastic::SinglePole(p) = model {
        check_no_exchange(p)?;
    }
    check_options(options)?;
    let (options, added) = with_refined_grid(model, options)?;
    let options = options.as_ref();
    let (inverse_mfp_per_m, probability, quantiles) = loss_rows(model, options)?;
    let energy = &options.energy_ev;
    let elf = model.optical_elf();
    let identity = model.model_identity();
    // the full model's rows read the DIIMFP from its loss grid (#256)
    let via = match model {
        PennInelastic::Full(_) => " (from the T-independent loss grid FullPenn::diimfp_grid)",
        _ => "",
    };
    let model_text = format!(
        "lindhard {} electron::inelastic::table: {identity} DIIMFP{via}, \
         W inverse CDF from an adaptive piecewise-linear density in ln W",
        env!("CARGO_PKG_VERSION")
    );
    let provenance = format!(
        "computed by lindhard {} from published formulas (see electron::inelastic docs); \
         model: {identity}; optical ELF: {} ({}); density tolerance {}; \
         {} energies {}..{} eV; {} probability points",
        env!("CARGO_PKG_VERSION"),
        elf.material(),
        elf.provenance(),
        options.density_tolerance,
        energy.len(),
        energy[0],
        energy[energy.len() - 1],
        probability.len()
    ) + &axis_provenance(model, options)
        + &refinement_provenance(options, added);
    Ok(CrossSectionTable::new(CrossSectionTableParts {
        model: model_text,
        material: material_identity(material),
        provenance,
        axis: SamplingAxis::InelasticEnergyLoss,
        energy_ev: energy.clone(),
        inverse_mfp_per_m,
        probability,
        quantiles,
    })?)
}

/// The model's kinetic energy `T` of each row of `options`: the axis energy
/// itself on [`EnergyAxis::ModelFermiLevel`], `E - E_F` (`E_F` the model's
/// Fermi energy) on [`EnergyAxis::BandBottom`].
fn model_energies_ev<M: LossModel>(penn: &M, options: &InelasticTableOptions) -> Vec<f64> {
    match options.axis {
        EnergyAxis::ModelFermiLevel => options.energy_ev.clone(),
        EnergyAxis::BandBottom => {
            let ef = penn.fermi_energy_ev();
            options.energy_ev.iter().map(|&e| e - ef).collect()
        }
    }
}

/// The provenance text of the axis: empty on the model's own axis (so tables
/// built before the axis option existed keep their bytes), the axis and the
/// Fermi energy otherwise.
fn axis_provenance<M: LossModel>(penn: &M, options: &InelasticTableOptions) -> String {
    match options.axis {
        EnergyAxis::ModelFermiLevel => String::new(),
        EnergyAxis::BandBottom => format!(
            "; energy axis: band bottom (row at E is the model's at T = E - {} eV)",
            penn.fermi_energy_ev()
        ),
    }
}

/// The tables of a [`ShellResolvedChannels`]: the valence channel's and one
/// per inner shell, in the order of
/// [`ShellResolvedChannels::inner_shells`]
/// ([`build_shell_channel_tables`]).
#[derive(Debug, Clone, PartialEq)]
pub struct ShellResolvedTables {
    /// The valence channel, as [`build_inelastic_table`] of
    /// [`ShellResolvedChannels::valence`].
    pub valence: CrossSectionTable,
    /// The inner-shell channels.
    pub shells: Vec<ShellChannelTable>,
}

/// Build the energy-loss tables of every channel of `channels` on the grid of
/// `options`: the valence table, and for inner shell `i` the table of
/// [`ShellResolvedChannels::shell_model`] (that shell's own optical ELF with
/// the valence model's settings), wrapped with the shell's identity, binding
/// energy and `binding_provenance` (the origin of the binding energies, for
/// example [`crate::electron::data::SubshellBindingTable::provenance`]).
///
/// Each table is built by [`build_inelastic_table`], one after the other in
/// channel order, so the result is bit-identical at any thread count. Its
/// rate is the channel's own `λ⁻¹_j(E)` and its rows the inverse CDF of that
/// channel's DIIMFP: the **channel-first** construction, in which a
/// transport picks a channel with probability `λ⁻¹_j / Σ λ⁻¹` and then draws
/// the loss from that channel's table (not the conditional channel choice of
/// [`ShellResolvedChannels::sample_channel`], which belongs with a loss drawn
/// from the total DIIMFP; the two are not mixed).
///
/// **Supported model.** Single-pole Penn per channel only, without the
/// exchange correction (rejected, as by [`build_inelastic_table`]); so the
/// binding energy does not enter a shell's DIIMFP, which is zero below the
/// first energy of the shell's ELF (at or above `B` by construction of
/// [`ShellResolvedChannels::new`]). Each shell row is floored at `B`, which
/// only removes the rounding of `exp(ln W)` at the edge. No full Penn or
/// Mermin valence model is wired to shell channels, and the valence ELF must
/// not already contain the shells (no partition rule is applied, see
/// [`super::inner_shell`](super::super::inner_shell)).
///
/// **Energy axis.** As for [`build_inelastic_table`], every table (valence
/// and shells) is on the axis of `options.axis`: kinetic energies above the
/// model's Fermi energy (zero by default) with
/// [`EnergyAxis::ModelFermiLevel`], above the band bottom with
/// [`EnergyAxis::BandBottom`], each row then the model's at `T = E - E_F`
/// with `E_F` the Fermi energy of the model (the shell models share the
/// valence model's settings, so its Fermi energy) (module docs, "Energy
/// axis"). See `crate::electron::transport`, "Energy reference of the
/// inelastic table", for how the transport reads them.
pub fn build_shell_channel_tables(
    channels: &ShellResolvedChannels,
    material: &Material,
    options: &InelasticTableOptions,
    binding_provenance: &str,
) -> Result<ShellResolvedTables, InelasticTableError> {
    let valence = build_inelastic_table(channels.valence(), material, options)?;
    let mut shells = Vec::with_capacity(channels.inner_shells().len());
    for (i, s) in channels.inner_shells().iter().enumerate() {
        let b = s.binding_energy_ev;
        let base = build_inelastic_table(channels.shell_model(i), material, options)?;
        let mut parts = base.parts().clone();
        let label = s.subshell.label();
        parts.model = format!("{} [inner-shell channel Z = {} {label}]", parts.model, s.z);
        parts.provenance = format!(
            "{}; inner-shell channel Z = {} {label}, binding energy {b} eV ({}), \
             losses floored at the binding energy",
            parts.provenance,
            s.z,
            binding_provenance.trim()
        );
        for row in &mut parts.quantiles {
            for w in row.iter_mut() {
                *w = w.max(b);
            }
        }
        let table = CrossSectionTable::new(parts)?;
        shells.push(ShellChannelTable::new(
            s.z,
            s.subshell,
            b,
            binding_provenance,
            table,
        )?);
    }
    Ok(ShellResolvedTables { valence, shells })
}

/// Inverse mean free paths, the (refined) probability grid and the quantile
/// rows of a table of `penn` on the grid of `options` (already checked).
type LossRows = (Vec<f64>, Vec<f64>, Vec<Vec<f64>>);

pub(super) fn loss_rows<M: LossModel>(
    penn: &M,
    options: &InelasticTableOptions,
) -> Result<LossRows, InelasticTableError> {
    // The model's energy `T` of each row (module docs, "Energy axis").
    let energy = &model_energies_ev(penn, options);
    // The DIIMFP grid covers only the rows that get losses: on the
    // band-bottom axis a row at or below the Fermi level has none (below),
    // and the grid would reject its energy.
    let grid = match options.axis {
        EnergyAxis::ModelFermiLevel => penn.diimfp_grid(energy)?,
        EnergyAxis::BandBottom => {
            let lossy: Vec<f64> = energy.iter().copied().filter(|&e| e > 0.0).collect();
            penn.diimfp_grid(&lossy)?
        }
    };
    // Rows in parallel, collected in grid order; the lowest failing energy
    // is the reported error whatever the thread count.
    type Row = Result<(f64, Option<LinearDensity>), InelasticTableError>;
    let results: Vec<Row> = energy
        .par_iter()
        .map(|&e| {
            let inv = row_rate(penn, options.axis, e)?;
            if inv > 0.0 {
                let d = loss_density(penn, grid.as_ref(), e, inv, options.density_tolerance)?;
                Ok((inv, Some(d)))
            } else {
                Ok((0.0, None))
            }
        })
        .collect();
    let mut inverse_mfp_per_m = Vec::with_capacity(energy.len());
    let mut rows: Vec<Option<LinearDensity>> = Vec::with_capacity(energy.len());
    for r in results {
        let (inv, d) = r?;
        inverse_mfp_per_m.push(inv);
        rows.push(d);
    }
    let probability = &options.probability;
    let quantiles: Vec<Vec<f64>> = rows
        .par_iter()
        .zip(energy.par_iter())
        .map(|(r, &e)| match r {
            Some(d) => probability
                .iter()
                .map(|&u| loss_quantile(d, e, u))
                .collect(),
            None => Vec::new(),
        })
        .collect();
    let (probability, quantiles) = match options.refine_tolerance {
        Some(tol) => refine_probability(&rows, energy, probability.clone(), quantiles, tol),
        None => (probability.clone(), quantiles),
    };
    Ok((inverse_mfp_per_m, probability, quantiles))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rate_refinement_needs_the_band_bottom_axis_and_a_valid_window() {
        let mut o = InelasticTableOptions::new(vec![5.0, 20.0])
            .with_rate_refinement(RateRefinement::default());
        assert!(check_options(&o)
            .unwrap_err()
            .to_string()
            .contains("band-bottom"));
        o.axis = EnergyAxis::BandBottom;
        assert!(check_options(&o).is_ok());
        for (window, tol) in [((3.0, 2.0), 0.01), ((0.0, 2.0), 0.01), ((2.0, 20.0), 0.0)] {
            let r = o.rate_refinement.as_mut().unwrap();
            r.window_ev = window;
            r.tolerance = tol;
            assert!(check_options(&o).is_err(), "{window:?} {tol}");
        }
    }
}
