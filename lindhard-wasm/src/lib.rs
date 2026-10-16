//! WebAssembly front end for `lindhard`: the range calculator behind the
//! static page in `web/`.
//!
//! [`range_profile`] is plain Rust (it builds an ordinary `lindhard` input,
//! runs the library's ion transport and reduces the stopped depths to `Rp`,
//! `ΔRp` and a depth histogram), so the native tests exercise exactly the code
//! the browser runs. The `#[wasm_bindgen]` items at the bottom only marshal
//! arguments and results.
//!
//! The library itself stays front-end agnostic: this crate depends on
//! `lindhard` and `wasm-bindgen`, and nothing else gains a WASM dependency.
//! It does not use `lindhard-cli`, whose run driver times itself with
//! `std::time::Instant`, which is unavailable on `wasm32-unknown-unknown`.

use lindhard::input::Input;
use lindhard::ion::bca::{Bca, BcaTally, Face, Particle};
use lindhard::ion::potential::Potential;
use lindhard::ion::scattering::ScatteringTable;
use lindhard::tally::Moments;
use wasm_bindgen::prelude::*;

/// Fallback displacement energy, eV, for a target element whose `E_d` the
/// element table does not record (the same illustrative 15 eV as
/// `examples/b_5keV_si.toml`). It only sets the damage bookkeeping; the
/// projected range does not depend on it for a light ion at low energy, and
/// the page says so.
const FALLBACK_E_D_EV: f64 = 15.0;

/// Number of bins of the returned depth histogram.
pub const PROFILE_BINS: usize = 100;

/// Largest run the calculator accepts (the page is a calculator, not a batch
/// runner).
pub const MAX_IONS: u32 = 1_000_000;

/// What to compute.
#[derive(Debug, Clone, PartialEq)]
pub struct Request {
    /// Beam element symbol, for example `"B"`.
    pub ion: String,
    /// Beam kinetic energy, eV.
    pub energy_ev: f64,
    /// Target element symbol (a pure element at its tabulated density).
    pub target: String,
    /// Number of ion histories.
    pub ions: u32,
    /// Seed of the per-ion random streams.
    pub seed: u32,
}

/// The result: the depth moments of the ions that stopped in the target and
/// their depth histogram.
#[derive(Debug, Clone, PartialEq)]
pub struct Profile {
    /// Projected range `Rp` (mean stopping depth), nm.
    pub rp_nm: f64,
    /// Straggle `ΔRp` (standard deviation of the stopping depth), nm.
    pub drp_nm: f64,
    /// Standard error of `rp_nm`, nm.
    pub rp_std_err_nm: f64,
    /// Standard error of `drp_nm`, nm.
    pub drp_std_err_nm: f64,
    /// Ions that came to rest in the target.
    pub stopped: u64,
    /// Ions run.
    pub histories: u64,
    /// Width of one histogram bin, nm.
    pub bin_nm: f64,
    /// Stopped ions per bin; bin `i` covers `[i, i + 1) * bin_nm`.
    pub counts: Vec<u64>,
    /// The library's non-fatal advice for this input, for example a beam
    /// energy outside the validity range of the stopping model.
    pub warnings: Vec<String>,
}

/// Records the stopping depth of every primary. `merge` appends in chunk
/// order, so the list (and everything computed from it) is independent of
/// scheduling.
#[derive(Debug, Default)]
struct DepthTally {
    histories: u64,
    depths_m: Vec<f64>,
}

impl BcaTally for DepthTally {
    fn begin_history(&mut self, _index: u64) {
        self.histories += 1;
    }

    fn stopped(&mut self, p: &Particle) {
        if p.is_primary() {
            self.depths_m.push(p.pos[0]);
        }
    }

    fn escaped(&mut self, _p: &Particle, _face: Face) {}

    fn merge(&mut self, other: Self) {
        self.histories += other.histories;
        self.depths_m.extend(other.depths_m);
    }
}

/// The TOML run description for `req`: amorphous target, ZBL screening,
/// Lindhard-Scharff electronic stopping, the choices of `examples/b_5keV_si.toml`.
fn input_toml(req: &Request) -> Result<String, String> {
    let symbol_ok =
        |s: &str| !s.is_empty() && s.len() <= 2 && s.chars().all(|c| c.is_ascii_alphabetic());
    if !symbol_ok(&req.ion) || !symbol_ok(&req.target) {
        return Err("species must be element symbols such as \"B\" or \"Si\"".into());
    }
    let target = lindhard::elements::element_by_symbol(&req.target)
        .ok_or_else(|| format!("unknown target element {:?}", req.target))?;
    let e_d = if target.default_displacement_ev.is_some() {
        String::new()
    } else {
        format!(
            "[physics.energies.{}]\ne_d_ev = {FALLBACK_E_D_EV:?}\n",
            req.target
        )
    };
    Ok(format!(
        "[beam]\nion = \"{ion}\"\nenergy_ev = {energy:?}\ntilt_deg = 7.0\n\n\
         [target]\nsubstrate = \"{target}\"\n\n\
         [physics]\npotential = \"zbl\"\nstopping = \"lindhard-scharff\"\n\
         free_path = \"constant\"\nprimary_cutoff_ev = 5.0\nrecoil_cutoff_ev = 2.0\n\n\
         {e_d}\n[run]\nions = {ions}\nseed = {seed}\n",
        ion = req.ion,
        energy = req.energy_ev,
        target = req.target,
        ions = req.ions,
        seed = req.seed,
    ))
}

/// Run `req` on the calling thread and reduce it to a [`Profile`].
///
/// Errors are the library's input errors (unknown element, energy out of
/// range, missing density) as text.
pub fn range_profile(req: &Request) -> Result<Profile, String> {
    if req.ions == 0 || req.ions > MAX_IONS {
        return Err(format!("ions must be between 1 and {MAX_IONS}"));
    }
    let input = Input::from_toml_str(&input_toml(req)?).map_err(|e| e.to_string())?;
    let r = input.resolve().map_err(|e| e.to_string())?;

    // The angle table depends on the screening function and length only, so
    // the Z pair is immaterial (the same construction as the CLI's driver).
    let z2 = r.layers[0].material.components()[0].z();
    let pot = Potential::new(r.screening, f64::from(r.beam.ion.z()), f64::from(z2))
        .with_length(r.screening_length);
    let table = ScatteringTable::build(&pot, &r.table_spec);
    let stopping = r.stopping_model();
    let bca =
        Bca::new(r.beam, &r.stack, r.config, &*stopping, &table).map_err(|e| e.to_string())?;
    let tally: DepthTally = bca.run(DepthTally::default).map_err(|e| e.to_string())?;

    let nm = 1e-9;
    let mut moments = Moments::new();
    for &x in &tally.depths_m {
        moments.push(x / nm);
    }
    let summary = moments
        .summary()
        .ok_or("fewer than two ions stopped in the target: no range to report")?;
    let deepest = tally.depths_m.iter().fold(0.0_f64, |a, &b| a.max(b)) / nm;
    let bin_nm = (deepest / PROFILE_BINS as f64).max(f64::MIN_POSITIVE) * (1.0 + 1e-9);
    let mut counts = vec![0u64; PROFILE_BINS];
    for &x in &tally.depths_m {
        let i = ((x / nm / bin_nm) as usize).min(PROFILE_BINS - 1);
        counts[i] += 1;
    }
    Ok(Profile {
        rp_nm: summary.mean,
        drp_nm: summary.std_dev,
        rp_std_err_nm: summary.mean_std_err,
        drp_std_err_nm: summary.std_dev_std_err,
        stopped: tally.depths_m.len() as u64,
        histories: tally.histories,
        bin_nm,
        counts,
        warnings: r.warnings,
    })
}

/// The result of [`calculate`], as JavaScript sees it.
#[wasm_bindgen]
#[derive(Debug)]
pub struct RangeResult(Profile);

#[wasm_bindgen]
impl RangeResult {
    /// Projected range `Rp`, nm.
    #[wasm_bindgen(getter)]
    pub fn rp_nm(&self) -> f64 {
        self.0.rp_nm
    }

    /// Straggle `ΔRp`, nm.
    #[wasm_bindgen(getter)]
    pub fn drp_nm(&self) -> f64 {
        self.0.drp_nm
    }

    /// Standard error of `Rp`, nm.
    #[wasm_bindgen(getter)]
    pub fn rp_std_err_nm(&self) -> f64 {
        self.0.rp_std_err_nm
    }

    /// Standard error of `ΔRp`, nm.
    #[wasm_bindgen(getter)]
    pub fn drp_std_err_nm(&self) -> f64 {
        self.0.drp_std_err_nm
    }

    /// Ions that stopped in the target.
    #[wasm_bindgen(getter)]
    pub fn stopped(&self) -> f64 {
        self.0.stopped as f64
    }

    /// Ions run.
    #[wasm_bindgen(getter)]
    pub fn histories(&self) -> f64 {
        self.0.histories as f64
    }

    /// Histogram bin width, nm.
    #[wasm_bindgen(getter)]
    pub fn bin_nm(&self) -> f64 {
        self.0.bin_nm
    }

    /// The library's warnings for this input, one per line (empty if none).
    #[wasm_bindgen(getter)]
    pub fn warnings(&self) -> String {
        self.0.warnings.join("\n")
    }

    /// Stopped ions per histogram bin.
    #[wasm_bindgen(getter)]
    pub fn counts(&self) -> Vec<f64> {
        self.0.counts.iter().map(|&c| c as f64).collect()
    }
}

/// Run the calculator: `ion` and `target` are element symbols, `energy_ev`
/// the beam energy in eV. Throws a JavaScript `Error` with the library's
/// message on invalid input.
#[wasm_bindgen]
pub fn calculate(
    ion: &str,
    energy_ev: f64,
    target: &str,
    ions: u32,
    seed: u32,
) -> Result<RangeResult, JsError> {
    let req = Request {
        ion: ion.to_owned(),
        energy_ev,
        target: target.to_owned(),
        ions,
        seed,
    };
    range_profile(&req)
        .map(RangeResult)
        .map_err(|e| JsError::new(&e))
}

/// The `lindhard` version this module was built from.
#[wasm_bindgen]
pub fn version() -> String {
    lindhard::VERSION.to_owned()
}
