//! A steppable simulation session: the full `lindhard` engine behind a handle
//! a Web Worker can drive in small chunks.
//!
//! [`Simulation`] is plain Rust, so the native tests exercise the code the
//! browser runs; the `#[wasm_bindgen]` wrapper in `lib.rs` only marshals
//! arguments. A session resolves one TOML run description (`docs/cli.md`) and
//! builds the engine with [`Resolved::bca`], so every engine option the input
//! configures (the beam divergence, the primary-energy spectrum) is applied
//! exactly as the CLI's driver applies it; then [`Simulation::step`] runs the
//! next histories with [`Bca::run_range`](lindhard::ion::bca::Bca::run_range)
//! and merges them into the accumulated tally, so the result does not depend
//! on how the run was cut into steps.
//!
//! # Chunk invariance
//!
//! History `i` always uses the random stream `stream(seed, i)`, and the
//! recorder keeps histories by their global index, so which histories are
//! simulated and recorded never depends on the step sizes. The floating-point
//! moments of the tally are combined chunk by chunk in index order; a tally
//! cut at a different place would be combined in a different order and differ
//! in the last bits. [`Simulation::step`] therefore always runs whole engine
//! chunks (`chunk_size` histories, 64 by default; the last chunk of the run
//! may be shorter) and merges them one at a time, so the report and the
//! trajectories are byte-identical for any sequence of step sizes. The price
//! is that a step rounds up to the next chunk boundary.
//!
//! # Not supported
//!
//! Crystal targets (`[[crystal]]`) and dynamic composition (`[dynamic]`) are
//! rejected with an error: the session transports an amorphous, static
//! target.

use lindhard::input::{Input, Resolved};
use lindhard::ion::bca::BcaTally;
use lindhard::ion::scattering::ScatteringTable;
use lindhard::ion::stopping::ElectronicStopping;
use lindhard::tally::{
    History, IonTally, Recorded, TrackEnd, TrajectoryConfig, TrajectoryRecorder,
};
use serde_json::{json, Value};

/// Largest run a session accepts (`run.ions` of the input).
pub const MAX_IONS: u32 = 1_000_000;

/// A run in progress. See the module documentation.
pub struct Simulation {
    resolved: Resolved,
    table: ScatteringTable,
    stopping: Box<dyn ElectronicStopping + Send + Sync>,
    /// An empty tally with the recorder's caps: what each engine chunk starts
    /// from.
    fresh: Recorded<IonTally>,
    /// Everything run so far.
    acc: Recorded<IonTally>,
    /// Histories run so far.
    done: u64,
    /// `run.ions`.
    total: u64,
}

impl std::fmt::Debug for Simulation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Simulation")
            .field("done", &self.done)
            .field("total", &self.total)
            .finish_non_exhaustive()
    }
}

impl Simulation {
    /// Resolve `input_toml` and prepare a run.
    ///
    /// The first `record_histories` histories (by index) are recorded as
    /// trajectories of at most `max_points` points each. Every error of the
    /// input (bad TOML, unknown element, energy out of range, an unsupported
    /// section) is returned as text.
    pub fn new(input_toml: &str, record_histories: u32, max_points: u32) -> Result<Self, String> {
        let input = Input::from_toml_str(input_toml).map_err(|e| e.to_string())?;
        if input.dynamic.is_some() {
            return Err("[dynamic] composition is not supported by the session".into());
        }
        if !input.crystal.is_empty() {
            return Err("[[crystal]] targets are not supported by the session".into());
        }
        let resolved = input.resolve().map_err(|e| e.to_string())?;
        let total = resolved.input.run.ions;
        if total == 0 || total > u64::from(MAX_IONS) {
            return Err(format!("run.ions must be between 1 and {MAX_IONS}"));
        }

        let table = resolved.scattering_table();
        let stopping = resolved.stopping_model();

        // A trial engine validates the beam, the configuration, the
        // materials and every engine option of the input now, so `step`
        // cannot fail on the input. `Bca` borrows the target, the stopping
        // model and the table, so each step rebuilds it with the same call
        // rather than storing it next to the data it borrows.
        let bca = resolved
            .bca(&*stopping, &table)
            .map_err(|e| e.to_string())?;
        let proto = IonTally::new(
            &resolved.stack,
            &bca.species_z(),
            resolved.ion_tally_config().map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        let acc = Recorded {
            inner: proto,
            recorder: TrajectoryRecorder::new(TrajectoryConfig {
                max_histories: u64::from(record_histories),
                max_points: max_points as usize,
            }),
        };
        Ok(Self {
            resolved,
            table,
            stopping,
            fresh: acc.clone(),
            acc,
            done: 0,
            total,
        })
    }

    /// Run at least the next `count` histories, rounded up to a whole
    /// engine chunk (or to the end of the run), and return the histories done
    /// so far. `count = 0` runs nothing. A transport failure is returned as
    /// text and leaves the session as it was.
    pub fn step(&mut self, count: u32) -> Result<u32, String> {
        if count == 0 || self.done() {
            return Ok(self.done_u32());
        }
        let chunk = self.resolved.config.chunk_size.max(1);
        let want = self.done + u64::from(count);
        let end = want.div_ceil(chunk).saturating_mul(chunk).min(self.total);
        let bca = self
            .resolved
            .bca(&*self.stopping, &self.table)
            .map_err(|e| e.to_string())?;
        // One engine chunk per call, merged in index order. All chunks are
        // run before any is merged, so a failure changes nothing.
        let mut parts = Vec::new();
        let mut next = self.done;
        while next < end {
            let stop = (next + chunk).min(end);
            let part = bca
                .run_range(next, stop - next, || self.fresh.clone())
                .map_err(|e| e.to_string())?;
            parts.push(part);
            next = stop;
        }
        for part in parts {
            self.acc.merge(part);
        }
        self.done = end;
        Ok(self.done_u32())
    }

    fn done_u32(&self) -> u32 {
        u32::try_from(self.done).unwrap_or(u32::MAX)
    }

    /// True once `run.ions` histories have run.
    pub fn done(&self) -> bool {
        self.done >= self.total
    }

    /// Histories run so far.
    pub fn histories_done(&self) -> u64 {
        self.done
    }

    /// The warnings of the resolver (for example a beam energy outside the
    /// validity range of the stopping model).
    pub fn warnings(&self) -> Vec<String> {
        self.resolved.warnings.clone()
    }

    /// The ion tally's report for everything run so far, as JSON. The dual
    /// Pearson fit is made if the input's `tally.dual_pearson` asks for it.
    pub fn report_json(&self) -> Result<String, String> {
        let report = self
            .acc
            .inner
            .report(self.resolved.input.tally.dual_pearson);
        serde_json::to_string(&report).map_err(|e| e.to_string())
    }

    /// The histories recorded since the last call, as JSON, and forget them.
    ///
    /// `{"histories": [{"index", "truncated", "tracks": [{"species",
    /// "generation", "parent", "end", "points": [[x, y, z, energy_ev], ...]}]}]}`
    /// with positions in metres. `end` is `"stopped"`, `"front"`, `"back"`,
    /// `"side"` or `null` (the point cap dropped it).
    pub fn take_trajectories_json(&mut self) -> String {
        let taken = std::mem::take(&mut self.acc.recorder.histories);
        json!({ "histories": taken.iter().map(history_json).collect::<Vec<_>>() }).to_string()
    }
}

fn history_json(h: &History) -> Value {
    let tracks: Vec<Value> = h
        .tracks
        .iter()
        .map(|t| {
            let end = match t.end {
                None => Value::Null,
                Some(TrackEnd::Stopped) => json!("stopped"),
                Some(TrackEnd::Escaped(f)) => json!(format!("{f:?}").to_lowercase()),
            };
            let points: Vec<Value> = t
                .points
                .iter()
                .map(|p| json!([p.pos[0], p.pos[1], p.pos[2], p.energy_ev]))
                .collect();
            json!({
                "species": t.species,
                "generation": t.generation,
                "parent": t.parent,
                "end": end,
                "points": points,
            })
        })
        .collect();
    json!({ "index": h.index, "truncated": h.truncated, "tracks": tracks })
}

/// The element table of `lindhard::elements` as JSON: an array of
/// `{"symbol", "z", "atomic_weight", "weight_is_mass_number",
/// "density_g_cm3", "default_displacement_ev", "default_surface_binding_ev"}`,
/// with `null` where the table records no value.
pub fn elements_json() -> String {
    let rows: Vec<Value> = lindhard::elements::ELEMENTS
        .iter()
        .map(|e| {
            json!({
                "symbol": e.symbol,
                "z": e.z,
                "atomic_weight": e.atomic_weight,
                "weight_is_mass_number": e.weight_is_mass_number,
                "density_g_cm3": e.density_g_cm3,
                "default_displacement_ev": e.default_displacement_ev,
                "default_surface_binding_ev": e.default_surface_binding_ev,
            })
        })
        .collect();
    Value::Array(rows).to_string()
}
