//! WebAssembly front end for `lindhard`: the steppable [`Simulation`] session
//! behind the static simulator page in `web/`, plus the element table
//! ([`elements_json`]) and the library version.
//!
//! [`Simulation`] is plain Rust (it resolves an ordinary `lindhard` input and
//! drives the library's ion transport one engine chunk at a time), so the
//! native tests exercise exactly the code the browser runs. The
//! `#[wasm_bindgen]` items below only marshal arguments and results.
//!
//! The library itself stays front-end agnostic: this crate depends on
//! `lindhard` and `wasm-bindgen`, and nothing else gains a WASM dependency.
//! It does not use `lindhard-cli`, whose run driver times itself with
//! `std::time::Instant`, which is unavailable on `wasm32-unknown-unknown`.

mod session;

pub use session::{elements_json, Simulation, MAX_IONS};

use wasm_bindgen::prelude::*;

/// The `lindhard` version this module was built from.
#[wasm_bindgen]
pub fn version() -> String {
    lindhard::VERSION.to_owned()
}

/// A steppable simulation session for the in-browser simulator: the full
/// engine behind a handle that a Web Worker drives in small chunks. The
/// methods marshal to [`Simulation`], where the behaviour is documented.
#[wasm_bindgen(js_name = Simulation)]
#[derive(Debug)]
pub struct WasmSimulation(Simulation);

#[wasm_bindgen(js_class = Simulation)]
impl WasmSimulation {
    /// Resolve the TOML run description `input_toml` and prepare the run.
    /// Throws a JavaScript `Error` with the library's message on invalid
    /// input.
    #[wasm_bindgen(constructor)]
    pub fn new(
        input_toml: &str,
        record_histories: u32,
        max_points: u32,
    ) -> Result<WasmSimulation, JsError> {
        Simulation::new(input_toml, record_histories, max_points)
            .map(WasmSimulation)
            .map_err(|e| JsError::new(&e))
    }

    /// Run the next `count` histories; returns the histories done so far.
    pub fn step(&mut self, count: u32) -> Result<u32, JsError> {
        self.0.step(count).map_err(|e| JsError::new(&e))
    }

    /// True once all `run.ions` histories have run.
    pub fn done(&self) -> bool {
        self.0.done()
    }

    /// The ion tally's report of everything run so far, as JSON.
    pub fn report_json(&self) -> Result<String, JsError> {
        self.0.report_json().map_err(|e| JsError::new(&e))
    }

    /// The trajectories recorded since the last call, as JSON; forgets them.
    pub fn take_trajectories_json(&mut self) -> String {
        self.0.take_trajectories_json()
    }

    /// The resolver's warnings.
    pub fn warnings(&self) -> Vec<String> {
        self.0.warnings()
    }
}

/// The element table of the library as JSON (see [`elements_json`]).
#[wasm_bindgen(js_name = elements_json)]
pub fn elements_json_js() -> String {
    session::elements_json()
}
