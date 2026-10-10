//! Screen-oxide characterization (issue #280): how an amorphous SiO2 overlayer
//! changes what reaches a crystalline Si substrate and the channeling tail
//! there. A numerical characterization of the implemented model, not an
//! experimental validation; see `docs/screen-oxide.md`.
//!
//! The ordinary test is a small conservation and determinism case. The sweep
//! is `#[ignore]`d and writes machine-readable results:
//!
//! ```text
//! LINDHARD_SCREEN_OXIDE_OUT=/path/to/results.json \
//!   cargo test --release -p lindhard --test crystal_screen_oxide -- --ignored --nocapture
//! ```
//!
//! `LINDHARD_SCREEN_OXIDE_N` overrides the primaries per configuration.
//!
//! Coordinates: depths are `x` along the inward surface normal. `x = 0` is
//! the entrance surface (the front of the oxide, or of the bare crystal). The
//! substrate interface is at `x = T`, the oxide thickness (`T = 0` for bare
//! silicon). "Below the interface" is `x - T`.
//!
//! Entry into the substrate is the first time a primary reaches `x = T`
//! moving inward, taken from the electronic-loss hook of the segment ending
//! at the boundary (the energy and direction there). For `T = 0` it is the
//! incident beam. Histories that stop in the oxide, or return through the
//! front before reaching the substrate, are counted separately and are not in
//! any substrate statistic, so a loss of transmitted histories is not read as
//! dechanneling.

use std::sync::OnceLock;

use lindhard::geometry::Stack;
use lindhard::ion::bca::{
    Bca, BcaConfig, BcaTally, Beam, CrystalTarget, ElectronicChannel, EnergyBudget, Face, Particle,
    Thermal,
};
use lindhard::ion::crystal::debye::THETA_D_SI;
use lindhard::ion::crystal::{Lattice, Orientation};
use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use lindhard::ion::stopping::Ion;
use lindhard::material::Material;
use serde::Serialize;

const NM: f64 = 1e-9;
/// Tolerance for "at the interface", m. The boundary point is computed by the
/// geometry, so the match is to rounding error; a collision this close to the
/// interface without crossing it is not a realistic concern.
const AT_INTERFACE: f64 = 1e-14;
const SEED: u64 = 1;
const TEMPERATURE_K: f64 = 300.0;

fn table() -> &'static ScatteringTable {
    static T: OnceLock<ScatteringTable> = OnceLock::new();
    T.get_or_init(|| {
        ScatteringTable::build(
            &Potential::new(Screening::ZblUniversal, 14.0, 14.0),
            &TableSpec {
                eps_min: 1e-6,
                eps_max: 1e4,
                beta_min: 1e-5,
                beta_max: 1e2,
                per_decade: 16,
            },
        )
    })
}

/// Silicon at the crystal density with `E_d` = 15 eV (as `crystal_bca.rs`).
fn si() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0)], None).unwrap();
    m.set_displacement_energy_ev(14, 15.0).unwrap();
    m
}

/// Amorphous SiO2 at 2200 kg/m3, `E_d` = 15 eV for both elements.
fn sio2() -> Material {
    let mut m = Material::from_atom_fractions(&[(14, 1.0), (8, 2.0)], Some(2200.0)).unwrap();
    for z in [14, 8] {
        m.set_displacement_energy_ev(z, 15.0).unwrap();
        if m.surface_binding_energy_ev(z).is_err() {
            m.set_surface_binding_energy_ev(z, 2.0).unwrap();
        }
    }
    m
}

fn crystal_si(tilt_deg: f64, twist_deg: f64) -> CrystalTarget {
    let lat = Lattice::silicon();
    let o = Orientation::new(
        &lat,
        [1, 0, 0],
        [0, 1, 0],
        tilt_deg.to_radians(),
        twist_deg.to_radians(),
        0.0,
    )
    .unwrap();
    CrystalTarget::new(lat, o).with_thermal(Thermal::new(TEMPERATURE_K, THETA_D_SI))
}

/// Per-history substrate-entry bookkeeping and sums, merged in chunk order.
#[derive(Debug, Clone, PartialEq)]
struct OxideTally {
    // Inputs.
    thickness_m: f64,
    beam_dir: [f64; 3],
    beam_energy_ev: f64,
    /// Depths below the interface (m) at which the tail is counted.
    tail_m: Vec<f64>,
    // Per-history state.
    entered: bool,
    // Counts.
    histories: u64,
    max_relative_residual: f64,
    stopped_in_oxide: u64,
    backscattered_before_entry: u64,
    entered_count: u64,
    entered_stopped: u64,
    entered_backscattered: u64,
    /// Entered the substrate, came back into the oxide and stopped there.
    entered_returned_to_oxide: u64,
    other_exits: u64,
    /// Entries that had to use the first hook inside the substrate instead
    /// of the boundary hook (should be 0).
    fallback_entries: u64,
    // Sums over entered histories.
    e_sum: f64,
    e_sq: f64,
    th_sum: f64,
    th_sq: f64,
    ax_sum: f64,
    ax_sq: f64,
    // Sums over primaries stopped in the substrate.
    x_sum: f64,
    x_sq: f64,
    y_sum: f64,
    y_sq: f64,
    tail_counts: Vec<u64>,
}

impl OxideTally {
    fn new(thickness_m: f64, b: &Beam, tail_m: &[f64]) -> Self {
        let (s, c) = b.polar_rad.sin_cos();
        let (sa, ca) = b.azimuth_rad.sin_cos();
        Self {
            thickness_m,
            beam_dir: [c, s * ca, s * sa],
            beam_energy_ev: b.energy_ev,
            tail_m: tail_m.to_vec(),
            entered: false,
            histories: 0,
            max_relative_residual: 0.0,
            stopped_in_oxide: 0,
            backscattered_before_entry: 0,
            entered_count: 0,
            entered_stopped: 0,
            entered_backscattered: 0,
            entered_returned_to_oxide: 0,
            other_exits: 0,
            fallback_entries: 0,
            e_sum: 0.0,
            e_sq: 0.0,
            th_sum: 0.0,
            th_sq: 0.0,
            ax_sum: 0.0,
            ax_sq: 0.0,
            x_sum: 0.0,
            x_sq: 0.0,
            y_sum: 0.0,
            y_sq: 0.0,
            tail_counts: vec![0; tail_m.len()],
        }
    }

    fn enter(&mut self, energy_ev: f64, dir: [f64; 3]) {
        self.entered = true;
        self.entered_count += 1;
        let th = dir[0].clamp(-1.0, 1.0).acos();
        let dot = dir[0] * self.beam_dir[0] + dir[1] * self.beam_dir[1] + dir[2] * self.beam_dir[2];
        let ax = dot.clamp(-1.0, 1.0).acos();
        self.e_sum += energy_ev;
        self.e_sq += energy_ev * energy_ev;
        self.th_sum += th;
        self.th_sq += th * th;
        self.ax_sum += ax;
        self.ax_sq += ax * ax;
    }
}

impl BcaTally for OxideTally {
    fn begin_history(&mut self, _index: u64) {
        self.entered = false;
        if self.thickness_m == 0.0 {
            self.enter(self.beam_energy_ev, self.beam_dir);
        }
    }

    fn electronic(
        &mut self,
        p: &Particle,
        _from: [f64; 3],
        _channel: ElectronicChannel,
        _energy_ev: f64,
    ) {
        if !p.is_primary() || self.entered {
            return;
        }
        if p.layer == 0 {
            if (p.pos[0] - self.thickness_m).abs() < AT_INTERFACE && p.dir[0] > 0.0 {
                self.enter(p.energy_ev, p.dir);
            }
        } else {
            // First hook already inside the substrate: boundary hook missed.
            self.fallback_entries += 1;
            self.enter(p.energy_ev, p.dir);
        }
    }

    fn stopped(&mut self, p: &Particle) {
        if !p.is_primary() {
            return;
        }
        if p.layer == 0 && self.thickness_m > 0.0 {
            if self.entered {
                self.entered_returned_to_oxide += 1;
            } else {
                self.stopped_in_oxide += 1;
            }
            return;
        }
        if !self.entered {
            self.fallback_entries += 1;
            self.enter(p.energy_ev, p.dir);
        }
        self.entered_stopped += 1;
        let x = p.pos[0];
        let y = x - self.thickness_m;
        self.x_sum += x;
        self.x_sq += x * x;
        self.y_sum += y;
        self.y_sq += y * y;
        for (c, &t) in self.tail_counts.iter_mut().zip(&self.tail_m) {
            if y > t {
                *c += 1;
            }
        }
    }

    fn escaped(&mut self, p: &Particle, face: Face) {
        if !p.is_primary() {
            return;
        }
        match (face, self.entered) {
            (Face::Front, false) => self.backscattered_before_entry += 1,
            (Face::Front, true) => self.entered_backscattered += 1,
            _ => self.other_exits += 1,
        }
    }

    fn end_history(&mut self, _index: u64, budget: &EnergyBudget) {
        self.histories += 1;
        let r = (budget.residual() / budget.incident).abs();
        self.max_relative_residual = self.max_relative_residual.max(r);
    }

    fn merge(&mut self, o: Self) {
        self.histories += o.histories;
        self.max_relative_residual = self.max_relative_residual.max(o.max_relative_residual);
        self.stopped_in_oxide += o.stopped_in_oxide;
        self.backscattered_before_entry += o.backscattered_before_entry;
        self.entered_count += o.entered_count;
        self.entered_stopped += o.entered_stopped;
        self.entered_backscattered += o.entered_backscattered;
        self.entered_returned_to_oxide += o.entered_returned_to_oxide;
        self.other_exits += o.other_exits;
        self.fallback_entries += o.fallback_entries;
        self.e_sum += o.e_sum;
        self.e_sq += o.e_sq;
        self.th_sum += o.th_sum;
        self.th_sq += o.th_sq;
        self.ax_sum += o.ax_sum;
        self.ax_sq += o.ax_sq;
        self.x_sum += o.x_sum;
        self.x_sq += o.x_sq;
        self.y_sum += o.y_sum;
        self.y_sq += o.y_sq;
        for (a, b) in self.tail_counts.iter_mut().zip(o.tail_counts) {
            *a += b;
        }
    }
}

/// A value and its one-standard-error statistical uncertainty.
#[derive(Debug, Clone, Copy, Serialize)]
struct Est {
    value: f64,
    stderr: f64,
}

/// Mean of `n` samples from their sum and sum of squares, with the error of
/// the mean (`sigma / sqrt(n)`).
fn mean_est(sum: f64, sq: f64, n: u64) -> Est {
    if n == 0 {
        return Est {
            value: f64::NAN,
            stderr: f64::NAN,
        };
    }
    let nf = n as f64;
    let m = sum / nf;
    let var = (sq / nf - m * m).max(0.0);
    Est {
        value: m,
        stderr: (var / nf).sqrt(),
    }
}

/// Binomial fraction `k / n` with `sqrt(f (1 - f) / n)`.
fn frac_est(k: u64, n: u64) -> Est {
    if n == 0 {
        return Est {
            value: f64::NAN,
            stderr: f64::NAN,
        };
    }
    let f = k as f64 / n as f64;
    Est {
        value: f,
        stderr: (f * (1.0 - f) / n as f64).sqrt(),
    }
}

#[derive(Debug, Serialize)]
struct TailPoint {
    /// Threshold depth below the interface, nm.
    depth_below_interface_nm: f64,
    /// Of primaries that reached the substrate and stopped in it.
    fraction_of_substrate_stopped: Est,
    /// Of all incident primaries (includes oxide stopping and reflection).
    fraction_of_incident: Est,
}

#[derive(Debug, Serialize)]
struct Outcome {
    histories: u64,
    max_relative_energy_residual: f64,
    /// Fractions of incident primaries.
    stopped_in_oxide: Est,
    backscattered_before_entry: Est,
    entered_substrate: Est,
    /// Of entered primaries.
    entered_then_backscattered: Est,
    entered_then_stopped: Est,
    /// Of entered primaries: came back into the oxide and stopped there.
    entered_then_returned_to_oxide: Est,
    other_exits: u64,
    fallback_entries: u64,
    entry_energy_ev: Est,
    /// Polar angle of the entry direction to the surface normal, degrees.
    entry_angle_to_normal_deg: Est,
    /// Angle between the entry direction and the incident beam, degrees.
    entry_angle_to_beam_deg: Est,
    /// Projected range of substrate-stopped primaries from the entrance
    /// surface, nm.
    rp_from_surface_nm: Est,
    /// The same relative to the substrate interface, nm.
    rp_from_interface_nm: Est,
    /// Standard deviation of the depth below the interface, nm.
    straggling_from_interface_nm: f64,
    tail: Vec<TailPoint>,
}

impl OxideTally {
    fn outcome(&self) -> Outcome {
        let n = self.histories;
        let ne = self.entered_count;
        let ns = self.entered_stopped;
        let deg = 180.0 / std::f64::consts::PI;
        let sc = |e: Est, k: f64| Est {
            value: e.value * k,
            stderr: e.stderr * k,
        };
        let y = mean_est(self.y_sum, self.y_sq, ns);
        let straggling = (self.y_sq / ns as f64 - (self.y_sum / ns as f64).powi(2))
            .max(0.0)
            .sqrt();
        Outcome {
            histories: n,
            max_relative_energy_residual: self.max_relative_residual,
            stopped_in_oxide: frac_est(self.stopped_in_oxide, n),
            backscattered_before_entry: frac_est(self.backscattered_before_entry, n),
            entered_substrate: frac_est(ne, n),
            entered_then_backscattered: frac_est(self.entered_backscattered, ne),
            entered_then_stopped: frac_est(ns, ne),
            entered_then_returned_to_oxide: frac_est(self.entered_returned_to_oxide, ne),
            other_exits: self.other_exits,
            fallback_entries: self.fallback_entries,
            entry_energy_ev: mean_est(self.e_sum, self.e_sq, ne),
            entry_angle_to_normal_deg: sc(mean_est(self.th_sum, self.th_sq, ne), deg),
            entry_angle_to_beam_deg: sc(mean_est(self.ax_sum, self.ax_sq, ne), deg),
            rp_from_surface_nm: sc(mean_est(self.x_sum, self.x_sq, ns), 1.0 / NM),
            rp_from_interface_nm: sc(y, 1.0 / NM),
            straggling_from_interface_nm: straggling / NM,
            tail: self
                .tail_m
                .iter()
                .zip(&self.tail_counts)
                .map(|(&t, &k)| TailPoint {
                    depth_below_interface_nm: t / NM,
                    fraction_of_substrate_stopped: frac_est(k, ns),
                    fraction_of_incident: frac_est(k, n),
                })
                .collect(),
        }
    }
}

struct Run {
    ion_z: u8,
    energy_ev: f64,
    tilt_deg: f64,
    twist_deg: f64,
    thickness_nm: f64,
    crystal: bool,
    count: u64,
}

/// Run one configuration; returns the tally and the crystal metadata JSON.
fn run(r: &Run, tail_m: &[f64], threads: Option<usize>) -> (OxideTally, serde_json::Value) {
    let t = r.thickness_nm * NM;
    let (st, region) = if t > 0.0 {
        (Stack::new(vec![(sio2(), t)], Some(si())).unwrap(), 1)
    } else {
        (Stack::semi_infinite(si()), 0)
    };
    let b = Beam {
        ion: Ion::new(r.ion_z).unwrap(),
        energy_ev: r.energy_ev,
        polar_rad: r.tilt_deg.to_radians(),
        azimuth_rad: r.twist_deg.to_radians(),
        count: r.count,
    };
    let ls = LindhardScharff::new();
    let mut cfg = BcaConfig::new(5.0, 2.0);
    cfg.seed = SEED;
    let mut bca = Bca::new(b, &st, cfg, &ls, table()).unwrap();
    if r.crystal {
        bca = bca
            .with_crystal(crystal_si(r.tilt_deg, r.twist_deg), &[region])
            .unwrap();
    }
    let meta = serde_json::to_value(bca.crystal_metadata()).unwrap();
    let go = || {
        let proto = OxideTally::new(t, &b, tail_m);
        bca.run(|| proto.clone()).unwrap()
    };
    let tally = match threads {
        None => go(),
        Some(n) => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(go),
    };
    (tally, meta)
}

/// Small case in every `cargo test`: energy is conserved, every primary is
/// accounted for exactly once, and the tally is bit-identical on 1, 2 and 8
/// threads, for the crystal and for its amorphous control.
#[test]
fn small_screen_oxide_case_conserves_energy_and_is_thread_independent() {
    for crystal in [true, false] {
        let r = Run {
            ion_z: 5,
            energy_ev: 5.0e3,
            tilt_deg: 45.0,
            twist_deg: 0.0,
            thickness_nm: 3.0,
            crystal,
            count: 150,
        };
        let tail = [10.0 * NM, 20.0 * NM];
        let (t1, _) = run(&r, &tail, Some(1));
        assert!(t1.max_relative_residual < 1e-9, "crystal={crystal}");
        assert_eq!(t1.histories, 150);
        assert_eq!(t1.fallback_entries, 0, "crystal={crystal}");
        assert_eq!(t1.other_exits, 0);
        // Every primary is in exactly one class.
        assert_eq!(
            t1.stopped_in_oxide
                + t1.backscattered_before_entry
                + t1.entered_stopped
                + t1.entered_backscattered
                + t1.entered_returned_to_oxide,
            t1.histories
        );
        assert_eq!(
            t1.entered_count,
            t1.entered_stopped + t1.entered_backscattered + t1.entered_returned_to_oxide
        );
        assert!(t1.entered_count > 100, "most B 5 keV ions pass 3 nm");
        // Entry energy is below the incident energy after the oxide loss.
        assert!(t1.e_sum / (t1.entered_count as f64) < 5.0e3);
        for n in [2, 8] {
            let (tn, _) = run(&r, &tail, Some(n));
            assert!(t1 == tn, "crystal={crystal}: {n} threads differ");
        }
    }
}

/// A bare-crystal run has the incident beam as its entry state.
#[test]
fn bare_entry_state_is_the_incident_beam() {
    let r = Run {
        ion_z: 5,
        energy_ev: 5.0e3,
        tilt_deg: 30.0,
        twist_deg: 17.0,
        thickness_nm: 0.0,
        crystal: false,
        count: 40,
    };
    let (t, _) = run(&r, &[], Some(1));
    assert_eq!(t.entered_count, 40);
    let o = t.outcome();
    assert!((o.entry_energy_ev.value - 5.0e3).abs() < 1e-6);
    assert!((o.entry_angle_to_normal_deg.value - 30.0).abs() < 1e-9);
    assert!(o.entry_angle_to_beam_deg.value.abs() < 1e-6);
}

#[derive(Serialize)]
struct Record {
    ion_z: u8,
    energy_ev: f64,
    direction: &'static str,
    tilt_deg: f64,
    twist_deg: f64,
    oxide_thickness_nm: f64,
    substrate: &'static str,
    outcome: Outcome,
    /// Crystal metadata of the run (null for amorphous controls).
    crystal_metadata: serde_json::Value,
}

/// The sweep. Prints a summary table and writes the JSON results to
/// `LINDHARD_SCREEN_OXIDE_OUT` if set. Asserts only conservation and
/// accounting: the measured baselines are the result, and no suppression is
/// presupposed.
#[test]
#[ignore = "statistical; run with --release -- --ignored"]
fn screen_oxide_sweep() {
    let n: u64 = std::env::var("LINDHARD_SCREEN_OXIDE_N")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4000);
    let thicknesses = [0.0, 1.0, 2.0, 5.0, 10.0];
    let dirs: [(&str, f64, f64); 2] = [("aligned_110", 45.0, 0.0), ("off_axis_30_17", 30.0, 17.0)];
    let ions: [(u8, f64); 2] = [(5, 5.0e3), (33, 3.0e4)];
    let mut records = Vec::new();
    for (z, e) in ions {
        for (label, tilt, twist) in dirs {
            let mk = |thickness_nm: f64, crystal: bool| Run {
                ion_z: z,
                energy_ev: e,
                tilt_deg: tilt,
                twist_deg: twist,
                thickness_nm,
                crystal,
                count: n,
            };
            // Reference range: the bare amorphous control, same ion and
            // direction. The tail thresholds are 2 and 3 times its Rp.
            let (bare_am, _) = run(&mk(0.0, false), &[], None);
            let rp = bare_am.x_sum / bare_am.entered_stopped as f64;
            let tail = [2.0 * rp, 3.0 * rp];
            for &t in &thicknesses {
                for crystal in [true, false] {
                    let (tally, meta) = run(&mk(t, crystal), &tail, None);
                    assert!(tally.max_relative_residual < 1e-9);
                    assert_eq!(tally.fallback_entries, 0);
                    assert_eq!(
                        tally.stopped_in_oxide
                            + tally.backscattered_before_entry
                            + tally.entered_stopped
                            + tally.entered_backscattered
                            + tally.entered_returned_to_oxide
                            + tally.other_exits,
                        tally.histories
                    );
                    let o = tally.outcome();
                    println!(
                        "Z={z} {label:>14} T={t:>4} nm {}: oxide-stop {:.3} refl {:.3} \
                         enter {:.3}  E_in {:.0} eV  th {:.1} deg  axis {:.1} deg  \
                         Rp(int) {:.1}+-{:.1} nm  tail>2Rp {:.4}+-{:.4}",
                        if crystal { "crystal  " } else { "amorphous" },
                        o.stopped_in_oxide.value,
                        o.backscattered_before_entry.value,
                        o.entered_substrate.value,
                        o.entry_energy_ev.value,
                        o.entry_angle_to_normal_deg.value,
                        o.entry_angle_to_beam_deg.value,
                        o.rp_from_interface_nm.value,
                        o.rp_from_interface_nm.stderr,
                        o.tail[0].fraction_of_substrate_stopped.value,
                        o.tail[0].fraction_of_substrate_stopped.stderr,
                    );
                    records.push(Record {
                        ion_z: z,
                        energy_ev: e,
                        direction: label,
                        tilt_deg: tilt,
                        twist_deg: twist,
                        oxide_thickness_nm: t,
                        substrate: if crystal {
                            "crystalline_si"
                        } else {
                            "amorphous_si_matched_density"
                        },
                        outcome: o,
                        crystal_metadata: meta,
                    });
                }
            }
        }
    }
    let doc = serde_json::json!({
        "schema": "lindhard.screen_oxide.v1",
        "issue": 280,
        "model": {
            "engine": "Bca, constant free path, crystal flight model in the substrate",
            "screening": "ZblUniversal",
            "stopping": "LindhardScharff",
            "electronic_loss": "NonLocal (default): no local loss; the crystal local/nonlocal ratio deviation of #250 is not addressed",
            "config": {"primary_cutoff_ev": 5.0, "recoil_cutoff_ev": 2.0, "seed": SEED, "follow_recoils": true},
            "primaries_per_configuration": n,
            "temperature_k": TEMPERATURE_K,
            "debye_temperature_k": THETA_D_SI,
            "oxide": "amorphous SiO2, 2200 kg/m3, E_d = 15 eV (Si, O), amorphous transport",
            "substrate_controls": "amorphous Si at the crystal density, same stack and seed",
            "coordinates": "x = 0 entrance surface; interface at x = T; *_from_interface = x - T",
            "tail_thresholds": "2 and 3 times the bare amorphous control Rp (same ion and direction), below the interface",
            "uncertainty": "one standard error: sigma/sqrt(n) for means, sqrt(f(1-f)/n) for fractions; per-history, seeds fixed and shared across configurations"
        },
        "records": records,
    });
    if let Ok(path) = std::env::var("LINDHARD_SCREEN_OXIDE_OUT") {
        std::fs::write(&path, serde_json::to_string(&doc).unwrap() + "\n").unwrap();
        println!("wrote {path}");
    }
}
