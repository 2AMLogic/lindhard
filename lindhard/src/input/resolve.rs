//! Validation and resolution: turns an [`Input`] into the engine's types.

use crate::elements::element_by_symbol;
use crate::geometry::Stack;
use crate::ion::bca::{
    BcaConfig, Beam, CrystalTarget, ElectronicLoss, MeanFreePath, Thermal, MAX_WEAK_COLLISIONS,
};
use crate::ion::crystal::Divergence;
use crate::ion::stopping::table::StoppingTable;
use crate::ion::stopping::Ion;
use crate::material::{EnergyKind, Material};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::path::Path;

use super::resolved::*;
use super::schema::*;
use super::{invalid, Input, InputError};

impl Input {
    /// Validate `[dynamic]`.
    fn check_dynamic(&self, d: &DynamicSpec) -> Result<(), InputError> {
        if !finite_pos(d.fluence_cm2) {
            return Err(invalid(
                "dynamic.fluence_cm2",
                "must be finite and positive",
            ));
        }
        if d.ions_per_step == 0 {
            return Err(invalid("dynamic.ions_per_step", "must be at least 1"));
        }
        match d.max_change {
            Some(c) if !finite_pos(c) => {
                return Err(invalid("dynamic.max_change", "must be finite and positive"))
            }
            Some(_) => {
                if d.min_ions_per_step == 0 || d.min_ions_per_step > d.ions_per_step {
                    return Err(invalid(
                        "dynamic.min_ions_per_step",
                        "must be between 1 and ions_per_step",
                    ));
                }
            }
            None => {
                if d.min_ions_per_step != default_min_ions() {
                    return Err(invalid(
                        "dynamic.min_ions_per_step",
                        "only used with max_change (adaptive steps)",
                    ));
                }
            }
        }
        if let Some(t) = d.slab_nm {
            if !finite_pos(t) {
                return Err(invalid("dynamic.slab_nm", "must be finite and positive"));
            }
        }
        match (d.relaxation, d.number_density_cm3) {
            (RelaxationChoice::FixedNumberDensity, None) => {
                return Err(invalid(
                    "dynamic.number_density_cm3",
                    "required with relaxation = \"fixed-number-density\"",
                ))
            }
            (RelaxationChoice::FixedNumberDensity, Some(n)) if !finite_pos(n) => {
                return Err(invalid(
                    "dynamic.number_density_cm3",
                    "must be finite and positive",
                ))
            }
            (RelaxationChoice::IdealMixing, Some(_)) => {
                return Err(invalid(
                    "dynamic.number_density_cm3",
                    "only used with relaxation = \"fixed-number-density\"",
                ))
            }
            _ => {}
        }
        if d.relaxation == RelaxationChoice::FixedNumberDensity && !d.atomic_volume_nm3.is_empty() {
            return Err(invalid(
                "dynamic.atomic_volume_nm3",
                "only used with relaxation = \"ideal-mixing\"",
            ));
        }
        for (sym, v) in &d.atomic_volume_nm3 {
            let field = format!("dynamic.atomic_volume_nm3.{sym}");
            if element_by_symbol(sym).is_none() {
                return Err(invalid(&field, format!("unknown element symbol {sym:?}")));
            }
            if !finite_pos(*v) {
                return Err(invalid(&field, "must be finite and positive"));
            }
        }
        for (sym, o) in &d.energies {
            let field = format!("dynamic.energies.{sym}");
            if element_by_symbol(sym).is_none() {
                return Err(invalid(&field, format!("unknown element symbol {sym:?}")));
            }
            for (key, v) in [
                ("e_d_ev", o.e_d_ev),
                ("e_b_ev", o.e_b_ev),
                ("e_s_ev", o.e_s_ev),
            ] {
                if v.is_some_and(|v| !(v.is_finite() && v >= 0.0)) {
                    return Err(invalid(
                        format!("{field}.{key}"),
                        "must be finite and non-negative",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Read and validate every `[stopping]` table file.
    fn load_stopping_tables(
        &self,
        base_dir: &Path,
    ) -> Result<Vec<LoadedStoppingTable>, InputError> {
        let mut loaded: Vec<LoadedStoppingTable> = Vec::new();
        for (i, path) in self.stopping.tables.iter().enumerate() {
            let field = format!("stopping.tables[{i}]");
            let joined = base_dir.join(path);
            let bytes = std::fs::read(&joined)
                .map_err(|e| invalid(&field, format!("cannot read {}: {e}", joined.display())))?;
            let text = std::str::from_utf8(&bytes).map_err(|_| {
                invalid(
                    &field,
                    format!("{} is not valid UTF-8 text", joined.display()),
                )
            })?;
            let table = StoppingTable::from_toml_str(text)
                .map_err(|e| invalid(&field, format!("{path}: {e}")))?;
            if let Some(prev) = loaded
                .iter()
                .position(|l| l.table.covers_pair(table.ion_z(), table.target_z()))
            {
                return Err(invalid(
                    &field,
                    format!(
                        "{path} declares the same ion/target pair (Z1={}, Z2={}) as \
                         stopping.tables[{prev}]; give each pair one table",
                        table.ion_z(),
                        table.target_z()
                    ),
                ));
            }
            let sha256 = Sha256::digest(&bytes)
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect();
            loaded.push(LoadedStoppingTable {
                path: path.clone(),
                resolved_path: std::fs::canonicalize(&joined).unwrap_or(joined),
                sha256,
                table,
            });
        }
        Ok(loaded)
    }

    /// Checks of the loaded tables against the rest of the run.
    fn check_stopping_tables(
        &self,
        tables: &[LoadedStoppingTable],
        ion: &Ion,
        layers: &[ResolvedLayer],
        warnings: &mut Vec<String>,
    ) -> Result<(), InputError> {
        if tables.is_empty() {
            return Ok(());
        }
        if self.physics.stopping == StoppingChoice::EquipartitionLsOr {
            return Err(invalid(
                "stopping.tables",
                "not usable with physics.stopping = \"equipartition-ls-or\", which carries \
                 its own Lindhard-Scharff/Oen-Robinson loss and would ignore the tables; \
                 use \"lindhard-scharff\" (or another model) as the fallback",
            ));
        }
        let e0 = self.beam.energy_ev;
        let cutoff = self.physics.primary_cutoff_ev;
        let recoil_cutoff = self.physics.recoil_cutoff_ev;
        let follow = self.physics.follow_recoils;
        let tol = crate::ion::stopping::table::MASS_RELATIVE_TOLERANCE;
        let target_elements: std::collections::BTreeSet<u8> = layers
            .iter()
            .flat_map(|ly| ly.material.components().iter().map(|c| c.z()))
            .collect();
        let symbol = |z: u8| crate::elements::element(z).map_or("?", |e| e.symbol);
        for (i, l) in tables.iter().enumerate() {
            let field = format!("stopping.tables[{i}]");
            let t = &l.table;
            let in_target = target_elements.contains(&t.target_z());
            // A table serves the beam ion of its element and, when recoils
            // are followed, every recoil of its element (`TableOverride`
            // keys on Z1, Z2 only).
            let serves_beam = t.ion_z() == ion.z();
            let serves_recoils = follow && target_elements.contains(&t.ion_z());
            if !in_target || !(serves_beam || serves_recoils) {
                let recoils = if follow {
                    ""
                } else {
                    " and physics.follow_recoils = false"
                };
                warnings.push(format!(
                    "{field}: the table ({}->{}, Z1={}, Z2={}) is for a pair that does not \
                     occur in this run (beam {}, target elements {}{recoils}); it is unused",
                    symbol(t.ion_z()),
                    symbol(t.target_z()),
                    t.ion_z(),
                    t.target_z(),
                    self.beam.ion,
                    target_elements
                        .iter()
                        .map(|&z| symbol(z))
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
                continue;
            }
            let (lo, hi) = t.energy_range_ev();
            if serves_recoils {
                // Recoils are built with the standard atomic weight
                // (`Ion::new`), so a table for another mass can never serve
                // them, and the run would stop at the first such recoil.
                let recoil = Ion::new(t.ion_z()).map_err(|e| invalid(&field, e.to_string()))?;
                if (recoil.mass_amu() - t.ion_mass_amu()).abs() > tol * t.ion_mass_amu() {
                    return Err(invalid(
                        &field,
                        format!(
                            "{}: table is for an ion of mass {} u, but with \
                             physics.follow_recoils = true it also serves the {} recoils, \
                             which have the standard atomic weight {} u; a table serves one \
                             mass per pair (use the standard weight in the table and the beam, \
                             drop the table, or set physics.follow_recoils = false)",
                            l.path,
                            t.ion_mass_amu(),
                            symbol(t.ion_z()),
                            recoil.mass_amu()
                        ),
                    ));
                }
                // Every recoil is followed down to recoil_cutoff_ev.
                if recoil_cutoff < lo {
                    return Err(invalid(
                        &field,
                        format!(
                            "{}: the table starts at {lo} eV, above physics.recoil_cutoff_ev = \
                             {recoil_cutoff} eV; with physics.follow_recoils = true it serves \
                             the {} recoils, which slow down to the recoil cutoff, so the run \
                             would stop outside the table range (extend the table down to the \
                             cutoff or raise physics.recoil_cutoff_ev)",
                            l.path,
                            symbol(t.ion_z()),
                        ),
                    ));
                }
                // The most a beam collision can hand a recoil of this element.
                let (m1, m2) = (ion.mass_amu(), recoil.mass_amu());
                let e_max = 4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2)) * e0;
                if hi < e_max {
                    warnings.push(format!(
                        "{field}: the table ends at {hi} eV, below the largest energy a \
                         {} recoil can receive from the beam ({e_max:.4e} eV); the run fails \
                         if a recoil starts above the table range",
                        symbol(t.ion_z())
                    ));
                }
            }
            if !serves_beam {
                continue;
            }
            // The beam ion's pair: mass and the energy range it travels.
            if (ion.mass_amu() - t.ion_mass_amu()).abs() > tol * t.ion_mass_amu() {
                return Err(invalid(
                    &field,
                    format!(
                        "{}: table is for an ion of mass {} u but the beam ion has mass {} u \
                         (set `ion_mass_amu` in the table, or beam.mass_amu)",
                        l.path,
                        t.ion_mass_amu(),
                        ion.mass_amu()
                    ),
                ));
            }
            if e0 < lo || e0 > hi {
                return Err(invalid(
                    &field,
                    format!(
                        "{}: beam energy {e0} eV is outside the table range [{lo}, {hi}] eV",
                        l.path
                    ),
                ));
            }
            if cutoff < lo {
                warnings.push(format!(
                    "{field}: the table starts at {lo} eV, above physics.primary_cutoff_ev = \
                     {cutoff} eV; the run fails if a projectile slows below the table range"
                ));
            }
        }
        Ok(())
    }

    fn material(&self, field: &str, r: &MaterialRef) -> Result<ResolvedLayer, InputError> {
        resolve_material(&self.materials, field, r)
    }

    /// Validate everything and build the engine's inputs. Errors name the
    /// offending key.
    pub fn resolve(&self) -> Result<Resolved, InputError> {
        self.resolve_in(Path::new("."))
    }

    /// Like [`Input::resolve`], with relative `[stopping]` table paths taken
    /// relative to `base_dir` (the directory of the input file).
    pub fn resolve_in(&self, base_dir: &Path) -> Result<Resolved, InputError> {
        self.resolve_with_sets(base_dir, TUNING_SETS)
    }

    /// [`Input::resolve_in`] against an explicit registry of tuning sets
    /// (tests use synthetic sets; production uses [`TUNING_SETS`]).
    pub fn resolve_with_sets(
        &self,
        base_dir: &Path,
        sets: &[TuningSet],
    ) -> Result<Resolved, InputError> {
        let mut warnings = Vec::new();

        // Beam.
        let b = &self.beam;
        let el = element_by_symbol(&b.ion).ok_or_else(|| {
            invalid(
                "beam.ion",
                format!(
                    "unknown element symbol {:?} (symbols are case-sensitive)",
                    b.ion
                ),
            )
        })?;
        let ion = match b.mass_amu {
            None => Ion::new(el.z),
            Some(m) => Ion::with_mass(el.z, m),
        }
        .map_err(|e| invalid("beam.mass_amu", e.to_string()))?;
        if !finite_pos(b.energy_ev) {
            return Err(invalid("beam.energy_ev", "must be finite and positive"));
        }
        if !(b.tilt_deg.is_finite() && (0.0..90.0).contains(&b.tilt_deg)) {
            return Err(invalid("beam.tilt_deg", "must be in [0, 90) degrees"));
        }
        if !b.azimuth_deg.is_finite() {
            return Err(invalid("beam.azimuth_deg", "must be finite"));
        }
        if self.run.ions == 0 {
            return Err(invalid("run.ions", "must be at least 1"));
        }
        if self.run.threads == Some(0) {
            return Err(invalid(
                "run.threads",
                "must be at least 1 (omit for all cores)",
            ));
        }
        let divergence = self.resolve_divergence()?;
        let beam = Beam {
            ion,
            energy_ev: b.energy_ev,
            polar_rad: b.tilt_deg.to_radians(),
            azimuth_rad: b.azimuth_deg.to_radians(),
            count: self.run.ions,
        };

        // Materials (all of them, so a broken unused entry is still reported).
        for (name, spec) in &self.materials {
            Material::try_from(spec.clone())
                .map_err(|e| invalid(format!("materials.{name}"), e.to_string()))?;
        }

        // Target.
        let t = &self.target;
        if t.layers.is_empty() && t.substrate.is_none() {
            return Err(invalid(
                "target",
                "needs at least one [[target.layers]] entry or a substrate",
            ));
        }
        let mut layers = Vec::with_capacity(t.layers.len() + 1);
        let mut fields = Vec::with_capacity(t.layers.len() + 1);
        let mut thick = Vec::with_capacity(t.layers.len());
        for (i, l) in t.layers.iter().enumerate() {
            let field = format!("target.layers[{i}]");
            if !finite_pos(l.thickness_nm) {
                return Err(invalid(
                    format!("{field}.thickness_nm"),
                    format!("{} nm must be finite and positive", l.thickness_nm),
                ));
            }
            layers.push(self.material(&format!("{field}.material"), &l.material)?);
            fields.push(format!("{field}.material"));
            thick.push(l.thickness_nm * 1e-9);
        }
        if let Some(s) = &t.substrate {
            layers.push(self.material("target.substrate", s)?);
            fields.push("target.substrate".to_string());
        }

        if let Some(d) = &self.dynamic {
            self.check_dynamic(d)?;
        }

        // Energy overrides.
        let p = &self.physics;
        for (sym, o) in &p.energies {
            let field = format!("physics.energies.{sym}");
            let el = element_by_symbol(sym)
                .ok_or_else(|| invalid(&field, format!("unknown element symbol {sym:?}")))?;
            let mut used = false;
            for l in &mut layers {
                if l.material.components().iter().all(|c| c.z() != el.z) {
                    continue;
                }
                used = true;
                let m = &mut l.material;
                let set = |r: Result<(), _>, key: &str| {
                    r.map_err(|e: crate::material::MaterialError| {
                        invalid(format!("{field}.{key}"), e.to_string())
                    })
                };
                if let Some(v) = o.e_d_ev {
                    set(m.set_displacement_energy_ev(el.z, v), "e_d_ev")?;
                }
                if let Some(v) = o.e_b_ev {
                    set(m.set_lattice_binding_energy_ev(el.z, v), "e_b_ev")?;
                }
                if let Some(v) = o.e_s_ev {
                    set(m.set_surface_binding_energy_ev(el.z, v), "e_s_ev")?;
                }
            }
            if !used {
                return Err(invalid(&field, format!("{sym} is not in any target layer")));
            }
        }
        for (l, field) in layers.iter().zip(&fields) {
            if let Some(&(z, kind)) = l.material.unset_energies().first() {
                let sym = crate::elements::element(z).expect("validated").symbol;
                let key = match kind {
                    EnergyKind::Displacement => "e_d_ev",
                    EnergyKind::LatticeBinding => "e_b_ev",
                    EnergyKind::SurfaceBinding => "e_s_ev",
                };
                return Err(invalid(
                    field.as_str(),
                    format!(
                        "{kind} of {sym} has no default; set `{key}` for {sym} in the \
                         material, or in [physics.energies.{sym}]"
                    ),
                ));
            }
        }

        let tuning = self.apply_tuning(sets, &mut layers, &mut warnings)?;

        let mut finite = Vec::with_capacity(thick.len());
        let mut it = layers.iter();
        for &th in &thick {
            finite.push((it.next().expect("one per layer").material.clone(), th));
        }
        let substrate = it.next().map(|l| l.material.clone());
        let stack = Stack::new(finite, substrate).map_err(|e| invalid("target", e.to_string()))?;

        // Physics.
        let mean_free_path = match (p.free_path, p.min_cm_angle_deg) {
            (FreePathChoice::Constant, None) => MeanFreePath::Constant,
            (FreePathChoice::Constant, Some(_)) => {
                return Err(invalid(
                    "physics.min_cm_angle_deg",
                    "only used with free_path = \"energy-dependent\"",
                ))
            }
            (FreePathChoice::EnergyDependent, None) => {
                return Err(invalid(
                    "physics.min_cm_angle_deg",
                    "required with free_path = \"energy-dependent\"",
                ))
            }
            (FreePathChoice::EnergyDependent, Some(a)) => {
                if !(a > 0.0 && a < 180.0) {
                    return Err(invalid(
                        "physics.min_cm_angle_deg",
                        format!("{a} must be in (0, 180) degrees"),
                    ));
                }
                MeanFreePath::EnergyDependent {
                    min_cm_angle_rad: a * PI / 180.0,
                }
            }
        };
        for (key, v) in [
            ("primary_cutoff_ev", p.primary_cutoff_ev),
            ("recoil_cutoff_ev", p.recoil_cutoff_ev),
        ] {
            if !finite_pos(v) {
                return Err(invalid(
                    format!("physics.{key}"),
                    format!("{v} must be finite and positive"),
                ));
            }
        }
        if !(p.primary_surface_binding_ev.is_finite() && p.primary_surface_binding_ev >= 0.0) {
            return Err(invalid(
                "physics.primary_surface_binding_ev",
                "must be finite and non-negative",
            ));
        }
        if p.primary_cutoff_ev >= b.energy_ev {
            return Err(invalid(
                "physics.primary_cutoff_ev",
                format!(
                    "{} eV is not below the beam energy {} eV",
                    p.primary_cutoff_ev, b.energy_ev
                ),
            ));
        }
        if p.follow_recoils {
            let min_es = layers
                .iter()
                .flat_map(|l| {
                    l.material
                        .components()
                        .iter()
                        .map(|c| l.material.surface_binding_energy_ev(c.z()).expect("set"))
                })
                .fold(f64::INFINITY, f64::min);
            if p.recoil_cutoff_ev > min_es {
                warnings.push(format!(
                    "physics.recoil_cutoff_ev = {} eV is above the smallest surface binding \
                     energy ({min_es} eV); sputtering will be underestimated",
                    p.recoil_cutoff_ev
                ));
            }
        }
        if p.weak_collisions > MAX_WEAK_COLLISIONS {
            return Err(invalid(
                "physics.weak_collisions",
                format!(
                    "{} must be at most {MAX_WEAK_COLLISIONS}",
                    p.weak_collisions
                ),
            ));
        }
        if p.weak_collisions > 0 && mean_free_path != MeanFreePath::Constant {
            return Err(invalid(
                "physics.weak_collisions",
                "only used with free_path = \"constant\"",
            ));
        }
        let mut config = BcaConfig::new(p.primary_cutoff_ev, p.recoil_cutoff_ev);
        config.mean_free_path = mean_free_path;
        config.weak_collisions = p.weak_collisions;
        config.electronic = match p.stopping {
            StoppingChoice::EquipartitionLsOr => ElectronicLoss::EquipartitionLsOr,
            _ => ElectronicLoss::NonLocal,
        };
        config.follow_recoils = p.follow_recoils;
        config.primary_surface_binding_ev = p.primary_surface_binding_ev;
        config.seed = self.run.seed;

        // Tally.
        if !finite_pos(self.tally.depth_bin_nm) {
            return Err(invalid("tally.depth_bin_nm", "must be finite and positive"));
        }
        if self.tally.depth_bins == 0 {
            return Err(invalid("tally.depth_bins", "must be at least 1"));
        }
        if !finite_pos(self.tally.lateral_bin_nm) {
            return Err(invalid(
                "tally.lateral_bin_nm",
                "must be finite and positive",
            ));
        }
        for (key, n) in [
            ("lateral_bins", self.tally.lateral_bins),
            ("escape_energy_bins", self.tally.escape_energy_bins),
            ("escape_polar_bins", self.tally.escape_polar_bins),
        ] {
            if n == 0 {
                return Err(invalid(format!("tally.{key}"), "must be at least 1"));
            }
        }
        if let Some(e) = self.tally.escape_energy_max_ev {
            if !finite_pos(e) {
                return Err(invalid(
                    "tally.escape_energy_max_ev",
                    "must be finite and positive",
                ));
            }
        }

        let crystals = self.resolve_crystals(&layers, &mut warnings)?;

        let screening = p.potential.screening();
        let screening_length = p
            .screening_length
            .map_or(screening.default_length(), LengthChoice::length);
        let stopping_tables = self.load_stopping_tables(base_dir)?;
        let fallback = stopping_model(p.stopping);
        let validity = fallback.validity(&ion);
        // The advisory range of the fallback model only matters for the pairs
        // it actually serves; a beam pair with a table has its range checked
        // as an error above.
        let beam_pair_uncovered = layers.iter().any(|l| {
            l.material.components().iter().any(|c| {
                !stopping_tables
                    .iter()
                    .any(|t| t.table.covers_pair(ion.z(), c.z()))
            })
        });
        if beam_pair_uncovered && !validity.contains(b.energy_ev) {
            warnings.push(format!(
                "beam energy {} eV is outside the advisory validity range of the {} \
                 stopping model for {} ([{}, {}] eV)",
                b.energy_ev,
                fallback.name(),
                b.ion,
                validity.min_energy_ev,
                validity.max_energy_ev
            ));
        }
        self.check_stopping_tables(&stopping_tables, &ion, &layers, &mut warnings)?;

        Ok(Resolved {
            input: self.echo(),
            beam,
            divergence,
            stack,
            layers,
            config,
            screening,
            screening_length,
            table_spec: TABLE_SPEC,
            stopping_tables,
            crystals,
            warnings,
            tuning,
        })
    }

    /// Validate `[beam.divergence]` and convert it to the engine's
    /// [`Divergence`]. Static ion runs only: a `[dynamic]` target is
    /// rejected, since the spread is not threaded through the fluence loop.
    fn resolve_divergence(&self) -> Result<Divergence, InputError> {
        let Some(spec) = &self.beam.divergence else {
            return Ok(Divergence::None);
        };
        if self.dynamic.is_some() {
            return Err(invalid(
                "beam.divergence",
                "not supported with a [dynamic] target; use a static run",
            ));
        }
        let (field, width_deg) = match *spec {
            DivergenceSpec::Gaussian { sigma_deg } => ("beam.divergence.sigma_deg", sigma_deg),
            DivergenceSpec::UniformCone { half_angle_deg } => {
                ("beam.divergence.half_angle_deg", half_angle_deg)
            }
        };
        if !(width_deg.is_finite() && (0.0..=MAX_DIVERGENCE_DEG).contains(&width_deg)) {
            return Err(invalid(
                field,
                format!(
                    "must be finite and in [0, {MAX_DIVERGENCE_DEG}] degrees (small-angle \
                     spreads only)"
                ),
            ));
        }
        let rad = width_deg.to_radians();
        Ok(match spec {
            DivergenceSpec::Gaussian { .. } => Divergence::Gaussian { sigma_rad: rad },
            DivergenceSpec::UniformCone { .. } => Divergence::UniformCone {
                half_angle_rad: rad,
            },
        })
    }

    /// Validate `[[crystal]]` and build the engine's crystal targets. The
    /// beam's tilt and azimuth are the orientation's tilt and twist. Every
    /// check the engine would make at [`crate::ion::bca::Bca::with_crystal`]
    /// that does not need an engine is made here, so `check` reports it.
    fn resolve_crystals(
        &self,
        layers: &[ResolvedLayer],
        warnings: &mut Vec<String>,
    ) -> Result<Vec<ResolvedCrystal>, InputError> {
        if self.crystal.is_empty() {
            return Ok(Vec::new());
        }
        if self.dynamic.is_some() {
            return Err(invalid(
                "crystal",
                "not supported with a [dynamic] target (the composition changes during the run \
                 while the lattice is fixed); use a static run",
            ));
        }
        let p = &self.physics;
        if p.free_path != FreePathChoice::Constant {
            return Err(invalid(
                "physics.free_path",
                "crystal regions need free_path = \"constant\" (the lattice replaces the free path)",
            ));
        }
        if p.weak_collisions > 0 {
            return Err(invalid(
                "physics.weak_collisions",
                "crystal regions need weak_collisions = 0 (the lattice replaces them)",
            ));
        }
        if p.tuning != NO_TUNING {
            return Err(invalid(
                "physics.tuning",
                "not supported together with [[crystal]] (the pilot sets were fitted for amorphous runs)",
            ));
        }
        if p.stopping == StoppingChoice::EquipartitionLsOr {
            warnings.push(
                "crystal: physics.stopping = \"equipartition-ls-or\" takes the local \
                 Oen-Robinson loss, whose constants are not verified against the paper \
                 (summary: physics.crystal[].electronic_constants_unverified); channeled \
                 ranges depend on them"
                    .to_string(),
            );
        }
        let beam = &self.beam;
        let mut owner: Vec<Option<usize>> = vec![None; layers.len()];
        let mut out = Vec::with_capacity(self.crystal.len());
        for (i, c) in self.crystal.iter().enumerate() {
            let f = |k: &str| format!("crystal[{i}].{k}");
            if c.layers.is_empty() {
                return Err(invalid(f("layers"), "needs at least one layer index"));
            }
            for &l in &c.layers {
                if l >= layers.len() {
                    return Err(invalid(
                        f("layers"),
                        format!(
                            "layer {l} does not exist: the target has {} layer(s), indices 0..{}",
                            layers.len(),
                            layers.len()
                        ),
                    ));
                }
                if let Some(prev) = owner[l] {
                    return Err(invalid(
                        f("layers"),
                        format!("layer {l} is already assigned by crystal[{prev}]"),
                    ));
                }
                owner[l] = Some(i);
            }
            if !c.wafer_rotation_deg.is_finite() {
                return Err(invalid(f("wafer_rotation_deg"), "must be finite"));
            }
            let lattice = c.preset.lattice();
            // Index triples, then the zone law, with the engine's own messages.
            let orientation = crate::ion::crystal::Orientation::new(
                &lattice,
                c.normal,
                c.reference,
                beam.tilt_deg.to_radians(),
                beam.azimuth_deg.to_radians(),
                c.wafer_rotation_deg.to_radians(),
            )
            .map_err(|e| {
                use crate::ion::crystal::CrystalError as E;
                let key = match &e {
                    E::ReferenceNotInPlane { .. } => "reference",
                    E::ZeroIndex(v) if *v == c.reference => "reference",
                    _ => "normal",
                };
                invalid(f(key), e.to_string())
            })?;
            let mut target = CrystalTarget::new(lattice, orientation);
            for (key, v, slot) in [
                ("p_max_nm", c.p_max_nm, &mut target.p_max_m),
                ("q_max_nm", c.q_max_nm, &mut target.q_max_m),
                (
                    "search_length_nm",
                    c.search_length_nm,
                    &mut target.search_length_m,
                ),
            ] {
                if let Some(v) = v {
                    if !finite_pos(v) {
                        return Err(invalid(f(key), "must be finite and positive"));
                    }
                    *slot = Some(v * 1e-9);
                }
            }
            if let Some(t) = &c.thermal {
                if !(t.temperature_k.is_finite() && t.temperature_k >= 0.0) {
                    return Err(invalid(
                        f("thermal.temperature_k"),
                        "must be finite and non-negative",
                    ));
                }
                let theta = t
                    .debye_temperature_k
                    .unwrap_or_else(|| c.preset.default_debye_temperature_k());
                if !finite_pos(theta) {
                    return Err(invalid(
                        f("thermal.debye_temperature_k"),
                        "must be finite and positive",
                    ));
                }
                target = target.with_thermal(Thermal::new(t.temperature_k, theta));
            }
            for &l in &c.layers {
                target.check_material(&layers[l].material).map_err(|m| {
                    invalid(
                        f("preset"),
                        format!(
                            "{} does not fit layer {l} ({}): {m}",
                            c.preset.name(),
                            layers[l].source
                        ),
                    )
                })?;
            }
            out.push(ResolvedCrystal {
                regions: c.layers.clone(),
                target,
            });
        }
        Ok(out)
    }

    /// Apply `[physics] tuning` to the resolved layers (see [`TuningSet`]).
    fn apply_tuning(
        &self,
        sets: &[TuningSet],
        layers: &mut [ResolvedLayer],
        warnings: &mut Vec<String>,
    ) -> Result<Option<TuningReport>, InputError> {
        let name = self.physics.tuning.as_str();
        if name == NO_TUNING {
            return Ok(None);
        }
        let f = "physics.tuning";
        let set = sets.iter().find(|s| s.name == name).ok_or_else(|| {
            let known: Vec<&str> = sets.iter().map(|s| s.name).collect();
            invalid(
                f,
                format!("unknown tuning set {name:?}; use \"none\" or one of {known:?}"),
            )
        })?;
        if self.dynamic.is_some() {
            return Err(invalid(
                f,
                "tuning is not supported with a [dynamic] target (pilot: static single-element layers)",
            ));
        }
        let (e_lo, e_hi) = set.energy_range_ev;
        if !(e_lo.is_finite() && e_hi.is_finite() && 0.0 < e_lo && e_lo <= e_hi) {
            return Err(invalid(
                f,
                format!(
                    "set {name:?}: energy range {e_lo}..{e_hi} eV is not a finite positive range"
                ),
            ));
        }
        if !set.ions.contains(&self.beam.ion.as_str()) {
            return Err(invalid(
                f,
                format!(
                    "set {name:?} was fitted for beam ion(s) {:?}, not {:?}; use \"none\"",
                    set.ions, self.beam.ion
                ),
            ));
        }
        let e = self.beam.energy_ev;
        if !(e_lo..=e_hi).contains(&e) {
            warnings.push(format!(
                "tuning set {name:?} was fitted at {e_lo}..{e_hi} eV; {e} eV is an extrapolation"
            ));
        }
        if self.beam.tilt_deg != 0.0 {
            warnings.push(format!(
                "tuning set {name:?} was fitted at normal incidence; tilt {} deg is an extrapolation",
                self.beam.tilt_deg
            ));
        }
        let mut factors = BTreeMap::new();
        for &(sym, k) in set.e_s_factors {
            let el = element_by_symbol(sym)
                .ok_or_else(|| invalid(f, format!("set {name:?} names unknown element {sym:?}")))?;
            if !(k.is_finite() && k > 0.0) {
                return Err(invalid(
                    f,
                    format!("set {name:?}: factor for {sym} is {k}; must be finite and positive"),
                ));
            }
            if factors.insert(el.z, k).is_some() {
                return Err(invalid(f, format!("set {name:?} lists {sym} twice")));
            }
        }
        for (i, l) in layers.iter().enumerate() {
            if l.material.components().len() != 1 {
                return Err(invalid(
                    f,
                    format!(
                        "tuning supports single-element layers only; layer {i} ({}) is a compound",
                        l.source
                    ),
                ));
            }
        }
        let mut components = Vec::with_capacity(layers.len());
        for (i, l) in layers.iter_mut().enumerate() {
            let z = l.material.components()[0].z();
            let sym = crate::elements::element(z).expect("validated").symbol;
            let original = l.material.surface_binding_energy_ev(z).expect("set");
            let factor = *factors.get(&z).ok_or_else(|| {
                invalid(
                    f,
                    format!(
                        "set {name:?} has no factor for {sym} (layer {i}); it supports {:?}",
                        set.e_s_factors.iter().map(|p| p.0).collect::<Vec<_>>()
                    ),
                )
            })?;
            let effective = original * factor;
            l.material
                .set_surface_binding_energy_ev(z, effective)
                .map_err(|e| invalid(f, format!("layer {i}: effective E_s of {sym}: {e}")))?;
            components.push(TunedEnergy {
                layer: i,
                element: sym.to_string(),
                e_s_original_ev: original,
                factor,
                e_s_effective_ev: effective,
            });
        }
        if components.iter().all(|c| c.factor == 1.0) {
            warnings.push(format!(
                "tuning set {name:?} has factor 1 for every element in the target; nothing changed"
            ));
        }
        Ok(Some(TuningReport {
            set: set.name.to_string(),
            version: set.version,
            provenance: set.provenance.to_string(),
            quantity: "surface-binding-energy",
            components,
        }))
    }
}

fn finite_pos(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
