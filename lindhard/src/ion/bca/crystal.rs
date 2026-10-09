//! The crystal flight model: collision partners taken from explicit lattice
//! sites instead of a random impact parameter (step 20b of the M2 crystal
//! plan), optionally displaced by thermal vibration (step 21b).
//!
//! # Choosing the model per region
//!
//! [`Geometry`] is unchanged. A [`Bca`] is amorphous everywhere unless
//! [`Bca::with_crystal`] names a [`CrystalTarget`] (lattice and orientation)
//! and the list of **regions** (layers of a stack, voxels of a grid) it
//! fills. A particle whose current region is in that list uses the crystal
//! model; in every other region, and for a [`Bca`] without `with_crystal`, the
//! amorphous model is used unchanged (and bit for bit: the amorphous code path
//! draws the same random numbers and does the same arithmetic as before; the
//! test `tests/crystal_amorphous_identity.rs` compares complete reports).
//! Several calls give several crystals (different lattices or orientations)
//! in different regions.
//!
//! The material of a crystal region still supplies the displacement and
//! binding energies, the stopping and the surface barrier; it must contain
//! every element of the lattice and its atom density must agree with the
//! lattice's to 5 %. The crystal model needs the constant-path convention and
//! no weak collisions (the lattice replaces both). Beam direction and surface
//! stay in the lab frame of the geometry; the orientation's lab-to-crystal
//! rotation ([`Orientation::to_crystal`]) is the only thing read from it (the
//! beam is the [`Beam`](super::Beam)'s, so give the `Beam` the orientation's tilt and
//! twist, as `Orientation::beam_lab` does).
//!
//! # Lattice frame of a history
//!
//! The crystal position of a lab point `r` is `R r + u`, with `R` the
//! orientation's rotation and `u` a translation drawn uniformly from the cube
//! `[0, a)^3` once per primary history (three uniform numbers per crystal,
//! drawn from a copy of the history's stream positioned at word 2^67, so the
//! transport draws are the same with or without crystals), shared by the
//! primary and all its recoils. A beam of finite
//! width samples all positions of the unit cell, so this is the static
//! equivalent of a random entry point; it also makes the surface termination
//! random. Lattice sites that lie outside the target (in front of the surface)
//! or in a region that is not part of this crystal are not partners.
//!
//! # Collision sequence (after DISPLATH, Tier A)
//!
//! The ion moves in straight segments. For a segment from position `r0` along
//! the unit direction `d` the partners are the lattice sites within `p_max`
//! of the line, ordered by the path distance `s = (r - r0) . d` of their
//! point of closest approach ([`crate::ion::crystal::LatticeSearch`]). The
//! scheme is the one of DISPLATH (<https://github.com/permissionx/DISPLATH>,
//! MIT, commit `7f461141c518305e6a5de8ce7e37dd314cd20ccb`, `src/dynamics.jl`:
//! `GetTargetsFromNeighbor` and `Collision!`, and `SimultaneousCriteria` in
//! `src/geometry.jl`), re-implemented here on [`LatticeSearch`] with
//! the lattice-site list of this crate, and with the project's own scattering
//! tables, kinematics and energy bookkeeping. DISPLATH is Copyright (c) 2024
//! 裴幂許Permission (MIT); the notice is in `THIRD_PARTY_LICENSES.md`.
//!
//! * **Free path.** The ion flies straight to the closest-approach point of
//!   the nearest partner (smallest `s`, ties broken by the search's documented
//!   order). There is no free-path length parameter: the path is the crystal's.
//!   The nonlocal electronic loss of the model chosen in the [`BcaConfig`](super::BcaConfig) is
//!   taken over that length at the energy at the start of the segment, as in
//!   the amorphous model. Under
//!   [`ElectronicLoss::EquipartitionLsOr`](super::ElectronicLoss::EquipartitionLsOr)
//!   that is half of the Lindhard-Scharff loss, and the other half is local
//!   (next bullets and "Local electronic loss"). When no partner lies within
//!   `search_length` the ion flies that distance and searches again (a
//!   channel). The geometry is asked about every segment, so surfaces and
//!   interfaces cut it exactly like an amorphous flight.
//! * **Simultaneous collisions.** Further sites join the nearest one if, for
//!   every partner `t` already accepted, `Δ = s - s_t` (≥ 0) satisfies
//!   `Δ <= q_max`, `p_t² + Δ² <= p_max²` and `p² + Δ² <= p_max²`
//!   (DISPLATH's `SimultaneousCriteria`: sites at similar path distance and
//!   inside a sphere of radius `p_max` around the collision point). DISPLATH
//!   sets `q_max` to the sum of the two atomic radii (about the
//!   nearest-neighbour distance); here the default is a small fraction of it,
//!   see the next section.
//! * **Collision.** Each partner `i` is scattered from the incoming state at
//!   its own impact parameter `p_i` (the distance from the site to the line)
//!   with the scattering table: transfer `T_i = γ E sin²(Θ_i/2)`, the recoil
//!   leaving along `sin(Θ_i/2) d + cos(Θ_i/2) n_i`, with `n_i` the unit vector
//!   from the line to the site (the half-angle form of
//!   [`kinematics::lab_angles`]: the recoil angle is `(π - Θ)/2`). The
//!   projectile takes the momentum that is left, `P = P0 - Σ sqrt(2 m_i T_i)
//!   r_i`, and the energies are scaled by the common factor
//!   `λ = (E - Σ Q_i) / (|P|²/2m + Σ T_i)` so that energy is conserved with
//!   the local electronic losses `Q_i` (DISPLATH's `Collision_!`). For one
//!   partner and no local loss this is the ordinary binary collision
//!   (`λ = 1`). The projectile stays at the closest-approach point of the
//!   nearest partner; the time-delay offsets of Robinson and Torrens (1974)
//!   are not applied (DISPLATH does not apply them either).
//! * **Recoils** start at their lattice site, with the usual displacement
//!   criterion ([`BcaConfig::follow_recoils`](super::BcaConfig::follow_recoils)) and follow the same model. The
//!   sites that were just hit are excluded from the next search (so is the
//!   site the recoil starts on): the ion cannot collide twice in a row with
//!   the same atom.
//!
//! # Choice of the search parameters
//!
//! Measured on this engine (B 5 keV into Si, 4000 ions, ZBL, Lindhard-Scharff,
//! 7 degrees tilt, 22 degrees twist, static lattice; ratios of Rp to the
//! amorphous result), with `nn` the nearest-neighbour distance (2.35 Å):
//!
//! * With `q_max = nn` (DISPLATH's rule) the result depends strongly on
//!   `p_max`: `p_max` = 1.53 Å gives Rp ratio 1.11, but `p_max` = 2.3 Å gives
//!   4.4. The merged collisions of distant sites cancel the kicks of a
//!   channel's walls, so almost every ion channels, at 20 and 30 degrees tilt
//!   as well. That is an artefact of merging over a sphere as large as the
//!   channel, not physics.
//! * With `q_max <= 0.1 nn` the result is independent of `q_max` and, for
//!   `p_max >= 1.5 * 1.53 Å = nn`, of `p_max` within the statistics (Rp ratio
//!   1.08-1.09 for `p_max` = 2.3, 3.1 Å and `q_max` = 0.01-0.1 nn; at
//!   `q_max` = 0.2 nn the artefact starts).
//!
//! So the defaults are `p_max = nn` and `q_max = 0.05 nn`: converged in both,
//! and the cost grows with `p_max`. These are properties of this engine's
//! implementation of the criterion, checked on silicon only; other lattices
//! need the same convergence check.
//!
//! # Measured checks (static lattice, 4000 ions, seed 1, Si, ZBL, Lindhard-Scharff)
//!
//! `tests/crystal_bca.rs` (the statistical ones are `#[ignore]`d; release
//! build). Ratios are crystal over amorphous:
//!
//! | Case | Rp | dRp | Notes |
//! |---|---|---|---|
//! | B 5 keV, 30° tilt, 17° twist | 0.97 | 1.09 | within 10 % / 15 % |
//! | As 30 keV, 30° tilt, 17° twist | 0.93 | 1.13 | within 10 % / 15 % |
//! | B 5 keV, 7° tilt, 22° twist | 1.06 | 1.24 | median 23.8 vs 23.8 nm, 99th percentile 76.2 vs 56.8 nm |
//! | As 30 keV, 7° tilt, 22° twist | 1.01 | 1.88 | median 22.2 vs 25.2 nm, 99th percentile 112.2 vs 50.8 nm |
//! | B 5 keV along <110> | **6.1** (121 vs 20 nm) | 3.2 | 91 % of ions deeper than twice the amorphous Rp (amorphous 5 %) |
//!
//! The random-direction check (Rp within 10 %, dRp within 15 % of the
//! amorphous result) is made at 30°/17° only. That beam is at least 8.0° from
//! every {100}, {110} and {111} plane. It passes within 1.3° of a {211} and
//! 0.9° of a {311} plane, which are higher-index planes with a lower atomic
//! areal density. The 7°/22° rows are not a random-direction check; see "The
//! 7°/22° criterion" below.
//!
//! # Measured checks with thermal vibration (same set-up, `THETA_D_SI`, zero-point term on)
//!
//! `tests/crystal_thermal.rs` (statistical ones `#[ignore]`d; release build).
//! "Tail" is the fraction of primaries deeper than twice the amorphous Rp
//! (amorphous: 0.034 at 7°/22°, 0.041 at 30°/17°), with its binomial error.
//!
//! | Case | Static | 0 K | 300 K | 600 K |
//! |---|---|---|---|---|
//! | B 5 keV along <110>: Rp, nm | 121.0 | 111.7 | 94.3 | 75.8 |
//! | B 5 keV along <110>: tail | 0.910 | 0.876 ± 0.005 | 0.833 ± 0.006 | 0.758 ± 0.007 |
//! | B 5 keV, 30°/17°: Rp / dRp ratio | 0.97 / 1.09 | 0.99 / 1.15 | 1.00 / 1.18 | 0.99 / 1.15 |
//! | As 30 keV, 30°/17°: Rp / dRp ratio | 0.93 / 1.13 | | 0.95 / 1.17 | |
//! | B 5 keV, 7°/22°: Rp / dRp ratio | 1.06 / 1.24 | 1.13 / 1.41 | 1.14 / 1.43 | 1.10 / 1.29 |
//! | As 30 keV, 7°/22°: Rp / dRp ratio | 1.01 / 1.88 | | 1.05 / 1.80 | |
//!
//! The static and 300 K entries of the 30°/17° and 7°/22° rows were measured
//! again for #225 (aarch64 macOS, commit `8513cfe`). The 0 K and 600 K
//! entries and the <110> rows are from PR #224. The crystal path is not bit
//! identical across platforms (`docs/architecture.md`). The re-measured
//! values differ from those of PR #224 by at most 0.02 in a dRp ratio. That
//! is within the seed-to-seed spread below.
//!
//! Along <110> the channeled tail falls with temperature, each step by more
//! than five binomial standard errors (the #181 criterion asks for a
//! monotonic decrease beyond the statistics). At the random direction the Rp
//! ratio stays within 10 % at 300 K.
//!
//! **Off-axis directions get a longer tail with vibration, not a shorter
//! one.** At 7°/22° (B) and 30°/17°, the tail beyond twice the amorphous Rp
//! grows from the static value as the amplitude rises from zero. It peaks
//! near `u1` ≈ 0.045-0.065 Å (0-300 K for Si) and falls again at higher
//! amplitudes. PR #224 measured, for B at 7°/22° with 4000 ions: tail 0.080
//! static, 0.075 at `u1` ≈ 0.0004 Å, 0.108 at 0 K, 0.115 at 300 K, 0.094 at
//! 600 K, and 0.076 at 1500 K. The dependence is smooth in the amplitude. It
//! is converged in the search parameters at 300 K: `q_max` = 0.01, 0.05 and
//! 0.1 nn and `p_max` = 1.5 nn spread the dRp ratio over 1.39-1.45 and change
//! Rp by under 1 %. So it is a property of the model as built, not of the
//! thermal code path.
//!
//! #225 measured how stable this is. The rise is significant at 16000 ions
//! (B 5 keV), under both electronic-loss treatments. The amorphous reference
//! uses the same loss as the crystal run
//! (`off_axis_tail_grows_with_vibration_under_both_losses`):
//!
//! | B 5 keV, 16000 ions | Amorphous tail | Static tail | 300 K tail | Rise | 300 K Rp / dRp ratio |
//! |---|---|---|---|---|---|
//! | 7°/22°, `NonLocal` | 0.0350 | 0.0806 ± 0.0022 | 0.1089 ± 0.0025 | 8.5 σ | 1.13 / 1.40 |
//! | 7°/22°, `EquipartitionLsOr` | 0.0395 | 0.0787 ± 0.0022 | 0.1059 ± 0.0025 | 8.2 σ | 1.12 / 1.39 |
//! | 30°/17°, `NonLocal` | 0.0452 | 0.0565 ± 0.0019 | 0.0714 ± 0.0021 | 5.2 σ | 1.00 / 1.20 |
//! | 30°/17°, `EquipartitionLsOr` | 0.0484 | 0.0523 ± 0.0018 | 0.0683 ± 0.0021 | 5.8 σ | 0.99 / 1.18 |
//!
//! With 4000 ions under `EquipartitionLsOr`, As 30 keV gives a tail of 0.057
//! static and 0.069 at 300 K at 7°/22° (amorphous 0.009). At 30°/17° it gives
//! 0.020 static and 0.028 at 300 K (amorphous 0.015).
//!
//! Over seeds 1 to 5 (4000 ions, `NonLocal`, 7°/22°), the tail is larger at
//! 300 K than static for every seed and both ions:
//!
//! * B: 0.075-0.085 static, 0.105-0.118 at 300 K.
//! * As: 0.057-0.060 static, 0.065-0.072 at 300 K.
//!
//! The dRp ratios have a seed-to-seed standard deviation of:
//!
//! * B: 0.018 static, 0.037 at 300 K.
//! * As: 0.052 static, 0.047 at 300 K.
//!
//! For As the dRp ratio does not rise with vibration (1.88 static, 1.80 at
//! 300 K). More ions reach the tail, but the deepest ones stop shallower:
//! the 99th percentile is 112 nm static and 98 nm at 300 K.
//!
//! The local Oen-Robinson half does not remove the rise. So the
//! amorphous-average nonlocal loss in channels is not its cause (#225, part
//! 3). Correlated vibration is absent: neighbouring atoms are displaced
//! independently ("Thermal vibration"). No source read for #225 gives the
//! sign of the effect of correlations on this tail, so its expected
//! direction is not stated.
//!
//! # Comparison with the literature (#225)
//!
//! Searched: measured depth profiles of B or As implanted into crystalline
//! (100) Si near 7° tilt at room temperature, and measurements against tilt,
//! twist or wafer temperature.
//!
//! **Not opened.** The leads from the issue could not be read:
//!
//! * Miyake et al., J. Electrochem. Soc. 130, 716 (1983),
//!   doi:10.1149/1.2119789. Paywalled at the publisher.
//! * Tian et al., Mikrochim. Acta (1992), doi:10.1007/BF01244469. Closed.
//! * Klein et al., IEEE Trans. Electron Devices 39 (1992),
//!   doi:10.1109/16.141226. Closed. It is a UT-MARLOWE paper, Tier C.
//!
//! Several further candidates were not readable either:
//!
//! * Semicond. Sci. Technol. 5 (1990), doi:10.1088/0268-1242/5/10/001,
//!   "Boron implants in <100> silicon at tilt angles of 0 degrees and 7
//!   degrees". The publisher served a bot check.
//! * J. Electrochem. Soc. papers on B channeling against tilt angle:
//!   doi:10.1149/1.2108785 and doi:10.1149/1.2085934. The publisher served a
//!   bot check.
//! * K. Nordlund, F. Djurabekova and G. Hobler, Phys. Rev. B 94, 214109
//!   (2016), doi:10.1103/PhysRevB.94.214109. The publisher refused access and
//!   the repository copy did not respond.
//!
//! None of these is cited for content here.
//!
//! **Read: measured profiles, not usable at matched conditions.** D. Cai,
//! N. Grønbech-Jensen, C. M. Snell and K. M. Beardmore, Phys. Rev. B 54,
//! 17147 (1996), arXiv:physics/9901056 (read in full). The paper shows SIMS
//! profiles of B in (100) Si:
//!
//! * at 15, 35 and 80 keV, at tilt 7° and rotation 30° (Fig. 3);
//! * at 5 keV, at tilt 7° and rotation 7° (Fig. 4; the SIMS data are from
//!   K. Gärtner, M. Nitschke and W. Eckstein, Nucl. Instrum. Methods B 83, 87
//!   (1993), which was not found openly readable).
//!
//! It also shows As at 8°/30° (Fig. 6). These profiles were not digitized,
//! for three reasons:
//!
//! 1. The paper does not give the in-plane reference direction of
//!    "rotation", so the orientation cannot be mapped through
//!    `docs/crystal-orientation.md`.
//! 2. It does not give the beam divergence, and it does not give the oxide
//!    or dose of the 5 keV measurement. Its 16 Å native oxide is a
//!    simulation setting.
//! 3. Each measured line overlaps the UT-MARLOWE curves (Tier C) in the same
//!    vector figure.
//!
//! A single room-temperature profile would also only show the model at
//! 300 K. It could not tell the static and the vibrating lattice apart.
//!
//! **Read: a published simulation of the same question.** M. I. Bratchenko,
//! A. S. Bakai and S. V. Dyuldya, J. Phys. Stud. 13, 1601 (2009),
//! doi:10.30970/jps.13.1601 (read in full). This is molecular dynamics with
//! the authors' own code, not a measurement. It simulates 15 keV B and As
//! into (001) Si at 7° tilt and 30° rotation, for a static lattice and at
//! 300 K, with uncorrelated Debye displacements. The paper reports:
//!
//! * In the static lattice, the deeply channeled states are almost empty.
//!   The authors attribute this to blocking, which forbids capture into
//!   them by a single strong collision.
//! * At 300 K "thermal vibrations facilitate the volume capture of ions into
//!   the stable channeling mode", that is, feeding-in.
//! * For As this lengthens the channeling tail substantially compared with
//!   the static lattice. For B the temperature effect on the tail length is
//!   "much weaker".
//! * "As long as a target is ordered there is no chance to obtain truly
//!   random-equivalent doping profiles".
//!
//! Its introduction also states, citing SIMS work not read here, that
//! off-axis SIMS profiles show long tails. It says the fraction of ions
//! beyond the amorphous Gaussian reaches about 20-30 % below 20 keV.
//!
//! The sign of the population effect here agrees with that paper: vibration
//! feeds ions into channels at an off-axis direction. But the energy (5 and
//! 30 keV here, against 15 keV) and the orientation (7°/22° here, against
//! 7°/30°) differ, and so does the measure. The paper discusses the tail
//! length, and for As at 30 keV our deepest percentile shortens. The paper
//! uses the same uncorrelated vibration model, so it does not test that
//! assumption. **The off-axis rise is therefore consistent with one
//! published simulation and not validated against measurement.** Validating
//! it needs a measured profile at matched, fully stated conditions, ideally
//! at two wafer temperatures (`docs/validation.md`, "Channeling (M2)").
//!
//! # The 7°/22° criterion (#225)
//!
//! #180 asked that dRp at 7°/22° be within 15 % of amorphous. That bound is
//! withdrawn. The beam lies 2.6° from a {100} plane and 2.7° from a {110}
//! plane (7.0° from <100>), so it is a near-planar direction, not a random
//! one. A tail of a few per cent of channeled ions is expected there, static
//! or vibrating, and the second moment weights it heavily. Without a
//! matched measurement, the 7°/22° runs are recorded-value regression checks
//! (`near_planar_7_22_static_lattice` in `tests/crystal_bca.rs` and
//! `near_planar_7_22_at_300_k` in `tests/crystal_thermal.rs`):
//!
//! * The dRp ratio stays within 0.12 (B) or 0.16 (As) of the table value.
//!   That is three seed-to-seed standard deviations, rounded up.
//! * The 90th percentile stays within 30 % of the amorphous one.
//! * Static lattice: Rp stays within 10 %.
//! * At 300 K: the tail exceeds the amorphous tail by more than five
//!   binomial standard errors.
//!
//! These are not agreement with experiment. They catch unintended changes of
//! the model.
//!
//! At 30°/17° the static lattice keeps both random-direction bounds. At
//! 300 K the dRp ratio, 1.18 (B) and 1.17 (As), is above the 15 % bound.
//! This is accepted as a property of the vibrating model, for three reasons:
//!
//! * The bulk of the profile still matches: Rp and the 90th percentile are
//!   within 10 % of amorphous.
//! * The excess is the thermally fed tail measured above, at 5 σ.
//! * The one published study read finds that thermal vibration causes this
//!   kind of feeding-in.
//!
//! The 300 K check asserts Rp and the 90th percentile within 10 %, and dRp
//! below its recorded value plus 0.12 (`random_direction_rp_at_300_k_matches_amorphous`).
//!
//! Unless marked `EquipartitionLsOr`, the ranges above use the nonlocal
//! Lindhard-Scharff loss only
//! ([`ElectronicLoss::NonLocal`](super::ElectronicLoss::NonLocal), the
//! default), so a channeled ion loses as much electronic energy per unit
//! path as in the amorphous target. With the local Oen-Robinson half (next
//! section) a channeled ion loses less per unit path than one in a random
//! direction. The off-axis tails under that mode are in the 16000-ion table
//! above.
//!
//! # Local electronic loss (Oen-Robinson)
//!
//! Under
//! [`ElectronicLoss::EquipartitionLsOr`](super::ElectronicLoss::EquipartitionLsOr)
//! half of the Lindhard-Scharff loss is nonlocal, taken along every flight
//! segment (above). The other half is local: every partner of a collision
//! step, so every lattice site that the search accepts within `p_max` of the
//! path, takes the Oen-Robinson loss
//! ([`OenRobinson::local_loss`](crate::ion::stopping::oen_robinson::OenRobinson::local_loss))
//! at its own distance of closest approach `r_min`. That distance comes from
//! the partner's impact parameter, which with thermal vibration is the
//! distance to the displaced position. The energy is the energy before the
//! collision step. The losses `Q_i` enter the energy-conserving scaling above.
//! O. S. Oen and M. T. Robinson, Nucl. Instrum. Methods 132, 647 (1976).
//!
//! What is not modelled:
//!
//! * **The `p_max` cutoff on the local half.** Sites farther than `p_max`
//!   from the path are not partners, so they take no local loss. The local
//!   half is therefore below half of Lindhard-Scharff wherever `p_max` is not
//!   large compared with the decay length `a / 0.3` of the Oen-Robinson loss.
//!   The amorphous model has the same cutoff. In a channel, where the ion
//!   stays farther than `p_max` from most rows, this cutoff is what lowers
//!   the loss.
//! * **Vibration beyond the displaced positions.** The local loss depends on
//!   the vibration amplitude only through the displaced site positions; there
//!   is no other thermal term.
//! * **Verified constants.** The Oen-Robinson constants are not verified
//!   against the paper
//!   ([`OR_CONSTANTS_UNVERIFIED`]).
//!   [`CrystalMetadata::electronic_constants_unverified`] records this in
//!   the run metadata.
//!
//! Measured by `tests/crystal_electronic.rs` (`#[ignore]`d; release build):
//! `R = E_local / E_nonlocal` of the primary (recoils not followed), counted
//! until its energy first falls below 80 % of the beam energy. Ar 20 keV into
//! Si, static lattice unless noted, seed 1, the crystal's `p_max` set to the
//! amorphous constant-path radius 1.532 Å. Errors are delta-method standard
//! errors.
//!
//! | Case | Histories | R | Notes |
//! |---|---|---|---|
//! | amorphous | 4000 | 0.8645 ± 0.0026 | the same at 30°/17° and 7°/22° |
//! | 30°/17° | 4000 | 0.9304 ± 0.0025 | 1.076 ± 0.004 × amorphous |
//! | 30°/17°, 300 K | 4000 | 0.9270 ± 0.0025 | 1.072 ± 0.004 × amorphous |
//! | 7°/22° | 4000 | 0.8744 ± 0.0030 | 1.011 ± 0.005 × amorphous |
//! | 30°/17°, `p_max` = nn | 4000 | 1.0033 ± 0.0026 | |
//! | along <110> | 8000 | 0.2448 ± 0.0022 | 0.263 × the 30°/17° crystal value (0.9316, same n) |
//! | along <100> | 8000 | 0.5600 ± 0.0010 | 0.601 × the 30°/17° crystal value |
//! | along <110> / <100>, 300 K | 8000 | 0.2548 / 0.5719 | 30°/17° at 300 K: 0.9283 |
//!
//! Along both axes the local share is far below the random direction (by
//! more than 150 standard errors), as channeling requires. **Gap (#250): in a
//! random direction the crystal ratio is 5-8 % above the amorphous one at the
//! same `p_max`**, outside the 5 % the validation asks for. The crystal
//! partners lean towards small impact parameters compared with the uniform
//! disc of the amorphous model. Vibration up to 3000 K and `q_max` hardly
//! change this. The cause is not established.
//!
//! # Thermal vibration (step 21b)
//!
//! With [`CrystalTarget::thermal`] set, each lattice site the particle meets
//! is displaced from its static position by a random vector `u` whose three
//! Cartesian components are independent Gaussians of standard deviation
//! `u1`, the one-dimensional RMS amplitude of the Debye model
//! ([`crate::ion::crystal::debye`]) for the site's species at the target
//! temperature [`Thermal::temperature_k`]. Displacements of different sites
//! are uncorrelated; a correlated mode is not implemented. The issue that
//! asked for this model (#181) names D. S. Gemmell, Rev. Mod. Phys. 46, 129
//! (1974), doi:10.1103/RevModPhys.46.129, as its reference; that paper could
//! not be opened, so nothing here is taken from it and none of its equations
//! is cited. The amplitude and its source are those of the `debye` module.
//!
//! * **Encounter.** A site's displacement is drawn when the site first
//!   appears in a search of the current *flight line* (a straight path of one
//!   particle in one crystal) and is kept while the line goes on: the next
//!   search segments of a channel see the same displaced atom. Any collision
//!   ends the line (and so does a change of direction at a surface or
//!   leaving the crystal), so every site met afterwards, by the projectile or
//!   by a recoil, gets a fresh draw. A site is therefore displaced once per
//!   encounter. Correlations in time (an atom met again later by the same
//!   cascade) are not kept.
//! * **What the displacement changes.** The impact parameter, the path
//!   distance (so the order of the partners and the simultaneity test), the
//!   direction `n_i` from the line to the atom, and so the deflection, the
//!   recoil direction and (under `EquipartitionLsOr`) the local electronic
//!   loss. A site's **region** (and so whether it is a partner at
//!   all) and the **starting point of its recoil** are its static lattice
//!   position, so that recoils, vacancies and tallies stay on lattice sites
//!   and a recoil's own site is recognised as such; the difference is
//!   `u1`, 0.065 Å for Si at 300 K (`THETA_D_SI`), against 2.35 Å between
//!   neighbours.
//! * **Search margin.** The static lattice is searched with the radius
//!   `p_max + m` over the path range `[-m, L + q_max + m]`, with
//!   `m =` [`THERMAL_MARGIN_U1`]` * u1` (largest `u1` of the crystal), and the
//!   static selection rules (`0 <= s <= L + q_max`, `p <= p_max`) are then
//!   applied to the displaced positions. A site whose displacement
//!   perpendicular to the path exceeds `m` can be missed: the probability is
//!   `exp(-m^2 / (2 u1^2)) = exp(-18)`, about `1.5e-8`, per site.
//! * **Random numbers.** The displacements come from the history's
//!   counter-based stream ([`crate::rng::stream`]), from a copy positioned at
//!   word 2^66, a segment that neither the transport draws nor the lattice
//!   translation (word 2^67) reach. Each displacement consumes exactly six
//!   words ([`sample_gaussian_3d`]), drawn in the documented order of the
//!   search results ([`crate::ion::crystal::search`], "Ordering"), and the
//!   primary and its recoils share the stream in the order they are followed.
//!   So a history's result depends only on `(seed, index)`, never on the
//!   thread count, and turning vibration on changes no other draw.
//! * **Static limit.** With every amplitude zero (`T = 0` and
//!   [`Thermal::include_zero_point`] `= false`, a test switch) the thermal
//!   code path gives the static model bit for bit
//!   (`tests/crystal_thermal.rs`). With the zero-point term (the physical
//!   case) `T = 0` still vibrates: `u1` = 0.045 Å for Si (0.088 Å at 600 K).
//! * **Not changed.** The lattice constant is the cited one at its own
//!   temperature (no thermal expansion). The electronic-loss model is the
//!   one of the static lattice: the nonlocal loss along a segment, and under
//!   `EquipartitionLsOr` the local loss at each partner's closest approach.
//!   The local loss sees the vibration only through the displaced position
//!   (the impact parameter, so `r_min`); it has no other thermal term.
//!
//! DISPLATH (Tier A, above; same commit) also perturbs lattice atoms by
//! independent Gaussian displacements of a Debye amplitude, but draws them
//! once per cell and cascade (`LatticeSiteCoordinatesBase!` and
//! `TemperatureToSigma` in `src/geometry.jl`); the per-encounter scheme here
//! is this crate's, as #181 asks, and no DISPLATH code was ported for it.
//!
//! # What is not here
//!
//! The lattice is **perfect** (static, or vibrating as above) and
//! stays perfect: a displaced atom leaves no vacancy, so a later particle can
//! still collide with its site (dynamic damage, dechanneling and
//! amorphization are later steps). The impact-parameter-dependent electronic
//! loss is the local Oen-Robinson half, taken only at sites within `p_max`
//! (see "Local electronic loss" for what that leaves out). Boundaries of a
//! periodic voxel axis do not continue the lattice. The equation-level treatment of Robinson and
//! Torrens, Phys. Rev. B 9, 5008 (1974) (MARLOWE, Tier C) was not available
//! to read when this was written; the choices above are DISPLATH's, not
//! quotations of that paper.
//!
//! [`Geometry`]: crate::geometry::Geometry
//! [`LatticeSearch`]: crate::ion::crystal::LatticeSearch

use crate::geometry::Flight;
use crate::ion::crystal::debye::{sample_gaussian_3d, ThermalVibration};
use crate::ion::crystal::search::{compare, path_metrics, unit_direction};
use crate::ion::crystal::{Lattice, LatticeSearch, Orientation};
use crate::ion::scattering::closest_approach;
use crate::ion::stopping::oen_robinson::OR_CONSTANTS_UNVERIFIED;
use crate::ion::stopping::Ion;
use crate::rng::ParticleRng;
use crate::units::J_PER_EV;

use super::{
    kinematics, Bca, BcaError, BcaTally, ElectronicChannel, EnergyBudget, MeanFreePath, Particle,
    Scratch,
};
pub(super) use crate::ion::crystal::Candidate;

/// A lattice site: its conventional cell and index in the cell.
pub(super) type SiteId = ([i64; 3], u8);

/// Search margin of the thermal model, in units of the largest
/// one-dimensional RMS displacement `u1` of the crystal's species (module
/// docs, "Thermal vibration").
pub const THERMAL_MARGIN_U1: f64 = 6.0;

/// How [`Thermal`] displacements are sampled, as recorded in
/// [`ThermalMetadata::sampling`].
pub const THERMAL_SAMPLING: &str = "per-encounter, uncorrelated: each lattice site met along a \
     straight flight line gets one isotropic Gaussian displacement (standard deviation u1 per \
     Cartesian axis, Debye model) drawn from the history's thermal stream; redrawn after every \
     collision";

/// Thermal vibration of the lattice atoms of a [`CrystalTarget`] (step 21b;
/// module docs, "Thermal vibration").
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Thermal {
    /// Target temperature, K (finite, `>= 0`). This sets the vibration
    /// amplitude only: the lattice constant stays the cited value at its own
    /// temperature ([`Lattice::temperature_k`]); no thermal expansion is
    /// applied.
    pub temperature_k: f64,
    /// Debye temperature, K, used for every species of the crystal with that
    /// species' own mass ([`crate::ion::crystal::debye`] module docs; the
    /// `THETA_D_*` constants there are cited defaults).
    pub debye_temperature_k: f64,
    /// Include the zero-point term of the Debye amplitude (`true` for any
    /// physical run). `false` keeps only the thermal part, so that `T = 0`
    /// gives zero amplitude: a **test switch** (the static-lattice identity
    /// test uses it), not a physical model.
    pub include_zero_point: bool,
}

impl Thermal {
    /// Vibration at `temperature_k` with Debye temperature
    /// `debye_temperature_k`, zero-point term included.
    pub fn new(temperature_k: f64, debye_temperature_k: f64) -> Self {
        Self {
            temperature_k,
            debye_temperature_k,
            include_zero_point: true,
        }
    }
}

/// Run metadata of one crystal of a [`Bca`] ([`Bca::crystal_metadata`]): what
/// a result needs to be reproduced from its own header.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct CrystalMetadata {
    /// Regions filled by this crystal.
    pub regions: Vec<usize>,
    /// Cubic lattice constant, m.
    pub lattice_constant_m: f64,
    /// Temperature of the cited lattice constant, K (not the target
    /// temperature; see [`Thermal::temperature_k`]).
    pub lattice_constant_temperature_k: f64,
    /// Miller indices of the surface normal and the in-plane reference.
    pub normal_hkl: [i32; 3],
    /// Indices of the in-plane reference direction.
    pub reference_uvw: [i32; 3],
    /// Tilt, twist and wafer rotation, rad.
    pub tilt_rad: f64,
    /// Twist, rad.
    pub twist_rad: f64,
    /// Wafer rotation, rad.
    pub wafer_rotation_rad: f64,
    /// Search parameters in use, m.
    pub p_max_m: f64,
    /// Simultaneous-collision window, m.
    pub q_max_m: f64,
    /// Search segment length, m.
    pub search_length_m: f64,
    /// Thermal vibration; `None` for a static lattice.
    pub thermal: Option<ThermalMetadata>,
    /// `true` when the run takes the local Oen-Robinson loss
    /// ([`ElectronicLoss::EquipartitionLsOr`](super::ElectronicLoss::EquipartitionLsOr))
    /// and the constants of that model are not verified against the paper
    /// ([`OR_CONSTANTS_UNVERIFIED`]). Channeled ranges computed with the local
    /// loss depend on those constants. `false` for a nonlocal-only run; a
    /// `false` value is left out of the serialised form.
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub electronic_constants_unverified: bool,
}

/// The thermal part of [`CrystalMetadata`].
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct ThermalMetadata {
    /// The input.
    pub input: Thermal,
    /// One-dimensional RMS displacement `u1` per species, as
    /// `(Z, mass in u, u1 in m)`, in order of first appearance in the
    /// lattice basis.
    pub rms_1d_m: Vec<(u8, f64, f64)>,
    /// Extra search radius and path margin, m ([`THERMAL_MARGIN_U1`] times
    /// the largest `u1`).
    pub search_margin_m: f64,
    /// The sampling rule ([`THERMAL_SAMPLING`]).
    pub sampling: &'static str,
}

/// Precomputed thermal data of one crystal.
#[derive(Debug, Clone)]
struct ThermalData {
    meta: ThermalMetadata,
    /// `u1` per conventional-cell site (index = [`Candidate::basis`]), m.
    sigma: Vec<f64>,
    margin: f64,
}

/// A site closer than this fraction of the lattice constant to the particle
/// is the particle's own starting site, not a partner.
const SELF_REL: f64 = 1e-6;

/// Default `q_max`, as a fraction of the nearest-neighbour distance (module
/// docs, "Choice of the search parameters").
pub const DEFAULT_Q_MAX_FRACTION: f64 = 0.05;

/// Largest relative density mismatch between lattice and material.
const DENSITY_TOLERANCE: f64 = 0.05;

/// A crystal for [`Bca::with_crystal`]: the lattice, the orientation of the
/// target in the lab frame and the three search parameters.
#[derive(Debug, Clone)]
pub struct CrystalTarget {
    /// The lattice (cubic; diamond or zincblende).
    pub lattice: Lattice,
    /// Wafer cut and angles. Only its lab-to-crystal rotation is used.
    pub orientation: Orientation,
    /// Largest impact parameter of a partner, m. Default `None`: the
    /// nearest-neighbour distance of the lattice. Steering by the atoms of a
    /// channel wall needs partners out to about that distance; smaller radii
    /// (DISPLATH's examples use `a/3` and the amorphous constant-path radius
    /// `(π N^(2/3))^(-1/2)`, `a/2/π^0.5` for silicon) leave the channeled
    /// fraction too small. See the module docs, "Choice of the search
    /// parameters".
    pub p_max_m: Option<f64>,
    /// Largest path distance between simultaneous partners, m. Default
    /// `None`: [`DEFAULT_Q_MAX_FRACTION`] of the nearest-neighbour distance.
    pub q_max_m: Option<f64>,
    /// Length of one search segment, m. Default `None`: the lattice
    /// constant. Affects speed only (a segment with no partner is followed by
    /// another), never the result.
    pub search_length_m: Option<f64>,
    /// Thermal vibration of the lattice atoms. Default `None`: a static,
    /// perfect lattice (the 20b model, unchanged bit for bit).
    pub thermal: Option<Thermal>,
}

impl CrystalTarget {
    /// A static crystal with the default search parameters.
    pub fn new(lattice: Lattice, orientation: Orientation) -> Self {
        Self {
            lattice,
            orientation,
            p_max_m: None,
            q_max_m: None,
            search_length_m: None,
            thermal: None,
        }
    }

    /// The same crystal with thermal vibration `thermal`.
    pub fn with_thermal(mut self, thermal: Thermal) -> Self {
        self.thermal = Some(thermal);
        self
    }
}

/// Validated, precomputed crystal data.
#[derive(Debug, Clone)]
pub(super) struct CrystalData {
    search: LatticeSearch,
    orientation: Orientation,
    a: f64,
    a_temperature_k: f64,
    p_max: f64,
    q_max: f64,
    search_len: f64,
    thermal: Option<ThermalData>,
}

impl CrystalData {
    pub(super) fn lattice_constant(&self) -> f64 {
        self.a
    }

    /// Whether this crystal's sites vibrate (a thermal stream is needed).
    pub(super) fn is_thermal(&self) -> bool {
        self.thermal.is_some()
    }

    /// Crystal-frame position of a lab point for the history translation
    /// `shift`.
    fn to_crystal(&self, lab: [f64; 3], shift: [f64; 3]) -> [f64; 3] {
        let c = self.orientation.to_crystal(lab);
        [c[0] + shift[0], c[1] + shift[1], c[2] + shift[2]]
    }

    /// Lab position of a crystal-frame point.
    fn to_lab(&self, c: [f64; 3], shift: [f64; 3]) -> [f64; 3] {
        self.orientation
            .to_lab([c[0] - shift[0], c[1] - shift[1], c[2] - shift[2]])
    }
}

/// One partner of the current collision step.
#[derive(Debug, Clone, Copy)]
pub(super) struct Partner {
    id: SiteId,
    /// Path distance of the closest approach, m.
    s: f64,
    /// Impact parameter, m.
    b: f64,
    /// Unit vector (lab) from the line to the site.
    n: [f64; 3],
    /// Site position (lab), m.
    site: [f64; 3],
    /// Region of the site and element index in its material.
    region: usize,
    j: usize,
    /// Transfer, recoil direction and local loss of the pair collision.
    t: f64,
    rdir: [f64; 3],
    q: f64,
}

/// Nearest-neighbour distance of `lattice`, m.
fn nearest_neighbour(lattice: &Lattice) -> f64 {
    let a = lattice.lattice_constant();
    let sites = lattice.conventional_cell_sites();
    let mut best = f64::INFINITY;
    for (_, f) in &sites {
        for (_, g) in &sites {
            for i in -1..=1 {
                for j in -1..=1 {
                    for k in -1..=1 {
                        let dx = a * (g[0] + f64::from(i) - f[0]);
                        let dy = a * (g[1] + f64::from(j) - f[1]);
                        let dz = a * (g[2] + f64::from(k) - f[2]);
                        let d = (dx * dx + dy * dy + dz * dz).sqrt();
                        if d > 1e-6 * a {
                            best = best.min(d);
                        }
                    }
                }
            }
        }
    }
    best
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl<'a> Bca<'a> {
    /// Use the crystal flight model in `regions` (see the module docs). The
    /// regions are indices into the geometry's regions: layer indices of a
    /// stack, flat voxel indices of a voxel grid.
    ///
    /// # Errors
    /// [`BcaError::InvalidConfig`] if the configuration is not the constant
    /// free path without weak collisions, a region does not exist or already
    /// belongs to a crystal, a search parameter is not finite and positive, a
    /// region's material lacks a lattice element, or its atom density differs
    /// from the lattice's by more than 5 %.
    pub fn with_crystal(
        mut self,
        target: CrystalTarget,
        regions: &[usize],
    ) -> Result<Self, BcaError> {
        let bad = |m: String| BcaError::InvalidConfig(m);
        if !matches!(self.config.mean_free_path, MeanFreePath::Constant)
            || self.config.weak_collisions > 0
        {
            return Err(bad(
                "crystal regions need the constant free path and no weak collisions".into(),
            ));
        }
        if regions.is_empty() {
            return Err(bad("a crystal needs at least one region".into()));
        }
        let lat = &target.lattice;
        let a = lat.lattice_constant();
        let n_lat = lat.atom_number_density();
        let zs: Vec<u8> = {
            let mut v: Vec<u8> = lat.conventional_cell_sites().iter().map(|s| s.0).collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        if self.region_crystal.is_empty() {
            self.region_crystal = vec![None; self.geometry.n_regions()];
        }
        let ci = self.crystals.len();
        for &r in regions {
            if r >= self.geometry.n_regions() {
                return Err(bad(format!("crystal region {r} does not exist")));
            }
            if self.region_crystal[r].is_some() {
                return Err(bad(format!("region {r} already belongs to a crystal")));
            }
            let lay = self.lay(r);
            for &z in &zs {
                if !lay.elems.iter().any(|e| e.z == z) {
                    return Err(bad(format!(
                        "region {r}: the material has no element Z={z} of the lattice"
                    )));
                }
            }
            if (lay.n / n_lat - 1.0).abs() > DENSITY_TOLERANCE {
                return Err(bad(format!(
                    "region {r}: material atom density {:.4e} m^-3 differs from the \
                     lattice's {n_lat:.4e} m^-3 by more than 5 %",
                    lay.n
                )));
            }
        }
        let nn = nearest_neighbour(lat);
        let p_max = target.p_max_m.unwrap_or(nn);
        let q_max = target.q_max_m.unwrap_or(DEFAULT_Q_MAX_FRACTION * nn);
        let search_len = target.search_length_m.unwrap_or(a);
        for (name, v) in [
            ("p_max_m", p_max),
            ("q_max_m", q_max),
            ("search_length_m", search_len),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(bad(format!("{name} = {v} must be finite and positive")));
            }
        }
        let thermal = match target.thermal {
            None => None,
            Some(th) => {
                // One Debye model per species, with the mass the engine uses
                // for that species (the material's element in the first
                // region; species data are shared by all regions).
                let lay = self.lay(regions[0]);
                let mut rms: Vec<(u8, f64, f64)> = Vec::new();
                for (z, _) in lat.conventional_cell_sites() {
                    if rms.iter().any(|r| r.0 == z) {
                        continue;
                    }
                    let elem = lay
                        .elems
                        .iter()
                        .find(|e| e.z == z)
                        .expect("lattice elements were checked above");
                    let mass = self.species[elem.species].ion.mass_amu();
                    let v = ThermalVibration::new(th.debye_temperature_k, mass, th.temperature_k)
                        .map_err(|e| bad(format!("thermal vibration: {e}")))?;
                    let u1_sq = if th.include_zero_point {
                        v.mean_square_1d()
                    } else {
                        v.thermal_only_mean_square_1d()
                    };
                    rms.push((z, mass, u1_sq.sqrt()));
                }
                let sigma: Vec<f64> = lat
                    .conventional_cell_sites()
                    .iter()
                    .map(|(z, _)| rms.iter().find(|r| r.0 == *z).map_or(0.0, |r| r.2))
                    .collect();
                let u1_max = sigma.iter().copied().fold(0.0, f64::max);
                let margin = THERMAL_MARGIN_U1 * u1_max;
                Some(ThermalData {
                    meta: ThermalMetadata {
                        input: th,
                        rms_1d_m: rms,
                        search_margin_m: margin,
                        sampling: THERMAL_SAMPLING,
                    },
                    sigma,
                    margin,
                })
            }
        };
        for &r in regions {
            self.region_crystal[r] = Some(ci);
        }
        self.crystals.push(CrystalData {
            search: LatticeSearch::new(lat),
            orientation: target.orientation,
            a,
            a_temperature_k: lat.temperature_k(),
            p_max,
            q_max,
            search_len,
            thermal,
        });
        Ok(self)
    }

    /// Run metadata of every crystal ([`Bca::with_crystal`] order): lattice,
    /// orientation, search parameters and the thermal model with its target
    /// temperature. Empty for an amorphous engine.
    pub fn crystal_metadata(&self) -> Vec<CrystalMetadata> {
        let local_or = self.config.electronic == super::ElectronicLoss::EquipartitionLsOr;
        self.crystals
            .iter()
            .enumerate()
            .map(|(ci, cr)| CrystalMetadata {
                regions: (0..self.region_crystal.len())
                    .filter(|&r| self.region_crystal[r] == Some(ci))
                    .collect(),
                lattice_constant_m: cr.a,
                lattice_constant_temperature_k: cr.a_temperature_k,
                normal_hkl: cr.orientation.normal_hkl(),
                reference_uvw: cr.orientation.reference_uvw(),
                tilt_rad: cr.orientation.tilt_rad(),
                twist_rad: cr.orientation.twist_rad(),
                wafer_rotation_rad: cr.orientation.wafer_rotation_rad(),
                p_max_m: cr.p_max,
                q_max_m: cr.q_max,
                search_length_m: cr.search_len,
                thermal: cr.thermal.as_ref().map(|t| t.meta.clone()),
                electronic_constants_unverified: local_or && OR_CONSTANTS_UNVERIFIED,
            })
            .collect()
    }

    /// The crystal of `region`, if it is a crystal region.
    #[inline]
    pub(super) fn crystal_at(&self, region: usize) -> Option<usize> {
        self.region_crystal.get(region).copied().flatten()
    }

    /// Region and element index of a lattice site that may be a partner:
    /// inside the target, in a region of crystal `ci`.
    fn site_region(&self, ci: usize, site_lab: [f64; 3], z: u8) -> Option<(usize, usize)> {
        let r = self.geometry.locate(site_lab)?;
        if self.crystal_at(r) != Some(ci) {
            return None;
        }
        let j = self.lay(r).elems.iter().position(|e| e.z == z)?;
        Some((r, j))
    }

    /// One step of a particle in crystal `ci`: a flight to the next
    /// collision (or to a surface or interface, or one search segment) and
    /// the collision there. Returns `true` if the particle is finished.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn crystal_step<T: BcaTally>(
        &self,
        ci: usize,
        shift: [f64; 3],
        p: &mut Particle,
        pending: &mut Vec<Particle>,
        budget: &mut EnergyBudget,
        tally: &mut T,
        scratch: &mut Scratch,
        trng: &mut Option<ParticleRng>,
    ) -> Result<bool, crate::ion::stopping::StoppingError> {
        let cr = &self.crystals[ci];
        let cutoff = self.species[p.species].cutoff_ev;
        let Scratch {
            cands,
            last,
            targets,
            homes,
            encounter,
            displaced,
            displaced_next,
            displaced_line,
            ..
        } = scratch;
        let o = cr.to_crystal(p.pos, shift);
        let d = cr.orientation.to_crystal(p.dir);
        let self_r2 = (SELF_REL * cr.a) * (SELF_REL * cr.a);
        match &cr.thermal {
            None => cr
                .search
                .search_into(o, d, cr.p_max, cr.search_len + cr.q_max, cands)
                .expect("particle state is finite"),
            Some(th) => {
                // Thermal vibration (module docs): search the static lattice
                // with a margin, displace every site found (once per
                // encounter: a site keeps its displacement while the flight
                // line is unchanged), then apply the static selection rules
                // to the displaced positions.
                let m = th.margin;
                let dn = unit_direction(d);
                let ob = [o[0] - m * dn[0], o[1] - m * dn[1], o[2] - m * dn[2]];
                cr.search
                    .search_into(
                        ob,
                        d,
                        cr.p_max + m,
                        cr.search_len + cr.q_max + 2.0 * m,
                        cands,
                    )
                    .expect("particle state is finite");
                let line = (ci, [d[0].to_bits(), d[1].to_bits(), d[2].to_bits()]);
                if *displaced_line != Some(line) {
                    displaced.clear();
                    *displaced_line = Some(line);
                }
                let rng = trng
                    .as_mut()
                    .expect("a thermal crystal has a thermal stream");
                displaced_next.clear();
                encounter.clear();
                let (pm2, len) = (cr.p_max * cr.p_max, cr.search_len + cr.q_max);
                for c in cands.iter() {
                    // The particle's own starting site (a recoil starts on its
                    // static site) is never a partner.
                    let (s0, p0sq) = path_metrics(o, dn, c.position);
                    let p0 = p0sq.sqrt();
                    if s0 * s0 + p0 * p0 < self_r2 {
                        continue;
                    }
                    let id = (c.cell, c.basis);
                    let u = match displaced.get(&id) {
                        Some(&u) => u,
                        None => sample_gaussian_3d(th.sigma[usize::from(c.basis)], rng),
                    };
                    displaced_next.insert(id, u);
                    let r = [
                        c.position[0] + u[0],
                        c.position[1] + u[1],
                        c.position[2] + u[2],
                    ];
                    let (s, p2) = path_metrics(o, dn, r);
                    if s >= 0.0 && s <= len && p2 <= pm2 {
                        encounter.push((
                            Candidate {
                                position: r,
                                s,
                                p: p2.sqrt(),
                                ..*c
                            },
                            c.position,
                        ));
                    }
                }
                // Sites no longer in the window are behind the particle on
                // this line and can never be found again.
                std::mem::swap(displaced, displaced_next);
                encounter.sort_by(|x, y| compare(&x.0, &y.0));
                cands.clear();
                homes.clear();
                for &(c, h) in encounter.iter() {
                    cands.push(c);
                    homes.push(h);
                }
            }
        }
        let thermal = cr.thermal.is_some();
        // Static (lattice-site) position of candidate `i`, crystal frame: the
        // site's region and the recoil's starting point.
        let home = |i: usize| if thermal { homes[i] } else { cands[i].position };

        // The nearest valid partner within one search segment.
        let usable = |i: usize| -> Option<(usize, usize)> {
            let c = &cands[i];
            if c.s * c.s + c.p * c.p < self_r2 || last.contains(&(c.cell, c.basis)) {
                return None;
            }
            self.site_region(ci, cr.to_lab(home(i), shift), c.z)
        };
        let hit = (0..cands.len()).find(|&i| cands[i].s <= cr.search_len && usable(i).is_some());
        let limit = hit.map_or(cr.search_len, |i| cands[i].s);

        match self.geometry.flight(p.layer, p.pos, p.dir, limit) {
            Flight::Event(ex) => {
                let from = p.pos;
                p.pos = ex.at;
                self.electronic_nonlocal(p, from, ex.distance, budget, tally)?;
                if p.energy_ev < cutoff {
                    return Ok(false);
                }
                return Ok(self.apply_event(p, ex.outcome, budget, tally));
            }
            Flight::Clear { region } => {
                let from = p.pos;
                Self::advance(p, limit);
                p.layer = region;
                self.electronic_nonlocal(p, from, limit, budget, tally)?;
            }
        }
        let Some(hit) = hit else {
            return Ok(false);
        };
        if p.energy_ev < cutoff || self.crystal_at(p.layer) != Some(ci) {
            return Ok(false);
        }

        // Partners: the nearest, then every later site that is simultaneous
        // with all partners accepted so far.
        targets.clear();
        let s0 = cands[hit].s;
        let sp2 = cr.p_max * cr.p_max;
        for (i, c) in cands.iter().enumerate().skip(hit) {
            if c.s - s0 > cr.q_max {
                break;
            }
            let Some((region, j)) = usable(i) else {
                continue;
            };
            let simultaneous = targets.iter().all(|t: &Partner| {
                let dl = c.s - t.s;
                dl <= cr.q_max && t.b * t.b + dl * dl <= sp2 && c.p * c.p + dl * dl <= sp2
            });
            if !simultaneous {
                continue;
            }
            // Unit vector from the line to the site, in the lab frame.
            let q = [
                c.position[0] - o[0] - c.s * d[0],
                c.position[1] - o[1] - c.s * d[1],
                c.position[2] - o[2] - c.s * d[2],
            ];
            let qn = dot(q, q).sqrt();
            let n_c = if qn > 1e-9 * cr.a {
                [q[0] / qn, q[1] / qn, q[2] / qn]
            } else {
                // On the line: any perpendicular does (the same rule as the
                // azimuth of an amorphous collision).
                kinematics::rotate_sc(d, (1.0, 0.0), (1.0, 0.0))
            };
            targets.push(Partner {
                id: (c.cell, c.basis),
                s: c.s,
                b: c.p,
                n: kinematics::normalize(cr.orientation.to_lab(n_c)),
                site: cr.to_lab(home(i), shift),
                region,
                j,
                t: 0.0,
                rdir: [0.0; 3],
                q: 0.0,
            });
        }
        last.clear();
        last.extend(targets.iter().map(|t| t.id));
        // The collision ends this flight line: every site met after it is a
        // new encounter with a fresh displacement.
        displaced.clear();
        *displaced_line = None;
        self.crystal_collide(p, targets, pending, budget, tally)?;
        Ok(false)
    }

    /// The simultaneous collision of `p` with `targets` (module docs).
    fn crystal_collide<T: BcaTally>(
        &self,
        p: &mut Particle,
        targets: &mut [Partner],
        pending: &mut Vec<Particle>,
        budget: &mut EnergyBudget,
        tally: &mut T,
    ) -> Result<(), crate::ion::stopping::StoppingError> {
        let ns = self.species.len();
        let e0 = p.energy_ev;
        let d0 = p.dir;
        let m_p = p.mass_amu;
        let beta_min = self.table.spec().beta_min;
        let p0 = (2.0 * m_p * e0).sqrt();
        let mut mom = [p0 * d0[0], p0 * d0[1], p0 * d0[2]];
        let (mut sum_t, mut sum_q) = (0.0, 0.0);
        for t in targets.iter_mut() {
            let elem = &self.lay(t.region).elems[t.j];
            let pair = self.pairs[p.species * ns + elem.species];
            let eps = e0 * pair.eps_per_ev;
            let beta = (t.b / pair.a).max(beta_min);
            let tan_half = self.half_angle_tan(eps, beta);
            let c = 1.0 / (1.0 + tan_half * tan_half).sqrt();
            let s = tan_half * c;
            t.t = (pair.gamma * e0 * s * s).min(e0);
            // Recoil leaves at (pi - Theta)/2 from the incoming direction,
            // towards the site: cos = sin(Theta/2), sin = cos(Theta/2).
            t.rdir = [
                s * d0[0] + c * t.n[0],
                s * d0[1] + c * t.n[1],
                s * d0[2] + c * t.n[2],
            ];
            if self.config.electronic == super::ElectronicLoss::EquipartitionLsOr {
                let r_min = closest_approach(self.screening, eps, beta) * pair.a;
                let ion: &Ion = &self.species[p.species].ion;
                let q = self.mix.0.local_loss(ion, elem.z, e0, r_min)? / J_PER_EV;
                t.q = q.min(e0 - sum_q).max(0.0);
                sum_q += t.q;
            }
            let pt = (2.0 * self.species[elem.species].ion.mass_amu() * t.t).sqrt();
            for (m, r) in mom.iter_mut().zip(t.rdir) {
                *m -= pt * r;
            }
            sum_t += t.t;
        }
        // Energy conservation with the local losses: scale the transfers and
        // the projectile's remaining energy by a common factor.
        let pm2 = dot(mom, mom);
        let e_p = pm2 / (2.0 * m_p);
        let need = (e0 - sum_q).max(0.0);
        let lambda = if e_p + sum_t > 0.0 {
            need / (e_p + sum_t)
        } else {
            1.0
        };
        if pm2 > 0.0 {
            p.dir = kinematics::normalize(mom);
        }
        p.energy_ev = lambda * e_p;
        for t in targets.iter() {
            if t.q > 0.0 {
                budget.electronic_local += t.q;
                tally.electronic(p, p.pos, ElectronicChannel::Local, t.q);
            }
        }
        for t in targets.iter() {
            let rdir = kinematics::normalize(t.rdir);
            self.emit_recoil(
                p.generation,
                t.site,
                t.region,
                t.j,
                lambda * t.t,
                || rdir,
                pending,
                budget,
                tally,
            );
        }
        Ok(())
    }
}
