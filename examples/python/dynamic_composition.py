# %% [markdown]
# # Fluence-dependent composition from Python
#
# Runs `examples/dynamic/as_1keV_si_film.toml` (1 keV As into a 10 nm Si film
# on Si) through `Run.run_dynamic()` and looks at the target as the fluence
# builds up. Cells are separated by `# %%` markers, so the file opens as a
# notebook in editors that understand them. Run it from the repository root.

# %%
import lindhard as lh

run = lh.Run.from_toml_file("examples/dynamic/as_1keV_si_film.toml")
res = run.run_dynamic(seed=1)  # releases the GIL; same numbers at any thread count
print(res)

# %% [markdown]
# ## Step series
#
# `res.steps` is a dict of arrays with one entry per step; entry 0 is the
# initial target. Fluence is in ions/cm^2, lengths in nm.

# %%
steps = res.steps
for i in range(res.n_steps):
    print(
        f"step {steps['step'][i]:3d}  fluence {steps['fluence_cm2'][i]:.3e} cm^-2  "
        f"slabs {steps['n_slabs'][i]}  sputter yield {steps['sputter_yield'][i]:.3f}"
    )

# %% [markdown]
# ## Composition snapshots
#
# `res.composition[i]` holds the slab profile after step `i`. The number of
# slabs can change between steps (with `erosion = true` emptied slabs are
# removed), so the snapshots are a list, not one rectangular array. Columns of
# `fraction` follow `res.species`.

# %%
first, last = res.composition[0], res.composition[-1]
print("species:", res.species)
print("initial slab fractions:\n", first["fraction"])
print("final slab fractions:\n", last["fraction"])
print("final slab centres (nm):", 0.5 * (last["front_nm"] + last["back_nm"]))

# %% [markdown]
# ## Writing the CLI's files
#
# `write()` produces the same `dynamic_summary.json`, `dynamic_steps.csv` and
# `dynamic_composition.csv` as `lindhard run`, publishing the summary last.

# %%
res.write("out/as_si_dynamic_py")
print(res.summary()["totals"])
