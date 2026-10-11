# %% [markdown]
# # An electron run from Python
#
# Runs `examples/electron/e_10keV_si.toml` (10 keV electrons into Si, with
# SYNTHETIC optical and band data; see the header of that file) through
# `ElectronRun.run()` and reads the escape spectra and the energy-deposition
# grid. Cells are separated by `# %%` markers. Run it from the repository
# root; the file's optical-ELF path resolves against the file's directory.

# %%
import lindhard as lh

run = lh.ElectronRun.from_toml_file("examples/electron/e_10keV_si.toml")
for w in run.validate():
    print("warning:", w)
# Releases the GIL while the tables are built and transport runs; same numbers
# at any thread count.
res = run.run(histories=200, seed=1)
print(res)

# %% [markdown]
# ## Escape spectra
#
# `res.escape_spectra["front"]` holds `Histogram`s: `energy` in eV and the
# polar angle of the slow (secondary) and fast (backscattered) electrons in
# degrees from the outward normal. `density` is per primary per unit.

# %%
energy = res.escape_spectra["front"]["energy"]
print("escaped electrons (front):", int(energy.counts.sum()))
print("peak of the spectrum at", energy.centers[energy.counts.argmax()], "eV")

# %% [markdown]
# ## Energy deposition
#
# `res.cylindrical` (with `[electron.tally.cylindrical]`; `None` otherwise) has
# `energy_ev` of shape `(n_r, n_depth)` summed over all histories and
# `ev_per_primary_per_nm3` normalised by histories and cell volume. Edges are
# in nm.

# %%
cyl = res.cylindrical
if cyl is not None:
    by_depth = cyl["energy_ev"].sum(axis=0) / res.histories
    peak = by_depth.argmax()
    print(
        "deposited per primary, peak depth bin:",
        cyl["depth_edges_nm"][peak],
        "-",
        cyl["depth_edges_nm"][peak + 1],
        "nm",
    )
y = res.summary()["results"]["yields"]
print("backscatter eta =", y["backscatter_eta"], " secondary delta =", y["secondary_delta"])
