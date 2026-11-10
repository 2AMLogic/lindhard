"""Python bindings for lindhard.

Configuration classes (:class:`Material`, :class:`Layer`, :class:`Target`,
:class:`Beam`, :class:`Physics`, :class:`Tally`) are assembled into a
:class:`Run`, which reads and writes the same TOML schema as the ``lindhard``
command (:meth:`Run.from_toml`, :meth:`Run.to_toml`). :meth:`Run.run` releases
the GIL while it transports the ions and returns a :class:`RunResult` whose
histograms are NumPy arrays and whose summary is a dict. A run with a
``[dynamic]`` section goes through :meth:`Run.run_dynamic`, which returns a
:class:`DynamicRunResult` (fluence-step series and slab composition). An
electron input (an ``[electron]`` table) is read with :class:`ElectronRun`,
whose :meth:`ElectronRun.run` returns an :class:`ElectronRunResult` (escape
spectra and energy-deposition grids as NumPy arrays, summary as a dict).
The dose of a multi-shot exposure from a fitted Gaussian-sum point-spread
function is :func:`exposure_dose` / :func:`exposure_dose_grid`, its large-area
limit :func:`area_dose` and its iso-dose contours :func:`iso_contours`.
A direct Monte Carlo evaluation of the same shot list, compared with the
convolution, is :func:`exposure_direct`.
"""

from ._lindhard import (
    FATE_NAMES,
    Beam,
    DynamicRunResult,
    ElectronRun,
    ElectronRunResult,
    Element,
    Histogram,
    InputError,
    Layer,
    LindhardError,
    Material,
    Physics,
    Run,
    RunError,
    RunResult,
    Tally,
    Target,
    __version__,
    area_dose,
    exposure_dose,
    exposure_dose_grid,
    exposure_direct,
    iso_contours,
    read_damage_export,
)

__all__ = [
    "FATE_NAMES",
    "Beam",
    "DynamicRunResult",
    "ElectronRun",
    "ElectronRunResult",
    "Element",
    "Histogram",
    "InputError",
    "Layer",
    "LindhardError",
    "Material",
    "Physics",
    "Run",
    "RunError",
    "RunResult",
    "Tally",
    "Target",
    "__version__",
    "area_dose",
    "exposure_dose",
    "exposure_direct",
    "exposure_dose_grid",
    "iso_contours",
    "read_damage_export",
]
