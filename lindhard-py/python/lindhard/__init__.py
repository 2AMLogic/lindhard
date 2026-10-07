"""Python bindings for lindhard.

Configuration classes (:class:`Material`, :class:`Layer`, :class:`Target`,
:class:`Beam`, :class:`Physics`, :class:`Tally`) are assembled into a
:class:`Run`, which reads and writes the same TOML schema as the ``lindhard``
command (:meth:`Run.from_toml`, :meth:`Run.to_toml`). :meth:`Run.run` releases
the GIL while it transports the ions and returns a :class:`RunResult` whose
histograms are NumPy arrays and whose summary is a dict.
"""

from ._lindhard import (
    FATE_NAMES,
    Beam,
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
)

__all__ = [
    "FATE_NAMES",
    "Beam",
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
]
