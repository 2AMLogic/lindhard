"""Type stubs for the compiled module ``lindhard._lindhard``."""

import os
from typing import Any, Mapping, Sequence

import numpy as np
import numpy.typing as npt

__version__: str

FATE_NAMES: tuple[str, str, str, str]
"""Names of the codes in ``RunResult.ions["fate"]``: stopped, backscattered, transmitted, lateral."""

class LindhardError(Exception):
    """Base class of every exception raised by lindhard."""

class InputError(LindhardError, ValueError):
    """The configuration is malformed or has an invalid value."""

class RunError(LindhardError, RuntimeError):
    """The configuration was accepted but the run failed."""

class Element:
    """One element of a material (give `symbol` or `z`, and one fraction kind)."""

    symbol: str | None
    z: int | None
    atom_fraction: float | None
    mass_fraction: float | None
    e_d_ev: float | None
    e_b_ev: float | None
    e_s_ev: float | None
    def __init__(
        self,
        symbol: str | None = None,
        *,
        z: int | None = None,
        atom_fraction: float | None = None,
        mass_fraction: float | None = None,
        e_d_ev: float | None = None,
        e_b_ev: float | None = None,
        e_s_ev: float | None = None,
    ) -> None: ...

class Material:
    """Elements with fractions and a mass density (required for compounds)."""

    elements: list[Element]
    density_g_cm3: float | None
    name: str | None
    def __init__(
        self,
        elements: Sequence[Element],
        density_g_cm3: float | None = None,
        name: str | None = None,
    ) -> None: ...

class Layer:
    """One finite layer: a material (name or `Material`) and a thickness in nm."""

    material: str | Material
    thickness_nm: float
    def __init__(self, material: str | Material, thickness_nm: float) -> None: ...

class Target:
    """Finite layers front to back, then an optional semi-infinite substrate."""

    layers: list[Layer]
    substrate: str | Material | None
    def __init__(
        self,
        layers: Sequence[Layer] | None = None,
        substrate: str | Material | None = None,
    ) -> None: ...

class Beam:
    """The incident beam."""

    ion: str
    energy_ev: float
    mass_amu: float | None
    tilt_deg: float
    azimuth_deg: float
    def __init__(
        self,
        ion: str,
        energy_ev: float,
        mass_amu: float | None = None,
        tilt_deg: float = 0.0,
        azimuth_deg: float = 0.0,
    ) -> None: ...

class Physics:
    """Model choices and cutoffs (the `[physics]` table)."""

    primary_cutoff_ev: float
    recoil_cutoff_ev: float
    potential: str
    screening_length: str | None
    stopping: str
    free_path: str
    min_cm_angle_deg: float | None
    weak_collisions: int
    follow_recoils: bool
    primary_surface_binding_ev: float
    energies: dict[str, dict[str, float]]
    def __init__(
        self,
        primary_cutoff_ev: float,
        recoil_cutoff_ev: float,
        *,
        potential: str = "zbl",
        screening_length: str | None = None,
        stopping: str = "lindhard-scharff",
        free_path: str = "constant",
        min_cm_angle_deg: float | None = None,
        weak_collisions: int = 0,
        follow_recoils: bool = True,
        primary_surface_binding_ev: float = 0.0,
        energies: Mapping[str, Mapping[str, float]] | None = None,
    ) -> None: ...

class Tally:
    """What the run records (the `[tally]` table); defaults are the CLI's."""

    depth_bin_nm: float
    depth_bins: int
    per_ion: bool
    lateral_bin_nm: float
    lateral_bins: int
    escape_energy_max_ev: float | None
    escape_energy_bins: int
    escape_polar_bins: int
    dual_pearson: bool
    def __init__(
        self,
        depth_bin_nm: float | None = None,
        depth_bins: int | None = None,
        per_ion: bool | None = None,
        lateral_bin_nm: float | None = None,
        lateral_bins: int | None = None,
        escape_energy_max_ev: float | None = None,
        escape_energy_bins: int | None = None,
        escape_polar_bins: int | None = None,
        dual_pearson: bool | None = None,
    ) -> None: ...

class Histogram:
    """Counts over uniform bins with edges and per-ion densities."""

    @property
    def edges(self) -> npt.NDArray[np.float64]: ...
    @property
    def centers(self) -> npt.NDArray[np.float64]: ...
    @property
    def counts(self) -> npt.NDArray[np.uint64]: ...
    @property
    def density(self) -> npt.NDArray[np.float64]: ...
    @property
    def underflow(self) -> int: ...
    @property
    def overflow(self) -> int: ...
    @property
    def unit(self) -> str: ...
    def __len__(self) -> int: ...

class RunResult:
    """The tallies of a finished run."""

    @property
    def histories(self) -> int: ...
    @property
    def depth(self) -> Histogram: ...
    @property
    def lateral_y(self) -> Histogram: ...
    @property
    def lateral_z(self) -> Histogram: ...
    @property
    def radial(self) -> Histogram: ...
    @property
    def vacancies(self) -> Histogram: ...
    @property
    def interstitials(self) -> Histogram: ...
    @property
    def replacements(self) -> Histogram: ...
    @property
    def escapes(self) -> list[dict[str, Any]]: ...
    @property
    def ions(self) -> dict[str, npt.NDArray[Any]] | None: ...
    def summary(self) -> dict[str, Any]: ...
    def summary_json(self) -> str: ...
    def write(self, out_dir: str | os.PathLike[str]) -> None: ...

class Run:
    """A complete run description: the whole input document of the CLI."""

    base_dir: str | os.PathLike[str] | None
    beam: Beam
    target: Target
    physics: Physics
    tally: Tally
    materials: dict[str, Material]
    stopping_tables: list[str]
    ions: int
    seed: int
    threads: int | None
    def __init__(
        self,
        beam: Beam,
        target: Target,
        physics: Physics,
        ions: int,
        seed: int,
        *,
        threads: int | None = None,
        tally: Tally | None = None,
        materials: Mapping[str, Material] | None = None,
        stopping_tables: Sequence[str] | None = None,
        base_dir: str | os.PathLike[str] | None = None,
    ) -> None: ...
    @staticmethod
    def from_toml(text: str, base_dir: str | os.PathLike[str] | None = None) -> Run: ...
    @staticmethod
    def from_toml_file(path: str | os.PathLike[str]) -> Run: ...
    def to_toml(self) -> str: ...
    def validate(self) -> list[str]: ...
    def run(
        self,
        *,
        ions: int | None = None,
        seed: int | None = None,
        threads: int | None = None,
    ) -> RunResult: ...
