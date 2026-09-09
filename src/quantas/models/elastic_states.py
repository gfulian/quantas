# -*- coding: utf-8 -*-

"""Backend-neutral elastic states sampled along a volume path."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray


class ElasticTensorKind(str, Enum):
    """Identify the physical convention of an elastic stiffness tensor."""

    RAW_ENERGY_STRAIN = "raw_energy_strain"
    WALLACE_HYDROSTATIC = "wallace_hydrostatic"
    INCREMENTAL_FULL_STRESS = "incremental_full_stress"
    UNKNOWN = "unknown"

    @property
    def is_incremental(self) -> bool:
        """Return whether the tensor is suitable for Christoffel acoustics."""
        return self in {
            ElasticTensorKind.WALLACE_HYDROSTATIC,
            ElasticTensorKind.INCREMENTAL_FULL_STRESS,
        }


class PressureSource(str, Enum):
    """Identify how the pressure associated with an elastic state was found."""

    OUTPUT_STRESS = "output_stress"
    APPLIED_PRESTRESS = "applied_prestress"
    MANUAL = "manual"
    ENERGY_EOS = "energy_eos"
    ENERGY_POLYNOMIAL = "energy_polynomial"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True, slots=True)
class PrestressProvenance:
    """Describe pressure and finite-stress treatment of one elastic tensor.

    Parameters
    ----------
    tensor_kind : ElasticTensorKind or str
        Physical convention of the stored stiffness matrix.
    pressure_gpa : float or None, optional
        Hydrostatic pressure in GPa. Compression is positive.
    pressure_source : PressureSource or str, optional
        Origin of ``pressure_gpa``.
    correction_method : str or None, optional
        Named finite-stress correction used to create an incremental tensor.
    correction_applied_by : str or None, optional
        Program or workflow that applied the correction.
    source_tensor_kind : ElasticTensorKind, str, or None, optional
        Convention of the tensor before the recorded correction.

    Raises
    ------
    ValueError
        If pressure or correction provenance is inconsistent.
    """

    tensor_kind: ElasticTensorKind | str
    pressure_gpa: float | None = None
    pressure_source: PressureSource | str = PressureSource.UNAVAILABLE
    correction_method: str | None = None
    correction_applied_by: str | None = None
    source_tensor_kind: ElasticTensorKind | str | None = None

    def __post_init__(self) -> None:
        """Normalize enums and reject ambiguous correction histories."""
        kind = ElasticTensorKind(self.tensor_kind)
        source = PressureSource(self.pressure_source)
        source_kind = (
            None
            if self.source_tensor_kind is None
            else ElasticTensorKind(self.source_tensor_kind)
        )
        pressure = self.pressure_gpa
        if pressure is not None:
            pressure = float(pressure)
            if not np.isfinite(pressure):
                raise ValueError("pressure_gpa must be finite when provided")
        if (pressure is None) != (source is PressureSource.UNAVAILABLE):
            raise ValueError(
                "pressure value and pressure source must be provided together"
            )
        method = _optional_text(self.correction_method)
        applied_by = _optional_text(self.correction_applied_by)
        if (method is None) != (applied_by is None):
            raise ValueError("correction method and applied-by provenance must coexist")
        if method is None and source_kind is not None:
            raise ValueError("source_tensor_kind requires a recorded correction")
        if method is not None:
            if not kind.is_incremental:
                raise ValueError(
                    "a corrected tensor must use an incremental tensor kind"
                )
            if source_kind is None:
                raise ValueError(
                    "a corrected tensor must record its source tensor kind"
                )
            if source_kind.is_incremental:
                raise ValueError(
                    "an incremental tensor cannot be corrected a second time"
                )
        object.__setattr__(self, "tensor_kind", kind)
        object.__setattr__(self, "pressure_gpa", pressure)
        object.__setattr__(self, "pressure_source", source)
        object.__setattr__(self, "correction_method", method)
        object.__setattr__(self, "correction_applied_by", applied_by)
        object.__setattr__(self, "source_tensor_kind", source_kind)

    def require_incremental(self) -> None:
        """Require a tensor convention suitable for Christoffel acoustics.

        Raises
        ------
        ValueError
            If the stored tensor is raw or its convention is unknown.
        """
        if not ElasticTensorKind(self.tensor_kind).is_incremental:
            raise ValueError(
                "Christoffel acoustics requires an explicitly incremental stiffness tensor"
            )


@dataclass(slots=True)
class ElasticState:
    """One backend-neutral elastic state at a primitive-cell volume.

    Parameters
    ----------
    volume : float
        Primitive-cell volume in angstrom cubed.
    density : float
        Density in kg m^-3.
    stiffness : array_like
        Symmetric stiffness matrix in GPa with shape ``(6, 6)``.
    prestress : PrestressProvenance
        Tensor convention and pressure provenance.
    energy : float or None, optional
        Static energy in the unit named by ``energy_unit``.
    energy_unit : str or None, optional
        Explicit unit for ``energy``.
    lattice : array_like or None, optional
        Primitive direct-lattice vectors by rows in angstrom.
    source : str, Path, or None, optional
        Native output or other data source.
    metadata : dict, optional
        Additional backend and structural provenance.

    Raises
    ------
    TypeError
        If ``prestress`` has an unsupported type.
    ValueError
        If a scalar, matrix, unit, or lattice is invalid.
    """

    volume: float
    density: float
    stiffness: ArrayLike
    prestress: PrestressProvenance
    energy: float | None = None
    energy_unit: str | None = None
    lattice: ArrayLike | None = None
    source: str | Path | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize values and validate one state."""
        self.volume = float(self.volume)
        self.density = float(self.density)
        if not np.isfinite(self.volume) or self.volume <= 0.0:
            raise ValueError("volume must be finite and positive")
        if not np.isfinite(self.density) or self.density <= 0.0:
            raise ValueError("density must be finite and positive")
        self.stiffness = _stiffness_matrix(self.stiffness)
        if not isinstance(self.prestress, PrestressProvenance):
            raise TypeError("prestress must be a PrestressProvenance instance")
        if self.energy is None:
            if self.energy_unit is not None:
                raise ValueError("energy_unit cannot be supplied without energy")
        else:
            self.energy = float(self.energy)
            if not np.isfinite(self.energy):
                raise ValueError("energy must be finite")
            self.energy_unit = _required_text(self.energy_unit, "energy_unit")
        if self.lattice is not None:
            lattice = np.asarray(self.lattice, dtype=np.float64)
            if lattice.shape != (3, 3) or not np.all(np.isfinite(lattice)):
                raise ValueError("lattice must be finite with shape (3, 3)")
            lattice_volume = abs(float(np.linalg.det(lattice)))
            if not np.isclose(lattice_volume, self.volume, rtol=2.0e-6, atol=1.0e-6):
                raise ValueError("lattice determinant and volume disagree")
            self.lattice = lattice.copy()
        self.source = None if self.source is None else str(self.source)
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class ElasticStateSeries:
    """Increasing volume series of compatible elastic states.

    Parameters
    ----------
    states : tuple of ElasticState
        States in strictly increasing primitive-cell volume order.
    reference_index : int
        Index of the reference state.
    orientation : str, optional
        Common Cartesian tensor-frame description.
    metadata : dict, optional
        Series-level provenance.

    Raises
    ------
    ValueError
        If the series is empty, unordered, duplicated, or inconsistent.
    """

    states: tuple[ElasticState, ...]
    reference_index: int
    orientation: str = "crystal"
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate ordering and cross-state conventions."""
        self.states = tuple(self.states)
        if not self.states:
            raise ValueError("elastic state series cannot be empty")
        if not 0 <= self.reference_index < len(self.states):
            raise ValueError("reference_index is outside the elastic series")
        if np.any(np.diff(self.volumes) <= 0.0):
            raise ValueError("elastic states must have unique increasing volumes")
        energy_units = {
            state.energy_unit for state in self.states if state.energy is not None
        }
        if len(energy_units) > 1:
            raise ValueError("elastic states use inconsistent energy units")
        self.orientation = _required_text(self.orientation, "orientation")
        self.metadata = dict(self.metadata)

    @property
    def nstates(self) -> int:
        """Return the number of volume states."""
        return len(self.states)

    @property
    def volumes(self) -> NDArray[np.float64]:
        """Return primitive-cell volumes in angstrom cubed."""
        return np.asarray([state.volume for state in self.states], dtype=np.float64)

    @property
    def densities(self) -> NDArray[np.float64]:
        """Return densities in kg m^-3."""
        return np.asarray([state.density for state in self.states], dtype=np.float64)

    @property
    def stiffness(self) -> NDArray[np.float64]:
        """Return stiffness matrices in GPa with shape ``(nstates, 6, 6)``."""
        return np.asarray([state.stiffness for state in self.states], dtype=np.float64)

    def require_incremental(self) -> None:
        """Require every state to be usable by Christoffel acoustics.

        Raises
        ------
        ValueError
            If any tensor convention is raw or unknown.
        """
        for index, state in enumerate(self.states):
            try:
                state.prestress.require_incremental()
            except ValueError as exc:
                raise ValueError(f"elastic state {index}: {exc}") from exc


def _optional_text(value: str | None) -> str | None:
    """Normalize optional non-empty text."""
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _stiffness_matrix(values: ArrayLike) -> NDArray[np.float64]:
    """Return a copied finite symmetric ``6 x 6`` stiffness matrix."""
    matrix = np.asarray(values, dtype=np.float64)
    if matrix.shape != (6, 6) or not np.all(np.isfinite(matrix)):
        raise ValueError("stiffness must be finite with shape (6, 6)")
    if not np.allclose(matrix, matrix.T, rtol=0.0, atol=1.0e-10):
        raise ValueError("stiffness matrix must be symmetric")
    return matrix.copy()


def _required_text(value: str | None, name: str) -> str:
    """Normalize required non-empty text."""
    text = _optional_text(value)
    if text is None:
        raise ValueError(f"{name} must be a non-empty string")
    return text


__all__ = [
    "ElasticState",
    "ElasticStateSeries",
    "ElasticTensorKind",
    "PressureSource",
    "PrestressProvenance",
]
