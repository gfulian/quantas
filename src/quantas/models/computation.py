# -*- coding: utf-8 -*-

"""Passive computational records shared by electronic-structure interfaces."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np
from numpy.typing import NDArray

from quantas.models.structures import CrystalStructure


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


class EnergyKind(str, Enum):
    """Semantic kind of an energy reported by an external code.

    The enumeration deliberately remains backend-neutral. Interface parsers
    should retain code-specific source information in :class:`EnergyRecord`
    metadata instead of introducing backend names here.
    """

    DFT = "dft"
    TOTAL = "total"
    REFERENCE = "reference"


@dataclass(slots=True)
class EnergyRecord:
    """One energy value together with its physical and provenance semantics.

    Parameters
    ----------
    value : float
        Energy value in the unit identified by ``unit``.
    unit : str
        Explicit unit label, for example ``"Ha"`` or ``"eV"``.
    kind : EnergyKind
        Semantic meaning of the energy.
    metadata : dict, optional
        Parser-specific provenance such as source marker or line index.
    """

    value: float
    unit: str
    kind: EnergyKind
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize the scalar value and validate the unit label."""
        self.value = float(self.value)
        self.unit = str(self.unit).strip()
        if not self.unit:
            raise ValueError("energy unit cannot be empty")


class RunTerminationStatus(str, Enum):
    """Termination state of one external-program output stream."""

    NORMAL = "normal"
    INCOMPLETE = "incomplete"
    UNKNOWN = "unknown"


@dataclass(slots=True)
class RunTermination:
    """Observed termination information for one external-program output.

    Parameters
    ----------
    status : RunTerminationStatus
        Parsed global termination state.
    line_index : int or None, optional
        Zero-based source line that established the state when available.
    metadata : dict, optional
        Interface-specific provenance associated with the marker.
    """

    status: RunTerminationStatus
    line_index: int | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


class SCFStatus(str, Enum):
    """Observed convergence state of one self-consistent-field block."""

    CONVERGED = "converged"
    TOO_MANY_CYCLES = "too_many_cycles"
    INCOMPLETE = "incomplete"
    UNKNOWN = "unknown"


@dataclass(slots=True)
class SCFResult:
    """Convergence history of one self-consistent-field calculation.

    Parameters
    ----------
    status : SCFStatus
        Observed termination state of the SCF block.
    cycle_numbers : array_like
        Integer SCF cycle labels in source order.
    energies : array_like
        Electronic energies associated with ``cycle_numbers``.
    delta_energies : array_like
        Energy differences printed for the corresponding cycles.
    energy_unit : str
        Unit shared by ``energies`` and ``delta_energies``.
    final_energy : EnergyRecord or None, optional
        Final energy reported for the block when independently identifiable.
    metadata : dict, optional
        Parser-specific provenance.

    Raises
    ------
    ValueError
        If convergence arrays have inconsistent shapes or the unit is empty.
    """

    status: SCFStatus
    cycle_numbers: IntArray
    energies: FloatArray
    delta_energies: FloatArray
    energy_unit: str
    final_energy: EnergyRecord | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays to Quantas dtypes and validate their shapes."""
        self.cycle_numbers = np.asarray(self.cycle_numbers, dtype=np.int64)
        self.energies = np.asarray(self.energies, dtype=np.float64)
        self.delta_energies = np.asarray(self.delta_energies, dtype=np.float64)
        self.energy_unit = str(self.energy_unit).strip()
        if not self.energy_unit:
            raise ValueError("SCF energy unit cannot be empty")
        if self.cycle_numbers.ndim != 1:
            raise ValueError("cycle_numbers must be one-dimensional")
        if self.energies.ndim != 1 or self.delta_energies.ndim != 1:
            raise ValueError("SCF energy arrays must be one-dimensional")
        if not (
            self.cycle_numbers.shape == self.energies.shape == self.delta_energies.shape
        ):
            raise ValueError("SCF convergence arrays must have identical shapes")

    @property
    def cycles(self) -> int:
        """Return the number of parsed SCF iterations.

        Returns
        -------
        int
            Number of entries in the convergence history.
        """
        return int(self.cycle_numbers.size)


class OptimizationStatus(str, Enum):
    """Observed convergence state of one geometry-optimization run."""

    CONVERGED = "converged"
    FAILED = "failed"
    INCOMPLETE = "incomplete"
    UNKNOWN = "unknown"


@dataclass(slots=True)
class OptimizationStep:
    """One geometry-optimization step reported by an external code.

    Parameters
    ----------
    index : int
        One-based optimization-point index reported by the source code.
    energy : EnergyRecord or None, optional
        Energy associated with this optimization point.
    delta_energy : float or None, optional
        Energy difference from the previous optimization point, in the same
        unit as ``energy`` when available.
    structure : CrystalStructure or None, optional
        Structure associated with the step when parsed.
    max_gradient : float or None, optional
        Maximum gradient convergence metric.
    rms_gradient : float or None, optional
        Root-mean-square gradient convergence metric.
    max_displacement : float or None, optional
        Maximum displacement convergence metric.
    rms_displacement : float or None, optional
        Root-mean-square displacement convergence metric.
    metadata : dict, optional
        Parser-specific provenance.
    """

    index: int
    energy: EnergyRecord | None = None
    delta_energy: float | None = None
    structure: CrystalStructure | None = None
    max_gradient: float | None = None
    rms_gradient: float | None = None
    max_displacement: float | None = None
    rms_displacement: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize scalar values and validate the step index."""
        self.index = int(self.index)
        if self.index < 1:
            raise ValueError("optimization step index must be positive")
        for name in (
            "delta_energy",
            "max_gradient",
            "rms_gradient",
            "max_displacement",
            "rms_displacement",
        ):
            value = getattr(self, name)
            if value is not None:
                setattr(self, name, float(value))


@dataclass(slots=True)
class OptimizationResult:
    """Convergence record for one geometry-optimization run.

    Parameters
    ----------
    status : OptimizationStatus
        Observed optimization termination state.
    steps : tuple of OptimizationStep
        Parsed optimization points in source order.
    final_energy : EnergyRecord or None, optional
        Final total energy reported by the optimization end marker.
    final_structure : CrystalStructure or None, optional
        Final optimized structure when independently parsed.
    metadata : dict, optional
        Parser-specific provenance.
    """

    status: OptimizationStatus
    steps: tuple[OptimizationStep, ...]
    final_energy: EnergyRecord | None = None
    final_structure: CrystalStructure | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def cycles(self) -> int:
        """Return the number of parsed optimization points.

        Returns
        -------
        int
            Number of optimization steps.
        """
        return len(self.steps)


@dataclass(slots=True)
class SourceProvenance:
    """Origin of one parsed observation supplied to Quantas.

    Parameters
    ----------
    interface : str
        Quantas interface that produced the observation, for example
        ``"crystal"`` or ``"vasp"``.
    source : str or None, optional
        Human-readable source path or identifier.
    record_index : int or None, optional
        Zero-based index of the record within a multi-record source.
    metadata : dict, optional
        Additional interface-specific provenance.

    Raises
    ------
    ValueError
        If the interface label is empty or ``record_index`` is negative.
    """

    interface: str
    source: str | None = None
    record_index: int | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize provenance labels and validate the optional record index."""
        self.interface = str(self.interface).strip()
        if not self.interface:
            raise ValueError("provenance interface cannot be empty")
        if self.source is not None:
            self.source = str(self.source)
        if self.record_index is not None:
            self.record_index = int(self.record_index)
            if self.record_index < 0:
                raise ValueError("record_index cannot be negative")

    def as_dict(self) -> dict[str, Any]:
        """Return a recursively serializable provenance mapping.

        Returns
        -------
        dict
            Interface, source, record index, and metadata.
        """
        return {
            "interface": self.interface,
            "source": self.source,
            "record_index": self.record_index,
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class StructureEnergyPoint:
    """One structure-energy observation with optional source provenance.

    Parameters
    ----------
    structure : CrystalStructure
        Periodic structure expressed in the canonical Quantas structural
        convention. Lattice vectors are therefore in angstrom.
    energy : EnergyRecord
        Static energy associated with ``structure``.
    provenance : SourceProvenance or None, optional
        Origin of the observation.
    metadata : dict, optional
        Additional backend-neutral information associated with the point.
    """

    structure: CrystalStructure
    energy: EnergyRecord
    provenance: SourceProvenance | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def volume(self) -> float:
        """Return the structure volume in cubic angstrom.

        Returns
        -------
        float
            Positive cell volume in cubic angstrom.
        """
        return self.structure.volume

    def as_dict(self) -> dict[str, Any]:
        """Return a recursively serializable observation mapping.

        Returns
        -------
        dict
            Structure, energy, provenance, and metadata.
        """
        data: dict[str, Any] = {
            "structure": self.structure.as_dict(),
            "energy": {
                "value": float(self.energy.value),
                "unit": self.energy.unit,
                "kind": self.energy.kind.value,
                "metadata": dict(self.energy.metadata),
            },
            "metadata": dict(self.metadata),
        }
        if self.provenance is not None:
            data["provenance"] = self.provenance.as_dict()
        return data


@dataclass(slots=True)
class StructureEnergySeries:
    """Analysis-ready sequence of compatible structure-energy observations.

    The series is intentionally independent of how its points were produced.
    A monolithic equation-of-state output and a set of independent
    fixed-volume optimizations can therefore be normalized to the same
    contract before entering an EOS, QHA, or workflow layer.

    Parameters
    ----------
    points : tuple of StructureEnergyPoint
        Ordered structure-energy observations.
    reference_index : int, optional
        Index of the reference observation within ``points``.
    metadata : dict, optional
        Dataset-level metadata and provenance.

    Raises
    ------
    ValueError
        If the series is empty, the reference index is invalid, or energies
        use inconsistent units or semantic kinds.
    """

    points: tuple[StructureEnergyPoint, ...]
    reference_index: int = 0
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize point storage and validate series-level compatibility."""
        self.points = tuple(self.points)
        if not self.points:
            raise ValueError("structure-energy series cannot be empty")
        self.reference_index = int(self.reference_index)
        if self.reference_index < 0 or self.reference_index >= len(self.points):
            raise ValueError("reference_index is outside the structure-energy series")

        unit = self.points[0].energy.unit.casefold()
        kind = self.points[0].energy.kind
        for point in self.points[1:]:
            if point.energy.unit.casefold() != unit:
                raise ValueError("structure-energy series uses inconsistent energy units")
            if point.energy.kind is not kind:
                raise ValueError("structure-energy series uses inconsistent energy kinds")

    @property
    def npoints(self) -> int:
        """Return the number of observations in the series.

        Returns
        -------
        int
            Number of stored points.
        """
        return len(self.points)

    @property
    def volumes(self) -> FloatArray:
        """Return point volumes in cubic angstrom.

        Returns
        -------
        ndarray
            One-dimensional float64 array of structure volumes.
        """
        return np.asarray([point.volume for point in self.points], dtype=np.float64)

    @property
    def energies(self) -> FloatArray:
        """Return point energies in their shared source unit.

        Returns
        -------
        ndarray
            One-dimensional float64 array of static energies.
        """
        return np.asarray(
            [point.energy.value for point in self.points],
            dtype=np.float64,
        )

    @property
    def energy_unit(self) -> str:
        """Return the shared energy unit.

        Returns
        -------
        str
            Unit label stored by the point energy records.
        """
        return self.points[0].energy.unit

    @property
    def energy_kind(self) -> EnergyKind:
        """Return the shared semantic energy kind.

        Returns
        -------
        EnergyKind
            Energy kind stored by every observation.
        """
        return self.points[0].energy.kind

    @property
    def volume_unit(self) -> str:
        """Return the canonical Quantas volume unit.

        Returns
        -------
        str
            ``"angstrom^3"`` because :class:`CrystalStructure` stores lattice
            vectors in angstrom.
        """
        return "angstrom^3"
