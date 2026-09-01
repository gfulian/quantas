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
