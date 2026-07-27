# -*- coding: utf-8 -*-

"""Passive input and context contracts for thermoelastic calibration."""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import numpy as np
from numpy.typing import NDArray
from quantas.core.physics.elasticity import validate_stiffness_matrix
from quantas.models import ResultData
from quantas.models.structures import CrystalStructure, SymmetryMetadata
from .types import FloatArray, QHAThermoelasticPayload, ThermoelasticMethod


@dataclass(slots=True)
class ElasticVolumePoint:
    """One hydrostatically pre-stressed elastic tensor at a sampled volume.

    Parameters
    ----------
    source : str or Path
        Source CRYSTAL output file.
    pressure : float
        Pressure used by CRYSTAL for the hydrostatic pre-stress correction, in
        GPa.
    stress_pressure : float
        Pressure calculated from the final unstrained stress tensor, in GPa.
    volume : float
        Primitive-cell volume in angstrom cubed.
    density : float
        Crystal density in kg m^-3.
    energy : float
        Static DFT energy in hartree.
    stiffness : array_like
        Symmetric ``(6, 6)`` Wallace stiffness matrix in GPa.
    lattice : array_like
        Final primitive direct-lattice vectors, stored by rows in angstrom.
    prestress_applied : bool, optional
        Whether the source explicitly used CRYSTAL's ``PRESSURE`` keyword.
    metadata : dict, optional
        Frame-normalization diagnostics and source provenance.

    Raises
    ------
    ValueError
        If scalar values or arrays are invalid.
    """

    source: str | Path
    pressure: float
    stress_pressure: float
    volume: float
    density: float
    energy: float
    stiffness: FloatArray
    lattice: FloatArray
    prestress_applied: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays and validate one elastic-volume point."""
        self.source = str(self.source)
        self.pressure = float(self.pressure)
        self.stress_pressure = float(self.stress_pressure)
        self.volume = float(self.volume)
        self.density = float(self.density)
        self.energy = float(self.energy)
        self.stiffness = validate_stiffness_matrix(self.stiffness, copy=True)
        self.lattice = np.asarray(self.lattice, dtype=np.float64)
        if self.lattice.shape != (3, 3):
            raise ValueError("lattice must have shape (3, 3)")
        if np.any(~np.isfinite(self.lattice)):
            raise ValueError("lattice must contain finite values")
        if not np.isfinite(self.pressure):
            raise ValueError("pressure must be finite")
        if not np.isfinite(self.stress_pressure):
            raise ValueError("stress_pressure must be finite")
        if not np.isfinite(self.volume) or self.volume <= 0.0:
            raise ValueError("volume must be finite and positive")
        if not np.isfinite(self.density) or self.density <= 0.0:
            raise ValueError("density must be finite and positive")
        if not np.isfinite(self.energy):
            raise ValueError("energy must be finite")
        lattice_volume = abs(float(np.linalg.det(self.lattice)))
        if not np.isclose(lattice_volume, self.volume, rtol=2.0e-6, atol=1.0e-6):
            raise ValueError(
                "lattice determinant and reported primitive-cell volume disagree"
            )
        self.metadata = dict(self.metadata)

    def as_dict(self) -> dict[str, Any]:
        """Return a recursively serializable point mapping."""
        return {
            "source": str(self.source),
            "pressure": float(self.pressure),
            "stress_pressure": float(self.stress_pressure),
            "volume": float(self.volume),
            "density": float(self.density),
            "energy": float(self.energy),
            "prestress_applied": bool(self.prestress_applied),
            "lattice": self.lattice.copy(),
            "stiffness": self.stiffness.copy(),
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class ElasticVolumeSeries:
    """Volume-dependent second-order elastic data in a common CRYSTAL frame.

    Parameters
    ----------
    points : tuple of ElasticVolumePoint
        Elastic data sorted by increasing primitive-cell volume.
    reference_structure : CrystalStructure
        Compact primitive structure used to define species, atom ordering, and
        the common CRYSTAL Cartesian frame.
    symmetry : SymmetryMetadata
        Crystallographic symmetry of the reference structure.
    elastic_symmetry : str
        Elastic crystal-system pattern detected from the stiffness matrices.
    reference_index : int
        Index of the reference point in the sorted series.
    orientation : str, optional
        Description of the tensor and lattice frame.
    metadata : dict, optional
        Additional provenance and validation metadata.

    Raises
    ------
    ValueError
        If the series is empty, unsorted, duplicated, or inconsistent.
    """

    points: tuple[ElasticVolumePoint, ...]
    reference_structure: CrystalStructure
    symmetry: SymmetryMetadata
    elastic_symmetry: str
    reference_index: int
    orientation: str = "crystal"
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate series ordering and reference consistency."""
        self.points = tuple(self.points)
        if not self.points:
            raise ValueError("elastic volume series cannot be empty")
        if self.reference_index < 0 or self.reference_index >= len(self.points):
            raise ValueError("reference_index is outside the elastic series")
        volumes = self.volumes
        if np.any(np.diff(volumes) <= 0.0):
            raise ValueError("elastic points must have unique increasing volumes")
        reference = self.points[self.reference_index]
        if not np.allclose(
            reference.lattice,
            self.reference_structure.lattice,
            rtol=0.0,
            atol=2.0e-7,
        ):
            raise ValueError(
                "reference_structure lattice does not match the reference point"
            )
        if not all(point.prestress_applied for point in self.points):
            raise ValueError(
                "all elastic points must include CRYSTAL hydrostatic pre-stress terms"
            )

    @property
    def npoints(self) -> int:
        """Return the number of sampled elastic volumes."""
        return len(self.points)

    @property
    def volumes(self) -> FloatArray:
        """Return sampled primitive-cell volumes in angstrom cubed."""
        return np.asarray([point.volume for point in self.points], dtype=np.float64)

    @property
    def pressures(self) -> FloatArray:
        """Return CRYSTAL elastic pre-stress pressures in GPa."""
        return np.asarray([point.pressure for point in self.points], dtype=np.float64)

    @property
    def stress_pressures(self) -> FloatArray:
        """Return final stress-tensor pressures in GPa."""
        return np.asarray(
            [point.stress_pressure for point in self.points],
            dtype=np.float64,
        )

    @property
    def densities(self) -> FloatArray:
        """Return sampled densities in kg m^-3."""
        return np.asarray([point.density for point in self.points], dtype=np.float64)

    @property
    def energies(self) -> FloatArray:
        """Return sampled static DFT energies in hartree."""
        return np.asarray([point.energy for point in self.points], dtype=np.float64)

    @property
    def stiffness(self) -> FloatArray:
        """Return sampled stiffness matrices with shape ``(npoints, 6, 6)``."""
        return np.asarray([point.stiffness for point in self.points], dtype=np.float64)

    @property
    def lattices(self) -> FloatArray:
        """Return sampled primitive lattice matrices."""
        return np.asarray([point.lattice for point in self.points], dtype=np.float64)

    @property
    def volume_bounds(self) -> tuple[float, float]:
        """Return minimum and maximum sampled elastic volumes."""
        return float(self.volumes[0]), float(self.volumes[-1])


@dataclass(slots=True)
class ThermoelasticInput:
    """Input contract for quasi-static thermoelastic calculations.

    Parameters
    ----------
    jobname : str
        Name or description of the workflow.
    elastic_series : ElasticVolumeSeries
        Volume-dependent CRYSTAL elastic data.
    method : {"quasistatic"}, optional
        Thermoelastic approximation represented by the input.
    source : str, Path, or None, optional
        YAML source path when the object was read from disk.
    metadata : dict, optional
        Schema and provenance metadata.
    """

    jobname: str
    elastic_series: ElasticVolumeSeries
    method: ThermoelasticMethod = "quasistatic"
    source: str | Path | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ThermoelasticContext:
    """Validated pairing of elastic-volume and QHA thermo-structural data.

    Parameters
    ----------
    input_data : ThermoelasticInput
        Parsed elastic-volume input.
    qha_result_data : ResultData
        Complete QHA result envelope, read from HDF5 or calculated from YAML.
    qha : object implementing the QHA thermoelastic protocol
        QHA payload exposing the thermo-structural fields required by the
        quasi-static coupling.
    extrapolation_mask : ndarray
        Boolean ``(nT, nP)`` mask marking QHA equilibrium volumes outside the
        elastic-volume interval.
    missing_qha_fields : tuple of str, optional
        QHA properties absent from the supplied record.  Equilibrium volume is
        mandatory; other fields are reported for staged implementation.
    metadata : dict, optional
        Context provenance and volume-coverage diagnostics.
    """

    input_data: ThermoelasticInput
    qha_result_data: ResultData
    qha: QHAThermoelasticPayload
    extrapolation_mask: NDArray[np.bool_]
    missing_qha_fields: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def has_complete_quasistatic_inputs(self) -> bool:
        """Return whether QHA data needed by the cold QSA fit are present.

        Returns
        -------
        bool
            ``True`` when the sampled static energy-volume path and the
            pressure-temperature equilibrium-volume map are available.
        """
        required = {
            "temperature",
            "pressure",
            "volume",
            "static_energy",
            "equilibrium_volume",
        }
        return not required.intersection(self.missing_qha_fields)

    @property
    def has_complete_adiabatic_inputs(self) -> bool:
        """Return whether anisotropic isothermal-to-adiabatic data exist.

        Returns
        -------
        bool
            ``True`` when temperature, volume, isochoric heat capacity, and
            the Cartesian thermal-expansion tensor are available.
        """
        required = {
            "isochoric_heat_capacity",
            "thermal_expansion_tensor",
            "equilibrium_volume",
            "temperature",
            "pressure",
        }
        return not required.intersection(self.missing_qha_fields)


__all__ = [
    "ElasticVolumePoint",
    "ElasticVolumeSeries",
    "ThermoelasticContext",
    "ThermoelasticInput",
]
