# -*- coding: utf-8 -*-

"""Resolve hydrostatic elastic pressures from static energy--volume data.

This module combines the backend-neutral ``P(V) = -dE/dV`` services with the
explicit volume-matching and pressure-assignment contracts used by elastic
volume series.  It does not apply any finite-pressure correction to stiffness
tensors; backend-specific transformations remain the responsibility of the
corresponding interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.physics.eos import (
    PressureEstimate,
    pressure_from_energy_eos,
    pressure_from_energy_polynomial,
)
from quantas.models.elastic_states import ElasticStateSeries, PressureSource
from quantas.models.volume_matching import VolumeMatch, match_sampled_volumes

from .prestress import assign_hydrostatic_pressures


_ENERGY_PRESSURE_SOURCES = {
    PressureSource.ENERGY_EOS,
    PressureSource.ENERGY_POLYNOMIAL,
}


@dataclass(slots=True)
class EnergyPressureResolution:
    """Result of assigning energy-derived pressures to an elastic series.

    Parameters
    ----------
    series : ElasticStateSeries
        Raw elastic series with one derived pressure attached to each state.
        Stiffness coefficients are unchanged.
    estimate : PressureEstimate
        Pressure estimate evaluated at the source energy-volume samples.
    pressures_gpa : ndarray
        Pressures assigned to the elastic states, in elastic-series order.
    matches : tuple of VolumeMatch
        Explicit elastic-volume to energy-volume associations.
    provenance : dict
        Backend-neutral fit provenance suitable for embedding in workflow
        metadata.
    """

    series: ElasticStateSeries
    estimate: PressureEstimate
    pressures_gpa: NDArray[np.float64]
    matches: tuple[VolumeMatch, ...]
    provenance: dict[str, Any]

    def __post_init__(self) -> None:
        """Normalize passive arrays and provenance containers."""
        self.pressures_gpa = np.asarray(self.pressures_gpa, dtype=np.float64).copy()
        self.matches = tuple(self.matches)
        self.provenance = dict(self.provenance)
        if self.pressures_gpa.shape != (self.series.nstates,):
            raise ValueError("pressures_gpa must contain one value per elastic state")
        if not np.all(np.isfinite(self.pressures_gpa)):
            raise ValueError("pressures_gpa must contain only finite values")


def resolve_energy_derived_pressures(
    series: ElasticStateSeries,
    energy_volumes: ArrayLike,
    energies: ArrayLike,
    *,
    pressure_source: PressureSource | str,
    source_dataset: str,
    energy_unit: str,
    volume_length_unit: str,
    volume_unit: str | None = None,
    eos: str = "BM3",
    polynomial_degree: int = 3,
    maxfev: int | None = None,
) -> EnergyPressureResolution:
    """Fit ``E(V)``, match volumes, and attach hydrostatic pressures.

    Parameters
    ----------
    series : ElasticStateSeries
        Increasing raw elastic-state series requiring external pressures.
    energy_volumes, energies : array_like
        Static energy-volume samples. Volumes must use the primitive-cell
        normalization of ``series``.
    pressure_source : PressureSource or str
        Energy-derived pressure method. Only ``energy_eos`` and
        ``energy_polynomial`` are accepted.
    source_dataset : str
        Stable provenance label or path identifying the energy dataset.
    energy_unit : str
        Unit of ``energies``.
    volume_length_unit : str
        Length unit whose cube defines ``energy_volumes``; for example,
        ``"angstrom"`` for values in cubic angstrom.
    volume_unit : str or None, optional
        Human-readable unit of the sampled volumes retained in provenance.
    eos : str, optional
        Integrated energy EOS used by ``energy_eos``.
    polynomial_degree : int, optional
        Polynomial degree used by ``energy_polynomial``.
    maxfev : int or None, optional
        Optional maximum number of EOS fitting function evaluations.

    Returns
    -------
    EnergyPressureResolution
        Assigned raw elastic series, fit diagnostics, explicit volume matches,
        and backend-neutral provenance.

    Raises
    ------
    TypeError
        If ``series`` is not an :class:`ElasticStateSeries`.
    ValueError
        If the pressure source, energy data, volume matching, fit, or pressure
        assignment is invalid.

    Notes
    -----
    This function intentionally stops before any finite-prestress or backend-specific
    finite-pressure correction.  Interfaces must transform the assigned raw
    tensors according to the convention emitted by the external code.
    """
    if not isinstance(series, ElasticStateSeries):
        raise TypeError("series must be an ElasticStateSeries")
    source = PressureSource(pressure_source)
    if source not in _ENERGY_PRESSURE_SOURCES:
        choices = ", ".join(sorted(item.value for item in _ENERGY_PRESSURE_SOURCES))
        raise ValueError(
            f"pressure_source must be energy-derived; choose one of {choices}"
        )
    dataset = str(source_dataset).strip()
    if not dataset:
        raise ValueError("source_dataset must be a non-empty string")

    volumes = np.asarray(energy_volumes, dtype=np.float64)
    energy = np.asarray(energies, dtype=np.float64)
    if (
        volumes.ndim != 1
        or energy.ndim != 1
        or volumes.shape != energy.shape
        or volumes.size < 3
        or not np.all(np.isfinite(volumes))
        or not np.all(np.isfinite(energy))
        or np.any(volumes <= 0.0)
    ):
        raise ValueError(
            "energy-derived pressure requires at least three finite aligned "
            "positive volume-energy points"
        )

    if source is PressureSource.ENERGY_EOS:
        estimate = pressure_from_energy_eos(
            volumes,
            energy,
            eos=eos,
            energy_unit=energy_unit,
            volume_unit=volume_length_unit,
            pressure_unit="GPa",
            maxfev=maxfev,
        )
    else:
        estimate = pressure_from_energy_polynomial(
            volumes,
            energy,
            degree=polynomial_degree,
            energy_unit=energy_unit,
            volume_unit=volume_length_unit,
            pressure_unit="GPa",
        )
    if not estimate.success:
        detail = estimate.fit.message or "fit did not return finite pressures"
        raise ValueError(f"{source.value} pressure fit failed: {detail}")

    matches = match_sampled_volumes(series.volumes, volumes)
    pressures = np.asarray(
        [estimate.pressure[match.source_index] for match in matches],
        dtype=np.float64,
    )
    provenance: dict[str, Any] = {
        "method": source.value,
        "relation": "P(V) = -dE/dV",
        "source_dataset": dataset,
        "energy_unit": str(energy_unit),
        "volume_length_unit": str(volume_length_unit),
        "pressure_unit": estimate.unit,
        "settings": dict(estimate.metadata),
        "fit": estimate.fit.as_dict(),
        "warnings": list(estimate.warnings),
        "volume_matches": [
            {
                "target_index": match.target_index,
                "source_index": match.source_index,
                "target_volume": match.target_volume,
                "source_volume": match.source_volume,
                "absolute_difference": match.absolute_difference,
                "relative_difference": match.relative_difference,
            }
            for match in matches
        ],
    }
    if volume_unit is not None:
        provenance["volume_unit"] = str(volume_unit)

    assigned = assign_hydrostatic_pressures(
        series,
        pressures,
        pressure_source=source,
        assignment_method=source.value,
        metadata=provenance,
    )
    return EnergyPressureResolution(
        series=assigned,
        estimate=estimate,
        pressures_gpa=pressures,
        matches=matches,
        provenance=provenance,
    )


__all__ = ["EnergyPressureResolution", "resolve_energy_derived_pressures"]
