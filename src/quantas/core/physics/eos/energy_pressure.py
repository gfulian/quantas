# -*- coding: utf-8 -*-

"""Pressure estimates obtained from static energy--volume data.

The functions in this module implement the thermodynamic relation
``P(V) = -dE/dV`` using either a polynomial representation of ``E(V)`` or an
integrated energy equation of state.  They are intentionally independent of
HA/QHA frontends so the same fitted pressure provenance can be reused by
elastic and equation-of-state workflows.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.math.fitting import FitQuality, FitResult, FitStatus, validate_xy
from quantas.core.math.polynomials import fit_polynomial
from quantas.core.physics.units import energy_to_pressure

from .energy import EnergyEOS


@dataclass(slots=True)
class PressureEstimate:
    """Pressure values and diagnostics for one energy--volume fit.

    Parameters
    ----------
    method : str
        Stable method identifier, currently ``"polynomial"`` or ``"eos"``.
    pressure : ndarray
        Pressures evaluated at the supplied volumes.
    fit : FitResult
        Complete diagnostics for the fitted energy representation.
    unit : str
        Unit of ``pressure``.
    warnings : list of str, optional
        Non-fatal fit or pressure-evaluation diagnostics.
    metadata : dict, optional
        Method-specific settings such as polynomial degree or EOS tag.
    """

    method: str
    pressure: NDArray[np.float64]
    fit: FitResult
    unit: str
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays and passive provenance containers."""
        self.method = str(self.method).strip()
        self.pressure = np.asarray(self.pressure, dtype=np.float64).copy()
        self.unit = str(self.unit).strip()
        self.warnings = [str(value) for value in self.warnings]
        self.metadata = dict(self.metadata)
        if not self.method:
            raise ValueError("pressure-estimate method must be non-empty")
        if self.pressure.ndim != 1:
            raise ValueError("pressure estimate must be one-dimensional")
        if not self.unit:
            raise ValueError("pressure-estimate unit must be non-empty")

    @property
    def success(self) -> bool:
        """Return whether the fit and all evaluated pressures are usable."""
        return bool(
            self.fit.success
            and self.pressure.size > 0
            and np.all(np.isfinite(self.pressure))
        )

    @property
    def pressure_min(self) -> float | None:
        """Return the minimum evaluated pressure, when available."""
        if self.pressure.size == 0:
            return None
        return float(np.nanmin(self.pressure))

    @property
    def pressure_max(self) -> float | None:
        """Return the maximum evaluated pressure, when available."""
        if self.pressure.size == 0:
            return None
        return float(np.nanmax(self.pressure))

    def as_dict(self) -> dict[str, Any]:
        """Return a YAML- and JSON-serializable representation."""
        return {
            "method": self.method,
            "success": self.success,
            "pressure": self.pressure.tolist(),
            "pressure_min": self.pressure_min,
            "pressure_max": self.pressure_max,
            "unit": self.unit,
            "fit": self.fit.as_dict(),
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


def pressure_from_energy_polynomial(
    volume: ArrayLike,
    energy: ArrayLike,
    *,
    degree: int,
    energy_unit: str,
    volume_unit: str,
    pressure_unit: str = "GPa",
) -> PressureEstimate:
    """Estimate ``-dE/dV`` from a polynomial energy--volume fit.

    The polynomial coordinate is centred and scaled over the sampled interval
    before fitting; its transform is retained in the fit metadata and the
    derivative is converted back to physical-volume coordinates.
    ``volume_unit`` is the length unit whose cube defines the input volumes,
    for example ``"angstrom"`` for values expressed in angstrom cubed.
    """
    volume_array, energy_array = _energy_volume_arrays(volume, energy)
    fit, fitted_polynomial = fit_polynomial(
        volume_array,
        energy_array,
        int(degree),
        scale_coordinate=True,
    )
    if not fit.success or fitted_polynomial is None:
        return PressureEstimate(
            "polynomial",
            np.asarray([], dtype=np.float64),
            fit,
            pressure_unit,
            [fit.message],
            {"degree": int(degree)},
        )

    pressure_energy_density = -fitted_polynomial.derivative(volume_array)
    pressure = np.asarray(
        energy_to_pressure(
            pressure_energy_density,
            energy_unit,
            volume_unit,
            pressure_unit,
        ),
        dtype=np.float64,
    )
    warnings: list[str] = []
    if fit.quality is FitQuality.POOR:
        warnings.append("the polynomial pressure estimate is based on a poor fit")
    warnings.extend(fit.warnings)
    return PressureEstimate(
        method="polynomial",
        pressure=pressure,
        fit=fit,
        unit=pressure_unit,
        warnings=warnings,
        metadata={"degree": int(degree)},
    )


def pressure_from_energy_eos(
    volume: ArrayLike,
    energy: ArrayLike,
    *,
    eos: str,
    energy_unit: str,
    volume_unit: str,
    pressure_unit: str = "GPa",
    maxfev: int | None = None,
) -> PressureEstimate:
    """Estimate ``-dE/dV`` from an integrated energy EOS fit.

    ``volume_unit`` is the length unit whose cube defines the input volumes.
    The returned metadata records the canonical EOS tag, family, and order.
    """
    volume_array, energy_array = _energy_volume_arrays(volume, energy)
    model = EnergyEOS()
    try:
        model_spec = model.model(eos)
    except ValueError as exc:
        return _failed_estimate(
            "eos",
            str(exc),
            unit=pressure_unit,
            metadata={"eos": eos},
        )

    fit = model.fit(model_spec, volume_array, energy_array, maxfev=maxfev)
    metadata = {
        "eos": model_spec.tag,
        "eos_family": model_spec.family.value,
        "eos_order": model_spec.order,
    }
    if not fit.success or fit.parameters is None:
        return PressureEstimate(
            "eos",
            np.asarray([], dtype=np.float64),
            fit,
            pressure_unit,
            [fit.message],
            metadata,
        )

    pressure_energy_density = model.pressure(
        model_spec,
        fit.parameters,
        volume_array,
    )
    pressure = np.asarray(
        energy_to_pressure(
            pressure_energy_density,
            energy_unit,
            volume_unit,
            pressure_unit,
        ),
        dtype=np.float64,
    )
    warnings: list[str] = []
    if fit.quality is FitQuality.POOR:
        warnings.append("the EOS pressure estimate is based on a poor fit")
    warnings.extend(fit.warnings)
    return PressureEstimate(
        method="eos",
        pressure=pressure,
        fit=fit,
        unit=pressure_unit,
        warnings=warnings,
        metadata=metadata,
    )


def _energy_volume_arrays(
    volume: ArrayLike,
    energy: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return validated one-dimensional energy--volume arrays."""
    volume_array, energy_array = validate_xy(
        np.asarray(volume, dtype=np.float64),
        np.asarray(energy, dtype=np.float64),
    )
    if volume_array.ndim != 1:
        raise ValueError("energy-volume pressure estimation requires volume values")
    if np.any(volume_array <= 0.0):
        raise ValueError("energy-volume pressure estimation requires positive volumes")
    return volume_array, energy_array


def _failed_estimate(
    method: str,
    message: str,
    *,
    unit: str,
    metadata: Mapping[str, Any] | None = None,
) -> PressureEstimate:
    """Create a failed pressure estimate with retained diagnostics."""
    details = dict(metadata or {})
    return PressureEstimate(
        method=method,
        pressure=np.asarray([], dtype=np.float64),
        fit=FitResult.failed(
            message,
            status=FitStatus.FAILED,
            metadata=details,
        ),
        unit=unit,
        warnings=[message],
        metadata=details,
    )


__all__ = [
    "PressureEstimate",
    "pressure_from_energy_eos",
    "pressure_from_energy_polynomial",
]
