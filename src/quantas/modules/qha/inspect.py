# -*- coding: utf-8 -*-

"""Inspection utilities for quasi-harmonic input data.

This module evaluates the static energy-volume dataset before a full QHA
calculation is started.  It provides pressure estimates from polynomial and
energy equation-of-state fits and returns structured diagnostics that can be
rendered by the command-line interface or by a graphical frontend.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np

from quantas.core.physics.units import energy_to_pressure
from quantas.core.physics.eos import EnergyEOS
from quantas.core.math.fitting import FitQuality, FitResult, FitStatus, validate_xy
from quantas.core.math.derivative import polynomial_derivative_from_coefficients
from quantas.core.math.polynomials import fit_polynomial_result as polynomial_fit
from quantas.modules.qha.models import QHAInput, QHAOptions

ArrayLike = np.ndarray | Sequence[float]


@dataclass(slots=True)
class PressureEstimate:
    """Pressure estimates associated with a single fitting method.

    Parameters
    ----------
    method : str
        Name of the method used to estimate the pressure values.
    pressure : ndarray
        Pressure values evaluated at the input volumes.
    fit : FitResult
        Fit diagnostics associated with the estimate.
    unit : str
        Pressure unit used for the returned values.
    warnings : list of str
        Non-fatal diagnostic messages.
    metadata : dict
        Additional method-specific information.
    """

    method: str
    pressure: np.ndarray
    fit: FitResult
    unit: str
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def success(self) -> bool:
        """Return whether the pressure estimate is usable.

        Returns
        -------
        bool
            ``True`` when the underlying fit succeeded and pressure values are
            finite.
        """
        return bool(self.fit.success and np.all(np.isfinite(self.pressure)))

    @property
    def pressure_min(self) -> float | None:
        """Return the minimum pressure value.

        Returns
        -------
        float or None
            Minimum pressure, or ``None`` when no values are available.
        """
        if self.pressure.size == 0:
            return None
        return float(np.nanmin(self.pressure))

    @property
    def pressure_max(self) -> float | None:
        """Return the maximum pressure value.

        Returns
        -------
        float or None
            Maximum pressure, or ``None`` when no values are available.
        """
        if self.pressure.size == 0:
            return None
        return float(np.nanmax(self.pressure))

    def as_dict(self) -> dict[str, Any]:
        """Return the pressure estimate as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the pressure estimate.
        """
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


@dataclass(slots=True)
class PressureVolumePreview:
    """Pressure-volume preview of a QHA input dataset.

    Parameters
    ----------
    volume : ndarray
        Input unit-cell volumes.
    energy : ndarray
        Input static energies.
    pressure_unit : str
        Pressure unit used by the pressure estimates.
    polynomial : PressureEstimate or None
        Pressure estimate from the polynomial fit.
    eos : PressureEstimate or None
        Pressure estimate from the energy equation-of-state fit.
    warnings : list of str
        Non-fatal diagnostic messages collected during inspection.
    metadata : dict
        Additional information associated with the preview.
    """

    volume: np.ndarray
    energy: np.ndarray
    pressure_unit: str
    polynomial: PressureEstimate | None = None
    eos: PressureEstimate | None = None
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def success(self) -> bool:
        """Return whether at least one pressure estimate succeeded.

        Returns
        -------
        bool
            ``True`` when a polynomial or EOS pressure estimate is usable.
        """
        return bool(
            (self.polynomial is not None and self.polynomial.success)
            or (self.eos is not None and self.eos.success)
        )

    def table_rows(self) -> list[dict[str, float | None]]:
        """Return volume-pressure values as neutral table rows.

        Returns
        -------
        list of dict
            Rows containing volume, energy, polynomial pressure and EOS
            pressure.  Missing estimates are represented by ``None``.
        """
        poly = None if self.polynomial is None else self.polynomial.pressure
        eos = None if self.eos is None else self.eos.pressure
        rows: list[dict[str, float | None]] = []
        for index, (volume, energy) in enumerate(
            zip(self.volume, self.energy, strict=True)
        ):
            rows.append(
                {
                    "volume": float(volume),
                    "energy": float(energy),
                    "pressure_polynomial": None if poly is None else float(poly[index]),
                    "pressure_eos": None if eos is None else float(eos[index]),
                }
            )
        return rows

    def as_dict(self) -> dict[str, Any]:
        """Return the preview as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the pressure-volume preview.
        """
        return {
            "success": self.success,
            "volume": self.volume.tolist(),
            "energy": self.energy.tolist(),
            "pressure_unit": self.pressure_unit,
            "polynomial": None
            if self.polynomial is None
            else self.polynomial.as_dict(),
            "eos": None if self.eos is None else self.eos.as_dict(),
            "rows": self.table_rows(),
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


def _failed_estimate(
    method: str,
    message: str,
    *,
    unit: str,
    metadata: Mapping[str, Any] | None = None,
) -> PressureEstimate:
    """Create a failed pressure estimate.

    Parameters
    ----------
    method : str
        Name of the method that failed.
    message : str
        Failure message.
    unit : str
        Pressure unit requested by the caller.
    metadata : mapping, optional
        Additional diagnostic information.

    Returns
    -------
    PressureEstimate
        Failed pressure estimate with an empty pressure array.
    """
    return PressureEstimate(
        method=method,
        pressure=np.asarray([], dtype=np.float64),
        fit=FitResult.failed(
            message, status=FitStatus.FAILED, metadata=dict(metadata or {})
        ),
        unit=unit,
        warnings=[message],
        metadata=dict(metadata or {}),
    )


def _polynomial_pressure(
    volume: np.ndarray,
    energy: np.ndarray,
    *,
    degree: int,
    energy_unit: str,
    volume_unit: str,
    pressure_unit: str,
) -> PressureEstimate:
    """Estimate pressure from a polynomial energy-volume fit.

    Parameters
    ----------
    volume : ndarray
        Unit-cell volumes.
    energy : ndarray
        Static energies.
    degree : int
        Polynomial degree.
    energy_unit : str
        Energy unit of the static energies.
    volume_unit : str
        Length unit defining the volume unit.
    pressure_unit : str
        Requested pressure unit.

    Returns
    -------
    PressureEstimate
        Polynomial pressure estimate.
    """
    fit = polynomial_fit(volume, energy, degree)
    if not fit.success or fit.parameters is None:
        return PressureEstimate(
            "polynomial",
            np.asarray([], dtype=np.float64),
            fit,
            pressure_unit,
            [fit.message],
        )

    pressure_energy_density = -polynomial_derivative_from_coefficients(
        fit.parameters,
        volume,
    )
    pressure = np.asarray(
        energy_to_pressure(
            pressure_energy_density, energy_unit, volume_unit, pressure_unit
        ),
        dtype=np.float64,
    )
    warnings_: list[str] = []
    if fit.quality is FitQuality.POOR:
        warnings_.append("the polynomial pressure estimate is based on a poor fit")
    warnings_.extend(fit.warnings)
    return PressureEstimate(
        method="polynomial",
        pressure=pressure,
        fit=fit,
        unit=pressure_unit,
        warnings=warnings_,
        metadata={"degree": int(degree)},
    )


def _eos_pressure(
    volume: np.ndarray,
    energy: np.ndarray,
    *,
    eos: str,
    energy_unit: str,
    volume_unit: str,
    pressure_unit: str,
    maxfev: int | None = None,
) -> PressureEstimate:
    """Estimate pressure from an energy equation-of-state fit.

    Parameters
    ----------
    volume : ndarray
        Unit-cell volumes.
    energy : ndarray
        Static energies.
    eos : str
        Equation-of-state model.
    energy_unit : str
        Energy unit of the static energies.
    volume_unit : str
        Length unit defining the volume unit.
    pressure_unit : str
        Requested pressure unit.
    maxfev : int, optional
        Maximum number of optimizer evaluations.

    Returns
    -------
    PressureEstimate
        EOS pressure estimate.
    """
    model = EnergyEOS()
    try:
        model_spec = model.model(eos)
    except ValueError as exc:
        return _failed_estimate(
            "eos", str(exc), unit=pressure_unit, metadata={"eos": eos}
        )

    fit = model.fit(model_spec, volume, energy, maxfev=maxfev)
    if not fit.success or fit.parameters is None:
        return PressureEstimate(
            "eos",
            np.asarray([], dtype=np.float64),
            fit,
            pressure_unit,
            [fit.message],
            {"eos": model_spec.tag},
        )

    pressure_energy_density = model.pressure(model_spec, fit.parameters, volume)
    pressure = np.asarray(
        energy_to_pressure(
            pressure_energy_density, energy_unit, volume_unit, pressure_unit
        ),
        dtype=np.float64,
    )
    warnings_: list[str] = []
    if fit.quality is FitQuality.POOR:
        warnings_.append("the EOS pressure estimate is based on a poor fit")
    warnings_.extend(fit.warnings)
    return PressureEstimate(
        method="eos",
        pressure=pressure,
        fit=fit,
        unit=pressure_unit,
        warnings=warnings_,
        metadata={
            "eos": model_spec.tag,
            "eos_family": model_spec.family.value,
            "eos_order": model_spec.order,
        },
    )


def pressure_volume_preview(
    qha_input: QHAInput,
    options: QHAOptions | None = None,
    *,
    include_polynomial: bool = True,
    include_eos: bool = True,
    polynomial_degree: int | None = None,
    eos: str | None = None,
    maxfev: int | None = None,
) -> PressureVolumePreview:
    """Estimate the pressure range sampled by a QHA input dataset.

    Parameters
    ----------
    qha_input : QHAInput
        Normalized QHA input data.
    options : QHAOptions, optional
        Calculation options providing units, polynomial degree and default EOS.
    include_polynomial : bool, optional
        If ``True``, include the polynomial pressure estimate.
    include_eos : bool, optional
        If ``True``, include the EOS pressure estimate.
    polynomial_degree : int, optional
        Polynomial degree.  If omitted, ``options.energy_degree`` is used.
    eos : str, optional
        EOS name.  If omitted, ``options.eos`` is used.
    maxfev : int, optional
        Maximum number of optimizer evaluations for the EOS fit.

    Returns
    -------
    PressureVolumePreview
        Structured pressure-volume preview with fit diagnostics.

    Raises
    ------
    ValueError
        If input volume and static energy arrays are missing or inconsistent.
    """
    options = QHAOptions() if options is None else options
    degree = (
        options.energy_degree if polynomial_degree is None else int(polynomial_degree)
    )
    eos_name = options.eos if eos is None else eos

    qha_input.validate_shapes()
    if qha_input.volume is None or qha_input.energy is None:
        raise ValueError("QHA input requires volume and static energy arrays")
    volume, energy = validate_xy(qha_input.volume, qha_input.energy)

    warnings_: list[str] = []
    polynomial = None
    eos_estimate = None

    if include_polynomial:
        polynomial = _polynomial_pressure(
            volume,
            energy,
            degree=degree,
            energy_unit=options.energy_unit,
            volume_unit=options.volume_unit,
            pressure_unit=options.pressure_unit,
        )
        warnings_.extend(polynomial.warnings)

    if include_eos:
        eos_estimate = _eos_pressure(
            volume,
            energy,
            eos=eos_name,
            energy_unit=options.energy_unit,
            volume_unit=options.volume_unit,
            pressure_unit=options.pressure_unit,
            maxfev=maxfev,
        )
        warnings_.extend(eos_estimate.warnings)

    return PressureVolumePreview(
        volume=volume,
        energy=energy,
        pressure_unit=options.pressure_unit,
        polynomial=polynomial,
        eos=eos_estimate,
        warnings=warnings_,
        metadata={
            "energy_unit": options.energy_unit,
            "volume_unit": options.volume_unit,
            "pressure_unit": options.pressure_unit,
            "polynomial_degree": degree,
            "eos": eos_name,
        },
    )
