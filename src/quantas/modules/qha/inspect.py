# -*- coding: utf-8 -*-

"""Inspection utilities for quasi-harmonic input data.

This module evaluates the static energy-volume dataset before a full QHA
calculation is started.  It provides pressure estimates from polynomial and
energy equation-of-state fits and returns structured diagnostics that can be
rendered by the command-line interface or by a graphical frontend.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np

from quantas.core.math.fitting import validate_xy
from quantas.core.physics.eos import (
    PressureEstimate,
    pressure_from_energy_eos,
    pressure_from_energy_polynomial,
)
from quantas.modules.qha.models import QHAInput, QHAOptions

ArrayLike = np.ndarray | Sequence[float]


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
        poly = (
            None
            if self.polynomial is None
            or not self.polynomial.success
            or self.polynomial.pressure.shape != self.volume.shape
            else self.polynomial.pressure
        )
        eos = (
            None
            if self.eos is None
            or not self.eos.success
            or self.eos.pressure.shape != self.volume.shape
            else self.eos.pressure
        )
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
        polynomial = pressure_from_energy_polynomial(
            volume,
            energy,
            degree=degree,
            energy_unit=options.energy_unit,
            volume_unit=options.volume_unit,
            pressure_unit=options.pressure_unit,
        )
        warnings_.extend(polynomial.warnings)

    if include_eos:
        eos_estimate = pressure_from_energy_eos(
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
