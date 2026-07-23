# -*- coding: utf-8 -*-

"""Dataset-coordinate classification for EOS workflows."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .models import (
    EOSCoordinateProfile,
    EOSCoordinateVariation,
    EOSDatasetClassification,
)

if TYPE_CHECKING:
    from .models import EOSDataset


_VALUE_COORDINATES = (
    "pressure",
    "temperature",
    "volume",
    "energy",
    "a",
    "b",
    "c",
    "alpha",
    "beta",
    "gamma",
)


def coordinate_profile(
    dataset: EOSDataset,
    name: str,
    *,
    mask: np.ndarray | None = None,
    relative_tolerance: float = 1.0e-12,
    absolute_tolerance: float | None = None,
) -> EOSCoordinateProfile:
    """Classify one measured coordinate as constant or variable.

    Parameters
    ----------
    dataset : EOSDataset
        Dataset containing the canonical measured coordinate.
    name : str
        Canonical non-uncertainty column name.
    mask : ndarray or None, optional
        Boolean observation selection.
    relative_tolerance : float, optional
        Relative constancy tolerance.
    absolute_tolerance : float or None, optional
        Absolute tolerance. If omitted, a scale-aware floating-point tolerance
        is derived from the selected values.

    Returns
    -------
    EOSCoordinateProfile
        Coordinate range and variation state.

    Raises
    ------
    KeyError
        If the coordinate is unavailable.
    ValueError
        If a sigma column, invalid tolerance, or invalid mask is supplied.
    """
    if name.startswith("sigma_"):
        raise ValueError(
            "EOS coordinate profiles require measured values, not sigma columns."
        )
    values = dataset.columns[name]
    selected = _selected_values(values, mask)
    if relative_tolerance < 0.0:
        raise ValueError("relative_tolerance cannot be negative")
    scale = max(1.0, float(np.max(np.abs(selected))))
    resolved_atol = float(
        64.0 * np.finfo(np.float64).eps * scale
        if absolute_tolerance is None
        else absolute_tolerance
    )
    if resolved_atol < 0.0:
        raise ValueError("absolute_tolerance cannot be negative")
    reference = float(np.mean(selected, dtype=np.float64))
    constant = bool(
        np.allclose(
            selected,
            reference,
            rtol=float(relative_tolerance),
            atol=resolved_atol,
        )
    )
    minimum = float(np.min(selected))
    maximum = float(np.max(selected))
    return EOSCoordinateProfile(
        name=name,
        variation=(
            EOSCoordinateVariation.CONSTANT
            if constant
            else EOSCoordinateVariation.VARIABLE
        ),
        minimum=minimum,
        maximum=maximum,
        span=maximum - minimum,
        reference_value=reference if constant else None,
        unit=dataset.units.get(name),
        npoints=int(selected.size),
        absolute_tolerance=resolved_atol,
        relative_tolerance=float(relative_tolerance),
    )


def classify_dataset(
    dataset: EOSDataset,
    *,
    mask: np.ndarray | None = None,
    relative_tolerance: float = 1.0e-12,
    absolute_tolerance: float | None = None,
) -> EOSDatasetClassification:
    """Classify all measured coordinates and identify reference conditions.

    Parameters
    ----------
    dataset : EOSDataset
        Source dataset.
    mask : ndarray or None, optional
        Boolean observation selection.
    relative_tolerance, absolute_tolerance : float or None, optional
        Constancy tolerances forwarded to :func:`coordinate_profile`.

    Returns
    -------
    EOSDatasetClassification
        Coordinate profiles and constant pressure/temperature conditions.
    """
    profiles = {
        name: coordinate_profile(
            dataset,
            name,
            mask=mask,
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )
        for name in _VALUE_COORDINATES
        if name in dataset.columns
    }
    variable = tuple(name for name, profile in profiles.items() if profile.is_variable)
    constant = tuple(name for name, profile in profiles.items() if profile.is_constant)
    pressure = profiles.get("pressure")
    temperature = profiles.get("temperature")
    return EOSDatasetClassification(
        profiles=profiles,
        variable_coordinates=variable,
        constant_coordinates=constant,
        is_isobaric=pressure is not None and pressure.is_constant,
        is_isothermal=temperature is not None and temperature.is_constant,
        reference_pressure=(
            pressure.reference_value
            if pressure is not None and pressure.is_constant
            else None
        ),
        reference_temperature=(
            temperature.reference_value
            if temperature is not None and temperature.is_constant
            else None
        ),
    )


def require_variable_coordinate(
    dataset: EOSDataset,
    name: str,
    *,
    purpose: str,
    mask: np.ndarray | None = None,
) -> EOSCoordinateProfile:
    """Return a profile or reject a constant fitting coordinate.

    Parameters
    ----------
    dataset : EOSDataset
        Source dataset.
    name : str
        Canonical measured coordinate.
    purpose : str
        Human-readable scientific operation used in the error message.
    mask : ndarray or None, optional
        Boolean observation selection.

    Returns
    -------
    EOSCoordinateProfile
        Variable coordinate profile.

    Raises
    ------
    ValueError
        If the selected coordinate is constant within numerical tolerance.
    """
    profile = coordinate_profile(dataset, name, mask=mask)
    if profile.is_constant:
        reference = profile.reference_value
        unit = f" {profile.unit}" if profile.unit else ""
        raise ValueError(
            f"The {name} column is constant within tolerance "
            f"({reference:.12g}{unit}) and cannot constrain {purpose}."
        )
    return profile


def _selected_values(
    values: np.ndarray,
    mask: np.ndarray | None,
) -> np.ndarray:
    """Return selected values after validating an optional boolean mask."""
    if mask is None:
        return values
    selected_mask = np.asarray(mask, dtype=bool)
    if selected_mask.ndim != 1 or selected_mask.shape != values.shape:
        raise ValueError("EOS coordinate mask must match the data shape.")
    if not np.any(selected_mask):
        raise ValueError("EOS coordinate mask must select at least one observation.")
    return values[selected_mask]


__all__ = [
    "classify_dataset",
    "coordinate_profile",
    "require_variable_coordinate",
]
