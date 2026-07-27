# -*- coding: utf-8 -*-

"""Result-aware public plot inventory for harmonic thermodynamics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from quantas.models import (
    PlotContextDescriptor,
    PlotInventory,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
)
from quantas.modules.ha.io.export import convert_property_values
from quantas.modules.ha.models import HAResult
from quantas.modules.ha.plot.spec import _property_array, _temperature_array


@dataclass(frozen=True, slots=True)
class _HAPropertyDefinition:
    """Static scientific metadata for one harmonic property."""

    key: str
    name: str
    symbol_math: str
    symbol_plain: str
    description: str
    category: str


_PROPERTY_DEFINITIONS = (
    _HAPropertyDefinition(
        "static_energy",
        "Static energy",
        "U_0",
        "U₀",
        "Electronic static energy at each sampled volume.",
        "energy",
    ),
    _HAPropertyDefinition(
        "zero_point_energy",
        "Zero-point energy",
        r"U_{\mathrm{ZP}}",
        "U_ZP",
        "Temperature-independent vibrational zero-point energy.",
        "energy",
    ),
    _HAPropertyDefinition(
        "thermal_energy",
        "Thermal energy",
        r"U_{\mathrm{th}}",
        "U_th",
        "Thermal vibrational contribution to the internal energy.",
        "energy",
    ),
    _HAPropertyDefinition(
        "internal_energy",
        "Internal energy",
        "U",
        "U",
        "Total static, zero-point, and thermal internal energy.",
        "energy",
    ),
    _HAPropertyDefinition(
        "entropy",
        "Entropy",
        "S",
        "S",
        "Harmonic vibrational entropy.",
        "entropy",
    ),
    _HAPropertyDefinition(
        "vibrational_free_energy",
        "Vibrational Helmholtz free energy",
        r"F_{\mathrm{vib}}",
        "F_vib",
        "Vibrational contribution to the Helmholtz free energy.",
        "energy",
    ),
    _HAPropertyDefinition(
        "free_energy",
        "Helmholtz free energy",
        "F",
        "F",
        "Total harmonic Helmholtz free energy.",
        "energy",
    ),
    _HAPropertyDefinition(
        "isochoric_heat_capacity",
        "Isochoric heat capacity",
        "C_V",
        "Cᵥ",
        "Heat capacity at constant volume.",
        "heat_capacity",
    ),
)


def describe_ha_plots(result: HAResult) -> PlotInventory:
    """Describe harmonic plots buildable from one result.

    The current HA builder emits temperature curves and includes every sampled
    volume as a separate series. Temperature and volume grids are therefore
    exposed as informative, exact result context rather than independent
    selections in this increment.
    """
    warnings: list[str] = []
    try:
        temperature = _temperature_array(result)
    except ValueError as exc:
        return PlotInventory(
            module="ha",
            properties=(),
            representations=(),
            warnings=(str(exc),),
        )

    available: list[tuple[_HAPropertyDefinition, str | None]] = []
    for definition in _PROPERTY_DEFINITIONS:
        try:
            values = _property_array(result, definition.key, temperature.size)
            _, unit = convert_property_values(
                value=values,
                attr=definition.key,
                units=result.metadata.get("units", {}),
                normalization=result.metadata.get("normalization"),
            )
        except (ValueError, KeyError) as exc:
            if getattr(result, definition.key) is not None:
                warnings.append(f"{definition.key}: {exc}")
            continue
        available.append((definition, _normalize_unit(unit)))

    property_keys = tuple(definition.key for definition, _ in available)
    if not property_keys:
        return PlotInventory(
            module="ha",
            properties=(),
            representations=(),
            warnings=tuple(warnings or ["No harmonic plot properties are available."]),
        )

    contexts = [
        PlotContextDescriptor(
            key="temperature_grid",
            name="Temperature grid",
            description="Exact temperature coordinates stored in the HA result.",
            values=tuple(float(value) for value in temperature),
            unit=_metadata_unit(result, "temperature", default="K"),
            selectable=False,
        )
    ]
    volume_values, volume_unit = _volume_context(result, property_keys)
    contexts.append(
        PlotContextDescriptor(
            key="sampled_volume",
            name="Sampled volume",
            description=(
                "Volumes represented as separate series by the current HA "
                "temperature-curve builder."
            ),
            values=volume_values,
            unit=volume_unit,
            selectable=False,
        )
    )

    representation = PlotRepresentationDescriptor(
        key="temperature_curves",
        name="Thermodynamic temperature curves",
        plot_kind="line",
        description=(
            "One curve for each sampled volume, with temperature on the "
            "independent axis."
        ),
        property_keys=property_keys,
        supported_contexts=("temperature_grid", "sampled_volume"),
        constraints=(
            "The current public builder emits all sampled volumes.",
            "Temperature-independent energies are broadcast across the temperature grid.",
        ),
    )
    properties = tuple(
        PlotPropertyDescriptor(
            key=definition.key,
            name=definition.name,
            symbol_math=definition.symbol_math,
            symbol_plain=definition.symbol_plain,
            unit=unit,
            description=definition.description,
            category=definition.category,
            representations=("temperature_curves",),
        )
        for definition, unit in available
    )
    return PlotInventory(
        module="ha",
        properties=properties,
        representations=(representation,),
        contexts=tuple(contexts),
        warnings=tuple(warnings),
    )


def _metadata_unit(result: HAResult, key: str, *, default: str) -> str | None:
    """Return one normalized unit from result metadata."""
    value = str(result.metadata.get("units", {}).get(key, default))
    return _normalize_unit(value)


def _normalize_unit(unit: str | None) -> str | None:
    """Represent unknown and dimensionless unit labels as ``None``."""
    if unit is None:
        return None
    value = str(unit).strip()
    if not value or value.lower() in {"unknown", "dimensionless", "none"}:
        return None
    return value


def _volume_context(
    result: HAResult,
    property_keys: tuple[str, ...],
) -> tuple[tuple[float | int, ...], str | None]:
    """Return exact sampled volumes or stable curve indices."""
    if result.volume is not None:
        volume = np.asarray(result.volume, dtype=np.float64)
        if volume.ndim == 1:
            return (
                tuple(float(value) for value in volume),
                _metadata_unit(result, "volume", default="unknown"),
            )

    nseries = 0
    temperature = _temperature_array(result)
    for key in property_keys:
        values = _property_array(result, key, temperature.size)
        nseries = max(nseries, int(values.shape[1]))
    return tuple(range(nseries)), None


__all__ = ["describe_ha_plots"]
