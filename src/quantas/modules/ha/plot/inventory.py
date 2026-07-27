# -*- coding: utf-8 -*-

"""Result-aware public plot inventory for harmonic thermodynamics."""

from __future__ import annotations

import numpy as np

from quantas.models import (
    PlotContextDescriptor,
    PlotInventory,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
)
from quantas.modules.ha.io.export import convert_property_values
from quantas.modules.ha.models import HAResult
from quantas.modules.ha.plot.labels import (
    HAPlotProperty,
    available_plot_properties,
)
from quantas.modules.ha.plot.spec import _property_array, _temperature_array


def describe_ha_plots(result: HAResult) -> PlotInventory:
    """Describe exact-grid HA sections and maps buildable from one result.

    HA properties are naturally defined on a temperature-volume grid.  The
    standard representation keeps temperature on the independent axis and one
    curve per sampled volume.  When physical volume coordinates match a
    property grid, the same stored data also support volume sections at exact
    stored temperatures and a volume-temperature contour map.  Discovery never
    implies interpolation.
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

    available: list[tuple[HAPlotProperty, str | None, tuple[str, ...]]] = []
    volume_capable_keys: list[str] = []
    contour_capable_keys: list[str] = []
    physical_volume: np.ndarray | None = None

    for property_info in available_plot_properties().values():
        try:
            values = _property_array(
                result,
                property_info.attribute,
                temperature.size,
            )
            _, unit = convert_property_values(
                value=values,
                attr=property_info.attribute,
                units=result.metadata.get("units", {}),
                normalization=result.metadata.get("normalization"),
            )
        except (ValueError, KeyError) as exc:
            if getattr(result, property_info.attribute) is not None:
                warnings.append(f"{property_info.attribute}: {exc}")
            continue

        representations = ["temperature_curves"]
        volume = _matched_volume(result, values.shape[1])
        if volume is not None:
            physical_volume = volume
            if volume.size >= 2:
                representations.append("volume_curves")
                volume_capable_keys.append(property_info.attribute)
            if temperature.size >= 2 and volume.size >= 2:
                representations.append("volume_temperature_contour")
                contour_capable_keys.append(property_info.attribute)
        available.append(
            (property_info, _normalize_unit(unit), tuple(representations))
        )

    property_keys = tuple(item[0].attribute for item in available)
    if not property_keys:
        return PlotInventory(
            module="ha",
            properties=(),
            representations=(),
            warnings=tuple(warnings or ["No harmonic plot properties are available."]),
        )

    contexts: list[PlotContextDescriptor] = [
        PlotContextDescriptor(
            key="curve_axis",
            name="Curve axis",
            description=(
                "Natural variable used as the independent axis of line sections."
            ),
            values=("temperature", "volume")
            if volume_capable_keys
            else ("temperature",),
            default="temperature",
        ),
        PlotContextDescriptor(
            key="temperature_grid",
            name="Temperature grid",
            description=(
                "Exact native temperature coordinates. They may be selected "
                "for volume sections without interpolation."
            ),
            values=tuple(float(value) for value in temperature),
            unit=_metadata_unit(result, "temperature", default="K"),
            selectable=True,
        ),
    ]

    if physical_volume is not None:
        volume_values: tuple[float | int, ...] = tuple(
            float(value) for value in physical_volume
        )
        volume_unit = _metadata_unit(result, "volume", default="unknown")
        volume_selectable = True
    else:
        volume_values, volume_unit = _volume_context(result, property_keys)
        volume_selectable = False
        if volume_values:
            warnings.append(
                "Physical sampled-volume coordinates are unavailable or do not "
                "match every property grid; only temperature curves are exposed "
                "for incompatible properties."
            )

    contexts.append(
        PlotContextDescriptor(
            key="sampled_volume",
            name="Sampled volume",
            description=(
                "Exact native volume coordinates used for volume sections and "
                "maps when available; otherwise stable curve indices."
            ),
            values=volume_values,
            unit=volume_unit,
            selectable=volume_selectable,
        )
    )

    representations: list[PlotRepresentationDescriptor] = [
        PlotRepresentationDescriptor(
            key="temperature_curves",
            name="Temperature sections",
            plot_kind="line",
            description=(
                "Property as a function of temperature, with one curve for each "
                "selected sampled volume."
            ),
            property_keys=property_keys,
            supported_contexts=("curve_axis", "sampled_volume"),
            constraints=(
                "Selected volumes must be exact values from the native grid.",
                "Temperature-independent energies are broadcast across the "
                "temperature grid.",
            ),
        )
    ]
    if volume_capable_keys:
        representations.append(
            PlotRepresentationDescriptor(
                key="volume_curves",
                name="Volume sections",
                plot_kind="line",
                description=(
                    "Property as a function of sampled volume, with one curve "
                    "for each selected stored temperature."
                ),
                property_keys=tuple(volume_capable_keys),
                supported_contexts=("curve_axis", "temperature_grid"),
                constraints=(
                    "Selected temperatures must be exact values from the native grid.",
                    "No interpolation or volume resampling is performed.",
                ),
            )
        )
    if contour_capable_keys:
        representations.append(
            PlotRepresentationDescriptor(
                key="volume_temperature_contour",
                name="Volume-temperature map",
                plot_kind="contour",
                description=(
                    "Property on the complete native temperature-volume grid."
                ),
                property_keys=tuple(contour_capable_keys),
                supported_contexts=("temperature_grid", "sampled_volume"),
                constraints=(
                    "At least two temperatures and two sampled volumes are required.",
                    "No interpolation or resampling is performed by Quantas.",
                ),
            )
        )

    properties = tuple(
        PlotPropertyDescriptor(
            key=property_info.attribute,
            name=property_info.name,
            symbol_math=property_info.symbol_math,
            symbol_plain=property_info.symbol_plain,
            unit=unit,
            description=property_info.description,
            category=property_info.category,
            representations=representations,
        )
        for property_info, unit, representations in available
    )
    return PlotInventory(
        module="ha",
        properties=properties,
        representations=tuple(representations),
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


def _matched_volume(result: HAResult, nseries: int) -> np.ndarray | None:
    """Return physical volumes when they match one property grid."""
    if result.volume is None:
        return None
    volume = np.asarray(result.volume, dtype=np.float64)
    if volume.ndim != 1 or volume.size != nseries:
        return None
    if np.unique(volume).size != volume.size:
        return None
    return volume


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
