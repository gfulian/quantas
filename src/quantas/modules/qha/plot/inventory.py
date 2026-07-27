# -*- coding: utf-8 -*-

"""Result-aware public plot inventory for quasi-harmonic thermodynamics."""

from __future__ import annotations

import numpy as np

from quantas.models import (
    PlotContextDescriptor,
    PlotInventory,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
)
from quantas.modules.qha.models import QHAResult
from quantas.modules.qha.plot.labels import available_plot_properties


def describe_qha_plots(result: QHAResult) -> PlotInventory:
    """Describe exact-grid QHA sections and maps buildable from one result.

    QHA scalar thermodynamic properties are represented on the native
    temperature-pressure grid.  The default line family uses temperature as
    the independent variable and pressure as the selected condition.  The
    complementary family uses pressure as the independent variable and exact
    stored temperatures as the selected condition.  Contour maps use the full
    native grid.  Discovery never implies interpolation.
    """
    warnings: list[str] = []
    try:
        temperature = _required_grid(result.temperature, "temperature")
        pressure = _required_grid(result.pressure, "pressure")
    except ValueError as exc:
        return PlotInventory(
            module="qha",
            properties=(),
            representations=(),
            warnings=(str(exc),),
        )

    line_representations = (
        ("temperature_curves", "pressure_curves")
        if pressure.size >= 2
        else ("temperature_curves",)
    )
    scalar_representations = (
        (*line_representations, "pressure_temperature_contour")
        if temperature.size >= 2 and pressure.size >= 2
        else line_representations
    )

    descriptors: list[PlotPropertyDescriptor] = []
    valid_attributes: set[str] = set()
    contour_keys: list[str] = []
    for prop in available_plot_properties().values():
        values = getattr(result, prop.attribute, None)
        if values is None:
            continue
        try:
            _as_grid(
                np.asarray(values, dtype=np.float64),
                ntemperature=temperature.size,
                npressure=pressure.size,
                name=prop.attribute,
            )
        except ValueError as exc:
            warnings.append(str(exc))
            continue

        key = prop.attribute
        valid_attributes.add(prop.attribute)
        if "pressure_temperature_contour" in scalar_representations:
            contour_keys.append(key)
        descriptors.append(
            PlotPropertyDescriptor(
                key=key,
                name=prop.description,
                symbol_math=_math_symbol(prop.symbol),
                symbol_plain=prop.symbol_plain,
                unit=_property_unit(result, prop.unit_key),
                description=_property_description(prop.attribute, prop.description),
                category=_property_category(prop.unit_key),
                representations=scalar_representations,
            )
        )

    heat_components = tuple(
        component
        for component, attribute in (
            ("isobaric", "isobaric_heat_capacity"),
            ("isochoric", "isochoric_heat_capacity"),
        )
        if attribute in valid_attributes
    )
    if heat_components:
        descriptors.append(
            PlotPropertyDescriptor(
                key="heat_capacities",
                name="Heat-capacity comparison",
                symbol_math="C_P, C_V",
                symbol_plain="Cₚ/Cᵥ",
                unit=_metadata_unit(
                    result,
                    "heat_capacity",
                    fallback=_metadata_unit(result, "entropy", fallback=None),
                ),
                description=(
                    "Combined isobaric and isochoric heat-capacity sections "
                    "when either component is available."
                ),
                category="comparison",
                components=heat_components,
                representations=line_representations,
            )
        )

    if not descriptors:
        return PlotInventory(
            module="qha",
            properties=(),
            representations=(),
            contexts=_contexts(result, temperature, pressure),
            warnings=tuple(warnings or ["No QHA plot properties are available."]),
        )

    all_line_keys = tuple(item.key for item in descriptors)
    representations: list[PlotRepresentationDescriptor] = [
        PlotRepresentationDescriptor(
            key="temperature_curves",
            name="Temperature sections",
            plot_kind="line",
            description=(
                "Property as a function of temperature, with one curve for "
                "each selected pressure."
            ),
            property_keys=all_line_keys,
            supported_contexts=("curve_axis", "pressure_grid"),
            constraints=(
                "Selected pressures must be exact values from the native grid.",
                "No interpolation is performed by Quantas.",
            ),
        )
    ]
    if pressure.size >= 2:
        representations.append(
            PlotRepresentationDescriptor(
                key="pressure_curves",
                name="Pressure sections",
                plot_kind="line",
                description=(
                    "Property as a function of pressure, with one curve for "
                    "each selected stored temperature."
                ),
                property_keys=all_line_keys,
                supported_contexts=("curve_axis", "temperature_grid"),
                constraints=(
                    "At least two stored pressures are required.",
                    "Selected temperatures must be exact values from the native grid.",
                    "No interpolation is performed by Quantas.",
                ),
            )
        )
    if contour_keys:
        representations.append(
            PlotRepresentationDescriptor(
                key="pressure_temperature_contour",
                name="Pressure-temperature map",
                plot_kind="contour",
                description=(
                    "One scalar property on the complete native pressure-"
                    "temperature grid."
                ),
                property_keys=tuple(contour_keys),
                supported_contexts=("temperature_grid", "pressure_grid"),
                constraints=(
                    "At least two temperatures and two pressures are required.",
                    "Combined heat-capacity comparison is a line representation only.",
                    "No interpolation or resampling is performed by Quantas.",
                ),
            )
        )

    return PlotInventory(
        module="qha",
        properties=tuple(descriptors),
        representations=tuple(representations),
        contexts=_contexts(result, temperature, pressure),
        warnings=tuple(warnings),
    )


def _contexts(
    result: QHAResult,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> tuple[PlotContextDescriptor, ...]:
    """Return exact native QHA coordinate and line-axis contexts."""
    return (
        PlotContextDescriptor(
            key="curve_axis",
            name="Curve axis",
            description="Natural variable used as the independent line axis.",
            values=("temperature", "pressure")
            if pressure.size >= 2
            else ("temperature",),
            default="temperature",
        ),
        PlotContextDescriptor(
            key="temperature_grid",
            name="Temperature grid",
            description=(
                "Exact native temperature coordinates available for pressure "
                "sections."
            ),
            values=tuple(float(value) for value in temperature),
            unit=_metadata_unit(result, "temperature", fallback="K"),
            selectable=True,
        ),
        PlotContextDescriptor(
            key="pressure_grid",
            name="Pressure grid",
            description=(
                "Exact native pressure coordinates available for temperature "
                "sections."
            ),
            values=tuple(float(value) for value in pressure),
            unit=_metadata_unit(result, "pressure", fallback="GPa"),
            selectable=True,
        ),
    )


def _required_grid(values: np.ndarray | None, name: str) -> np.ndarray:
    """Return one non-empty, unique, one-dimensional native coordinate grid."""
    if values is None:
        raise ValueError(f"QHA {name} grid is not available")
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"QHA {name} grid must be a non-empty one-dimensional array")
    if np.unique(array).size != array.size:
        raise ValueError(f"QHA {name} grid contains duplicate coordinates")
    return array


def _as_grid(
    values: np.ndarray,
    *,
    ntemperature: int,
    npressure: int,
    name: str,
) -> np.ndarray:
    """Validate one property as a native ``(nT, nP)`` grid."""
    if values.ndim == 1 and values.shape == (ntemperature,):
        return np.repeat(values[:, np.newaxis], npressure, axis=1)
    if values.ndim == 2 and values.shape == (ntemperature, npressure):
        return values
    raise ValueError(
        f"QHA property {name} has shape {values.shape}, expected "
        f"({ntemperature}, {npressure})"
    )


def _math_symbol(symbol: str) -> str:
    """Remove renderer delimiters from one established QHA symbol."""
    value = symbol.strip()
    if value.startswith("$") and value.endswith("$"):
        return value[1:-1]
    return value


def _property_category(unit_key: str) -> str:
    """Return one stable scientific category from the property unit family."""
    mapping = {
        "volume": "structural",
        "pressure": "elastic",
        "energy": "energy",
        "entropy": "entropy",
        "heat_capacity": "heat_capacity",
        "temperature_inverse": "thermal_expansion",
        "dimensionless": "dimensionless",
    }
    return mapping.get(unit_key, "thermodynamic")


def _property_description(attribute: str, description: str) -> str:
    """Add the explicit scaling convention for volumetric expansion."""
    if attribute == "thermal_expansion":
        return f"{description} Prepared plot values are scaled by 10^5."
    return description


def _property_unit(result: QHAResult, unit_key: str) -> str | None:
    """Resolve one frontend-neutral unit from native result metadata."""
    if unit_key == "volume":
        return _metadata_unit(result, "volume", fallback="A^3")
    if unit_key == "pressure":
        return _metadata_unit(result, "pressure", fallback="GPa")
    if unit_key == "energy":
        return _metadata_unit(result, "energy", fallback=None)
    if unit_key == "entropy":
        return _metadata_unit(
            result,
            "entropy",
            fallback=_metadata_unit(result, "heat_capacity", fallback=None),
        )
    if unit_key == "heat_capacity":
        return _metadata_unit(result, "heat_capacity", fallback=None)
    if unit_key == "temperature_inverse":
        temperature = _metadata_unit(result, "temperature", fallback="K") or "K"
        return f"10^-5 {temperature}^-1"
    if unit_key == "dimensionless":
        return None
    return _metadata_unit(result, unit_key, fallback=None)


def _metadata_unit(
    result: QHAResult,
    key: str,
    *,
    fallback: str | None,
) -> str | None:
    """Return a normalized native unit label."""
    units = result.metadata.get("units", {})
    value = units.get(key, fallback) if isinstance(units, dict) else fallback
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized or normalized.lower() in {"unknown", "none", "dimensionless"}:
        return None
    return normalized


__all__ = ["describe_qha_plots"]
