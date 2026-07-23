# -*- coding: utf-8 -*-

"""Shared helpers for frontend-neutral thermoelastic plot builders."""

from __future__ import annotations

from dataclasses import dataclass
from statistics import NormalDist
from typing import Iterable

import numpy as np
from numpy.typing import NDArray

from quantas.models import PlotMask, PlotSpan, ResultData
from quantas.modules.thermoelasticity.models import (
    ThermoelasticProfileResult,
    ThermoelasticResult,
)
from quantas.modules.thermoelasticity.plot.components import component_indices
from quantas.modules.thermoelasticity.plot.options import (
    ThermoelasticPlotStyleOptions,
)


FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class ResolvedThermoelasticStyle:
    """Concrete common style values derived from a thermoelastic preset.

    Parameters
    ----------
    preset : str
        Selected high-level style profile.
    show_title, grid : bool
        Concrete title and grid choices.
    line_width, marker_size : float
        Concrete curve dimensions.
    marker_edge_color : str or None
        Concrete marker border color.
    marker_edge_width : float
        Concrete marker border width.
    errorbar_width, errorbar_capsize : float
        Concrete error-bar dimensions.
    """

    preset: str
    show_title: bool
    grid: bool
    line_width: float
    marker_size: float
    marker_edge_color: str | None
    marker_edge_width: float
    errorbar_width: float
    errorbar_capsize: float


def extract_thermoelastic_result(
    result: ResultData | ThermoelasticResult,
) -> ThermoelasticResult:
    """Return the thermoelastic payload from a result envelope or payload.

    Parameters
    ----------
    result : ResultData or ThermoelasticResult
        Native result envelope or direct passive payload.

    Returns
    -------
    ThermoelasticResult
        Validated thermoelastic payload.

    Raises
    ------
    ValueError
        If a generic result envelope does not contain the expected payload.
    """
    if isinstance(result, ThermoelasticResult):
        return result
    payload = result.results.get("thermoelasticity")
    if not isinstance(payload, ThermoelasticResult):
        raise ValueError("result does not contain a thermoelasticity payload")
    return payload


def resolve_style(options: ThermoelasticPlotStyleOptions) -> ResolvedThermoelasticStyle:
    """Resolve nullable common options against one style preset."""
    if options.preset == "analysis":
        default_show_title = True
        default_grid = True
        default_line_width = 1.4
        default_marker_size = 5.0
    elif options.preset == "publication":
        default_show_title = False
        default_grid = False
        default_line_width = 1.2
        default_marker_size = 3.5
    else:
        default_show_title = True
        default_grid = True
        default_line_width = 1.4
        default_marker_size = 5.0
    return ResolvedThermoelasticStyle(
        preset=options.preset,
        show_title=(
            default_show_title if options.show_title is None else options.show_title
        ),
        grid=default_grid if options.grid is None else options.grid,
        line_width=(
            default_line_width
            if options.line_width is None
            else float(options.line_width)
        ),
        marker_size=(
            default_marker_size
            if options.marker_size is None
            else float(options.marker_size)
        ),
        marker_edge_color=options.marker_edge_color,
        marker_edge_width=float(options.marker_edge_width),
        errorbar_width=float(options.errorbar_width),
        errorbar_capsize=float(options.errorbar_capsize),
    )


def component_grid(
    result: ThermoelasticResult,
    label: str,
    *,
    tensor_condition: str = "isothermal",
) -> tuple[FloatArray, FloatArray]:
    """Return one component and one-sigma uncertainty on the P-T grid.

    Parameters
    ----------
    result : ThermoelasticResult
        Result containing reconstructed tensors.
    label : str
        Canonical component label.

    Returns
    -------
    tuple of ndarray
        Values and uncertainties with shape ``(nT, nP)``.

    Raises
    ------
    ValueError
        If reconstructed tensors are absent.
    """
    if tensor_condition == "isothermal":
        stiffness = result.stiffness_isothermal
        sigma_stiffness = result.sigma_stiffness_isothermal
    elif tensor_condition == "adiabatic":
        stiffness = result.stiffness_adiabatic
        sigma_stiffness = result.sigma_stiffness_adiabatic
    else:
        raise ValueError("tensor_condition must be 'isothermal' or 'adiabatic'")
    if stiffness is None:
        raise ValueError(
            f"thermoelastic result does not contain {tensor_condition} tensors"
        )
    first, second = component_indices(label)
    values = np.asarray(stiffness[..., first, second], dtype=np.float64)
    if sigma_stiffness is None:
        uncertainty = np.full(values.shape, np.nan, dtype=np.float64)
    else:
        uncertainty = np.asarray(
            sigma_stiffness[..., first, second],
            dtype=np.float64,
        )
    return values.copy(), uncertainty.copy()


def profile_component(
    profile: ThermoelasticProfileResult,
    label: str,
    *,
    tensor_condition: str = "isothermal",
) -> tuple[FloatArray, FloatArray]:
    """Return one component and one-sigma uncertainty along a depth profile."""
    stiffness: np.ndarray | None
    sigma_stiffness: np.ndarray | None
    if tensor_condition == "isothermal":
        stiffness = profile.stiffness_isothermal
        sigma_stiffness = profile.sigma_stiffness_isothermal
    elif tensor_condition == "adiabatic":
        stiffness = profile.stiffness_adiabatic
        sigma_stiffness = profile.sigma_stiffness_adiabatic
    else:
        raise ValueError("tensor_condition must be 'isothermal' or 'adiabatic'")
    if stiffness is None:
        raise ValueError(f"profile does not contain {tensor_condition} tensors")
    first, second = component_indices(label)
    values = np.asarray(stiffness[..., first, second], dtype=np.float64)
    if sigma_stiffness is None:
        uncertainty = np.full(values.shape, np.nan, dtype=np.float64)
    else:
        uncertainty = np.asarray(sigma_stiffness[..., first, second], dtype=np.float64)
    return values.copy(), uncertainty.copy()


def confidence_multiplier(probability: float) -> float:
    """Return the two-sided standard-normal multiplier for a probability."""
    return float(NormalDist().inv_cdf(0.5 + 0.5 * float(probability)))


def extrapolation_masks(
    result: ThermoelasticResult,
) -> list[PlotMask]:
    """Build distinct QHA and elastic extrapolation masks for a P-T map."""
    masks: list[PlotMask] = []
    qha = np.asarray(result.qha_extrapolation_mask, dtype=np.bool_)
    elastic = np.asarray(result.extrapolation_mask, dtype=np.bool_)
    if np.any(qha):
        masks.append(
            PlotMask(
                key="qha_extrapolation",
                label="QHA coordinate extrapolation",
                x=result.pressure,
                y=result.temperature,
                mask=qha,
                hatch="///",
                color="none",
            )
        )
    if np.any(elastic):
        masks.append(
            PlotMask(
                key="elastic_extrapolation",
                label="Elastic-volume extrapolation",
                x=result.pressure,
                y=result.temperature,
                mask=elastic,
                hatch="\\\\\\",
                color="none",
            )
        )
    return masks


def profile_extrapolation_spans(
    profile: ThermoelasticProfileResult,
) -> list[PlotSpan]:
    """Return contiguous depth intervals requiring extrapolation."""
    spans: list[PlotSpan] = []
    for key, label, mask, hatch in (
        (
            "qha_extrapolation",
            "QHA coordinate extrapolation",
            profile.qha_extrapolation_mask,
            "///",
        ),
        (
            "elastic_extrapolation",
            "Elastic-volume extrapolation",
            profile.elastic_extrapolation_mask,
            "\\\\\\",
        ),
    ):
        for index, (start, end) in enumerate(
            contiguous_intervals(profile.depth, np.asarray(mask, dtype=np.bool_))
        ):
            spans.append(
                PlotSpan(
                    key=f"{key}_{index}",
                    label=label if index == 0 else "_nolegend_",
                    axis="y",
                    start=start,
                    end=end,
                    color="0.82",
                    alpha=0.18,
                    hatch=hatch,
                    metadata={"extrapolation_kind": key},
                )
            )
    return spans


def contiguous_intervals(
    coordinates: Iterable[float],
    mask: Iterable[bool],
) -> tuple[tuple[float, float], ...]:
    """Return cell-edge intervals covering contiguous true mask values.

    Parameters
    ----------
    coordinates : iterable of float
        Strictly ordered sample centers.
    mask : iterable of bool
        Boolean state for each center.

    Returns
    -------
    tuple of tuple
        Closed intervals in the coordinate units.
    """
    values = np.asarray(tuple(coordinates), dtype=np.float64)
    flags = np.asarray(tuple(mask), dtype=np.bool_)
    if values.ndim != 1 or flags.shape != values.shape:
        raise ValueError("coordinates and mask must be aligned one-dimensional data")
    if values.size == 0 or not np.any(flags):
        return ()
    edges = _coordinate_edges(values)
    intervals: list[tuple[float, float]] = []
    start: int | None = None
    for index, active in enumerate(flags):
        if active and start is None:
            start = index
        if start is not None and (not active or index == flags.size - 1):
            stop = index if active and index == flags.size - 1 else index - 1
            intervals.append((float(edges[start]), float(edges[stop + 1])))
            start = None
    return tuple(intervals)


def _coordinate_edges(coordinates: FloatArray) -> FloatArray:
    """Return cell edges for ordered one-dimensional centers."""
    if coordinates.size == 1:
        return np.asarray(
            [coordinates[0] - 0.5, coordinates[0] + 0.5], dtype=np.float64
        )
    midpoint = 0.5 * (coordinates[:-1] + coordinates[1:])
    first = coordinates[0] - (midpoint[0] - coordinates[0])
    last = coordinates[-1] + (coordinates[-1] - midpoint[-1])
    return np.concatenate(([first], midpoint, [last])).astype(np.float64)


__all__ = [
    "ResolvedThermoelasticStyle",
    "component_grid",
    "confidence_multiplier",
    "contiguous_intervals",
    "extrapolation_masks",
    "extract_thermoelastic_result",
    "profile_component",
    "profile_extrapolation_spans",
    "resolve_style",
]
