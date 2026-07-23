# -*- coding: utf-8 -*-

"""Extrapolation diagnostics for thermoelastic reports."""

from __future__ import annotations
from typing import Any
import numpy as np
from quantas.models import ReportTable
from quantas.modules.thermoelasticity.models import ThermoelasticResult
from .common import (
    _maximum_relative_uncertainty,
    _metadata_float,
    _metadata_range_text,
    _outside_distance,
)


def _profile_component_ranges_table(
    name: str,
    profile: Any,
    result: ThermoelasticResult,
) -> ReportTable:
    rows: list[list[Any]] = []
    for index, label in enumerate(result.independent_labels):
        component = profile.independent_stiffness[:, index]
        sigma_component = profile.sigma_independent_stiffness[:, index]
        rows.append(
            [
                label,
                float(np.nanmin(component)),
                float(np.nanmax(component)),
                float(np.nanmax(sigma_component)),
                _maximum_relative_uncertainty(component, sigma_component),
            ]
        )
    return ReportTable(
        f"Profile {name}: elastic-component ranges",
        [
            "Component",
            "Minimum (GPa)",
            "Maximum (GPa)",
            "Max σ (GPa)",
            "Max σ/|C| (%)",
        ],
        rows,
    )


def _extrapolation_policy_table(
    result: ThermoelasticResult,
    policy: str,
) -> ReportTable:
    descriptions = {
        "fail": "Abort before persisting any state that requires extrapolation.",
        "warn": "Evaluate extrapolated states and persist explicit warnings and masks.",
        "allow": (
            "Evaluate extrapolated states without workflow warnings; "
            "masks remain persisted."
        ),
    }
    return ReportTable(
        "Extrapolation debug context",
        ["Diagnostic", "Value"],
        [
            ["Policy", policy],
            ["Policy action", descriptions.get(policy, "Unknown policy")],
            [
                "QHA interpolation",
                result.metadata.get("qha_volume_interpolation", "not recorded"),
            ],
            [
                "Source T range (K)",
                _metadata_range_text(
                    result,
                    "source_grid_temperature_min_K",
                    "source_grid_temperature_max_K",
                ),
            ],
            [
                "Source P range (GPa)",
                _metadata_range_text(
                    result,
                    "source_grid_pressure_min_GPa",
                    "source_grid_pressure_max_GPa",
                ),
            ],
            [
                "Elastic V range (Å³)",
                _metadata_range_text(
                    result,
                    "elastic_volume_min_A3",
                    "elastic_volume_max_A3",
                ),
            ],
        ],
        metadata={
            "notes": [
                "QHA-coordinate extrapolation uses endpoint-slope "
                "piecewise-linear continuation of V(P,T).",
                "Elastic extrapolation evaluates the fitted finite-strain "
                "Cij(V) model outside the sampled SOEC volume interval.",
            ],
            "debug_only": True,
        },
    )


def _grid_extrapolation_tables(
    result: ThermoelasticResult,
) -> tuple[ReportTable, ...]:
    """Return compact domain and uncertainty diagnostics for grid extrapolation."""
    qha_mask = np.asarray(result.qha_extrapolation_mask, dtype=np.bool_)
    elastic_mask = np.asarray(result.extrapolation_mask, dtype=np.bool_)
    mask = qha_mask | elastic_mask
    indices = np.argwhere(mask)
    if indices.size == 0:
        return (
            ReportTable(
                "Grid extrapolation diagnostics",
                ["Status", "Value"],
                [
                    ["Extrapolated states", 0],
                    ["Assessment", "All grid states are interpolated"],
                ],
                metadata={"debug_only": True},
            ),
        )
    rows: list[list[Any]] = []
    for it, ip in indices[:500]:
        rows.append(
            _extrapolation_row(
                result,
                pressure=float(result.pressure[ip]),
                temperature=float(result.temperature[it]),
                volume=float(result.equilibrium_volume[it, ip]),
                qha=bool(qha_mask[it, ip]),
                elastic=bool(elastic_mask[it, ip]),
                sigma_values=(
                    None
                    if result.sigma_independent_stiffness is None
                    else result.sigma_independent_stiffness[it, ip]
                ),
                component_values=(
                    None
                    if result.independent_stiffness is None
                    else result.independent_stiffness[it, ip]
                ),
                identifier=f"T{it + 1:04d}_P{ip + 1:04d}",
                depth=None,
            )
        )
    note = (
        f"Showing {min(indices.shape[0], 500)} of {indices.shape[0]} "
        "extrapolated grid states; all masks and uncertainties are persisted "
        "in HDF5."
    )
    return _split_extrapolation_tables(
        "Grid extrapolation",
        rows,
        notes=[note],
    )


def _profile_extrapolation_tables(
    name: str,
    profile: Any,
    result: ThermoelasticResult,
) -> tuple[ReportTable, ...]:
    """Return compact domain and uncertainty diagnostics for one profile."""
    mask = np.asarray(profile.qha_extrapolation_mask, dtype=np.bool_) | np.asarray(
        profile.elastic_extrapolation_mask, dtype=np.bool_
    )
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return (
            ReportTable(
                f"Profile {name}: extrapolation diagnostics",
                ["Status", "Value"],
                [
                    ["Extrapolated states", 0],
                    ["Assessment", "All profile states are interpolated"],
                ],
                metadata={"debug_only": True},
            ),
        )
    rows: list[list[Any]] = []
    for index in indices[:500]:
        rows.append(
            _extrapolation_row(
                result,
                pressure=float(profile.pressure[index]),
                temperature=float(profile.temperature[index]),
                volume=float(profile.volume[index]),
                qha=bool(profile.qha_extrapolation_mask[index]),
                elastic=bool(profile.elastic_extrapolation_mask[index]),
                sigma_values=profile.sigma_independent_stiffness[index],
                component_values=profile.independent_stiffness[index],
                identifier=f"{name}:{index + 1:04d}",
                depth=float(profile.depth[index]),
            )
        )
    note = (
        f"Showing {min(indices.size, 500)} of {indices.size} extrapolated "
        "profile states."
    )
    return _split_extrapolation_tables(
        f"Profile {name}: extrapolation",
        rows,
        notes=[note],
    )


def _split_extrapolation_tables(
    title: str,
    rows: list[list[Any]],
    *,
    notes: list[str],
) -> tuple[ReportTable, ReportTable]:
    """Split wide extrapolation diagnostics into two readable tables."""
    domains = [[row[index] for index in (0, 1, 2, 3, 4, 5, 8)] for row in rows]
    distances = [[row[index] for index in (0, 6, 7, 9, 10, 11, 12)] for row in rows]
    metadata = {"debug_only": True, "notes": notes}
    return (
        ReportTable(
            f"{title} states",
            [
                "State",
                "Depth (km)",
                "P (GPa)",
                "T (K)",
                "V (Å³)",
                "QHA extrap.",
                "Elastic extrap.",
            ],
            domains,
            metadata=metadata,
        ),
        ReportTable(
            f"{title} distances and uncertainty",
            [
                "State",
                "ΔP outside (GPa)",
                "ΔT outside (K)",
                "ΔV outside (Å³)",
                "ΔV/span (%)",
                "Max σ(C) (GPa)",
                "Max σ/|C| (%)",
            ],
            distances,
            metadata=metadata,
        ),
    )


def _extrapolation_row(
    result: ThermoelasticResult,
    *,
    pressure: float,
    temperature: float,
    volume: float,
    qha: bool,
    elastic: bool,
    sigma_values: Any,
    component_values: Any,
    identifier: str,
    depth: float | None,
) -> list[Any]:
    p_distance = _outside_distance(
        pressure,
        result.metadata.get("source_grid_pressure_min_GPa"),
        result.metadata.get("source_grid_pressure_max_GPa"),
    )
    t_distance = _outside_distance(
        temperature,
        result.metadata.get("source_grid_temperature_min_K"),
        result.metadata.get("source_grid_temperature_max_K"),
    )
    v_lower = _metadata_float(result, "elastic_volume_min_A3")
    v_upper = _metadata_float(result, "elastic_volume_max_A3")
    v_distance = _outside_distance(volume, v_lower, v_upper)
    v_span = None if v_lower is None or v_upper is None else v_upper - v_lower
    relative_v = (
        None
        if v_distance is None or v_span is None or v_span <= 0.0
        else 100.0 * abs(v_distance) / v_span
    )
    sigma_array = (
        None if sigma_values is None else np.asarray(sigma_values, dtype=np.float64)
    )
    component_array = (
        None
        if component_values is None
        else np.asarray(component_values, dtype=np.float64)
    )
    return [
        identifier,
        depth,
        pressure,
        temperature,
        volume,
        qha,
        p_distance,
        t_distance,
        elastic,
        v_distance,
        relative_v,
        None if sigma_array is None else float(np.nanmax(sigma_array)),
        _maximum_relative_uncertainty(component_array, sigma_array),
    ]


__all__ = [
    "_extrapolation_policy_table",
    "_grid_extrapolation_tables",
    "_profile_extrapolation_tables",
]
