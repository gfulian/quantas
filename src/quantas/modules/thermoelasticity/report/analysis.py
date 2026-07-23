# -*- coding: utf-8 -*-

"""Grid- and profile-analysis thermoelastic reports."""

from __future__ import annotations

from typing import Any

import numpy as np
from quantas.models import ReportTable
from quantas.modules.thermoelasticity.models import (
    ThermoelasticReportLevel,
    ThermoelasticResult,
)
from .common import (
    _maximum_relative_uncertainty,
    _minimum_stability_eigenvalue,
    _range_text,
)
from .extrapolation import (
    _extrapolation_policy_table,
    _profile_component_ranges_table,
    _grid_extrapolation_tables,
    _profile_extrapolation_tables,
)


def build_thermoelastic_analysis_report(
    result: ThermoelasticResult,
    *,
    level: ThermoelasticReportLevel = "standard",
    extrapolation_policy: str | None = None,
) -> tuple[ReportTable, ...]:
    """Build post-fit pressure-temperature analysis reports.

    Unlike :func:`build_thermoelastic_report`, this function does not repeat
    static-EOS and component-fit diagnostics.  It reports only the requested
    reconstructed states, uncertainty ranges, geological profiles, and—at
    debug level—the quantitative distance of extrapolated states from the
    archived QHA and elastic domains.

    Parameters
    ----------
    result : ThermoelasticResult
        Reconstructed pressure-temperature result.
    level : {"standard", "extended", "debug"}, optional
        Requested amount of analysis information.
    extrapolation_policy : str or None, optional
        Effective policy used for the reconstruction.  When omitted, the
        value stored in result metadata is used.

    Returns
    -------
    tuple of ReportTable
        Ordered frontend-neutral report tables.
    """
    if level not in {"standard", "extended", "debug"}:
        raise ValueError("invalid thermoelastic report level")
    policy = (
        str(result.metadata.get("extrapolation_policy", "warn"))
        if extrapolation_policy is None
        else str(extrapolation_policy)
    )
    tables: list[ReportTable] = [_analysis_summary_table(result, policy)]
    profile_summary = _analysis_profile_summary_table(result)
    if profile_summary is not None:
        tables.append(profile_summary)
    if level in {"extended", "debug"}:
        grid_ranges = _analysis_component_ranges_table(result)
        if grid_ranges is not None:
            tables.append(grid_ranges)
        for name, profile in result.profiles.items():
            tables.append(_profile_component_ranges_table(name, profile, result))
    if level == "debug":
        tables.append(_extrapolation_policy_table(result, policy))
        tables.extend(_grid_extrapolation_tables(result))
        for name, profile in result.profiles.items():
            tables.extend(_profile_extrapolation_tables(name, profile, result))
    return tuple(tables)


def _analysis_summary_table(
    result: ThermoelasticResult,
    policy: str,
) -> ReportTable:
    grid_tensors = (
        0
        if result.stiffness_isothermal is None
        else int(result.temperature.size * result.pressure.size)
    )
    qha_count = int(
        np.count_nonzero(np.asarray(result.qha_extrapolation_mask, dtype=np.bool_))
    )
    elastic_count = int(np.count_nonzero(result.extrapolation_mask))
    rows: list[list[Any]] = [
        ["Analysis kind", result.metadata.get("analysis_kind", "P-T reconstruction")],
        ["Grid reconstructed", result.stiffness_isothermal is not None],
        ["Grid tensors", grid_tensors],
        ["Temperature range (K)", _range_text(result.temperature)],
        ["Pressure range (GPa)", _range_text(result.pressure)],
        ["Volume range (Å³)", _range_text(result.equilibrium_volume)],
        ["Density range (kg m⁻³)", _range_text(result.density)],
        ["Depth profiles", len(result.profiles)],
        ["QHA-coordinate extrapolated states", qha_count],
        ["Elastic-volume extrapolated states", elastic_count],
        ["Extrapolation policy", policy],
        [
            "Mechanically stable states",
            None
            if result.stability is None
            else int(np.count_nonzero(result.stability.stable_mask)),
        ],
        [
            "Mechanically unstable states",
            None
            if result.stability is None
            else int(np.count_nonzero(result.stability.unstable_mask)),
        ],
        [
            "Stability indeterminate states",
            None
            if result.stability is None
            else int(np.count_nonzero(result.stability.indeterminate_mask)),
        ],
        [
            "Minimum stiffness eigenvalue (GPa)",
            _minimum_stability_eigenvalue(result.stability),
        ],
        [
            "Isothermal tensor condition",
            "available" if result.stiffness_isothermal is not None else "absent",
        ],
        [
            "Adiabatic tensor condition",
            "absent" if result.stiffness_adiabatic is None else "available",
        ],
        [
            "Adiabatic valid grid states",
            None
            if result.adiabatic_valid_mask is None
            else int(np.count_nonzero(result.adiabatic_valid_mask)),
        ],
        [
            "Adiabatic invalid grid states",
            None
            if result.adiabatic_valid_mask is None
            else int(
                result.adiabatic_valid_mask.size
                - np.count_nonzero(result.adiabatic_valid_mask)
            ),
        ],
        [
            "Maximum |C^S-C^T| (GPa)",
            None
            if result.adiabatic_correction is None
            else float(np.nanmax(np.abs(result.adiabatic_correction))),
        ],
        ["Approximation", "quasi-static"],
    ]
    return ReportTable(
        "Pressure-temperature analysis summary",
        ["Property", "Value"],
        rows,
        metadata={
            "notes": [
                "Fit diagnostics are intentionally omitted here; they belong "
                "to the thermoelastic run report.",
                "Extended output summarizes ranges rather than printing "
                "every grid state.",
            ]
        },
    )


def _analysis_profile_summary_table(
    result: ThermoelasticResult,
) -> ReportTable | None:
    if not result.profiles:
        return None
    rows: list[list[Any]] = []
    for name, profile in result.profiles.items():
        rows.append(
            [
                name,
                profile.depth.size,
                _range_text(profile.depth),
                _range_text(profile.pressure),
                _range_text(profile.temperature),
                _range_text(profile.volume),
                int(
                    np.count_nonzero(
                        profile.qha_extrapolation_mask
                        | profile.elastic_extrapolation_mask
                    )
                ),
                None
                if profile.stability is None
                else int(np.count_nonzero(profile.stability.unstable_mask)),
                profile.metadata.get("kind", "user"),
            ]
        )
    return ReportTable(
        "Geological profile summary",
        [
            "Profile",
            "Points",
            "Depth (km)",
            "Pressure (GPa)",
            "Temperature (K)",
            "Volume (Å³)",
            "Extrapolated",
            "Unstable",
            "Kind",
        ],
        rows,
    )


def _analysis_component_ranges_table(
    result: ThermoelasticResult,
) -> ReportTable | None:
    values = result.independent_stiffness
    if values is None:
        return None
    sigma = result.sigma_independent_stiffness
    rows: list[list[Any]] = []
    for index, label in enumerate(result.independent_labels):
        component = values[..., index]
        sigma_component = None if sigma is None else sigma[..., index]
        rows.append(
            [
                label,
                float(np.nanmin(component)),
                float(np.nanmax(component)),
                None if sigma_component is None else float(np.nanmax(sigma_component)),
                _maximum_relative_uncertainty(component, sigma_component),
            ]
        )
    return ReportTable(
        "P-T elastic-component ranges",
        [
            "Component",
            "Minimum (GPa)",
            "Maximum (GPa)",
            "Max σ (GPa)",
            "Max σ/|C| (%)",
        ],
        rows,
        metadata={
            "column_formats": [None, ".4f", ".4f", ".4f", ".4f"],
            "notes": [
                "Ranges summarize the reconstructed grid; individual tensors "
                "remain in HDF5 and are selected explicitly by export."
            ],
        },
    )


__all__ = ["build_thermoelastic_analysis_report"]
