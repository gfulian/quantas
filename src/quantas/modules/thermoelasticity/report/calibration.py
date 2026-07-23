# -*- coding: utf-8 -*-

"""Calibration-stage thermoelastic reports."""

from __future__ import annotations
from typing import Any, Mapping
import numpy as np
from quantas.models import ReportTable
from quantas.modules.thermoelasticity.models import (
    ElasticComponentFit,
    ThermoelasticReportLevel,
    ThermoelasticResult,
)
from .common import _minimum_stability_eigenvalue


def build_thermoelastic_report(
    result: ThermoelasticResult,
    *,
    level: ThermoelasticReportLevel = "standard",
) -> tuple[ReportTable, ...]:
    """Build standard, extended, or debug thermoelastic report tables.

    Parameters
    ----------
    result : ThermoelasticResult
        Completed or diagnostic-only QSA result.
    level : {"standard", "extended", "debug"}, optional
        Requested amount of rendered information. HDF5 persistence is
        unaffected and always contains complete diagnostics.

    Returns
    -------
    tuple of ReportTable
        Ordered neutral tables suitable for Rich, plain text, or a GUI.
    """
    if level not in {"standard", "extended", "debug"}:
        raise ValueError("invalid thermoelastic report level")
    tables: list[ReportTable] = [
        _reference_eos_table(result),
        _component_summary_table(result),
        _reconstruction_summary_table(result),
    ]
    if level in {"extended", "debug"}:
        for component in result.component_fits.values():
            tables.append(_component_points_table(component))
            covariance = _component_covariance_table(component)
            if covariance is not None:
                tables.append(covariance)
    if level == "debug":
        for component in result.component_fits.values():
            tables.append(_component_debug_table(component))
            trace = _component_trace_table(component)
            if trace is not None:
                tables.append(trace)
    return tuple(tables)


def _reference_eos_table(result: ThermoelasticResult) -> ReportTable:
    reference = result.reference_eos
    fit = reference.fit
    parameter_errors: dict[str, float | None] = {
        "V0": None,
        "K0": None,
        "Kprime": None,
    }
    if reference.covariance is not None:
        diagonal = np.sqrt(np.clip(np.diag(reference.covariance), 0.0, None))
        parameter_errors = dict(zip(parameter_errors, diagonal, strict=True))
    rows: list[list[Any]] = [
        ["EOS", reference.eos, None],
        ["V0", reference.reference_volume, parameter_errors["V0"]],
        ["K0", reference.bulk_modulus, parameter_errors["K0"]],
        ["Kprime", reference.bulk_modulus_derivative, parameter_errors["Kprime"]],
        ["Ksecond", reference.bulk_modulus_second_derivative, None],
        ["Fit success", fit.success, None],
        ["Fit status", fit.status.value, None],
        ["Fit quality", fit.quality.value, None],
        ["R squared", fit.r_squared, None],
        ["RMSE", fit.rmse, None],
        ["Maximum residual", fit.max_abs_error, None],
    ]
    return ReportTable(
        "Static reference EOS",
        ["Quantity", "Value", "E.S.D."],
        rows,
        metadata={
            "units": {"V0": "angstrom^3", "K0": "GPa", "Ksecond": "GPa^-1"},
            "notes": ["V0, K0 and Kprime are fixed for all elastic-component fits."],
            "column_formats": [None, ".4f", ".4f"],
        },
    )


def _component_summary_table(result: ThermoelasticResult) -> ReportTable:
    rows: list[list[Any]] = []
    for label, record in result.component_fits.items():
        fit = record.fit
        if record.zero_by_tolerance:
            rows.append(
                [
                    label,
                    "zero",
                    True,
                    "not_applicable",
                    "component_retained_as_exact_zero",
                    0.0,
                    0.0,
                    None,
                    None,
                    None,
                    0.0,
                    0.0,
                    0.0,
                    None,
                    None,
                    None,
                ]
            )
            continue
        parameters = None if fit is None else fit.parameters
        errors = None if fit is None else fit.errors
        metadata = {} if fit is None else fit.metadata
        quality = record.quality
        rows.append(
            [
                label,
                "missing" if fit is None else fit.status.value,
                False if fit is None else fit.success,
                None if quality is None else quality.level,
                None if quality is None else ", ".join(quality.issues),
                None if parameters is None else float(parameters[0]),
                None if errors is None else float(errors[0]),
                None if parameters is None else float(parameters[1]),
                None if errors is None else float(errors[1]),
                None if fit is None else fit.r_squared,
                None if fit is None else metadata.get("unweighted_chi_square_GPa2"),
                None if fit is None else fit.max_abs_error,
                None
                if fit is None
                else 100.0 * float(metadata.get("maximum_relative_error", np.nan)),
                None if quality is None else quality.eulerian_strain_span,
                None if quality is None else quality.design_condition_number,
                None
                if quality is None
                else quality.maximum_leave_one_out_parameter_change,
            ]
        )
    return ReportTable(
        "Quasi-static elastic-component fits",
        [
            "Component",
            "Status",
            "Success",
            "Support",
            "Support issues",
            "C0 (GPa)",
            "σ(C0)",
            "Cprime",
            "σ(Cprime)",
            "R²",
            "χ² / RSS (GPa²)",
            "Max |Δ| (GPa)",
            "Max |Δrel| (%)",
            "ΔfE",
            "Design κ₂",
            "Max LOO Δparameter",
        ],
        rows,
        metadata={
            "column_formats": [
                None,
                None,
                None,
                None,
                None,
                ".4f",
                ".4f",
                ".4f",
                ".4f",
                ".6f",
                ".4e",
                ".4f",
                ".4f",
                ".6f",
                ".4e",
                ".4f",
            ],
            "notes": [
                "χ² is the unweighted residual sum of squares because "
                "CRYSTAL SOEC uncertainties are unavailable.",
                "All complete fit diagnostics and covariances are persisted "
                "in HDF5 independently of report level.",
                "Support diagnostics classify data coverage and sensitivity; "
                "they do not alter the fitted parameters.",
            ],
        },
    )


def _reconstruction_summary_table(result: ThermoelasticResult) -> ReportTable:
    qha_mask = np.asarray(result.qha_extrapolation_mask, dtype=np.bool_)
    elastic_mask = np.asarray(result.extrapolation_mask, dtype=np.bool_)
    rows: list[list[Any]] = [
        ["Completed", result.completed],
        ["Temperature points", result.temperature.size],
        ["Pressure points", result.pressure.size],
        [
            "Grid tensors",
            0
            if result.stiffness_isothermal is None
            else result.temperature.size * result.pressure.size,
        ],
        [
            "Grid reconstruction",
            "deferred" if result.stiffness_isothermal is None else "completed",
        ],
        ["Independent components", len(result.independent_labels)],
        [
            "QHA-coordinate extrapolated points",
            int(np.count_nonzero(qha_mask)),
        ],
        [
            "Elastic extrapolated points",
            int(np.count_nonzero(elastic_mask)),
        ],
        ["Depth profiles", len(result.profiles)],
        [
            "Mechanically stable grid states",
            None
            if result.stability is None
            else int(np.count_nonzero(result.stability.stable_mask)),
        ],
        [
            "Mechanically unstable grid states",
            None
            if result.stability is None
            else int(np.count_nonzero(result.stability.unstable_mask)),
        ],
        [
            "Mechanical stability indeterminate",
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
    for name, profile in result.profiles.items():
        rows.extend(
            [
                [f"Profile {name}: points", profile.depth.size],
                [
                    f"Profile {name}: extrapolated",
                    int(
                        np.count_nonzero(
                            profile.qha_extrapolation_mask
                            | profile.elastic_extrapolation_mask
                        )
                    ),
                ],
                [
                    f"Profile {name}: unstable",
                    None
                    if profile.stability is None
                    else int(np.count_nonzero(profile.stability.unstable_mask)),
                ],
                [
                    f"Profile {name}: adiabatic valid",
                    None
                    if profile.adiabatic_valid_mask is None
                    else int(np.count_nonzero(profile.adiabatic_valid_mask)),
                ],
            ]
        )
    return ReportTable(
        "Thermoelastic reconstruction",
        ["Property", "Value"],
        rows,
    )


def _component_points_table(component: ElasticComponentFit) -> ReportTable:
    leverage = (
        np.zeros(component.volumes.shape, dtype=np.float64)
        if component.quality is None
        else component.quality.leverage
    )
    rows = [
        [
            float(volume),
            float(pressure),
            float(observed),
            float(fitted) if np.isfinite(fitted) else None,
            float(residual) if np.isfinite(residual) else None,
            100.0 * float(relative) if np.isfinite(relative) else None,
            float(spread),
            float(point_leverage),
        ]
        for (
            volume,
            pressure,
            observed,
            fitted,
            residual,
            relative,
            spread,
            point_leverage,
        ) in zip(
            component.volumes,
            component.pressures,
            component.observed,
            component.fitted,
            component.residuals,
            component.relative_residuals,
            component.symmetry_spread,
            leverage,
            strict=True,
        )
    ]
    return ReportTable(
        f"{component.label} observed and fitted values",
        [
            "V (Å³)",
            "P (GPa)",
            "Observed (GPa)",
            "Fitted (GPa)",
            "Residual (GPa)",
            "Relative residual (%)",
            "Symmetry spread (GPa)",
            "Leverage",
        ],
        rows,
        metadata={
            "column_formats": [None, ".4f", ".4f", ".4f", ".4f", ".4f", ".4f", ".4f"]
        },
    )


def _component_covariance_table(
    component: ElasticComponentFit,
) -> ReportTable | None:
    fit = component.fit
    if fit is None or fit.covariance is None:
        return None
    correlation = None if fit.diagnostics is None else fit.diagnostics.correlation
    names = fit.parameter_names or ("C0", "Cprime")
    rows: list[list[Any]] = []
    for i, first in enumerate(names):
        for j, second in enumerate(names):
            rows.append(
                [
                    first,
                    second,
                    float(fit.covariance[i, j]),
                    None if correlation is None else float(correlation[i, j]),
                ]
            )
    return ReportTable(
        f"{component.label} parameter covariance",
        ["Parameter i", "Parameter j", "Covariance", "Correlation"],
        rows,
    )


def _component_debug_table(component: ElasticComponentFit) -> ReportTable:
    fit = component.fit
    rows: list[list[Any]]
    if fit is None:
        rows = [["Status", "exact zero"], ["Message", "No optimizer required"]]
    else:
        diagnostics = fit.diagnostics
        metadata: Mapping[str, Any] = fit.metadata
        diagnostic_metadata: Mapping[str, Any] = (
            {} if diagnostics is None else diagnostics.metadata
        )
        rows = [
            ["Status", fit.status.value],
            ["Quality", fit.quality.value],
            ["Message", fit.message],
            ["Stop reason", "" if diagnostics is None else diagnostics.stop_reason],
            ["Objective", "" if diagnostics is None else diagnostics.objective],
            ["Iterations", None if diagnostics is None else diagnostics.n_iterations],
            ["Evaluations", None if diagnostics is None else diagnostics.n_evaluations],
            [
                "Jacobian rank",
                None if diagnostics is None else diagnostics.jacobian_rank,
            ],
            [
                "Condition number",
                None if diagnostics is None else diagnostics.condition_number,
            ],
            ["Warnings", "; ".join(fit.warnings)],
            ["Termination category", diagnostic_metadata.get("termination_category")],
            ["Backend", diagnostic_metadata.get("backend")],
            ["Trace truncated", diagnostic_metadata.get("evaluation_trace_truncated")],
            ["Maximum symmetry spread", metadata.get("maximum_symmetry_spread_GPa")],
            [
                "Scientific support",
                None if component.quality is None else component.quality.level,
            ],
            [
                "Support issues",
                None
                if component.quality is None
                else "; ".join(component.quality.issues),
            ],
            [
                "Eulerian strain span",
                None
                if component.quality is None
                else component.quality.eulerian_strain_span,
            ],
            [
                "Reference volume bracketed",
                None
                if component.quality is None
                else component.quality.reference_volume_bracketed,
            ],
            [
                "Design condition number",
                None
                if component.quality is None
                else component.quality.design_condition_number,
            ],
            [
                "Maximum leave-one-out parameter change",
                None
                if component.quality is None
                else component.quality.maximum_leave_one_out_parameter_change,
            ],
            [
                "Finite-strain-order parameter change",
                None
                if component.quality is None
                else component.quality.maximum_order_parameter_change,
            ],
        ]
    return ReportTable(
        f"{component.label} solver diagnostics",
        ["Diagnostic", "Value"],
        rows,
    )


def _component_trace_table(component: ElasticComponentFit) -> ReportTable | None:
    fit = component.fit
    if fit is None or fit.diagnostics is None:
        return None
    trace = fit.diagnostics.metadata.get("evaluation_trace")
    if not isinstance(trace, list) or not trace:
        return None
    rows: list[list[Any]] = []
    for index, record in enumerate(trace):
        if not isinstance(record, Mapping):
            continue
        rows.append(
            [
                index + 1,
                record.get("parameters"),
                record.get("finite"),
                record.get("minimum"),
                record.get("maximum"),
                record.get("error"),
            ]
        )
    if not rows:
        return None
    return ReportTable(
        f"{component.label} model-evaluation trace",
        ["Evaluation", "Parameters", "Finite", "Minimum", "Maximum", "Error"],
        rows,
        metadata={"debug_only": True},
    )


__all__ = ["build_thermoelastic_report"]
