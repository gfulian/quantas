# -*- coding: utf-8 -*-

"""Frontend-neutral tables for EOS batch input and fitting results."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

import numpy as np

from quantas.models import ReportTable

from .batch import EOSBatchJob, EOSBatchJobResult, EOSBatchPlan, EOSBatchResult
from .models import EOSDataset, EOSFitResult
from .presentation import (
    domain_label,
    format_unit,
    model_label,
    parameter_label,
    solver_label,
    target_label,
)


class EOSReportDetail(str, Enum):
    """Amount of EOS result detail requested from a renderer."""

    SHORT = "short"
    EXTENDED = "extended"


@dataclass(frozen=True, slots=True)
class EOSReportOptions:
    """Configure frontend-neutral EOS batch report generation."""

    detail: EOSReportDetail = EOSReportDetail.SHORT
    show_uncertainties: bool = False
    max_data_rows: int | None = None
    debug: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "detail", EOSReportDetail(self.detail))
        if self.max_data_rows is not None and int(self.max_data_rows) <= 0:
            raise ValueError("max_data_rows must be positive or None")
        if self.max_data_rows is not None:
            object.__setattr__(self, "max_data_rows", int(self.max_data_rows))
        object.__setattr__(self, "debug", bool(self.debug))


def build_eos_batch_preamble(
    dataset: EOSDataset,
    plan: EOSBatchPlan,
    archive_path: Any,
    options: EOSReportOptions | None = None,
) -> tuple[ReportTable, ...]:
    """Build input and requested-configuration tables before execution.

    Parameters
    ----------
    dataset : EOSDataset
        Parsed and normalized input data.
    plan : EOSBatchPlan
        Declarative ordered batch plan.
    archive_path : path-like
        Native HDF5 destination shown in the input summary.
    options : EOSReportOptions or None, optional
        Reporting detail and input-row selection.

    Returns
    -------
    tuple of ReportTable
        Input identity, data, optional uncertainties, plan summary, and one
        normalized request table per job.
    """
    resolved = EOSReportOptions() if options is None else options
    tables: list[ReportTable] = [
        eos_input_summary_table(dataset, archive_path),
        eos_data_table(dataset, max_rows=resolved.max_data_rows),
    ]
    if dataset.groups is not None or dataset.excluded_npoints:
        tables.append(eos_data_selection_table(dataset))
    if resolved.show_uncertainties:
        tables.append(eos_uncertainty_table(dataset, max_rows=resolved.max_data_rows))
    tables.append(eos_plan_table_from_plan(plan))
    tables.extend(
        eos_requested_fit_table(job, index)
        for index, job in enumerate(plan.jobs, start=1)
    )
    return tuple(tables)


def build_eos_batch_result_tables(
    result: EOSBatchResult,
    options: EOSReportOptions | None = None,
) -> tuple[ReportTable, ...]:
    """Build detailed results followed by compact batch summaries."""
    resolved = EOSReportOptions() if options is None else options
    tables: list[ReportTable] = []
    for job in result.jobs:
        tables.extend(
            eos_job_tables(
                job,
                dataset=result.dataset,
                detail=resolved.detail,
                debug=resolved.debug,
            )
        )
    if result.jobs:
        tables.extend(
            (
                eos_batch_fit_summary_table(result),
                eos_batch_parameter_summary_table(result),
            )
        )
    return tuple(tables)


def eos_batch_fit_summary_table(result: EOSBatchResult) -> ReportTable:
    """Return a compact one-row-per-job summary sorted by scientific identity."""
    domain_order = {"pv": 0, "vt": 1, "pvt": 2, "ev": 3}
    jobs = sorted(
        result.jobs,
        key=lambda job: (
            domain_order.get(job.request.domain.value, 99),
            job.request.target,
            model_label(job.request.model),
            solver_label(_solver_method(job.request)),
            job.job_id,
        ),
    )
    rows: list[list[Any]] = []
    for job in jobs:
        fit = job.result.fit
        diagnostics = fit.diagnostics
        rows.append(
            [
                domain_label(job.request.domain),
                target_label(job.request.target),
                model_label(job.request.model),
                solver_label(_solver_method(job.request)),
                job.job_id,
                fit.status.value,
                job.accepted,
                fit.n_points,
                fit.n_parameters,
                fit.rmse,
                None if diagnostics is None else diagnostics.reduced_chi_square,
                fit.max_abs_error,
                None if diagnostics is None else diagnostics.condition_number,
            ]
        )
    return ReportTable(
        "EOS batch fit summary",
        [
            "Domain",
            "Quantity",
            "Formulation",
            "Solver",
            "Job",
            "Status",
            "Accepted",
            "Points",
            "Free",
            "RMSE",
            "Reduced χ²",
            "Max residual",
            "Condition number",
        ],
        rows,
        metadata={
            "column_formats": [
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                "eos_statistic",
                "eos_statistic",
                "eos_statistic",
                "eos_statistic",
            ],
            "column_alignments": [
                "left",
                "left",
                "left",
                "left",
                "left",
                "left",
                "center",
                "right",
                "right",
                "right",
                "right",
                "right",
                "right",
            ],
            "notes": ["Rows are ordered by domain, quantity, formulation, and solver."],
        },
    )


def eos_batch_parameter_summary_table(result: EOSBatchResult) -> ReportTable:
    """Return a paper-ready long-form summary of all reported parameters."""
    domain_order = {"pv": 0, "vt": 1, "pvt": 2, "ev": 3}
    jobs = sorted(
        result.jobs,
        key=lambda job: (
            domain_order.get(job.request.domain.value, 99),
            job.request.target,
            model_label(job.request.model),
            solver_label(_solver_method(job.request)),
            job.job_id,
        ),
    )
    rows: list[list[Any]] = []
    for job in jobs:
        fit = job.result.fit
        if fit.parameters is None:
            continue
        definitions = fit.metadata.get("parameter_map", {}).get("definitions", [])
        units = {
            item.get("name"): format_unit(item.get("unit"))
            for item in definitions
            if isinstance(item, dict)
        }
        errors = (
            fit.errors
            if fit.errors is not None
            else np.full(len(fit.parameter_names), np.nan)
        )
        for index, name in enumerate(fit.parameter_names):
            state = (
                fit.parameter_states[index].value
                if index < len(fit.parameter_states)
                else "unknown"
            )
            esd = None if not np.isfinite(errors[index]) else float(errors[index])
            rows.append(
                [
                    domain_label(job.request.domain),
                    target_label(job.request.target),
                    model_label(job.request.model),
                    solver_label(_solver_method(job.request)),
                    job.job_id,
                    job.accepted,
                    parameter_label(name),
                    state,
                    float(fit.parameters[index]),
                    esd,
                    units.get(name) or "—",
                ]
            )
    return ReportTable(
        "EOS batch parameter summary",
        [
            "Domain",
            "Quantity",
            "Formulation",
            "Solver",
            "Job",
            "Accepted",
            "Parameter",
            "State",
            "Value",
            "E.S.D.",
            "Unit",
        ],
        rows,
        metadata={
            "column_formats": [
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                "eos_parameter",
                "eos_uncertainty",
                None,
            ],
            "column_alignments": [
                "left",
                "left",
                "left",
                "left",
                "left",
                "center",
                "left",
                "left",
                "right",
                "right",
                "left",
            ],
            "notes": [
                "Machine-readable parameter names remain unchanged in HDF5 and CSV metadata."
            ],
        },
    )


def build_eos_batch_report(
    result: EOSBatchResult,
    options: EOSReportOptions | None = None,
) -> tuple[ReportTable, ...]:
    """Build the complete ordered report for one completed batch."""
    return (
        *build_eos_batch_preamble(
            result.dataset, result.plan, result.archive_path, options
        ),
        *build_eos_batch_result_tables(result, options),
    )


def eos_input_summary_table(dataset: EOSDataset, archive_path: Any) -> ReportTable:
    """Return input identity, size, units, and uncertainty availability."""
    measured = [
        name
        for name in ("pressure", "temperature", "volume", "a", "b", "c")
        if dataset.has(name)
    ]
    uncertain = [name for name in measured if dataset.has(f"sigma_{name}")]
    rows: list[list[Any]] = [
        ["Title", dataset.jobname],
        ["Source", str(dataset.source or "in-memory")],
        ["Observations", dataset.npoints],
        ["Selected by default", dataset.selected_npoints],
        ["Excluded by default", dataset.excluded_npoints],
        ["Data groups", len(dataset.group_ids) or "none"],
        [
            "Crystal system",
            "unknown"
            if dataset.crystal_system is None
            else dataset.crystal_system.value,
        ],
        ["Available quantities", ", ".join(measured) or "none"],
        ["Standard uncertainties", ", ".join(uncertain) or "none"],
        ["Output archive", str(archive_path)],
    ]
    return ReportTable("EOS input summary", ["Property", "Value"], rows)


def eos_data_table(dataset: EOSDataset, *, max_rows: int | None = None) -> ReportTable:
    """Return only input quantities actually present in the dataset."""
    available = [
        name
        for name in ("pressure", "temperature", "volume", "a", "b", "c")
        if dataset.has(name)
    ]
    count = dataset.npoints if max_rows is None else min(dataset.npoints, max_rows)
    include_selection = dataset.groups is not None or dataset.excluded_npoints > 0
    columns = [target_label(name) for name in available]
    if include_selection:
        columns.extend(["Group", "Use"])
    rows: list[list[Any]] = []
    default_mask = dataset.selection_mask()
    for index in range(count):
        row: list[Any] = [_column_value(dataset, name, index) for name in available]
        if include_selection:
            group = None if dataset.groups is None else int(dataset.groups[index])
            row.extend([group, bool(default_mask[index])])
        rows.append(row)
    uncertainty_names = [name for name in available if dataset.has(f"sigma_{name}")]
    notes = [
        "Standard uncertainties available for: "
        + (
            ", ".join(target_label(name) for name in uncertainty_names)
            if uncertainty_names
            else "none"
        ),
        f"Default selection: {dataset.selected_npoints} included, "
        f"{dataset.excluded_npoints} excluded",
    ]
    if count < dataset.npoints:
        notes.append(f"Showing {count} of {dataset.npoints} observations")
    format_by_name = {
        "pressure": "eos_pressure",
        "temperature": "eos_temperature",
        "volume": "eos_structural",
        "a": "eos_structural",
        "b": "eos_structural",
        "c": "eos_structural",
    }
    formats: list[str | None] = [format_by_name[name] for name in available]
    alignments = ["right"] * len(available)
    units = [format_unit(dataset.units.get(name, "")) or "" for name in available]
    if include_selection:
        formats.extend([None, None])
        alignments.extend(["right", "center"])
        units.extend(["", ""])
    return ReportTable(
        "EOS input data",
        columns,
        rows,
        metadata={
            "column_units": units,
            "column_formats": formats,
            "column_alignments": alignments,
            "notes": notes,
        },
    )


def eos_data_selection_table(
    dataset: EOSDataset,
    mask: np.ndarray | None = None,
) -> ReportTable:
    """Return non-destructive selection counts overall and by data group."""
    selected = dataset.selection_mask(mask)
    rows: list[list[Any]] = [
        [
            "All data",
            dataset.npoints,
            int(np.count_nonzero(selected)),
            int(np.count_nonzero(~selected)),
        ]
    ]
    for item in dataset.group_summary(selected):
        rows.append(
            [
                f"Group {item['group']}",
                item["total"],
                item["selected"],
                item["excluded"],
            ]
        )
    return ReportTable(
        "EOS data selection",
        ["Subset", "Available", "Selected", "Excluded"],
        rows,
        metadata={
            "column_alignments": ["left", "right", "right", "right"],
            "notes": [
                "Exclusions are non-destructive; every observation remains in the HDF5 archive."
            ],
        },
    )


def eos_uncertainty_table(
    dataset: EOSDataset, *, max_rows: int | None = None
) -> ReportTable:
    """Return standard uncertainties for measured quantities that exist."""
    available = [
        name
        for name in ("pressure", "temperature", "volume", "a", "b", "c")
        if dataset.has(f"sigma_{name}")
    ]
    count = dataset.npoints if max_rows is None else min(dataset.npoints, max_rows)
    rows = [
        [_column_value(dataset, f"sigma_{name}", index) for name in available]
        for index in range(count)
    ]
    notes = []
    if count < dataset.npoints:
        notes.append(f"Showing {count} of {dataset.npoints} observations")
    formats = [
        "eos_pressure_uncertainty"
        if name == "pressure"
        else "eos_temperature_uncertainty"
        if name == "temperature"
        else "eos_uncertainty"
        for name in available
    ]
    return ReportTable(
        "EOS input standard uncertainties",
        [f"σ({target_label(name)})" for name in available],
        rows,
        metadata={
            "column_units": [
                format_unit(dataset.units.get(f"sigma_{name}", "")) or ""
                for name in available
            ],
            "column_formats": formats,
            "column_alignments": ["right"] * len(available),
            "notes": notes,
        },
    )


def eos_plan_table(result: EOSBatchResult) -> ReportTable:
    """Return one row per declarative batch job in a completed result."""
    return eos_plan_table_from_plan(result.plan)


def eos_plan_table_from_plan(plan: EOSBatchPlan) -> ReportTable:
    """Return one human-readable row per declarative batch job."""
    rows = []
    for index, job in enumerate(plan.jobs, start=1):
        rows.append(
            [
                index,
                job.job_id or job.request.request_id or f"job-{index:06d}",
                domain_label(job.request.domain),
                target_label(job.request.target),
                model_label(job.request.model),
                solver_label(_solver_method(job.request)),
                None
                if job.request.mask is None
                else int(np.count_nonzero(job.request.mask)),
                job.accept,
                job.replace_accepted,
            ]
        )
    return ReportTable(
        "EOS batch plan",
        [
            "#",
            "Job",
            "Domain",
            "Quantity",
            "Formulation",
            "Solver",
            "Selected",
            "Accept",
            "Replace",
        ],
        rows,
        metadata={
            "column_alignments": [
                "right",
                "left",
                "left",
                "left",
                "left",
                "left",
                "right",
                "center",
                "center",
            ]
        },
    )


def eos_requested_fit_table(job: EOSBatchJob, index: int) -> ReportTable:
    """Return normalized model, solver, and constraint settings for one job."""
    request = job.request
    options = request.options.as_dict()
    constraints = request.as_dict().get("constraints", [])
    constraint_text = (
        "; ".join(
            f"{parameter_label(item['name'])}="
            f"{item.get('value') if item['state'] == 'fixed' else item.get('initial_value')} "
            f"({item['state']})"
            for item in constraints
        )
        or "none"
    )
    solver_options = options.get("solver_options", {})
    rows: list[list[Any]] = [
        ["Job", job.job_id or request.request_id or f"job-{index:06d}"],
        ["Domain", domain_label(request.domain)],
        ["Quantity", target_label(request.target)],
        ["Formulation", model_label(request.model)],
        ["Solver", solver_label(_solver_method(request))],
        ["Solver options", _compact_mapping(solver_options)],
        ["Parameter constraints", constraint_text],
        ["Data selection", _compact_mapping(request.metadata.get("selection", {}))],
        ["Accept successful result", job.accept],
        ["Replace accepted result", job.replace_accepted],
    ]
    return ReportTable(
        f"EOS requested fit — {job.job_id or index}",
        ["Property", "Value"],
        rows,
    )


def eos_job_tables(
    job: EOSBatchJobResult,
    *,
    dataset: EOSDataset,
    detail: EOSReportDetail = EOSReportDetail.SHORT,
    debug: bool = False,
) -> tuple[ReportTable, ...]:
    """Return short or extended tables for one fit result."""
    tables: list[ReportTable] = [
        eos_job_configuration_table(job),
        eos_parameter_table(job.result),
        eos_diagnostics_table(job, debug=debug),
    ]
    if debug:
        tables.extend(eos_solver_debug_tables(job))
    if EOSReportDetail(detail) is EOSReportDetail.EXTENDED and job.result.fit.success:
        tables.append(eos_observed_calculated_table(job.result, dataset))
        if job.result.fit.covariance is not None:
            tables.append(
                _matrix_table(
                    "Parameter covariance",
                    job.result.fit.parameter_names,
                    job.result.fit.covariance,
                )
            )
        correlation = (
            None
            if job.result.fit.diagnostics is None
            else job.result.fit.diagnostics.correlation
        )
        if correlation is not None:
            tables.append(
                _matrix_table(
                    "Parameter correlation", job.result.fit.parameter_names, correlation
                )
            )
    return tuple(tables)


def eos_job_configuration_table(job: EOSBatchJobResult) -> ReportTable:
    """Return normalized human-readable model and solver configuration."""
    request = job.request
    rows: list[list[Any]] = [
        ["Job", job.job_id],
        ["Domain", domain_label(request.domain)],
        ["Quantity", target_label(request.target)],
        ["Formulation", model_label(request.model)],
        ["Solver", solver_label(_solver_method(request))],
        [
            "Selected observations",
            "all" if request.mask is None else int(np.count_nonzero(request.mask)),
        ],
        ["Record ID", job.record_id],
        ["Accepted", job.accepted],
        ["Status", job.result.fit.status.value],
    ]
    return ReportTable(
        f"EOS fit configuration — {job.job_id}", ["Property", "Value"], rows
    )


def eos_parameter_table(result: EOSFitResult) -> ReportTable:
    """Return parameter values with scientific symbols and explicit units."""
    fit = result.fit
    parameter_map = fit.metadata.get("parameter_map", {})
    definitions = (
        parameter_map.get("definitions", []) if isinstance(parameter_map, dict) else []
    )
    initial = {
        item.get("name"): item.get("initial_value")
        for item in definitions
        if isinstance(item, dict)
    }
    bounds = {
        item.get("name"): item.get("bounds")
        for item in definitions
        if isinstance(item, dict)
    }
    units = {
        item.get("name"): format_unit(item.get("unit"))
        for item in definitions
        if isinstance(item, dict)
    }
    rows: list[list[Any]] = []
    parameters = (
        fit.parameters
        if fit.parameters is not None
        else np.full(len(fit.parameter_names), np.nan)
    )
    errors = (
        fit.errors
        if fit.errors is not None
        else np.full(len(fit.parameter_names), np.nan)
    )
    for index, name in enumerate(fit.parameter_names):
        value = float(parameters[index])
        start = initial.get(name)
        shift = None if start is None else value - float(start)
        esd = None if not np.isfinite(errors[index]) else float(errors[index])
        state = (
            fit.parameter_states[index].value
            if index < len(fit.parameter_states)
            else "unknown"
        )
        rows.append(
            [
                parameter_label(name),
                state,
                start,
                value,
                shift,
                esd,
                units.get(name) or "—",
                bounds.get(name),
            ]
        )
    if not rows:
        rows.append(["—", "—", None, None, None, None, "—", None])
    return ReportTable(
        "EOS parameters",
        [
            "Parameter",
            "State",
            "Initial",
            "Value",
            "Final - initial",
            "E.S.D.",
            "Unit",
            "Bounds",
        ],
        rows,
        metadata={
            "column_formats": [
                None,
                None,
                "eos_parameter",
                "eos_parameter",
                "eos_parameter",
                "eos_uncertainty",
                None,
                None,
            ],
            "column_alignments": [
                "left",
                "left",
                "right",
                "right",
                "right",
                "right",
                "left",
                "left",
            ],
        },
    )


def eos_diagnostics_table(
    job: EOSBatchJobResult, *, debug: bool = False
) -> ReportTable:
    """Return objective and objective diagnostics without subjective judgement."""
    fit = job.result.fit
    diagnostics = fit.diagnostics
    values: list[list[Any]] = [
        ["Success", fit.success],
        ["Status", fit.status.value],
        ["Message", fit.message],
        ["Points", fit.n_points],
        ["Free parameters", fit.n_parameters],
        ["Degrees of freedom", fit.dof],
        ["RMSE", fit.rmse],
        ["MAE", fit.mae],
        ["Maximum absolute residual", fit.max_abs_error],
        [
            "Reduced chi-square",
            None if diagnostics is None else diagnostics.reduced_chi_square,
        ],
        [
            "Condition number",
            None if diagnostics is None else diagnostics.condition_number,
        ],
        ["Stop reason", None if diagnostics is None else diagnostics.stop_reason],
    ]
    notes = list(job.result.warnings) + list(fit.warnings)
    if diagnostics is not None:
        notes.extend(diagnostics.warnings)
    if not fit.success:
        notes.append(
            "Detailed solver diagnostics follow."
            if debug
            else "Rerun the EOS batch with -v debug for detailed solver diagnostics."
        )
    return ReportTable(
        f"EOS diagnostics — {job.job_id}",
        ["Metric", "Value"],
        values,
        metadata={
            "column_formats": [None, "eos_statistic"],
            "notes": tuple(dict.fromkeys(note for note in notes if note)),
        },
    )


def eos_solver_debug_tables(job: EOSBatchJobResult) -> tuple[ReportTable, ...]:
    """Return detailed backend-neutral diagnostics for one EOS solver result.

    Parameters
    ----------
    job : EOSBatchJobResult
        Completed or failed EOS batch job.

    Returns
    -------
    tuple of ReportTable
        Summary, numerical ranges, parameter state, and available iteration or
        model-evaluation traces.
    """
    fit = job.result.fit
    diagnostics = fit.diagnostics
    if diagnostics is None:
        return (
            ReportTable(
                f"EOS solver debug — {job.job_id}",
                ["Property", "Value"],
                [["Diagnostics", "No solver diagnostics were returned"]],
            ),
        )
    metadata = diagnostics.metadata
    tables: list[ReportTable] = [
        _solver_debug_summary_table(job, diagnostics, metadata)
    ]
    range_table = _solver_debug_ranges_table(job, metadata)
    if range_table is not None:
        tables.append(range_table)
    parameter_table = _solver_debug_parameter_table(job, metadata)
    if parameter_table is not None:
        tables.append(parameter_table)
    history = metadata.get("effective_variance_history")
    if history is None:
        history = metadata.get("reweighting_history")
    if isinstance(history, list) and history:
        tables.append(_effective_variance_history_table(job, history))
    trace = metadata.get("evaluation_trace")
    if isinstance(trace, list) and trace:
        tables.append(_solver_evaluation_trace_table(job, trace, metadata))
    return tuple(tables)


def _solver_debug_summary_table(
    job: EOSBatchJobResult, diagnostics: Any, metadata: dict[str, Any]
) -> ReportTable:
    """Return termination, backend, and work-counter diagnostics."""
    fit = job.result.fit
    rows = [
        ["Method", None if fit.method is None else fit.method.value],
        ["Backend", metadata.get("backend")],
        ["Backend version", metadata.get("backend_version")],
        ["Termination category", metadata.get("termination_category")],
        ["Backend status code", metadata.get("backend_status_code")],
        ["Exception type", metadata.get("exception_type")],
        ["Exception message", metadata.get("exception_message")],
        ["Objective", diagnostics.objective],
        ["Weighted objective", diagnostics.weighted],
        ["Iterations", diagnostics.n_iterations],
        ["Function/model evaluations", diagnostics.n_evaluations],
        ["Jacobian evaluations", metadata.get("jacobian_evaluations")],
        ["Jacobian rank", diagnostics.jacobian_rank],
        ["Condition number", diagnostics.condition_number],
        ["Stop reason", diagnostics.stop_reason],
        ["Trace requested", metadata.get("detailed_trace_requested")],
        ["Trace truncated", metadata.get("evaluation_trace_truncated")],
        [
            "Period-two oscillation",
            metadata.get("period_two_oscillation_detected"),
        ],
        ["Last valid iteration", metadata.get("last_valid_iteration")],
        ["Last parameter change / tolerance", metadata.get("last_parameter_change")],
        ["Last sigma change / tolerance", metadata.get("last_sigma_change")],
        ["Last inner stop reason", metadata.get("last_inner_stop_reason")],
    ]
    return ReportTable(
        f"EOS solver debug summary — {job.job_id}",
        ["Property", "Value"],
        rows,
        metadata={"column_formats": [None, "eos_statistic"]},
    )


def _solver_debug_ranges_table(
    job: EOSBatchJobResult, metadata: dict[str, Any]
) -> ReportTable | None:
    """Return compact ranges for data and uncertainty arrays."""
    rows: list[list[Any]] = []
    for key, label in (
        ("x_summary", "Explanatory coordinates"),
        ("y_summary", "Observed response"),
        ("sigma_summary", "Objective sigma"),
        ("sigma_x_summary", "Independent-variable sigma"),
        ("sigma_y_summary", "Dependent-variable sigma"),
    ):
        values = metadata.get(key)
        if not isinstance(values, dict):
            continue
        rows.append(
            [
                label,
                values.get("shape"),
                values.get("minimum"),
                values.get("maximum"),
                values.get("median"),
                values.get("norm"),
                values.get("positive_dynamic_range"),
            ]
        )
    if not rows:
        return None
    return ReportTable(
        f"EOS solver numerical ranges — {job.job_id}",
        ["Quantity", "Shape", "Minimum", "Maximum", "Median", "Norm", "Max/min"],
        rows,
        metadata={
            "column_formats": [None, None, *(["eos_statistic"] * 5)],
            "column_alignments": ["left", "left", *(["right"] * 5)],
        },
    )


def _solver_debug_parameter_table(
    job: EOSBatchJobResult, metadata: dict[str, Any]
) -> ReportTable | None:
    """Return initial, last evaluated, and bound values for free parameters."""
    fit = job.result.fit
    initial = metadata.get("initial_parameters")
    lower = metadata.get("lower_bounds")
    upper = metadata.get("upper_bounds")
    last = metadata.get("last_valid_free_parameters")
    if last is None:
        evaluation = metadata.get("last_evaluation")
        if isinstance(evaluation, dict):
            last = evaluation.get("parameters")
    if last is None and fit.optimizer_parameters is not None:
        last = fit.optimizer_parameters.tolist()
    if not isinstance(initial, list):
        return None
    size = len(initial)
    names = _free_parameter_names(fit, metadata, size)
    lower_values = (
        lower if isinstance(lower, list) and len(lower) == size else [None] * size
    )
    upper_values = (
        upper if isinstance(upper, list) and len(upper) == size else [None] * size
    )
    last_values = (
        last if isinstance(last, list) and len(last) == size else [None] * size
    )
    rows = []
    for index, name in enumerate(names):
        start = initial[index]
        end = last_values[index]
        shift = None if end is None else float(end) - float(start)
        rows.append([name, start, end, shift, lower_values[index], upper_values[index]])
    return ReportTable(
        f"EOS solver parameter path — {job.job_id}",
        ["Parameter", "Initial", "Last evaluated", "Shift", "Lower", "Upper"],
        rows,
        metadata={
            "column_formats": [None, *(["eos_parameter"] * 5)],
            "column_alignments": ["left", *(["right"] * 5)],
        },
    )


def _free_parameter_names(fit: Any, metadata: dict[str, Any], size: int) -> list[str]:
    """Return free parameter names matching one optimizer vector."""
    explicit = metadata.get("free_parameter_names")
    if isinstance(explicit, list) and len(explicit) == size:
        return [str(name) for name in explicit]
    parameter_map = fit.metadata.get("parameter_map")
    if isinstance(parameter_map, dict):
        definitions = parameter_map.get("definitions")
        if isinstance(definitions, list):
            names = [
                str(item.get("name"))
                for item in definitions
                if isinstance(item, dict) and item.get("state") == "free"
            ]
            if len(names) == size:
                return names
    if len(fit.parameter_names) == size:
        return list(fit.parameter_names)
    return [f"p{index}" for index in range(size)]


def _effective_variance_history_table(
    job: EOSBatchJobResult, history: list[dict[str, Any]]
) -> ReportTable:
    """Return one row per effective-variance reweighting cycle."""
    rows = [
        [
            item.get("iteration"),
            item.get("parameter_change"),
            item.get("sigma_change"),
            item.get("chi_square"),
            item.get("rmse"),
            item.get("maximum_absolute_residual"),
            item.get("sigma_minimum"),
            item.get("sigma_maximum"),
            item.get("sigma_dynamic_range"),
            item.get("inner_evaluations"),
            item.get("inner_stop_reason"),
        ]
        for item in history
        if isinstance(item, dict)
    ]
    return ReportTable(
        f"Effective-variance iteration history — {job.job_id}",
        [
            "Cycle",
            "Δ parameters / tol",
            "Δ sigma / tol",
            "Chi-square",
            "RMSE",
            "Max residual",
            "Sigma min",
            "Sigma max",
            "Sigma max/min",
            "Inner evals",
            "Inner stop reason",
        ],
        rows,
        metadata={
            "column_formats": [
                None,
                *(["eos_statistic"] * 8),
                None,
                None,
            ],
            "column_alignments": ["right"] * 10 + ["left"],
        },
    )


def _solver_evaluation_trace_table(
    job: EOSBatchJobResult, trace: list[dict[str, Any]], metadata: dict[str, Any]
) -> ReportTable:
    """Return the bounded model-evaluation trace retained by OLS or WLS."""
    rows = [
        [
            item.get("evaluation"),
            item.get("objective"),
            item.get("rmse"),
            item.get("maximum_absolute_residual"),
            _compact_vector(item.get("parameters")),
        ]
        for item in trace
        if isinstance(item, dict)
    ]
    notes = []
    if metadata.get("evaluation_trace_truncated"):
        notes.append(
            "The trace contains the first and most-recent evaluations; "
            "intermediate evaluations were truncated."
        )
    return ReportTable(
        f"Solver model-evaluation trace — {job.job_id}",
        ["Evaluation", "Objective", "RMSE", "Max residual", "Parameters"],
        rows,
        metadata={
            "column_formats": [
                None,
                "eos_statistic",
                "eos_statistic",
                "eos_statistic",
                None,
            ],
            "column_alignments": ["right", "right", "right", "right", "left"],
            "notes": notes,
        },
    )


def _compact_vector(values: Any) -> str:
    """Return compact deterministic text for a parameter vector."""
    if not isinstance(values, list):
        return "—"
    return "[" + ", ".join(f"{float(value):.8g}" for value in values) + "]"


def eos_observed_calculated_table(
    result: EOSFitResult,
    dataset: EOSDataset,
) -> ReportTable:
    """Return observed, calculated, and residual data in domain orientation."""
    fit = result.fit
    fitted = np.asarray(fit.fitted if fit.fitted is not None else [], dtype=float)
    residuals = np.asarray(
        fit.residuals if fit.residuals is not None else [], dtype=float
    )
    request = result.request
    mask = dataset.selection_mask(request.mask)
    if request.domain.value == "pv":
        coordinate = np.asarray(dataset.column(request.target)[mask], dtype=float)
        pressure = np.asarray(dataset.column("pressure")[mask], dtype=float)
        columns = [
            request.target,
            "Pressure observed",
            "Pressure calculated",
            "Residual",
        ]
        rows = [
            [coordinate[i], pressure[i], fitted[i], residuals[i]]
            for i in range(fitted.size)
        ]
        formats = [
            "eos_structural",
            "eos_pressure",
            "eos_pressure",
            "eos_residual",
        ]
    elif request.domain.value == "vt":
        temperature = np.asarray(dataset.column("temperature")[mask], dtype=float)
        observed = np.asarray(dataset.column(request.target)[mask], dtype=float)
        calculated = np.asarray(
            result.predictions.get(request.target, fitted), dtype=float
        )[mask]
        columns = [
            "Temperature",
            f"{request.target} observed",
            f"{request.target} calculated",
            "Residual",
        ]
        rows = [
            [temperature[i], observed[i], calculated[i], observed[i] - calculated[i]]
            for i in range(temperature.size)
        ]
        formats = [
            "eos_temperature",
            "eos_structural",
            "eos_structural",
            "eos_residual",
        ]
    else:
        volume = np.asarray(dataset.column("volume")[mask], dtype=float)
        temperature = np.asarray(dataset.column("temperature")[mask], dtype=float)
        pressure = np.asarray(dataset.column("pressure")[mask], dtype=float)
        columns = [
            "Volume",
            "Temperature",
            "Pressure observed",
            "Pressure calculated",
            "Residual",
        ]
        rows = [
            [volume[i], temperature[i], pressure[i], fitted[i], residuals[i]]
            for i in range(fitted.size)
        ]
        formats = [
            "eos_structural",
            "eos_temperature",
            "eos_pressure",
            "eos_pressure",
            "eos_residual",
        ]
    return ReportTable(
        "Observed and calculated EOS data",
        columns,
        rows,
        metadata={
            "column_formats": formats,
            "column_alignments": ["right"] * len(columns),
        },
    )


def _prediction_or_index(result: EOSFitResult, name: str, size: int) -> np.ndarray:
    values = result.predictions.get(name)
    if values is None or np.asarray(values).size != size:
        return np.arange(size, dtype=float)
    return np.asarray(values, dtype=float)


def _matrix_table(
    title: str, names: tuple[str, ...], matrix: np.ndarray
) -> ReportTable:
    rows = [
        [name, *[float(value) for value in row]]
        for name, row in zip(names, np.asarray(matrix), strict=True)
    ]
    numeric_format = (
        "eos_covariance" if title == "Parameter covariance" else "eos_correlation"
    )
    return ReportTable(
        title,
        ["Parameter", *names],
        rows,
        metadata={
            "column_formats": [None, *([numeric_format] * len(names))],
            "column_alignments": ["left", *(["right"] * len(names))],
        },
    )


def _solver_method(request: Any) -> str:
    """Return the normalized method name from one EOS request."""
    solver = request.options.solver_options
    if solver is not None:
        return solver.method.value
    method = request.options.method
    return getattr(method, "value", str(method))


def _compact_mapping(values: Any) -> str:
    """Return deterministic compact text for nested solver-option mappings."""
    if not isinstance(values, dict):
        return str(values)
    parts = []
    for key in sorted(values):
        value = values[key]
        if isinstance(value, dict):
            rendered = _compact_mapping(value)
        else:
            rendered = str(value)
        parts.append(f"{key}={rendered}")
    return ", ".join(parts)


def _column_value(dataset: EOSDataset, name: str, index: int) -> float | None:
    if not dataset.has(name):
        return None
    return float(dataset.column(name)[index])


__all__ = [
    "EOSReportDetail",
    "EOSReportOptions",
    "build_eos_batch_preamble",
    "build_eos_batch_report",
    "build_eos_batch_result_tables",
    "eos_batch_fit_summary_table",
    "eos_batch_parameter_summary_table",
    "eos_data_selection_table",
    "eos_data_table",
    "eos_solver_debug_tables",
    "eos_uncertainty_table",
]
