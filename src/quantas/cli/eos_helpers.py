# -*- coding: utf-8 -*-

"""Private construction helpers shared by EOS Click commands."""

from __future__ import annotations

from typing import Any, cast

import click
import numpy as np
from click.core import ParameterSource

from quantas.cli.contracts import ReportVerbosity
from quantas.cli.output import CLIOutput
from quantas.core.events import Event, EventLevel
from quantas.core.math.fitting import (
    CovarianceScaling,
    EffectiveVarianceOptions,
    ODRDifferenceScheme,
    OLSOptions,
    OrthogonalDistanceOptions,
    SOLVER_DEBUG_METADATA_KEY,
    WLSOptions,
)
from quantas.core.physics.eos import (
    MGDNormalization,
    PVTModel,
    ThermalPressureFamily,
    parse_eos_model,
    parse_temperature_eos_model,
    parse_thermal_pressure_model,
)
from quantas.api.eos import (
    BatchPlan as EOSBatchPlan,
    FitDomain as EOSFitDomain,
    FitOptions as EOSFitOptions,
    FitRequest as EOSFitRequest,
    ParameterConstraint,
    ReportDetail as EOSReportDetail,
    ReportOptions as EOSReportOptions,
    validate_request,
)


def _combined_grid(
    values: tuple[float, ...], range_spec: str | None, name: str
) -> np.ndarray | None:
    explicit = np.asarray(values, dtype=np.float64)
    ranged = None if range_spec is None else _parse_numeric_range(range_spec, name)
    if explicit.size == 0 and ranged is None:
        return None
    if explicit.size and ranged is not None:
        return np.concatenate((explicit, ranged))
    return explicit if explicit.size else ranged


def _parse_numeric_range(specification: str, name: str) -> np.ndarray:
    parts = str(specification).split(":")
    if len(parts) != 3:
        raise click.BadParameter(
            f"{name} range must use START:STOP:STEP", param_hint=f"--{name}-range"
        )
    try:
        start, stop, step = (float(value) for value in parts)
    except ValueError as exc:
        raise click.BadParameter(
            f"{name} range contains a non-numeric value",
            param_hint=f"--{name}-range",
        ) from exc
    if not all(np.isfinite([start, stop, step])) or step == 0.0:
        raise click.BadParameter(
            f"{name} range values must be finite and STEP non-zero",
            param_hint=f"--{name}-range",
        )
    if (stop - start) * step < 0.0:
        raise click.BadParameter(
            f"{name} range STEP points away from STOP",
            param_hint=f"--{name}-range",
        )
    count = int(np.floor((stop - start) / step + 1.0e-12)) + 1
    if count <= 0 or count > 1_000_000:
        raise click.BadParameter(
            f"{name} range produces an invalid number of points",
            param_hint=f"--{name}-range",
        )
    grid = start + step * np.arange(count, dtype=np.float64)
    tolerance = 1.0e-12 * max(abs(start), abs(stop), 1.0)
    if abs(grid[-1] - stop) > tolerance and (stop - grid[-1]) * step > 0.0:
        grid = np.append(grid, stop)
    return grid


def _combine_pvt_coordinates(
    primary: np.ndarray,
    temperature: np.ndarray,
    *,
    pairwise: bool,
) -> tuple[np.ndarray, np.ndarray]:
    if pairwise:
        if (
            primary.size != temperature.size
            and primary.size != 1
            and temperature.size != 1
        ):
            raise click.UsageError(
                "--pairwise requires equal lengths or one scalar coordinate"
            )
        first, second = np.broadcast_arrays(primary, temperature)
        return first.reshape(-1), second.reshape(-1)
    first, second = np.meshgrid(primary, temperature, indexing="ij")
    return first.reshape(-1), second.reshape(-1)


_SPEC_SCIENTIFIC_PARAMETERS = (
    "pressure_unit",
    "length_unit",
    "temperature_unit",
    "domain",
    "targets",
    "pv_eos",
    "pv_order",
    "vt_eos",
    "vt_variant",
    "pvt_coupling",
    "thermal_pressure_model",
    "atoms_per_cell",
    "formula",
    "formula_units_per_cell",
    "fixed_parameters",
    "initial_parameters",
    "parameter_bounds",
    "solver",
    "max_iterations",
    "ftol",
    "xtol",
    "gtol",
    "covariance_scaling",
    "inner_max_iterations",
    "odr_difference",
    "odr_ndigit",
    "failure_policy",
)


def _reject_spec_cli_conflicts(ctx: click.Context) -> None:
    conflicts = [
        name
        for name in _SPEC_SCIENTIFIC_PARAMETERS
        if ctx.get_parameter_source(name) is ParameterSource.COMMANDLINE
    ]
    if conflicts:
        rendered = ", ".join("--" + name.replace("_", "-") for name in conflicts)
        raise click.UsageError(
            "--spec is the authority for scientific batch settings and cannot "
            f"be combined with: {rendered}"
        )


def _merge_spec_report_options(
    ctx: click.Context,
    spec_options: EOSReportOptions,
    *,
    verbosity: ReportVerbosity,
    show_uncertainties: bool,
    max_data_rows: int | None,
) -> EOSReportOptions:
    detail = spec_options.detail
    debug = spec_options.debug
    if ctx.get_parameter_source("verbosity") is ParameterSource.COMMANDLINE:
        detail = (
            EOSReportDetail.EXTENDED
            if verbosity.includes_extended
            else EOSReportDetail.SHORT
        )
        debug = verbosity.includes_debug
    show = spec_options.show_uncertainties
    if (
        ctx.get_parameter_source("show_uncertainties") is ParameterSource.COMMANDLINE
        and show_uncertainties
    ):
        show = True
    maximum = spec_options.max_data_rows
    if ctx.get_parameter_source("max_data_rows") is ParameterSource.COMMANDLINE:
        maximum = max_data_rows
    return EOSReportOptions(
        detail=detail,
        show_uncertainties=show,
        max_data_rows=maximum,
        debug=debug,
    )


def _configure_plan_debug(plan: EOSBatchPlan, enabled: bool) -> None:
    """Enable detailed numerical tracing for every resolved batch request.

    Parameters
    ----------
    plan : EOSBatchPlan
        Resolved frontend-neutral batch plan.
    enabled : bool
        Whether solver debug traces should be collected.
    """
    if not enabled:
        return
    for job in plan.jobs:
        job.request.options.metadata = {
            **job.request.options.metadata,
            SOLVER_DEBUG_METADATA_KEY: True,
        }


def _validate_resolved_plan(dataset: Any, plan: EOSBatchPlan) -> None:
    for index, job in enumerate(plan.jobs, start=1):
        try:
            validate_request(dataset, job.request)
        except Exception as exc:
            job_id = job.job_id or job.request.request_id or f"job-{index:06d}"
            raise click.UsageError(
                f"resolved EOS job {job_id!r} is invalid: {exc}"
            ) from exc


class _EOSCLIEventObserver:
    """Render EOS workflow events without interpreting scientific results."""

    def __init__(self, output: CLIOutput) -> None:
        self.output = output

    def __call__(self, event: Event) -> None:
        if event.level is EventLevel.PROGRESS:
            return
        if event.level is EventLevel.RESULT and event.data.get("kind") == "job_result":
            return
        self.output.message(event.message, level=event.level)


def _resolve_targets(
    dataset: Any, domain: EOSFitDomain, targets: tuple[str, ...]
) -> tuple[str, ...]:
    requested = tuple(dict.fromkeys(value.lower() for value in targets))
    if "all" not in requested:
        resolved = requested
    elif domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
        resolved = ("volume",)
    else:
        resolved = tuple(
            name for name in ("volume", "a", "b", "c") if dataset.has(name)
        )
    if domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE and resolved != ("volume",):
        raise click.UsageError("P-V-T fitting currently supports only --fit volume")
    if not resolved:
        raise click.UsageError("no requested fit target is available in the dataset")
    missing = [name for name in resolved if not dataset.has(name)]
    if missing:
        raise click.UsageError(
            "dataset does not contain target(s): " + ", ".join(missing)
        )
    return resolved


def _build_request(
    target: str,
    domain: EOSFitDomain,
    *,
    pv_eos: str,
    pv_order: int,
    vt_eos: str,
    vt_variant: str | None,
    pvt_coupling: str,
    thermal_pressure_model: str,
    atoms_per_cell: float | None,
    formula: str | None,
    formula_units_per_cell: float | None,
    solver_options: Any,
    constraint_specs: dict[str, list[tuple[str | None, str, Any]]],
) -> EOSFitRequest:
    has_mgd_normalization = any(
        value is not None for value in (atoms_per_cell, formula, formula_units_per_cell)
    )
    if domain is not EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE and (
        thermal_pressure_model.lower() != "holland-powell-einstein"
        or has_mgd_normalization
    ):
        raise click.BadParameter(
            "thermal-pressure model and MGD normalization options require --domain pvt"
        )
    if domain is EOSFitDomain.PRESSURE_VOLUME:
        model: Any = _parse_pressure_model(pv_eos, pv_order)
    elif domain is EOSFitDomain.VOLUME_TEMPERATURE:
        model = parse_temperature_eos_model(vt_eos, vt_variant)
    else:
        thermal_pressure = None
        normalization = None
        if pvt_coupling.lower() == "thermal-pressure":
            thermal = None
            thermal_pressure = parse_thermal_pressure_model(thermal_pressure_model)
            if (
                thermal_pressure.family_name
                is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE
            ):
                normalization = MGDNormalization.cell(
                    atoms_per_cell=atoms_per_cell,
                    formula=formula,
                    formula_units_per_cell=formula_units_per_cell,
                )
            elif has_mgd_normalization:
                raise click.BadParameter(
                    "MGD normalization options require --thermal-pressure-model mgd "
                    "or mgd:q-compromise"
                )
        else:
            thermal = parse_temperature_eos_model(vt_eos, vt_variant)
            if (
                thermal_pressure_model.lower() != "holland-powell-einstein"
                or has_mgd_normalization
            ):
                raise click.BadParameter(
                    "thermal-pressure model and MGD normalization options require "
                    "--pvt-coupling thermal-pressure"
                )
        model = PVTModel(
            pressure_model=_parse_pressure_model(pv_eos, pv_order),
            coupling=pvt_coupling,
            temperature_model=thermal,
            thermal_pressure_model=thermal_pressure,
            mgd_normalization=normalization,
        )
    constraints = _constraints_for_target(target, constraint_specs)
    return EOSFitRequest(
        model=model,
        target=target,
        domain=domain,
        constraints=constraints,
        options=EOSFitOptions(solver_options=solver_options),
        request_id=f"{domain.value}-{target}",
        metadata={"source": "cli"},
    )


def _parse_pressure_model(family: str, order: int) -> Any:
    normalized = str(family).strip().lower()
    if normalized in {"m", "murnaghan"}:
        return parse_eos_model(family)
    return parse_eos_model(family, order=order)


def _build_solver_options(
    solver: str,
    **values: Any,
) -> Any:
    method = solver.lower()
    covariance = values["covariance_scaling"]
    covariance_value = (
        None if covariance is None else CovarianceScaling(covariance.replace("-", "_"))
    )
    common = {
        "covariance_scaling": covariance_value,
        "max_iterations": values["max_iterations"],
        "ftol": values["ftol"],
        "xtol": values["xtol"],
        "gtol": values["gtol"],
    }
    if method == "ols":
        return OLSOptions(**common)
    if method == "wls":
        return WLSOptions(**common)
    if method == "effective-variance":
        return EffectiveVarianceOptions(
            **common,
            inner_max_iterations=values["inner_max_iterations"],
        )
    return OrthogonalDistanceOptions(
        **common,
        difference_scheme=ODRDifferenceScheme(values["odr_difference"].lower()),
        ndigit=values["odr_ndigit"],
    )


def _parse_constraint_specs(
    fixed: tuple[str, ...],
    initial: tuple[str, ...],
    bounds: tuple[str, ...],
) -> dict[str, list[tuple[str | None, str, Any]]]:
    specs: dict[str, list[tuple[str | None, str, Any]]] = {
        "fixed": [],
        "initial": [],
        "bounds": [],
    }
    for value in fixed:
        target, name, parsed = _parse_assignment(value)
        specs["fixed"].append((target, name, parsed))
    for value in initial:
        target, name, parsed = _parse_assignment(value)
        specs["initial"].append((target, name, parsed))
    for value in bounds:
        target, name, text = _parse_assignment(value, parse_value=False)
        parts = str(text).split(":")
        if len(parts) != 2:
            raise click.BadParameter("bounds must use NAME=LOW:HIGH")
        low = _optional_float(parts[0])
        high = _optional_float(parts[1])
        specs["bounds"].append((target, name, (low, high)))
    return specs


def _parse_assignment(
    value: str, *, parse_value: bool = True
) -> tuple[str | None, str, Any]:
    if "=" not in value:
        raise click.BadParameter("parameter assignments must use NAME=VALUE")
    left, right = value.split("=", 1)
    if ":" in left:
        target, name = left.split(":", 1)
        target = target.strip().lower()
        if target not in {"volume", "a", "b", "c"}:
            raise click.BadParameter(f"unknown target prefix {target!r}")
    else:
        target, name = None, left
    name = _canonical_parameter_name(name)
    if not right.strip():
        raise click.BadParameter("parameter assignment has an empty value")
    parsed: Any = float(right) if parse_value else right.strip()
    return target, name, parsed


def _optional_float(value: str) -> float | None:
    text = value.strip().lower()
    if text in {"", "none", "null", "-inf", "+inf", "inf"}:
        return None
    return float(value)


def _canonical_parameter_name(value: str) -> str:
    text = value.strip()
    key = text.lower().replace("_", "").replace("'", "p")
    aliases = {
        "b0": "K0",
        "k0": "K0",
        "kp": "KP",
        "kprime": "KP",
        "kpp": "KPP",
        "kdoubleprime": "KPP",
        "v0": "V0",
        "l0": "L0",
        "m0": "M0",
        "mp": "MP",
        "mpp": "MPP",
        "tref": "temperature_ref",
        "temperatureref": "temperature_ref",
        "thetae": "theta_e",
        "thetad": "theta_d0",
        "thetad0": "theta_d0",
        "debyetemperature": "theta_d0",
        "gamma": "gamma0",
        "gamma0": "gamma0",
        "gruneisen": "gamma0",
        "q": "q",
        "thetasat": "theta_sat",
        "dk0dt": "dK0_dT",
    }
    return aliases.get(key, text)


def _constraints_for_target(
    target: str,
    specs: dict[str, list[tuple[str | None, str, Any]]],
) -> tuple[ParameterConstraint, ...]:
    combined: dict[str, dict[str, Any]] = {}
    for kind, entries in specs.items():
        for prefix, name, value in entries:
            if prefix is not None and prefix != target:
                continue
            state = combined.setdefault(name, {})
            if kind == "fixed":
                state["fixed_value"] = float(value)
                state["initial_value"] = float(value)
            elif kind == "initial":
                state["initial_value"] = float(value)
            else:
                state["bounds"] = value
    definitions: list[ParameterConstraint] = []
    for name, values in combined.items():
        if "fixed_value" in values:
            definitions.append(
                cast(
                    ParameterConstraint,
                    ParameterConstraint.fixed(name, values["fixed_value"]),
                )
            )
            continue
        initial = values.get("initial_value")
        bounds = values.get("bounds")
        if initial is None:
            raise click.BadParameter(
                f"free constraint for {name!r} requires --initial {name}=VALUE"
            )
        lower = (
            float("-inf") if bounds is None or bounds[0] is None else float(bounds[0])
        )
        upper = (
            float("inf") if bounds is None or bounds[1] is None else float(bounds[1])
        )
        definitions.append(
            cast(
                ParameterConstraint,
                ParameterConstraint.free(
                    name,
                    float(initial),
                    lower_bound=lower,
                    upper_bound=upper,
                ),
            )
        )
    return tuple(definitions)
