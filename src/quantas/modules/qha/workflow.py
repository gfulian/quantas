# -*- coding: utf-8 -*-

"""Workflow control for quasi-harmonic approximation calculations.

This module coordinates pressure-temperature grid calculations, local fit
recovery, diagnostic collection, progress reporting, and controlled early-stop
conditions.  It does not render messages and does not depend on command-line or
graphical user-interface objects.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from quantas.core.physics.units import energy_to_pressure
from quantas.modules.qha.analysis import (
    apply_minimum_to_result,
    initialize_result,
    prepare_free_energy_grid,
    pressure_energy_density,
)
from quantas.modules.qha.core.minimization import (
    VolumeMinimumResult,
    VolumeRangeStatus,
    evaluate_fitted_eos_at_pressure,
    evaluate_fitted_polynomial_at_pressure,
    fit_eos_free_energy_model,
    fit_polynomial_free_energy_model,
)
from quantas.core.math.fitting import FitQuality
from quantas.modules.qha.models import (
    QHAFailedPoint,
    QHAFitRecord,
    QHAInput,
    QHAOptions,
    QHAResult,
)

QHAWorkflowEventKind = Literal[
    "started",
    "point_started",
    "fit_record",
    "warning",
    "point_completed",
    "point_failed",
    "stopped",
    "completed",
]

QHAWorkflowCallback = Callable[["QHAWorkflowEvent"], None]


@dataclass(slots=True)
class QHAWorkflowEvent:
    """Event emitted by a QHA workflow controller.

    Parameters
    ----------
    kind : str
        Type of workflow event.
    message : str
        Human-readable event description.
    progress : float or None, optional
        Fractional workflow progress between 0 and 1.
    temperature : float or None, optional
        Temperature associated with the event.
    pressure : float or None, optional
        Pressure associated with the event.
    data : dict, optional
        Structured event payload.
    """

    kind: QHAWorkflowEventKind
    message: str
    progress: float | None = None
    temperature: float | None = None
    pressure: float | None = None
    data: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class QHAWorkflowState:
    """State accumulated during a QHA pressure-temperature workflow.

    Parameters
    ----------
    total_points : int, optional
        Number of pressure-temperature points requested by the calculation.
    processed_points : int, optional
        Number of pressure-temperature points already processed.
    successful_points : int, optional
        Number of points that produced a valid result.
    failed_points : int, optional
        Number of points that failed during local analysis.
    warning_points : int, optional
        Number of points that produced diagnostic warnings.
    consecutive_failures : int, optional
        Number of consecutive failed pressure-temperature points.
    stopped : bool, optional
        Whether the workflow requested an early stop.
    stop_message : str, optional
        Explanation associated with an early stop.
    metadata : dict, optional
        Additional workflow details.
    """

    total_points: int = 0
    processed_points: int = 0
    successful_points: int = 0
    failed_points: int = 0
    warning_points: int = 0
    consecutive_failures: int = 0
    stopped: bool = False
    stop_message: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def progress(self) -> float:
        """Return the completed fraction of the requested workflow.

        Returns
        -------
        float
            Fractional progress between 0 and 1.
        """
        if self.total_points <= 0:
            return 0.0
        return min(1.0, self.processed_points / self.total_points)

    def as_dict(self) -> dict[str, Any]:
        """Return workflow counters as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary containing workflow counters and stop information.
        """
        return {
            "total_points": self.total_points,
            "processed_points": self.processed_points,
            "successful_points": self.successful_points,
            "failed_points": self.failed_points,
            "warning_points": self.warning_points,
            "consecutive_failures": self.consecutive_failures,
            "stopped": self.stopped,
            "stop_message": self.stop_message,
            "progress": self.progress,
            "metadata": dict(self.metadata),
        }


def _notify(
    callback: QHAWorkflowCallback | None,
    kind: QHAWorkflowEventKind,
    message: str,
    *,
    state: QHAWorkflowState,
    temperature: float | None = None,
    pressure: float | None = None,
    data: dict[str, Any] | None = None,
) -> None:
    """Send a workflow event to an optional callback.

    Parameters
    ----------
    callback : callable or None
        Function receiving workflow events.
    kind : str
        Event kind.
    message : str
        Event message.
    state : QHAWorkflowState
        Current workflow state.
    temperature : float or None, optional
        Temperature associated with the event.
    pressure : float or None, optional
        Pressure associated with the event.
    data : dict or None, optional
        Structured event payload.
    """
    if callback is None:
        return
    payload = {
        "current": int(state.processed_points),
        "total": int(state.total_points),
        "label": "pressure-temperature grid",
        **dict(data or {}),
    }
    callback(
        QHAWorkflowEvent(
            kind=kind,
            message=message,
            progress=state.progress,
            temperature=temperature,
            pressure=pressure,
            data=payload,
        )
    )


def _register_failure(
    result: QHAResult,
    state: QHAWorkflowState,
    *,
    temperature: float,
    pressure: float,
    stage: str,
    message: str,
    diagnostics: dict[str, Any] | None,
    options: QHAOptions,
    callback: QHAWorkflowCallback | None = None,
    count_consecutive: bool = True,
) -> None:
    """Record a failed point and apply the selected failure policy.

    Parameters
    ----------
    result : QHAResult
        Result container to update.
    state : QHAWorkflowState
        Mutable workflow state.
    temperature : float
        Temperature associated with the failure.
    pressure : float
        Pressure associated with the failure.
    stage : str
        Workflow stage where the failure occurred.
    message : str
        Failure message.
    diagnostics : dict or None
        Additional failure details.
    options : QHAOptions
        Workflow options controlling failure policy.
    callback : callable or None, optional
        Optional workflow-event callback.
    count_consecutive : bool, optional
        If ``True``, count this point toward the consecutive-failure policy.
        Shared temperature-fit failures set this to ``False`` and update the
        counter once after the complete pressure row has been recorded.

    Raises
    ------
    RuntimeError
        If ``options.fit_failure_policy`` is ``"raise"``.
    """
    payload = dict(diagnostics or {})
    result.add_failed_point(
        QHAFailedPoint(
            temperature=float(temperature),
            pressure=float(pressure),
            stage=stage,
            message=message,
            diagnostics=payload,
        )
    )
    state.failed_points += 1
    if count_consecutive:
        state.consecutive_failures += 1

    _notify(
        callback,
        "point_failed",
        message,
        state=state,
        temperature=float(temperature),
        pressure=float(pressure),
        data={"stage": stage, "diagnostics": payload},
    )

    if count_consecutive:
        _apply_failure_policy(
            result,
            state,
            temperature=temperature,
            pressure=pressure,
            stage=stage,
            message=message,
            options=options,
            callback=callback,
            scope="point",
        )


def _apply_failure_policy(
    result: QHAResult,
    state: QHAWorkflowState,
    *,
    temperature: float,
    pressure: float | None,
    stage: str,
    message: str,
    options: QHAOptions,
    callback: QHAWorkflowCallback | None,
    scope: str,
) -> None:
    """Apply the configured failure policy to the current failure counter.

    Parameters
    ----------
    result : QHAResult
        Result container updated when the workflow stops.
    state : QHAWorkflowState
        Mutable workflow counters.
    temperature : float
        Temperature associated with the failure.
    pressure : float or None
        Pressure associated with a point failure, or ``None`` for a shared
        temperature-fit failure.
    stage : str
        Failure stage.
    message : str
        Human-readable failure message.
    options : QHAOptions
        Failure-policy options.
    callback : callable or None
        Optional workflow callback.
    scope : str
        Counter scope, either ``"point"`` or ``"temperature"``.

    Raises
    ------
    RuntimeError
        If the selected policy is ``"raise"``.
    """
    if options.fit_failure_policy == "raise":
        raise RuntimeError(message)
    if (
        options.fit_failure_policy != "stop"
        or state.consecutive_failures < options.max_consecutive_failures
    ):
        return

    qualifier = "failed temperature fits" if scope == "temperature" else "failed fits"
    state.stopped = True
    state.stop_message = (
        f"stopped after {state.consecutive_failures} consecutive {qualifier}"
    )
    result.completed = False
    result.metadata.setdefault("stop", {})
    stop_metadata: dict[str, Any] = {
        "message": state.stop_message,
        "temperature": float(temperature),
        "max_consecutive_failures": int(options.max_consecutive_failures),
        "failure_scope": scope,
    }
    if pressure is not None:
        stop_metadata["pressure"] = float(pressure)
    result.metadata["stop"].update(stop_metadata)
    _notify(
        callback,
        "stopped",
        state.stop_message,
        state=state,
        temperature=float(temperature),
        pressure=None if pressure is None else float(pressure),
        data={"stage": stage, "scope": scope},
    )


def run_volume_minimization_workflow(
    input_data: QHAInput,
    options: QHAOptions,
    *,
    free_energy: np.ndarray | None = None,
    local_free_energy_evaluator: Callable[[np.ndarray, float, int], np.ndarray]
    | None = None,
    callback: QHAWorkflowCallback | None = None,
) -> QHAResult:
    """Run volume minimization over the pressure-temperature grid.

    Parameters
    ----------
    input_data : QHAInput
        Normalized input data containing sampled volumes and static energies.
    options : QHAOptions
        Workflow options controlling pressure and temperature grids,
        minimization method, diagnostics and failure policies.
    free_energy : ndarray or None, optional
        Free-energy data sampled on the input volume grid. If omitted, the
        static energies stored in ``input_data`` are used at all temperatures.
    local_free_energy_evaluator : callable or None, optional
        Function receiving local volumes, temperature and temperature index.
        It is used by the polynomial local-grid derivative method to rebuild
        Helmholtz free energies around each equilibrium volume.
    callback : callable or None, optional
        Function receiving structured workflow events.

    Returns
    -------
    QHAResult
        Result containing valid points, failed points, fit records, workflow
        counters and partial data when an early stop occurs.

    Raises
    ------
    ValueError
        If input arrays or options are inconsistent.
    RuntimeError
        If the selected failure policy requests exceptions on local failures.
    """
    result = initialize_result(input_data, options)
    volumes = np.asarray(input_data.volume, dtype=np.float64)
    temperatures = np.asarray(result.temperature, dtype=np.float64)
    pressures = np.asarray(result.pressure, dtype=np.float64)
    source_energy = input_data.energy if free_energy is None else free_energy
    if source_energy is None:
        raise ValueError("QHA workflow requires static or supplied free-energy data")
    energy_grid = prepare_free_energy_grid(
        source_energy,
        ntemperatures=temperatures.size,
        nvolumes=volumes.size,
    )

    state = QHAWorkflowState(total_points=int(temperatures.size * pressures.size))
    result.metadata["workflow"] = state.as_dict()
    _notify(callback, "started", "QHA workflow started", state=state)

    if options.minimization == "eos":
        result.metadata["eos_workflow"] = {
            "eos": options.eos,
            "fit_count": 0,
            "state_count": 0,
            "fit_scope": "temperature",
            "pressure_evaluation": "single_fitted_eos",
            "uncertainty_method": (
                options.uncertainty_method if options.estimate_uncertainties else "none"
            ),
        }
    else:
        result.metadata["polynomial_workflow"] = {
            "fit_count": 0,
            "state_count": 0,
            "fit_scope": "temperature",
            "pressure_evaluation": "single_fitted_polynomial",
            "derivative_method": options.polynomial_derivative_method,
            "local_grid_points": int(options.polynomial_grid_points),
            "local_grid_separation_percent": float(options.polynomial_grid_separation),
            "local_free_energy_source": (
                "not_used"
                if options.polynomial_derivative_method == "analytic"
                else (
                    "external_evaluator"
                    if local_free_energy_evaluator is not None
                    else "fitted_polynomial"
                )
            ),
        }

    for it, temperature in enumerate(temperatures):
        fitted_eos = None
        eos_fit_rejected = False
        eos_uncertainty_method = "none"
        fitted_polynomial = None
        polynomial_fit_rejected = False
        shared_failure: VolumeMinimumResult | None = None

        if options.minimization == "eos":
            fitted_eos = fit_eos_free_energy_model(
                volumes,
                energy_grid[it],
                eos=options.eos,
                energy_unit=options.energy_unit,
                volume_unit=options.volume_unit,
                pressure_unit=options.pressure_unit,
            )
            result.metadata["eos_workflow"]["fit_count"] += 1

            if options.estimate_uncertainties:
                eos_uncertainty_method = options.uncertainty_method
                if eos_uncertainty_method == "bootstrap":
                    fitted_eos.warnings.append(
                        "bootstrap uncertainty propagation is not available for "
                        "EOS pressure states; uncertainties were not calculated"
                    )
                    eos_uncertainty_method = "none"
                elif (
                    eos_uncertainty_method != "none"
                    and fitted_eos.model is not None
                    and fitted_eos.model.covariance is None
                ):
                    eos_uncertainty_method = "none"

            if options.store_fit_diagnostics:
                record = QHAFitRecord(
                    quantity="F",
                    method="eos",
                    temperature=float(temperature),
                    pressure=None,
                    fit=fitted_eos.fit,
                    success=fitted_eos.success,
                    message=fitted_eos.message,
                    metadata=fitted_eos.as_dict(),
                )
                result.add_fit_record(record)
                _notify(
                    callback,
                    "fit_record",
                    record.message or "EOS fit diagnostics recorded",
                    state=state,
                    temperature=float(temperature),
                    pressure=None,
                    data=record.as_dict(),
                )

            if fitted_eos.warnings:
                state.warning_points += 1
                _notify(
                    callback,
                    "warning",
                    "the temperature EOS fit produced diagnostic warnings",
                    state=state,
                    temperature=float(temperature),
                    pressure=None,
                    data={
                        "warnings": list(fitted_eos.warnings),
                        "eos_fit": fitted_eos.as_dict(),
                    },
                )

            eos_fit_rejected = (
                fitted_eos.success
                and fitted_eos.fit.quality is FitQuality.POOR
                and options.fit_quality_policy == "stop"
            )
            if not fitted_eos.success:
                shared_failure = VolumeMinimumResult.failed(
                    "the free-energy EOS fit failed",
                    method="eos",
                    fit=fitted_eos.fit,
                    metadata={
                        "failure_stage": "eos_fit",
                        "eos_fit": fitted_eos.as_dict(),
                    },
                )
            elif eos_fit_rejected:
                shared_failure = VolumeMinimumResult.failed(
                    "the temperature EOS fit did not satisfy the selected quality policy",
                    method="eos",
                    fit=fitted_eos.fit,
                    metadata={
                        "failure_stage": "fit_quality",
                        "eos_fit": fitted_eos.as_dict(),
                    },
                )
        else:
            fitted_polynomial = fit_polynomial_free_energy_model(
                volumes,
                energy_grid[it],
                degree=options.free_energy_degree,
            )
            result.metadata["polynomial_workflow"]["fit_count"] += 1

            if options.store_fit_diagnostics:
                record = QHAFitRecord(
                    quantity="F",
                    method="polynomial",
                    temperature=float(temperature),
                    pressure=None,
                    fit=fitted_polynomial.fit,
                    success=fitted_polynomial.success,
                    message=fitted_polynomial.message,
                    metadata=fitted_polynomial.as_dict(),
                )
                result.add_fit_record(record)
                _notify(
                    callback,
                    "fit_record",
                    record.message or "polynomial fit diagnostics recorded",
                    state=state,
                    temperature=float(temperature),
                    pressure=None,
                    data=record.as_dict(),
                )

            if fitted_polynomial.warnings:
                state.warning_points += 1
                _notify(
                    callback,
                    "warning",
                    "the temperature polynomial fit produced diagnostic warnings",
                    state=state,
                    temperature=float(temperature),
                    pressure=None,
                    data={
                        "warnings": list(fitted_polynomial.warnings),
                        "polynomial_fit": fitted_polynomial.as_dict(),
                    },
                )

            polynomial_fit_rejected = (
                fitted_polynomial.success
                and fitted_polynomial.fit.quality is FitQuality.POOR
                and options.fit_quality_policy == "stop"
            )
            if not fitted_polynomial.success:
                shared_failure = VolumeMinimumResult.failed(
                    "the free-energy polynomial fit failed at this temperature",
                    method="polynomial",
                    fit=fitted_polynomial.fit,
                    metadata={
                        "failure_stage": "polynomial_fit",
                        "polynomial_fit": fitted_polynomial.as_dict(),
                    },
                )
            elif polynomial_fit_rejected:
                shared_failure = VolumeMinimumResult.failed(
                    "the temperature polynomial fit did not satisfy the selected quality policy",
                    method="polynomial",
                    fit=fitted_polynomial.fit,
                    metadata={
                        "failure_stage": "fit_quality",
                        "polynomial_fit": fitted_polynomial.as_dict(),
                    },
                )

        if shared_failure is not None:
            stage = str(shared_failure.metadata.get("failure_stage", "temperature_fit"))
            for ip, pressure in enumerate(pressures):
                _notify(
                    callback,
                    "point_started",
                    "processing pressure-temperature point",
                    state=state,
                    temperature=float(temperature),
                    pressure=float(pressure),
                    data={"temperature_index": it, "pressure_index": ip},
                )
                state.processed_points += 1
                _register_failure(
                    result,
                    state,
                    temperature=float(temperature),
                    pressure=float(pressure),
                    stage=stage,
                    message=shared_failure.message,
                    diagnostics=shared_failure.as_dict(),
                    options=options,
                    callback=callback,
                    count_consecutive=False,
                )

            state.consecutive_failures += 1
            _apply_failure_policy(
                result,
                state,
                temperature=float(temperature),
                pressure=None,
                stage=stage,
                message=shared_failure.message,
                options=options,
                callback=callback,
                scope="temperature",
            )
            result.metadata["workflow"] = state.as_dict()
            if state.stopped:
                return result
            continue

        for ip, pressure in enumerate(pressures):
            _notify(
                callback,
                "point_started",
                "processing pressure-temperature point",
                state=state,
                temperature=float(temperature),
                pressure=float(pressure),
                data={"temperature_index": it, "pressure_index": ip},
            )

            if options.minimization == "eos":
                if fitted_eos is None or not fitted_eos.success:
                    minimum = VolumeMinimumResult.failed(
                        "the free-energy EOS fit failed",
                        method="eos",
                        fit=None if fitted_eos is None else fitted_eos.fit,
                        metadata={
                            "failure_stage": "eos_fit",
                            "eos_fit": (
                                None if fitted_eos is None else fitted_eos.as_dict()
                            ),
                        },
                    )
                elif eos_fit_rejected:
                    minimum = VolumeMinimumResult.failed(
                        "the temperature EOS fit did not satisfy the selected quality policy",
                        method="eos",
                        fit=fitted_eos.fit,
                        metadata={
                            "failure_stage": "fit_quality",
                            "eos_fit": fitted_eos.as_dict(),
                        },
                    )
                else:
                    random_state = options.uncertainty_seed
                    if random_state is not None:
                        random_state = int(random_state) + it
                    minimum = evaluate_fitted_eos_at_pressure(
                        fitted_eos,
                        float(pressure),
                        uncertainty_method=eos_uncertainty_method,
                        relative_step=options.uncertainty_relative_step,
                        confidence_level=options.uncertainty_confidence_level,
                        monte_carlo_samples=options.uncertainty_samples,
                        random_state=random_state,
                        minimum_success_fraction=(
                            options.uncertainty_minimum_success_fraction
                        ),
                    )
                    result.metadata["eos_workflow"]["state_count"] += 1
            else:
                if fitted_polynomial is None or not fitted_polynomial.success:
                    minimum = VolumeMinimumResult.failed(
                        "the free-energy polynomial fit failed at this temperature",
                        method="polynomial",
                        fit=(
                            None if fitted_polynomial is None else fitted_polynomial.fit
                        ),
                        metadata={
                            "failure_stage": "polynomial_fit",
                            "polynomial_fit": (
                                None
                                if fitted_polynomial is None
                                else fitted_polynomial.as_dict()
                            ),
                        },
                    )
                elif polynomial_fit_rejected:
                    minimum = VolumeMinimumResult.failed(
                        "the temperature polynomial fit did not satisfy the selected quality policy",
                        method="polynomial",
                        fit=fitted_polynomial.fit,
                        metadata={
                            "failure_stage": "fit_quality",
                            "polynomial_fit": fitted_polynomial.as_dict(),
                        },
                    )
                else:
                    local_evaluator = None
                    if local_free_energy_evaluator is not None:

                        def local_evaluator(local_volume):
                            return local_free_energy_evaluator(
                                np.asarray(local_volume, dtype=np.float64),
                                float(temperature),
                                int(it),
                            )

                    minimum = evaluate_fitted_polynomial_at_pressure(
                        fitted_polynomial,
                        pressure_energy_density(float(pressure), options),
                        derivative_method=options.polynomial_derivative_method,
                        local_free_energy=local_evaluator,
                        local_grid_points=options.polynomial_grid_points,
                        local_grid_separation=options.polynomial_grid_separation,
                        local_degree=options.energy_degree,
                    )
                    if minimum.bulk_modulus is not None:
                        minimum.bulk_modulus = float(
                            energy_to_pressure(
                                minimum.bulk_modulus,
                                options.energy_unit,
                                options.volume_unit,
                                options.pressure_unit,
                            )
                        )
                    result.metadata["polynomial_workflow"]["state_count"] += 1

            state.processed_points += 1

            if not minimum.success:
                _register_failure(
                    result,
                    state,
                    temperature=float(temperature),
                    pressure=float(pressure),
                    stage=str(
                        minimum.metadata.get("failure_stage", "volume_minimization")
                    ),
                    message=minimum.message,
                    diagnostics=minimum.as_dict(),
                    options=options,
                    callback=callback,
                )
                result.metadata["workflow"] = state.as_dict()
                if state.stopped:
                    return result
                continue

            if (
                minimum.range_status is VolumeRangeStatus.OUTSIDE
                and options.extrapolation_policy == "fail"
            ):
                _register_failure(
                    result,
                    state,
                    temperature=float(temperature),
                    pressure=float(pressure),
                    stage="volume_extrapolation",
                    message="the equilibrium volume is outside the sampled interval",
                    diagnostics=minimum.as_dict(),
                    options=options,
                    callback=callback,
                )
                result.metadata["workflow"] = state.as_dict()
                if state.stopped:
                    return result
                continue

            if minimum.warnings:
                state.warning_points += 1
                _notify(
                    callback,
                    "warning",
                    "the pressure-temperature state produced diagnostic warnings",
                    state=state,
                    temperature=float(temperature),
                    pressure=float(pressure),
                    data={
                        "warnings": list(minimum.warnings),
                        "minimum": minimum.as_dict(),
                    },
                )

            if (
                options.minimization != "eos"
                and minimum.warnings
                and options.fit_quality_policy == "stop"
            ):
                _register_failure(
                    result,
                    state,
                    temperature=float(temperature),
                    pressure=float(pressure),
                    stage="fit_quality",
                    message="the local fit produced diagnostic warnings",
                    diagnostics=minimum.as_dict(),
                    options=options,
                    callback=callback,
                )
                result.metadata["workflow"] = state.as_dict()
                if state.stopped:
                    return result
                continue

            if (
                options.minimization == "poly"
                and options.polynomial_derivative_method == "local_grid"
                and options.store_fit_diagnostics
                and options.debug
            ):
                local_record = QHAFitRecord(
                    quantity="KT/Kp",
                    method="polynomial_local_grid",
                    temperature=float(temperature),
                    pressure=float(pressure),
                    fit=None,
                    success=True,
                    message="local polynomial thermoelastic derivatives evaluated",
                    metadata=minimum.as_dict(),
                )
                result.add_fit_record(local_record)
                _notify(
                    callback,
                    "fit_record",
                    local_record.message,
                    state=state,
                    temperature=float(temperature),
                    pressure=float(pressure),
                    data=local_record.as_dict(),
                )

            apply_minimum_to_result(
                result,
                minimum,
                temperature_index=it,
                pressure_index=ip,
            )
            state.successful_points += 1
            state.consecutive_failures = 0
            result.metadata["workflow"] = state.as_dict()
            _notify(
                callback,
                "point_completed",
                "pressure-temperature point completed",
                state=state,
                temperature=float(temperature),
                pressure=float(pressure),
                data={
                    "current": int(state.processed_points),
                    "total": int(state.total_points),
                    "label": "pressure-temperature grid",
                    "minimum": minimum.as_dict(),
                    "volume": minimum.volume,
                    "bulk_modulus": minimum.bulk_modulus,
                    "bulk_modulus_derivative": minimum.bulk_modulus_derivative,
                    "sigma_volume": minimum.sigma_volume,
                    "sigma_bulk_modulus": minimum.sigma_bulk_modulus,
                    "sigma_bulk_modulus_derivative": (
                        minimum.sigma_bulk_modulus_derivative
                    ),
                    "range_status": (
                        None
                        if minimum.range_status is None
                        else minimum.range_status.value
                    ),
                },
            )

    result.metadata["workflow"] = state.as_dict()
    _notify(callback, "completed", "QHA workflow completed", state=state)
    return result
