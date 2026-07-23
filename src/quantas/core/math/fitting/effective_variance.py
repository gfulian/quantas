# -*- coding: utf-8 -*-

r"""Iteratively reweighted fitting with the effective-variance method.

For observations with response :math:`y_i` and one or more explanatory
coordinates :math:`x_{j,i}`, the method uses

.. math::

    \sigma_{\mathrm{eff},i}^2 = \sigma_{y,i}^2 +
    \sum_j\left(
    \frac{\partial f}{\partial x_j}\sigma_{x_j,i}
    \right)^2.

The familiar one-coordinate expression is recovered when the sum contains
only one term.

Because the model derivative normally depends on fitted parameters, Quantas
updates the effective standard uncertainties between nonlinear weighted
least-squares cycles until both parameters and weights converge.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from .diagnostics import standardized_residuals
from .least_squares import LeastSquaresFitter
from .model import BaseFitModel, MappedFitModel
from .observations import FitObservations
from .options import (
    EffectiveVarianceOptions,
    FitMethod,
    FitOptions,
    WLSOptions,
)
from .parameters import ParameterMap
from .result import (
    FitDiagnostics,
    FitQuality,
    FitResult,
    FitStatus,
    resolve_mapped_fit_result,
)

from .solver_debug import (
    array_summary,
    detect_period_two_oscillation,
    problem_summary,
    solver_debug_enabled,
    termination_metadata,
)


@dataclass(slots=True)
class _PreparedEffectiveVariance:
    """Validated inputs for iterative effective-variance fitting."""

    selected: FitObservations
    sigma_x: np.ndarray
    sigma_y: np.ndarray
    mapped_model: MappedFitModel
    parameter_map: ParameterMap
    initial: np.ndarray
    bounds: tuple[np.ndarray, np.ndarray]
    settings: EffectiveVarianceOptions
    metadata: dict[str, Any]


@dataclass(slots=True)
class _EffectiveVarianceExecution:
    """Converged reduced fit and final uncertainty components."""

    result: FitResult
    sigma_effective: np.ndarray
    derivative: np.ndarray
    projected_sigma_x: np.ndarray
    history: list[dict[str, Any]]


class EffectiveVarianceFitter:
    """Fit models by iterative effective-variance reweighting.

    Parameters
    ----------
    least_squares : LeastSquaresFitter or None, optional
        Weighted least-squares backend used for every reweighting cycle.
        Supplying a backend is primarily useful for testing or controlled
        dependency injection.
    """

    def __init__(self, least_squares: LeastSquaresFitter | None = None) -> None:
        self._least_squares = least_squares or LeastSquaresFitter()

    def fit_problem(
        self,
        model: BaseFitModel,
        observations: FitObservations,
        parameters: ParameterMap,
        options: FitOptions,
    ) -> FitResult:
        """Fit one complete effective-variance problem.

        Parameters
        ----------
        model : BaseFitModel
            Full physical model. An analytical ``derivative_x`` implementation
            is preferred; the base contract supplies a numerical fallback.
        observations : FitObservations
            Coordinates, observations, standard uncertainties, and mask. At
            least one of ``sigma_x`` or ``sigma_y`` must be available. Matrix
            explanatory coordinates contribute one projected variance term
            per coordinate.
        parameters : ParameterMap
            Complete FREE/FIXED/IMPLIED/DERIVED parameter mapping.
        options : FitOptions
            Effective-variance and optimizer controls.

        Returns
        -------
        FitResult
            Complete resolved result with final effective uncertainties and
            iteration history in :attr:`FitDiagnostics.metadata`.
        """
        try:
            prepared = _prepare_problem(model, observations, parameters, options)
        except (TypeError, ValueError) as exc:
            return _invalid_result(model, observations, parameters, options, exc)

        execution = _run_reweighting(prepared, self._least_squares)
        if isinstance(execution, FitResult):
            return execution
        reduced = _effective_variance_result(prepared, execution)
        return resolve_mapped_fit_result(reduced, parameters)


def _prepare_problem(
    model: BaseFitModel,
    observations: FitObservations,
    parameters: ParameterMap,
    options: FitOptions,
) -> _PreparedEffectiveVariance:
    """Validate and normalize an effective-variance problem."""
    settings = _effective_variance_options(options)
    selected = observations.selected()
    sigma_x, sigma_y = _uncertainty_components(selected)
    mapped_model = MappedFitModel(model, parameters)
    initial = parameters.initial_free_values()
    if selected.size <= initial.size:
        raise ValueError(
            "effective-variance fitting requires more observations than free parameters"
        )
    return _PreparedEffectiveVariance(
        selected=selected,
        sigma_x=sigma_x,
        sigma_y=sigma_y,
        mapped_model=mapped_model,
        parameter_map=parameters,
        initial=initial,
        bounds=parameters.free_bounds(),
        settings=settings,
        metadata={**model.metadata(), **options.metadata},
    )


def _run_reweighting(
    prepared: _PreparedEffectiveVariance,
    least_squares: LeastSquaresFitter,
) -> _EffectiveVarianceExecution | FitResult:
    """Execute weighted least-squares cycles until weights converge."""
    try:
        sigma_effective, derivative, projected = _effective_sigma(
            prepared.mapped_model,
            prepared.selected.x,
            prepared.initial,
            prepared.sigma_x,
            prepared.sigma_y,
        )
    except ValueError as exc:
        return _prepared_failure(
            prepared,
            str(exc),
            status=FitStatus.INVALID_INPUT,
            error=exc,
            termination_category="invalid_input",
        )

    current = prepared.initial
    history: list[dict[str, Any]] = []
    last_inner: FitResult | None = None
    max_iterations = int(prepared.settings.max_iterations or 25)
    for iteration in range(1, max_iterations + 1):
        inner = _weighted_cycle(
            prepared,
            least_squares,
            current,
            sigma_effective,
        )
        if not inner.success or inner.parameters is None:
            return _failed_cycle_result(
                prepared, inner, iteration, history, last_inner, sigma_effective
            )
        last_inner = inner
        try:
            updated_sigma, updated_derivative, updated_projected = _effective_sigma(
                prepared.mapped_model,
                prepared.selected.x,
                inner.parameters,
                prepared.sigma_x,
                prepared.sigma_y,
            )
        except ValueError as exc:
            return _prepared_failure(
                prepared,
                f"effective-variance update failed: {exc}",
                history=history,
                last_result=inner,
                last_sigma=sigma_effective,
                error=exc,
                termination_category="invalid_update",
            )
        parameter_change, sigma_change = _record_iteration(
            prepared.settings,
            history,
            iteration,
            current,
            inner.parameters,
            sigma_effective,
            updated_sigma,
            inner,
        )
        if parameter_change <= 1.0 and sigma_change <= 1.0:
            return _EffectiveVarianceExecution(
                result=inner,
                sigma_effective=updated_sigma,
                derivative=updated_derivative,
                projected_sigma_x=updated_projected,
                history=history,
            )
        current = inner.parameters.copy()
        sigma_effective = updated_sigma
        derivative = updated_derivative
        projected = updated_projected

    del derivative, projected
    return _prepared_failure(
        prepared,
        "effective-variance weights did not converge within "
        f"{max_iterations} iterations",
        history=history,
        last_result=last_inner,
        last_sigma=sigma_effective,
        termination_category="iteration_limit",
    )


def _weighted_cycle(
    prepared: _PreparedEffectiveVariance,
    least_squares: LeastSquaresFitter,
    initial: np.ndarray,
    sigma_effective: np.ndarray,
) -> FitResult:
    """Run one pressure-like weighted nonlinear least-squares cycle."""
    return least_squares.fit(
        prepared.mapped_model,
        prepared.selected.x,
        prepared.selected.y,
        initial_parameters=initial,
        bounds=prepared.bounds,
        sigma=sigma_effective,
        options=_inner_options(prepared.settings),
    )


def _record_iteration(
    settings: EffectiveVarianceOptions,
    history: list[dict[str, Any]],
    iteration: int,
    previous_parameters: np.ndarray,
    current_parameters: np.ndarray,
    previous_sigma: np.ndarray,
    current_sigma: np.ndarray,
    result: FitResult,
) -> tuple[float, float]:
    """Append one convergence record and return normalized changes."""
    parameter_change = _scaled_change(
        previous_parameters,
        current_parameters,
        settings.parameter_atol,
        settings.parameter_rtol,
    )
    sigma_change = _scaled_change(
        previous_sigma,
        current_sigma,
        settings.sigma_atol,
        settings.sigma_rtol,
    )
    record: dict[str, Any] = {
        "iteration": iteration,
        "parameter_change": parameter_change,
        "sigma_change": sigma_change,
        "chi_square": _chi_square(result.residuals, previous_sigma),
        "rmse": result.rmse,
        "maximum_absolute_residual": result.max_abs_error,
        "parameters": np.asarray(current_parameters, dtype=np.float64).tolist(),
        "sigma_minimum": float(np.min(current_sigma)),
        "sigma_maximum": float(np.max(current_sigma)),
        "sigma_dynamic_range": float(np.max(current_sigma) / np.min(current_sigma)),
    }
    if result.diagnostics is not None:
        record["inner_evaluations"] = result.diagnostics.n_evaluations
        record["inner_stop_reason"] = result.diagnostics.stop_reason
    history.append(record)
    return parameter_change, sigma_change


def _effective_variance_options(options: FitOptions) -> EffectiveVarianceOptions:
    """Normalize common controls for effective-variance regression."""
    if isinstance(options, EffectiveVarianceOptions):
        return options
    return EffectiveVarianceOptions(
        covariance_scaling=options.covariance_policy,
        max_iterations=options.max_iterations,
        ftol=options.ftol,
        xtol=options.xtol,
        gtol=options.gtol,
        metadata=options.metadata,
    )


def _inner_options(settings: EffectiveVarianceOptions) -> WLSOptions:
    """Return typed controls for one inner weighted least-squares cycle."""
    return WLSOptions(
        covariance_scaling=settings.covariance_policy,
        max_iterations=settings.inner_max_iterations,
        ftol=settings.ftol,
        xtol=settings.xtol,
        gtol=settings.gtol,
        metadata=settings.metadata,
    )


def _uncertainty_components(
    observations: FitObservations,
) -> tuple[np.ndarray, np.ndarray]:
    """Return non-negative x and y uncertainty components."""
    if observations.sigma_x is None and observations.sigma_y is None:
        raise ValueError("effective variance requires sigma_x, sigma_y, or both")
    sigma_x = (
        np.zeros(observations.x.shape, dtype=np.float64)
        if observations.sigma_x is None
        else np.asarray(observations.sigma_x, dtype=np.float64).copy()
    )
    sigma_y = (
        np.zeros(observations.y.shape, dtype=np.float64)
        if observations.sigma_y is None
        else np.asarray(observations.sigma_y, dtype=np.float64).copy()
    )
    if np.any(sigma_x < 0.0) or np.any(sigma_y < 0.0):
        raise ValueError("effective-variance uncertainties cannot be negative")
    return sigma_x, sigma_y


def _effective_sigma(
    model: BaseFitModel,
    x: np.ndarray,
    parameters: np.ndarray,
    sigma_x: np.ndarray,
    sigma_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate effective sigma and its projected x component."""
    derivative = np.asarray(model.derivative_x(x, parameters), dtype=np.float64)
    if derivative.shape != x.shape or not np.all(np.isfinite(derivative)):
        raise ValueError("model derivative must be finite and match x")
    projected = derivative * sigma_x
    projected_variance = (
        projected**2 if projected.ndim == 1 else np.sum(projected**2, axis=0)
    )
    variance = sigma_y**2 + projected_variance
    if not np.all(np.isfinite(variance)) or np.any(variance <= 0.0):
        raise ValueError(
            "effective variance must be finite and strictly positive at every point"
        )
    return np.sqrt(variance), derivative, projected


def _scaled_change(
    previous: np.ndarray,
    current: np.ndarray,
    atol: float,
    rtol: float,
) -> float:
    """Return the maximum change normalized by absolute/relative tolerance."""
    scale = atol + rtol * np.maximum(np.abs(previous), np.abs(current))
    scale = np.maximum(scale, np.finfo(np.float64).tiny)
    return float(np.max(np.abs(current - previous) / scale))


def _chi_square(residuals: np.ndarray | None, sigma: np.ndarray) -> float | None:
    """Return chi-square for residuals and effective uncertainties."""
    if residuals is None:
        return None
    return float(np.sum((np.asarray(residuals, dtype=np.float64) / sigma) ** 2))


def _effective_variance_result(
    prepared: _PreparedEffectiveVariance,
    execution: _EffectiveVarianceExecution,
) -> FitResult:
    """Relabel a converged WLS cycle as effective variance."""
    result = execution.result
    residuals = np.asarray(result.residuals, dtype=np.float64)
    chi_square = float(np.sum((residuals / execution.sigma_effective) ** 2))
    reduced_chi_square = chi_square / result.dof if result.dof > 0 else None
    previous = result.diagnostics or FitDiagnostics()
    diagnostics = _effective_diagnostics(
        prepared,
        execution,
        previous,
        residuals,
        chi_square,
        reduced_chi_square,
    )
    return replace(
        result,
        method=FitMethod.EFFECTIVE_VARIANCE,
        diagnostics=diagnostics,
        metadata={
            **result.metadata,
            "fit_method": FitMethod.EFFECTIVE_VARIANCE.value,
        },
    )


def _effective_diagnostics(
    prepared: _PreparedEffectiveVariance,
    execution: _EffectiveVarianceExecution,
    previous: FitDiagnostics,
    residuals: np.ndarray,
    chi_square: float,
    reduced_chi_square: float | None,
) -> FitDiagnostics:
    """Build method-neutral diagnostics for effective variance."""
    settings = prepared.settings
    projected = execution.projected_sigma_x
    return FitDiagnostics(
        objective=(
            "sum(((y_observed-y_calculated)/sigma_effective)^2), "
            "sigma_effective^2=sigma_y^2+sum_j((df/dx_j*sigma_x_j)^2)"
        ),
        weighted=True,
        chi_square=chi_square,
        reduced_chi_square=reduced_chi_square,
        correlation=previous.correlation,
        standardized_residuals=standardized_residuals(
            residuals, execution.sigma_effective
        ),
        jacobian_rank=previous.jacobian_rank,
        condition_number=previous.condition_number,
        parameters_at_bounds=previous.parameters_at_bounds,
        n_iterations=len(execution.history),
        n_evaluations=previous.n_evaluations,
        stop_reason="effective parameters and weights converged",
        warnings=previous.warnings,
        metadata={
            **previous.metadata,
            "covariance_scaling": settings.covariance_policy.value,
            "covariance_scale_factor": previous.metadata.get(
                "covariance_scale_factor", 1.0
            ),
            "effective_sigma": execution.sigma_effective.tolist(),
            "sigma_x": prepared.sigma_x.tolist(),
            "sigma_y": prepared.sigma_y.tolist(),
            "derivative_x": execution.derivative.tolist(),
            "projected_sigma_x": projected.tolist(),
            "variance_x_component": (projected**2).tolist(),
            "variance_y_component": (prepared.sigma_y**2).tolist(),
            **problem_summary(
                prepared.selected.x,
                prepared.selected.y,
                prepared.initial,
                prepared.bounds,
                execution.sigma_effective,
            ),
            "sigma_x_summary": array_summary(prepared.sigma_x, positive_ratio=True),
            "sigma_y_summary": array_summary(prepared.sigma_y, positive_ratio=True),
            "free_parameter_names": list(prepared.parameter_map.free_names),
            "detailed_trace_requested": solver_debug_enabled(settings.metadata),
            **termination_metadata(
                category="converged",
                backend="quantas.effective_variance+scipy.optimize.curve_fit",
            ),
            "reweighting_history": execution.history,
            "reweighting_iterations": len(execution.history),
            "period_two_oscillation_detected": detect_period_two_oscillation(
                execution.history
            ),
            "parameter_rtol": settings.parameter_rtol,
            "parameter_atol": settings.parameter_atol,
            "sigma_rtol": settings.sigma_rtol,
            "sigma_atol": settings.sigma_atol,
        },
    )


def _invalid_result(
    model: BaseFitModel,
    observations: FitObservations,
    parameters: ParameterMap,
    options: FitOptions,
    error: Exception,
) -> FitResult:
    """Return an invalid-input result with complete mapping labels."""
    metadata = {
        **model.metadata(),
        **options.metadata,
        **termination_metadata(
            category="invalid_input",
            backend="quantas.effective_variance",
            exception=error,
        ),
    }
    diagnostics = FitDiagnostics(
        objective="iteratively reweighted effective-variance objective",
        weighted=True,
        stop_reason=str(error),
        metadata=metadata,
    )
    return FitResult.failed(
        str(error),
        status=FitStatus.INVALID_INPUT,
        n_points=observations.selected_size,
        n_parameters=parameters.n_free,
        metadata=metadata,
        method=FitMethod.EFFECTIVE_VARIANCE,
        parameter_names=parameters.names,
        parameter_states=parameters.states,
        diagnostics=diagnostics,
    )


def _prepared_failure(
    prepared: _PreparedEffectiveVariance,
    message: str,
    *,
    status: FitStatus = FitStatus.FAILED,
    history: list[dict[str, Any]] | None = None,
    last_result: FitResult | None = None,
    last_sigma: np.ndarray | None = None,
    error: Exception | None = None,
    termination_category: str = "nonconvergence",
) -> FitResult:
    """Return a failure tied to one prepared problem with rich diagnostics."""
    records = [] if history is None else list(history)
    sigma = prepared.sigma_y if last_sigma is None else np.asarray(last_sigma)
    metadata: dict[str, Any] = {
        **prepared.metadata,
        **problem_summary(
            prepared.selected.x,
            prepared.selected.y,
            prepared.initial,
            prepared.bounds,
            sigma,
        ),
        "sigma_x_summary": array_summary(prepared.sigma_x, positive_ratio=True),
        "sigma_y_summary": array_summary(prepared.sigma_y, positive_ratio=True),
        "free_parameter_names": list(prepared.parameter_map.free_names),
        "detailed_trace_requested": solver_debug_enabled(prepared.settings.metadata),
        **termination_metadata(
            category=termination_category,
            backend="quantas.effective_variance+scipy.optimize.curve_fit",
            exception=error,
        ),
        "effective_variance_history": records,
        "effective_variance_converged": False,
        "period_two_oscillation_detected": detect_period_two_oscillation(records),
    }
    if records:
        metadata["last_valid_iteration"] = int(records[-1]["iteration"])
        metadata["last_parameter_change"] = records[-1].get("parameter_change")
        metadata["last_sigma_change"] = records[-1].get("sigma_change")
    if last_result is not None and last_result.parameters is not None:
        free = np.asarray(last_result.parameters, dtype=np.float64)
        metadata["last_valid_free_parameters"] = free.tolist()
        metadata["last_valid_parameters"] = prepared.parameter_map.expand(
            free
        ).values.tolist()
        metadata["last_valid_parameter_names"] = list(prepared.parameter_map.names)
        metadata["last_valid_rmse"] = last_result.rmse
        metadata["last_valid_maximum_absolute_residual"] = last_result.max_abs_error
        if last_result.diagnostics is not None:
            metadata["last_inner_stop_reason"] = last_result.diagnostics.stop_reason
            metadata["last_inner_evaluations"] = last_result.diagnostics.n_evaluations
            inner_metadata = last_result.diagnostics.metadata
            for key in (
                "first_evaluation",
                "last_evaluation",
                "evaluation_trace",
                "evaluation_trace_truncated",
                "recorded_model_evaluations",
            ):
                if key in inner_metadata:
                    metadata[key] = inner_metadata[key]
    total_evaluations = sum(int(item.get("inner_evaluations") or 0) for item in records)
    metadata["total_inner_evaluations"] = total_evaluations
    diagnostics = FitDiagnostics(
        objective=(
            "sum(((y_observed-y_calculated)/sigma_effective)^2), "
            "sigma_effective^2=sigma_y^2+sum_j((df/dx_j*sigma_x_j)^2)"
        ),
        weighted=True,
        chi_square=(None if not records else records[-1].get("chi_square")),
        n_iterations=len(records),
        n_evaluations=total_evaluations or None,
        stop_reason=message,
        metadata=metadata,
    )
    return FitResult.failed(
        message,
        status=status,
        n_points=prepared.selected.size,
        n_parameters=prepared.parameter_map.n_free,
        metadata=metadata,
        method=FitMethod.EFFECTIVE_VARIANCE,
        parameter_names=prepared.parameter_map.names,
        parameter_states=prepared.parameter_map.states,
        diagnostics=diagnostics,
    )


def _failed_cycle_result(
    prepared: _PreparedEffectiveVariance,
    result: FitResult,
    iteration: int,
    history: list[dict[str, Any]],
    last_result: FitResult | None,
    sigma_effective: np.ndarray,
) -> FitResult:
    """Return a stable failure when an inner WLS cycle fails."""
    message = f"effective-variance cycle {iteration} failed: {result.message}"
    failed_metadata = (
        {} if result.diagnostics is None else dict(result.diagnostics.metadata)
    )
    metadata = {
        **prepared.metadata,
        **problem_summary(
            prepared.selected.x,
            prepared.selected.y,
            prepared.initial,
            prepared.bounds,
            sigma_effective,
        ),
        **termination_metadata(
            category="inner_solver_failure",
            backend="quantas.effective_variance+scipy.optimize.curve_fit",
        ),
        "effective_variance_history": list(history),
        "effective_variance_converged": False,
        "free_parameter_names": list(prepared.parameter_map.free_names),
        "detailed_trace_requested": solver_debug_enabled(prepared.settings.metadata),
        "failed_iteration": iteration,
        "inner_failure": failed_metadata,
    }
    if last_result is not None and last_result.parameters is not None:
        free = np.asarray(last_result.parameters, dtype=np.float64)
        metadata["last_valid_free_parameters"] = free.tolist()
        metadata["last_valid_parameters"] = prepared.parameter_map.expand(
            free
        ).values.tolist()
        metadata["last_valid_parameter_names"] = list(prepared.parameter_map.names)
    diagnostics = FitDiagnostics(
        objective="iteratively reweighted effective-variance objective",
        weighted=True,
        n_iterations=len(history),
        n_evaluations=(
            None if result.diagnostics is None else result.diagnostics.n_evaluations
        ),
        stop_reason=message,
        metadata=metadata,
    )
    return FitResult(
        success=False,
        status=result.status,
        quality=FitQuality.FAILED,
        message=message,
        n_points=result.n_points,
        n_parameters=prepared.parameter_map.n_free,
        dof=result.dof,
        warnings=result.warnings,
        metadata=metadata,
        method=FitMethod.EFFECTIVE_VARIANCE,
        parameter_names=prepared.parameter_map.names,
        parameter_states=prepared.parameter_map.states,
        diagnostics=diagnostics,
    )


__all__ = ["EffectiveVarianceFitter"]
