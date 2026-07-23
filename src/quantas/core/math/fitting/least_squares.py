# -*- coding: utf-8 -*-

"""Nonlinear ordinary and weighted least squares based on SciPy."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any
import warnings

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from .diagnostics import (
    assess_fit_quality,
    covariance_correlation,
    covariance_errors,
    covariance_scale_factor,
    parameters_at_bounds,
    residual_metrics,
    standardized_residuals,
    validate_xy,
)
from .model import (
    ArrayLike,
    BaseFitModel,
    CallableFitModel,
    MappedFitModel,
    ModelFunction,
)
from .observations import FitObservations
from .options import (
    CovarianceScaling,
    FitMethod,
    FitOptions,
    LeastSquaresOptions,
    OLSOptions,
    WLSOptions,
)
from .parameters import ParameterMap, ParameterState
from .result import (
    FitDiagnostics,
    FitQuality,
    FitResult,
    FitStatus,
    resolve_mapped_fit_result,
)
from .solver_debug import (
    ModelEvaluationRecorder,
    problem_summary,
    solver_debug_enabled,
    termination_metadata,
)


@dataclass(slots=True)
class _PreparedLeastSquares:
    """Validated internal inputs for one least-squares execution."""

    x: np.ndarray
    y: np.ndarray
    initial: np.ndarray
    bounds: tuple[np.ndarray, np.ndarray]
    sigma: np.ndarray | None
    settings: LeastSquaresOptions
    metadata: dict[str, Any]
    parameter_names: tuple[str, ...]
    parameter_states: tuple[ParameterState, ...]

    @property
    def n_points(self) -> int:
        """Return the number of observations."""
        return int(self.y.size)

    @property
    def n_parameters(self) -> int:
        """Return the number of optimized parameters."""
        return int(self.initial.size)


@dataclass(slots=True)
class _LeastSquaresExecution:
    """Backend output and trace metadata for one SciPy execution."""

    parameters: np.ndarray
    covariance: np.ndarray
    warning_messages: list[str]
    status_code: int | None
    stop_reason: str
    n_evaluations: int | None
    recorder: ModelEvaluationRecorder


class _LeastSquaresExecutionError(RuntimeError):
    """SciPy execution failure retaining the completed evaluation trace."""

    def __init__(self, error: Exception, recorder: ModelEvaluationRecorder) -> None:
        super().__init__(str(error))
        self.error = error
        self.recorder = recorder


class _InvalidFitProblem(ValueError):
    """Input validation failure with available fit dimensions."""

    def __init__(
        self,
        message: str,
        *,
        n_points: int = 0,
        n_parameters: int = 0,
    ) -> None:
        super().__init__(message)
        self.n_points = n_points
        self.n_parameters = n_parameters


class LeastSquaresFitter:
    """Fit model objects with ordinary or weighted nonlinear least squares."""

    def fit_problem(
        self,
        model: BaseFitModel,
        observations: FitObservations,
        parameters: ParameterMap,
        options: FitOptions,
    ) -> FitResult:
        """Fit one complete backend-neutral numerical problem.

        Parameters
        ----------
        model : BaseFitModel
            Full model equation. Its parameter order must match the non-derived
            entries in ``parameters``.
        observations : FitObservations
            Coordinates, observations, uncertainties, and selection mask.
        parameters : ParameterMap
            Complete FREE/FIXED/IMPLIED/DERIVED parameter mapping.
        options : FitOptions
            General controls. This fitter accepts only explicit OLS or WLS.

        Returns
        -------
        FitResult
            Complete resolved parameters and method-neutral diagnostics.
        """
        try:
            selected = observations.selected()
            settings = _least_squares_options(options)
            sigma = (
                selected.require_positive_sigma("y")
                if settings.method is FitMethod.WLS
                else None
            )
            mapped_model = MappedFitModel(model, parameters)
        except (TypeError, ValueError) as exc:
            diagnostic_metadata = {
                **model.metadata(),
                **options.metadata,
                "free_parameter_names": list(parameters.free_names),
                "detailed_trace_requested": solver_debug_enabled(options.metadata),
                **termination_metadata(
                    category="invalid_input",
                    backend="scipy.optimize.curve_fit",
                    exception=exc,
                ),
            }
            diagnostics = FitDiagnostics(
                objective=(
                    "sum((observed-calculated)^2)"
                    if options.method is FitMethod.OLS
                    else "sum(((observed-calculated)/sigma_y)^2)"
                ),
                weighted=options.method is FitMethod.WLS,
                stop_reason=str(exc),
                metadata=diagnostic_metadata,
            )
            return FitResult.failed(
                str(exc),
                status=FitStatus.INVALID_INPUT,
                n_points=observations.selected_size,
                n_parameters=parameters.n_free,
                metadata=diagnostic_metadata,
                method=options.method,
                parameter_names=parameters.names,
                parameter_states=parameters.states,
                diagnostics=diagnostics,
            )

        reduced = self.fit(
            mapped_model,
            selected.x,
            selected.y,
            sigma=sigma,
            options=settings,
        )
        if not reduced.success:
            return reduced
        return resolve_mapped_fit_result(reduced, parameters)

    def fit(
        self,
        model: BaseFitModel,
        x: ArrayLike,
        y: ArrayLike,
        *,
        initial_parameters: ArrayLike | None = None,
        bounds: tuple[ArrayLike, ArrayLike] | None = None,
        sigma: ArrayLike | None = None,
        options: LeastSquaresOptions | None = None,
    ) -> FitResult:
        """Fit a model to observations.

        Parameters
        ----------
        model : BaseFitModel
            Model equation and parameter definition.
        x, y : array-like
            Coordinates and observations.
        initial_parameters : sequence of float, optional
            Initial estimate overriding ``model.initial_guess``.
        bounds : tuple of array-like, optional
            Bounds overriding ``model.bounds``.
        sigma : array-like or None, optional
            Dependent-variable standard uncertainties for explicit weighted
            least squares. Ordinary least squares rejects supplied sigma.
        options : LeastSquaresOptions, optional
            Optimizer and metadata settings.

        Returns
        -------
        FitResult
            Structured optimized parameters and diagnostics.
        """
        settings = options or OLSOptions()
        metadata = {**model.metadata(), **settings.metadata}
        parameter_names = tuple(model.parameter_names)
        parameter_states = tuple(ParameterState.FREE for _ in parameter_names)
        try:
            prepared = _prepare_problem(
                model,
                x,
                y,
                initial_parameters=initial_parameters,
                bounds=bounds,
                sigma=sigma,
                settings=settings,
                metadata=metadata,
                parameter_names=parameter_names,
                parameter_states=parameter_states,
            )
        except _InvalidFitProblem as exc:
            diagnostic_metadata = {
                **metadata,
                "free_parameter_names": list(parameter_names),
                "detailed_trace_requested": solver_debug_enabled(settings.metadata),
                **termination_metadata(
                    category="invalid_input",
                    backend="scipy.optimize.curve_fit",
                    exception=exc,
                ),
            }
            diagnostics = FitDiagnostics(
                objective=(
                    "sum((observed-calculated)^2)"
                    if settings.method is FitMethod.OLS
                    else "sum(((observed-calculated)/sigma_y)^2)"
                ),
                weighted=settings.method is FitMethod.WLS,
                stop_reason=str(exc),
                metadata=diagnostic_metadata,
            )
            return _failed_result(
                str(exc),
                settings,
                diagnostic_metadata,
                parameter_names,
                parameter_states,
                status=FitStatus.INVALID_INPUT,
                n_points=exc.n_points,
                n_parameters=exc.n_parameters,
                diagnostics=diagnostics,
            )

        try:
            execution = _execute_fit(model, prepared)
        except _LeastSquaresExecutionError as exc:
            diagnostics = _least_squares_failure_diagnostics(
                prepared,
                exc.recorder,
                exc.error,
            )
            return _failed_result(
                str(exc.error),
                settings,
                {**metadata, **diagnostics.metadata},
                parameter_names,
                parameter_states,
                n_points=prepared.n_points,
                n_parameters=prepared.n_parameters,
                diagnostics=diagnostics,
            )
        if not np.all(np.isfinite(execution.parameters)):
            error = ValueError("the optimizer returned non-finite parameters")
            diagnostics = _least_squares_failure_diagnostics(
                prepared, execution.recorder, error
            )
            return _failed_result(
                str(error),
                settings,
                {**metadata, **diagnostics.metadata},
                parameter_names,
                parameter_states,
                n_points=prepared.n_points,
                n_parameters=prepared.n_parameters,
                diagnostics=diagnostics,
            )
        return _build_result(model, prepared, execution)


def _least_squares_options(options: FitOptions) -> LeastSquaresOptions:
    """Translate common controls to typed least-squares options."""
    if isinstance(options, LeastSquaresOptions):
        return options
    if options.method not in {FitMethod.OLS, FitMethod.WLS}:
        raise ValueError(
            "least-squares solver supports only ordinary or weighted least squares"
        )
    option_type = OLSOptions if options.method is FitMethod.OLS else WLSOptions
    return option_type(
        covariance_scaling=options.covariance_policy,
        max_iterations=options.max_iterations,
        ftol=options.ftol,
        xtol=options.xtol,
        gtol=options.gtol,
        metadata=options.metadata,
    )


def _prepare_problem(
    model: BaseFitModel,
    x: ArrayLike,
    y: ArrayLike,
    *,
    initial_parameters: ArrayLike | None,
    bounds: tuple[ArrayLike | float, ArrayLike | float] | None,
    sigma: ArrayLike | None,
    settings: LeastSquaresOptions,
    metadata: dict[str, Any],
    parameter_names: tuple[str, ...],
    parameter_states: tuple[ParameterState, ...],
) -> _PreparedLeastSquares:
    """Validate and normalize one least-squares problem."""
    del metadata, parameter_states  # retained in the caller's failure context
    try:
        x_array, y_array = validate_xy(x, y)
    except ValueError as exc:
        raise _InvalidFitProblem(str(exc)) from exc
    try:
        initial = np.asarray(
            model.initial_guess(x_array, y_array)
            if initial_parameters is None
            else initial_parameters,
            dtype=np.float64,
        )
    except Exception as exc:  # noqa: BLE001 - model diagnostics
        raise _InvalidFitProblem(
            f"could not determine initial fit parameters: {exc}",
            n_points=int(y_array.size),
            n_parameters=len(parameter_names),
        ) from exc
    if initial.ndim != 1 or initial.size == 0:
        raise _InvalidFitProblem(
            "at least one initial parameter is required",
            n_points=int(y_array.size),
        )
    n_parameters = int(initial.size)
    if n_parameters != len(parameter_names):
        raise _InvalidFitProblem(
            "initial parameters do not match the model definition",
            n_points=int(y_array.size),
            n_parameters=n_parameters,
        )
    if y_array.size < n_parameters:
        raise _InvalidFitProblem(
            "the number of fit points is smaller than the number of parameters",
            n_points=int(y_array.size),
            n_parameters=n_parameters,
        )
    try:
        lower, upper = model.bounds(x_array, y_array) if bounds is None else bounds
        fit_bounds = (
            np.broadcast_to(np.asarray(lower, dtype=np.float64), initial.shape).copy(),
            np.broadcast_to(np.asarray(upper, dtype=np.float64), initial.shape).copy(),
        )
    except (TypeError, ValueError) as exc:
        raise _InvalidFitProblem(
            f"invalid parameter bounds: {exc}",
            n_points=int(y_array.size),
            n_parameters=n_parameters,
        ) from exc
    fit_sigma = _prepare_sigma(sigma, y_array, settings, n_parameters)
    return _PreparedLeastSquares(
        x=x_array,
        y=y_array,
        initial=initial,
        bounds=fit_bounds,
        sigma=fit_sigma,
        settings=settings,
        metadata={**model.metadata(), **settings.metadata},
        parameter_names=parameter_names,
        parameter_states=tuple(ParameterState.FREE for _ in parameter_names),
    )


def _prepare_sigma(
    sigma: ArrayLike | None,
    y: np.ndarray,
    settings: LeastSquaresOptions,
    n_parameters: int,
) -> np.ndarray | None:
    """Validate explicit ordinary or weighted least-squares uncertainty use."""
    if settings.method is FitMethod.OLS:
        if sigma is not None:
            raise _InvalidFitProblem(
                "ordinary least squares does not accept sigma; select WLS explicitly",
                n_points=int(y.size),
                n_parameters=n_parameters,
            )
        return None
    if sigma is None:
        raise _InvalidFitProblem(
            "weighted least squares requires dependent-variable sigma",
            n_points=int(y.size),
            n_parameters=n_parameters,
        )
    fit_sigma = np.asarray(sigma, dtype=np.float64)
    if (
        fit_sigma.ndim != 1
        or fit_sigma.shape != y.shape
        or not np.all(np.isfinite(fit_sigma))
        or np.any(fit_sigma <= 0.0)
    ):
        raise _InvalidFitProblem(
            "weighted least-squares sigma must be finite, positive, and match "
            "observations",
            n_points=int(y.size),
            n_parameters=n_parameters,
        )
    return fit_sigma.copy()


def _execute_fit(
    model: BaseFitModel,
    prepared: _PreparedLeastSquares,
) -> _LeastSquaresExecution:
    """Execute SciPy ``curve_fit`` and retain backend diagnostics."""
    fit_kwargs: dict[str, Any] = {}
    if prepared.settings.max_iterations is not None:
        fit_kwargs["maxfev"] = prepared.settings.max_iterations
    for name in ("ftol", "xtol", "gtol"):
        value = getattr(prepared.settings, name)
        if value is not None:
            fit_kwargs[name] = value

    recorder = ModelEvaluationRecorder(
        prepared.y,
        prepared.sigma,
        detailed=solver_debug_enabled(prepared.settings.metadata),
    )
    curve_function = model.curve_function()

    def traced_function(x_values: np.ndarray, *parameters: float) -> np.ndarray:
        fitted = np.asarray(curve_function(x_values, *parameters), dtype=np.float64)
        recorder.evaluate(parameters, fitted)
        return fitted

    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", OptimizeWarning)
            parameters, covariance, info, message, status_code = curve_fit(
                traced_function,
                prepared.x,
                prepared.y,
                p0=prepared.initial,
                bounds=prepared.bounds,
                sigma=prepared.sigma,
                absolute_sigma=(prepared.sigma is not None),
                full_output=True,
                **fit_kwargs,
            )
    except Exception as exc:  # noqa: BLE001 - preserve numerical diagnostics
        raise _LeastSquaresExecutionError(exc, recorder) from exc
    return _LeastSquaresExecution(
        parameters=np.asarray(parameters, dtype=np.float64),
        covariance=np.asarray(covariance, dtype=np.float64),
        warning_messages=[str(item.message) for item in caught],
        status_code=None if status_code is None else int(status_code),
        stop_reason=str(message),
        n_evaluations=(
            None
            if not isinstance(info, dict) or info.get("nfev") is None
            else int(info["nfev"])
        ),
        recorder=recorder,
    )


def _build_result(
    model: BaseFitModel,
    prepared: _PreparedLeastSquares,
    execution: _LeastSquaresExecution,
) -> FitResult:
    """Evaluate predictions and build a method-neutral result contract."""
    parameters = execution.parameters
    covariance = execution.covariance
    warning_messages = list(execution.warning_messages)
    fitted = np.asarray(
        model.curve_function()(prepared.x, *parameters), dtype=np.float64
    )
    residuals = prepared.y - fitted
    metrics = residual_metrics(prepared.y, fitted)
    dof = max(prepared.n_points - prepared.n_parameters, 0)
    standardized = (
        None
        if prepared.sigma is None
        else standardized_residuals(residuals, prepared.sigma)
    )
    chi_square = None if standardized is None else float(np.sum(standardized**2))
    reduced_chi_square = None if chi_square is None or dof == 0 else chi_square / dof
    covariance_scale = covariance_scale_factor(
        prepared.settings.covariance_policy,
        reduced_chi_square,
        weighted=prepared.sigma is not None,
    )
    covariance = np.asarray(covariance, dtype=np.float64) * covariance_scale
    if prepared.sigma is not None and covariance_scale > 1.0 + 1.0e-12:
        warning_messages.append(
            "parameter covariance inflated by reduced chi-square "
            f"factor {covariance_scale:.6g}"
        )
    errors, condition_number, covariance_usable = covariance_errors(covariance)
    quality = assess_fit_quality(
        success=True,
        covariance_usable=covariance_usable,
        warnings_=warning_messages,
        r_squared=metrics["r_squared"],
        condition_number=condition_number,
    )
    message = (
        "fit converged"
        if quality is FitQuality.GOOD
        else "fit converged with diagnostic warnings"
    )
    at_bounds = parameters_at_bounds(parameters, prepared.bounds)
    diagnostics = FitDiagnostics(
        objective=(
            "sum((observed-calculated)^2)"
            if prepared.sigma is None
            else "sum(((observed-calculated)/sigma_y)^2)"
        ),
        weighted=prepared.sigma is not None,
        chi_square=chi_square,
        reduced_chi_square=reduced_chi_square,
        correlation=(covariance_correlation(covariance) if covariance_usable else None),
        standardized_residuals=standardized,
        parameters_at_bounds=at_bounds,
        n_evaluations=execution.n_evaluations,
        stop_reason=execution.stop_reason or message,
        warnings=warning_messages,
        metadata={
            "covariance_scaling": (
                "residual_variance"
                if prepared.sigma is None
                else prepared.settings.covariance_policy.value
            ),
            "covariance_scale_factor": covariance_scale,
            "parameter_error_scale_factor": float(np.sqrt(covariance_scale)),
            "absolute_sigma": (
                False
                if prepared.sigma is None
                else (
                    True
                    if prepared.settings.covariance_policy is CovarianceScaling.ABSOLUTE
                    else (
                        False
                        if prepared.settings.covariance_policy
                        is CovarianceScaling.REDUCED_CHI_SQUARE
                        else None
                    )
                )
            ),
            **problem_summary(
                prepared.x,
                prepared.y,
                prepared.initial,
                prepared.bounds,
                prepared.sigma,
            ),
            "free_parameter_names": list(prepared.parameter_names),
            "detailed_trace_requested": solver_debug_enabled(
                prepared.settings.metadata
            ),
            **termination_metadata(
                category="converged",
                backend="scipy.optimize.curve_fit",
                status_code=execution.status_code,
            ),
            **execution.recorder.as_metadata(),
        },
    )
    return FitResult(
        success=True,
        status=FitStatus.SUCCESS,
        quality=quality,
        message=message,
        parameters=parameters,
        covariance=covariance,
        errors=errors,
        fitted=fitted,
        residuals=residuals,
        rmse=metrics["rmse"],
        mae=metrics["mae"],
        max_abs_error=metrics["max_abs_error"],
        r_squared=metrics["r_squared"],
        n_points=prepared.n_points,
        n_parameters=prepared.n_parameters,
        dof=dof,
        condition_number=condition_number,
        warnings=warning_messages,
        metadata=prepared.metadata,
        method=prepared.settings.method,
        parameter_names=prepared.parameter_names,
        parameter_states=prepared.parameter_states,
        optimizer_parameters=parameters,
        optimizer_covariance=covariance,
        optimizer_errors=errors,
        diagnostics=diagnostics,
    )


def _least_squares_failure_diagnostics(
    prepared: _PreparedLeastSquares,
    recorder: ModelEvaluationRecorder,
    error: Exception,
) -> FitDiagnostics:
    """Build diagnostics for a failed OLS or WLS backend execution."""
    message = str(error)
    lowered = message.lower()
    category = (
        "iteration_limit"
        if "maximum number" in lowered or "maxfev" in lowered
        else "backend_exception"
    )
    metadata = {
        **problem_summary(
            prepared.x, prepared.y, prepared.initial, prepared.bounds, prepared.sigma
        ),
        "free_parameter_names": list(prepared.parameter_names),
        "detailed_trace_requested": solver_debug_enabled(prepared.settings.metadata),
        **termination_metadata(
            category=category,
            backend="scipy.optimize.curve_fit",
            exception=error,
        ),
        **recorder.as_metadata(),
    }
    return FitDiagnostics(
        objective=(
            "sum((observed-calculated)^2)"
            if prepared.sigma is None
            else "sum(((observed-calculated)/sigma_y)^2)"
        ),
        weighted=prepared.sigma is not None,
        n_evaluations=recorder.n_evaluations or None,
        stop_reason=message,
        metadata=metadata,
    )


def _failed_result(
    message: str,
    settings: LeastSquaresOptions,
    metadata: Mapping[str, Any],
    parameter_names: tuple[str, ...],
    parameter_states: tuple[ParameterState, ...],
    *,
    status: FitStatus = FitStatus.FAILED,
    n_points: int = 0,
    n_parameters: int = 0,
    diagnostics: FitDiagnostics | None = None,
) -> FitResult:
    """Build a failed result while retaining method and parameter metadata."""
    return FitResult.failed(
        message,
        status=status,
        n_points=n_points,
        n_parameters=n_parameters,
        metadata=metadata,
        method=settings.method,
        parameter_names=parameter_names,
        parameter_states=parameter_states,
        diagnostics=diagnostics,
    )


def least_squares_fit(
    function: ModelFunction,
    x: ArrayLike,
    y: ArrayLike,
    *,
    p0: Sequence[float],
    bounds: tuple[ArrayLike | float, ArrayLike | float] = (-np.inf, np.inf),
    max_iterations: int | None = None,
    sigma: ArrayLike | None = None,
    method: FitMethod = FitMethod.OLS,
    absolute_sigma: bool | None = None,
    covariance_scaling: CovarianceScaling | str | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> FitResult:
    """Fit a callable function and return structured diagnostics.

    Parameters
    ----------
    function : callable
        Model function with signature ``f(x, *parameters)``.
    x, y : array-like
        Coordinates and observations.
    p0 : sequence of float
        Initial parameter estimate.
    bounds : tuple of array-like, optional
        Lower and upper parameter bounds.
    max_iterations : int, optional
        Maximum number of model evaluations for the SciPy backend.
    sigma : array-like or None, optional
        Dependent-variable uncertainties used only with explicit WLS.
    method : FitMethod, optional
        Ordinary or weighted least squares.
    absolute_sigma : bool or None, optional
        Compatibility alias mapping ``True`` to absolute covariance and
        ``False`` to reduced-chi-square rescaling.
    covariance_scaling : CovarianceScaling, str, or None, optional
        Explicit covariance policy for weighted fits.
    metadata : mapping, optional
        Additional information stored in the result.

    Returns
    -------
    FitResult
        Optimized parameters and diagnostics.
    """
    extra_metadata = dict(metadata or {})
    parameter_order = (
        tuple(str(item) for item in extra_metadata["parameter_order"])
        if "parameter_order" in extra_metadata
        else None
    )
    model_name = str(extra_metadata.get("model", "callable"))
    try:
        model = CallableFitModel(
            function,
            p0,
            name=model_name,
            parameter_names=parameter_order,
            bounds=bounds,
        )
    except (TypeError, ValueError) as exc:
        return FitResult.failed(
            str(exc),
            status=FitStatus.INVALID_INPUT,
            n_parameters=len(tuple(p0)),
            metadata=extra_metadata,
            method=FitMethod(method),
        )
    return LeastSquaresFitter().fit(
        model,
        x,
        y,
        sigma=sigma,
        options=(
            OLSOptions(
                covariance_scaling=(
                    None
                    if covariance_scaling is None
                    else CovarianceScaling(covariance_scaling)
                ),
                absolute_sigma=absolute_sigma,
                max_iterations=max_iterations,
                metadata=extra_metadata,
            )
            if FitMethod(method) is FitMethod.OLS
            else WLSOptions(
                covariance_scaling=(
                    None
                    if covariance_scaling is None
                    else CovarianceScaling(covariance_scaling)
                ),
                absolute_sigma=absolute_sigma,
                max_iterations=max_iterations,
                metadata=extra_metadata,
            )
        ),
    )
