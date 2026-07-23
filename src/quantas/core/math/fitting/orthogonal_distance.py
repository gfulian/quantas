# -*- coding: utf-8 -*-

"""Backend-neutral weighted orthogonal distance regression.

The public classes in this module do not expose :mod:`odrpack` objects.  The
required runtime dependency is imported lazily only when an ODR calculation is
requested, keeping unrelated fitting workflows decoupled from backend objects
while using the modern ODRPACK95 implementation for scientific ODR.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import importlib
import importlib.util
from typing import Any, Protocol

import numpy as np
from numpy.typing import NDArray

from .diagnostics import (
    assess_fit_quality,
    covariance_correlation,
    covariance_errors,
    covariance_scale_factor,
    parameters_at_bounds,
    residual_metrics,
)
from .model import BaseFitModel, MappedFitModel
from .observations import FitObservations
from .options import (
    FitMethod,
    FitOptions,
    ODRDifferenceScheme,
    OrthogonalDistanceOptions,
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
    array_summary,
    problem_summary,
    solver_debug_enabled,
    termination_metadata,
)


class ODRBackendUnavailableError(ImportError):
    """Raised when the optional orthogonal-regression backend is unavailable."""


@dataclass(frozen=True, slots=True)
class ODRBackendResult:
    """Backend-neutral output required from an ODR implementation.

    Parameters
    ----------
    parameters : ndarray
        Optimized free-parameter vector.
    x_corrections, y_corrections : ndarray
        Additive corrections satisfying ``x_adjusted = x_observed + delta``
        and ``y_fitted = y_observed + epsilon``.
    adjusted_x, fitted_y : ndarray
        Corrected explanatory coordinates and corresponding fitted responses.
    covariance : ndarray
        Absolute free-parameter covariance before optional residual-variance
        rescaling.
    residual_variance : float
        Weighted sum of squares divided by residual degrees of freedom.
    sum_square, sum_square_x, sum_square_y : float
        Total weighted ODR objective and its explanatory/response components.
    success : bool
        Whether the backend declared convergence.
    status_code : int
        Backend termination code.
    stop_reason : str
        Human-readable termination reason.
    n_iterations, n_function_evaluations, n_jacobian_evaluations : int
        Backend work counters.
    rank_deficiency : int
        Number of parameter directions dropped because of rank deficiency.
    inverse_condition_number : float
        Backend inverse condition number at termination.
    backend_name, backend_version : str
        Backend provenance.
    """

    parameters: NDArray[np.float64]
    x_corrections: NDArray[np.float64]
    y_corrections: NDArray[np.float64]
    adjusted_x: NDArray[np.float64]
    fitted_y: NDArray[np.float64]
    covariance: NDArray[np.float64]
    residual_variance: float
    sum_square: float
    sum_square_x: float
    sum_square_y: float
    success: bool
    status_code: int
    stop_reason: str
    n_iterations: int
    n_function_evaluations: int
    n_jacobian_evaluations: int
    rank_deficiency: int
    inverse_condition_number: float
    backend_name: str
    backend_version: str


class ODRBackend(Protocol):
    """Protocol implemented by optional orthogonal-regression backends."""

    @property
    def name(self) -> str:
        """Return a stable backend name."""
        ...

    @property
    def version(self) -> str:
        """Return backend version provenance."""
        ...

    def fit(
        self,
        model: Callable[
            [NDArray[np.float64], NDArray[np.float64]], NDArray[np.float64]
        ],
        x: NDArray[np.float64],
        y: NDArray[np.float64],
        initial_parameters: NDArray[np.float64],
        *,
        sigma_x: NDArray[np.float64],
        sigma_y: NDArray[np.float64],
        bounds: tuple[NDArray[np.float64], NDArray[np.float64]],
        options: OrthogonalDistanceOptions,
    ) -> ODRBackendResult:
        """Execute one explicit weighted ODR problem."""
        ...


class ODRPackBackend:
    """Adapt the required :mod:`odrpack` runtime package to Quantas contracts.

    Instances are created through :meth:`load`, which performs the lazy import.
    The adapter always requests an explicit ODR problem, supplies inverse
    variance weights, applies parameter bounds, suppresses backend reports, and
    returns only Quantas-owned passive data.
    """

    def __init__(self, odr_fit: Callable[..., Any], version: str) -> None:
        self._odr_fit = odr_fit
        self._version = str(version)

    @classmethod
    def load(cls) -> ODRPackBackend:
        """Import and validate the Python ODRPACK95 binding.

        Returns
        -------
        ODRPackBackend
            Loaded backend adapter.

        Raises
        ------
        ODRBackendUnavailableError
            If :mod:`odrpack` is not installed or does not expose ``odr_fit``.
        """
        try:
            module = importlib.import_module("odrpack")
        except ModuleNotFoundError as exc:
            if exc.name == "odrpack":
                raise ODRBackendUnavailableError(
                    "orthogonal distance regression requires the Quantas runtime "
                    "dependency 'odrpack'; reinstall Quantas"
                ) from exc
            raise
        function = getattr(module, "odr_fit", None)
        if not callable(function):
            raise ODRBackendUnavailableError(
                "the installed 'odrpack' package does not provide odr_fit"
            )
        return cls(function, str(getattr(module, "__version__", "unknown")))

    @property
    def name(self) -> str:
        """Return the stable backend name."""
        return "odrpack95"

    @property
    def version(self) -> str:
        """Return the Python binding version."""
        return self._version

    def fit(
        self,
        model: Callable[
            [NDArray[np.float64], NDArray[np.float64]], NDArray[np.float64]
        ],
        x: NDArray[np.float64],
        y: NDArray[np.float64],
        initial_parameters: NDArray[np.float64],
        *,
        sigma_x: NDArray[np.float64],
        sigma_y: NDArray[np.float64],
        bounds: tuple[NDArray[np.float64], NDArray[np.float64]],
        options: OrthogonalDistanceOptions,
    ) -> ODRBackendResult:
        """Execute explicit weighted ODR through :func:`odrpack.odr_fit`."""
        result = self._odr_fit(
            model,
            x,
            y,
            initial_parameters,
            weight_x=np.reciprocal(np.square(sigma_x)),
            weight_y=np.reciprocal(np.square(sigma_y)),
            bounds=bounds,
            task="explicit-ODR",
            delta0=options.initial_x_corrections,
            diff_scheme=ODRDifferenceScheme(options.difference_scheme).value,
            report="none",
            maxit=int(options.max_iterations or 50),
            ndigit=options.ndigit,
            taufac=options.trust_region_factor,
            sstol=options.ftol,
            partol=options.xtol,
            step_beta=options.parameter_steps,
            step_delta=options.x_steps,
            scale_beta=options.parameter_scales,
            scale_delta=options.x_scales,
        )
        return ODRBackendResult(
            parameters=np.asarray(result.beta, dtype=np.float64).copy(),
            x_corrections=np.asarray(result.delta, dtype=np.float64).copy(),
            y_corrections=np.asarray(result.eps, dtype=np.float64).copy(),
            adjusted_x=np.asarray(result.xplusd, dtype=np.float64).copy(),
            fitted_y=np.asarray(result.yest, dtype=np.float64).copy(),
            covariance=np.asarray(result.cov_beta, dtype=np.float64).copy(),
            residual_variance=float(result.res_var),
            sum_square=float(result.sum_square),
            sum_square_x=float(result.sum_square_delta),
            sum_square_y=float(result.sum_square_eps),
            success=bool(result.success),
            status_code=int(result.info),
            stop_reason=str(result.stopreason),
            n_iterations=int(result.niter),
            n_function_evaluations=int(result.nfev),
            n_jacobian_evaluations=int(result.njev),
            rank_deficiency=int(result.irank),
            inverse_condition_number=float(result.inv_condnum),
            backend_name=self.name,
            backend_version=self.version,
        )


def odr_backend_available() -> bool:
    """Return whether the :mod:`odrpack` runtime backend can be discovered.

    The function does not import the extension module and therefore is suitable
    for capability checks in future CLI and GUI frontends.
    """
    return importlib.util.find_spec("odrpack") is not None


class OrthogonalDistanceFitter:
    """Fit explicit models by weighted orthogonal distance regression.

    Both selected ``sigma_x`` and ``sigma_y`` vectors are mandatory. The
    objective is the sum of squared corrections normalized by their stated
    standard uncertainties. Fixed and implied physical parameters are removed
    before backend execution through :class:`ParameterMap` and restored after
    the fit with covariance propagation.

    Parameters
    ----------
    backend : ODRBackend or None, optional
        Injected backend, primarily for testing. ``None`` loads the required
        ODRPACK95 runtime adapter lazily on the first fit.
    """

    def __init__(self, backend: ODRBackend | None = None) -> None:
        self._backend = backend

    def fit_problem(
        self,
        model: BaseFitModel,
        observations: FitObservations,
        parameters: ParameterMap,
        options: FitOptions,
    ) -> FitResult:
        """Fit one backend-neutral errors-in-variables problem.

        Parameters
        ----------
        model : BaseFitModel
            Full physical model in non-derived parameter order.
        observations : FitObservations
            Coordinates, responses, both standard-uncertainty vectors, and mask.
        parameters : ParameterMap
            FREE/FIXED/IMPLIED/DERIVED parameter contract.
        options : FitOptions
            General or ODR-specific numerical controls.

        Returns
        -------
        FitResult
            Complete parameters, covariance, coordinate corrections, and ODR
            diagnostics. If the required runtime backend cannot be loaded,
            ``success`` is false and the message reports the installation error.
        """
        metadata = {**model.metadata(), **options.metadata}
        try:
            selected = observations.selected()
            settings = _orthogonal_distance_options(options)
            sigma_x = selected.require_positive_sigma("x")
            sigma_y = selected.require_positive_sigma("y")
            if selected.size <= parameters.n_free:
                raise ValueError(
                    "ODR requires more selected observations than free parameters"
                )
            mapped_model = MappedFitModel(model, parameters)
            initial = parameters.initial_free_values()
            bounds = parameters.free_bounds()
            _validate_option_shapes(
                settings, selected.x.shape, selected.size, parameters.n_free
            )
            backend = self._backend or ODRPackBackend.load()
        except ODRBackendUnavailableError as exc:
            diagnostics = FitDiagnostics(
                objective="weighted orthogonal distance objective",
                weighted=True,
                stop_reason=str(exc),
                metadata={
                    **termination_metadata(
                        category="backend_unavailable",
                        backend="odrpack95",
                        exception=exc,
                    ),
                    "backend_available": False,
                    "free_parameter_names": list(parameters.free_names),
                    "detailed_trace_requested": solver_debug_enabled(options.metadata),
                },
            )
            return FitResult.failed(
                str(exc),
                n_points=observations.selected_size,
                n_parameters=parameters.n_free,
                metadata={**metadata, **diagnostics.metadata},
                method=FitMethod.ODR,
                parameter_names=parameters.names,
                parameter_states=parameters.states,
                diagnostics=diagnostics,
            )
        except (TypeError, ValueError) as exc:
            diagnostics = FitDiagnostics(
                objective="weighted orthogonal distance objective",
                weighted=True,
                stop_reason=str(exc),
                metadata={
                    "free_parameter_names": list(parameters.free_names),
                    "detailed_trace_requested": solver_debug_enabled(options.metadata),
                    **termination_metadata(
                        category="invalid_input",
                        backend="odrpack95",
                        exception=exc,
                    ),
                },
            )
            return FitResult.failed(
                str(exc),
                status=FitStatus.INVALID_INPUT,
                n_points=observations.selected_size,
                n_parameters=parameters.n_free,
                metadata={**metadata, **diagnostics.metadata},
                method=FitMethod.ODR,
                parameter_names=parameters.names,
                parameter_states=parameters.states,
                diagnostics=diagnostics,
            )

        try:
            backend_result = backend.fit(
                lambda x, beta: mapped_model.evaluate(x, beta),
                selected.x,
                selected.y,
                initial,
                sigma_x=sigma_x,
                sigma_y=sigma_y,
                bounds=bounds,
                options=settings,
            )
            reduced = _build_odr_result(
                mapped_model,
                selected,
                sigma_x,
                sigma_y,
                bounds,
                settings,
                backend_result,
            )
        except Exception as exc:  # noqa: BLE001 - preserve backend diagnostics
            diagnostic_metadata = {
                **problem_summary(selected.x, selected.y, initial, bounds, sigma_y),
                "sigma_x_summary": array_summary(sigma_x, positive_ratio=True),
                "sigma_y_summary": array_summary(sigma_y, positive_ratio=True),
                "free_parameter_names": list(parameters.free_names),
                "detailed_trace_requested": solver_debug_enabled(settings.metadata),
                **termination_metadata(
                    category="backend_exception",
                    backend=backend.name,
                    backend_version=backend.version,
                    exception=exc,
                ),
            }
            diagnostics = FitDiagnostics(
                objective="weighted orthogonal distance objective",
                weighted=True,
                stop_reason=str(exc),
                metadata=diagnostic_metadata,
            )
            return FitResult.failed(
                f"orthogonal distance regression failed: {exc}",
                n_points=selected.size,
                n_parameters=parameters.n_free,
                metadata={**metadata, **diagnostic_metadata},
                method=FitMethod.ODR,
                parameter_names=parameters.names,
                parameter_states=parameters.states,
                diagnostics=diagnostics,
            )
        if not reduced.success:
            return reduced
        return resolve_mapped_fit_result(reduced, parameters)


def _orthogonal_distance_options(options: FitOptions) -> OrthogonalDistanceOptions:
    """Normalize common controls for orthogonal distance regression."""
    if isinstance(options, OrthogonalDistanceOptions):
        return options
    return OrthogonalDistanceOptions(
        covariance_scaling=options.covariance_policy,
        max_iterations=options.max_iterations,
        ftol=options.ftol,
        xtol=options.xtol,
        gtol=options.gtol,
        metadata=options.metadata,
    )


def _validate_option_shapes(
    options: OrthogonalDistanceOptions,
    x_shape: tuple[int, ...],
    n_points: int,
    n_parameters: int,
) -> None:
    """Validate ODR option arrays after problem dimensions are known."""
    expected_shapes = {
        "initial_x_corrections": x_shape,
        "x_steps": x_shape,
        "x_scales": x_shape,
        "parameter_steps": (n_parameters,),
        "parameter_scales": (n_parameters,),
    }
    for name, shape in expected_shapes.items():
        value = getattr(options, name)
        if value is None:
            continue
        if value.shape != shape:
            # Preserve the historical one-dimensional convenience for one x.
            if not (len(x_shape) == 1 and value.shape == (n_points,)):
                raise ValueError(f"{name} must have shape {shape}")


@dataclass(frozen=True, slots=True)
class _ODRAnalysis:
    """Post-processed numerical quantities for one successful ODR result."""

    fitted: NDArray[np.float64]
    residuals: NDArray[np.float64]
    metrics: dict[str, float]
    standardized_x: NDArray[np.float64]
    standardized_y: NDArray[np.float64]
    standardized_orthogonal: NDArray[np.float64]
    chi_square: float
    reduced_chi_square: float | None
    covariance_scale: float
    covariance: NDArray[np.float64]
    errors: NDArray[np.float64] | None
    covariance_condition: float | None
    covariance_usable: bool
    backend_condition: float
    warnings: tuple[str, ...]
    parameters_at_bounds: tuple[bool, ...]


def _build_odr_result(
    model: MappedFitModel,
    observations: FitObservations,
    sigma_x: NDArray[np.float64],
    sigma_y: NDArray[np.float64],
    bounds: tuple[NDArray[np.float64], NDArray[np.float64]],
    options: OrthogonalDistanceOptions,
    backend: ODRBackendResult,
) -> FitResult:
    """Convert one backend output to the Quantas fit-result contract."""
    n_points = observations.size
    n_parameters = len(model.parameter_names)
    if not backend.success:
        return _failed_backend_result(
            model,
            observations,
            sigma_x,
            sigma_y,
            bounds,
            options,
            backend,
            n_points,
            n_parameters,
        )
    _validate_backend_result(backend, observations.x.shape, n_points, n_parameters)
    analysis = _analyze_odr_result(
        observations,
        sigma_x,
        sigma_y,
        bounds,
        options,
        backend,
    )
    quality = assess_fit_quality(
        success=True,
        covariance_usable=analysis.covariance_usable,
        warnings_=analysis.warnings,
        r_squared=analysis.metrics["r_squared"],
        condition_number=analysis.backend_condition,
    )
    message = (
        "orthogonal distance regression converged"
        if quality is FitQuality.GOOD
        else "orthogonal distance regression converged with diagnostic warnings"
    )
    diagnostics = _odr_diagnostics(
        model,
        observations,
        sigma_x,
        sigma_y,
        bounds,
        options,
        backend,
        analysis,
        n_parameters=n_parameters,
    )
    dof = max(n_points - n_parameters, 0)
    return FitResult(
        success=True,
        status=FitStatus.SUCCESS,
        quality=quality,
        message=message,
        parameters=backend.parameters,
        covariance=analysis.covariance,
        errors=analysis.errors,
        fitted=analysis.fitted,
        residuals=analysis.residuals,
        rmse=analysis.metrics["rmse"],
        mae=analysis.metrics["mae"],
        max_abs_error=analysis.metrics["max_abs_error"],
        r_squared=analysis.metrics["r_squared"],
        n_points=n_points,
        n_parameters=n_parameters,
        dof=dof,
        condition_number=analysis.covariance_condition,
        warnings=list(analysis.warnings),
        metadata={
            **model.metadata(),
            **options.metadata,
            "backend": backend.backend_name,
            "backend_version": backend.backend_version,
            "adjusted_x": backend.adjusted_x.tolist(),
        },
        method=FitMethod.ODR,
        parameter_names=model.parameter_names,
        parameter_states=tuple(ParameterState.FREE for _ in model.parameter_names),
        optimizer_parameters=backend.parameters,
        optimizer_covariance=analysis.covariance,
        optimizer_errors=analysis.errors,
        diagnostics=diagnostics,
    )


def _failed_backend_result(
    model: MappedFitModel,
    observations: FitObservations,
    sigma_x: NDArray[np.float64],
    sigma_y: NDArray[np.float64],
    bounds: tuple[NDArray[np.float64], NDArray[np.float64]],
    options: OrthogonalDistanceOptions,
    backend: ODRBackendResult,
    n_points: int,
    n_parameters: int,
) -> FitResult:
    """Return a structured non-convergence result from the ODR backend."""
    metadata = {
        **model.metadata(),
        **options.metadata,
        **problem_summary(
            observations.x,
            observations.y,
            model.initial_guess(observations.x, observations.y),
            bounds,
            sigma_y,
        ),
        "sigma_x_summary": array_summary(sigma_x, positive_ratio=True),
        "sigma_y_summary": array_summary(sigma_y, positive_ratio=True),
        "free_parameter_names": list(model.parameter_names),
        "detailed_trace_requested": solver_debug_enabled(options.metadata),
        **termination_metadata(
            category=(
                "iteration_limit"
                if backend.n_iterations >= int(options.max_iterations or 50)
                else "backend_nonconvergence"
            ),
            backend=backend.backend_name,
            backend_version=backend.backend_version,
            status_code=backend.status_code,
        ),
        "last_valid_parameters": backend.parameters.tolist(),
        "last_x_corrections": backend.x_corrections.tolist(),
        "last_y_corrections": backend.y_corrections.tolist(),
        "sum_square": backend.sum_square,
        "sum_square_x": backend.sum_square_x,
        "sum_square_y": backend.sum_square_y,
        "residual_variance": backend.residual_variance,
        "rank_deficiency": backend.rank_deficiency,
        "inverse_condition_number": backend.inverse_condition_number,
        "jacobian_evaluations": backend.n_jacobian_evaluations,
    }
    diagnostics = FitDiagnostics(
        objective=(
            "sum_i(sum_j((delta_x_j,i/sigma_x_j,i)^2) + "
            "((f(x+delta_x,beta)-y)/sigma_y)^2)"
        ),
        weighted=True,
        chi_square=backend.sum_square,
        condition_number=(
            float(np.inf)
            if backend.inverse_condition_number <= 0.0
            else float(1.0 / backend.inverse_condition_number)
        ),
        jacobian_rank=max(n_parameters - backend.rank_deficiency, 0),
        n_iterations=backend.n_iterations,
        n_evaluations=backend.n_function_evaluations,
        stop_reason=backend.stop_reason,
        metadata=metadata,
    )
    return FitResult.failed(
        f"ODR backend did not converge: {backend.stop_reason}",
        n_points=n_points,
        n_parameters=n_parameters,
        metadata=metadata,
        method=FitMethod.ODR,
        parameter_names=model.parameter_names,
        parameter_states=tuple(ParameterState.FREE for _ in model.parameter_names),
        diagnostics=diagnostics,
    )


def _analyze_odr_result(
    observations: FitObservations,
    sigma_x: NDArray[np.float64],
    sigma_y: NDArray[np.float64],
    bounds: tuple[NDArray[np.float64], NDArray[np.float64]],
    options: OrthogonalDistanceOptions,
    backend: ODRBackendResult,
) -> _ODRAnalysis:
    """Calculate residual, covariance, conditioning, and warning quantities."""
    fitted = np.asarray(backend.fitted_y, dtype=np.float64)
    residuals = observations.y - fitted
    metrics = residual_metrics(observations.y, fitted)
    standardized_x = backend.x_corrections / sigma_x
    standardized_y = backend.y_corrections / sigma_y
    x_component = (
        np.square(standardized_x)
        if standardized_x.ndim == 1
        else np.sum(np.square(standardized_x), axis=0)
    )
    standardized_orthogonal = np.sign(residuals) * np.sqrt(
        x_component + np.square(standardized_y)
    )
    dof = max(observations.size - backend.parameters.size, 0)
    chi_square = float(backend.sum_square)
    reduced = None if dof == 0 else chi_square / float(dof)
    scale = covariance_scale_factor(
        options.covariance_policy,
        reduced,
        weighted=True,
    )
    covariance = np.asarray(backend.covariance, dtype=np.float64) * scale
    warnings = _odr_warnings(scale, backend.rank_deficiency)
    errors, condition, usable = covariance_errors(covariance)
    backend_condition = (
        float(np.inf)
        if backend.inverse_condition_number <= 0.0
        else float(1.0 / backend.inverse_condition_number)
    )
    return _ODRAnalysis(
        fitted=fitted,
        residuals=residuals,
        metrics=metrics,
        standardized_x=standardized_x,
        standardized_y=standardized_y,
        standardized_orthogonal=standardized_orthogonal,
        chi_square=chi_square,
        reduced_chi_square=reduced,
        covariance_scale=scale,
        covariance=covariance,
        errors=errors,
        covariance_condition=condition,
        covariance_usable=usable,
        backend_condition=backend_condition,
        warnings=warnings,
        parameters_at_bounds=parameters_at_bounds(backend.parameters, bounds),
    )


def _odr_warnings(scale: float, rank_deficiency: int) -> tuple[str, ...]:
    """Return stable scientific warnings for a successful ODR result."""
    messages: list[str] = []
    if scale > 1.0 + 1.0e-12:
        messages.append(
            f"parameter covariance inflated by reduced chi-square factor {scale:.6g}"
        )
    if rank_deficiency > 0:
        messages.append(f"ODR Jacobian has rank deficiency {rank_deficiency}")
    return tuple(messages)


def _odr_diagnostics(
    model: MappedFitModel,
    observations: FitObservations,
    sigma_x: NDArray[np.float64],
    sigma_y: NDArray[np.float64],
    bounds: tuple[NDArray[np.float64], NDArray[np.float64]],
    options: OrthogonalDistanceOptions,
    backend: ODRBackendResult,
    analysis: _ODRAnalysis,
    *,
    n_parameters: int,
) -> FitDiagnostics:
    """Build the common diagnostic contract from analyzed ODR quantities."""
    return FitDiagnostics(
        objective=(
            "sum_i(sum_j((delta_x_j,i/sigma_x_j,i)^2) + "
            "((f(x+delta_x,beta)-y)/sigma_y)^2)"
        ),
        weighted=True,
        chi_square=analysis.chi_square,
        reduced_chi_square=analysis.reduced_chi_square,
        correlation=(
            covariance_correlation(analysis.covariance)
            if analysis.covariance_usable
            else None
        ),
        standardized_residuals=analysis.standardized_orthogonal,
        x_corrections=backend.x_corrections,
        y_corrections=backend.y_corrections,
        jacobian_rank=max(n_parameters - backend.rank_deficiency, 0),
        condition_number=analysis.backend_condition,
        parameters_at_bounds=analysis.parameters_at_bounds,
        n_iterations=backend.n_iterations,
        n_evaluations=backend.n_function_evaluations,
        stop_reason=backend.stop_reason,
        warnings=list(analysis.warnings),
        metadata=_odr_diagnostic_metadata(
            model,
            observations,
            sigma_x,
            sigma_y,
            bounds,
            options,
            backend,
            analysis,
        ),
    )


def _odr_diagnostic_metadata(
    model: MappedFitModel,
    observations: FitObservations,
    sigma_x: NDArray[np.float64],
    sigma_y: NDArray[np.float64],
    bounds: tuple[NDArray[np.float64], NDArray[np.float64]],
    options: OrthogonalDistanceOptions,
    backend: ODRBackendResult,
    analysis: _ODRAnalysis,
) -> dict[str, Any]:
    """Return serialization-ready backend and objective provenance."""
    return {
        **problem_summary(
            observations.x,
            observations.y,
            model.initial_guess(observations.x, observations.y),
            bounds,
            sigma_y,
        ),
        "sigma_x_summary": array_summary(sigma_x, positive_ratio=True),
        "sigma_y_summary": array_summary(sigma_y, positive_ratio=True),
        "free_parameter_names": list(model.parameter_names),
        "detailed_trace_requested": solver_debug_enabled(options.metadata),
        **termination_metadata(
            category="converged",
            backend=backend.backend_name,
            backend_version=backend.backend_version,
            status_code=backend.status_code,
        ),
        "rank_deficiency": backend.rank_deficiency,
        "inverse_condition_number": backend.inverse_condition_number,
        "jacobian_evaluations": backend.n_jacobian_evaluations,
        "sum_square_x": backend.sum_square_x,
        "sum_square_y": backend.sum_square_y,
        "residual_variance": backend.residual_variance,
        "standardized_x_corrections": analysis.standardized_x.tolist(),
        "standardized_y_corrections": analysis.standardized_y.tolist(),
        "x_correction_convention": "adjusted_x-observed_x",
        "y_correction_convention": "fitted_y-observed_y",
        "covariance_scaling": options.covariance_policy.value,
        "covariance_scale_factor": analysis.covariance_scale,
        "parameter_error_scale_factor": float(np.sqrt(analysis.covariance_scale)),
        "difference_scheme": ODRDifferenceScheme(options.difference_scheme).value,
        "x_jacobian": "ODRPACK95 finite differences",
        "parameter_jacobian": "ODRPACK95 finite differences",
    }


def _validate_backend_result(
    result: ODRBackendResult,
    x_shape: tuple[int, ...],
    n_points: int,
    n_parameters: int,
) -> None:
    """Reject malformed or non-finite backend outputs."""
    exact_shapes = {
        "parameters": (result.parameters, (n_parameters,)),
        "x_corrections": (result.x_corrections, x_shape),
        "y_corrections": (result.y_corrections, (n_points,)),
        "adjusted_x": (result.adjusted_x, x_shape),
        "fitted_y": (result.fitted_y, (n_points,)),
    }
    for name, (value, shape) in exact_shapes.items():
        array = np.asarray(value, dtype=np.float64)
        if array.shape != shape or not np.all(np.isfinite(array)):
            raise ValueError(f"ODR backend returned invalid {name}")
    covariance = np.asarray(result.covariance, dtype=np.float64)
    if covariance.shape != (n_parameters, n_parameters):
        raise ValueError("ODR backend returned covariance with invalid shape")
    if not np.all(np.isfinite(covariance)):
        raise ValueError("ODR backend returned non-finite covariance")
    for name in (
        "residual_variance",
        "sum_square",
        "sum_square_x",
        "sum_square_y",
        "inverse_condition_number",
    ):
        scalar_value = float(getattr(result, name))
        if not np.isfinite(scalar_value) or scalar_value < 0.0:
            raise ValueError(f"ODR backend returned invalid {name}")


__all__ = [
    "ODRBackend",
    "ODRBackendResult",
    "ODRBackendUnavailableError",
    "ODRPackBackend",
    "OrthogonalDistanceFitter",
    "odr_backend_available",
]
