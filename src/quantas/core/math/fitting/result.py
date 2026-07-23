# -*- coding: utf-8 -*-

"""Structured results returned by numerical fitting algorithms."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping

import numpy as np

from .options import FitMethod
from .parameters import ParameterMap, ParameterState


class FitStatus(str, Enum):
    """Execution status of a fit operation.

    Attributes
    ----------
    SUCCESS
        The optimizer converged and returned finite parameters.
    FAILED
        The optimizer did not return a usable solution.
    INVALID_INPUT
        The input data are not suitable for fitting.
    """

    SUCCESS = "success"
    FAILED = "failed"
    INVALID_INPUT = "invalid_input"


class FitQuality(str, Enum):
    """Qualitative assessment of a fitted model.

    Attributes
    ----------
    GOOD
        The fit converged and the diagnostic metrics do not indicate problems.
    POOR
        The fit converged but at least one diagnostic metric is suspicious.
    FAILED
        The fit did not converge or the result is not usable.
    """

    GOOD = "good"
    POOR = "poor"
    FAILED = "failed"


@dataclass(slots=True)
class FitDiagnostics:
    """Method-neutral numerical and statistical diagnostics.

    Parameters
    ----------
    objective : str, optional
        Exact residual or likelihood objective minimized by the solver.
    weighted : bool, optional
        Whether observation uncertainties or weights entered the objective.
    chi_square, reduced_chi_square : float or None, optional
        Weighted residual criteria when defined by the method.
    correlation : ndarray or None, optional
        Parameter correlation matrix in complete reporting order.
    standardized_residuals : ndarray or None, optional
        Residuals divided by the standard deviation used in the objective.
    x_corrections, y_corrections : ndarray or None, optional
        Coordinate corrections returned by errors-in-variables methods.
    jacobian_rank : int or None, optional
        Numerical rank of the local Jacobian or information matrix.
    condition_number : float or None, optional
        Conditioning of the local fitted problem. This is distinct from the
        legacy ``FitResult.condition_number`` field, which may describe the
        covariance returned by an older fitter.
    parameters_at_bounds : tuple of bool, optional
        Flags in complete reporting order.
    n_iterations, n_evaluations : int or None, optional
        Solver work counters.
    stop_reason : str, optional
        Backend-neutral convergence or termination explanation.
    warnings : list of str, optional
        Non-fatal statistical or numerical warnings.
    metadata : dict, optional
        Method-specific passive details.
    """

    objective: str = ""
    weighted: bool = False
    chi_square: float | None = None
    reduced_chi_square: float | None = None
    correlation: np.ndarray | None = None
    standardized_residuals: np.ndarray | None = None
    x_corrections: np.ndarray | None = None
    y_corrections: np.ndarray | None = None
    jacobian_rank: int | None = None
    condition_number: float | None = None
    parameters_at_bounds: tuple[bool, ...] = ()
    n_iterations: int | None = None
    n_evaluations: int | None = None
    stop_reason: str = ""
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize optional arrays and passive containers."""
        for name in (
            "correlation",
            "standardized_residuals",
            "x_corrections",
            "y_corrections",
        ):
            value = getattr(self, name)
            if value is not None:
                setattr(self, name, np.asarray(value, dtype=np.float64).copy())
        self.parameters_at_bounds = tuple(
            bool(value) for value in self.parameters_at_bounds
        )
        self.warnings = [str(value) for value in self.warnings]
        self.metadata = dict(self.metadata)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready diagnostic mapping."""
        return {
            "objective": self.objective,
            "weighted": self.weighted,
            "chi_square": self.chi_square,
            "reduced_chi_square": self.reduced_chi_square,
            "correlation": _array_list(self.correlation),
            "standardized_residuals": _array_list(self.standardized_residuals),
            "x_corrections": _array_list(self.x_corrections),
            "y_corrections": _array_list(self.y_corrections),
            "jacobian_rank": self.jacobian_rank,
            "condition_number": self.condition_number,
            "parameters_at_bounds": list(self.parameters_at_bounds),
            "n_iterations": self.n_iterations,
            "n_evaluations": self.n_evaluations,
            "stop_reason": self.stop_reason,
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class FitResult:
    """Store optimized parameters and numerical fit diagnostics.

    ``parameters`` and ``covariance`` remain the complete model-facing fields
    used by existing Quantas workflows. New mapped solvers also retain their
    reduced optimizer quantities separately so fixed and implied parameters are
    never misrepresented as directly refined variables.

    Parameters
    ----------
    success : bool
        Whether the fit returned usable optimized parameters.
    status : FitStatus
        Execution status of the fit operation.
    quality : FitQuality
        Qualitative assessment derived from the diagnostic metrics.
    message : str, optional
        Human-readable status message.
    parameters : ndarray or None, optional
        Complete optimized or resolved model parameters.
    covariance : ndarray or None, optional
        Complete parameter covariance matrix.
    errors : ndarray or None, optional
        Complete one-standard-deviation parameter uncertainties.
    fitted, residuals : ndarray or None, optional
        Model predictions and observed-minus-predicted residuals.
    method : FitMethod or None, optional
        Regression strategy.
    parameter_names : tuple of str, optional
        Complete reporting order.
    parameter_states : tuple of ParameterState, optional
        FREE/FIXED/IMPLIED/DERIVED state for every complete parameter.
    optimizer_parameters, optimizer_covariance, optimizer_errors : ndarray or
        None, optional
        Reduced free-parameter quantities returned by the numerical backend.
    diagnostics : FitDiagnostics or None, optional
        Extended method-neutral diagnostics.
    """

    success: bool
    status: FitStatus
    quality: FitQuality
    message: str = ""
    parameters: np.ndarray | None = None
    covariance: np.ndarray | None = None
    errors: np.ndarray | None = None
    fitted: np.ndarray | None = None
    residuals: np.ndarray | None = None
    rmse: float | None = None
    mae: float | None = None
    max_abs_error: float | None = None
    r_squared: float | None = None
    n_points: int = 0
    n_parameters: int = 0
    dof: int = 0
    condition_number: float | None = None
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    method: FitMethod | None = None
    parameter_names: tuple[str, ...] = ()
    parameter_states: tuple[ParameterState, ...] = ()
    optimizer_parameters: np.ndarray | None = None
    optimizer_covariance: np.ndarray | None = None
    optimizer_errors: np.ndarray | None = None
    diagnostics: FitDiagnostics | None = None

    def __post_init__(self) -> None:
        """Normalize enums, arrays, and passive containers."""
        self.status = FitStatus(self.status)
        self.quality = FitQuality(self.quality)
        if self.method is not None:
            self.method = FitMethod(self.method)
        self.parameter_names = tuple(str(name) for name in self.parameter_names)
        self.parameter_states = tuple(
            ParameterState(state) for state in self.parameter_states
        )
        for name in (
            "parameters",
            "covariance",
            "errors",
            "fitted",
            "residuals",
            "optimizer_parameters",
            "optimizer_covariance",
            "optimizer_errors",
        ):
            value = getattr(self, name)
            if value is not None:
                setattr(self, name, np.asarray(value, dtype=np.float64).copy())
        self.warnings = [str(value) for value in self.warnings]
        self.metadata = dict(self.metadata)

    @classmethod
    def failed(
        cls,
        message: str,
        *,
        status: FitStatus = FitStatus.FAILED,
        n_points: int = 0,
        n_parameters: int = 0,
        metadata: Mapping[str, Any] | None = None,
        method: FitMethod | None = None,
        parameter_names: tuple[str, ...] = (),
        parameter_states: tuple[ParameterState, ...] = (),
        diagnostics: FitDiagnostics | None = None,
    ) -> FitResult:
        """Create a failed fit result."""
        return cls(
            success=False,
            status=status,
            quality=FitQuality.FAILED,
            message=message,
            n_points=n_points,
            n_parameters=n_parameters,
            dof=max(n_points - n_parameters, 0),
            metadata=dict(metadata or {}),
            method=method,
            parameter_names=parameter_names,
            parameter_states=parameter_states,
            diagnostics=diagnostics,
        )

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable dictionary representation."""
        return {
            "success": self.success,
            "status": self.status.value,
            "quality": self.quality.value,
            "message": self.message,
            "parameters": _array_list(self.parameters),
            "covariance": _array_list(self.covariance),
            "errors": _array_list(self.errors),
            "fitted": _array_list(self.fitted),
            "residuals": _array_list(self.residuals),
            "rmse": self.rmse,
            "mae": self.mae,
            "max_abs_error": self.max_abs_error,
            "r_squared": self.r_squared,
            "n_points": self.n_points,
            "n_parameters": self.n_parameters,
            "dof": self.dof,
            "condition_number": self.condition_number,
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
            "method": None if self.method is None else self.method.value,
            "parameter_names": list(self.parameter_names),
            "parameter_states": [state.value for state in self.parameter_states],
            "optimizer_parameters": _array_list(self.optimizer_parameters),
            "optimizer_covariance": _array_list(self.optimizer_covariance),
            "optimizer_errors": _array_list(self.optimizer_errors),
            "diagnostics": None
            if self.diagnostics is None
            else self.diagnostics.as_dict(),
        }


def fit_result_from_dict(data: Mapping[str, Any]) -> FitResult:
    """Reconstruct a :class:`FitResult` from :meth:`FitResult.as_dict`.

    Parameters
    ----------
    data : mapping
        Serialization-ready fit-result mapping.

    Returns
    -------
    FitResult
        Reconstructed result including complete diagnostics and solver traces.

    Raises
    ------
    ValueError
        If required status or quality fields are missing or invalid.
    """
    diagnostics_data = data.get("diagnostics")
    diagnostics = None
    if isinstance(diagnostics_data, Mapping):
        diagnostics = FitDiagnostics(
            objective=str(diagnostics_data.get("objective", "")),
            weighted=bool(diagnostics_data.get("weighted", False)),
            chi_square=_optional_float(diagnostics_data.get("chi_square")),
            reduced_chi_square=_optional_float(
                diagnostics_data.get("reduced_chi_square")
            ),
            correlation=_optional_array(diagnostics_data.get("correlation")),
            standardized_residuals=_optional_array(
                diagnostics_data.get("standardized_residuals")
            ),
            x_corrections=_optional_array(diagnostics_data.get("x_corrections")),
            y_corrections=_optional_array(diagnostics_data.get("y_corrections")),
            jacobian_rank=_optional_int(diagnostics_data.get("jacobian_rank")),
            condition_number=_optional_float(diagnostics_data.get("condition_number")),
            parameters_at_bounds=tuple(
                bool(value)
                for value in diagnostics_data.get("parameters_at_bounds", ())
            ),
            n_iterations=_optional_int(diagnostics_data.get("n_iterations")),
            n_evaluations=_optional_int(diagnostics_data.get("n_evaluations")),
            stop_reason=str(diagnostics_data.get("stop_reason", "")),
            warnings=[str(value) for value in diagnostics_data.get("warnings", ())],
            metadata=dict(diagnostics_data.get("metadata", {})),
        )

    method_value = data.get("method")
    method = None if method_value is None else FitMethod(str(method_value))
    return FitResult(
        success=bool(data.get("success", False)),
        status=FitStatus(str(data.get("status", FitStatus.FAILED.value))),
        quality=FitQuality(str(data.get("quality", FitQuality.FAILED.value))),
        message=str(data.get("message", "")),
        parameters=_optional_array(data.get("parameters")),
        covariance=_optional_array(data.get("covariance")),
        errors=_optional_array(data.get("errors")),
        fitted=_optional_array(data.get("fitted")),
        residuals=_optional_array(data.get("residuals")),
        rmse=_optional_float(data.get("rmse")),
        mae=_optional_float(data.get("mae")),
        max_abs_error=_optional_float(data.get("max_abs_error")),
        r_squared=_optional_float(data.get("r_squared")),
        n_points=int(data.get("n_points", 0)),
        n_parameters=int(data.get("n_parameters", 0)),
        dof=int(data.get("dof", 0)),
        condition_number=_optional_float(data.get("condition_number")),
        warnings=[str(value) for value in data.get("warnings", ())],
        metadata=dict(data.get("metadata", {})),
        method=method,
        parameter_names=tuple(str(value) for value in data.get("parameter_names", ())),
        parameter_states=tuple(
            ParameterState(str(value)) for value in data.get("parameter_states", ())
        ),
        optimizer_parameters=_optional_array(data.get("optimizer_parameters")),
        optimizer_covariance=_optional_array(data.get("optimizer_covariance")),
        optimizer_errors=_optional_array(data.get("optimizer_errors")),
        diagnostics=diagnostics,
    )


def _optional_array(value: Any) -> np.ndarray | None:
    """Convert a serialized optional numerical array."""
    if value is None:
        return None
    return np.asarray(value, dtype=np.float64)


def _optional_float(value: Any) -> float | None:
    """Convert a serialized optional floating scalar."""
    return None if value is None else float(value)


def _optional_int(value: Any) -> int | None:
    """Convert a serialized optional integer scalar."""
    return None if value is None else int(value)


def _array_list(value: np.ndarray | None) -> list[Any] | None:
    """Convert an optional NumPy array to nested Python lists."""
    return None if value is None else value.tolist()


def resolve_mapped_fit_result(
    result: FitResult,
    parameter_map: ParameterMap,
) -> FitResult:
    """Transform a successful reduced fit result to complete parameters.

    Parameters
    ----------
    result : FitResult
        Result produced for the free optimizer vector.
    parameter_map : ParameterMap
        Mapping used to evaluate the model.

    Returns
    -------
    FitResult
        New result with complete parameters, propagated covariance and errors,
        while retaining the original optimizer quantities.

    Raises
    ------
    ValueError
        If the result is unsuccessful or lacks a compatible free vector.
    """
    if not result.success or result.parameters is None:
        raise ValueError("only successful fit results can be resolved")
    optimizer_parameters = (
        result.parameters
        if result.optimizer_parameters is None
        else result.optimizer_parameters
    )
    if optimizer_parameters.size != parameter_map.n_free:
        raise ValueError("fit result does not match the parameter map")
    complete = parameter_map.expand(optimizer_parameters)
    optimizer_covariance = (
        result.covariance
        if result.optimizer_covariance is None
        else result.optimizer_covariance
    )
    complete_covariance: np.ndarray | None = None
    complete_errors: np.ndarray | None = None
    correlation: np.ndarray | None = None
    if optimizer_covariance is not None:
        complete_covariance = parameter_map.propagate_covariance(
            optimizer_covariance,
            optimizer_parameters,
        )
        complete_errors = np.sqrt(np.maximum(np.diag(complete_covariance), 0.0))
        correlation = _covariance_correlation(complete_covariance)

    diagnostics = result.diagnostics
    if diagnostics is not None:
        diagnostics = FitDiagnostics(
            objective=diagnostics.objective,
            weighted=diagnostics.weighted,
            chi_square=diagnostics.chi_square,
            reduced_chi_square=diagnostics.reduced_chi_square,
            correlation=correlation,
            standardized_residuals=diagnostics.standardized_residuals,
            x_corrections=diagnostics.x_corrections,
            y_corrections=diagnostics.y_corrections,
            jacobian_rank=diagnostics.jacobian_rank,
            condition_number=diagnostics.condition_number,
            parameters_at_bounds=tuple(
                False
                if state is not ParameterState.FREE
                else diagnostics.parameters_at_bounds[
                    parameter_map.free_names.index(name)
                ]
                if diagnostics.parameters_at_bounds
                else False
                for name, state in zip(
                    parameter_map.names, parameter_map.states, strict=True
                )
            ),
            n_iterations=diagnostics.n_iterations,
            n_evaluations=diagnostics.n_evaluations,
            stop_reason=diagnostics.stop_reason,
            warnings=diagnostics.warnings,
            metadata=diagnostics.metadata,
        )

    optimizer_errors = (
        result.errors if result.optimizer_errors is None else result.optimizer_errors
    )
    return FitResult(
        success=result.success,
        status=result.status,
        quality=result.quality,
        message=result.message,
        parameters=complete.values,
        covariance=complete_covariance,
        errors=complete_errors,
        fitted=result.fitted,
        residuals=result.residuals,
        rmse=result.rmse,
        mae=result.mae,
        max_abs_error=result.max_abs_error,
        r_squared=result.r_squared,
        n_points=result.n_points,
        n_parameters=parameter_map.n_free,
        dof=result.dof,
        condition_number=result.condition_number,
        warnings=result.warnings,
        metadata={**result.metadata, "parameter_map": parameter_map.as_dict()},
        method=result.method,
        parameter_names=parameter_map.names,
        parameter_states=parameter_map.states,
        optimizer_parameters=optimizer_parameters,
        optimizer_covariance=optimizer_covariance,
        optimizer_errors=optimizer_errors,
        diagnostics=diagnostics,
    )


def _covariance_correlation(covariance: np.ndarray) -> np.ndarray:
    """Return a correlation matrix without importing diagnostic modules."""
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(errors, errors)
    correlation = np.zeros_like(covariance)
    np.divide(covariance, denominator, out=correlation, where=denominator > 0.0)
    diagonal = np.diag_indices_from(correlation)
    correlation[diagonal] = (errors > 0.0).astype(np.float64)
    return np.clip(correlation, -1.0, 1.0)
