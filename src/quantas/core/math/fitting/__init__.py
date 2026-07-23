# -*- coding: utf-8 -*-

"""General fitting models, algorithms, results and diagnostics."""

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
from .effective_variance import EffectiveVarianceFitter
from .least_squares import LeastSquaresFitter, least_squares_fit
from .orthogonal_distance import (
    ODRBackend,
    ODRBackendResult,
    ODRBackendUnavailableError,
    ODRPackBackend,
    OrthogonalDistanceFitter,
    odr_backend_available,
)
from .model import BaseFitModel, CallableFitModel, FittedModel, MappedFitModel
from .observations import FitObservations
from .options import (
    CovarianceScaling,
    EffectiveVarianceOptions,
    FitMethod,
    FitOptions,
    LeastSquaresOptions,
    OLSOptions,
    ODROptions,
    SolverOptions,
    WLSOptions,
    default_solver_options,
    ODRDifferenceScheme,
    OrthogonalDistanceOptions,
)
from .parameters import (
    ParameterDefinition,
    ParameterMap,
    ParameterSet,
    ParameterState,
)
from .result import (
    FitDiagnostics,
    FitQuality,
    FitResult,
    FitStatus,
    fit_result_from_dict,
    resolve_mapped_fit_result,
)
from .solver import FitSolver
from .solver_debug import SOLVER_DEBUG_METADATA_KEY

__all__ = [
    "BaseFitModel",
    "CallableFitModel",
    "CovarianceScaling",
    "EffectiveVarianceFitter",
    "EffectiveVarianceOptions",
    "FitDiagnostics",
    "FitMethod",
    "FitObservations",
    "FitOptions",
    "FitQuality",
    "FitResult",
    "FitSolver",
    "FitStatus",
    "fit_result_from_dict",
    "SOLVER_DEBUG_METADATA_KEY",
    "FittedModel",
    "LeastSquaresFitter",
    "LeastSquaresOptions",
    "OLSOptions",
    "ODROptions",
    "SolverOptions",
    "WLSOptions",
    "default_solver_options",
    "MappedFitModel",
    "ParameterDefinition",
    "ParameterMap",
    "ParameterSet",
    "ParameterState",
    "assess_fit_quality",
    "covariance_correlation",
    "covariance_errors",
    "least_squares_fit",
    "residual_metrics",
    "resolve_mapped_fit_result",
    "standardized_residuals",
    "ODRBackend",
    "ODRBackendResult",
    "ODRBackendUnavailableError",
    "ODRDifferenceScheme",
    "ODRPackBackend",
    "OrthogonalDistanceFitter",
    "OrthogonalDistanceOptions",
    "covariance_scale_factor",
    "odr_backend_available",
    "parameters_at_bounds",
    "validate_xy",
]
