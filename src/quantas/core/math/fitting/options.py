# -*- coding: utf-8 -*-

"""Passive, typed configuration objects for numerical fitting solvers."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, TypeAlias

import numpy as np


class FitMethod(str, Enum):
    """Supported general regression strategies.

    Attributes
    ----------
    OLS
        Ordinary nonlinear least squares using vertical residuals.
    WLS
        Weighted nonlinear least squares using dependent-variable standard
        uncertainties.
    EFFECTIVE_VARIANCE
        Iteratively reweighted least squares including projected uncertainty
        from the independent variable.
    ODR
        Orthogonal distance regression or errors-in-variables fitting.
    """

    OLS = "ordinary_least_squares"
    WLS = "weighted_least_squares"
    EFFECTIVE_VARIANCE = "effective_variance"
    ODR = "orthogonal_distance_regression"


class ODRDifferenceScheme(str, Enum):
    """Finite-difference scheme used by an ODR backend.

    Attributes
    ----------
    FORWARD
        Forward finite differences. This is the native ODRPACK95 default.
    CENTRAL
        Central finite differences. Quantas uses this as the scientific
        default because it is generally more accurate, at the cost of extra
        model evaluations.
    """

    FORWARD = "forward"
    CENTRAL = "central"


class CovarianceScaling(str, Enum):
    """Policy applied to covariance from weighted regression.

    Attributes
    ----------
    ABSOLUTE
        Preserve the absolute scale of supplied standard uncertainties.
    REDUCED_CHI_SQUARE
        Multiply covariance by the reduced chi-square. Parameter standard
        errors may therefore increase or decrease.
    INFLATE_ONLY
        Multiply covariance by ``max(1, reduced_chi_square)``. This reproduces
        the uncertainty-rescaling convention used by EosFit7.
    """

    ABSOLUTE = "absolute"
    REDUCED_CHI_SQUARE = "reduced_chi_square"
    INFLATE_ONLY = "inflate_only"


@dataclass(slots=True)
class FitOptions:
    """Base contract shared by all numerical solver options.

    This class supplies the names and semantics that are common to every
    Quantas fitting solver. Public workflow APIs should normally receive one
    of the typed subclasses such as :class:`OLSOptions`, :class:`WLSOptions`,
    :class:`EffectiveVarianceOptions`, or
    :class:`OrthogonalDistanceOptions`.

    Parameters
    ----------
    method : FitMethod, optional
        Regression strategy represented by the options object.
    covariance_scaling : CovarianceScaling or None, optional
        Policy for covariance returned by a weighted fit. ``None`` selects
        :attr:`CovarianceScaling.ABSOLUTE` unless ``absolute_sigma`` is used.
    absolute_sigma : bool or None, optional
        Compatibility alias for the former two-state covariance contract.
        ``True`` maps to :attr:`CovarianceScaling.ABSOLUTE`; ``False`` maps to
        :attr:`CovarianceScaling.REDUCED_CHI_SQUARE`.
    max_iterations : int or None, optional
        Maximum iterations allowed by the represented solver. For SciPy
        nonlinear least squares this value is passed as the maximum model
        evaluation count because that is the native backend stopping control.
        For effective variance it is the number of outer reweighting cycles;
        for ODR it is the ODR iteration limit.
    ftol, xtol, gtol : float or None, optional
        Function, parameter, and gradient convergence tolerances where
        supported by the represented solver.
    metadata : dict, optional
        Passive information copied into the numerical result.

    Raises
    ------
    ValueError
        If iteration limits, tolerances, or covariance controls are invalid.
    """

    method: FitMethod = FitMethod.OLS
    covariance_scaling: CovarianceScaling | None = None
    absolute_sigma: bool | None = None
    max_iterations: int | None = None
    ftol: float | None = None
    xtol: float | None = None
    gtol: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    _covariance_scaling_explicit: bool = field(default=False, init=False, repr=False)

    def __post_init__(self) -> None:
        """Normalize enum values, mappings, and common numerical controls."""
        self.method = FitMethod(self.method)
        self._covariance_scaling_explicit = (
            self.covariance_scaling is not None or self.absolute_sigma is not None
        )
        self.covariance_scaling = _resolve_covariance_scaling(
            self.covariance_scaling,
            self.absolute_sigma,
        )
        if self.max_iterations is not None and self.max_iterations <= 0:
            raise ValueError("max_iterations must be positive")
        for name in ("ftol", "xtol", "gtol"):
            value = getattr(self, name)
            if value is not None and (not np.isfinite(value) or value <= 0.0):
                raise ValueError(f"{name} must be a positive finite value")
        self.metadata = dict(self.metadata)

    @property
    def covariance_scaling_explicit(self) -> bool:
        """Return whether the caller explicitly selected a covariance policy."""
        return bool(self._covariance_scaling_explicit)

    @property
    def covariance_policy(self) -> CovarianceScaling:
        """Return the normalized covariance-scaling policy.

        Returns
        -------
        CovarianceScaling
            Non-optional policy after dataclass validation.
        """
        policy = self.covariance_scaling
        if policy is None:  # defensive guard for static type checkers
            return CovarianceScaling.ABSOLUTE
        return CovarianceScaling(policy)

    def with_metadata(self, metadata: dict[str, Any]) -> FitOptions:
        """Return an independent options object with merged metadata.

        Parameters
        ----------
        metadata : dict
            Values merged after the existing options metadata.

        Returns
        -------
        FitOptions
            Deep copy preserving the concrete solver-options type.
        """
        clone = deepcopy(self)
        clone.metadata = {**self.metadata, **metadata}
        return clone

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready description of the solver options.

        Returns
        -------
        dict
            Common fields plus solver-specific fields supplied by subclasses.
        """
        payload: dict[str, Any] = {
            "type": type(self).__name__,
            "method": self.method.value,
            "covariance_scaling": self.covariance_policy.value,
            "max_iterations": self.max_iterations,
            "ftol": self.ftol,
            "xtol": self.xtol,
            "gtol": self.gtol,
            "metadata": dict(self.metadata),
        }
        payload.update(self._specific_as_dict())
        return payload

    def _specific_as_dict(self) -> dict[str, Any]:
        """Return subclass-specific serialization fields."""
        return {}


@dataclass(slots=True)
class LeastSquaresOptions(FitOptions):
    """Base typed options for SciPy nonlinear least-squares solvers.

    New public code should use :class:`OLSOptions` or :class:`WLSOptions`.
    This common base remains useful to the general mathematical implementation.

    Raises
    ------
    ValueError
        If a method other than OLS or WLS is requested.
    """

    def __post_init__(self) -> None:
        """Validate least-squares-specific settings."""
        FitOptions.__post_init__(self)
        if self.method not in {FitMethod.OLS, FitMethod.WLS}:
            raise ValueError(
                "LeastSquaresOptions supports only ordinary or weighted least squares"
            )


@dataclass(slots=True)
class OLSOptions(LeastSquaresOptions):
    """Configure ordinary nonlinear least squares.

    The solver minimizes unweighted vertical residuals. The ``method`` is
    fixed by the class and therefore cannot disagree with the selected option
    contract.
    """

    method: FitMethod = FitMethod.OLS

    def __post_init__(self) -> None:
        """Validate that this contract represents OLS."""
        LeastSquaresOptions.__post_init__(self)
        if self.method is not FitMethod.OLS:
            raise ValueError("OLSOptions requires method='ordinary_least_squares'")


@dataclass(slots=True)
class WLSOptions(LeastSquaresOptions):
    """Configure weighted nonlinear least squares.

    The solver minimizes vertical residuals divided by dependent-variable
    standard uncertainties. The ``method`` is fixed by the class.
    """

    method: FitMethod = FitMethod.WLS

    def __post_init__(self) -> None:
        """Validate that this contract represents WLS."""
        LeastSquaresOptions.__post_init__(self)
        if self.method is not FitMethod.WLS:
            raise ValueError("WLSOptions requires method='weighted_least_squares'")


@dataclass(slots=True)
class EffectiveVarianceOptions(FitOptions):
    """Configure iterative effective-variance regression.

    Effective variance projects independent-variable uncertainty into the
    dependent-variable direction through the local model derivative and
    updates the effective standard uncertainties after every nonlinear WLS
    cycle.

    Parameters
    ----------
    inner_max_iterations : int or None, optional
        Maximum model-evaluation count for each inner WLS cycle. ``None``
        leaves the SciPy backend default unchanged.
    parameter_rtol, parameter_atol : float, optional
        Relative and absolute convergence tolerances for the free-parameter
        vector between reweighting cycles.
    sigma_rtol, sigma_atol : float, optional
        Relative and absolute convergence tolerances for the effective standard
        uncertainty vector.
    max_iterations : int or None, optional
        Maximum number of outer reweighting cycles. ``None`` selects 25.
    covariance_scaling, absolute_sigma, ftol, xtol, gtol, metadata
        Common controls inherited from :class:`FitOptions`. The tolerances
        ``ftol``, ``xtol``, and ``gtol`` are passed to each inner WLS cycle.

    Raises
    ------
    ValueError
        If convergence tolerances are invalid.
    """

    method: FitMethod = FitMethod.EFFECTIVE_VARIANCE
    inner_max_iterations: int | None = None
    parameter_rtol: float = 1.0e-10
    parameter_atol: float = 1.0e-12
    sigma_rtol: float = 1.0e-10
    sigma_atol: float = 1.0e-12

    def __post_init__(self) -> None:
        """Validate effective-variance-specific controls."""
        FitOptions.__post_init__(self)
        if self.method is not FitMethod.EFFECTIVE_VARIANCE:
            raise ValueError(
                "EffectiveVarianceOptions requires method='effective_variance'"
            )
        if self.inner_max_iterations is not None and self.inner_max_iterations <= 0:
            raise ValueError("inner_max_iterations must be positive")
        for name in (
            "parameter_rtol",
            "parameter_atol",
            "sigma_rtol",
            "sigma_atol",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be a finite non-negative value")
            setattr(self, name, value)
        if self.max_iterations is None:
            self.max_iterations = 25

    def _specific_as_dict(self) -> dict[str, Any]:
        """Return effective-variance convergence controls."""
        return {
            "inner_max_iterations": self.inner_max_iterations,
            "parameter_rtol": self.parameter_rtol,
            "parameter_atol": self.parameter_atol,
            "sigma_rtol": self.sigma_rtol,
            "sigma_atol": self.sigma_atol,
        }


@dataclass(slots=True)
class OrthogonalDistanceOptions(FitOptions):
    """Configure weighted orthogonal distance regression.

    The options are backend-neutral even though the first Quantas adapter uses
    ODRPACK95 through the required :mod:`odrpack` runtime package. The fitter owns the
    regression task, weights, bounds, and silence policy; callers cannot
    request backend terminal reports from the numerical core.

    Parameters
    ----------
    difference_scheme : ODRDifferenceScheme or str, optional
        Finite-difference scheme for derivatives not supplied analytically.
        Quantas defaults to central differences.
    ndigit : int or None, optional
        Number of reliable decimal digits in model evaluations. ``None`` asks
        the backend to estimate it.
    trust_region_factor : float or None, optional
        Initial trust-region radius factor in the interval ``(0, 1]``.
    initial_x_corrections : ndarray or None, optional
        Starting corrections to the selected explanatory coordinates.
    parameter_steps, x_steps : ndarray or None, optional
        Positive relative finite-difference step sizes.
    parameter_scales, x_scales : ndarray or None, optional
        Positive scaling values used to improve numerical conditioning.
    max_iterations : int or None, optional
        Maximum ODR iterations. ``None`` selects 50.
    ftol, xtol : float or None, optional
        Sum-of-squares and parameter convergence tolerances mapped to the
        ODRPACK95 controls ``sstol`` and ``partol``.
    gtol : float or None, optional
        Unsupported for ODR and therefore required to remain ``None``.
    covariance_scaling, absolute_sigma, metadata
        Common controls inherited from :class:`FitOptions`.

    Raises
    ------
    ValueError
        If tolerances, arrays, or ODR controls are invalid.
    """

    method: FitMethod = FitMethod.ODR
    difference_scheme: ODRDifferenceScheme | str = ODRDifferenceScheme.CENTRAL
    ndigit: int | None = None
    trust_region_factor: float | None = None
    initial_x_corrections: np.ndarray | None = None
    parameter_steps: np.ndarray | None = None
    x_steps: np.ndarray | None = None
    parameter_scales: np.ndarray | None = None
    x_scales: np.ndarray | None = None

    def __post_init__(self) -> None:
        """Validate ODR-specific controls without importing the backend."""
        FitOptions.__post_init__(self)
        if self.method is not FitMethod.ODR:
            raise ValueError(
                "OrthogonalDistanceOptions requires "
                "method='orthogonal_distance_regression'"
            )
        self.difference_scheme = ODRDifferenceScheme(self.difference_scheme)
        if self.gtol is not None:
            raise ValueError("orthogonal distance regression does not support gtol")
        if self.max_iterations is None:
            self.max_iterations = 50
        if self.ndigit is not None and self.ndigit <= 0:
            raise ValueError("ndigit must be positive")
        if self.trust_region_factor is not None:
            value = float(self.trust_region_factor)
            if not np.isfinite(value) or not 0.0 < value <= 1.0:
                raise ValueError("trust_region_factor must be finite and lie in (0, 1]")
            self.trust_region_factor = value
        for name in (
            "initial_x_corrections",
            "parameter_steps",
            "x_steps",
            "parameter_scales",
            "x_scales",
        ):
            value = getattr(self, name)
            if value is None:
                continue
            array = np.asarray(value, dtype=np.float64)
            if (
                array.ndim not in {1, 2}
                or array.size == 0
                or not np.all(np.isfinite(array))
            ):
                raise ValueError(
                    f"{name} must be a non-empty finite vector or coordinate matrix"
                )
            if name in {
                "parameter_steps",
                "x_steps",
                "parameter_scales",
                "x_scales",
            } and np.any(array <= 0.0):
                raise ValueError(f"{name} must contain strictly positive values")
            setattr(self, name, array.copy())

    def _specific_as_dict(self) -> dict[str, Any]:
        """Return ODR-specific backend-neutral controls."""
        return {
            "difference_scheme": ODRDifferenceScheme(self.difference_scheme).value,
            "ndigit": self.ndigit,
            "trust_region_factor": self.trust_region_factor,
            "initial_x_corrections": _optional_array(self.initial_x_corrections),
            "parameter_steps": _optional_array(self.parameter_steps),
            "x_steps": _optional_array(self.x_steps),
            "parameter_scales": _optional_array(self.parameter_scales),
            "x_scales": _optional_array(self.x_scales),
        }


# Concise public spelling used by workflow APIs and documentation.
ODROptions = OrthogonalDistanceOptions

# Public type alias for any concrete solver-options contract.
SolverOptions: TypeAlias = (
    OLSOptions | WLSOptions | EffectiveVarianceOptions | OrthogonalDistanceOptions
)


def default_solver_options(method: FitMethod | str) -> SolverOptions:
    """Return default typed options for a fitting method.

    Parameters
    ----------
    method : FitMethod or str
        Requested statistical strategy.

    Returns
    -------
    SolverOptions
        Fresh method-specific options object.
    """
    resolved = FitMethod(method)
    if resolved is FitMethod.OLS:
        return OLSOptions()
    if resolved is FitMethod.WLS:
        return WLSOptions()
    if resolved is FitMethod.EFFECTIVE_VARIANCE:
        return EffectiveVarianceOptions()
    if resolved is FitMethod.ODR:
        return OrthogonalDistanceOptions()
    raise ValueError(f"unsupported fitting method: {resolved.value}")


def _optional_array(value: np.ndarray | None) -> list[Any] | None:
    """Convert an optional NumPy vector or matrix to nested lists."""
    if value is None:
        return None
    return np.asarray(value, dtype=np.float64).tolist()


def _resolve_covariance_scaling(
    covariance_scaling: CovarianceScaling | str | None,
    absolute_sigma: bool | None,
) -> CovarianceScaling:
    """Resolve the modern covariance policy and legacy boolean alias.

    Parameters
    ----------
    covariance_scaling : CovarianceScaling, str, or None
        Explicit modern policy.
    absolute_sigma : bool or None
        Legacy alias.

    Returns
    -------
    CovarianceScaling
        Normalized covariance policy.

    Raises
    ------
    ValueError
        If both controls are supplied with conflicting meanings.
    """
    policy = (
        None if covariance_scaling is None else CovarianceScaling(covariance_scaling)
    )
    if absolute_sigma is None:
        return policy or CovarianceScaling.ABSOLUTE
    alias = (
        CovarianceScaling.ABSOLUTE
        if bool(absolute_sigma)
        else CovarianceScaling.REDUCED_CHI_SQUARE
    )
    if policy is not None and policy is not alias:
        raise ValueError(
            "absolute_sigma conflicts with the requested covariance_scaling policy"
        )
    return alias
