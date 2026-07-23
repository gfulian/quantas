# -*- coding: utf-8 -*-

"""Fit request and result contracts for EOS workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, cast

import numpy as np

from quantas.core.math.fitting import (
    CovarianceScaling,
    FitMethod,
    FitOptions as SolverFitOptionsBase,
    FitResult,
    OLSOptions,
    ParameterDefinition,
    SolverOptions,
    default_solver_options,
)
from quantas.core.physics.eos import (
    EOSModel,
    PVTModel,
    TemperatureEOSModel,
    parse_eos_model,
    parse_temperature_eos_model,
)

from ._model_utils import as_float64_vector
from .dataset_models import EOS_TARGET_NAMES


class EOSFitDomain(str, Enum):
    """Scientific domain addressed by one EOS fitting request."""

    PRESSURE_VOLUME = "pv"
    ENERGY_VOLUME = "ev"
    VOLUME_TEMPERATURE = "vt"
    PRESSURE_VOLUME_TEMPERATURE = "pvt"


@dataclass(frozen=True, slots=True)
class ParameterConstraint(ParameterDefinition):
    """EOS-facing parameter constraint using the general fitting contract."""


@dataclass(slots=True)
class EOSFitOptions:
    """Configure one EOS fit independently of its frontend.

    The EOS workflow owns only workflow-level controls. Numerical settings are
    carried by one typed ``solver_options`` object, so every regression method
    exposes an explicit, validated contract rather than an untyped mapping.

    Parameters
    ----------
    solver_options : SolverOptions or None, optional
        Typed numerical configuration. Use
        :class:`~quantas.api.eos.OLSOptions`,
        :class:`~quantas.api.eos.WLSOptions`,
        :class:`~quantas.api.eos.EffectiveVarianceOptions`, or
        :class:`~quantas.api.eos.ODROptions`. If omitted, defaults are generated
        from ``method``; when both are omitted, OLS is selected.
    method : FitMethod, str, or None, optional
        Convenience selector and consistency guard. Internally Quantas always
        stores typed ``solver_options``. When both arguments are supplied, the
        methods must agree.
    allow_extrapolation : bool, optional
        Whether downstream prediction requests may leave the observed domain.
    metadata : dict, optional
        Passive workflow provenance. It is merged into the solver metadata
        without mutating the caller's options object.

    Notes
    -----
    For weighted EOS methods, an unspecified covariance policy is normalized
    to :attr:`~quantas.api.eos.CovarianceScaling.INFLATE_ONLY`,
    reproducing the EosFit7 uncertainty convention. A policy explicitly set
    on ``solver_options`` is always preserved.

    Raises
    ------
    TypeError
        If ``solver_options`` is not a typed Quantas fitting-options object.
    ValueError
        If ``method`` disagrees with the typed solver-options method.
    """

    solver_options: SolverOptions | None = None
    method: FitMethod | str | None = None
    allow_extrapolation: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Resolve and detach the typed numerical solver configuration."""
        requested = None if self.method is None else FitMethod(self.method)
        solver: SolverOptions
        if self.solver_options is None:
            solver = default_solver_options(requested or FitMethod.OLS)
        else:
            if not isinstance(self.solver_options, SolverFitOptionsBase):
                raise TypeError(
                    "solver_options must be a typed Quantas fitting-options object"
                )
            solver = cast(SolverOptions, self.solver_options.with_metadata({}))
        if requested is not None and requested is not solver.method:
            raise ValueError(
                "EOSFitOptions method conflicts with solver_options: "
                f"{requested.value!r} != {solver.method.value!r}"
            )
        if (
            solver.method is not FitMethod.OLS
            and not solver.covariance_scaling_explicit
        ):
            solver.covariance_scaling = CovarianceScaling.INFLATE_ONLY
        self.solver_options = solver
        self.method = solver.method
        self.allow_extrapolation = bool(self.allow_extrapolation)
        self.metadata = dict(self.metadata)

    def to_fit_options(
        self,
        extra_metadata: dict[str, Any] | None = None,
    ) -> SolverOptions:
        """Return detached mathematical options with merged provenance.

        Parameters
        ----------
        extra_metadata : dict or None, optional
            Workflow metadata appended after :attr:`metadata`.

        Returns
        -------
        SolverOptions
            Concrete solver-options object ready for the numerical service.
        """
        solver = self.solver_options
        if solver is None:  # defensive guard after dataclass validation
            solver = OLSOptions()
        merged = {**self.metadata, **dict(extra_metadata or {})}
        return cast(SolverOptions, solver.with_metadata(merged))

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready workflow options mapping.

        Returns
        -------
        dict
            Workflow controls with nested typed solver settings.
        """
        solver = self.solver_options
        if solver is None:  # defensive guard after dataclass validation
            solver = OLSOptions()
        return {
            "method": solver.method.value,
            "solver_options": solver.as_dict(),
            "allow_extrapolation": self.allow_extrapolation,
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class EOSFitRequest:
    """Declarative request for one EOS fit.

    Parameters
    ----------
    model : EOSModel, TemperatureEOSModel, or str
        Pressure/energy EOS or volume-temperature model, selected according
        to ``domain``.
    target : str, optional
        Structural or energetic quantity selected from an :class:`EOSDataset`.
    domain : EOSFitDomain, optional
        Scientific relationship being fitted.
    constraints : tuple of ParameterConstraint, optional
        User overrides for initial values, fixed values, and bounds. Empty
        tuples request model defaults from the future EOS fitting service.
    options : EOSFitOptions, optional
        Explicit statistical and optimizer controls.
    mask : ndarray or None, optional
        Observation selection. It is validated against the dataset only when
        the request is executed.
    request_id : str or None, optional
        Stable batch or session identifier.
    metadata : dict, optional
        Passive provenance.

    Raises
    ------
    ValueError
        If the model/domain/target combination or constraints are invalid.
    """

    model: EOSModel | TemperatureEOSModel | PVTModel | str
    target: str = "volume"
    domain: EOSFitDomain = EOSFitDomain.PRESSURE_VOLUME
    constraints: tuple[ParameterConstraint, ...] = ()
    options: EOSFitOptions = field(default_factory=EOSFitOptions)
    mask: np.ndarray | None = None
    request_id: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize and validate one frontend-neutral request."""
        self.domain = EOSFitDomain(self.domain)
        if self.domain is EOSFitDomain.VOLUME_TEMPERATURE:
            if isinstance(self.model, (EOSModel, PVTModel)):
                raise ValueError("V-T requests require a temperature EOS model")
            self.model = parse_temperature_eos_model(self.model)
        elif self.domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            if not isinstance(self.model, PVTModel):
                raise ValueError("P-V-T requests require a compositional PVTModel")
        else:
            if isinstance(self.model, (TemperatureEOSModel, PVTModel)):
                raise ValueError("P-V and E-V requests require a pressure EOS model")
            self.model = parse_eos_model(self.model)
        self.target = str(self.target).lower()
        if self.target not in EOS_TARGET_NAMES:
            raise ValueError(f"unsupported EOS fitting target: {self.target}")
        if self.domain is EOSFitDomain.PRESSURE_VOLUME and self.target not in {
            "volume",
            "a",
            "b",
            "c",
        }:
            raise ValueError(
                "P-V EOS requests support volume or linear cell parameters only"
            )
        if self.domain is EOSFitDomain.VOLUME_TEMPERATURE and self.target not in {
            "volume",
            "a",
            "b",
            "c",
        }:
            raise ValueError(
                "V-T EOS requests support volume or linear cell parameters only"
            )
        if self.domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            if self.target != "volume":
                raise ValueError("P-V-T requests currently require target='volume'")
        if self.domain is EOSFitDomain.ENERGY_VOLUME:
            if self.target != "energy":
                raise ValueError("E-V requests require target='energy'")
            if not isinstance(self.model, EOSModel) or not self.model.supports_energy:
                raise ValueError(f"{self.model.tag} has no integrated E-V form")
        constraints = tuple(self.constraints)
        names = [constraint.name for constraint in constraints]
        if len(set(names)) != len(names):
            raise ValueError("EOS parameter constraints must have unique names")
        self.constraints = constraints
        if not isinstance(self.options, EOSFitOptions):
            raise TypeError("options must be an EOSFitOptions instance")
        if self.mask is not None:
            mask = np.asarray(self.mask, dtype=np.bool_)
            if mask.ndim != 1 or mask.size == 0 or not np.any(mask):
                raise ValueError("EOS request mask must be a non-empty boolean vector")
            self.mask = mask.copy()
        self.request_id = None if self.request_id is None else str(self.request_id)
        self.metadata = dict(self.metadata)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready fit request."""
        return {
            "request_id": self.request_id,
            "domain": self.domain.value,
            "target": self.target,
            "model": (
                {
                    "family": self.model.family.value,
                    "variant": (
                        None if self.model.variant is None else self.model.variant.value
                    ),
                    "tag": self.model.tag,
                }
                if isinstance(self.model, TemperatureEOSModel)
                else (
                    self.model.as_dict()
                    if isinstance(self.model, PVTModel)
                    else parse_eos_model(self.model).as_dict()
                )
            ),
            "constraints": [item.as_dict() for item in self.constraints],
            "options": self.options.as_dict(),
            "mask": None if self.mask is None else self.mask.tolist(),
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class EOSFitResult:
    """EOS-domain result wrapping the general numerical fit result."""

    request: EOSFitRequest
    fit: FitResult
    predictions: dict[str, np.ndarray] = field(default_factory=dict)
    derived: dict[str, float] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize prediction arrays and passive result mappings."""
        normalized: dict[str, np.ndarray] = {}
        for name, values in self.predictions.items():
            normalized[str(name)] = as_float64_vector(values, name=str(name))
        self.predictions = normalized
        self.derived = {str(name): float(value) for name, value in self.derived.items()}
        if not all(np.isfinite(value) for value in self.derived.values()):
            raise ValueError("derived EOS result values must be finite")
        self.warnings = [str(value) for value in self.warnings]
        self.metadata = dict(self.metadata)

    @property
    def parameter_values(self) -> dict[str, float]:
        """Return complete physical fit parameters keyed by stable name.

        Returns
        -------
        dict
            Empty for unsuccessful fits; otherwise the complete FREE, FIXED,
            and IMPLIED parameter values in reporting order.
        """
        if self.fit.parameters is None:
            return {}
        return {
            name: float(value)
            for name, value in zip(
                self.fit.parameter_names, self.fit.parameters, strict=True
            )
        }

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready EOS result."""
        return {
            "request": self.request.as_dict(),
            "fit": self.fit.as_dict(),
            "predictions": {
                name: values.tolist() for name, values in self.predictions.items()
            },
            "derived": dict(self.derived),
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }

__all__ = [
    "EOSFitDomain",
    "EOSFitOptions",
    "EOSFitRequest",
    "EOSFitResult",
    "ParameterConstraint",
]
