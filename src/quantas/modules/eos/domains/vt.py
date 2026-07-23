# -*- coding: utf-8 -*-

"""Volume-temperature fit models and parameter preparation for EOS.

The objects in this module adapt the analytical :class:`TemperatureEOS` core
models to the frontend-neutral fitting contracts.  They contain no CLI,
rendering, persistence, or terminal dependencies.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from numpy.polynomial import Polynomial

from quantas.core.math.fitting import (
    BaseFitModel,
    ParameterDefinition,
    ParameterMap,
    ParameterState,
)
from quantas.core.physics.eos import (
    TemperatureEOS,
    TemperatureEOSFamily,
    TemperatureEOSModel,
    TemperatureEOSVariant,
    parse_temperature_eos_model,
)

from ..models import ParameterConstraint

_POSITIVE_LOWER_BOUND = float(np.nextafter(0.0, 1.0))
_DEFAULT_REFERENCE_TEMPERATURE = 298.15
_DEFAULT_EINSTEIN_TEMPERATURE = 500.0
_DEFAULT_KP = 4.0
_DEFAULT_SATURATION_TEMPERATURE = 300.0


class TemperatureEOSFitModel(BaseFitModel):
    r"""Adapt an analytical V--T model to the general fitting contract.

    Parameters
    ----------
    model : TemperatureEOSModel or str
        Thermal-expansion family and variant.
    target : {"volume", "a", "b", "c"}
        Structural quantity represented by the model.
    axial : bool, optional
        When ``True`` the dependent solver coordinate is the auxiliary
        quantity :math:`q=x^3`, not the original cell-axis length.
    initial_parameters : mapping, optional
        Complete initial values in :attr:`parameter_names` order.  The normal
        EOS workflow uses a :class:`~quantas.core.math.fitting.ParameterMap`
        and therefore obtains initial free values from that mapping.

    Notes
    -----
    For an axial fit, ``L0`` is the public physical reference length while
    evaluation uses :math:`q_0=L_0^3`. Other coefficients remain the exact
    parameters of the auxiliary :math:`q(T)=x(T)^3` equation. The physical
    linear expansion coefficient is :math:`\alpha_x=\alpha_q/3`; coefficients
    are not silently divided by three when such a conversion would be
    model-dependent. Because ``L0`` is optimized directly, its covariance is
    already reported in the physical linear parameter space.
    """

    def __init__(
        self,
        model: TemperatureEOSModel | str,
        target: str,
        *,
        axial: bool = False,
        initial_parameters: Mapping[str, float] | None = None,
    ) -> None:
        normalized_target = str(target).lower()
        if normalized_target not in {"volume", "a", "b", "c"}:
            raise ValueError("V-T target must be volume, a, b, or c")
        if axial and normalized_target == "volume":
            raise ValueError("axial V-T models require target a, b, or c")
        if not axial and normalized_target != "volume":
            raise ValueError("non-axial V-T models require target volume")
        self.temperature_model = parse_temperature_eos_model(model)
        self.target = normalized_target
        self.axial = bool(axial)
        self._temperature_eos = TemperatureEOS()
        self._parameter_names = temperature_parameter_names(
            self.temperature_model,
            axial=self.axial,
        )
        defaults = _default_parameter_values(
            self.temperature_model,
            axial=self.axial,
        )
        if initial_parameters is not None:
            defaults.update(
                {
                    _normalize_temperature_parameter_name(
                        name,
                        axial=self.axial,
                    ): float(value)
                    for name, value in initial_parameters.items()
                }
            )
        self._initial = np.asarray(
            [defaults[name] for name in self._parameter_names], dtype=np.float64
        )

    @property
    def name(self) -> str:
        """Return a stable physical-model identifier."""
        prefix = "axial_temperature_eos" if self.axial else "temperature_eos"
        return f"{prefix}:{self.target}:{self.temperature_model.tag}"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return complete thermal-model parameter order."""
        return self._parameter_names

    def evaluate(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Evaluate the structural quantity at supplied temperatures."""
        mapping = self.parameter_mapping(parameters)
        return self._temperature_eos.value(self.temperature_model, mapping, x)

    def derivative_x(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        r"""Return the analytical derivative :math:`\partial X/\partial T`."""
        mapping = self.parameter_mapping(parameters)
        return self._temperature_eos.derivative(self.temperature_model, mapping, x)

    def expansion_coefficient(
        self,
        temperature: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return the exact expansion coefficient of the solver quantity."""
        mapping = self.parameter_mapping(parameters)
        return self._temperature_eos.expansion_coefficient(
            self.temperature_model, mapping, temperature
        )

    def linear_expansion_coefficient(
        self,
        temperature: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return physical linear expansion for an axial :math:`q=x^3` fit.

        Raises
        ------
        ValueError
            If called for a volumetric model.
        """
        if not self.axial:
            raise ValueError("linear expansion is defined only for axial V-T fits")
        return TemperatureEOS.linear_expansion_coefficient(
            self.expansion_coefficient(temperature, parameters)
        )

    def parameter_mapping(
        self,
        parameters: np.ndarray | Sequence[float],
    ) -> dict[str, float]:
        """Return a core-compatible parameter mapping."""
        values = np.asarray(parameters, dtype=np.float64)
        if values.ndim != 1 or values.size != len(self._parameter_names):
            raise ValueError(
                "temperature EOS model received an incompatible parameter vector"
            )
        public = {
            name: float(value)
            for name, value in zip(self._parameter_names, values, strict=True)
        }
        reference_name = "L0" if self.axial else "V0"
        reference = public.pop(reference_name)
        return {
            "V0": reference**3 if self.axial else reference,
            **public,
        }

    def initial_guess(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return configured complete initial values after validation."""
        _validate_temperature_data(x, y)
        return self._initial.copy()

    def bounds(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return physical default bounds in complete parameter order."""
        _validate_temperature_data(x, y)
        lower = np.asarray(
            [_default_parameter_bounds(name)[0] for name in self._parameter_names],
            dtype=np.float64,
        )
        upper = np.asarray(
            [_default_parameter_bounds(name)[1] for name in self._parameter_names],
            dtype=np.float64,
        )
        return lower, upper

    def metadata(self) -> dict[str, Any]:
        """Return thermal model and coordinate metadata."""
        return {
            **super().metadata(),
            "temperature_eos_model": {
                "family": self.temperature_model.family.value,
                "variant": (
                    None
                    if self.temperature_model.variant is None
                    else self.temperature_model.variant.value
                ),
                "tag": self.temperature_model.tag,
            },
            "relationship": f"{self.target}(temperature)",
            "linear_eos": self.axial,
            "linear_target": self.target if self.axial else None,
            "coordinate_transform": "q=x^3" if self.axial else None,
            "reference_parameter": "L0" if self.axial else "V0",
            "reference_mapping": "V0=L0^3" if self.axial else None,
            "coefficient_space": (
                "auxiliary_cubed_length" if self.axial else "structural_quantity"
            ),
        }


def temperature_parameter_names(
    model: TemperatureEOSModel | str,
    *,
    axial: bool = False,
) -> tuple[str, ...]:
    """Return complete reporting parameters for one thermal model.

    ``V0`` is used for a volume target and physical ``L0`` for an axial
    target. The sequence also includes fixed reference quantities and
    variant-implied coefficients so results preserve FREE/FIXED/IMPLIED
    provenance.
    """
    spec = parse_temperature_eos_model(model)
    reference = "L0" if axial else "V0"
    if spec.family is TemperatureEOSFamily.BERMAN:
        return (reference, "temperature_ref", "alpha0", "alpha1")
    if spec.family is TemperatureEOSFamily.FEI:
        return (reference, "temperature_ref", "alpha0", "alpha1", "alpha2")
    if spec.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL:
        return (reference, "temperature_ref", "alpha0", "alpha1")
    if spec.family is TemperatureEOSFamily.SALJE:
        return (reference, "temperature_ref", "p1", "theta_sat")
    return (reference, "temperature_ref", "alpha_ref", "theta_e", "kp")


def estimate_temperature_parameters(
    model: TemperatureEOSModel | str,
    temperature: np.ndarray | Sequence[float],
    value: np.ndarray | Sequence[float],
    *,
    reference_temperature: float | None = None,
    axial: bool = False,
) -> dict[str, float]:
    """Estimate complete thermal-model parameters from V--T observations.

    Parameters
    ----------
    model : TemperatureEOSModel or str
        Thermal model and variant.
    temperature, value : array-like
        Temperatures in kelvin and positive structural quantities. For axial
        fitting, ``value`` must already be the auxiliary coordinate ``x**3``;
        the returned reference estimate is the physical ``L0``.
    reference_temperature : float or None, optional
        Fixed reference temperature.  The default is 298.15 K, except for
        Salje which is intrinsically referenced to 0 K.

    Returns
    -------
    dict
        Complete initial values in :func:`temperature_parameter_names` order.
    """
    spec = parse_temperature_eos_model(model)
    t, y = _validate_temperature_data(temperature, value)
    tref = _reference_temperature(spec, reference_temperature)
    if spec.family is TemperatureEOSFamily.SALJE:
        internal = _estimate_salje(t, y)
        reference = float(internal.pop("V0"))
        return {
            "L0" if axial else "V0": (
                float(np.cbrt(reference)) if axial else reference
            ),
            **internal,
        }

    value_ref = _estimate_reference_value(t, y, tref)
    if spec.family is TemperatureEOSFamily.BERMAN:
        estimate = _estimate_berman(spec, t, y, value_ref, tref)
    elif spec.family is TemperatureEOSFamily.FEI:
        estimate = _estimate_fei(spec, t, y, value_ref, tref)
    elif spec.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL:
        estimate = _estimate_modified_holland_powell(spec, t, y, value_ref, tref)
    else:
        estimate = _estimate_khp(t, y, value_ref, tref)
    return {
        "L0" if axial else "V0": (float(np.cbrt(value_ref)) if axial else value_ref),
        "temperature_ref": tref,
        **estimate,
    }


def build_temperature_parameter_map(
    model: TemperatureEOSModel | str,
    temperature: np.ndarray | Sequence[float],
    value: np.ndarray | Sequence[float],
    constraints: Sequence[ParameterConstraint] = (),
    *,
    value_unit: str = "angstrom^3",
    reference_temperature: float | None = None,
    axial: bool = False,
) -> ParameterMap:
    """Build a reduced/full parameter mapping for one V--T fit.

    Parameters
    ----------
    model : TemperatureEOSModel or str
        Thermal model and variant.
    temperature, value : array-like
        Selected observations.  ``value`` is the measured volume or the
        auxiliary cubed length for an axial fit.
    constraints : sequence of ParameterConstraint, optional
        User overrides for initial values, fixed values, and bounds.
    value_unit : str, optional
        Unit of the public reference parameter: volume for a volumetric fit or
        length for an axial fit.
    reference_temperature : float or None, optional
        Default fixed reference temperature when no explicit constraint is
        supplied.
    axial : bool, optional
        Whether coefficients belong to the auxiliary cubed-length equation.

    Returns
    -------
    ParameterMap
        Complete FREE/FIXED/IMPLIED parameter contract.
    """
    spec = parse_temperature_eos_model(model)
    overrides = _temperature_constraint_overrides(constraints, axial=axial)
    allowed_names = set(temperature_parameter_names(spec, axial=axial))
    unknown = sorted(set(overrides) - allowed_names)
    if unknown:
        raise ValueError(
            "unknown temperature-EOS parameter constraint(s): " + ", ".join(unknown)
        )
    requested_tref = reference_temperature
    tref_override = overrides.get("temperature_ref")
    if tref_override is not None:
        if tref_override.state is not ParameterState.FIXED:
            raise ValueError(
                "temperature_ref is a reference condition and must be fixed"
            )
        if tref_override.value is None:
            raise ValueError("fixed temperature_ref requires a value")
        requested_tref = float(tref_override.value)
    estimates = estimate_temperature_parameters(
        spec,
        temperature,
        value,
        reference_temperature=requested_tref,
        axial=axial,
    )
    definitions = tuple(
        _temperature_parameter_definition(
            spec,
            name,
            estimates,
            overrides.get(name),
            value_unit=value_unit,
            axial=axial,
        )
        for name in temperature_parameter_names(spec, axial=axial)
    )
    return ParameterMap(definitions, resolver=_temperature_resolver(spec))


def _temperature_parameter_definition(
    model: TemperatureEOSModel,
    name: str,
    estimates: Mapping[str, float],
    override: ParameterConstraint | None,
    *,
    value_unit: str,
    axial: bool,
) -> ParameterDefinition:
    """Build one fitted, fixed, or implied thermal parameter."""
    unit = _temperature_parameter_unit(
        model,
        name,
        value_unit,
        axial=axial,
    )
    description = _temperature_parameter_description(name, axial=axial)
    implied = _is_implied_temperature_parameter(model, name)
    if implied:
        if override is not None and override.state is not ParameterState.IMPLIED:
            raise ValueError(
                f"parameter {name} is implied by {model.tag} and cannot be "
                f"declared {override.state.value}"
            )
        return ParameterDefinition.implied(
            name,
            value=None if override is None else override.value,
            lower_bound=-np.inf if override is None else override.lower_bound,
            upper_bound=np.inf if override is None else override.upper_bound,
            unit=unit if override is None else override.unit or unit,
            description=description,
            metadata={"source": "temperature_model_variant"},
        )

    default_fixed = _is_default_fixed_temperature_parameter(model, name)
    if override is None and default_fixed:
        return ParameterDefinition.fixed(
            name,
            estimates[name],
            lower_bound=_default_parameter_bounds(name)[0],
            upper_bound=_default_parameter_bounds(name)[1],
            unit=unit,
            description=description,
            metadata={"source": "model_default"},
        )
    if override is None:
        lower, upper = _default_parameter_bounds(name)
        return ParameterDefinition.free(
            name,
            estimates[name],
            lower_bound=lower,
            upper_bound=upper,
            unit=unit,
            description=description,
            metadata={"initial_source": "analytical_vt_estimate"},
        )
    if override.state is ParameterState.FREE:
        initial = (
            estimates[name]
            if override.initial_value is None
            else override.initial_value
        )
        lower, upper = _default_parameter_bounds(name)
        return ParameterDefinition.free(
            name,
            initial,
            lower_bound=(
                lower if np.isneginf(override.lower_bound) else override.lower_bound
            ),
            upper_bound=(
                upper if np.isposinf(override.upper_bound) else override.upper_bound
            ),
            unit=override.unit or unit,
            description=override.description or description,
            metadata={**override.metadata, "initial_source": "user"},
        )
    if override.state is ParameterState.FIXED:
        if override.value is None:
            raise ValueError(f"fixed parameter {name} requires a value")
        lower, upper = _default_parameter_bounds(name)
        return ParameterDefinition.fixed(
            name,
            override.value,
            lower_bound=(
                lower if np.isneginf(override.lower_bound) else override.lower_bound
            ),
            upper_bound=(
                upper if np.isposinf(override.upper_bound) else override.upper_bound
            ),
            unit=override.unit or unit,
            description=override.description or description,
            metadata={**override.metadata, "source": "user_fixed"},
        )
    raise ValueError(f"parameter {name} cannot be declared {override.state.value}")


def _temperature_resolver(model: TemperatureEOSModel):
    """Return a resolver for variant-implied thermal coefficients."""

    def resolver(values: Mapping[str, float]) -> Mapping[str, float]:
        resolved: dict[str, float] = {}
        if (
            model.family is TemperatureEOSFamily.BERMAN
            and model.variant is TemperatureEOSVariant.LINEAR
        ):
            resolved["alpha1"] = 0.0
        elif (
            model.family is TemperatureEOSFamily.FEI
            and model.variant is TemperatureEOSVariant.LINEAR
        ):
            resolved["alpha2"] = 0.0
        elif (
            model.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL
            and model.variant is TemperatureEOSVariant.SIMPLIFIED
        ):
            resolved["alpha1"] = 0.0
        return resolved

    return resolver


def _is_implied_temperature_parameter(
    model: TemperatureEOSModel,
    name: str,
) -> bool:
    return (
        (
            name == "alpha1"
            and model.family is TemperatureEOSFamily.BERMAN
            and model.variant is TemperatureEOSVariant.LINEAR
        )
        or (
            name == "alpha2"
            and model.family is TemperatureEOSFamily.FEI
            and model.variant is TemperatureEOSVariant.LINEAR
        )
        or (
            name == "alpha1"
            and model.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL
            and model.variant is TemperatureEOSVariant.SIMPLIFIED
        )
    )


def _is_default_fixed_temperature_parameter(
    model: TemperatureEOSModel,
    name: str,
) -> bool:
    if name == "temperature_ref":
        return True
    return model.family is TemperatureEOSFamily.KROLL_HOLLAND_POWELL and name in {
        "theta_e",
        "kp",
    }


def _temperature_constraint_overrides(
    constraints: Sequence[ParameterConstraint],
    *,
    axial: bool,
) -> dict[str, ParameterConstraint]:
    """Normalize thermal parameter aliases and reject duplicates."""
    overrides: dict[str, ParameterConstraint] = {}
    for constraint in constraints:
        name = _normalize_temperature_parameter_name(
            constraint.name,
            axial=axial,
        )
        if name in overrides:
            raise ValueError(f"duplicate temperature-EOS parameter constraint: {name}")
        overrides[name] = (
            constraint
            if constraint.name == name
            else ParameterConstraint(
                name=name,
                state=constraint.state,
                initial_value=constraint.initial_value,
                value=constraint.value,
                lower_bound=constraint.lower_bound,
                upper_bound=constraint.upper_bound,
                unit=constraint.unit,
                description=constraint.description,
                metadata=constraint.metadata,
            )
        )
    return overrides


def _normalize_temperature_parameter_name(name: str, *, axial: bool) -> str:
    key = str(name).strip().lower().replace("-", "_")
    aliases = {
        "v0": "V0",
        "l0": "L0",
        "tref": "temperature_ref",
        "t_ref": "temperature_ref",
        "alpharef": "alpha_ref",
        "alpha_ref": "alpha_ref",
        "thetae": "theta_e",
        "theta_e": "theta_e",
        "thetasat": "theta_sat",
        "theta_sat": "theta_sat",
        "kprime": "kp",
        "kp": "kp",
    }
    normalized = aliases.get(key, key)
    if normalized == "V0" and axial:
        raise ValueError("axial V-T fits use L0, not V0")
    if normalized == "L0" and not axial:
        raise ValueError("volumetric V-T fits use V0, not L0")
    if key in {"value_ref", "valueref", "vref", "xref", "qref", "x0"}:
        raise ValueError("V-T constraints use V0 for volume or L0 for an axis")
    return normalized


def _default_parameter_values(
    model: TemperatureEOSModel,
    *,
    axial: bool,
) -> dict[str, float]:
    values = {
        "L0" if axial else "V0": 1.0,
        "temperature_ref": 0.0
        if model.family is TemperatureEOSFamily.SALJE
        else _DEFAULT_REFERENCE_TEMPERATURE,
        "alpha0": 1.0e-5,
        "alpha1": 0.0,
        "alpha2": 0.0,
        "alpha_ref": 1.0e-5,
        "p1": 1.0e-5,
        "theta_sat": _DEFAULT_SATURATION_TEMPERATURE,
        "theta_e": _DEFAULT_EINSTEIN_TEMPERATURE,
        "kp": _DEFAULT_KP,
    }
    return values


def _default_parameter_bounds(name: str) -> tuple[float, float]:
    if name in {"V0", "L0"}:
        return _POSITIVE_LOWER_BOUND, np.inf
    if name == "temperature_ref":
        return 0.0, np.inf
    if name in {"theta_sat", "theta_e", "kp"}:
        return _POSITIVE_LOWER_BOUND, np.inf
    return -np.inf, np.inf


def _temperature_parameter_unit(
    model: TemperatureEOSModel,
    name: str,
    value_unit: str,
    *,
    axial: bool,
) -> str:
    if name == "alpha1":
        return (
            "K^-1/2"
            if model.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL
            else "K^-2"
        )
    return {
        "V0": str(value_unit),
        "L0": str(value_unit),
        "temperature_ref": "K",
        "alpha0": "K^-1",
        "alpha2": "K",
        "alpha_ref": "K^-1",
        "p1": (f"{value_unit} K^-1" if axial else f"({value_unit})^(1/3) K^-1"),
        "theta_sat": "K",
        "theta_e": "K",
        "kp": "1",
    }[name]


def _temperature_parameter_description(name: str, *, axial: bool) -> str:
    return {
        "V0": "reference volume",
        "L0": "reference cell-axis length",
        "temperature_ref": "fixed reference temperature",
        "alpha0": "model alpha0 coefficient",
        "alpha1": "model alpha1 coefficient",
        "alpha2": "Fei inverse-square coefficient",
        "alpha_ref": "exact expansion coefficient at the reference temperature",
        "p1": "Salje expansion parameter",
        "theta_sat": "Salje saturation temperature",
        "theta_e": "Einstein temperature",
        "kp": "reference pressure derivative of the bulk modulus",
    }[name]


def _reference_temperature(
    model: TemperatureEOSModel,
    requested: float | None,
) -> float:
    if model.family is TemperatureEOSFamily.SALJE:
        if requested is not None and float(requested) != 0.0:
            raise ValueError("Salje requires temperature_ref = 0 K")
        return 0.0
    value = _DEFAULT_REFERENCE_TEMPERATURE if requested is None else float(requested)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("reference temperature must be finite and positive")
    return value


def _estimate_reference_value(t: np.ndarray, y: np.ndarray, tref: float) -> float:
    degree = min(3, int(t.size) - 1)
    polynomial = Polynomial.fit(t, y, degree)
    value = float(polynomial(tref))
    if not np.isfinite(value) or value <= 0.0:
        value = float(y[int(np.argmin(np.abs(t - tref)))])
    if value <= 0.0 or not np.isfinite(value):
        raise ValueError("V-T data did not produce a positive reference quantity")
    return value


def _estimate_berman(
    model: TemperatureEOSModel,
    t: np.ndarray,
    y: np.ndarray,
    value_ref: float,
    tref: float,
) -> dict[str, float]:
    dt = t - tref
    response = y / value_ref - 1.0
    columns = [dt]
    if model.variant is TemperatureEOSVariant.QUADRATIC:
        columns.append(0.5 * dt**2)
    coefficients = _linear_coefficients(columns, response)
    alpha0 = float(coefficients[0])
    alpha1 = 0.0 if len(coefficients) == 1 else float(coefficients[1])
    return {"alpha0": alpha0, "alpha1": alpha1}


def _estimate_fei(
    model: TemperatureEOSModel,
    t: np.ndarray,
    y: np.ndarray,
    value_ref: float,
    tref: float,
) -> dict[str, float]:
    if np.any(t <= 0.0):
        raise ValueError("Fei fitting requires temperatures above 0 K")
    response = np.log(y / value_ref)
    columns = [t - tref, 0.5 * (t**2 - tref**2)]
    if model.variant is TemperatureEOSVariant.INVERSE_SQUARE:
        columns.append(-(1.0 / t - 1.0 / tref))
    coefficients = _linear_coefficients(columns, response)
    return {
        "alpha0": float(coefficients[0]),
        "alpha1": float(coefficients[1]),
        "alpha2": 0.0 if len(coefficients) == 2 else float(coefficients[2]),
    }


def _estimate_modified_holland_powell(
    model: TemperatureEOSModel,
    t: np.ndarray,
    y: np.ndarray,
    value_ref: float,
    tref: float,
) -> dict[str, float]:
    if np.any(t <= 0.0):
        raise ValueError("modified Holland-Powell fitting requires T > 0 K")
    root_difference = np.sqrt(t) - np.sqrt(tref)
    response = y / value_ref - 1.0
    alpha0_basis = (t - tref) - 20.0 * root_difference
    columns = [alpha0_basis]
    if model.variant is TemperatureEOSVariant.GENERAL:
        columns.append(-2.0 * root_difference)
    coefficients = _linear_coefficients(columns, response)
    return {
        "alpha0": float(coefficients[0]),
        "alpha1": 0.0 if len(coefficients) == 1 else float(coefficients[1]),
    }


def _estimate_salje(t: np.ndarray, y: np.ndarray) -> dict[str, float]:
    positive = t > 0.0
    if not np.any(positive):
        raise ValueError("Salje fitting requires at least one temperature above 0 K")
    tp = t[positive]
    root = np.cbrt(y[positive])
    theta = float(np.clip(np.median(tp), 50.0, 500.0))
    basis = theta * (1.0 / np.tanh(theta / tp) - 1.0)
    design = np.column_stack((np.ones_like(basis), basis))
    coefficients, *_ = np.linalg.lstsq(design, root, rcond=None)
    root_ref = float(coefficients[0])
    p1 = float(coefficients[1])
    if not np.isfinite(root_ref) or root_ref <= 0.0:
        root_ref = float(np.min(root))
    return {
        "V0": root_ref**3,
        "temperature_ref": 0.0,
        "p1": p1,
        "theta_sat": theta,
    }


def _estimate_khp(
    t: np.ndarray,
    y: np.ndarray,
    value_ref: float,
    tref: float,
) -> dict[str, float]:
    degree = min(3, int(t.size) - 1)
    polynomial = Polynomial.fit(t, y, degree)
    derivative = float(polynomial.deriv()(tref))
    alpha_ref = derivative / value_ref
    if not np.isfinite(alpha_ref):
        alpha_ref = 0.0
    return {
        "alpha_ref": alpha_ref,
        "theta_e": _DEFAULT_EINSTEIN_TEMPERATURE,
        "kp": _DEFAULT_KP,
    }


def _linear_coefficients(
    columns: Sequence[np.ndarray],
    response: np.ndarray,
) -> np.ndarray:
    design = np.column_stack(columns)
    coefficients, *_ = np.linalg.lstsq(design, response, rcond=None)
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("V-T data did not produce finite initial coefficients")
    return np.asarray(coefficients, dtype=np.float64)


def _validate_temperature_data(
    temperature: np.ndarray | Sequence[float],
    value: np.ndarray | Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(temperature, dtype=np.float64)
    y = np.asarray(value, dtype=np.float64)
    if t.ndim != 1 or y.ndim != 1:
        raise ValueError("V-T arrays must be one-dimensional")
    if t.shape != y.shape or t.size < 2:
        raise ValueError("V-T arrays must have equal length and at least two points")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(y)):
        raise ValueError("V-T arrays must contain finite values")
    if np.any(t < 0.0):
        raise ValueError("temperature cannot be below absolute zero")
    if np.any(y <= 0.0):
        raise ValueError("V-T structural quantities must be strictly positive")
    return t.copy(), y.copy()


__all__ = [
    "TemperatureEOSFitModel",
    "build_temperature_parameter_map",
    "estimate_temperature_parameters",
    "temperature_parameter_names",
]
