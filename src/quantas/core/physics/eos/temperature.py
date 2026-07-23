# -*- coding: utf-8 -*-

"""Analytical volume-temperature and thermal-expansion models.

This module implements the Berman, Fei, modified Holland--Powell, Salje,
and Kroll--Holland--Powell formulations audited against Angel,
Gonzalez-Platas, and Alvaro (2014).  The implementation is frontend-neutral
and contains no fitting, persistence, plotting, or terminal dependencies.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
import math
import re
from typing import TypeAlias

import numpy as np

ArrayLike: TypeAlias = np.ndarray | float | Sequence[float]
ParameterLike: TypeAlias = "TemperatureEOSParameters | Mapping[str, float]"


class TemperatureEOSFamily(str, Enum):
    """Supported analytical volume-temperature model families."""

    BERMAN = "berman"
    FEI = "fei"
    MODIFIED_HOLLAND_POWELL = "modified-holland-powell"
    SALJE = "salje"
    KROLL_HOLLAND_POWELL = "kroll-holland-powell"


class TemperatureEOSVariant(str, Enum):
    """Supported variants of analytical volume-temperature models."""

    LINEAR = "linear"
    QUADRATIC = "quadratic"
    INVERSE_SQUARE = "inverse-square"
    SIMPLIFIED = "simplified"
    GENERAL = "general"
    STANDARD = "standard"


_ALLOWED_VARIANTS: dict[TemperatureEOSFamily, tuple[TemperatureEOSVariant, ...]] = {
    TemperatureEOSFamily.BERMAN: (
        TemperatureEOSVariant.LINEAR,
        TemperatureEOSVariant.QUADRATIC,
    ),
    TemperatureEOSFamily.FEI: (
        TemperatureEOSVariant.LINEAR,
        TemperatureEOSVariant.INVERSE_SQUARE,
    ),
    TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL: (
        TemperatureEOSVariant.SIMPLIFIED,
        TemperatureEOSVariant.GENERAL,
    ),
    TemperatureEOSFamily.SALJE: (TemperatureEOSVariant.STANDARD,),
    TemperatureEOSFamily.KROLL_HOLLAND_POWELL: (TemperatureEOSVariant.STANDARD,),
}

_DEFAULT_VARIANT = {
    TemperatureEOSFamily.BERMAN: TemperatureEOSVariant.QUADRATIC,
    TemperatureEOSFamily.FEI: TemperatureEOSVariant.INVERSE_SQUARE,
    TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL: TemperatureEOSVariant.GENERAL,
    TemperatureEOSFamily.SALJE: TemperatureEOSVariant.STANDARD,
    TemperatureEOSFamily.KROLL_HOLLAND_POWELL: TemperatureEOSVariant.STANDARD,
}

_ALIASES: dict[str, TemperatureEOSFamily] = {
    "berman": TemperatureEOSFamily.BERMAN,
    "b": TemperatureEOSFamily.BERMAN,
    "fei": TemperatureEOSFamily.FEI,
    "f": TemperatureEOSFamily.FEI,
    "modifiedhollandpowell": TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL,
    "modified-holland-powell": TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL,
    "mhp": TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL,
    "salje": TemperatureEOSFamily.SALJE,
    "s": TemperatureEOSFamily.SALJE,
    "krollhollandpowell": TemperatureEOSFamily.KROLL_HOLLAND_POWELL,
    "kroll-holland-powell": TemperatureEOSFamily.KROLL_HOLLAND_POWELL,
    "khp": TemperatureEOSFamily.KROLL_HOLLAND_POWELL,
}


@dataclass(frozen=True, slots=True)
class TemperatureEOSModel:
    """Describe one volume-temperature model family and variant.

    Parameters
    ----------
    family : TemperatureEOSFamily
        Analytical model family.
    variant : TemperatureEOSVariant, optional
        Family-specific variant.  The audited default is used when omitted.
    """

    family: TemperatureEOSFamily
    variant: TemperatureEOSVariant | None = None

    def __post_init__(self) -> None:
        """Resolve and validate the model variant."""
        variant = self.variant or _DEFAULT_VARIANT[self.family]
        if variant not in _ALLOWED_VARIANTS[self.family]:
            allowed = ", ".join(value.value for value in _ALLOWED_VARIANTS[self.family])
            raise ValueError(
                f"unsupported variant {variant.value!r} for {self.family.value}; "
                f"available variants are {allowed}"
            )
        object.__setattr__(self, "variant", variant)

    @property
    def tag(self) -> str:
        """Return a stable compact model tag."""
        prefixes = {
            TemperatureEOSFamily.BERMAN: "BERMAN",
            TemperatureEOSFamily.FEI: "FEI",
            TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL: "MHP",
            TemperatureEOSFamily.SALJE: "SALJE",
            TemperatureEOSFamily.KROLL_HOLLAND_POWELL: "KHP",
        }
        variant = self.variant
        assert variant is not None
        return f"{prefixes[self.family]}:{variant.value}"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return parameters required by the selected formulation."""
        if self.family is TemperatureEOSFamily.BERMAN:
            return (
                ("V0", "temperature_ref", "alpha0")
                if self.variant is TemperatureEOSVariant.LINEAR
                else ("V0", "temperature_ref", "alpha0", "alpha1")
            )
        if self.family is TemperatureEOSFamily.FEI:
            return (
                ("V0", "temperature_ref", "alpha0", "alpha1")
                if self.variant is TemperatureEOSVariant.LINEAR
                else ("V0", "temperature_ref", "alpha0", "alpha1", "alpha2")
            )
        if self.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL:
            return (
                ("V0", "temperature_ref", "alpha0")
                if self.variant is TemperatureEOSVariant.SIMPLIFIED
                else ("V0", "temperature_ref", "alpha0", "alpha1")
            )
        if self.family is TemperatureEOSFamily.SALJE:
            return ("V0", "p1", "theta_sat")
        return ("V0", "temperature_ref", "alpha_ref", "theta_e", "kp")


@dataclass(frozen=True, slots=True)
class TemperatureEOSParameters:
    """Physical parameters of one volume-temperature formulation.

    Parameters
    ----------
    V0 : float
        Reference volume-like quantity. It may be an absolute volume, a
        normalized volume ratio, or the auxiliary ``L0**3`` supplied by an
        axial fitting adapter.
    temperature_ref : float, optional
        Reference temperature in kelvin.  Salje is intrinsically referenced to
        0 K and ignores this field after validating it.
    alpha0, alpha1, alpha2 : float, optional
        Family-specific thermal-expansion coefficients.
    alpha_ref : float, optional
        Exact expansion coefficient at ``temperature_ref`` for KHP.
    p1 : float, optional
        Salje expansion parameter.
    theta_sat : float, optional
        Salje saturation temperature in kelvin.
    theta_e : float, optional
        Einstein temperature in kelvin for KHP.
    kp : float, optional
        Reference pressure derivative of the bulk modulus for KHP.
    """

    V0: float
    temperature_ref: float = 298.15
    alpha0: float = 0.0
    alpha1: float = 0.0
    alpha2: float = 0.0
    alpha_ref: float = 0.0
    p1: float = 0.0
    theta_sat: float = 0.0
    theta_e: float = 0.0
    kp: float = 4.0

    def __post_init__(self) -> None:
        """Validate universal physical conditions."""
        values = np.asarray(
            [
                self.V0,
                self.temperature_ref,
                self.alpha0,
                self.alpha1,
                self.alpha2,
                self.alpha_ref,
                self.p1,
                self.theta_sat,
                self.theta_e,
                self.kp,
            ],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("temperature-EOS parameters must be finite")
        if self.V0 <= 0.0:
            raise ValueError("V0 must be strictly positive")
        if self.temperature_ref < 0.0:
            raise ValueError("temperature_ref cannot be negative")


def parse_temperature_eos_model(
    model: str | TemperatureEOSFamily | TemperatureEOSModel,
    variant: str | TemperatureEOSVariant | None = None,
) -> TemperatureEOSModel:
    """Return a canonical volume-temperature model specification.

    Parameters
    ----------
    model : str, TemperatureEOSFamily, or TemperatureEOSModel
        Family name, alias, or existing specification.
    variant : str or TemperatureEOSVariant, optional
        Family-specific variant.

    Returns
    -------
    TemperatureEOSModel
        Canonical immutable model specification.

    Raises
    ------
    ValueError
        If the family or variant is unsupported.
    """
    if isinstance(model, TemperatureEOSModel):
        if variant is not None:
            raise ValueError("variant cannot be supplied with TemperatureEOSModel")
        return model
    embedded_variant: str | None = None
    family: TemperatureEOSFamily | None
    if isinstance(model, TemperatureEOSFamily):
        family = model
    else:
        model_text = str(model).strip().lower()
        if ":" in model_text:
            family_text, embedded_variant = model_text.split(":", 1)
            if variant is not None:
                raise ValueError(
                    "variant cannot be supplied both in the model tag and separately"
                )
            model_text = family_text
        family = _ALIASES.get(re.sub(r"[\s_]", "", model_text))
    if family is None:
        raise ValueError(f"unknown temperature EOS: {model!r}")
    parsed_variant: TemperatureEOSVariant | None
    variant_value = embedded_variant if embedded_variant is not None else variant
    if variant_value is None:
        parsed_variant = None
    elif isinstance(variant_value, TemperatureEOSVariant):
        parsed_variant = variant_value
    else:
        key = str(variant_value).strip().lower().replace("_", "-")
        try:
            parsed_variant = TemperatureEOSVariant(key)
        except ValueError as exc:
            raise ValueError(
                f"unknown temperature-EOS variant: {variant_value!r}"
            ) from exc
    return TemperatureEOSModel(family, parsed_variant)


def available_temperature_eos_models() -> tuple[TemperatureEOSModel, ...]:
    """Return all supported volume-temperature family/variant combinations."""
    return tuple(
        TemperatureEOSModel(family, variant)
        for family in TemperatureEOSFamily
        for variant in _ALLOWED_VARIANTS[family]
    )


class TemperatureEOS:
    r"""Evaluate analytical structural-quantity--temperature equations.

    For a structural quantity :math:`X(T)`, the exact thermal expansion
    coefficient returned by :meth:`expansion_coefficient` is

    .. math::

        \alpha(T)=\frac{1}{X(T)}\left(\frac{\partial X}{\partial T}\right)_P.

    ``X`` may represent a physical volume, an auxiliary cubed cell length, or
    a normalized structural ratio. Public volumetric parameters use ``V0``.
    An axial fitting adapter exposes the physical reference length ``L0`` and
    maps it internally to :math:`q_0=L_0^3`; the physical linear coefficient is
    :math:`\alpha_x=\alpha_q/3`.
    """

    def model(
        self,
        model: str | TemperatureEOSFamily | TemperatureEOSModel,
        variant: str | TemperatureEOSVariant | None = None,
    ) -> TemperatureEOSModel:
        """Return the canonical model specification."""
        return parse_temperature_eos_model(model, variant)

    def value(
        self,
        model: str | TemperatureEOSFamily | TemperatureEOSModel,
        parameters: ParameterLike,
        temperature: ArrayLike,
        *,
        variant: str | TemperatureEOSVariant | None = None,
    ) -> np.ndarray:
        """Evaluate the structural quantity at one or more temperatures.

        Raises
        ------
        ValueError
            If temperatures or parameters lie outside the mathematical domain.
        """
        spec = self.model(model, variant)
        pars = self._resolve_parameters(spec, parameters)
        temp = self._validate_temperature(temperature)
        value, _ = self._evaluate(spec, pars, temp)
        self._validate_output(value, "structural quantity")
        return value

    def expansion_coefficient(
        self,
        model: str | TemperatureEOSFamily | TemperatureEOSModel,
        parameters: ParameterLike,
        temperature: ArrayLike,
        *,
        variant: str | TemperatureEOSVariant | None = None,
    ) -> np.ndarray:
        """Evaluate the exact volumetric or auxiliary expansion coefficient."""
        spec = self.model(model, variant)
        pars = self._resolve_parameters(spec, parameters)
        temp = self._validate_temperature(temperature)
        value, alpha = self._evaluate(spec, pars, temp)
        self._validate_output(value, "structural quantity")
        self._validate_output(alpha, "thermal expansion coefficient", positive=False)
        return alpha

    def derivative(
        self,
        model: str | TemperatureEOSFamily | TemperatureEOSModel,
        parameters: ParameterLike,
        temperature: ArrayLike,
        *,
        variant: str | TemperatureEOSVariant | None = None,
    ) -> np.ndarray:
        r"""Evaluate :math:`\mathrm dX/\mathrm dT = \alpha X`."""
        spec = self.model(model, variant)
        pars = self._resolve_parameters(spec, parameters)
        temp = self._validate_temperature(temperature)
        value, alpha = self._evaluate(spec, pars, temp)
        result = value * alpha
        self._validate_output(result, "temperature derivative", positive=False)
        return result

    @staticmethod
    def linear_expansion_coefficient(auxiliary_alpha: ArrayLike) -> np.ndarray:
        """Convert the expansion of :math:`q=x^3` to linear expansion."""
        values = np.asarray(auxiliary_alpha, dtype=np.float64)
        if not np.all(np.isfinite(values)):
            raise ValueError("auxiliary expansion coefficients must be finite")
        return values / 3.0

    def _evaluate(
        self,
        model: TemperatureEOSModel,
        pars: TemperatureEOSParameters,
        temperature: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        if model.family is TemperatureEOSFamily.BERMAN:
            return self._berman(model, pars, temperature)
        if model.family is TemperatureEOSFamily.FEI:
            return self._fei(model, pars, temperature)
        if model.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL:
            return self._modified_holland_powell(model, pars, temperature)
        if model.family is TemperatureEOSFamily.SALJE:
            return self._salje(pars, temperature)
        return self._kroll_holland_powell(pars, temperature)

    @staticmethod
    def _berman(
        model: TemperatureEOSModel,
        p: TemperatureEOSParameters,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        alpha1 = 0.0 if model.variant is TemperatureEOSVariant.LINEAR else p.alpha1
        dt = t - p.temperature_ref
        ratio = 1.0 + p.alpha0 * dt + 0.5 * alpha1 * dt**2
        if np.any(ratio <= 0.0):
            raise ValueError("Berman model predicts a non-positive structural quantity")
        value = p.V0 * ratio
        alpha = (p.alpha0 + alpha1 * dt) / ratio
        return value, alpha

    @staticmethod
    def _fei(
        model: TemperatureEOSModel,
        p: TemperatureEOSParameters,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        alpha2 = 0.0 if model.variant is TemperatureEOSVariant.LINEAR else p.alpha2
        if np.any(t <= 0.0) or p.temperature_ref <= 0.0:
            raise ValueError("Fei model requires T > 0 K and T_ref > 0 K")
        exponent = (
            p.alpha0 * (t - p.temperature_ref)
            + 0.5 * p.alpha1 * (t**2 - p.temperature_ref**2)
            - alpha2 * (1.0 / t - 1.0 / p.temperature_ref)
        )
        value = p.V0 * np.exp(exponent)
        alpha = p.alpha0 + p.alpha1 * t + alpha2 / t**2
        return value, alpha

    @staticmethod
    def _modified_holland_powell(
        model: TemperatureEOSModel,
        p: TemperatureEOSParameters,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        if np.any(t <= 0.0) or p.temperature_ref <= 0.0:
            raise ValueError("modified Holland-Powell model requires T > 0 K")
        alpha1 = 0.0 if model.variant is TemperatureEOSVariant.SIMPLIFIED else p.alpha1
        coefficient = 10.0 * p.alpha0 + alpha1
        ratio = (
            1.0
            + p.alpha0 * (t - p.temperature_ref)
            - 2.0 * coefficient * (np.sqrt(t) - math.sqrt(p.temperature_ref))
        )
        if np.any(ratio <= 0.0):
            raise ValueError(
                "modified Holland-Powell model predicts a non-positive structural quantity"
            )
        value = p.V0 * ratio
        alpha = (p.alpha0 - coefficient / np.sqrt(t)) / ratio
        return value, alpha

    @staticmethod
    def _salje(
        p: TemperatureEOSParameters,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        if p.temperature_ref != 0.0:
            raise ValueError("Salje model is intrinsically referenced to T_ref = 0 K")
        if p.theta_sat <= 0.0:
            raise ValueError("theta_sat must be strictly positive for Salje")
        root_ref = p.V0 ** (1.0 / 3.0)
        value = np.empty_like(t)
        alpha = np.empty_like(t)
        zero = t == 0.0
        positive = ~zero
        value[zero] = p.V0
        alpha[zero] = 0.0
        if np.any(positive):
            tp = t[positive]
            z = p.theta_sat / tp
            coth = 1.0 / np.tanh(z)
            root = root_ref + p.p1 * p.theta_sat * (coth - 1.0)
            if np.any(root <= 0.0):
                raise ValueError("Salje model predicts a non-positive cube root")
            value[positive] = root**3
            # Stable csch²(z): use expm1 form for large z to avoid sinh overflow.
            csch_sq = np.empty_like(z)
            moderate = z < 350.0
            csch_sq[moderate] = 1.0 / np.sinh(z[moderate]) ** 2
            ez = np.exp(-2.0 * z[~moderate])
            csch_sq[~moderate] = 4.0 * ez / (1.0 - ez) ** 2
            alpha[positive] = 3.0 * p.p1 * p.theta_sat**2 * csch_sq / (root * tp**2)
        return value, alpha

    @staticmethod
    def _kroll_holland_powell(
        p: TemperatureEOSParameters,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        if p.temperature_ref <= 0.0:
            raise ValueError("KHP requires T_ref > 0 K")
        if p.theta_e <= 0.0:
            raise ValueError("theta_e must be strictly positive for KHP")
        if p.kp in {0.0, -1.0, -2.0}:
            raise ValueError("KHP is singular for kp = 0, -1, or -2")
        x_ref = p.theta_e / p.temperature_ref
        xi_ref = TemperatureEOS._einstein_xi(x_ref)
        occupation = TemperatureEOS._einstein_occupation(p.theta_e, t)
        occupation_ref = float(
            TemperatureEOS._einstein_occupation(
                p.theta_e, np.asarray([p.temperature_ref], dtype=np.float64)
            )[0]
        )
        A = p.alpha_ref * p.theta_e / xi_ref * (occupation - occupation_ref)
        C = p.kp * (p.kp + 2.0) / (p.kp + 1.0)
        exponent = -1.0 / (p.kp * (p.kp + 2.0))
        u = 1.0 - C * A
        if np.any(u <= 0.0):
            raise ValueError("KHP model lies outside its real-valued domain")
        ratio = -p.kp + (1.0 + p.kp) * u**exponent
        if np.any(ratio <= 0.0):
            raise ValueError("KHP model predicts a non-positive structural quantity")
        value = p.V0 * ratio
        Aprime = (
            p.alpha_ref
            * p.theta_e
            / xi_ref
            * TemperatureEOS._occupation_derivative(p.theta_e, t)
        )
        alpha = u ** (exponent - 1.0) * Aprime / ratio
        return value, alpha

    @staticmethod
    def _einstein_occupation(theta: float, t: np.ndarray) -> np.ndarray:
        result = np.zeros_like(t)
        positive = t > 0.0
        if np.any(positive):
            x = theta / t[positive]
            values = np.zeros_like(x)
            safe = x < 700.0
            values[safe] = 1.0 / np.expm1(x[safe])
            result[positive] = values
        return result

    @staticmethod
    def _occupation_derivative(theta: float, t: np.ndarray) -> np.ndarray:
        result = np.zeros_like(t)
        positive = t > 0.0
        if np.any(positive):
            tp = t[positive]
            x = theta / tp
            values = np.zeros_like(x)
            safe = x < 350.0
            exp_x = np.exp(x[safe])
            values[safe] = theta / tp[safe] ** 2 * exp_x / np.expm1(x[safe]) ** 2
            # exp(x)/expm1(x)^2 ~ exp(-x) for large x.
            values[~safe] = theta / tp[~safe] ** 2 * np.exp(-x[~safe])
            result[positive] = values
        return result

    @staticmethod
    def _einstein_xi(x: float) -> float:
        if x <= 0.0 or not math.isfinite(x):
            raise ValueError("Einstein argument must be finite and positive")
        if x < 350.0:
            return x * x * math.exp(x) / math.expm1(x) ** 2
        return x * x * math.exp(-x)

    @staticmethod
    def _resolve_parameters(
        model: TemperatureEOSModel,
        parameters: ParameterLike,
    ) -> TemperatureEOSParameters:
        if isinstance(parameters, TemperatureEOSParameters):
            pars = parameters
        else:
            known = set(TemperatureEOSParameters.__dataclass_fields__)
            unknown = set(parameters) - known
            if unknown:
                raise ValueError(
                    "unknown temperature-EOS parameter(s): "
                    + ", ".join(sorted(unknown))
                )
            missing = [name for name in model.parameter_names if name not in parameters]
            if missing:
                raise ValueError(
                    "missing temperature-EOS parameters: " + ", ".join(missing)
                )
            try:
                pars = TemperatureEOSParameters(
                    **{key: float(value) for key, value in parameters.items()}
                )
            except TypeError as exc:
                raise ValueError("invalid temperature-EOS parameter mapping") from exc
        return pars

    @staticmethod
    def _validate_temperature(temperature: ArrayLike) -> np.ndarray:
        values = np.asarray(temperature, dtype=np.float64)
        if values.ndim == 0:
            values = values.reshape(1)
        if values.ndim != 1:
            raise ValueError("temperature must be a scalar or one-dimensional array")
        if values.size == 0:
            raise ValueError("temperature cannot be empty")
        if not np.all(np.isfinite(values)):
            raise ValueError("temperature values must be finite")
        if np.any(values < 0.0):
            raise ValueError("temperature cannot be below absolute zero")
        return values

    @staticmethod
    def _validate_output(
        values: np.ndarray,
        name: str,
        *,
        positive: bool = True,
    ) -> None:
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} contains non-finite values")
        if positive and np.any(values <= 0.0):
            raise ValueError(f"{name} must remain strictly positive")


__all__ = [
    "TemperatureEOS",
    "TemperatureEOSFamily",
    "TemperatureEOSModel",
    "TemperatureEOSParameters",
    "TemperatureEOSVariant",
    "available_temperature_eos_models",
    "parse_temperature_eos_model",
]
