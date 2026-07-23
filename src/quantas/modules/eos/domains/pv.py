# -*- coding: utf-8 -*-

"""Pressure-volume fit models and parameter preparation for EOS."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
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
    EOSModel,
    PressureEOS,
    implied_kp,
    implied_kpp,
    parse_eos_model,
)

from ..models import ParameterConstraint

_PRESSURE_PARAMETER_ORDER = ("K0", "KP", "KPP", "V0")
_AXIAL_PARAMETER_ORDER = ("M0", "MP", "MPP", "L0")
_POSITIVE_LOWER_BOUND = float(np.nextafter(0.0, 1.0))


class PressureEOSFitModel(BaseFitModel):
    """Adapt the shared pressure EOS core to the generic fitting contract.

    Parameters
    ----------
    model : EOSModel or str
        EOS family and order evaluated by the shared physical core.
    initial_parameters : mapping, optional
        Complete initial physical values used only by the unmapped model API.
        The normal EOS workflow supplies a :class:`ParameterMap` and
        therefore obtains its reduced initial vector from that mapping.
    """

    def __init__(
        self,
        model: EOSModel | str,
        *,
        initial_parameters: Mapping[str, float] | None = None,
    ) -> None:
        self.eos_model = parse_eos_model(model)
        self._pressure = PressureEOS()
        defaults = {"K0": 1.0, "KP": 4.0, "KPP": 0.0, "V0": 1.0}
        if initial_parameters is not None:
            defaults.update(
                {
                    str(name).upper(): float(value)
                    for name, value in initial_parameters.items()
                }
            )
        self._initial = np.asarray(
            [defaults[name] for name in _PRESSURE_PARAMETER_ORDER],
            dtype=np.float64,
        )

    @property
    def name(self) -> str:
        """Return a stable physical-model identifier."""
        return f"pressure_eos:{self.eos_model.tag}"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return complete physical pressure-EOS parameter order."""
        return _PRESSURE_PARAMETER_ORDER

    def evaluate(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Evaluate pressure at the supplied volumes.

        Parameters
        ----------
        x : array-like
            Positive volumes.
        parameters : array-like
            Complete ``K0, KP, KPP, V0`` physical parameter vector.

        Returns
        -------
        ndarray
            Calculated pressures as a ``float64`` array.
        """
        values = np.asarray(parameters, dtype=np.float64)
        if values.ndim != 1 or values.size != len(_PRESSURE_PARAMETER_ORDER):
            raise ValueError("pressure EOS model requires K0, KP, KPP, V0")
        physical = dict(zip(_PRESSURE_PARAMETER_ORDER, values, strict=True))
        return self._pressure.pressure(self.eos_model, physical, x)

    def derivative_x(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        r"""Evaluate the analytical pressure derivative with respect to volume.

        For every isothermal EOS implemented by the shared physical core,

        .. math::

            \frac{\partial P}{\partial V}=-\frac{K_T(V)}{V}.

        Parameters
        ----------
        x : array-like
            Positive volumes.
        parameters : array-like
            Complete ``K0, KP, KPP, V0`` physical parameter vector.

        Returns
        -------
        ndarray
            Analytical :math:`\partial P/\partial V` values.
        """
        volume = np.asarray(x, dtype=np.float64)
        values = np.asarray(parameters, dtype=np.float64)
        if values.ndim != 1 or values.size != len(_PRESSURE_PARAMETER_ORDER):
            raise ValueError("pressure EOS model requires K0, KP, KPP, V0")
        physical = dict(zip(_PRESSURE_PARAMETER_ORDER, values, strict=True))
        bulk = self._pressure.bulk_modulus(self.eos_model, physical, volume)
        return -np.asarray(bulk, dtype=np.float64) / volume

    def initial_guess(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return configured complete initial values after data validation."""
        volume, pressure = _validate_pressure_data(x, y)
        del volume, pressure
        return self._initial.copy()

    def bounds(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return minimally restrictive physical bounds for the full model."""
        volume, pressure = _validate_pressure_data(x, y)
        del volume, pressure
        lower: np.ndarray = np.asarray(
            [_POSITIVE_LOWER_BOUND, -np.inf, -np.inf, _POSITIVE_LOWER_BOUND],
            dtype=np.float64,
        )
        upper: np.ndarray = np.full(4, np.inf, dtype=np.float64)
        return lower, upper

    def metadata(self) -> dict[str, Any]:
        """Return EOS family, order, and physical parameter metadata."""
        return {
            **super().metadata(),
            "eos_model": self.eos_model.as_dict(),
            "relationship": "pressure(volume)",
        }


class AxialEOSFitModel(BaseFitModel):
    r"""Adapt a linear-axis EOS to the generic fitting contract.

    Quantas follows the finite-strain convention described by Angel et al.
    (2014): a measured length :math:`x` is transformed to the auxiliary
    volume-like coordinate :math:`q=x^3`, and the selected volumetric EOS is
    evaluated in :math:`q`.  Public parameters remain linear quantities:

    .. math::

        K_0 = M_0/3,\quad K'_0 = M'_0/3,\quad
        K''_0 = M''_0/3,\quad V_0 = X_0^3.

    Parameters
    ----------
    model : EOSModel or str
        EOS family and order evaluated by the shared physical core.
    target : {"a", "b", "c"}
        Linear cell parameter represented by the model.
    initial_parameters : mapping, optional
        Complete initial ``M0, MP, MPP, L0`` values used by the unmapped
        model API.

    Raises
    ------
    ValueError
        If ``target`` is not a supported linear cell parameter.
    """

    def __init__(
        self,
        model: EOSModel | str,
        target: str,
        *,
        initial_parameters: Mapping[str, float] | None = None,
    ) -> None:
        normalized_target = str(target).lower()
        if normalized_target not in {"a", "b", "c"}:
            raise ValueError("axial EOS target must be 'a', 'b', or 'c'")
        self.eos_model = parse_eos_model(model)
        self.target = normalized_target
        self._pressure = PressureEOS()
        defaults = {"M0": 3.0, "MP": 12.0, "MPP": 0.0, "L0": 1.0}
        if initial_parameters is not None:
            defaults.update(
                {
                    str(name).upper(): float(value)
                    for name, value in initial_parameters.items()
                }
            )
        self._initial = np.asarray(
            [defaults[name] for name in _AXIAL_PARAMETER_ORDER],
            dtype=np.float64,
        )

    @property
    def name(self) -> str:
        """Return a stable physical-model identifier."""
        return f"axial_pressure_eos:{self.target}:{self.eos_model.tag}"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return complete linear parameter reporting order."""
        return _AXIAL_PARAMETER_ORDER

    def evaluate(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Evaluate pressure at cubed-length coordinates.

        Parameters
        ----------
        x : array-like
            Positive auxiliary coordinates :math:`q=x_{axis}^3`.
        parameters : array-like
            Complete ``M0, MP, MPP, L0`` linear parameter vector.

        Returns
        -------
        ndarray
            Calculated pressures as a ``float64`` array.
        """
        coordinate = np.asarray(x, dtype=np.float64)
        if np.any(coordinate <= 0.0):
            raise ValueError("cubed axial coordinates must be strictly positive")
        physical = axial_to_volume_parameters(parameters)
        return self._pressure.pressure(self.eos_model, physical, coordinate)

    def derivative_x(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        r"""Return :math:`\partial P/\partial q` for :math:`q=x^3`."""
        coordinate = np.asarray(x, dtype=np.float64)
        physical = axial_to_volume_parameters(parameters)
        bulk = self._pressure.bulk_modulus(self.eos_model, physical, coordinate)
        return -np.asarray(bulk, dtype=np.float64) / coordinate

    def derivative_length(
        self,
        length: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        r"""Return the analytical derivative :math:`\partial P/\partial x`.

        The chain rule gives

        .. math::

            \frac{\partial P}{\partial x}
            = 3x^2 \frac{\partial P}{\partial q}
            = -\frac{M(x)}{x}.
        """
        axis = np.asarray(length, dtype=np.float64)
        if np.any(axis <= 0.0):
            raise ValueError("axial lengths must be strictly positive")
        return 3.0 * axis**2 * self.derivative_x(axis**3, parameters)

    def linear_modulus(
        self,
        length: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return the instantaneous linear modulus along the fitted axis."""
        axis = np.asarray(length, dtype=np.float64)
        physical = axial_to_volume_parameters(parameters)
        bulk = self._pressure.bulk_modulus(self.eos_model, physical, axis**3)
        return 3.0 * np.asarray(bulk, dtype=np.float64)

    def initial_guess(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return configured complete initial values after validation."""
        coordinate, pressure = _validate_pressure_data(x, y)
        del coordinate, pressure
        return self._initial.copy()

    def bounds(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return minimally restrictive bounds for linear parameters."""
        coordinate, pressure = _validate_pressure_data(x, y)
        del coordinate, pressure
        lower = np.asarray(
            [_POSITIVE_LOWER_BOUND, -np.inf, -np.inf, _POSITIVE_LOWER_BOUND],
            dtype=np.float64,
        )
        return lower, np.full(4, np.inf, dtype=np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return EOS, target, and coordinate-transformation metadata."""
        return {
            **super().metadata(),
            "eos_model": self.eos_model.as_dict(),
            "relationship": f"pressure({self.target})",
            "linear_target": self.target,
            "coordinate_transform": "q=x^3",
            "parameter_transform": {
                "K0": "M0/3",
                "KP": "MP/3",
                "KPP": "MPP/3",
                "V0": "L0^3",
            },
        }


def axial_to_volume_parameters(
    parameters: Mapping[str, float] | np.ndarray | Sequence[float],
) -> dict[str, float]:
    """Convert complete linear parameters to auxiliary volumetric values.

    Parameters
    ----------
    parameters : mapping or array-like
        ``M0, MP, MPP, L0`` linear quantities.

    Returns
    -------
    dict
        ``K0, KP, KPP, V0`` values consumed by the shared pressure EOS core.
    """
    if isinstance(parameters, Mapping):
        values = {name: float(parameters[name]) for name in _AXIAL_PARAMETER_ORDER}
    else:
        array = np.asarray(parameters, dtype=np.float64)
        if array.ndim != 1 or array.size != len(_AXIAL_PARAMETER_ORDER):
            raise ValueError("axial EOS model requires M0, MP, MPP, L0")
        values = dict(zip(_AXIAL_PARAMETER_ORDER, array, strict=True))
    if values["M0"] <= 0.0 or values["L0"] <= 0.0:
        raise ValueError("M0 and L0 must be strictly positive")
    return {
        "K0": values["M0"] / 3.0,
        "KP": values["MP"] / 3.0,
        "KPP": values["MPP"] / 3.0,
        "V0": values["L0"] ** 3,
    }


def estimate_pressure_parameters(
    volume: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
) -> dict[str, float]:
    """Estimate a complete pressure-EOS parameter set from P-V data.

    The estimate follows the historical Quantas/EosFit strategy while using
    scaled :class:`numpy.polynomial.Polynomial` fits for numerical stability.
    A low-order polynomial in ``V(P)`` provides the zero-pressure volume and a
    polynomial in ``P(V)`` provides the local slope used for ``K0``.

    Parameters
    ----------
    volume, pressure : array-like
        Positive volumes and corresponding pressures.

    Returns
    -------
    dict
        Initial values for ``K0``, ``KP``, ``KPP``, and ``V0``.

    Raises
    ------
    ValueError
        If fewer than three valid data points are provided or an estimate
        cannot be made finite and physical.
    """
    volume_array, pressure_array = _validate_pressure_data(volume, pressure)
    if volume_array.size < 3:
        raise ValueError("at least three P-V points are required for an EOS estimate")

    degree = min(3, int(volume_array.size) - 1)
    volume_of_pressure = Polynomial.fit(pressure_array, volume_array, degree)
    v0 = float(volume_of_pressure(0.0))
    if not np.isfinite(v0) or v0 <= 0.0:
        nonnegative = volume_array[pressure_array >= 0.0]
        source = nonnegative if nonnegative.size else volume_array
        v0 = float(np.max(source))

    pressure_of_volume = Polynomial.fit(volume_array, pressure_array, degree)
    slope = float(pressure_of_volume.deriv()(v0))
    k0 = float(-v0 * slope)
    if not np.isfinite(k0) or k0 <= 0.0:
        order = np.argsort(np.abs(pressure_array))
        count = min(max(3, degree + 1), volume_array.size)
        selected = order[:count]
        local = Polynomial.fit(volume_array[selected], pressure_array[selected], 1)
        k0 = float(-v0 * local.deriv()(v0))
    if not np.isfinite(k0) or k0 <= 0.0:
        raise ValueError("P-V data did not produce a positive initial bulk modulus")

    kp = 4.0
    return {"K0": k0, "KP": kp, "KPP": 0.0, "V0": v0}


def estimate_axial_parameters(
    length: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
) -> dict[str, float]:
    """Estimate complete linear-EOS parameters from pressure-length data.

    The measured length is cubed before calling the shared P-V estimator.
    Auxiliary volumetric estimates are then converted to linear quantities.

    Parameters
    ----------
    length, pressure : array-like
        Positive cell-axis lengths and corresponding pressures.

    Returns
    -------
    dict
        Initial values for ``M0``, ``MP``, ``MPP``, and ``L0``.
    """
    axis, pressure_array = _validate_axial_data(length, pressure)
    auxiliary = estimate_pressure_parameters(axis**3, pressure_array)
    return {
        "M0": 3.0 * auxiliary["K0"],
        "MP": 3.0 * auxiliary["KP"],
        "MPP": 3.0 * auxiliary["KPP"],
        "L0": float(np.cbrt(auxiliary["V0"])),
    }


def build_axial_parameter_map(
    model: EOSModel | str,
    length: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
    constraints: Sequence[ParameterConstraint] = (),
    *,
    pressure_unit: str = "GPa",
    length_unit: str = "angstrom",
) -> ParameterMap:
    """Build a reduced/full mapping for one linear-axis EOS fit.

    Parameters
    ----------
    model : EOSModel or str
        EOS family and order applied to the cubed length.
    length, pressure : array-like
        Selected pressure-axis observations used for initial estimates.
    constraints : sequence of ParameterConstraint, optional
        Overrides using the canonical linear names ``M0``, ``MP``, ``MPP``,
        and ``L0``.
    pressure_unit, length_unit : str, optional
        Units retained for reporting.

    Returns
    -------
    ParameterMap
        Mapping in stable ``M0, MP, MPP, L0`` reporting order.
    """
    eos_model = parse_eos_model(model)
    estimates = _complete_axial_initial_estimates(eos_model, length, pressure)
    overrides = _axial_constraint_overrides(constraints)
    definitions = tuple(
        _axial_parameter_definition(
            eos_model,
            name,
            estimates,
            overrides.get(name),
            pressure_unit=pressure_unit,
            length_unit=length_unit,
        )
        for name in _AXIAL_PARAMETER_ORDER
    )
    return ParameterMap(definitions, resolver=_axial_resolver(eos_model))


def build_pressure_parameter_map(
    model: EOSModel | str,
    volume: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
    constraints: Sequence[ParameterConstraint] = (),
    *,
    pressure_unit: str = "GPa",
    volume_unit: str = "angstrom^3",
) -> ParameterMap:
    """Build the reduced/full parameter mapping for one pressure-EOS fit.

    Parameters
    ----------
    model : EOSModel or str
        EOS family and order.
    volume, pressure : array-like
        Selected P-V observations used for initial estimates.
    constraints : sequence of ParameterConstraint, optional
        Explicit user overrides. Parameters implied by EOS order cannot be
        made free or fixed, while fitted parameters may be free or fixed.
    pressure_unit, volume_unit : str, optional
        Units retained for K0/KPP and V0 reporting. The equations themselves
        are unit-consistent and do not perform silent conversion.

    Returns
    -------
    ParameterMap
        Mapping with reporting order ``K0, KP, KPP, V0``.

    Raises
    ------
    ValueError
        If constraints are unknown, duplicated, or conflict with EOS order.
    """
    eos_model = parse_eos_model(model)
    estimates = _complete_initial_estimates(eos_model, volume, pressure)
    overrides = _constraint_overrides(constraints)
    definitions = tuple(
        _pressure_parameter_definition(
            eos_model,
            name,
            estimates,
            overrides.get(name),
            pressure_unit=pressure_unit,
            volume_unit=volume_unit,
        )
        for name in _PRESSURE_PARAMETER_ORDER
    )
    return ParameterMap(definitions, resolver=_pressure_resolver(eos_model))


def _complete_axial_initial_estimates(
    model: EOSModel,
    length: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
) -> dict[str, float]:
    """Return linear initial values consistent with EOS order."""
    estimates = estimate_axial_parameters(length, pressure)
    if model.order == 2:
        kp_value = implied_kp(model)
        if kp_value is None:
            raise ValueError(f"{model.tag} does not imply KP")
        estimates["MP"] = 3.0 * kp_value
    estimates["MPP"] = 3.0 * implied_kpp(
        model,
        estimates["M0"] / 3.0,
        estimates["MP"] / 3.0,
    )
    return estimates


def _axial_parameter_definition(
    model: EOSModel,
    name: str,
    estimates: Mapping[str, float],
    override: ParameterConstraint | None,
    *,
    pressure_unit: str,
    length_unit: str,
) -> ParameterDefinition:
    """Build one fitted, fixed, or implied linear-EOS parameter."""
    source_name = {"M0": "K0", "MP": "KP", "MPP": "KPP", "L0": "V0"}[name]
    source = model.parameter_sources[source_name]
    unit = _axial_parameter_unit(name, pressure_unit, length_unit)
    if source == "implied":
        return _axial_implied_definition(model, name, override, unit)
    return _axial_fitted_definition(name, estimates[name], override, unit)


def _axial_implied_definition(
    model: EOSModel,
    name: str,
    override: ParameterConstraint | None,
    unit: str,
) -> ParameterDefinition:
    """Build one linear parameter imposed by EOS family/order."""
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
        description=f"linear parameter implied by {model.tag}",
        metadata={"source": "eos_order", "linear": True},
    )


def _axial_fitted_definition(
    name: str,
    estimate: float,
    override: ParameterConstraint | None,
    unit: str,
) -> ParameterDefinition:
    """Build one normally fitted linear parameter."""
    default_lower = _POSITIVE_LOWER_BOUND if name in {"M0", "L0"} else -np.inf
    if override is None:
        return ParameterDefinition.free(
            name,
            estimate,
            lower_bound=default_lower,
            unit=unit,
            description=_axial_parameter_description(name),
            metadata={"initial_source": "cubed_length_polynomial_estimate"},
        )
    if override.state is ParameterState.FREE:
        initial = estimate if override.initial_value is None else override.initial_value
        return ParameterDefinition.free(
            name,
            initial,
            lower_bound=override.lower_bound,
            upper_bound=override.upper_bound,
            unit=override.unit or unit,
            description=override.description or _axial_parameter_description(name),
            metadata={**override.metadata, "initial_source": "user"},
        )
    if override.state is ParameterState.FIXED:
        if override.value is None:
            raise ValueError(f"fixed parameter {name} requires a value")
        return ParameterDefinition.fixed(
            name,
            override.value,
            lower_bound=override.lower_bound,
            upper_bound=override.upper_bound,
            unit=override.unit or unit,
            description=override.description or _axial_parameter_description(name),
            metadata={**override.metadata, "source": "user_fixed"},
        )
    raise ValueError(
        f"fitted parameter {name} cannot be declared {override.state.value}"
    )


def _axial_resolver(model: EOSModel):
    """Return a resolver for linear parameters implied by one EOS model."""
    sources = model.parameter_sources

    def resolver(values: Mapping[str, float]) -> Mapping[str, float]:
        m0 = float(values["M0"])
        if sources["KP"] == "implied":
            implied = implied_kp(model)
            if implied is None:
                raise ValueError(f"{model.tag} does not imply KP")
            mp = 3.0 * float(implied)
        else:
            mp = float(values["MP"])
        resolved: dict[str, float] = {}
        if sources["KP"] == "implied":
            resolved["MP"] = mp
        if sources["KPP"] == "implied":
            resolved["MPP"] = 3.0 * implied_kpp(model, m0 / 3.0, mp / 3.0)
        return resolved

    return resolver


def _complete_initial_estimates(
    model: EOSModel,
    volume: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
) -> dict[str, float]:
    """Return initial values consistent with the selected EOS order."""
    estimates = estimate_pressure_parameters(volume, pressure)
    if model.order == 2:
        kp_value = implied_kp(model)
        if kp_value is None:
            raise ValueError(f"{model.tag} does not imply KP")
        estimates["KP"] = kp_value
    estimates["KPP"] = implied_kpp(model, estimates["K0"], estimates["KP"])
    return estimates


def _pressure_parameter_definition(
    model: EOSModel,
    name: str,
    estimates: Mapping[str, float],
    override: ParameterConstraint | None,
    *,
    pressure_unit: str,
    volume_unit: str,
) -> ParameterDefinition:
    """Build one fitted, fixed, or implied pressure-EOS parameter."""
    source = model.parameter_sources[name]
    unit = _parameter_unit(name, pressure_unit, volume_unit)
    if source == "implied":
        return _implied_definition(model, name, override, unit)
    return _fitted_definition(name, estimates[name], override, unit)


def _implied_definition(
    model: EOSModel,
    name: str,
    override: ParameterConstraint | None,
    unit: str,
) -> ParameterDefinition:
    """Build one parameter imposed by EOS family/order conventions."""
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
        description=f"parameter implied by {model.tag}",
        metadata={"source": "eos_order"},
    )


def _fitted_definition(
    name: str,
    estimate: float,
    override: ParameterConstraint | None,
    unit: str,
) -> ParameterDefinition:
    """Build one normally fitted parameter with an optional user override."""
    default_lower = _POSITIVE_LOWER_BOUND if name in {"K0", "V0"} else -np.inf
    if override is None:
        return ParameterDefinition.free(
            name,
            estimate,
            lower_bound=default_lower,
            unit=unit,
            description=_parameter_description(name),
            metadata={"initial_source": "polynomial_estimate"},
        )
    if override.state is ParameterState.FREE:
        initial = estimate if override.initial_value is None else override.initial_value
        return ParameterDefinition.free(
            name,
            initial,
            lower_bound=override.lower_bound,
            upper_bound=override.upper_bound,
            unit=override.unit or unit,
            description=override.description or _parameter_description(name),
            metadata={**override.metadata, "initial_source": "user"},
        )
    if override.state is ParameterState.FIXED:
        if override.value is None:
            raise ValueError(f"fixed parameter {name} requires a value")
        return ParameterDefinition.fixed(
            name,
            override.value,
            lower_bound=override.lower_bound,
            upper_bound=override.upper_bound,
            unit=override.unit or unit,
            description=override.description or _parameter_description(name),
            metadata={**override.metadata, "source": "user_fixed"},
        )
    raise ValueError(
        f"fitted parameter {name} cannot be declared {override.state.value}"
    )


def _pressure_resolver(model: EOSModel):
    """Return a resolver for parameters implied by one EOS model."""
    sources = model.parameter_sources

    def resolver(values: Mapping[str, float]) -> Mapping[str, float]:
        k0 = float(values["K0"])
        if sources["KP"] == "implied":
            implied = implied_kp(model)
            if implied is None:
                raise ValueError(f"{model.tag} does not imply KP")
            kp = float(implied)
        else:
            kp = float(values["KP"])
        resolved: dict[str, float] = {}
        if sources["KP"] == "implied":
            resolved["KP"] = kp
        if sources["KPP"] == "implied":
            resolved["KPP"] = implied_kpp(model, k0, kp)
        return resolved

    return resolver


def _axial_constraint_overrides(
    constraints: Sequence[ParameterConstraint],
) -> dict[str, ParameterConstraint]:
    """Normalize linear-axis constraints by canonical parameter name."""
    overrides: dict[str, ParameterConstraint] = {}
    for constraint in constraints:
        name = constraint.name.upper()
        if name not in _AXIAL_PARAMETER_ORDER:
            raise ValueError(f"unknown axial-EOS parameter constraint: {name}")
        if name in overrides:
            raise ValueError(f"duplicate axial-EOS parameter constraint: {name}")
        overrides[name] = (
            constraint if constraint.name == name else replace(constraint, name=name)
        )
    return overrides


def _axial_parameter_unit(name: str, pressure_unit: str, length_unit: str) -> str:
    """Return units for linear EOS parameters."""
    return {
        "M0": str(pressure_unit),
        "MP": "1",
        "MPP": f"({pressure_unit})^-1",
        "L0": str(length_unit),
    }[name]


def _axial_parameter_description(name: str) -> str:
    """Return a technical linear-parameter description."""
    return {
        "M0": "reference linear modulus, -x(dP/dx)",
        "MP": "first pressure derivative of the linear modulus",
        "MPP": "second pressure derivative of the linear modulus",
        "L0": "reference zero-pressure cell-axis length",
    }[name]


def _constraint_overrides(
    constraints: Sequence[ParameterConstraint],
) -> dict[str, ParameterConstraint]:
    """Normalize user constraints by canonical parameter name."""
    overrides: dict[str, ParameterConstraint] = {}
    for constraint in constraints:
        name = constraint.name.upper()
        if name not in _PRESSURE_PARAMETER_ORDER:
            raise ValueError(f"unknown pressure-EOS parameter constraint: {name}")
        if name in overrides:
            raise ValueError(f"duplicate pressure-EOS parameter constraint: {name}")
        overrides[name] = (
            constraint if constraint.name == name else replace(constraint, name=name)
        )
    return overrides


def _parameter_unit(name: str, pressure_unit: str, volume_unit: str) -> str:
    """Return units consistent with the selected P-V observations."""
    units = {
        "K0": str(pressure_unit),
        "KP": "1",
        "KPP": f"({pressure_unit})^-1",
        "V0": str(volume_unit),
    }
    return units[name]


def _parameter_description(name: str) -> str:
    """Return a technical parameter description."""
    descriptions = {
        "K0": "reference isothermal bulk modulus",
        "KP": "first pressure derivative of the bulk modulus",
        "KPP": "second pressure derivative of the bulk modulus",
        "V0": "reference zero-pressure volume",
    }
    return descriptions[name]


def _validate_pressure_data(
    volume: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return finite one-dimensional P-V arrays with positive volumes."""
    volume_array = np.asarray(volume, dtype=np.float64)
    pressure_array = np.asarray(pressure, dtype=np.float64)
    if volume_array.ndim != 1 or pressure_array.ndim != 1:
        raise ValueError("P-V arrays must be one-dimensional")
    if volume_array.shape != pressure_array.shape or volume_array.size < 2:
        raise ValueError("P-V arrays must have equal length and at least two points")
    if not np.all(np.isfinite(volume_array)) or not np.all(np.isfinite(pressure_array)):
        raise ValueError("P-V arrays must contain finite values")
    if np.any(volume_array <= 0.0):
        raise ValueError("P-V volumes must be strictly positive")
    return volume_array.copy(), pressure_array.copy()


def _validate_axial_data(
    length: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return finite one-dimensional pressure-length arrays."""
    axis = np.asarray(length, dtype=np.float64)
    pressure_array = np.asarray(pressure, dtype=np.float64)
    if axis.ndim != 1 or pressure_array.ndim != 1:
        raise ValueError("pressure-axis arrays must be one-dimensional")
    if axis.shape != pressure_array.shape or axis.size < 2:
        raise ValueError(
            "pressure-axis arrays must have equal length and at least two points"
        )
    if not np.all(np.isfinite(axis)) or not np.all(np.isfinite(pressure_array)):
        raise ValueError("pressure-axis arrays must contain finite values")
    if np.any(axis <= 0.0):
        raise ValueError("cell-axis lengths must be strictly positive")
    return axis.copy(), pressure_array.copy()


__all__ = [
    "AxialEOSFitModel",
    "PressureEOSFitModel",
    "axial_to_volume_parameters",
    "build_axial_parameter_map",
    "build_pressure_parameter_map",
    "estimate_axial_parameters",
    "estimate_pressure_parameters",
]
