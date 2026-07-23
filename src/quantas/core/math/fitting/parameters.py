# -*- coding: utf-8 -*-

"""Parameter definitions and reduced/full mappings for numerical fitting."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, TypeVar

import numpy as np
from numpy.typing import ArrayLike, NDArray


ParameterDefinitionT = TypeVar("ParameterDefinitionT", bound="ParameterDefinition")


class ParameterState(str, Enum):
    """Role of one parameter in a fitting problem.

    Attributes
    ----------
    FREE
        Optimized directly by the numerical solver.
    FIXED
        Held at a user-specified value.
    IMPLIED
        Required by the model but calculated from other parameters and model
        conventions.
    DERIVED
        Reported after fitting but not passed to the model equation.
    """

    FREE = "free"
    FIXED = "fixed"
    IMPLIED = "implied"
    DERIVED = "derived"


@dataclass(frozen=True, slots=True)
class ParameterDefinition:
    """Describe one physical or mathematical fit parameter.

    Parameters
    ----------
    name : str
        Stable parameter identifier.
    state : ParameterState
        Role of the parameter in optimization and reporting.
    initial_value : float or None, optional
        Initial estimate. It is required for free parameters.
    value : float or None, optional
        Constant value for fixed parameters, or an optional constant fallback
        for implied or derived parameters.
    lower_bound, upper_bound : float, optional
        Admissible physical or numerical bounds.
    unit : str or None, optional
        Unit label used for reporting.
    description : str, optional
        Technical description of the parameter.
    metadata : mapping, optional
        Additional passive metadata.

    Raises
    ------
    ValueError
        If the definition is incomplete, non-finite, or has invalid bounds.
    """

    name: str
    state: ParameterState = ParameterState.FREE
    initial_value: float | None = None
    value: float | None = None
    lower_bound: float = -np.inf
    upper_bound: float = np.inf
    unit: str | None = None
    description: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize and validate the parameter definition."""
        normalized_name = str(self.name).strip()
        if not normalized_name:
            raise ValueError("parameter name cannot be empty")
        object.__setattr__(self, "name", normalized_name)
        object.__setattr__(self, "state", ParameterState(self.state))
        lower = float(self.lower_bound)
        upper = float(self.upper_bound)
        if np.isnan(lower) or np.isnan(upper) or lower >= upper:
            raise ValueError(
                f"parameter '{normalized_name}' requires lower_bound < upper_bound"
            )
        object.__setattr__(self, "lower_bound", lower)
        object.__setattr__(self, "upper_bound", upper)

        initial = _optional_finite(self.initial_value, "initial_value", normalized_name)
        value = _optional_finite(self.value, "value", normalized_name)
        object.__setattr__(self, "initial_value", initial)
        object.__setattr__(self, "value", value)

        if self.state is ParameterState.FREE:
            if initial is None:
                raise ValueError(
                    f"free parameter '{normalized_name}' requires an initial value"
                )
            if not lower <= initial <= upper:
                raise ValueError(
                    f"initial value for parameter '{normalized_name}' is outside bounds"
                )
            if value is not None:
                raise ValueError(
                    f"free parameter '{normalized_name}' cannot define a fixed value"
                )
        elif self.state is ParameterState.FIXED:
            if value is None:
                raise ValueError(
                    f"fixed parameter '{normalized_name}' requires a value"
                )
            if not lower <= value <= upper:
                raise ValueError(
                    f"fixed value for parameter '{normalized_name}' is outside bounds"
                )
        elif initial is not None:
            raise ValueError(
                f"{self.state.value} parameter '{normalized_name}' cannot define "
                "an initial value"
            )

        object.__setattr__(self, "unit", None if self.unit is None else str(self.unit))
        object.__setattr__(self, "description", str(self.description))
        object.__setattr__(self, "metadata", dict(self.metadata))

    @classmethod
    def free(
        cls: type[ParameterDefinitionT],
        name: str,
        initial_value: float,
        *,
        lower_bound: float = -np.inf,
        upper_bound: float = np.inf,
        unit: str | None = None,
        description: str = "",
        metadata: Mapping[str, Any] | None = None,
    ) -> ParameterDefinitionT:
        """Construct a free parameter definition."""
        return cls(
            name=name,
            state=ParameterState.FREE,
            initial_value=initial_value,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            unit=unit,
            description=description,
            metadata=dict(metadata or {}),
        )

    @classmethod
    def fixed(
        cls: type[ParameterDefinitionT],
        name: str,
        value: float,
        *,
        lower_bound: float = -np.inf,
        upper_bound: float = np.inf,
        unit: str | None = None,
        description: str = "",
        metadata: Mapping[str, Any] | None = None,
    ) -> ParameterDefinitionT:
        """Construct a fixed parameter definition."""
        return cls(
            name=name,
            state=ParameterState.FIXED,
            value=value,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            unit=unit,
            description=description,
            metadata=dict(metadata or {}),
        )

    @classmethod
    def implied(
        cls: type[ParameterDefinitionT],
        name: str,
        *,
        value: float | None = None,
        lower_bound: float = -np.inf,
        upper_bound: float = np.inf,
        unit: str | None = None,
        description: str = "",
        metadata: Mapping[str, Any] | None = None,
    ) -> ParameterDefinitionT:
        """Construct an implied parameter definition."""
        return cls(
            name=name,
            state=ParameterState.IMPLIED,
            value=value,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            unit=unit,
            description=description,
            metadata=dict(metadata or {}),
        )

    @classmethod
    def derived(
        cls: type[ParameterDefinitionT],
        name: str,
        *,
        value: float | None = None,
        lower_bound: float = -np.inf,
        upper_bound: float = np.inf,
        unit: str | None = None,
        description: str = "",
        metadata: Mapping[str, Any] | None = None,
    ) -> ParameterDefinitionT:
        """Construct a derived parameter definition."""
        return cls(
            name=name,
            state=ParameterState.DERIVED,
            value=value,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            unit=unit,
            description=description,
            metadata=dict(metadata or {}),
        )

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable parameter definition."""
        return {
            "name": self.name,
            "state": self.state.value,
            "initial_value": self.initial_value,
            "value": self.value,
            "lower_bound": self.lower_bound,
            "upper_bound": self.upper_bound,
            "unit": self.unit,
            "description": self.description,
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True, slots=True)
class ParameterSet:
    """Resolved complete parameter values in stable reporting order.

    Parameters
    ----------
    names : tuple of str
        Parameter names.
    values : array-like
        Resolved parameter values.
    states : tuple of ParameterState
        Parameter roles corresponding to ``names``.
    units : tuple of str or None
        Unit labels corresponding to ``names``.
    """

    names: tuple[str, ...]
    values: NDArray[np.float64]
    states: tuple[ParameterState, ...]
    units: tuple[str | None, ...]

    def __post_init__(self) -> None:
        """Normalize the resolved vector and validate aligned metadata."""
        names = tuple(str(name) for name in self.names)
        values = np.asarray(self.values, dtype=np.float64)
        states = tuple(ParameterState(state) for state in self.states)
        units = tuple(None if unit is None else str(unit) for unit in self.units)
        if values.ndim != 1 or values.size != len(names):
            raise ValueError("resolved parameter values must match parameter names")
        if len(states) != len(names) or len(units) != len(names):
            raise ValueError("resolved parameter metadata must have matching length")
        if not np.all(np.isfinite(values)):
            raise ValueError("resolved parameter values must be finite")
        object.__setattr__(self, "names", names)
        object.__setattr__(self, "values", values.copy())
        object.__setattr__(self, "states", states)
        object.__setattr__(self, "units", units)

    def as_mapping(self) -> dict[str, float]:
        """Return resolved values keyed by parameter name."""
        return {
            name: float(value)
            for name, value in zip(self.names, self.values, strict=True)
        }

    def model_values(self) -> NDArray[np.float64]:
        """Return parameters required by the model, excluding derived values."""
        return np.asarray(
            [
                value
                for value, state in zip(self.values, self.states, strict=True)
                if state is not ParameterState.DERIVED
            ],
            dtype=np.float64,
        )

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable resolved parameter set."""
        return {
            "names": list(self.names),
            "values": self.values.tolist(),
            "states": [state.value for state in self.states],
            "units": list(self.units),
        }


ParameterResolver = Callable[[Mapping[str, float]], Mapping[str, float]]


class ParameterMap:
    """Map reduced optimizer parameters to complete physical parameters.

    The mapping is general numerical infrastructure. A domain-specific model
    supplies passive definitions and, where needed, a resolver for implied and
    derived values. Solvers only see the reduced free vector and never optimize
    fixed or implied quantities.

    Parameters
    ----------
    definitions : sequence of ParameterDefinition
        Complete reporting-order parameter definitions.
    resolver : callable or None, optional
        Function that receives the current free and fixed values and returns
        values for declared implied or derived parameters.

    Raises
    ------
    ValueError
        If definitions are empty, duplicate names are present, or a mapping
        cannot be resolved consistently.
    """

    def __init__(
        self,
        definitions: Sequence[ParameterDefinition],
        *,
        resolver: ParameterResolver | None = None,
    ) -> None:
        normalized = tuple(definitions)
        if not normalized:
            raise ValueError("parameter map requires at least one definition")
        names = tuple(definition.name for definition in normalized)
        if len(set(names)) != len(names):
            raise ValueError("parameter names must be unique")
        if not any(
            definition.state is ParameterState.FREE for definition in normalized
        ):
            raise ValueError("parameter map requires at least one free parameter")
        self._definitions = normalized
        self._resolver = resolver
        self._index = {name: index for index, name in enumerate(names)}
        self._free_indices = tuple(
            index
            for index, definition in enumerate(normalized)
            if definition.state is ParameterState.FREE
        )

    @property
    def definitions(self) -> tuple[ParameterDefinition, ...]:
        """Return complete parameter definitions."""
        return self._definitions

    @property
    def names(self) -> tuple[str, ...]:
        """Return complete parameter names in reporting order."""
        return tuple(definition.name for definition in self._definitions)

    @property
    def states(self) -> tuple[ParameterState, ...]:
        """Return complete parameter states in reporting order."""
        return tuple(definition.state for definition in self._definitions)

    @property
    def free_names(self) -> tuple[str, ...]:
        """Return names optimized by a numerical solver."""
        return tuple(self._definitions[index].name for index in self._free_indices)

    @property
    def n_parameters(self) -> int:
        """Return the number of complete reported parameters."""
        return len(self._definitions)

    @property
    def n_free(self) -> int:
        """Return the number of optimized parameters."""
        return len(self._free_indices)

    def initial_free_values(self) -> NDArray[np.float64]:
        """Return the optimizer initial vector."""
        return np.asarray(
            [self._definitions[index].initial_value for index in self._free_indices],
            dtype=np.float64,
        )

    def free_bounds(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return lower and upper bounds in optimizer order."""
        lower = np.asarray(
            [self._definitions[index].lower_bound for index in self._free_indices],
            dtype=np.float64,
        )
        upper = np.asarray(
            [self._definitions[index].upper_bound for index in self._free_indices],
            dtype=np.float64,
        )
        return lower, upper

    def expand(self, free_parameters: ArrayLike) -> ParameterSet:
        """Resolve a free optimizer vector to complete physical parameters.

        Parameters
        ----------
        free_parameters : array-like
            Values in ``free_names`` order.

        Returns
        -------
        ParameterSet
            Complete values in reporting order.

        Raises
        ------
        ValueError
            If the free vector is invalid or implied values cannot be resolved.
        """
        free = np.asarray(free_parameters, dtype=np.float64)
        if free.ndim != 1 or free.size != self.n_free:
            raise ValueError(f"free parameter vector must contain {self.n_free} values")
        if not np.all(np.isfinite(free)):
            raise ValueError("free parameter values must be finite")

        values: dict[str, float] = {}
        for definition, value in zip(
            (self._definitions[index] for index in self._free_indices),
            free,
            strict=True,
        ):
            numeric = float(value)
            _validate_within_bounds(definition, numeric)
            values[definition.name] = numeric
        for definition in self._definitions:
            if definition.state is ParameterState.FIXED:
                assert definition.value is not None
                values[definition.name] = definition.value

        if self._resolver is not None:
            resolved = dict(self._resolver(dict(values)))
            for name, value in resolved.items():
                if name not in self._index:
                    raise ValueError(f"resolver returned unknown parameter '{name}'")
                definition = self._definitions[self._index[name]]
                if definition.state not in {
                    ParameterState.IMPLIED,
                    ParameterState.DERIVED,
                }:
                    raise ValueError(
                        f"resolver cannot overwrite {definition.state.value} "
                        f"parameter '{name}'"
                    )
                numeric = float(value)
                if not np.isfinite(numeric):
                    raise ValueError(f"resolved parameter '{name}' must be finite")
                values[name] = numeric

        for definition in self._definitions:
            if definition.name not in values and definition.value is not None:
                values[definition.name] = definition.value
            if definition.name not in values:
                raise ValueError(
                    f"no value was resolved for {definition.state.value} parameter "
                    f"'{definition.name}'"
                )
            _validate_within_bounds(definition, values[definition.name])

        return ParameterSet(
            names=self.names,
            values=np.asarray([values[name] for name in self.names], dtype=np.float64),
            states=self.states,
            units=tuple(definition.unit for definition in self._definitions),
        )

    def reduce(
        self, parameters: Mapping[str, float] | ArrayLike
    ) -> NDArray[np.float64]:
        """Extract optimizer values from a complete parameter representation."""
        if isinstance(parameters, Mapping):
            try:
                values = [float(parameters[name]) for name in self.free_names]
            except KeyError as exc:
                raise ValueError(
                    f"complete parameter mapping is missing '{exc.args[0]}'"
                ) from exc
            free = np.asarray(values, dtype=np.float64)
        else:
            complete = np.asarray(parameters, dtype=np.float64)
            if complete.ndim != 1 or complete.size != self.n_parameters:
                raise ValueError(
                    f"complete parameter vector must contain {self.n_parameters} values"
                )
            free = complete[np.asarray(self._free_indices, dtype=int)]
        if not np.all(np.isfinite(free)):
            raise ValueError("reduced parameter values must be finite")
        return free.copy()

    def resolved_jacobian(
        self,
        free_parameters: ArrayLike,
        *,
        relative_step: float | None = None,
    ) -> NDArray[np.float64]:
        """Return the Jacobian of complete values with respect to free values.

        Free rows are set analytically to the identity and fixed rows to zero.
        Implied and derived rows are evaluated with bound-aware numerical
        differences through the resolver.
        """
        free = np.asarray(free_parameters, dtype=np.float64)
        base = self.expand(free).values
        lower, upper = self.free_bounds()
        step_scale = (
            float(relative_step)
            if relative_step is not None
            else float(np.cbrt(np.finfo(np.float64).eps))
        )
        if not np.isfinite(step_scale) or step_scale <= 0.0:
            raise ValueError("relative_step must be a positive finite value")
        jacobian = np.zeros((self.n_parameters, self.n_free), dtype=np.float64)

        for free_position, parameter_index in enumerate(self._free_indices):
            jacobian[parameter_index, free_position] = 1.0
            step = step_scale * max(1.0, abs(float(free[free_position])))
            can_lower = free[free_position] - step >= lower[free_position]
            can_upper = free[free_position] + step <= upper[free_position]
            if can_lower and can_upper:
                plus = free.copy()
                minus = free.copy()
                plus[free_position] += step
                minus[free_position] -= step
                derivative = (self.expand(plus).values - self.expand(minus).values) / (
                    2.0 * step
                )
            elif can_upper:
                plus = free.copy()
                plus[free_position] += step
                derivative = (self.expand(plus).values - base) / step
            elif can_lower:
                minus = free.copy()
                minus[free_position] -= step
                derivative = (base - self.expand(minus).values) / step
            else:
                raise ValueError(
                    f"cannot differentiate parameter '{self.free_names[free_position]}' "
                    "inside its bounds"
                )
            for row, definition in enumerate(self._definitions):
                if definition.state in {
                    ParameterState.IMPLIED,
                    ParameterState.DERIVED,
                }:
                    jacobian[row, free_position] = derivative[row]
        return jacobian

    def propagate_covariance(
        self,
        free_covariance: ArrayLike,
        free_parameters: ArrayLike,
    ) -> NDArray[np.float64]:
        """Propagate optimizer covariance to complete physical parameters."""
        covariance = np.asarray(free_covariance, dtype=np.float64)
        expected = (self.n_free, self.n_free)
        if covariance.shape != expected:
            raise ValueError(f"free covariance must have shape {expected}")
        if not np.all(np.isfinite(covariance)):
            raise ValueError("free covariance must contain finite values")
        jacobian = self.resolved_jacobian(free_parameters)
        propagated = jacobian @ covariance @ jacobian.T
        return np.asarray(0.5 * (propagated + propagated.T), dtype=np.float64)

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable mapping definition."""
        return {
            "definitions": [definition.as_dict() for definition in self._definitions],
            "free_names": list(self.free_names),
        }


def _optional_finite(value: float | None, label: str, name: str) -> float | None:
    """Normalize an optional finite scalar."""
    if value is None:
        return None
    numeric = float(value)
    if not np.isfinite(numeric):
        raise ValueError(f"parameter '{name}' {label} must be finite")
    return numeric


def _validate_within_bounds(definition: ParameterDefinition, value: float) -> None:
    """Validate one resolved value against its declared bounds."""
    if not definition.lower_bound <= value <= definition.upper_bound:
        raise ValueError(
            f"parameter '{definition.name}' value {value} is outside "
            f"[{definition.lower_bound}, {definition.upper_bound}]"
        )
