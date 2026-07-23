# -*- coding: utf-8 -*-

"""Model contracts used by the general fitting infrastructure."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable, Sequence
from typing import Any, TypeAlias

import numpy as np

from .diagnostics import validate_xy
from .parameters import ParameterMap, ParameterState
from .result import FitResult

ArrayLike: TypeAlias = np.ndarray | Sequence[float]
ModelFunction: TypeAlias = Callable[..., np.ndarray]


class BaseFitModel(ABC):
    """Define a mathematical or physical model that can be fitted to data.

    Subclasses provide the model equation, parameter names and initial
    estimates. Optimizer selection and diagnostic handling are delegated to a
    fitter object.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the model name.

        Returns
        -------
        str
            Stable model identifier.
        """

    @property
    @abstractmethod
    def parameter_names(self) -> tuple[str, ...]:
        """Return parameter names in fitting order.

        Returns
        -------
        tuple of str
            Ordered parameter names.
        """

    @abstractmethod
    def evaluate(self, x: ArrayLike, parameters: ArrayLike) -> np.ndarray:
        """Evaluate the model.

        Parameters
        ----------
        x : array-like
            Independent coordinates.
        parameters : array-like
            Model parameters in ``parameter_names`` order.

        Returns
        -------
        ndarray
            Model values.
        """

    @abstractmethod
    def initial_guess(self, x: ArrayLike, y: ArrayLike) -> np.ndarray:
        """Estimate initial parameters from observed data.

        Parameters
        ----------
        x, y : array-like
            Coordinates and observations.

        Returns
        -------
        ndarray
            Initial parameter estimate.
        """

    def bounds(self, x: ArrayLike, y: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        """Return lower and upper parameter bounds.

        Parameters
        ----------
        x, y : array-like
            Coordinates and observations.

        Returns
        -------
        tuple of ndarray
            Lower and upper bounds. The default is unbounded.
        """
        n_parameters = len(self.parameter_names)
        return (
            np.full(n_parameters, -np.inf, dtype=np.float64),
            np.full(n_parameters, np.inf, dtype=np.float64),
        )

    def derivative_x(
        self,
        x: ArrayLike,
        parameters: ArrayLike,
    ) -> np.ndarray:
        """Evaluate the derivative of the model with respect to ``x``.

        The default implementation uses centered finite differences with a
        scale-aware step. Physical models should override this method when an
        analytical derivative is available, particularly when the derivative
        enters an uncertainty model such as effective variance.

        Parameters
        ----------
        x : array-like
            Independent coordinates.
        parameters : array-like
            Model parameters in ``parameter_names`` order.

        Returns
        -------
        ndarray
            Values of :math:`\\partial y/\\partial x` at each coordinate.

        Raises
        ------
        ValueError
            If the numerical derivative is non-finite.
        """
        coordinates = np.asarray(x, dtype=np.float64)
        values = np.asarray(parameters, dtype=np.float64)
        if coordinates.ndim not in {1, 2}:
            raise ValueError(
                "x must be a vector or coordinate matrix for derivative evaluation"
            )
        if coordinates.ndim == 1:
            step = np.cbrt(np.finfo(np.float64).eps) * np.maximum(
                np.abs(coordinates), 1.0
            )
            upper = self.evaluate(coordinates + step, values)
            lower = self.evaluate(coordinates - step, values)
            derivative = (upper - lower) / (2.0 * step)
        else:
            derivative = np.empty_like(coordinates)
            for axis in range(coordinates.shape[0]):
                step = np.cbrt(np.finfo(np.float64).eps) * np.maximum(
                    np.abs(coordinates[axis]), 1.0
                )
                upper_coordinates = coordinates.copy()
                lower_coordinates = coordinates.copy()
                upper_coordinates[axis] += step
                lower_coordinates[axis] -= step
                upper = self.evaluate(upper_coordinates, values)
                lower = self.evaluate(lower_coordinates, values)
                derivative[axis] = (upper - lower) / (2.0 * step)
        if derivative.shape != coordinates.shape or not np.all(np.isfinite(derivative)):
            raise ValueError("model derivative with respect to x is not finite")
        return np.asarray(derivative, dtype=np.float64)

    def curve_function(self) -> ModelFunction:
        """Return a ``curve_fit`` compatible callable.

        Returns
        -------
        callable
            Function with signature ``f(x, *parameters)``.
        """

        def function(x: np.ndarray, *parameters: float) -> np.ndarray:
            return self.evaluate(x, np.asarray(parameters, dtype=np.float64))

        return function

    def metadata(self) -> dict[str, Any]:
        """Return metadata stored with fit results.

        Returns
        -------
        dict
            Model name and parameter order.
        """
        return {"model": self.name, "parameter_order": list(self.parameter_names)}


class CallableFitModel(BaseFitModel):
    """Adapt a callable function to the general fit-model contract.

    Parameters
    ----------
    function : callable
        Function with signature ``f(x, *parameters)``.
    initial_parameters : sequence of float
        Default initial parameters.
    name : str, optional
        Model identifier.
    parameter_names : sequence of str, optional
        Ordered parameter labels. Generic labels are generated when omitted.
    bounds : tuple of array-like, optional
        Default lower and upper bounds.
    """

    def __init__(
        self,
        function: ModelFunction,
        initial_parameters: Sequence[float],
        *,
        name: str = "callable",
        parameter_names: Sequence[str] | None = None,
        bounds: tuple[ArrayLike | float, ArrayLike | float] = (-np.inf, np.inf),
    ) -> None:
        self._function = function
        self._initial_parameters = np.asarray(initial_parameters, dtype=np.float64)
        if self._initial_parameters.ndim != 1 or self._initial_parameters.size == 0:
            raise ValueError(
                "initial_parameters must be a non-empty one-dimensional array"
            )
        self._name = str(name)
        if parameter_names is None:
            self._parameter_names = tuple(
                f"p{i}" for i in range(self._initial_parameters.size)
            )
        else:
            labels = tuple(str(item) for item in parameter_names)
            if len(labels) != self._initial_parameters.size:
                raise ValueError("parameter_names must match initial_parameters")
            self._parameter_names = labels
        lower, upper = bounds
        self._bounds = (
            np.broadcast_to(
                np.asarray(lower, dtype=np.float64), self._initial_parameters.shape
            ).copy(),
            np.broadcast_to(
                np.asarray(upper, dtype=np.float64), self._initial_parameters.shape
            ).copy(),
        )

    @property
    def name(self) -> str:
        """Return the model name.

        Returns
        -------
        str
            Model identifier.
        """
        return self._name

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return parameter names.

        Returns
        -------
        tuple of str
            Ordered labels.
        """
        return self._parameter_names

    def evaluate(self, x: ArrayLike, parameters: ArrayLike) -> np.ndarray:
        """Evaluate the wrapped callable.

        Parameters
        ----------
        x : array-like
            Independent coordinates.
        parameters : array-like
            Model parameters.

        Returns
        -------
        ndarray
            Model values.
        """
        values = np.asarray(parameters, dtype=np.float64)
        return np.asarray(
            self._function(np.asarray(x, dtype=np.float64), *values), dtype=np.float64
        )

    def initial_guess(self, x: ArrayLike, y: ArrayLike) -> np.ndarray:
        """Return the configured initial parameters.

        Parameters
        ----------
        x, y : array-like
            Coordinates and observations, validated for consistency.

        Returns
        -------
        ndarray
            Initial parameter estimate.
        """
        validate_xy(x, y)
        return self._initial_parameters.copy()

    def bounds(self, x: ArrayLike, y: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        """Return configured parameter bounds.

        Parameters
        ----------
        x, y : array-like
            Coordinates and observations, validated for consistency.

        Returns
        -------
        tuple of ndarray
            Lower and upper bounds.
        """
        validate_xy(x, y)
        return self._bounds[0].copy(), self._bounds[1].copy()

    def curve_function(self) -> ModelFunction:
        """Return the wrapped callable.

        Returns
        -------
        callable
            Original ``curve_fit`` compatible function.
        """
        return self._function


class FittedModel:
    """Combine a fit model with optimized parameters.

    Parameters
    ----------
    model : BaseFitModel
        Mathematical or physical model.
    fit_result : FitResult
        Successful fit result containing optimized parameters.

    Raises
    ------
    ValueError
        If the fit result does not contain usable parameters.
    """

    def __init__(self, model: BaseFitModel, fit_result: FitResult) -> None:
        if not fit_result.success or fit_result.parameters is None:
            raise ValueError("fit_result must contain successful optimized parameters")
        if fit_result.parameters.size != len(model.parameter_names):
            raise ValueError("fit parameters do not match the model definition")
        self.model = model
        self.fit_result = fit_result

    @property
    def parameters(self) -> np.ndarray:
        """Return optimized parameters.

        Returns
        -------
        ndarray
            Copy of the fitted parameter vector.
        """
        assert self.fit_result.parameters is not None
        return self.fit_result.parameters.copy()

    def evaluate(self, x: ArrayLike) -> np.ndarray:
        """Evaluate the fitted model.

        Parameters
        ----------
        x : array-like
            Independent coordinates.

        Returns
        -------
        ndarray
            Fitted model values.
        """
        return self.model.evaluate(x, self.parameters)


class MappedFitModel(BaseFitModel):
    """Expose only free parameters of a complete physical model to a solver.

    Parameters
    ----------
    model : BaseFitModel
        Model whose parameter order corresponds to all non-derived definitions
        in ``parameter_map``.
    parameter_map : ParameterMap
        Reduced-to-complete parameter mapping.

    Raises
    ------
    ValueError
        If model parameter names do not match the non-derived mapping order.
    """

    def __init__(self, model: BaseFitModel, parameter_map: ParameterMap) -> None:
        model_names = tuple(
            definition.name
            for definition in parameter_map.definitions
            if definition.state is not ParameterState.DERIVED
        )
        if tuple(model.parameter_names) != model_names:
            raise ValueError(
                "model parameter names must match non-derived ParameterMap order"
            )
        self.model = model
        self.parameter_map = parameter_map

    @property
    def name(self) -> str:
        """Return the wrapped model name."""
        return self.model.name

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return only the names optimized by the solver."""
        return self.parameter_map.free_names

    def evaluate(self, x: ArrayLike, parameters: ArrayLike) -> np.ndarray:
        """Resolve free values and evaluate the complete physical model."""
        resolved = self.parameter_map.expand(parameters)
        return self.model.evaluate(x, resolved.model_values())

    def derivative_x(
        self,
        x: ArrayLike,
        parameters: ArrayLike,
    ) -> np.ndarray:
        """Resolve free values and evaluate the wrapped x derivative."""
        resolved = self.parameter_map.expand(parameters)
        return self.model.derivative_x(x, resolved.model_values())

    def initial_guess(self, x: ArrayLike, y: ArrayLike) -> np.ndarray:
        """Return the configured reduced initial vector."""
        validate_xy(x, y)
        return self.parameter_map.initial_free_values()

    def bounds(self, x: ArrayLike, y: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        """Return configured reduced parameter bounds."""
        validate_xy(x, y)
        return self.parameter_map.free_bounds()

    def metadata(self) -> dict[str, Any]:
        """Return wrapped model and complete mapping metadata."""
        metadata = self.model.metadata()
        metadata["optimizer_parameter_order"] = list(self.parameter_names)
        metadata["parameter_map"] = self.parameter_map.as_dict()
        return metadata
