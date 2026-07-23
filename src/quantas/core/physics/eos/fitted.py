# -*- coding: utf-8 -*-

"""Reusable fitted energy equations of state and uncertainty propagation."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, TypeAlias

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

from .energy import EnergyEOS
from .parameters import (
    EOSParameters,
    resolve_energy_parameters,
    resolved_energy_parameter_covariance,
)
from .pressure import PressureEOS
from .spec import EOSModel
from .state import EOSState, EOSStateUncertainty

ArrayLike: TypeAlias = np.ndarray | float | Sequence[float]
RandomState: TypeAlias = int | np.random.Generator | None
_STATE_ORDER = ("volume", "bulk_modulus", "bulk_modulus_derivative")


class FittedEnergyEOS:
    """Evaluate a fitted integrated EOS over an arbitrary pressure grid.

    The stored parameter vector contains only the free parameters of the
    selected EOS order.  Fixed and implied parameters are resolved every time
    a state is evaluated, including during covariance or Monte Carlo
    propagation.
    """

    def __init__(
        self,
        eos: str | EOSModel,
        parameters: ArrayLike,
        *,
        sampled_volumes: ArrayLike | None = None,
        covariance: np.ndarray | None = None,
    ) -> None:
        """Initialize a fitted equation-of-state model."""
        self._energy_eos = EnergyEOS()
        self._pressure_eos = PressureEOS()
        self._model = self._energy_eos.model(eos)
        self._parameter_names = self._model.energy_parameter_names
        self._parameters = self._validate_parameters(parameters)
        self._sampled_volumes = self._validate_sampled_volumes(sampled_volumes)
        self._covariance = self._validate_covariance(covariance)

    @property
    def eos(self) -> str:
        """Return the canonical EOS family name."""
        return self._model.family.value

    @property
    def eos_tag(self) -> str:
        """Return the compact family-and-order tag."""
        return self._model.tag

    @property
    def eos_order(self) -> int | None:
        """Return the EOS order."""
        return self._model.order

    @property
    def model(self) -> EOSModel:
        """Return the immutable EOS model specification."""
        return self._model

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return free parameter names in fitting order."""
        return self._parameter_names

    @property
    def parameters(self) -> np.ndarray:
        """Return a copy of the free fitted parameters."""
        return self._parameters.copy()

    @property
    def resolved_parameters(self) -> EOSParameters:
        """Return the complete physical parameter set."""
        return resolve_energy_parameters(self._model, self._parameters)

    @property
    def covariance(self) -> np.ndarray | None:
        """Return a copy of the free-parameter covariance matrix."""
        return None if self._covariance is None else self._covariance.copy()

    @property
    def resolved_covariance(self) -> np.ndarray | None:
        """Return covariance of ``E0, K0, KP, KPP, V0``.

        Returns
        -------
        ndarray or None
            Complete physical-parameter covariance.  Parameters fixed by the
            selected EOS order have zero variance, while uncertainties of
            implied parameters are propagated from the free fit covariance.
        """
        if self._covariance is None:
            return None
        return resolved_energy_parameter_covariance(
            self._model,
            self._parameters,
            self._covariance,
        )

    @property
    def resolved_errors(self) -> np.ndarray | None:
        """Return one-sigma errors of ``E0, K0, KP, KPP, V0``."""
        covariance = self.resolved_covariance
        if covariance is None:
            return None
        return np.sqrt(np.clip(np.diag(covariance), 0.0, None))

    @property
    def sampled_volume_range(self) -> tuple[float, float] | None:
        """Return the sampled volume interval."""
        if self._sampled_volumes is None:
            return None
        return float(self._sampled_volumes[0]), float(self._sampled_volumes[-1])

    @property
    def reference_volume(self) -> float:
        """Return the zero-pressure reference volume."""
        return float(self.resolved_parameters.V0)

    @property
    def reference_bulk_modulus(self) -> float:
        """Return the zero-pressure bulk modulus."""
        return float(self.resolved_parameters.K0)

    @property
    def reference_bulk_modulus_derivative(self) -> float:
        """Return the first pressure derivative at the reference state."""
        return float(self.resolved_parameters.KP)

    @property
    def reference_bulk_modulus_second_derivative(self) -> float:
        """Return the second pressure derivative at the reference state."""
        return float(self.resolved_parameters.KPP)

    def pressure(self, volume: ArrayLike) -> np.ndarray:
        """Evaluate pressure at one or more volumes."""
        return self._pressure_eos.pressure(
            self._model, self.resolved_parameters, volume
        )

    def bulk_modulus(self, volume: ArrayLike) -> np.ndarray:
        """Evaluate the isothermal bulk modulus."""
        return self._pressure_eos.bulk_modulus(
            self._model, self.resolved_parameters, volume
        )

    def bulk_modulus_derivative(self, volume: ArrayLike) -> np.ndarray:
        """Evaluate the first pressure derivative of the bulk modulus."""
        return self._pressure_eos.bulk_modulus_derivative(
            self._model, self.resolved_parameters, volume
        )

    def bulk_modulus_second_derivative(self, volume: ArrayLike) -> np.ndarray:
        """Evaluate the second pressure derivative of the bulk modulus."""
        return self._pressure_eos.bulk_modulus_second_derivative(
            self._model, self.resolved_parameters, volume
        )

    def volume_at_pressure(
        self,
        pressure: float,
        *,
        bounds: tuple[float, float] | None = None,
        xtol: float = 1.0e-12,
        rtol: float = 1.0e-12,
        maxiter: int = 200,
    ) -> float:
        """Return the positive volume satisfying ``P(V)=pressure``."""
        return self._volume_at_pressure_for_parameters(
            pressure,
            self._parameters,
            bounds=bounds,
            xtol=xtol,
            rtol=rtol,
            maxiter=maxiter,
        )

    def state_at_pressure(
        self,
        pressure: float,
        *,
        bounds: tuple[float, float] | None = None,
        uncertainty_method: str = "none",
        relative_step: float = 1.0e-5,
        confidence_level: float = 0.95,
        monte_carlo_samples: int = 10000,
        random_state: RandomState = None,
        minimum_success_fraction: float = 0.8,
    ) -> EOSState:
        """Evaluate volume, bulk modulus and ``K'`` at a target pressure."""
        vector = self._state_vector_at_pressure(
            pressure, self._parameters, bounds=bounds
        )
        method = self._canonical_uncertainty_method(uncertainty_method)
        uncertainty = None
        if method == "covariance":
            uncertainty = self._covariance_uncertainty(
                pressure,
                vector,
                bounds=bounds,
                relative_step=relative_step,
                confidence_level=confidence_level,
            )
        elif method == "montecarlo":
            uncertainty = self._monte_carlo_uncertainty(
                pressure,
                bounds=bounds,
                confidence_level=confidence_level,
                n_samples=monte_carlo_samples,
                random_state=random_state,
                minimum_success_fraction=minimum_success_fraction,
            )
        volume, bulk_modulus, derivative = vector
        return EOSState(
            pressure=float(pressure),
            volume=float(volume),
            bulk_modulus=float(bulk_modulus),
            bulk_modulus_derivative=float(derivative),
            extrapolated=self._is_extrapolated(float(volume)),
            uncertainty=uncertainty,
            metadata={
                "eos": self.eos,
                "eos_tag": self.eos_tag,
                "eos_order": self.eos_order,
            },
        )

    def states_at_pressures(
        self,
        pressures: ArrayLike,
        *,
        uncertainty_method: str = "none",
        relative_step: float = 1.0e-5,
        confidence_level: float = 0.95,
        monte_carlo_samples: int = 10000,
        random_state: RandomState = None,
        minimum_success_fraction: float = 0.8,
    ) -> list[EOSState]:
        """Evaluate EOS states at several pressures."""
        values = np.asarray(pressures, dtype=np.float64)
        if values.ndim != 1:
            raise ValueError("pressures must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("pressures must be finite")
        return [
            self.state_at_pressure(
                float(value),
                uncertainty_method=uncertainty_method,
                relative_step=relative_step,
                confidence_level=confidence_level,
                monte_carlo_samples=monte_carlo_samples,
                random_state=random_state,
                minimum_success_fraction=minimum_success_fraction,
            )
            for value in values
        ]

    def state_jacobian_at_pressure(
        self,
        pressure: float,
        *,
        bounds: tuple[float, float] | None = None,
        relative_step: float = 1.0e-5,
    ) -> np.ndarray:
        """Return the state Jacobian with respect to free EOS parameters."""
        if not np.isfinite(relative_step) or relative_step <= 0.0:
            raise ValueError("relative_step must be finite and positive")
        jacobian = np.zeros((3, len(self._parameter_names)), dtype=np.float64)
        for index, name in enumerate(self._parameter_names):
            if name == "E0":
                continue
            step = self._parameter_step(index, relative_step)
            plus_two = self._parameters.copy()
            plus_one = self._parameters.copy()
            minus_one = self._parameters.copy()
            minus_two = self._parameters.copy()
            plus_two[index] += 2.0 * step
            plus_one[index] += step
            minus_one[index] -= step
            minus_two[index] -= 2.0 * step
            if name in {"K0", "V0"} and minus_two[index] <= 0.0:
                f0 = self._state_vector_at_pressure(
                    pressure, self._parameters, bounds=bounds
                )
                f1 = self._state_vector_at_pressure(pressure, plus_one, bounds=bounds)
                f2 = self._state_vector_at_pressure(pressure, plus_two, bounds=bounds)
                jacobian[:, index] = (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * step)
                continue
            f2p = self._state_vector_at_pressure(pressure, plus_two, bounds=bounds)
            f1p = self._state_vector_at_pressure(pressure, plus_one, bounds=bounds)
            f1m = self._state_vector_at_pressure(pressure, minus_one, bounds=bounds)
            f2m = self._state_vector_at_pressure(pressure, minus_two, bounds=bounds)
            jacobian[:, index] = (-f2p + 8.0 * f1p - 8.0 * f1m + f2m) / (12.0 * step)
        return jacobian

    def state_covariance_at_pressure(
        self,
        pressure: float,
        *,
        bounds: tuple[float, float] | None = None,
        relative_step: float = 1.0e-5,
    ) -> np.ndarray:
        """Propagate the free-parameter covariance to one pressure state."""
        covariance = self._require_covariance()
        jacobian = self.state_jacobian_at_pressure(
            pressure, bounds=bounds, relative_step=relative_step
        )
        result = jacobian @ covariance @ jacobian.T
        return self._stabilize_covariance(result, name="state covariance")

    def _validate_parameters(self, parameters: ArrayLike) -> np.ndarray:
        values = np.asarray(parameters, dtype=np.float64)
        if values.shape != (len(self._parameter_names),):
            names = ", ".join(self._parameter_names)
            raise ValueError(f"parameters for {self.eos_tag} must contain {names}")
        if not np.all(np.isfinite(values)):
            raise ValueError("EOS parameters must be finite")
        resolve_energy_parameters(self._model, values)
        return values.copy()

    @staticmethod
    def _validate_sampled_volumes(
        sampled_volumes: ArrayLike | None,
    ) -> np.ndarray | None:
        if sampled_volumes is None:
            return None
        values = np.asarray(sampled_volumes, dtype=np.float64)
        if values.ndim != 1 or values.size < 2:
            raise ValueError(
                "sampled volumes must be a one-dimensional array with at least two values"
            )
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("sampled volumes must be finite and positive")
        return np.sort(values)

    def _validate_covariance(self, covariance: np.ndarray | None) -> np.ndarray | None:
        if covariance is None:
            return None
        values = np.asarray(covariance, dtype=np.float64)
        expected = (len(self._parameter_names), len(self._parameter_names))
        if values.shape != expected:
            raise ValueError(f"covariance must have shape {expected}")
        if not np.all(np.isfinite(values)):
            raise ValueError("covariance must contain only finite values")
        scale = max(float(np.max(np.abs(values))), 1.0)
        if not np.allclose(values, values.T, rtol=1.0e-8, atol=1.0e-12 * scale):
            raise ValueError("covariance must be symmetric")
        return self._stabilize_covariance(0.5 * (values + values.T), name="covariance")

    @staticmethod
    def _stabilize_covariance(covariance: np.ndarray, *, name: str) -> np.ndarray:
        values = np.asarray(covariance, dtype=np.float64)
        values = 0.5 * (values + values.T)
        eigenvalues, eigenvectors = np.linalg.eigh(values)
        scale = max(
            float(np.max(np.abs(eigenvalues))), float(np.finfo(np.float64).tiny)
        )
        tolerance = max(1.0e-14, 1.0e-10 * scale)
        if float(np.min(eigenvalues)) < -tolerance:
            raise ValueError(f"{name} must be positive semidefinite")
        clipped = np.clip(eigenvalues, 0.0, None)
        stabilized = (eigenvectors * clipped) @ eigenvectors.T
        return 0.5 * (stabilized + stabilized.T)

    def _require_covariance(self) -> np.ndarray:
        if self._covariance is None:
            raise ValueError(
                "parameter covariance is required for uncertainty propagation"
            )
        return self._covariance

    @staticmethod
    def _canonical_uncertainty_method(method: str) -> str:
        aliases = {
            "": "none",
            "none": "none",
            "off": "none",
            "covariance": "covariance",
            "linear": "covariance",
            "delta": "covariance",
            "montecarlo": "montecarlo",
            "monte-carlo": "montecarlo",
            "mc": "montecarlo",
        }
        key = str(method).strip().lower()
        try:
            return aliases[key]
        except KeyError as exc:
            raise ValueError(
                "uncertainty_method must be none, covariance, or montecarlo"
            ) from exc

    def _covariance_uncertainty(
        self,
        pressure: float,
        state: np.ndarray,
        *,
        bounds: tuple[float, float] | None,
        relative_step: float,
        confidence_level: float,
    ) -> EOSStateUncertainty:
        self._validate_confidence_level(confidence_level)
        covariance = self.state_covariance_at_pressure(
            pressure, bounds=bounds, relative_step=relative_step
        )
        standard = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
        z_value = float(norm.ppf(0.5 + 0.5 * confidence_level))
        intervals = np.column_stack(
            (state - z_value * standard, state + z_value * standard)
        )
        return EOSStateUncertainty(
            method="covariance",
            mean=state,
            covariance=covariance,
            correlation=self._covariance_to_correlation(covariance),
            confidence_intervals=intervals,
            metadata={
                "confidence_level": float(confidence_level),
                "relative_step": float(relative_step),
                "parameter_order": list(self._parameter_names),
                "state_order": list(_STATE_ORDER),
                "eos_tag": self.eos_tag,
            },
        )

    def _monte_carlo_uncertainty(
        self,
        pressure: float,
        *,
        bounds: tuple[float, float] | None,
        confidence_level: float,
        n_samples: int,
        random_state: RandomState,
        minimum_success_fraction: float,
    ) -> EOSStateUncertainty:
        covariance = self._require_covariance()
        self._validate_confidence_level(confidence_level)
        if int(n_samples) != n_samples or n_samples < 2:
            raise ValueError("monte_carlo_samples must be an integer greater than one")
        if (
            not np.isfinite(minimum_success_fraction)
            or minimum_success_fraction <= 0.0
            or minimum_success_fraction > 1.0
        ):
            raise ValueError("minimum_success_fraction must lie in (0, 1]")
        rng = (
            random_state
            if isinstance(random_state, np.random.Generator)
            else np.random.default_rng(random_state)
        )
        samples = rng.multivariate_normal(
            self._parameters,
            covariance,
            size=int(n_samples),
            check_valid="raise",
        )
        states: list[np.ndarray] = []
        nonphysical = 0
        failed = 0
        for parameters in samples:
            try:
                resolve_energy_parameters(self._model, parameters)
            except ValueError:
                nonphysical += 1
                continue
            try:
                state = self._state_vector_at_pressure(
                    pressure, parameters, bounds=bounds
                )
            except (ValueError, FloatingPointError, OverflowError):
                failed += 1
                continue
            if not np.all(np.isfinite(state)):
                failed += 1
                continue
            states.append(state)
        n_successful = len(states)
        required = max(2, int(np.ceil(float(n_samples) * minimum_success_fraction)))
        if n_successful < required:
            raise ValueError(
                "Monte Carlo uncertainty propagation produced too few physical EOS states "
                f"({n_successful}/{n_samples}; required {required})"
            )
        values = np.asarray(states, dtype=np.float64)
        mean = np.mean(values, axis=0)
        covariance_state = np.cov(values, rowvar=False, ddof=1)
        covariance_state = self._stabilize_covariance(
            covariance_state, name="Monte Carlo state covariance"
        )
        tail = 0.5 * (1.0 - confidence_level)
        intervals = np.quantile(values, [tail, 1.0 - tail], axis=0).T
        return EOSStateUncertainty(
            method="montecarlo",
            mean=mean,
            covariance=covariance_state,
            correlation=self._covariance_to_correlation(covariance_state),
            confidence_intervals=intervals,
            n_samples=int(n_samples),
            n_successful=n_successful,
            metadata={
                "confidence_level": float(confidence_level),
                "minimum_success_fraction": float(minimum_success_fraction),
                "nonphysical_parameter_samples": nonphysical,
                "failed_state_samples": failed,
                "parameter_order": list(self._parameter_names),
                "state_order": list(_STATE_ORDER),
                "eos_tag": self.eos_tag,
            },
        )

    @staticmethod
    def _covariance_to_correlation(covariance: np.ndarray) -> np.ndarray:
        standard = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
        denominator = np.outer(standard, standard)
        correlation = np.zeros_like(covariance, dtype=np.float64)
        np.divide(covariance, denominator, out=correlation, where=denominator > 0.0)
        correlation[np.diag_indices_from(correlation)] = (standard > 0.0).astype(float)
        return np.clip(0.5 * (correlation + correlation.T), -1.0, 1.0)

    @staticmethod
    def _validate_confidence_level(confidence_level: float) -> None:
        if (
            not np.isfinite(confidence_level)
            or confidence_level <= 0.0
            or confidence_level >= 1.0
        ):
            raise ValueError("confidence_level must lie in (0, 1)")

    def _parameter_step(self, index: int, relative_step: float) -> float:
        value = float(self._parameters[index])
        scale = max(abs(value), 1.0)
        step = max(
            relative_step * scale,
            np.sqrt(np.finfo(np.float64).eps) * scale,
        )
        if self._parameter_names[index] in {"K0", "V0"}:
            step = min(step, 0.2 * value)
        return float(step)

    def _state_vector_at_pressure(
        self,
        pressure: float,
        parameters: np.ndarray,
        *,
        bounds: tuple[float, float] | None,
    ) -> np.ndarray:
        values = self._validate_parameters(parameters)
        resolved = resolve_energy_parameters(self._model, values)
        volume = self._volume_at_pressure_for_parameters(
            pressure, values, bounds=bounds
        )
        bulk_modulus = float(
            np.asarray(self._pressure_eos.bulk_modulus(self._model, resolved, volume))
        )
        derivative = float(
            np.asarray(
                self._pressure_eos.bulk_modulus_derivative(
                    self._model, resolved, volume
                )
            )
        )
        state = np.asarray([volume, bulk_modulus, derivative], dtype=np.float64)
        if not np.all(np.isfinite(state)):
            raise ValueError("EOS state contains non-finite values")
        return state

    def _volume_at_pressure_for_parameters(
        self,
        pressure: float,
        parameters: np.ndarray,
        *,
        bounds: tuple[float, float] | None,
        xtol: float = 1.0e-12,
        rtol: float = 1.0e-12,
        maxiter: int = 200,
    ) -> float:
        target = float(pressure)
        if not np.isfinite(target):
            raise ValueError("pressure must be finite")
        values = self._validate_parameters(parameters)
        resolved = resolve_energy_parameters(self._model, values)
        lower, upper = self._initial_bounds(bounds, resolved.V0)

        def residual(volume: float) -> float:
            value = self._pressure_eos.pressure(self._model, resolved, volume)
            return float(np.asarray(value)) - target

        lower, upper = self._expand_bracket(residual, lower, upper)
        return float(
            brentq(residual, lower, upper, xtol=xtol, rtol=rtol, maxiter=maxiter)
        )

    def _initial_bounds(
        self,
        bounds: tuple[float, float] | None,
        reference_volume: float,
    ) -> tuple[float, float]:
        if bounds is not None:
            lower, upper = map(float, bounds)
        elif self._sampled_volumes is not None:
            lower = float(self._sampled_volumes[0])
            upper = float(self._sampled_volumes[-1])
        else:
            lower = 0.75 * reference_volume
            upper = 1.25 * reference_volume
        if not np.isfinite(lower) or not np.isfinite(upper):
            raise ValueError("volume bounds must be finite")
        if lower <= 0.0 or upper <= lower:
            raise ValueError("volume bounds must satisfy 0 < lower < upper")
        return lower, upper

    @staticmethod
    def _expand_bracket(
        residual: Any,
        lower: float,
        upper: float,
        *,
        max_expansions: int = 100,
    ) -> tuple[float, float]:
        f_lower = float(residual(lower))
        f_upper = float(residual(upper))
        if not np.isfinite(f_lower) or not np.isfinite(f_upper):
            raise ValueError(
                "EOS returned non-finite pressure while bracketing the volume"
            )
        if f_lower == 0.0:
            return lower, lower * (1.0 + 1.0e-12)
        if f_upper == 0.0:
            return upper * (1.0 - 1.0e-12), upper
        minimum = float(np.finfo(np.float64).tiny)
        for _ in range(max_expansions):
            if f_lower * f_upper < 0.0:
                return lower, upper
            lower = max(lower / 1.5, minimum)
            upper *= 1.5
            f_lower = float(residual(lower))
            f_upper = float(residual(upper))
            if not np.isfinite(f_lower) or not np.isfinite(f_upper):
                continue
        raise ValueError("the requested pressure state could not be bracketed")

    def _is_extrapolated(self, volume: float) -> bool:
        if self._sampled_volumes is None:
            return False
        return bool(
            volume < float(self._sampled_volumes[0])
            or volume > float(self._sampled_volumes[-1])
        )


__all__ = ["FittedEnergyEOS"]
