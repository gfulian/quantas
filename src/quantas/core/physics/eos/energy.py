# -*- coding: utf-8 -*-

"""Volume-integrated equations of state and energy fitting models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Callable, TypeAlias

import numpy as np

from quantas.core.math.fitting import (
    BaseFitModel,
    FitQuality,
    FitResult,
    FitStatus,
    LeastSquaresFitter,
    OLSOptions,
    validate_xy,
)

from .parameters import (
    EOSParameters,
    free_energy_parameters,
    resolve_energy_parameters,
    resolved_energy_parameter_covariance,
)
from .pressure import PressureEOS
from .spec import EOSFamily, EOSModel, parse_eos_model

ArrayLike: TypeAlias = np.ndarray | float | Sequence[float]
ParameterArrayLike: TypeAlias = np.ndarray | Sequence[float]
ParameterLike: TypeAlias = ParameterArrayLike | Mapping[str, float] | EOSParameters


class EnergyEOS:
    r"""Evaluate and fit volume-integrated equations of state.

    For every supported model the pressure counterpart is obtained from

    .. math::

        P(V)=-\frac{\mathrm dE}{\mathrm dV}.

    The free fit parameters depend on the EOS order.  Parameters omitted from
    a lower-order fit are resolved from the corresponding truncation rules.
    """

    def model(self, eos: str | EOSModel, order: int | None = None) -> EOSModel:
        """Return a canonical family-and-order specification."""
        model = parse_eos_model(eos, order)
        if not model.supports_energy:
            raise ValueError(f"{model.tag} has no implemented energy-volume form")
        return model

    def canonical_name(self, eos: str | EOSModel, order: int | None = None) -> str:
        """Return the canonical EOS family name."""
        return self.model(eos, order).family.value

    def canonical_tag(self, eos: str | EOSModel, order: int | None = None) -> str:
        """Return the compact canonical EOS tag."""
        return self.model(eos, order).tag

    _canonical_name = canonical_name

    def function(
        self,
        eos: str | EOSModel,
        order: int | None = None,
    ) -> Callable[..., np.ndarray]:
        """Return a ``curve_fit`` compatible energy function."""
        return EnergyEOSFitModel(self, self.model(eos, order)).curve_function()

    def evaluate(
        self,
        eos: str | EOSModel,
        volume: ArrayLike,
        parameters: ParameterLike,
        *,
        order: int | None = None,
    ) -> np.ndarray:
        """Evaluate an integrated energy EOS."""
        model = self.model(eos, order)
        pars = resolve_energy_parameters(model, parameters)
        values = self._validate_volume(volume)
        if model.family is EOSFamily.MURNAGHAN:
            return self._murnaghan_energy(values, pars)
        if model.family is EOSFamily.BIRCH_MURNAGHAN:
            if model.order == 3:
                return self._birch_murnaghan_third_order_energy(values, pars)
            return self._strain_energy(values, pars, EOSFamily.BIRCH_MURNAGHAN)
        if model.family is EOSFamily.NATURAL_STRAIN:
            return self._strain_energy(values, pars, EOSFamily.NATURAL_STRAIN)
        if model.family is EOSFamily.VINET:
            return self._vinet_energy(values, pars, model.order)
        raise ValueError(f"unknown energy EOS: {eos!r}")

    def fit(
        self,
        eos: str | EOSModel,
        volume: ParameterArrayLike,
        energy: ParameterArrayLike,
        *,
        order: int | None = None,
        p0: Sequence[float] | Mapping[str, float] | EOSParameters | None = None,
        maxfev: int | None = None,
    ) -> FitResult:
        """Fit energy-volume data with a selected EOS family and order."""
        try:
            model_spec = self.model(eos, order)
        except ValueError as exc:
            return FitResult.failed(str(exc), status=FitStatus.INVALID_INPUT)
        metadata = self._fit_metadata(model_spec)
        try:
            volume_array, energy_array = validate_xy(volume, energy)
        except ValueError as exc:
            return FitResult.failed(
                str(exc), status=FitStatus.INVALID_INPUT, metadata=metadata
            )

        fit_model = EnergyEOSFitModel(self, model_spec)
        try:
            if p0 is None:
                initial = fit_model.initial_guess(volume_array, energy_array)
            else:
                resolved = resolve_energy_parameters(model_spec, p0)
                initial = free_energy_parameters(model_spec, resolved)
        except Exception as exc:  # noqa: BLE001 - initial estimate diagnostics
            return FitResult.failed(
                f"could not estimate initial EOS parameters: {exc}",
                n_points=int(volume_array.size),
                n_parameters=len(model_spec.energy_parameter_names),
                metadata=metadata,
            )

        result = LeastSquaresFitter().fit(
            fit_model,
            volume_array,
            energy_array,
            initial_parameters=initial,
            options=OLSOptions(
                max_iterations=maxfev,
                ftol=1.0e-15,
                xtol=1.0e-15,
                metadata=metadata,
            ),
        )
        if result.success and result.parameters is not None:
            resolved = resolve_energy_parameters(model_spec, result.parameters)
            resolved_metadata: dict[str, object] = {
                "resolved_parameter_order": ["E0", "K0", "KP", "KPP", "V0"],
                "resolved_parameters": [
                    resolved.E0,
                    resolved.K0,
                    resolved.KP,
                    resolved.KPP,
                    resolved.V0,
                ],
                "parameter_sources": {
                    "E0": "fitted",
                    **model_spec.parameter_sources,
                },
            }
            if result.covariance is not None:
                covariance = np.asarray(result.covariance, dtype=np.float64)
                covariance_status: str | None = None
                covariance_warning: str | None = None
                if result.dof == 0:
                    covariance_status = "unavailable_zero_degrees_of_freedom"
                    covariance_warning = (
                        f"the {model_spec.tag} fit is exactly determined: "
                        f"{result.n_points} observations and "
                        f"{result.n_parameters} free parameters give zero "
                        "residual degrees of freedom; parameter covariance and "
                        "statistical uncertainties are unavailable"
                    )
                elif not np.all(np.isfinite(covariance)):
                    covariance_status = "unavailable_non_finite"
                    covariance_warning = (
                        "the optimizer returned a non-finite EOS covariance; "
                        "parameter uncertainties are unavailable"
                    )

                if covariance_status is None:
                    resolved_covariance = resolved_energy_parameter_covariance(
                        model_spec,
                        result.parameters,
                        covariance,
                    )
                    resolved_metadata["resolved_covariance"] = (
                        resolved_covariance.tolist()
                    )
                    resolved_metadata["resolved_errors"] = np.sqrt(
                        np.clip(np.diag(resolved_covariance), 0.0, None)
                    ).tolist()
                else:
                    assert covariance_warning is not None
                    result.covariance = None
                    result.errors = None
                    result.optimizer_covariance = None
                    result.optimizer_errors = None
                    result.condition_number = None
                    result.quality = FitQuality.POOR
                    result.message = "fit converged with diagnostic warnings"
                    if covariance_warning not in result.warnings:
                        result.warnings.append(covariance_warning)
                    if result.diagnostics is not None:
                        result.diagnostics.correlation = None
                        if covariance_warning not in result.diagnostics.warnings:
                            result.diagnostics.warnings.append(covariance_warning)
                        result.diagnostics.metadata["covariance_status"] = (
                            covariance_status
                        )
                    resolved_metadata["covariance_status"] = covariance_status
            result.metadata.update(resolved_metadata)
        return result

    def fit_parameters(
        self,
        eos: str | EOSModel,
        volume: ParameterArrayLike,
        energy: ParameterArrayLike,
        *,
        order: int | None = None,
        p0: Sequence[float] | Mapping[str, float] | EOSParameters | None = None,
        maxfev: int | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return optimized free parameters and one-sigma errors."""
        result = self.fit(eos, volume, energy, order=order, p0=p0, maxfev=maxfev)
        if not result.success or result.parameters is None:
            raise RuntimeError(result.message)
        errors = (
            result.errors
            if result.errors is not None
            else np.full(result.parameters.shape, np.nan)
        )
        return result.parameters, errors

    def guess(
        self, volume: ParameterArrayLike, energy: ParameterArrayLike
    ) -> EOSParameters:
        """Estimate a complete physical parameter set from a polynomial fit."""
        volume_array, energy_array = validate_xy(volume, energy)
        if volume_array.size < 4:
            raise ValueError(
                "at least four points are required for the EOS initial estimate"
            )
        degree = 4 if volume_array.size >= 5 else 2
        coefficients = np.polyfit(volume_array, energy_array, degree)
        first = np.polyder(coefficients, 1)
        second = np.polyder(coefficients, 2)
        roots = np.roots(first)
        real_roots = np.real(roots[np.isclose(np.imag(roots), 0.0)])
        if real_roots.size == 0:
            raise ValueError(
                "the polynomial estimate did not produce a real stationary point"
            )
        curvatures = np.asarray([np.polyval(second, root) for root in real_roots])
        minima = real_roots[curvatures > 0.0]
        candidates = minima if minima.size else real_roots
        v0 = float(
            candidates[
                np.argmin(np.abs(candidates - volume_array[np.argmin(energy_array)]))
            ]
        )
        e0 = float(np.polyval(coefficients, v0))
        k0 = float(v0 * np.polyval(second, v0))
        if not np.isfinite(k0) or k0 <= 0.0:
            k0 = float(abs(k0)) if np.isfinite(k0) and k0 != 0.0 else 1.0
        return EOSParameters(E0=e0, K0=k0, KP=4.0, KPP=0.0, V0=v0)

    def pressure(
        self,
        eos: str | EOSModel,
        parameters: ParameterLike,
        volume: ArrayLike,
        *,
        order: int | None = None,
    ) -> np.ndarray:
        """Evaluate the matching pressure EOS from energy-fit parameters."""
        model = self.model(eos, order)
        resolved = resolve_energy_parameters(model, parameters)
        return PressureEOS().pressure(model, resolved, volume)

    def murnaghan(
        self,
        volume: ArrayLike,
        E0: float,
        K0: float,
        KP: float,
        V0: float,
    ) -> np.ndarray:
        r"""Return the volume-integrated Murnaghan energy.

        For :math:`K'_0\ne1`,

        .. math::

            E(V)=E_0+\frac{K_0V}{K'_0}
            \left[\frac{(V_0/V)^{K'_0}}{K'_0-1}+1\right]
            -\frac{K_0V_0}{K'_0-1}.

        Parameters
        ----------
        volume : array-like
            Positive volume values.
        E0, K0, KP, V0 : float
            Reference energy, bulk modulus, first pressure derivative and
            reference volume.

        Returns
        -------
        ndarray
            Energy values at ``volume``.
        """
        return self.evaluate("M", volume, [E0, K0, KP, V0])

    def birchmurnaghan(
        self,
        volume: ArrayLike,
        E0: float,
        K0: float,
        KP: float,
        V0: float,
    ) -> np.ndarray:
        r"""Return the third-order Birch-Murnaghan energy.

        With :math:`\eta=(V_0/V)^{1/3}`,

        .. math::

            E(V)=E_0+\frac{9K_0V_0}{16}\left[
            K'_0(\eta^2-1)^3+(\eta^2-1)^2(6-4\eta^2)\right].

        Parameters
        ----------
        volume : array-like
            Positive volume values.
        E0, K0, KP, V0 : float
            Reference energy, bulk modulus, first pressure derivative and
            reference volume.

        Returns
        -------
        ndarray
            Energy values at ``volume``.
        """
        return self.evaluate("BM3", volume, [E0, K0, KP, V0])

    def poirier_tarantola(
        self,
        volume: ArrayLike,
        E0: float,
        K0: float,
        KP: float,
        V0: float,
    ) -> np.ndarray:
        r"""Return the third-order natural-strain energy.

        With :math:`\rho=\ln(V_0/V)`,

        .. math::

            E(V)=E_0+\frac{K_0V_0\rho^2}{6}
            \left[3+\rho(K'_0-2)\right].

        Parameters
        ----------
        volume : array-like
            Positive volume values.
        E0, K0, KP, V0 : float
            Reference energy, bulk modulus, first pressure derivative and
            reference volume.

        Returns
        -------
        ndarray
            Energy values at ``volume``.
        """
        return self.evaluate("PT3", volume, [E0, K0, KP, V0])

    pouriertarantola = poirier_tarantola

    def vinet(
        self,
        volume: ArrayLike,
        E0: float,
        K0: float,
        KP: float,
        V0: float,
    ) -> np.ndarray:
        r"""Return the third-order Vinet energy.

        With :math:`x=(V/V_0)^{1/3}`,

        .. math::

            E(V)=E_0+\frac{2K_0V_0}{(K'_0-1)^2}
            \left\{2-\left[5+3K'_0(x-1)-3x\right]
            \exp\left[-\frac{3}{2}(K'_0-1)(x-1)\right]\right\}.

        Parameters
        ----------
        volume : array-like
            Positive volume values.
        E0, K0, KP, V0 : float
            Reference energy, bulk modulus, first pressure derivative and
            reference volume.

        Returns
        -------
        ndarray
            Energy values at ``volume``.
        """
        return self.evaluate("V3", volume, [E0, K0, KP, V0])

    @staticmethod
    def _fit_metadata(model: EOSModel) -> dict[str, object]:
        return {
            "model": model.family.value,
            "eos_tag": model.tag,
            "eos_order": model.order,
            "parameter_order": list(model.energy_parameter_names),
        }

    @staticmethod
    def _murnaghan_energy(volume: np.ndarray, pars: EOSParameters) -> np.ndarray:
        if pars.E0 is None:
            raise ValueError("energy EOS parameters require E0")
        if np.isclose(pars.KP, 1.0):
            ratio = volume / pars.V0
            return pars.E0 + pars.K0 * (volume - pars.V0 - pars.V0 * np.log(ratio))
        return (
            pars.E0
            + pars.K0
            * volume
            / pars.KP
            * (((pars.V0 / volume) ** pars.KP) / (pars.KP - 1.0) + 1.0)
            - pars.V0 * pars.K0 / (pars.KP - 1.0)
        )

    @staticmethod
    def _birch_murnaghan_third_order_energy(
        volume: np.ndarray,
        pars: EOSParameters,
    ) -> np.ndarray:
        """Return the third-order Birch-Murnaghan energy.

        The algebra is equivalent to the general finite-strain expansion.
        This operation order is retained to keep BM3 numerical results exactly
        reproducible after introducing order-selectable EOS models.
        """
        if pars.E0 is None:
            raise ValueError("energy EOS parameters require E0")
        eta = (pars.V0 / volume) ** (1.0 / 3.0)
        eta2_minus_one = eta**2 - 1.0
        return pars.E0 + (9.0 * pars.K0 * pars.V0 / 16.0) * (
            pars.KP * eta2_minus_one**3 + eta2_minus_one**2 * (6.0 - 4.0 * eta**2)
        )

    @staticmethod
    def _strain_energy(
        volume: np.ndarray,
        pars: EOSParameters,
        family: EOSFamily,
    ) -> np.ndarray:
        if pars.E0 is None:
            raise ValueError("energy EOS parameters require E0")
        if family is EOSFamily.BIRCH_MURNAGHAN:
            strain = 0.5 * ((pars.V0 / volume) ** (2.0 / 3.0) - 1.0)
            a = 1.5 * (pars.KP - 4.0)
            b = 1.5 * (
                pars.K0 * pars.KPP + (pars.KP - 4.0) * (pars.KP - 3.0) + 35.0 / 9.0
            )
        elif family is EOSFamily.NATURAL_STRAIN:
            strain = np.log(pars.V0 / volume) / 3.0
            delta = pars.KP - 2.0
            a = 1.5 * delta
            b = 1.5 * (1.0 + pars.K0 * pars.KPP + delta + delta**2)
        else:
            raise ValueError(f"unsupported strain family: {family.value}")
        return pars.E0 + 9.0 * pars.K0 * pars.V0 * (
            strain**2 / 2.0 + a * strain**3 / 3.0 + b * strain**4 / 4.0
        )

    @staticmethod
    def _vinet_energy(
        volume: np.ndarray,
        pars: EOSParameters,
        order: int | None,
    ) -> np.ndarray:
        if pars.E0 is None:
            raise ValueError("energy EOS parameters require E0")
        x = (volume / pars.V0) ** (1.0 / 3.0)
        if order == 2:
            return pars.E0 + 4.5 * pars.K0 * pars.V0 * (x - 1.0) ** 2
        return pars.E0 + 2.0 * pars.K0 * pars.V0 / (pars.KP - 1.0) ** 2 * (
            2.0
            - (5.0 + 3.0 * pars.KP * (x - 1.0) - 3.0 * x)
            * np.exp(-1.5 * (pars.KP - 1.0) * (x - 1.0))
        )

    @staticmethod
    def _validate_volume(volume: ArrayLike) -> np.ndarray:
        values = np.asarray(volume, dtype=np.float64)
        if not np.all(np.isfinite(values)):
            raise ValueError("volume values must be finite")
        if np.any(values <= 0.0):
            raise ValueError("volume values must be positive")
        return values


class EnergyEOSFitModel(BaseFitModel):
    """Adapt an integrated energy EOS to the general fitting contract."""

    def __init__(self, eos: EnergyEOS, model: str | EOSModel) -> None:
        self._eos = eos
        self._model = eos.model(model)

    @property
    def name(self) -> str:
        """Return the compact EOS tag."""
        return self._model.tag

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return the free parameter names for this order."""
        return self._model.energy_parameter_names

    def evaluate(
        self, x: ParameterArrayLike, parameters: ParameterArrayLike
    ) -> np.ndarray:
        """Evaluate the integrated EOS."""
        return self._eos.evaluate(self._model, x, parameters)

    def initial_guess(self, x: ParameterArrayLike, y: ParameterArrayLike) -> np.ndarray:
        """Return an order-aware initial parameter vector."""
        full = self._eos.guess(x, y)
        if self._model.order == 2:
            if self._model.family is EOSFamily.BIRCH_MURNAGHAN:
                full = EOSParameters(
                    E0=full.E0, K0=full.K0, KP=4.0, KPP=0.0, V0=full.V0
                )
            elif self._model.family is EOSFamily.NATURAL_STRAIN:
                full = EOSParameters(
                    E0=full.E0, K0=full.K0, KP=2.0, KPP=0.0, V0=full.V0
                )
            elif self._model.family is EOSFamily.VINET:
                full = EOSParameters(
                    E0=full.E0, K0=full.K0, KP=1.0, KPP=0.0, V0=full.V0
                )
        return free_energy_parameters(self._model, full)

    def metadata(self) -> dict[str, object]:
        """Return model metadata for fit results."""
        return EnergyEOS._fit_metadata(self._model)


__all__ = ["EnergyEOS", "EnergyEOSFitModel"]
