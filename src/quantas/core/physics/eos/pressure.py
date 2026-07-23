# -*- coding: utf-8 -*-

"""Analytical pressure-volume equations of state.

The module evaluates pressure, isothermal bulk modulus and its pressure
derivatives from a common family-and-order specification.  Parameters fixed or
implied by a lower-order truncation are resolved before the equations are
evaluated.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TypeAlias

import numpy as np

from .parameters import EOSParameters, resolve_pressure_parameters
from .spec import EOSFamily, EOSModel, parse_eos_model

ArrayLike: TypeAlias = np.ndarray | float | Sequence[float]
ParameterArrayLike: TypeAlias = np.ndarray | Sequence[float]
ParameterLike: TypeAlias = ParameterArrayLike | Mapping[str, float] | EOSParameters


class PressureEOS:
    r"""Pressure-volume equations of state and analytical derivatives.

    At a selected model and order, the local isothermal bulk modulus and its
    pressure derivatives are

    .. math::

        K_T(V) = -V\frac{\mathrm dP}{\mathrm dV},\qquad
        K'_T(V) = \frac{\mathrm dK_T}{\mathrm dP},\qquad
        K''_T(V) = \frac{\mathrm d^2K_T}{\mathrm dP^2}.
    """

    def model(self, eos: str | EOSModel, order: int | None = None) -> EOSModel:
        """Return the canonical family-and-order specification."""
        return parse_eos_model(eos, order)

    def canonical_name(self, eos: str | EOSModel, order: int | None = None) -> str:
        """Return the canonical EOS family name."""
        return self.model(eos, order).family.value

    def canonical_tag(self, eos: str | EOSModel, order: int | None = None) -> str:
        """Return the compact canonical EOS tag."""
        return self.model(eos, order).tag

    def pressure(
        self,
        eos: str | EOSModel,
        parameters: ParameterLike,
        volume: ArrayLike,
        *,
        order: int | None = None,
    ) -> np.ndarray:
        """Evaluate pressure at one or more volumes."""
        model = self.model(eos, order)
        pars = resolve_pressure_parameters(model, parameters)
        values = self._validate_volume(volume)
        if model.family is EOSFamily.MURNAGHAN:
            return self._murnaghan_pressure(values, pars)
        if model.family is EOSFamily.BIRCH_MURNAGHAN:
            if model.order == 3:
                return self._birch_murnaghan_third_order_pressure(values, pars)
            pressure, _, _, _ = self._birch_derivatives(values, pars)
            return pressure
        if model.family is EOSFamily.NATURAL_STRAIN:
            pressure, _, _, _ = self._natural_derivatives(values, pars)
            return pressure
        if model.family is EOSFamily.VINET:
            pressure, _, _, _ = self._vinet_derivatives(values, pars)
            return pressure
        if model.family is EOSFamily.TAIT:
            return self._tait_pressure(values, pars)
        raise ValueError(f"unknown pressure EOS: {eos!r}")

    def bulk_modulus(
        self,
        eos: str | EOSModel,
        parameters: ParameterLike,
        volume: ArrayLike,
        *,
        order: int | None = None,
    ) -> np.ndarray:
        """Evaluate the isothermal bulk modulus."""
        model = self.model(eos, order)
        pars = resolve_pressure_parameters(model, parameters)
        values = self._validate_volume(volume)
        if model.family is EOSFamily.MURNAGHAN:
            return pars.K0 * (pars.V0 / values) ** pars.KP
        if model.family is EOSFamily.BIRCH_MURNAGHAN:
            if model.order == 3:
                return self._birch_murnaghan_third_order_bulk_modulus(values, pars)
            _, first, _, _ = self._birch_derivatives(values, pars)
            strain = 0.5 * ((pars.V0 / values) ** (2.0 / 3.0) - 1.0)
            return (1.0 + 2.0 * strain) * first / 3.0
        if model.family is EOSFamily.NATURAL_STRAIN:
            _, first, _, _ = self._natural_derivatives(values, pars)
            return first / 3.0
        if model.family is EOSFamily.VINET:
            _, first, _, _ = self._vinet_derivatives(values, pars)
            x = (values / pars.V0) ** (1.0 / 3.0)
            return -x * first / 3.0
        if model.family is EOSFamily.TAIT:
            pressure = self._tait_pressure(values, pars)
            a, b, c = self._tait_coefficients(pars)
            ratio = values / pars.V0
            return pars.K0 * ratio * (1.0 + b * pressure) ** (c + 1.0)
        raise ValueError(f"unknown pressure EOS: {eos!r}")

    def bulk_modulus_derivative(
        self,
        eos: str | EOSModel,
        parameters: ParameterLike,
        volume: ArrayLike,
        *,
        order: int | None = None,
    ) -> np.ndarray:
        """Evaluate the first pressure derivative of the bulk modulus."""
        model = self.model(eos, order)
        pars = resolve_pressure_parameters(model, parameters)
        values = self._validate_volume(volume)
        if model.family is EOSFamily.MURNAGHAN:
            return np.full_like(values, pars.KP, dtype=np.float64)
        if model.family is EOSFamily.BIRCH_MURNAGHAN:
            if model.order == 3:
                return self._birch_murnaghan_third_order_bulk_derivative(values, pars)
            _, first, second, _ = self._birch_derivatives(values, pars)
            strain = 0.5 * ((pars.V0 / values) ** (2.0 / 3.0) - 1.0)
            g = 1.0 + 2.0 * strain
            return 2.0 / 3.0 + g * second / (3.0 * first)
        if model.family is EOSFamily.NATURAL_STRAIN:
            _, first, second, _ = self._natural_derivatives(values, pars)
            return second / (3.0 * first)
        if model.family is EOSFamily.VINET:
            _, first, second, _ = self._vinet_derivatives(values, pars)
            x = (values / pars.V0) ** (1.0 / 3.0)
            return -(1.0 + x * second / first) / 3.0
        if model.family is EOSFamily.TAIT:
            pressure = self._tait_pressure(values, pars)
            a, b, c = self._tait_coefficients(pars)
            return (pars.KP + 1.0) * ((1.0 - a) * (1.0 + b * pressure) ** c + a) - 1.0
        raise ValueError(f"unknown pressure EOS: {eos!r}")

    def bulk_modulus_second_derivative(
        self,
        eos: str | EOSModel,
        parameters: ParameterLike,
        volume: ArrayLike,
        *,
        order: int | None = None,
    ) -> np.ndarray:
        """Evaluate the second pressure derivative of the bulk modulus."""
        model = self.model(eos, order)
        pars = resolve_pressure_parameters(model, parameters)
        values = self._validate_volume(volume)
        if model.family is EOSFamily.MURNAGHAN:
            return np.zeros_like(values, dtype=np.float64)
        if model.family is EOSFamily.BIRCH_MURNAGHAN:
            _, first, second, third = self._birch_derivatives(values, pars)
            strain = 0.5 * ((pars.V0 / values) ** (2.0 / 3.0) - 1.0)
            g = 1.0 + 2.0 * strain
            derivative = (
                (2.0 * second + g * third) / first - g * (second / first) ** 2
            ) / 3.0
            return derivative / first
        if model.family is EOSFamily.NATURAL_STRAIN:
            _, first, second, third = self._natural_derivatives(values, pars)
            return (third * first - second**2) / (3.0 * first**3)
        if model.family is EOSFamily.VINET:
            _, first, second, third = self._vinet_derivatives(values, pars)
            x = (values / pars.V0) ** (1.0 / 3.0)
            ratio = second / first
            derivative = -(ratio + x * (third / first - ratio**2)) / 3.0
            return derivative / first
        if model.family is EOSFamily.TAIT:
            pressure = self._tait_pressure(values, pars)
            a, b, c = self._tait_coefficients(pars)
            return (
                (pars.KP + 1.0) * (1.0 - a) * c * b * (1.0 + b * pressure) ** (c - 1.0)
            )
        raise ValueError(f"unknown pressure EOS: {eos!r}")

    @staticmethod
    def murnaghan(volume: ArrayLike, K0: float, KP: float, V0: float) -> np.ndarray:
        r"""Return pressure from the Murnaghan equation.

        .. math::

            P(V)=\frac{K_0}{K'_0}
            \left[\left(\frac{V_0}{V}\right)^{K'_0}-1\right].
        """
        pars = EOSParameters(K0=K0, KP=KP, KPP=0.0, V0=V0)
        return PressureEOS._murnaghan_pressure(
            np.asarray(volume, dtype=np.float64), pars
        )

    @staticmethod
    def birch_murnaghan(
        volume: ArrayLike,
        K0: float,
        KP: float,
        V0: float,
    ) -> np.ndarray:
        r"""Return pressure from the third-order Birch-Murnaghan EOS.

        With :math:`f_E=[(V_0/V)^{2/3}-1]/2`,

        .. math::

            P(V)=3K_0f_E(1+2f_E)^{5/2}
            \left[1+\frac{3}{2}(K'_0-4)f_E\right].

        Parameters
        ----------
        volume : array-like
            Positive volume values.
        K0, KP, V0 : float
            Reference bulk modulus, first pressure derivative and reference
            volume.

        Returns
        -------
        ndarray
            Pressure values in the same scale as ``K0``.
        """
        model = parse_eos_model("BM3")
        pars = resolve_pressure_parameters(model, [K0, KP, V0])
        return PressureEOS().pressure(model, pars, volume)

    @staticmethod
    def poirier_tarantola(
        volume: ArrayLike,
        K0: float,
        KP: float,
        V0: float,
    ) -> np.ndarray:
        r"""Return pressure from the third-order natural-strain EOS.

        With :math:`f_N=\ln(V_0/V)/3`,

        .. math::

            P(V)=3K_0\frac{V_0}{V}f_N
            \left[1+\frac{3}{2}(K'_0-2)f_N\right].

        Parameters
        ----------
        volume : array-like
            Positive volume values.
        K0, KP, V0 : float
            Reference bulk modulus, first pressure derivative and reference
            volume.

        Returns
        -------
        ndarray
            Pressure values in the same scale as ``K0``.
        """
        model = parse_eos_model("PT3")
        pars = resolve_pressure_parameters(model, [K0, KP, V0])
        return PressureEOS().pressure(model, pars, volume)

    @staticmethod
    def vinet(
        volume: ArrayLike,
        K0: float,
        KP: float,
        V0: float,
    ) -> np.ndarray:
        r"""Return pressure from the third-order Vinet EOS.

        With :math:`x=(V/V_0)^{1/3}`,

        .. math::

            P(V)=3K_0\frac{1-x}{x^2}
            \exp\left[\frac{3}{2}(K'_0-1)(1-x)\right].

        Parameters
        ----------
        volume : array-like
            Positive volume values.
        K0, KP, V0 : float
            Reference bulk modulus, first pressure derivative and reference
            volume.

        Returns
        -------
        ndarray
            Pressure values in the same scale as ``K0``.
        """
        model = parse_eos_model("V3")
        pars = resolve_pressure_parameters(model, [K0, KP, V0])
        return PressureEOS().pressure(model, pars, volume)

    @staticmethod
    def _birch_murnaghan_third_order_pressure(
        volume: np.ndarray,
        pars: EOSParameters,
    ) -> np.ndarray:
        """Return BM3 pressure with the reproducible third-order expression."""
        strain = 0.5 * ((pars.V0 / volume) ** (2.0 / 3.0) - 1.0)
        return (
            3.0
            * pars.K0
            * strain
            * (1.0 + 2.0 * strain) ** 2.5
            * (1.0 + 1.5 * (pars.KP - 4.0) * strain)
        )

    @staticmethod
    def _birch_murnaghan_third_order_bulk_modulus(
        volume: np.ndarray,
        pars: EOSParameters,
    ) -> np.ndarray:
        """Return the BM3 bulk modulus with the reproducible expression."""
        x = (pars.V0 / volume) ** (1.0 / 3.0)
        y = x**2 - 1.0
        a = 0.75 * (pars.KP - 4.0)
        h = 5.0 * y * (1.0 + a * y) + 2.0 * x**2 * (1.0 + 2.0 * a * y)
        dpdx = 1.5 * pars.K0 * x**4 * h
        return x * dpdx / 3.0

    @staticmethod
    def _birch_murnaghan_third_order_bulk_derivative(
        volume: np.ndarray,
        pars: EOSParameters,
    ) -> np.ndarray:
        """Return BM3 ``dK/dP`` with the reproducible expression."""
        x = (pars.V0 / volume) ** (1.0 / 3.0)
        y = x**2 - 1.0
        a = 0.75 * (pars.KP - 4.0)
        h = 5.0 * y * (1.0 + a * y) + 2.0 * x**2 * (1.0 + 2.0 * a * y)
        dhdx = 14.0 * x * (1.0 + 2.0 * a * y) + 8.0 * a * x**3
        dpdx = 1.5 * pars.K0 * x**4 * h
        d2pdx2 = 1.5 * pars.K0 * (4.0 * x**3 * h + x**4 * dhdx)
        return (1.0 + x * d2pdx2 / dpdx) / 3.0

    @staticmethod
    def _murnaghan_pressure(volume: np.ndarray, pars: EOSParameters) -> np.ndarray:
        return (pars.K0 / pars.KP) * ((pars.V0 / volume) ** pars.KP - 1.0)

    @staticmethod
    def _strain_polynomial(
        pars: EOSParameters, family: EOSFamily
    ) -> tuple[float, float]:
        if family is EOSFamily.BIRCH_MURNAGHAN:
            a = 1.5 * (pars.KP - 4.0)
            b = 1.5 * (
                pars.K0 * pars.KPP + (pars.KP - 4.0) * (pars.KP - 3.0) + 35.0 / 9.0
            )
            return a, b
        if family is EOSFamily.NATURAL_STRAIN:
            delta = pars.KP - 2.0
            a = 1.5 * delta
            b = 1.5 * (1.0 + pars.K0 * pars.KPP + delta + delta**2)
            return a, b
        raise ValueError(f"unsupported strain family: {family.value}")

    @classmethod
    def _birch_derivatives(
        cls,
        volume: np.ndarray,
        pars: EOSParameters,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        f = 0.5 * ((pars.V0 / volume) ** (2.0 / 3.0) - 1.0)
        g = 1.0 + 2.0 * f
        a, b = cls._strain_polynomial(pars, EOSFamily.BIRCH_MURNAGHAN)
        h = f + a * f**2 + b * f**3
        h1 = 1.0 + 2.0 * a * f + 3.0 * b * f**2
        h2 = 2.0 * a + 6.0 * b * f
        h3 = np.full_like(f, 6.0 * b)
        prefactor = 3.0 * pars.K0
        pressure = prefactor * g**2.5 * h
        first = prefactor * g**1.5 * (5.0 * h + g * h1)
        second = prefactor * g**0.5 * (15.0 * h + 10.0 * g * h1 + g**2 * h2)
        third = (
            prefactor
            * g**-0.5
            * (15.0 * h + 45.0 * g * h1 + 15.0 * g**2 * h2 + g**3 * h3)
        )
        return pressure, first, second, third

    @classmethod
    def _natural_derivatives(
        cls,
        volume: np.ndarray,
        pars: EOSParameters,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        f = np.log(pars.V0 / volume) / 3.0
        a, b = cls._strain_polynomial(pars, EOSFamily.NATURAL_STRAIN)
        h = f + a * f**2 + b * f**3
        h1 = 1.0 + 2.0 * a * f + 3.0 * b * f**2
        h2 = 2.0 * a + 6.0 * b * f
        h3 = np.full_like(f, 6.0 * b)
        prefactor = 3.0 * pars.K0 * np.exp(3.0 * f)
        pressure = prefactor * h
        first = prefactor * (3.0 * h + h1)
        second = prefactor * (9.0 * h + 6.0 * h1 + h2)
        third = prefactor * (27.0 * h + 27.0 * h1 + 9.0 * h2 + h3)
        return pressure, first, second, third

    @staticmethod
    def _vinet_derivatives(
        volume: np.ndarray,
        pars: EOSParameters,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        x = (volume / pars.V0) ** (1.0 / 3.0)
        a = 1.5 * (pars.KP - 1.0)
        exponential = np.exp(a * (1.0 - x))
        u = (1.0 - x) / x**2
        u1 = (x - 2.0) / x**3
        u2 = (6.0 - 2.0 * x) / x**4
        u3 = (6.0 * x - 24.0) / x**5
        prefactor = 3.0 * pars.K0 * exponential
        pressure = prefactor * u
        first = prefactor * (u1 - a * u)
        second = prefactor * (u2 - 2.0 * a * u1 + a**2 * u)
        third = prefactor * (u3 - 3.0 * a * u2 + 3.0 * a**2 * u1 - a**3 * u)
        return pressure, first, second, third

    @staticmethod
    def _tait_coefficients(pars: EOSParameters) -> tuple[float, float, float]:
        denominator = 1.0 + pars.KP + pars.K0 * pars.KPP
        a = (1.0 + pars.KP) / denominator
        b = pars.KP / pars.K0 - pars.KPP / (1.0 + pars.KP)
        c = denominator / (pars.KP**2 + pars.KP - pars.K0 * pars.KPP)
        if not np.all(np.isfinite([a, b, c])) or a == 0.0 or b == 0.0 or c == 0.0:
            raise ValueError("Tait parameters produce a singular equation")
        return float(a), float(b), float(c)

    @classmethod
    def _tait_pressure(cls, volume: np.ndarray, pars: EOSParameters) -> np.ndarray:
        a, b, c = cls._tait_coefficients(pars)
        argument = ((volume / pars.V0) + a - 1.0) / a
        if np.any(argument <= 0.0):
            raise ValueError("Tait EOS is undefined for the requested volume")
        return (argument ** (-1.0 / c) - 1.0) / b

    @staticmethod
    def _validate_volume(volume: ArrayLike) -> np.ndarray:
        values = np.asarray(volume, dtype=np.float64)
        if not np.all(np.isfinite(values)):
            raise ValueError("volume values must be finite")
        if np.any(values <= 0.0):
            raise ValueError("volume values must be positive")
        return values


__all__ = ["PressureEOS"]
