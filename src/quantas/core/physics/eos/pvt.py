# -*- coding: utf-8 -*-

r"""Pressure-volume-temperature equation-of-state coupling.

The classes in this module compose an isothermal pressure-volume EOS with a
thermal model without duplicating either family of equations.  Three coupling
prescriptions are implemented:

``linear-bulk-modulus``
    The zero-pressure bulk modulus varies linearly with temperature,

    .. math::

        K_{0T}=K_{00}+\left(\frac{\partial K_0}{\partial T}\right)
        (T-T_{\mathrm{ref}}).

``anderson-gruneisen``
    The zero-pressure bulk modulus follows

    .. math::

        K_{0T}=K_{00}\left(\frac{V_{00}}{V_{0T}}\right)^\delta.

``thermal-pressure``
    The total pressure is the sum of a reference isotherm and a selected
    Holland--Powell Einstein or Mie--Grüneisen--Debye thermal-pressure term,

    .. math::

        P(V,T)=P(V,T_{\mathrm{ref}})+P_{\mathrm{th}}(T).

The first two prescriptions use any audited :class:`TemperatureEOS` to obtain
:math:`V_{0T}`.  The thermal-pressure formulation is a separate coupling model
and does not consume an independent V--T equation.  In all prescriptions the
pressure derivative parameters :math:`K'_0` and :math:`K''_0` are held
temperature-independent unless a future, explicitly named coupling extends
that assumption.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
import math
from typing import TypeAlias

import numpy as np
from scipy.optimize import brentq

from .parameters import EOSParameters, resolve_pressure_parameters
from .pressure import PressureEOS
from .spec import EOSModel, parse_eos_model
from .thermal_pressure import (
    MGDNormalization,
    MGDParameters,
    MGDThermalPressure,
    ThermalPressureFamily,
    ThermalPressureModel,
    parse_thermal_pressure_model,
)
from .temperature import (
    TemperatureEOS,
    TemperatureEOSModel,
    TemperatureEOSParameters,
    parse_temperature_eos_model,
)

ArrayLike: TypeAlias = np.ndarray | float | Sequence[float]
PressureParameterLike: TypeAlias = (
    EOSParameters | Mapping[str, float] | np.ndarray | Sequence[float]
)
TemperatureParameterLike: TypeAlias = TemperatureEOSParameters | Mapping[str, float]
_TINY = float(np.finfo(np.float64).tiny)


class PVTCouplingFamily(str, Enum):
    """Supported pressure-volume-temperature coupling prescriptions."""

    LINEAR_BULK_MODULUS = "linear-bulk-modulus"
    ANDERSON_GRUNEISEN = "anderson-gruneisen"
    THERMAL_PRESSURE = "thermal-pressure"


_COUPLING_ALIASES: dict[str, PVTCouplingFamily] = {
    "linear": PVTCouplingFamily.LINEAR_BULK_MODULUS,
    "linearbulkmodulus": PVTCouplingFamily.LINEAR_BULK_MODULUS,
    "lineark": PVTCouplingFamily.LINEAR_BULK_MODULUS,
    "dkdt": PVTCouplingFamily.LINEAR_BULK_MODULUS,
    "andersongruneisen": PVTCouplingFamily.ANDERSON_GRUNEISEN,
    "anderson-gruneisen": PVTCouplingFamily.ANDERSON_GRUNEISEN,
    "ag": PVTCouplingFamily.ANDERSON_GRUNEISEN,
    "delta": PVTCouplingFamily.ANDERSON_GRUNEISEN,
    "thermalpressure": PVTCouplingFamily.THERMAL_PRESSURE,
    "thermal-pressure": PVTCouplingFamily.THERMAL_PRESSURE,
    "pth": PVTCouplingFamily.THERMAL_PRESSURE,
}


@dataclass(frozen=True, slots=True)
class PVTModel:
    """Describe one compositional pressure-volume-temperature model.

    Parameters
    ----------
    pressure_model : EOSModel or str
        Reference isothermal pressure-volume EOS.
    coupling : PVTCouplingFamily or str
        Temperature-coupling prescription.
    temperature_model : TemperatureEOSModel, str, or None, optional
        V--T model used to calculate ``V0(T)`` for the linear-bulk-modulus and
        Anderson--Gruneisen prescriptions. It must be omitted for the
        thermal-pressure prescription.
    thermal_pressure_model : ThermalPressureModel, str, or None, optional
        Oscillator model used by the thermal-pressure prescription. Omission
        preserves the historical Holland--Powell Einstein default.
    mgd_normalization : MGDNormalization or None, optional
        Cell or molar atom-count normalization. It is required only for MGD.

    Raises
    ------
    ValueError
        If the temperature-model requirement is inconsistent with the
        coupling family.
    """

    pressure_model: EOSModel | str
    coupling: PVTCouplingFamily | str
    temperature_model: TemperatureEOSModel | str | None = None
    thermal_pressure_model: ThermalPressureModel | str | None = None
    mgd_normalization: MGDNormalization | None = None

    def __post_init__(self) -> None:
        """Normalize component models and validate the composition."""
        pressure = parse_eos_model(self.pressure_model)
        coupling = parse_pvt_coupling(self.coupling)
        thermal = (
            None
            if self.temperature_model is None
            else parse_temperature_eos_model(self.temperature_model)
        )
        thermal_pressure = (
            None
            if self.thermal_pressure_model is None
            else parse_thermal_pressure_model(self.thermal_pressure_model)
        )
        normalization = self.mgd_normalization
        if coupling is PVTCouplingFamily.THERMAL_PRESSURE:
            if thermal is not None:
                raise ValueError(
                    "thermal-pressure coupling defines its own thermal model; "
                    "temperature_model must be omitted"
                )
            if thermal_pressure is None:
                thermal_pressure = ThermalPressureModel(
                    ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN
                )
            if thermal_pressure.requires_mgd_normalization:
                if normalization is None:
                    raise ValueError(
                        "Mie-Gruneisen-Debye thermal pressure requires "
                        "mgd_normalization"
                    )
            elif normalization is not None:
                raise ValueError(
                    "mgd_normalization is only valid for Mie-Gruneisen-Debye "
                    "thermal pressure"
                )
        else:
            if thermal is None:
                raise ValueError(
                    f"{coupling.value} coupling requires a temperature_model"
                )
            if thermal_pressure is not None or normalization is not None:
                raise ValueError(
                    "thermal_pressure_model and mgd_normalization are only valid "
                    "for thermal-pressure coupling"
                )
        object.__setattr__(self, "pressure_model", pressure)
        object.__setattr__(self, "coupling", coupling)
        object.__setattr__(self, "temperature_model", thermal)
        object.__setattr__(self, "thermal_pressure_model", thermal_pressure)
        object.__setattr__(self, "mgd_normalization", normalization)

    @property
    def pressure_spec(self) -> EOSModel:
        """Return the canonical reference isothermal EOS specification."""
        value = self.pressure_model
        assert isinstance(value, EOSModel)
        return value

    @property
    def coupling_family(self) -> PVTCouplingFamily:
        """Return the canonical temperature-coupling family."""
        value = self.coupling
        assert isinstance(value, PVTCouplingFamily)
        return value

    @property
    def thermal_spec(self) -> TemperatureEOSModel | None:
        """Return the canonical V--T model, when the coupling uses one."""
        value = self.temperature_model
        assert value is None or isinstance(value, TemperatureEOSModel)
        return value

    @property
    def thermal_pressure_spec(self) -> ThermalPressureModel | None:
        """Return the selected thermal-pressure oscillator model."""
        value = self.thermal_pressure_model
        assert value is None or isinstance(value, ThermalPressureModel)
        return value

    @property
    def mgd_normalization_spec(self) -> MGDNormalization | None:
        """Return the MGD normalization, when applicable."""
        value = self.mgd_normalization
        assert value is None or isinstance(value, MGDNormalization)
        return value

    @property
    def tag(self) -> str:
        """Return a stable composite model tag."""
        thermal = self.thermal_spec
        thermal_tag = "none" if thermal is None else thermal.tag
        thermal_pressure = self.thermal_pressure_spec
        if (
            thermal_pressure is not None
            and thermal_pressure.family_name
            is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE
        ):
            thermal_tag = thermal_pressure.tag
        return f"{self.pressure_spec.tag}+{thermal_tag}+{self.coupling_family.value}"

    def as_dict(self) -> dict[str, object]:
        """Return a serialization-ready model description."""
        thermal = self.thermal_spec
        description: dict[str, object] = {
            "pressure_model": self.pressure_spec.as_dict(),
            "temperature_model": (
                None
                if thermal is None
                else {
                    "family": thermal.family.value,
                    "variant": (thermal.variant.value if thermal.variant else None),
                    "tag": thermal.tag,
                }
            ),
            "coupling": self.coupling_family.value,
            "tag": self.tag,
        }
        thermal_pressure = self.thermal_pressure_spec
        if (
            thermal_pressure is not None
            and thermal_pressure.family_name
            is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE
        ):
            normalization = self.mgd_normalization_spec
            assert normalization is not None
            description["thermal_pressure_model"] = thermal_pressure.as_dict()
            description["mgd_normalization"] = normalization.as_dict()
        return description


def parse_pvt_coupling(
    coupling: str | PVTCouplingFamily,
) -> PVTCouplingFamily:
    """Return a canonical P--V--T coupling family.

    Parameters
    ----------
    coupling : str or PVTCouplingFamily
        Canonical name or documented alias.

    Returns
    -------
    PVTCouplingFamily
        Canonical coupling family.

    Raises
    ------
    ValueError
        If the coupling name is unknown.
    """
    if isinstance(coupling, PVTCouplingFamily):
        return coupling
    text = str(coupling).strip().lower()
    key = "".join(character for character in text if character.isalnum())
    resolved = _COUPLING_ALIASES.get(text) or _COUPLING_ALIASES.get(key)
    if resolved is None:
        raise ValueError(f"unknown P-V-T coupling: {coupling!r}")
    return resolved


def available_pvt_couplings() -> tuple[PVTCouplingFamily, ...]:
    """Return all implemented P--V--T coupling prescriptions."""
    return tuple(PVTCouplingFamily)


class PVTEOS:
    r"""Evaluate compositional pressure-volume-temperature equations.

    Notes
    -----
    For the ``linear-bulk-modulus`` and ``anderson-gruneisen`` couplings,
    ``temperature_parameters['V0']`` must equal the reference ``V0``
    of the pressure EOS, and both component models must share the same
    reference temperature.  This constraint prevents two inconsistent
    reference volumes from entering a nominally single P--V--T model.
    """

    def __init__(self) -> None:
        self._pressure = PressureEOS()
        self._temperature = TemperatureEOS()
        self._mgd = MGDThermalPressure()

    def pressure(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        temperature_parameters: TemperatureParameterLike | None,
        coupling_parameters: Mapping[str, float],
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Evaluate pressure at paired volumes and temperatures.

        Parameters
        ----------
        model : PVTModel
            Composite reference EOS, thermal model, and coupling prescription.
        pressure_parameters : mapping, sequence, or EOSParameters
            Reference ``K0, KP, KPP, V0`` parameters.
        temperature_parameters : mapping, TemperatureEOSParameters, or None
            Parameters of the selected V--T model.  Required for non-thermal-
            pressure couplings and omitted for thermal pressure.
        coupling_parameters : mapping
            ``dK0_dT`` for linear coupling, ``delta`` for Anderson--Grüneisen,
            or ``temperature_ref, alpha_ref, theta_e`` for thermal pressure.
        volume, temperature : array-like
            Broadcast-compatible volumes and temperatures.

        Returns
        -------
        ndarray
            Pressure values in the unit used by ``K0``.
        """
        volume_array, temperature_array = self._broadcast_state(volume, temperature)
        pressure_spec = model.pressure_spec
        coupling = model.coupling_family
        assert isinstance(pressure_spec, EOSModel)
        assert isinstance(coupling, PVTCouplingFamily)
        reference = resolve_pressure_parameters(pressure_spec, pressure_parameters)
        if coupling is PVTCouplingFamily.THERMAL_PRESSURE:
            p_ref = self._pressure.pressure(pressure_spec, reference, volume_array)
            return p_ref + self._thermal_pressure_for_model(
                model,
                reference,
                coupling_parameters,
                volume_array,
                temperature_array,
            )
        thermal_spec = model.thermal_spec
        assert isinstance(thermal_spec, TemperatureEOSModel)
        thermal = self._resolve_temperature_parameters(
            thermal_spec, temperature_parameters
        )
        self._validate_reference_consistency(reference, thermal)
        v0_t = self._temperature.value(thermal_spec, thermal, temperature_array)
        k0_t = self.zero_pressure_bulk_modulus(
            model,
            reference,
            thermal,
            coupling_parameters,
            temperature_array,
        )
        parameters: dict[str, np.ndarray | float] = {
            "K0": k0_t,
            "KP": reference.KP,
            "KPP": reference.KPP,
            "V0": v0_t,
        }
        return self._pressure_vectorized(pressure_spec, parameters, volume_array)

    def reference_volume(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        temperature_parameters: TemperatureParameterLike | None,
        coupling_parameters: Mapping[str, float],
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return the zero-pressure volume at temperature ``T``."""
        temp = self._validate_temperature(temperature)
        pressure_spec = model.pressure_spec
        coupling = model.coupling_family
        assert isinstance(pressure_spec, EOSModel)
        assert isinstance(coupling, PVTCouplingFamily)
        reference = resolve_pressure_parameters(pressure_spec, pressure_parameters)
        if coupling is not PVTCouplingFamily.THERMAL_PRESSURE:
            thermal_spec = model.thermal_spec
            assert isinstance(thermal_spec, TemperatureEOSModel)
            thermal = self._resolve_temperature_parameters(
                thermal_spec, temperature_parameters
            )
            self._validate_reference_consistency(reference, thermal)
            return self._temperature.value(thermal_spec, thermal, temp)
        return np.asarray(
            [
                self._solve_volume_scalar(
                    model,
                    reference,
                    None,
                    coupling_parameters,
                    0.0,
                    float(item),
                )
                for item in temp.reshape(-1)
            ],
            dtype=np.float64,
        ).reshape(temp.shape)

    def zero_pressure_bulk_modulus(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        temperature_parameters: TemperatureParameterLike | None,
        coupling_parameters: Mapping[str, float],
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return :math:`K_{0T}` at one or more temperatures."""
        temp = self._validate_temperature(temperature)
        pressure_spec = model.pressure_spec
        coupling = model.coupling_family
        assert isinstance(pressure_spec, EOSModel)
        assert isinstance(coupling, PVTCouplingFamily)
        reference = resolve_pressure_parameters(pressure_spec, pressure_parameters)
        if coupling is PVTCouplingFamily.THERMAL_PRESSURE:
            v0_t = self.reference_volume(
                model,
                reference,
                None,
                coupling_parameters,
                temp,
            )
            return self.bulk_modulus(
                model,
                reference,
                None,
                coupling_parameters,
                v0_t,
                temp,
            )
        thermal_spec = model.thermal_spec
        assert isinstance(thermal_spec, TemperatureEOSModel)
        thermal = self._resolve_temperature_parameters(
            thermal_spec, temperature_parameters
        )
        self._validate_reference_consistency(reference, thermal)
        if coupling is PVTCouplingFamily.LINEAR_BULK_MODULUS:
            derivative = self._required_finite(coupling_parameters, "dK0_dT")
            values = reference.K0 + derivative * (temp - thermal.temperature_ref)
        else:
            delta = self._required_finite(coupling_parameters, "delta")
            v0_t = self._temperature.value(thermal_spec, thermal, temp)
            values = reference.K0 * (reference.V0 / v0_t) ** delta
        invalid = np.logical_or(values <= 0.0, ~np.isfinite(values))
        if np.any(invalid):
            index = int(np.flatnonzero(invalid)[0])
            value = float(np.ravel(values)[index])
            temperature_value = float(np.ravel(temp)[index])
            if coupling is PVTCouplingFamily.LINEAR_BULK_MODULUS:
                derivative = self._required_finite(coupling_parameters, "dK0_dT")
                raise ValueError(
                    "P-V-T coupling is invalid at "
                    f"T={temperature_value:.12g} K: "
                    f"K0(T)={value:.12g} GPa from K0={reference.K0:.12g} GPa, "
                    f"dK0_dT={derivative:.12g} GPa K^-1, and "
                    f"Tref={thermal.temperature_ref:.12g} K"
                )
            raise ValueError(
                "P-V-T coupling is invalid at "
                f"T={temperature_value:.12g} K: K0(T)={value:.12g} GPa"
            )
        return np.asarray(values, dtype=np.float64)

    def zero_pressure_bulk_modulus_temperature_derivative(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        temperature_parameters: TemperatureParameterLike | None,
        coupling_parameters: Mapping[str, float],
        temperature: ArrayLike,
    ) -> np.ndarray:
        r"""Return :math:`\partial K_{0T}/\partial T` at zero pressure."""
        temp = self._validate_temperature(temperature)
        coupling = model.coupling_family
        assert isinstance(coupling, PVTCouplingFamily)
        pressure_spec = model.pressure_spec
        assert isinstance(pressure_spec, EOSModel)
        reference = resolve_pressure_parameters(pressure_spec, pressure_parameters)
        if coupling is PVTCouplingFamily.LINEAR_BULK_MODULUS:
            value = self._required_finite(coupling_parameters, "dK0_dT")
            return np.full(temp.shape, value, dtype=np.float64)
        if coupling is PVTCouplingFamily.ANDERSON_GRUNEISEN:
            thermal_spec = model.thermal_spec
            assert isinstance(thermal_spec, TemperatureEOSModel)
            thermal = self._resolve_temperature_parameters(
                thermal_spec, temperature_parameters
            )
            self._validate_reference_consistency(reference, thermal)
            delta = self._required_finite(coupling_parameters, "delta")
            k0_t = self.zero_pressure_bulk_modulus(
                model, reference, thermal, coupling_parameters, temp
            )
            alpha = self._temperature.expansion_coefficient(thermal_spec, thermal, temp)
            return -delta * alpha * k0_t
        # Thermal pressure: evaluate dK0/dT from the analytical chain rule
        # dK/dT = (dK/dP)*(dP_ref/dT)|V0 + (dK/dV)*(dV0/dT).  A stable
        # centered derivative of the already analytical K0(T) is preferable to
        # duplicating family-specific third derivatives here.
        step = np.cbrt(np.finfo(np.float64).eps) * np.maximum(temp, 1.0)
        lower = np.maximum(temp - step, 0.0)
        upper = temp + step
        k_upper = self.zero_pressure_bulk_modulus(
            model, reference, None, coupling_parameters, upper
        )
        k_lower = self.zero_pressure_bulk_modulus(
            model, reference, None, coupling_parameters, lower
        )
        return (k_upper - k_lower) / (upper - lower)

    def expansion_coefficient_zero_pressure(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        temperature_parameters: TemperatureParameterLike | None,
        coupling_parameters: Mapping[str, float],
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return the zero-pressure volume expansion coefficient."""
        temp = self._validate_temperature(temperature)
        coupling = model.coupling_family
        assert isinstance(coupling, PVTCouplingFamily)
        if coupling is not PVTCouplingFamily.THERMAL_PRESSURE:
            thermal_spec = model.thermal_spec
            assert isinstance(thermal_spec, TemperatureEOSModel)
            thermal = self._resolve_temperature_parameters(
                thermal_spec, temperature_parameters
            )
            return self._temperature.expansion_coefficient(thermal_spec, thermal, temp)
        pressure_spec = model.pressure_spec
        assert isinstance(pressure_spec, EOSModel)
        reference = resolve_pressure_parameters(pressure_spec, pressure_parameters)
        k0_t = self.zero_pressure_bulk_modulus(
            model, reference, None, coupling_parameters, temp
        )
        v0_t = self.reference_volume(model, reference, None, coupling_parameters, temp)
        return (
            self._thermal_pressure_temperature_derivative_for_model(
                model,
                reference,
                coupling_parameters,
                v0_t,
                temp,
            )
            / k0_t
        )

    def bulk_modulus(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        temperature_parameters: TemperatureParameterLike | None,
        coupling_parameters: Mapping[str, float],
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return the isothermal bulk modulus at ``V`` and ``T``."""
        volume_array, temperature_array = self._broadcast_state(volume, temperature)
        pressure_spec = model.pressure_spec
        coupling = model.coupling_family
        assert isinstance(pressure_spec, EOSModel)
        assert isinstance(coupling, PVTCouplingFamily)
        reference = resolve_pressure_parameters(pressure_spec, pressure_parameters)
        if coupling is PVTCouplingFamily.THERMAL_PRESSURE:
            static = self._pressure.bulk_modulus(pressure_spec, reference, volume_array)
            thermal_pressure = model.thermal_pressure_spec
            assert isinstance(thermal_pressure, ThermalPressureModel)
            if (
                thermal_pressure.family_name
                is ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN
            ):
                return static
            normalization = model.mgd_normalization_spec
            assert normalization is not None
            mgd_parameters = self._resolve_mgd_parameters(
                thermal_pressure, coupling_parameters
            )
            return static + self._mgd.bulk_modulus_contribution(
                thermal_pressure,
                mgd_parameters,
                normalization,
                reference.V0,
                volume_array,
                temperature_array,
            )
        thermal_spec = model.thermal_spec
        assert isinstance(thermal_spec, TemperatureEOSModel)
        thermal = self._resolve_temperature_parameters(
            thermal_spec, temperature_parameters
        )
        self._validate_reference_consistency(reference, thermal)
        v0_t = self._temperature.value(thermal_spec, thermal, temperature_array)
        k0_t = self.zero_pressure_bulk_modulus(
            model, reference, thermal, coupling_parameters, temperature_array
        )
        return self._bulk_modulus_vectorized(
            pressure_spec,
            {"K0": k0_t, "KP": reference.KP, "KPP": reference.KPP, "V0": v0_t},
            volume_array,
        )

    def volume(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        temperature_parameters: TemperatureParameterLike | None,
        coupling_parameters: Mapping[str, float],
        pressure: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Invert the P--V--T model and return volume at paired ``P, T``."""
        pressure_array, temperature_array = self._broadcast_state(pressure, temperature)
        pressure_spec = model.pressure_spec
        assert isinstance(pressure_spec, EOSModel)
        reference = resolve_pressure_parameters(pressure_spec, pressure_parameters)
        thermal = None
        if model.coupling_family is not PVTCouplingFamily.THERMAL_PRESSURE:
            thermal_spec = model.thermal_spec
            assert isinstance(thermal_spec, TemperatureEOSModel)
            thermal = self._resolve_temperature_parameters(
                thermal_spec, temperature_parameters
            )
            self._validate_reference_consistency(reference, thermal)
        values = [
            self._solve_volume_scalar(
                model,
                reference,
                thermal,
                coupling_parameters,
                float(p),
                float(t),
            )
            for p, t in zip(
                pressure_array.reshape(-1),
                temperature_array.reshape(-1),
                strict=True,
            )
        ]
        return np.asarray(values, dtype=np.float64).reshape(pressure_array.shape)

    def _thermal_pressure_for_model(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        coupling_parameters: Mapping[str, float],
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return the selected thermal-pressure contribution."""
        thermal_pressure = model.thermal_pressure_spec
        assert isinstance(thermal_pressure, ThermalPressureModel)
        if (
            thermal_pressure.family_name
            is ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN
        ):
            return self.thermal_pressure(
                pressure_parameters, coupling_parameters, temperature
            )
        normalization = model.mgd_normalization_spec
        assert normalization is not None
        reference = resolve_pressure_parameters(
            model.pressure_spec, pressure_parameters
        )
        parameters = self._resolve_mgd_parameters(thermal_pressure, coupling_parameters)
        return self._mgd.pressure(
            thermal_pressure,
            parameters,
            normalization,
            reference.V0,
            volume,
            temperature,
        )

    def _thermal_pressure_temperature_derivative_for_model(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        coupling_parameters: Mapping[str, float],
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return the selected thermal-pressure temperature derivative."""
        thermal_pressure = model.thermal_pressure_spec
        assert isinstance(thermal_pressure, ThermalPressureModel)
        if (
            thermal_pressure.family_name
            is ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN
        ):
            return self.thermal_pressure_temperature_derivative(
                pressure_parameters, coupling_parameters, temperature
            )
        normalization = model.mgd_normalization_spec
        assert normalization is not None
        reference = resolve_pressure_parameters(
            model.pressure_spec, pressure_parameters
        )
        parameters = self._resolve_mgd_parameters(thermal_pressure, coupling_parameters)
        return self._mgd.temperature_derivative(
            thermal_pressure,
            parameters,
            normalization,
            reference.V0,
            volume,
            temperature,
        )

    def thermal_pressure_contribution(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        coupling_parameters: Mapping[str, float],
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return the selected thermal-pressure contribution.

        Parameters
        ----------
        model : PVTModel
            Composite model using the thermal-pressure coupling.
        pressure_parameters : mapping, sequence, or EOSParameters
            Reference pressure-EOS parameters.
        coupling_parameters : mapping
            Parameters of the selected Einstein or MGD thermal-pressure model.
        volume, temperature : array-like
            Broadcast-compatible state coordinates.

        Returns
        -------
        ndarray
            Thermal-pressure contribution in the pressure unit of ``K0``.

        Raises
        ------
        ValueError
            If ``model`` does not use thermal pressure or its parameters are
            incomplete.
        """
        if model.coupling_family is not PVTCouplingFamily.THERMAL_PRESSURE:
            raise ValueError(
                "thermal_pressure_contribution requires thermal-pressure coupling"
            )
        volume_array, temperature_array = self._broadcast_state(volume, temperature)
        return self._thermal_pressure_for_model(
            model,
            pressure_parameters,
            coupling_parameters,
            volume_array,
            temperature_array,
        )

    def thermal_pressure_temperature_derivative_for_model(
        self,
        model: PVTModel,
        pressure_parameters: PressureParameterLike,
        coupling_parameters: Mapping[str, float],
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        r"""Return :math:`(\partial P_{th}/\partial T)_V`.

        Parameters
        ----------
        model : PVTModel
            Composite model using the thermal-pressure coupling.
        pressure_parameters : mapping, sequence, or EOSParameters
            Reference pressure-EOS parameters.
        coupling_parameters : mapping
            Parameters of the selected Einstein or MGD thermal-pressure model.
        volume, temperature : array-like
            Broadcast-compatible state coordinates.

        Returns
        -------
        ndarray
            Temperature derivative in pressure per kelvin.
        """
        if model.coupling_family is not PVTCouplingFamily.THERMAL_PRESSURE:
            raise ValueError(
                "thermal_pressure_temperature_derivative_for_model requires "
                "thermal-pressure coupling"
            )
        volume_array, temperature_array = self._broadcast_state(volume, temperature)
        return self._thermal_pressure_temperature_derivative_for_model(
            model,
            pressure_parameters,
            coupling_parameters,
            volume_array,
            temperature_array,
        )

    def thermal_pressure(
        self,
        pressure_parameters: PressureParameterLike,
        coupling_parameters: Mapping[str, float],
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return the Holland--Powell Einstein thermal pressure."""
        reference_k0 = self._reference_bulk_modulus(pressure_parameters)
        temp = self._validate_temperature(temperature)
        temperature_ref = self._required_positive(
            coupling_parameters, "temperature_ref"
        )
        alpha_ref = self._required_finite(coupling_parameters, "alpha_ref")
        theta_e = self._required_positive(coupling_parameters, "theta_e")
        xi_ref = self._einstein_xi(theta_e / temperature_ref)
        occupation = self._einstein_occupation(theta_e, temp)
        occupation_ref = float(
            self._einstein_occupation(
                theta_e, np.asarray([temperature_ref], dtype=np.float64)
            )[0]
        )
        return (
            alpha_ref * reference_k0 * theta_e / xi_ref * (occupation - occupation_ref)
        )

    def thermal_pressure_temperature_derivative(
        self,
        pressure_parameters: PressureParameterLike,
        coupling_parameters: Mapping[str, float],
        temperature: ArrayLike,
    ) -> np.ndarray:
        """Return the analytical derivative ``dP_th/dT``."""
        reference_k0 = self._reference_bulk_modulus(pressure_parameters)
        temp = self._validate_temperature(temperature)
        temperature_ref = self._required_positive(
            coupling_parameters, "temperature_ref"
        )
        alpha_ref = self._required_finite(coupling_parameters, "alpha_ref")
        theta_e = self._required_positive(coupling_parameters, "theta_e")
        xi_ref = self._einstein_xi(theta_e / temperature_ref)
        result = np.zeros_like(temp)
        positive = temp > 0.0
        if np.any(positive):
            result[positive] = (
                alpha_ref
                * reference_k0
                * self._einstein_xi(theta_e / temp[positive])
                / xi_ref
            )
        return result

    def _solve_volume_scalar(
        self,
        model: PVTModel,
        reference: EOSParameters,
        thermal: TemperatureEOSParameters | None,
        coupling_parameters: Mapping[str, float],
        pressure: float,
        temperature: float,
    ) -> float:
        """Solve one scalar P--V--T inversion with adaptive bracketing."""
        if not np.isfinite(pressure):
            raise ValueError("pressure values must be finite")
        temp_array = np.asarray([temperature], dtype=np.float64)
        if model.coupling_family is PVTCouplingFamily.THERMAL_PRESSURE:
            centre = reference.V0
        else:
            thermal_spec = model.thermal_spec
            assert isinstance(thermal_spec, TemperatureEOSModel)
            assert thermal is not None
            centre = float(
                self._temperature.value(thermal_spec, thermal, temp_array)[0]
            )

        def residual(volume: float) -> float:
            return float(
                self.pressure(
                    model,
                    reference,
                    thermal,
                    coupling_parameters,
                    np.asarray([volume], dtype=np.float64),
                    temp_array,
                )[0]
                - pressure
            )

        f_centre = residual(centre)
        if f_centre == 0.0:
            return centre
        lower = centre
        upper = centre
        f_lower = f_centre
        f_upper = f_centre
        lower_active = True
        upper_active = True
        expansion = 1.25
        for _ in range(96):
            if lower_active:
                candidate = float(max(lower / expansion, _TINY))
                if candidate == lower:
                    lower_active = False
                else:
                    try:
                        f_candidate = residual(candidate)
                    except ValueError:
                        lower_active = False
                    else:
                        if f_candidate == 0.0:
                            return candidate
                        if f_candidate * f_lower < 0.0:
                            return float(
                                brentq(
                                    residual,
                                    candidate,
                                    lower,
                                    xtol=1.0e-12,
                                    rtol=1.0e-12,
                                )
                            )
                        lower = candidate
                        f_lower = f_candidate
            if upper_active:
                candidate = float(upper * expansion)
                if not np.isfinite(candidate):
                    upper_active = False
                else:
                    try:
                        f_candidate = residual(candidate)
                    except ValueError:
                        upper_active = False
                    else:
                        if f_candidate == 0.0:
                            return candidate
                        if f_upper * f_candidate < 0.0:
                            return float(
                                brentq(
                                    residual,
                                    upper,
                                    candidate,
                                    xtol=1.0e-12,
                                    rtol=1.0e-12,
                                )
                            )
                        upper = candidate
                        f_upper = f_candidate
            if not lower_active and not upper_active:
                break
        raise ValueError(
            f"could not bracket P-V-T volume root at P={pressure:g}, T={temperature:g}"
        )

    @classmethod
    def _resolve_mgd_parameters(
        cls,
        model: ThermalPressureModel,
        parameters: Mapping[str, float],
    ) -> MGDParameters:
        """Build validated MGD parameters from a coupling mapping."""
        q = None
        if model.mgd_variant is not None and model.mgd_variant.value == "full":
            q = cls._required_finite(parameters, "q")
        return MGDParameters(
            temperature_ref=cls._required_nonnegative(parameters, "temperature_ref"),
            theta_d0=cls._required_positive(parameters, "theta_d0"),
            gamma0=cls._required_finite(parameters, "gamma0"),
            q=q,
        )

    @staticmethod
    def _resolve_temperature_parameters(
        model: TemperatureEOSModel,
        parameters: TemperatureParameterLike | None,
    ) -> TemperatureEOSParameters:
        if parameters is None:
            raise ValueError("temperature parameters are required by this coupling")
        if isinstance(parameters, TemperatureEOSParameters):
            return parameters
        required = set(model.parameter_names)
        supplied = {str(name): float(value) for name, value in parameters.items()}
        missing = sorted(required.difference(supplied))
        if missing:
            raise ValueError(
                "missing temperature-EOS parameters: " + ", ".join(missing)
            )
        unknown = sorted(
            set(supplied).difference(TemperatureEOSParameters.__dataclass_fields__)
        )
        if unknown:
            raise ValueError(
                "unknown temperature-EOS parameters: " + ", ".join(unknown)
            )
        return TemperatureEOSParameters(**supplied)

    @staticmethod
    def _validate_reference_consistency(
        pressure: EOSParameters,
        thermal: TemperatureEOSParameters,
    ) -> None:
        if not math.isclose(
            pressure.V0,
            thermal.V0,
            rel_tol=1.0e-10,
            abs_tol=1.0e-12,
        ):
            raise ValueError(
                "pressure V0 and temperature V0 must describe the same reference volume"
            )

    def _pressure_vectorized(
        self,
        model: EOSModel,
        parameters: Mapping[str, np.ndarray | float],
        volume: np.ndarray,
    ) -> np.ndarray:
        """Evaluate a pressure EOS with state-dependent ``K0`` and ``V0``."""
        result = np.empty_like(volume)
        for index in np.ndindex(volume.shape):
            local = self._parameters_at_index(parameters, index)
            result[index] = self._pressure.pressure(model, local, float(volume[index]))
        return result

    def _bulk_modulus_vectorized(
        self,
        model: EOSModel,
        parameters: Mapping[str, np.ndarray | float],
        volume: np.ndarray,
    ) -> np.ndarray:
        """Evaluate bulk modulus with state-dependent ``K0`` and ``V0``."""
        result = np.empty_like(volume)
        for index in np.ndindex(volume.shape):
            local = self._parameters_at_index(parameters, index)
            result[index] = self._pressure.bulk_modulus(
                model, local, float(volume[index])
            )
        return result

    @staticmethod
    def _parameters_at_index(
        parameters: Mapping[str, np.ndarray | float],
        index: tuple[int, ...],
    ) -> dict[str, float]:
        """Return scalar EOS parameters for one broadcast state index."""
        local: dict[str, float] = {}
        for name, value in parameters.items():
            array = np.asarray(value, dtype=np.float64)
            local[name] = float(array[index]) if array.ndim else float(array)
        return local

    @staticmethod
    def _broadcast_state(
        first: ArrayLike,
        second: ArrayLike,
    ) -> tuple[np.ndarray, np.ndarray]:
        first_array = np.asarray(first, dtype=np.float64)
        second_array = np.asarray(second, dtype=np.float64)
        first_array, second_array = np.broadcast_arrays(first_array, second_array)
        if not np.all(np.isfinite(first_array)) or not np.all(
            np.isfinite(second_array)
        ):
            raise ValueError("P-V-T state coordinates must be finite")
        if np.any(second_array < 0.0):
            raise ValueError("temperature cannot be below absolute zero")
        return first_array.copy(), second_array.copy()

    @staticmethod
    def _validate_temperature(temperature: ArrayLike) -> np.ndarray:
        values = np.asarray(temperature, dtype=np.float64)
        if not np.all(np.isfinite(values)):
            raise ValueError("temperature values must be finite")
        if np.any(values < 0.0):
            raise ValueError("temperature cannot be below absolute zero")
        return values.copy()

    @staticmethod
    def _reference_bulk_modulus(
        parameters: PressureParameterLike,
    ) -> float:
        """Return and validate the reference ``K0`` without assuming an EOS."""
        if isinstance(parameters, EOSParameters):
            value = parameters.K0
        elif isinstance(parameters, Mapping):
            normalized = {
                str(name).upper(): float(item) for name, item in parameters.items()
            }
            if "K0" not in normalized:
                raise ValueError("pressure parameters must define K0")
            value = normalized["K0"]
        else:
            array = np.asarray(parameters, dtype=np.float64)
            if array.ndim != 1 or array.size < 1:
                raise ValueError(
                    "pressure parameters must form a one-dimensional array"
                )
            value = float(array[0])
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError("reference K0 must be finite and positive")
        return float(value)

    @staticmethod
    def _required_finite(parameters: Mapping[str, float], name: str) -> float:
        if name not in parameters:
            raise ValueError(f"missing P-V-T coupling parameter: {name}")
        value = float(parameters[name])
        if not np.isfinite(value):
            raise ValueError(f"P-V-T coupling parameter {name} must be finite")
        return value

    @classmethod
    def _required_nonnegative(cls, parameters: Mapping[str, float], name: str) -> float:
        """Return a required finite non-negative parameter."""
        value = cls._required_finite(parameters, name)
        if value < 0.0:
            raise ValueError(f"coupling parameter {name!r} cannot be negative")
        return value

    @classmethod
    def _required_positive(cls, parameters: Mapping[str, float], name: str) -> float:
        value = cls._required_finite(parameters, name)
        if value <= 0.0:
            raise ValueError(f"P-V-T coupling parameter {name} must be positive")
        return value

    @staticmethod
    def _einstein_occupation(theta: float, temperature: np.ndarray) -> np.ndarray:
        result = np.zeros_like(temperature)
        positive = temperature > 0.0
        if np.any(positive):
            argument = theta / temperature[positive]
            moderate = argument < 700.0
            values = np.zeros_like(argument)
            values[moderate] = 1.0 / np.expm1(argument[moderate])
            result[positive] = values
        return result

    @staticmethod
    def _einstein_xi(argument: np.ndarray | float) -> np.ndarray | float:
        values = np.asarray(argument, dtype=np.float64)
        result = np.empty_like(values)
        moderate = values < 350.0
        if np.any(moderate):
            expm1 = np.expm1(values[moderate])
            result[moderate] = values[moderate] ** 2 * (expm1 + 1.0) / expm1**2
        if np.any(~moderate):
            result[~moderate] = values[~moderate] ** 2 * np.exp(-values[~moderate])
        if np.isscalar(argument):
            return float(result)
        return result


__all__ = [
    "PVTCouplingFamily",
    "PVTEOS",
    "PVTModel",
    "available_pvt_couplings",
    "parse_pvt_coupling",
]
