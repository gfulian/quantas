# -*- coding: utf-8 -*-

"""Fitting adapters and parameter contracts for P--V--T equations."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from quantas.core.math.fitting import (
    BaseFitModel,
    ParameterDefinition,
    ParameterMap,
    ParameterState,
)
from quantas.core.physics.eos import (
    EOSFamily,
    EOSModel,
    PVTCouplingFamily,
    PVTEOS,
    PVTModel,
    TemperatureEOSFamily,
    TemperatureEOSModel,
    TemperatureEOSVariant,
    ThermalPressureFamily,
    MGDVariant,
    implied_kp,
    implied_kpp,
)

from .pv import estimate_pressure_parameters
from ..models import ParameterConstraint
from .vt import (
    estimate_temperature_parameters,
    temperature_parameter_names,
)

_POSITIVE = float(np.finfo(np.float64).tiny)
_DEFAULT_TREF = 298.15
_DEFAULT_THETA = 500.0


class PVTEOSFitModel(BaseFitModel):
    r"""Adapt a compositional P--V--T equation to the general fit contract.

    The explanatory coordinate is a two-row matrix with volume in the first
    row and temperature in the second.  The response is pressure.  Thus the
    model implements

    .. math::

        P_i = P(V_i,T_i;\boldsymbol\beta).

    Parameters
    ----------
    model : PVTModel
        Reference pressure EOS, optional V--T model, and coupling family.
    initial_parameters : mapping, optional
        Complete initial parameter values used by the unmapped API.
    """

    def __init__(
        self,
        model: PVTModel,
        *,
        initial_parameters: Mapping[str, float] | None = None,
    ) -> None:
        if not isinstance(model, PVTModel):
            raise TypeError("P-V-T fitting requires a PVTModel")
        self.pvt_model = model
        self._core = PVTEOS()
        defaults = _default_values(model)
        if initial_parameters is not None:
            defaults.update(
                {str(name): float(value) for name, value in initial_parameters.items()}
            )
        self._initial = np.asarray(
            [defaults[name] for name in pvt_parameter_names(model)],
            dtype=np.float64,
        )

    @property
    def name(self) -> str:
        """Return a stable composite model identifier."""
        return f"pvt_eos:{self.pvt_model.tag}"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return complete P--V--T parameter reporting order."""
        return pvt_parameter_names(self.pvt_model)

    @property
    def core(self) -> PVTEOS:
        """Return the pure compositional P--V--T evaluator."""
        return self._core

    def evaluate(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Evaluate pressure at paired volume and temperature coordinates."""
        volume, temperature = _split_coordinates(x)
        pressure, thermal, coupling = self.split_parameters(parameters)
        return self._core.pressure(
            self.pvt_model,
            pressure,
            thermal,
            coupling,
            volume,
            temperature,
        )

    def derivative_x(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        r"""Return derivatives with respect to volume and temperature.

        The volume derivative is analytical,

        .. math::

            \left(\frac{\partial P}{\partial V}\right)_T=-\frac{K_T}{V}.

        The temperature derivative is analytical for thermal pressure and a
        stable scale-aware numerical derivative for the other compositional
        couplings.  It represents :math:`(\partial P/\partial T)_V`, the
        quantity required to project temperature uncertainty into pressure.
        """
        volume, temperature = _split_coordinates(x)
        pressure, thermal, coupling = self.split_parameters(parameters)
        bulk = self._core.bulk_modulus(
            self.pvt_model,
            pressure,
            thermal,
            coupling,
            volume,
            temperature,
        )
        d_pressure_d_volume = -bulk / volume
        if self.pvt_model.coupling_family is PVTCouplingFamily.THERMAL_PRESSURE:
            d_pressure_d_temperature = (
                self._core.thermal_pressure_temperature_derivative_for_model(
                    self.pvt_model,
                    pressure,
                    coupling,
                    volume,
                    temperature,
                )
            )
        else:
            step = np.cbrt(np.finfo(np.float64).eps) * np.maximum(temperature, 1.0)
            lower_temperature = np.maximum(temperature - step, 0.0)
            upper_temperature = temperature + step
            upper = self._core.pressure(
                self.pvt_model,
                pressure,
                thermal,
                coupling,
                volume,
                upper_temperature,
            )
            lower = self._core.pressure(
                self.pvt_model,
                pressure,
                thermal,
                coupling,
                volume,
                lower_temperature,
            )
            d_pressure_d_temperature = (upper - lower) / (
                upper_temperature - lower_temperature
            )
        result = np.vstack([d_pressure_d_volume, d_pressure_d_temperature])
        if not np.all(np.isfinite(result)):
            raise ValueError("P-V-T coordinate derivatives are not finite")
        return result

    def initial_guess(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return configured complete initial parameters after validation."""
        volume, temperature = _split_coordinates(x)
        pressure = np.asarray(y, dtype=np.float64)
        _validate_pvt_data(volume, temperature, pressure)
        return self._initial.copy()

    def bounds(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return minimally restrictive physical bounds."""
        self.initial_guess(x, y)
        lower: list[float] = []
        upper: list[float] = []
        for name in self.parameter_names:
            low, high = _default_bounds(name)
            lower.append(low)
            upper.append(high)
        return np.asarray(lower), np.asarray(upper)

    def split_parameters(
        self,
        parameters: np.ndarray | Sequence[float],
    ) -> tuple[dict[str, float], dict[str, float] | None, dict[str, float]]:
        """Split a complete vector into pressure, thermal, and coupling parts."""
        array = np.asarray(parameters, dtype=np.float64)
        names = self.parameter_names
        if array.ndim != 1 or array.size != len(names):
            raise ValueError("P-V-T parameter vector has an invalid shape")
        values = dict(zip(names, array, strict=True))
        pressure = {
            "K0": values["K0"],
            "KP": values["KP"],
            "KPP": values["KPP"],
            "V0": values["V0"],
        }
        thermal: dict[str, float] | None
        thermal_model = self.pvt_model.thermal_spec
        if thermal_model is None:
            thermal = None
        else:
            assert isinstance(thermal_model, TemperatureEOSModel)
            thermal = {
                "V0": values["V0"],
                "temperature_ref": values["temperature_ref"],
            }
            for name in temperature_parameter_names(thermal_model):
                if name in {"V0", "temperature_ref", "kp"}:
                    continue
                thermal[name] = values[name]
            if thermal_model.family is TemperatureEOSFamily.KROLL_HOLLAND_POWELL:
                thermal["kp"] = values["KP"]
        coupling_family = self.pvt_model.coupling_family
        assert isinstance(coupling_family, PVTCouplingFamily)
        if coupling_family is PVTCouplingFamily.LINEAR_BULK_MODULUS:
            coupling = {"dK0_dT": values["dK0_dT"]}
        elif coupling_family is PVTCouplingFamily.ANDERSON_GRUNEISEN:
            coupling = {"delta": values["delta"]}
        else:
            thermal_pressure = self.pvt_model.thermal_pressure_spec
            assert thermal_pressure is not None
            if (
                thermal_pressure.family_name
                is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE
            ):
                coupling = {
                    "temperature_ref": values["temperature_ref"],
                    "theta_d0": values["theta_d0"],
                    "gamma0": values["gamma0"],
                }
                if thermal_pressure.mgd_variant is MGDVariant.FULL:
                    coupling["q"] = values["q"]
            else:
                coupling = {
                    "temperature_ref": values["temperature_ref"],
                    "alpha_ref": values["alpha_ref"],
                    "theta_e": values["theta_e"],
                }
        return pressure, thermal, coupling

    def metadata(self) -> dict[str, Any]:
        """Return model composition and coordinate semantics."""
        return {
            **super().metadata(),
            "pvt_model": self.pvt_model.as_dict(),
            "relationship": "pressure(volume,temperature)",
            "coordinate_order": ["volume", "temperature"],
            "temperature_dependent_parameters": ["V0(T)", "K0(T)"],
            "temperature_invariant_parameters": ["KP", "KPP"],
        }


def pvt_parameter_names(model: PVTModel) -> tuple[str, ...]:
    """Return complete reporting order for one P--V--T model."""
    names: list[str] = ["K0", "KP", "KPP", "V0", "temperature_ref"]
    thermal = model.thermal_spec
    if thermal is not None:
        assert isinstance(thermal, TemperatureEOSModel)
        names.extend(
            name
            for name in temperature_parameter_names(thermal)
            if name not in {"V0", "temperature_ref", "kp"}
        )
    coupling = model.coupling_family
    assert isinstance(coupling, PVTCouplingFamily)
    if coupling is PVTCouplingFamily.LINEAR_BULK_MODULUS:
        names.append("dK0_dT")
    elif coupling is PVTCouplingFamily.ANDERSON_GRUNEISEN:
        names.append("delta")
    else:
        thermal_pressure = model.thermal_pressure_spec
        assert thermal_pressure is not None
        if thermal_pressure.family_name is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE:
            names.extend(["theta_d0", "gamma0"])
            if thermal_pressure.mgd_variant is MGDVariant.FULL:
                names.append("q")
        else:
            names.extend(["alpha_ref", "theta_e"])
    return tuple(dict.fromkeys(names))


def estimate_pvt_parameters(
    model: PVTModel,
    volume: np.ndarray | Sequence[float],
    temperature: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
    *,
    reference_temperature: float | None = None,
) -> dict[str, float]:
    """Estimate complete parameters from a P--V--T table.

    The reference P--V estimate uses observations closest to the reference
    temperature.  The V--T estimate uses observations at the smallest absolute
    pressures.  These subsets provide initialization only; the subsequent fit
    uses every selected observation.
    """
    v, t, p = _validate_pvt_data(volume, temperature, pressure)
    thermal = model.thermal_spec
    tref = (
        0.0
        if thermal is not None and thermal.family is TemperatureEOSFamily.SALJE
        else float(
            _DEFAULT_TREF if reference_temperature is None else reference_temperature
        )
    )
    minimum = min(v.size, max(4, int(np.ceil(0.25 * v.size))))
    isotherm_indices = np.argsort(np.abs(t - tref))[:minimum]
    try:
        pressure_estimate = estimate_pressure_parameters(
            v[isotherm_indices], p[isotherm_indices]
        )
    except ValueError:
        pressure_estimate = estimate_pressure_parameters(v, p)

    estimates = {**pressure_estimate, "temperature_ref": tref}
    if thermal is not None:
        isobar_indices = np.argsort(np.abs(p))[:minimum]
        try:
            thermal_estimate = estimate_temperature_parameters(
                thermal,
                t[isobar_indices],
                v[isobar_indices],
                reference_temperature=tref,
            )
        except (ValueError, np.linalg.LinAlgError):
            thermal_estimate = estimate_temperature_parameters(
                thermal,
                t,
                v,
                reference_temperature=tref,
            )
        estimates["V0"] = float(thermal_estimate["V0"])
        for name, value in thermal_estimate.items():
            if name not in {"V0", "temperature_ref", "kp"}:
                estimates[name] = float(value)
    coupling = model.coupling_family
    assert isinstance(coupling, PVTCouplingFamily)
    if coupling is PVTCouplingFamily.LINEAR_BULK_MODULUS:
        estimates["dK0_dT"] = -1.0e-4 * estimates["K0"]
    elif coupling is PVTCouplingFamily.ANDERSON_GRUNEISEN:
        estimates["delta"] = max(float(estimates["KP"]), 1.0)
    else:
        thermal_pressure = model.thermal_pressure_spec
        assert thermal_pressure is not None
        if thermal_pressure.family_name is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE:
            estimates.setdefault("theta_d0", _DEFAULT_THETA)
            estimates.setdefault("gamma0", 1.0)
            if thermal_pressure.mgd_variant is MGDVariant.FULL:
                estimates.setdefault("q", 1.0)
        else:
            estimates.setdefault("alpha_ref", 3.0e-5)
            estimates.setdefault("theta_e", _DEFAULT_THETA)
    return {name: float(estimates[name]) for name in pvt_parameter_names(model)}


def build_pvt_parameter_map(
    model: PVTModel,
    volume: np.ndarray | Sequence[float],
    temperature: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
    constraints: Sequence[ParameterConstraint] = (),
    *,
    pressure_unit: str = "GPa",
    volume_unit: str = "angstrom^3",
    reference_temperature: float | None = None,
) -> ParameterMap:
    """Build FREE/FIXED/IMPLIED parameters for a global P--V--T fit."""
    estimates = estimate_pvt_parameters(
        model,
        volume,
        temperature,
        pressure,
        reference_temperature=reference_temperature,
    )
    overrides = _constraint_overrides(constraints)
    allowed = set(pvt_parameter_names(model))
    unknown = sorted(set(overrides) - allowed)
    if unknown:
        raise ValueError("unknown P-V-T parameter constraint(s): " + ", ".join(unknown))
    definitions = tuple(
        _definition(
            model,
            name,
            estimates,
            overrides.get(name),
            pressure_unit=pressure_unit,
            volume_unit=volume_unit,
        )
        for name in pvt_parameter_names(model)
    )
    return ParameterMap(definitions, resolver=_resolver(model))


def _definition(
    model: PVTModel,
    name: str,
    estimates: Mapping[str, float],
    override: ParameterConstraint | None,
    *,
    pressure_unit: str,
    volume_unit: str,
) -> ParameterDefinition:
    implied = _is_implied(model, name)
    default_fixed = _is_default_fixed(model, name)
    unit = _unit(name, pressure_unit, volume_unit)
    description = _description(name)
    if implied:
        if override is not None and override.state is not ParameterState.IMPLIED:
            raise ValueError(
                f"parameter {name} is implied by {model.tag} and cannot be "
                f"declared {override.state.value}"
            )
        return ParameterDefinition.implied(
            name,
            unit=unit,
            description=description,
            metadata={
                "source": "model_composition",
                "default_state": ParameterState.IMPLIED.value,
            },
        )
    lower, upper = _default_bounds(name)
    if override is None and default_fixed:
        return ParameterDefinition.fixed(
            name,
            estimates[name],
            lower_bound=lower,
            upper_bound=upper,
            unit=unit,
            description=description,
            metadata={
                "source": "model_default",
                "default_state": ParameterState.FIXED.value,
            },
        )
    if override is None:
        return ParameterDefinition.free(
            name,
            estimates[name],
            lower_bound=lower,
            upper_bound=upper,
            unit=unit,
            description=description,
            metadata={
                "initial_source": "pvt_subset_estimate",
                "default_state": ParameterState.FREE.value,
            },
        )
    selected_lower = (
        lower if np.isneginf(override.lower_bound) else override.lower_bound
    )
    selected_upper = (
        upper if np.isposinf(override.upper_bound) else override.upper_bound
    )
    if override.state is ParameterState.FREE:
        initial = (
            estimates[name]
            if override.initial_value is None
            else override.initial_value
        )
        return ParameterDefinition.free(
            name,
            initial,
            lower_bound=selected_lower,
            upper_bound=selected_upper,
            unit=override.unit or unit,
            description=override.description or description,
            metadata={
                **override.metadata,
                "initial_source": "user",
                "default_state": (
                    ParameterState.FIXED.value
                    if default_fixed
                    else ParameterState.FREE.value
                ),
            },
        )
    if override.state is ParameterState.FIXED:
        if override.value is None:
            raise ValueError(f"fixed parameter {name} requires a value")
        return ParameterDefinition.fixed(
            name,
            override.value,
            lower_bound=selected_lower,
            upper_bound=selected_upper,
            unit=override.unit or unit,
            description=override.description or description,
            metadata={
                **override.metadata,
                "source": "user_fixed",
                "default_state": (
                    ParameterState.FIXED.value
                    if default_fixed
                    else ParameterState.FREE.value
                ),
            },
        )
    raise ValueError(f"parameter {name} cannot be declared {override.state.value}")


def _resolver(model: PVTModel):
    """Return implied pressure and thermal parameter resolution."""
    pressure_model = model.pressure_spec
    assert isinstance(pressure_model, EOSModel)
    thermal_model = model.thermal_spec

    def resolver(values: Mapping[str, float]) -> Mapping[str, float]:
        resolved: dict[str, float] = {}
        if pressure_model.order == 2:
            kp = implied_kp(pressure_model)
            if kp is None:
                raise ValueError(f"{pressure_model.tag} does not imply KP")
            resolved["KP"] = kp
        current_kp = resolved.get("KP", values.get("KP"))
        if current_kp is None:
            raise ValueError("P-V-T parameter map cannot resolve KP")
        if "KPP" in pvt_parameter_names(model) and _pressure_kpp_implied(
            pressure_model
        ):
            resolved["KPP"] = implied_kpp(
                pressure_model,
                float(values["K0"]),
                float(current_kp),
            )
        if thermal_model is not None:
            if (
                thermal_model.family is TemperatureEOSFamily.BERMAN
                and thermal_model.variant is TemperatureEOSVariant.LINEAR
            ):
                resolved["alpha1"] = 0.0
            elif (
                thermal_model.family is TemperatureEOSFamily.FEI
                and thermal_model.variant is TemperatureEOSVariant.LINEAR
            ):
                resolved["alpha2"] = 0.0
            elif (
                thermal_model.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL
                and thermal_model.variant is TemperatureEOSVariant.SIMPLIFIED
            ):
                resolved["alpha1"] = 0.0
        return resolved

    return resolver


def _is_implied(model: PVTModel, name: str) -> bool:
    pressure_model = model.pressure_spec
    assert isinstance(pressure_model, EOSModel)
    if name == "KP" and pressure_model.order == 2:
        return True
    if name == "KPP" and _pressure_kpp_implied(pressure_model):
        return True
    thermal = model.thermal_spec
    if thermal is None:
        return False
    return (
        (
            name == "alpha1"
            and thermal.family is TemperatureEOSFamily.BERMAN
            and thermal.variant is TemperatureEOSVariant.LINEAR
        )
        or (
            name == "alpha2"
            and thermal.family is TemperatureEOSFamily.FEI
            and thermal.variant is TemperatureEOSVariant.LINEAR
        )
        or (
            name == "alpha1"
            and thermal.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL
            and thermal.variant is TemperatureEOSVariant.SIMPLIFIED
        )
    )


def _pressure_kpp_implied(model: EOSModel) -> bool:
    """Return whether KPP is fixed by family/order rather than refined."""
    return model.family is EOSFamily.MURNAGHAN or model.order in {2, 3}


def _is_default_fixed(model: PVTModel, name: str) -> bool:
    if name == "temperature_ref":
        return True
    thermal = model.thermal_spec
    if (
        name in {"theta_e", "theta_sat"}
        and thermal is not None
        and thermal.family
        in {
            TemperatureEOSFamily.KROLL_HOLLAND_POWELL,
            TemperatureEOSFamily.SALJE,
        }
    ):
        return True
    if model.coupling_family is not PVTCouplingFamily.THERMAL_PRESSURE:
        return False
    thermal_pressure = model.thermal_pressure_spec
    assert thermal_pressure is not None
    if thermal_pressure.family_name is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE:
        return name in {"theta_d0", "q"}
    return name == "theta_e"


def _default_values(model: PVTModel) -> dict[str, float]:
    values = {
        "K0": 100.0,
        "KP": 4.0,
        "KPP": 0.0,
        "V0": 100.0,
        "temperature_ref": 0.0
        if model.thermal_spec is not None
        and model.thermal_spec.family is TemperatureEOSFamily.SALJE
        else _DEFAULT_TREF,
        "alpha0": 1.0e-5,
        "alpha1": 0.0,
        "alpha2": 0.0,
        "alpha_ref": 1.0e-5,
        "p1": 1.0e-5,
        "theta_sat": 300.0,
        "theta_e": _DEFAULT_THETA,
        "theta_d0": _DEFAULT_THETA,
        "gamma0": 1.0,
        "q": 1.0,
        "dK0_dT": -0.01,
        "delta": 4.0,
    }
    return values


def _default_bounds(name: str) -> tuple[float, float]:
    if name in {
        "K0",
        "V0",
        "theta_e",
        "theta_d0",
        "theta_sat",
        "delta",
    }:
        return _POSITIVE, np.inf
    if name == "temperature_ref":
        return 0.0, np.inf
    return -np.inf, np.inf


def _unit(name: str, pressure_unit: str, volume_unit: str) -> str | None:
    return {
        "K0": pressure_unit,
        "KP": "1",
        "KPP": f"{pressure_unit}^-1",
        "V0": volume_unit,
        "temperature_ref": "K",
        "alpha0": "K^-1",
        "alpha1": "K^-2",
        "alpha2": "K",
        "alpha_ref": "K^-1",
        "p1": f"{volume_unit}^(1/3) K^-1",
        "theta_sat": "K",
        "theta_e": "K",
        "theta_d0": "K",
        "gamma0": "1",
        "q": "1",
        "dK0_dT": f"{pressure_unit} K^-1",
        "delta": "1",
    }.get(name)


def _description(name: str) -> str:
    return {
        "K0": "reference-temperature zero-pressure bulk modulus",
        "KP": "reference pressure derivative of the bulk modulus",
        "KPP": "reference second pressure derivative of the bulk modulus",
        "V0": "zero-pressure volume at the reference temperature",
        "temperature_ref": "fixed reference temperature",
        "alpha0": "thermal expansion coefficient parameter alpha0",
        "alpha1": "thermal expansion coefficient parameter alpha1",
        "alpha2": "Fei inverse-square thermal coefficient alpha2",
        "alpha_ref": "exact thermal expansion at the reference temperature",
        "p1": "Salje expansion parameter",
        "theta_sat": "Salje saturation temperature",
        "theta_e": "Einstein temperature",
        "theta_d0": "Debye temperature at the reference volume",
        "gamma0": "thermal Gruneisen parameter at the reference volume",
        "q": "logarithmic volume exponent of the Gruneisen parameter",
        "dK0_dT": "linear temperature derivative of K0",
        "delta": "Anderson-Gruneisen parameter",
    }.get(name, name)


def _constraint_overrides(
    constraints: Sequence[ParameterConstraint],
) -> dict[str, ParameterConstraint]:
    aliases = {
        "v0": "V0",
        "k0": "K0",
        "kp": "KP",
        "kprime": "KP",
        "kpp": "KPP",
        "tref": "temperature_ref",
        "d_k0_dt": "dK0_dT",
        "dk0dt": "dK0_dT",
        "dkdt": "dK0_dT",
        "thetae": "theta_e",
        "thetad": "theta_d0",
        "thetad0": "theta_d0",
        "debyetemperature": "theta_d0",
        "gamma": "gamma0",
        "gamma_0": "gamma0",
        "gruneisen": "gamma0",
        "thetasat": "theta_sat",
        "alpharef": "alpha_ref",
    }
    overrides: dict[str, ParameterConstraint] = {}
    for constraint in constraints:
        key = str(constraint.name).strip().replace("-", "_")
        normalized = aliases.get(key.lower(), key)
        if normalized in overrides:
            raise ValueError(f"duplicate P-V-T parameter constraint: {normalized}")
        overrides[normalized] = (
            constraint
            if constraint.name == normalized
            else ParameterConstraint(
                name=normalized,
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


def _split_coordinates(
    x: np.ndarray | Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    coordinates = np.asarray(x, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[0] != 2:
        raise ValueError(
            "P-V-T coordinates must have shape (2, n): volume, temperature"
        )
    volume = coordinates[0]
    temperature = coordinates[1]
    if np.any(volume <= 0.0):
        raise ValueError("P-V-T volumes must be strictly positive")
    if np.any(temperature < 0.0):
        raise ValueError("P-V-T temperatures cannot be negative")
    return volume, temperature


def _validate_pvt_data(
    volume: np.ndarray | Sequence[float],
    temperature: np.ndarray | Sequence[float],
    pressure: np.ndarray | Sequence[float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    v = np.asarray(volume, dtype=np.float64)
    t = np.asarray(temperature, dtype=np.float64)
    p = np.asarray(pressure, dtype=np.float64)
    if v.ndim != 1 or t.ndim != 1 or p.ndim != 1:
        raise ValueError("P-V-T data arrays must be one-dimensional")
    if v.shape != t.shape or v.shape != p.shape or v.size < 4:
        raise ValueError("P-V-T data require equal arrays with at least four points")
    if (
        not np.all(np.isfinite(v))
        or not np.all(np.isfinite(t))
        or not np.all(np.isfinite(p))
    ):
        raise ValueError("P-V-T data must be finite")
    if np.any(v <= 0.0) or np.any(t < 0.0):
        raise ValueError("P-V-T volumes must be positive and temperatures non-negative")
    return v.copy(), t.copy(), p.copy()


__all__ = [
    "PVTEOSFitModel",
    "build_pvt_parameter_map",
    "estimate_pvt_parameters",
    "pvt_parameter_names",
]
