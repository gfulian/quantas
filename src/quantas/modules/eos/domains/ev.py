# -*- coding: utf-8 -*-

"""Energy-volume fit models and parameter preparation for EOS workflows.

This module adapts the integrated :mod:`quantas.core.physics.eos` equations to
public EOS workflow units.  The numerical core remains unit-consistent and
works with an energy-density bulk modulus; the workflow reports the physical
parameters conventionally used by solid-state calculations: ``E0`` in energy
units, ``V0`` in volume units, ``K0`` in pressure units, ``KP`` dimensionless,
and ``KPP`` in inverse-pressure units.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import Any

import numpy as np

from quantas.core.math.fitting import (
    BaseFitModel,
    ParameterDefinition,
    ParameterMap,
    ParameterState,
)
from quantas.core.physics.eos import (
    EOSModel,
    EnergyEOS,
    EnergyEOSFitModel as CoreEnergyEOSFitModel,
    implied_kp,
    implied_kpp,
    parse_eos_model,
    resolve_energy_parameters,
)
from quantas.core.physics.units import energy_to_pressure, pressure_to_energy

from ..models import ParameterConstraint

_ENERGY_PARAMETER_ORDER = ("E0", "K0", "KP", "KPP", "V0")
_POSITIVE_LOWER_BOUND = float(np.nextafter(0.0, 1.0))


class EnergyEOSFitModel(BaseFitModel):
    """Adapt an integrated energy EOS to the public EOS fitting contract.

    Parameters
    ----------
    model : EOSModel or str
        Integrated EOS family and order.
    energy_unit : str, optional
        Public energy unit.  The normalized EOS workflow currently uses Ha.
    volume_unit : str, optional
        Public volume unit.  The normalized theoretical workflow currently
        uses ``angstrom^3``.
    pressure_unit : str, optional
        Public pressure unit used for ``K0`` and ``KPP``.

    Notes
    -----
    The shared :class:`~quantas.core.physics.eos.EnergyEOS` evaluates energy
    with ``K0`` expressed as energy per volume.  Conversion between that
    natural representation and pressure units is intentionally confined to
    this adapter.
    """

    def __init__(
        self,
        model: EOSModel | str,
        *,
        energy_unit: str = "Ha",
        volume_unit: str = "angstrom^3",
        pressure_unit: str = "GPa",
    ) -> None:
        self.eos_model = parse_eos_model(model)
        if not self.eos_model.supports_energy:
            raise ValueError(f"{self.eos_model.tag} has no integrated E-V form")
        self.energy_unit = str(energy_unit)
        self.volume_unit = str(volume_unit)
        self.pressure_unit = str(pressure_unit)
        if self.energy_unit != "Ha" or self.volume_unit != "angstrom^3":
            raise ValueError(
                "public E-V fitting currently requires normalized Ha and angstrom^3"
            )
        self._energy = EnergyEOS()

    @property
    def name(self) -> str:
        """Return a stable physical-model identifier."""
        return f"energy_eos:{self.eos_model.tag}"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return complete public physical-parameter order."""
        return _ENERGY_PARAMETER_ORDER

    def evaluate(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Evaluate energy at the supplied volumes.

        Parameters
        ----------
        x : array-like
            Positive volumes in ``volume_unit``.
        parameters : array-like
            Complete ``E0, K0, KP, KPP, V0`` vector in public units.

        Returns
        -------
        ndarray
            Calculated energies in ``energy_unit``.
        """
        volume = _validate_volume(x)
        physical = self.core_parameters(parameters)
        return np.asarray(
            self._energy.evaluate(self.eos_model, volume, physical),
            dtype=np.float64,
        )

    def derivative_x(
        self,
        x: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        r"""Return analytical :math:`\partial E/\partial V` in energy density.

        The derivative follows directly from :math:`P=-\partial E/\partial V`.
        It is returned in the natural ``energy_unit / volume_unit`` numerical
        units required by the generic effective-variance and ODR services.

        Parameters
        ----------
        x : np.ndarray | Sequence[float]
            Positive volumes in ``volume_unit``.
        parameters : np.ndarray | Sequence[float]
            Complete ``E0, K0, KP, KPP, V0`` vector in public workflow units.

        Returns
        -------
        np.ndarray
            Analytical :math:`\partial E/\partial V` in energy density.
        """
        volume = _validate_volume(x)
        physical = self.core_parameters(parameters)
        pressure_density = self._energy.pressure(
            self.eos_model,
            physical,
            volume,
        )
        return -np.asarray(pressure_density, dtype=np.float64)

    def pressure(
        self,
        volume: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return pressure derived from the fitted energy EOS in public units.

        Parameters
        ----------
        volume : np.ndarray | Sequence[float]
            Positive volumes in ``volume_unit``.
        parameters : np.ndarray | Sequence[float]
            Complete ``E0, K0, KP, KPP, V0`` vector in public workflow units.

        Returns
        -------
        np.ndarray
            Pressure derived from the fitted energy EOS in public units.
        """
        values = _validate_volume(volume)
        density = self._energy.pressure(
            self.eos_model,
            self.core_parameters(parameters),
            values,
        )
        return np.asarray(
            energy_to_pressure(
                density,
                self.energy_unit,
                "angstrom",
                self.pressure_unit,
            ),
            dtype=np.float64,
        )

    def bulk_modulus(
        self,
        volume: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return the instantaneous bulk modulus in public pressure units.

        Parameters
        ----------
        volume : np.ndarray | Sequence[float]
            Positive volumes in ``volume_unit``.
        parameters : np.ndarray | Sequence[float]
            Complete ``E0, K0, KP, KPP, V0`` vector in public workflow units.

        Returns
        -------
        np.ndarray
            The instantaneous bulk modulus in public pressure units.
        """
        values = _validate_volume(volume)
        core = self.core_parameters(parameters)
        density = self._energy.pressure(self.eos_model, core, values)
        del density  # pressure evaluation validates the same physical state
        from quantas.core.physics.eos import PressureEOS

        bulk_density = PressureEOS().bulk_modulus(self.eos_model, core, values)
        return np.asarray(
            energy_to_pressure(
                bulk_density,
                self.energy_unit,
                "angstrom",
                self.pressure_unit,
            ),
            dtype=np.float64,
        )

    def bulk_modulus_derivative(
        self,
        volume: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return :math:`K'(V)` in dimensionless form.

        Parameters
        ----------
        volume : np.ndarray | Sequence[float]
            Positive volumes in ``volume_unit``.
        parameters : np.ndarray | Sequence[float]
            Complete ``E0, K0, KP, KPP, V0`` vector in public workflow units.

        Returns
        -------
        np.ndarray
            :math:`K'(V)` in dimensionless form.
        """
        from quantas.core.physics.eos import PressureEOS

        values = _validate_volume(volume)
        return np.asarray(
            PressureEOS().bulk_modulus_derivative(
                self.eos_model,
                self.core_parameters(parameters),
                values,
            ),
            dtype=np.float64,
        )

    def bulk_modulus_second_derivative(
        self,
        volume: np.ndarray | Sequence[float],
        parameters: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return :math:`K''(V)` in inverse public pressure units.

        Parameters
        ----------
        volume : np.ndarray | Sequence[float]
            Positive volumes in ``volume_unit``.
        parameters : np.ndarray | Sequence[float]
            Complete ``E0, K0, KP, KPP, V0`` vector in public workflow units.

        Returns
        -------
        np.ndarray
            :math:`K''(V)` in inverse public pressure units.
        """
        from quantas.core.physics.eos import PressureEOS

        values = _validate_volume(volume)
        internal = np.asarray(
            PressureEOS().bulk_modulus_second_derivative(
                self.eos_model,
                self.core_parameters(parameters),
                values,
            ),
            dtype=np.float64,
        )
        factor = _pressure_per_energy_density(
            self.energy_unit,
            self.pressure_unit,
        )
        return internal / factor

    def initial_guess(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> np.ndarray:
        """Return a complete public initial physical parameter vector.

        Parameters
        ----------
        x : np.ndarray | Sequence[float]
            Positive sampled volumes in ``volume_unit``.
        y : np.ndarray | Sequence[float]
            Sampled energies in ``energy_unit`` aligned with ``x``.

        Returns
        -------
        np.ndarray
            A complete public initial physical parameter vector.
        """
        estimates = estimate_energy_parameters(
            self.eos_model,
            x,
            y,
            energy_unit=self.energy_unit,
            pressure_unit=self.pressure_unit,
        )
        return np.asarray(
            [estimates[name] for name in _ENERGY_PARAMETER_ORDER],
            dtype=np.float64,
        )

    def bounds(
        self,
        x: np.ndarray | Sequence[float],
        y: np.ndarray | Sequence[float],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return minimally restrictive bounds for public parameters.

        Parameters
        ----------
        x : np.ndarray | Sequence[float]
            Positive sampled volumes in ``volume_unit``.
        y : np.ndarray | Sequence[float]
            Sampled energies in ``energy_unit`` aligned with ``x``.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Minimally restrictive bounds for public parameters.
        """
        _validate_energy_data(x, y)
        lower = np.asarray(
            [-np.inf, _POSITIVE_LOWER_BOUND, -np.inf, -np.inf, _POSITIVE_LOWER_BOUND],
            dtype=np.float64,
        )
        return lower, np.full(5, np.inf, dtype=np.float64)

    def core_parameters(
        self,
        parameters: np.ndarray | Sequence[float] | Mapping[str, float],
    ) -> dict[str, float]:
        """Convert complete public parameters to core energy-density units.

        Parameters
        ----------
        parameters : np.ndarray | Sequence[float] | Mapping[str, float]
            Complete public ``E0, K0, KP, KPP, V0`` parameters.

        Returns
        -------
        dict[str, float]
            Core parameter mapping with ``K0`` expressed as energy density and
            ``KPP`` expressed in the reciprocal of that energy-density unit.
        """
        public = _parameter_mapping(parameters)
        factor = _pressure_per_energy_density(
            self.energy_unit,
            self.pressure_unit,
        )
        return {
            "E0": public["E0"],
            "K0": float(
                pressure_to_energy(
                    public["K0"],
                    self.energy_unit,
                    "angstrom",
                    self.pressure_unit,
                )
            ),
            "KP": public["KP"],
            "KPP": public["KPP"] * factor,
            "V0": public["V0"],
        }

    def metadata(self) -> dict[str, Any]:
        """Return model, relationship, and public-unit metadata.

        Returns
        -------
        dict[str, Any]
            Model, relationship, and public-unit metadata.
        """
        return {
            **super().metadata(),
            "eos_model": self.eos_model.as_dict(),
            "relationship": "energy(volume)",
            "pressure_relation": "P(V)=-dE/dV",
            "parameter_units": {
                "E0": self.energy_unit,
                "K0": self.pressure_unit,
                "KP": "1",
                "KPP": f"{self.pressure_unit}^-1",
                "V0": self.volume_unit,
            },
        }


def build_energy_parameter_map(
    model: EOSModel | str,
    volume: np.ndarray | Sequence[float],
    energy: np.ndarray | Sequence[float],
    constraints: Sequence[ParameterConstraint] = (),
    *,
    energy_unit: str = "Ha",
    pressure_unit: str = "GPa",
    volume_unit: str = "angstrom^3",
) -> ParameterMap:
    """Build the reduced/full parameter mapping for one E-V fit.

    Parameters
    ----------
    model : EOSModel or str
        Integrated EOS family and order.
    volume, energy : array-like
        Selected observations used for initial estimates.
    constraints : sequence of ParameterConstraint, optional
        User overrides in public physical units.
    energy_unit, pressure_unit, volume_unit : str, optional
        Public parameter units.

    Returns
    -------
    ParameterMap
        Mapping with reporting order ``E0, K0, KP, KPP, V0``.

    Raises
    ------
    ValueError
        If the supplied data or workflow state violates the documented contract.
    """
    eos_model = parse_eos_model(model)
    if not eos_model.supports_energy:
        raise ValueError(f"{eos_model.tag} has no integrated E-V form")
    estimates = estimate_energy_parameters(
        eos_model,
        volume,
        energy,
        energy_unit=energy_unit,
        pressure_unit=pressure_unit,
    )
    overrides = _constraint_overrides(constraints)
    definitions = tuple(
        _energy_parameter_definition(
            eos_model,
            name,
            estimates,
            overrides.get(name),
            energy_unit=energy_unit,
            pressure_unit=pressure_unit,
            volume_unit=volume_unit,
        )
        for name in _ENERGY_PARAMETER_ORDER
    )
    return ParameterMap(definitions, resolver=_energy_resolver(eos_model))


def estimate_energy_parameters(
    model: EOSModel | str,
    volume: np.ndarray | Sequence[float],
    energy: np.ndarray | Sequence[float],
    *,
    energy_unit: str = "Ha",
    pressure_unit: str = "GPa",
) -> dict[str, float]:
    """Estimate complete public E-V parameters from sampled data.

    Parameters
    ----------
    model : EOSModel | str
        Integrated energy-EOS family and order.
    volume : np.ndarray | Sequence[float]
        Positive sampled volumes in angstrom cubed.
    energy : np.ndarray | Sequence[float]
        Sampled electronic energies in ``energy_unit`` aligned with ``volume``.
    energy_unit : str
        Energy unit used for public energy values.
    pressure_unit : str
        Pressure unit used for public pressure values.

    Returns
    -------
    dict[str, float]
        Complete ``E0, K0, KP, KPP, V0`` estimate in public workflow units.

    Raises
    ------
    RuntimeError
        If no stable initial estimate can be constructed from the sampled data.
    """
    eos_model = parse_eos_model(model)
    volume_values, energy_values = _validate_energy_data(volume, energy)
    core = EnergyEOS()
    adapter = CoreEnergyEOSFitModel(core, eos_model)
    free = adapter.initial_guess(volume_values, energy_values)
    resolved = resolve_energy_parameters(eos_model, free)
    if resolved.E0 is None:
        raise RuntimeError("resolved E-V parameters do not contain E0")
    factor = _pressure_per_energy_density(energy_unit, pressure_unit)
    return {
        "E0": float(resolved.E0),
        "K0": float(
            energy_to_pressure(
                resolved.K0,
                energy_unit,
                "angstrom",
                pressure_unit,
            )
        ),
        "KP": float(resolved.KP),
        "KPP": float(resolved.KPP / factor),
        "V0": float(resolved.V0),
    }


def _energy_parameter_definition(
    model: EOSModel,
    name: str,
    estimates: Mapping[str, float],
    override: ParameterConstraint | None,
    *,
    energy_unit: str,
    pressure_unit: str,
    volume_unit: str,
) -> ParameterDefinition:
    """Build one fitted, fixed, or EOS-implied E-V parameter."""
    source = "fitted" if name == "E0" else model.parameter_sources[name]
    unit = _parameter_unit(name, energy_unit, pressure_unit, volume_unit)
    if source == "implied":
        return _implied_definition(model, name, override, unit)
    return _fitted_definition(name, estimates[name], override, unit)


def _implied_definition(
    model: EOSModel,
    name: str,
    override: ParameterConstraint | None,
    unit: str,
) -> ParameterDefinition:
    """Build one physical parameter imposed by EOS family/order."""
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
        description=f"energy-EOS parameter implied by {model.tag}",
        metadata={"source": "eos_order"},
    )


def _fitted_definition(
    name: str,
    estimate: float,
    override: ParameterConstraint | None,
    unit: str,
) -> ParameterDefinition:
    """Build one normally fitted E-V parameter with an optional override."""
    default_lower = _POSITIVE_LOWER_BOUND if name in {"K0", "V0"} else -np.inf
    if override is None:
        return ParameterDefinition.free(
            name,
            estimate,
            lower_bound=default_lower,
            unit=unit,
            description=_parameter_description(name),
            metadata={"initial_source": "integrated_eos_estimate"},
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


def _energy_resolver(model: EOSModel):
    """Return a resolver for public parameters implied by one EOS model."""
    sources = model.parameter_sources

    def resolver(values: Mapping[str, float]) -> Mapping[str, float]:
        k0 = float(values["K0"])
        if sources["KP"] == "implied":
            kp_value = implied_kp(model)
            if kp_value is None:
                raise ValueError(f"{model.tag} does not imply KP")
            kp = float(kp_value)
        else:
            kp = float(values["KP"])
        resolved: dict[str, float] = {}
        if sources["KP"] == "implied":
            resolved["KP"] = kp
        if sources["KPP"] == "implied":
            resolved["KPP"] = implied_kpp(model, k0, kp)
        return resolved

    return resolver


def _constraint_overrides(
    constraints: Sequence[ParameterConstraint],
) -> dict[str, ParameterConstraint]:
    """Normalize E-V user constraints by canonical parameter name."""
    overrides: dict[str, ParameterConstraint] = {}
    for constraint in constraints:
        name = constraint.name.upper()
        if name not in _ENERGY_PARAMETER_ORDER:
            raise ValueError(f"unknown energy-EOS parameter constraint: {name}")
        if name in overrides:
            raise ValueError(f"duplicate energy-EOS parameter constraint: {name}")
        overrides[name] = (
            constraint if constraint.name == name else replace(constraint, name=name)
        )
    return overrides


def _parameter_mapping(
    parameters: np.ndarray | Sequence[float] | Mapping[str, float],
) -> dict[str, float]:
    """Return complete public E-V parameters as a validated mapping."""
    if isinstance(parameters, Mapping):
        try:
            values = {name: float(parameters[name]) for name in _ENERGY_PARAMETER_ORDER}
        except KeyError as exc:
            raise ValueError(
                f"energy EOS parameters are missing {exc.args[0]!r}"
            ) from exc
    else:
        array = np.asarray(parameters, dtype=np.float64)
        if array.ndim != 1 or array.size != len(_ENERGY_PARAMETER_ORDER):
            raise ValueError("energy EOS model requires E0, K0, KP, KPP, V0")
        values = dict(zip(_ENERGY_PARAMETER_ORDER, map(float, array), strict=True))
    if not np.all(np.isfinite(list(values.values()))):
        raise ValueError("energy EOS parameters must be finite")
    if values["K0"] <= 0.0 or values["V0"] <= 0.0:
        raise ValueError("energy EOS requires positive K0 and V0")
    return values


def _validate_energy_data(
    volume: np.ndarray | Sequence[float],
    energy: np.ndarray | Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return finite one-dimensional E-V arrays with positive volume."""
    x = _validate_volume(volume)
    y = np.asarray(energy, dtype=np.float64)
    if y.ndim != 1 or y.size != x.size or y.size == 0:
        raise ValueError("energy-volume data must be equal-length non-empty vectors")
    if not np.all(np.isfinite(y)):
        raise ValueError("energy values must be finite")
    return x, y


def _validate_volume(volume: np.ndarray | Sequence[float]) -> np.ndarray:
    """Return a finite positive one-dimensional volume vector."""
    values = np.asarray(volume, dtype=np.float64)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("volume must be a non-empty finite vector")
    if np.any(values <= 0.0):
        raise ValueError("volume must be strictly positive")
    return values


def _pressure_per_energy_density(energy_unit: str, pressure_unit: str) -> float:
    """Return public pressure units represented by one energy-density unit."""
    return float(
        energy_to_pressure(1.0, energy_unit, "angstrom", pressure_unit)
    )


def _parameter_unit(
    name: str,
    energy_unit: str,
    pressure_unit: str,
    volume_unit: str,
) -> str:
    """Return the public unit of one E-V parameter."""
    return {
        "E0": str(energy_unit),
        "K0": str(pressure_unit),
        "KP": "1",
        "KPP": f"{pressure_unit}^-1",
        "V0": str(volume_unit),
    }[name]


def _parameter_description(name: str) -> str:
    """Return a technical description for one E-V parameter."""
    return {
        "E0": "energy at the zero-pressure reference volume",
        "K0": "reference isothermal bulk modulus",
        "KP": "first pressure derivative of the bulk modulus",
        "KPP": "second pressure derivative of the bulk modulus",
        "V0": "zero-pressure reference volume",
    }[name]


__all__ = [
    "EnergyEOSFitModel",
    "build_energy_parameter_map",
    "estimate_energy_parameters",
]
