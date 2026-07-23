# -*- coding: utf-8 -*-

r"""Thermal-pressure model and normalization contracts.

This module defines the passive scientific contracts and the float64
numerical evaluator used by pressure-volume-temperature equations that supply
their own thermal-pressure term.  The Holland--Powell Einstein evaluator
remains in :mod:`quantas.core.physics.eos.pvt`; the Mie--Grüneisen--Debye
(MGD) implementation here is independent of fitting, persistence, rendering,
and frontend objects.

Two MGD variants are represented:

``full``
    The thermal Grüneisen parameter follows

    .. math::

        \gamma(V)=\gamma_0\left(\frac{V}{V_0}\right)^q,

    and the Debye temperature follows from
    :math:`\gamma=-\partial\ln\theta_D/\partial\ln V`.

``q-compromise``
    The Debye temperature and :math:`\gamma/V` are both held constant.  This
    is the explicitly named EosFit approximation; it is not represented by a
    hidden numerical value of ``q`` because its two assumptions correspond to
    incompatible values of ``q`` in the full MGD law.

MGD thermal pressure may be normalized either to a crystallographic cell or to
one mole of formula units.  The native Quantas convention is a cell volume in
Angstrom cubed together with the number of atoms in that cell.  The molar
formula-unit convention is retained as an equivalent external representation.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
import math
from typing import TypeAlias
import unicodedata

import numpy as np

from quantas.core.chemistry import parse_formula
from quantas.core.physics.units import N, kB

ArrayLike: TypeAlias = np.ndarray | float | Sequence[float]

_CELL_PRESSURE_FACTOR = 3.0 * kB * 1.0e21
_MOLAR_PRESSURE_FACTOR = 3.0 * N * kB * 1.0e-3
_DEBYE_SERIES_LIMIT = 2.0
_DEBYE_ASYMPTOTIC_LIMIT = 50.0
_DEBYE_INFINITY_INTEGRAL = math.pi**4 / 15.0
_DEBYE_SERIES_COEFFICIENTS: tuple[tuple[int, float], ...] = (
    (0, 1.0),
    (1, -3.0 / 8.0),
    (2, 1.0 / 20.0),
    (4, -1.0 / 1680.0),
    (6, 1.0 / 90720.0),
    (8, -1.0 / 4435200.0),
    (10, 4.817713151046484e-09),
    (12, -1.056838027737503e-10),
    (14, 2.361624093650241e-12),
    (16, -5.352126783667251e-14),
    (18, 1.226580293753982e-15),
    (20, -2.836785258988787e-17),
    (22, 6.610803394032284e-19),
    (24, -1.5504960762013934e-20),
    (26, 3.656593489271868e-22),
    (28, -8.664694284229898e-24),
    (30, 2.0617749566706239e-25),
    (32, -4.924106287604752e-27),
    (34, 1.1798695748228652e-28),
    (36, -2.8353807235887045e-30),
    (38, 6.831756773484189e-32),
    (40, -1.6500156388609074e-33),
)


class ThermalPressureFamily(str, Enum):
    """Named thermal-pressure oscillator families.

    ``HOLLAND_POWELL_EINSTEIN`` is the currently executable Quantas
    thermal-pressure formulation. ``MIE_GRUNEISEN_DEBYE`` identifies the MGD
    contracts verified during the EOS thermal-pressure review.
    """

    HOLLAND_POWELL_EINSTEIN = "holland-powell-einstein"
    MIE_GRUNEISEN_DEBYE = "mie-gruneisen-debye"


class MGDVariant(str, Enum):
    """Supported Mie--Grüneisen--Debye model variants."""

    FULL = "full"
    Q_COMPROMISE = "q-compromise"


class MGDVolumeBasis(str, Enum):
    """Volume normalization used by an MGD parameter set."""

    CELL = "cell"
    MOLAR_FORMULA_UNIT = "molar-formula-unit"


@dataclass(frozen=True, slots=True)
class ThermalPressureModel:
    """Describe one thermal-pressure oscillator model.

    Parameters
    ----------
    family : ThermalPressureFamily or str
        Oscillator family.
    variant : MGDVariant, str, or None, optional
        MGD variant. It is required for MGD and must be omitted for the
        Holland--Powell Einstein family. The shorthand ``"mgd"`` resolves to
        the full MGD variant.

    Raises
    ------
    ValueError
        If the family and variant are inconsistent.
    """

    family: ThermalPressureFamily | str
    variant: MGDVariant | str | None = None

    def __post_init__(self) -> None:
        """Normalize the family and validate its variant."""
        family = _parse_thermal_pressure_family(self.family)
        if family is ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN:
            if self.variant is not None:
                raise ValueError(
                    "Holland-Powell Einstein thermal pressure has no variant"
                )
            variant = None
        else:
            variant = (
                MGDVariant.FULL
                if self.variant is None
                else _parse_mgd_variant(self.variant)
            )
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "variant", variant)

    @property
    def family_name(self) -> ThermalPressureFamily:
        """Return the canonical oscillator family."""
        value = self.family
        assert isinstance(value, ThermalPressureFamily)
        return value

    @property
    def mgd_variant(self) -> MGDVariant | None:
        """Return the canonical MGD variant, when applicable."""
        value = self.variant
        assert value is None or isinstance(value, MGDVariant)
        return value

    @property
    def tag(self) -> str:
        """Return a stable serialization tag."""
        variant = self.mgd_variant
        if variant is None:
            return self.family_name.value
        return f"{self.family_name.value}:{variant.value}"

    @property
    def requires_mgd_normalization(self) -> bool:
        """Return whether atom-count and volume-basis metadata are required."""
        return self.family_name is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE

    def as_dict(self) -> dict[str, str | None]:
        """Return a serialization-ready model description."""
        variant = self.mgd_variant
        return {
            "family": self.family_name.value,
            "variant": None if variant is None else variant.value,
            "tag": self.tag,
        }


@dataclass(frozen=True, slots=True)
class MGDNormalization:
    """Describe the atom-count normalization of MGD thermal pressure.

    Parameters
    ----------
    volume_basis : MGDVolumeBasis or str
        ``"cell"`` means that volumes are crystallographic cell volumes in
        Angstrom cubed. ``"molar-formula-unit"`` means molar volumes in cubic
        centimetres per mole of formula units.
    atoms_per_unit : float
        Number of atoms in the volume unit. This is atoms per cell for cell
        volumes and atoms per formula unit for molar formula-unit volumes.
    formula : str or None, optional
        Optional chemical formula retained for provenance.
    formula_units_per_cell : float or None, optional
        Number of formula units in the cell. It is meaningful only for the
        cell basis and is retained for validation and provenance.

    Raises
    ------
    ValueError
        If counts are non-positive, inconsistent with the formula, or invalid
        for the selected volume basis.
    """

    volume_basis: MGDVolumeBasis | str
    atoms_per_unit: float
    formula: str | None = None
    formula_units_per_cell: float | None = None

    def __post_init__(self) -> None:
        """Normalize and validate the volume basis and atom counts."""
        basis = _parse_volume_basis(self.volume_basis)
        atoms = _positive_finite(self.atoms_per_unit, "atoms_per_unit")
        formula = None if self.formula is None else self.formula.strip()
        if formula == "":
            raise ValueError("formula must be non-empty when supplied")
        z = self.formula_units_per_cell
        if z is not None:
            z = _positive_finite(z, "formula_units_per_cell")
        if basis is MGDVolumeBasis.MOLAR_FORMULA_UNIT and z is not None:
            raise ValueError(
                "formula_units_per_cell is only valid for cell-volume normalization"
            )
        if basis is MGDVolumeBasis.CELL and formula is not None and z is None:
            raise ValueError("cell formula metadata requires formula_units_per_cell")
        if basis is MGDVolumeBasis.CELL and z is not None and formula is None:
            raise ValueError("formula_units_per_cell requires a formula")
        if formula is not None:
            atoms_formula = _atoms_in_formula(formula)
            expected = atoms_formula * z if z is not None else atoms_formula
            if not math.isclose(atoms, expected, rel_tol=1.0e-12, abs_tol=1.0e-12):
                label = "cell" if basis is MGDVolumeBasis.CELL else "formula unit"
                raise ValueError(
                    "atoms_per_unit is inconsistent with formula metadata "
                    f"for the {label}"
                )
        object.__setattr__(self, "volume_basis", basis)
        object.__setattr__(self, "atoms_per_unit", atoms)
        object.__setattr__(self, "formula", formula)
        object.__setattr__(self, "formula_units_per_cell", z)

    @classmethod
    def cell(
        cls,
        *,
        atoms_per_cell: float | None = None,
        formula: str | None = None,
        formula_units_per_cell: float | None = None,
    ) -> MGDNormalization:
        """Build cell-volume normalization from atoms or ``formula + Z``.

        Parameters
        ----------
        atoms_per_cell : float or None, optional
            Explicit number of atoms in the crystallographic cell.
        formula : str or None, optional
            Formula of one formula unit.
        formula_units_per_cell : float or None, optional
            Number ``Z`` of formula units in the cell.

        Returns
        -------
        MGDNormalization
            Validated cell-volume normalization.

        Raises
        ------
        ValueError
            If neither an explicit atom count nor a complete ``formula + Z``
            specification is supplied, or if redundant values disagree.
        """
        inferred: float | None = None
        if formula is not None and formula_units_per_cell is not None:
            inferred = _atoms_in_formula(formula) * _positive_finite(
                formula_units_per_cell, "formula_units_per_cell"
            )
        if atoms_per_cell is None:
            if inferred is None:
                raise ValueError(
                    "cell MGD normalization requires atoms_per_cell or formula + "
                    "formula_units_per_cell"
                )
            atoms_per_cell = inferred
        return cls(
            volume_basis=MGDVolumeBasis.CELL,
            atoms_per_unit=atoms_per_cell,
            formula=formula,
            formula_units_per_cell=formula_units_per_cell,
        )

    @classmethod
    def molar_formula_unit(
        cls,
        *,
        atoms_per_formula_unit: float | None = None,
        formula: str | None = None,
    ) -> MGDNormalization:
        """Build molar formula-unit normalization.

        Parameters
        ----------
        atoms_per_formula_unit : float or None, optional
            Explicit atom count in one formula unit.
        formula : str or None, optional
            Formula from which the atom count may be inferred.

        Returns
        -------
        MGDNormalization
            Validated molar formula-unit normalization.

        Raises
        ------
        ValueError
            If neither a formula nor an explicit atom count is supplied, or if
            redundant values disagree.
        """
        if atoms_per_formula_unit is None:
            if formula is None:
                raise ValueError(
                    "molar MGD normalization requires atoms_per_formula_unit or formula"
                )
            atoms_per_formula_unit = _atoms_in_formula(formula)
        return cls(
            volume_basis=MGDVolumeBasis.MOLAR_FORMULA_UNIT,
            atoms_per_unit=atoms_per_formula_unit,
            formula=formula,
        )

    @property
    def basis(self) -> MGDVolumeBasis:
        """Return the canonical volume basis."""
        value = self.volume_basis
        assert isinstance(value, MGDVolumeBasis)
        return value

    @property
    def canonical_volume_unit(self) -> str:
        """Return the canonical volume unit required by this basis."""
        if self.basis is MGDVolumeBasis.CELL:
            return "angstrom^3"
        return "cm^3 mol^-1"

    @property
    def atoms_per_cell(self) -> float | None:
        """Return atoms per cell when using cell-volume normalization."""
        if self.basis is MGDVolumeBasis.CELL:
            return self.atoms_per_unit
        return None

    @property
    def atoms_per_formula_unit(self) -> float | None:
        """Return atoms per formula unit for molar normalization."""
        if self.basis is MGDVolumeBasis.MOLAR_FORMULA_UNIT:
            return self.atoms_per_unit
        if self.formula is not None:
            return _atoms_in_formula(self.formula)
        return None

    def as_dict(self) -> dict[str, str | float | None]:
        """Return a serialization-ready normalization description."""
        return {
            "volume_basis": self.basis.value,
            "volume_unit": self.canonical_volume_unit,
            "atoms_per_unit": self.atoms_per_unit,
            "formula": self.formula,
            "formula_units_per_cell": self.formula_units_per_cell,
        }


@dataclass(frozen=True, slots=True)
class MGDParameters:
    """Reference parameters of an MGD thermal-pressure model.

    Parameters
    ----------
    temperature_ref : float
        Reference-isotherm temperature in kelvin. Zero kelvin is allowed.
    theta_d0 : float
        Debye temperature at the reference volume, in kelvin.
    gamma0 : float
        Thermal Grüneisen parameter at the reference volume. Negative values
        are permitted for materials with negative average mode-Grüneisen
        response.
    q : float or None, optional
        Logarithmic volume exponent of ``gamma`` for the full MGD model. It
        must be omitted for the q-compromise model.

    Raises
    ------
    ValueError
        If a value is non-finite or outside its physical domain.
    """

    temperature_ref: float
    theta_d0: float
    gamma0: float
    q: float | None = None

    def __post_init__(self) -> None:
        """Validate reference-temperature and MGD parameter domains."""
        temperature_ref = _nonnegative_finite(self.temperature_ref, "temperature_ref")
        theta_d0 = _positive_finite(self.theta_d0, "theta_d0")
        gamma0 = _finite(self.gamma0, "gamma0")
        q = None if self.q is None else _finite(self.q, "q")
        object.__setattr__(self, "temperature_ref", temperature_ref)
        object.__setattr__(self, "theta_d0", theta_d0)
        object.__setattr__(self, "gamma0", gamma0)
        object.__setattr__(self, "q", q)

    def validate_for(self, model: ThermalPressureModel | str) -> None:
        """Validate parameter presence against an MGD model contract.

        Parameters
        ----------
        model : ThermalPressureModel or str
            MGD model specification.

        Raises
        ------
        ValueError
            If the model is not MGD, if ``q`` is absent for the full model, or
            if ``q`` is supplied for the q-compromise model.
        """
        specification = parse_thermal_pressure_model(model)
        if specification.family_name is not ThermalPressureFamily.MIE_GRUNEISEN_DEBYE:
            raise ValueError("MGD parameters require a Mie-Gruneisen-Debye model")
        variant = specification.mgd_variant
        assert isinstance(variant, MGDVariant)
        if variant is MGDVariant.FULL and self.q is None:
            raise ValueError("full MGD parameters require q")
        if variant is MGDVariant.Q_COMPROMISE and self.q is not None:
            raise ValueError("q-compromise MGD parameters must omit q")

    def as_dict(self) -> dict[str, float | None]:
        """Return a serialization-ready parameter mapping."""
        return {
            "temperature_ref": self.temperature_ref,
            "theta_d0": self.theta_d0,
            "gamma0": self.gamma0,
            "q": self.q,
        }


def debye_function_3(argument: ArrayLike) -> np.ndarray:
    r"""Evaluate the third-order Debye function in float64.

    The adopted convention is

    .. math::

        D_3(x)=\frac{3}{x^3}\int_0^x\frac{t^3}{\exp(t)-1}\,\mathrm dt,
        \qquad D_3(0)=1.

    A Bernoulli-series branch is used for small and moderate arguments.  For
    larger arguments the integral is evaluated from its exact infinite limit
    minus an exponentially convergent tail.  This avoids numerical quadrature
    in repeated P--V--T calculations.

    Parameters
    ----------
    argument : array-like
        Non-negative dimensionless Debye arguments.

    Returns
    -------
    ndarray
        Debye-function values with the input shape.

    Raises
    ------
    ValueError
        If any argument is non-finite or negative.
    """
    values = np.asarray(argument, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError("Debye arguments must be finite")
    if np.any(values < 0.0):
        raise ValueError("Debye arguments cannot be negative")
    result = np.empty_like(values, dtype=np.float64)
    flat_values = values.reshape(-1)
    flat_result = result.reshape(-1)
    for index, item in enumerate(flat_values):
        flat_result[index] = _debye_function_3_scalar(float(item))
    return result


class MGDThermalPressure:
    r"""Evaluate Mie--Grüneisen--Debye thermal-pressure properties.

    The evaluator implements both the full volume-dependent MGD formulation
    and the explicitly named EosFit ``q-compromise`` approximation.  It owns no
    workflow state and emits no Quantas events.
    """

    def gamma(
        self,
        model: str | ThermalPressureModel,
        parameters: MGDParameters,
        reference_volume: float,
        volume: ArrayLike,
    ) -> np.ndarray:
        """Return the thermal Grüneisen parameter at one or more volumes.

        Parameters
        ----------
        model : str or ThermalPressureModel
            Full or q-compromise MGD model.
        parameters : MGDParameters
            Reference MGD parameters.
        reference_volume : float
            Reference cell or molar volume in the normalization basis.
        volume : array-like
            Positive volumes in the same unit as ``reference_volume``.

        Returns
        -------
        ndarray
            Dimensionless Grüneisen parameters.
        """
        specification, reference, values = self._validated_volume_inputs(
            model, parameters, reference_volume, volume
        )
        ratio = values / reference
        if specification.mgd_variant is MGDVariant.Q_COMPROMISE:
            result = parameters.gamma0 * ratio
        else:
            assert parameters.q is not None
            with np.errstate(over="ignore", invalid="ignore"):
                result = parameters.gamma0 * np.exp(parameters.q * np.log(ratio))
        if not np.all(np.isfinite(result)):
            raise ValueError("MGD model predicts a non-finite gamma(V)")
        return np.asarray(result, dtype=np.float64)

    def debye_temperature(
        self,
        model: str | ThermalPressureModel,
        parameters: MGDParameters,
        reference_volume: float,
        volume: ArrayLike,
    ) -> np.ndarray:
        """Return the volume-dependent Debye temperature in kelvin.

        Parameters
        ----------
        model : str or ThermalPressureModel
            Full or q-compromise MGD model.
        parameters : MGDParameters
            Reference MGD parameters.
        reference_volume : float
            Reference cell or molar volume.
        volume : array-like
            Positive volumes in the same unit as ``reference_volume``.

        Returns
        -------
        ndarray
            Positive Debye temperatures in kelvin.
        """
        specification, reference, values = self._validated_volume_inputs(
            model, parameters, reference_volume, volume
        )
        if specification.mgd_variant is MGDVariant.Q_COMPROMISE:
            return np.full_like(values, parameters.theta_d0, dtype=np.float64)
        assert parameters.q is not None
        logarithm = np.log(values / reference)
        exponent = self._full_theta_exponent(logarithm, parameters.gamma0, parameters.q)
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            result = parameters.theta_d0 * np.exp(exponent)
        if np.any(result <= 0.0) or not np.all(np.isfinite(result)):
            raise ValueError(
                "MGD model predicts a non-positive or non-finite theta_D(V)"
            )
        return np.asarray(result, dtype=np.float64)

    def pressure(
        self,
        model: str | ThermalPressureModel,
        parameters: MGDParameters,
        normalization: MGDNormalization,
        reference_volume: float,
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        r"""Return MGD thermal pressure relative to the reference isotherm.

        Parameters
        ----------
        model : str or ThermalPressureModel
            Full or q-compromise MGD model.
        parameters : MGDParameters
            Reference MGD parameters.
        normalization : MGDNormalization
            Atom count and volume basis.
        reference_volume : float
            Reference volume in Angstrom cubed for cell normalization or
            cubic centimetres per mole for molar normalization.
        volume, temperature : array-like
            Broadcast-compatible positive volumes and non-negative
            temperatures.

        Returns
        -------
        ndarray
            Thermal pressure in GPa.
        """
        specification, reference, volumes, temperatures = self._validated_state(
            model, parameters, reference_volume, volume, temperature
        )
        gamma = self.gamma(specification, parameters, reference, volumes)
        theta = self.debye_temperature(specification, parameters, reference, volumes)
        thermal_difference = _debye_thermal_term(temperatures, theta) - (
            _debye_thermal_term(parameters.temperature_ref, theta)
        )
        if specification.mgd_variant is MGDVariant.Q_COMPROMISE:
            gamma_over_volume = np.full_like(
                volumes, parameters.gamma0 / reference, dtype=np.float64
            )
        else:
            gamma_over_volume = gamma / volumes
        result = (
            self._normalization_factor(normalization)
            * gamma_over_volume
            * thermal_difference
        )
        result = np.asarray(result, dtype=np.float64)
        result[temperatures == parameters.temperature_ref] = 0.0
        if not np.all(np.isfinite(result)):
            raise ValueError("MGD model predicts a non-finite thermal pressure")
        return result

    def temperature_derivative(
        self,
        model: str | ThermalPressureModel,
        parameters: MGDParameters,
        normalization: MGDNormalization,
        reference_volume: float,
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        r"""Return :math:`(\partial P_{th}/\partial T)_V` in GPa K\ :sup:`-1`.

        Parameters are identical to :meth:`pressure`.
        """
        specification, reference, volumes, temperatures = self._validated_state(
            model, parameters, reference_volume, volume, temperature
        )
        gamma = self.gamma(specification, parameters, reference, volumes)
        theta = self.debye_temperature(specification, parameters, reference, volumes)
        if specification.mgd_variant is MGDVariant.Q_COMPROMISE:
            gamma_over_volume = np.full_like(
                volumes, parameters.gamma0 / reference, dtype=np.float64
            )
        else:
            gamma_over_volume = gamma / volumes
        result = (
            self._normalization_factor(normalization)
            * gamma_over_volume
            * _debye_thermal_term_temperature_derivative(temperatures, theta)
        )
        if not np.all(np.isfinite(result)):
            raise ValueError("MGD model predicts a non-finite dP_th/dT")
        return np.asarray(result, dtype=np.float64)

    def volume_derivative(
        self,
        model: str | ThermalPressureModel,
        parameters: MGDParameters,
        normalization: MGDNormalization,
        reference_volume: float,
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        r"""Return :math:`(\partial P_{th}/\partial V)_T`.

        The returned unit is GPa per unit of the selected volume basis:
        GPa Angstrom\ :sup:`-3` for cell volumes and GPa per
        cm\ :sup:`3` mol\ :sup:`-1` for molar volumes.

        Parameters are identical to :meth:`pressure`.
        """
        specification, reference, volumes, temperatures = self._validated_state(
            model, parameters, reference_volume, volume, temperature
        )
        if specification.mgd_variant is MGDVariant.Q_COMPROMISE:
            return np.zeros_like(volumes, dtype=np.float64)
        assert parameters.q is not None
        gamma = self.gamma(specification, parameters, reference, volumes)
        theta = self.debye_temperature(specification, parameters, reference, volumes)
        thermal_difference = _debye_thermal_term(temperatures, theta) - (
            _debye_thermal_term(parameters.temperature_ref, theta)
        )
        theta_difference = _debye_thermal_term_theta_derivative(
            temperatures, theta
        ) - _debye_thermal_term_theta_derivative(parameters.temperature_ref, theta)
        result = (
            self._normalization_factor(normalization)
            * gamma
            / volumes**2
            * (
                (parameters.q - 1.0) * thermal_difference
                - gamma * theta * theta_difference
            )
        )
        result = np.asarray(result, dtype=np.float64)
        result[temperatures == parameters.temperature_ref] = 0.0
        if not np.all(np.isfinite(result)):
            raise ValueError("MGD model predicts a non-finite dP_th/dV")
        return result

    def bulk_modulus_contribution(
        self,
        model: str | ThermalPressureModel,
        parameters: MGDParameters,
        normalization: MGDNormalization,
        reference_volume: float,
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> np.ndarray:
        r"""Return the MGD contribution :math:`-V(\partial P_{th}/\partial V)_T`.

        Parameters are identical to :meth:`pressure`.

        Returns
        -------
        ndarray
            Thermal contribution to the isothermal bulk modulus in GPa.
        """
        volumes, temperatures = np.broadcast_arrays(
            np.asarray(volume, dtype=np.float64),
            np.asarray(temperature, dtype=np.float64),
        )
        derivative = self.volume_derivative(
            model,
            parameters,
            normalization,
            reference_volume,
            volumes,
            temperatures,
        )
        return np.asarray(-volumes * derivative, dtype=np.float64)

    @staticmethod
    def _full_theta_exponent(
        logarithmic_volume_ratio: np.ndarray,
        gamma0: float,
        q: float,
    ) -> np.ndarray:
        """Return a cancellation-safe full-MGD theta exponent."""
        if q == 0.0:
            return -gamma0 * logarithmic_volume_ratio
        scaled = q * logarithmic_volume_ratio
        ratio = np.ones_like(scaled, dtype=np.float64)
        nonzero = scaled != 0.0
        ratio[nonzero] = np.expm1(scaled[nonzero]) / scaled[nonzero]
        return -gamma0 * logarithmic_volume_ratio * ratio

    @staticmethod
    def _normalization_factor(normalization: MGDNormalization) -> float:
        """Return the cell- or molar-basis pressure prefactor."""
        if normalization.basis is MGDVolumeBasis.CELL:
            return _CELL_PRESSURE_FACTOR * normalization.atoms_per_unit
        return _MOLAR_PRESSURE_FACTOR * normalization.atoms_per_unit

    @staticmethod
    def _validated_volume_inputs(
        model: str | ThermalPressureModel,
        parameters: MGDParameters,
        reference_volume: float,
        volume: ArrayLike,
    ) -> tuple[ThermalPressureModel, float, np.ndarray]:
        """Validate an MGD model, parameters, and volume coordinates."""
        specification = parse_thermal_pressure_model(model)
        parameters.validate_for(specification)
        reference = _positive_finite(reference_volume, "reference_volume")
        values = np.asarray(volume, dtype=np.float64)
        if not np.all(np.isfinite(values)):
            raise ValueError("volume values must be finite")
        if np.any(values <= 0.0):
            raise ValueError("volume values must be positive")
        return specification, reference, values.copy()

    @classmethod
    def _validated_state(
        cls,
        model: str | ThermalPressureModel,
        parameters: MGDParameters,
        reference_volume: float,
        volume: ArrayLike,
        temperature: ArrayLike,
    ) -> tuple[ThermalPressureModel, float, np.ndarray, np.ndarray]:
        """Validate and broadcast MGD state coordinates."""
        specification, reference, volumes = cls._validated_volume_inputs(
            model, parameters, reference_volume, volume
        )
        temperatures = np.asarray(temperature, dtype=np.float64)
        volumes, temperatures = np.broadcast_arrays(volumes, temperatures)
        if not np.all(np.isfinite(temperatures)):
            raise ValueError("temperature values must be finite")
        if np.any(temperatures < 0.0):
            raise ValueError("temperature cannot be below absolute zero")
        return (
            specification,
            reference,
            volumes.copy(),
            temperatures.copy(),
        )


def _debye_function_3_scalar(argument: float) -> float:
    """Evaluate one non-negative Debye argument."""
    if argument == 0.0:
        return 1.0
    if argument <= _DEBYE_SERIES_LIMIT:
        return math.fsum(
            coefficient * argument**power
            for power, coefficient in _DEBYE_SERIES_COEFFICIENTS
        )
    if argument >= _DEBYE_ASYMPTOTIC_LIMIT:
        return math.pi**4 / (5.0 * argument**3)
    tail = 0.0
    index = 1
    while True:
        value = float(index)
        exponential = math.exp(-value * argument)
        term = exponential * (
            argument**3 / value
            + 3.0 * argument**2 / value**2
            + 6.0 * argument / value**3
            + 6.0 / value**4
        )
        tail += term
        if term <= np.finfo(np.float64).eps * max(abs(tail), 1.0):
            break
        index += 1
    integral = _DEBYE_INFINITY_INTEGRAL - tail
    return 3.0 * integral / argument**3


def _debye_function_3_derivative(argument: np.ndarray) -> np.ndarray:
    """Return the analytical derivative of the third-order Debye function."""
    values = np.asarray(argument, dtype=np.float64)
    result = np.empty_like(values)
    small = values <= _DEBYE_SERIES_LIMIT
    if np.any(small):
        x = values[small]
        result[small] = np.asarray(
            [
                math.fsum(
                    power * coefficient * float(item) ** (power - 1)
                    for power, coefficient in _DEBYE_SERIES_COEFFICIENTS
                    if power > 0
                )
                for item in x.reshape(-1)
            ],
            dtype=np.float64,
        ).reshape(x.shape)
    if np.any(~small):
        x = values[~small]
        debye = debye_function_3(x)
        occupation = np.zeros_like(x)
        moderate = x < _DEBYE_ASYMPTOTIC_LIMIT
        if np.any(moderate):
            occupation[moderate] = 1.0 / np.expm1(x[moderate])
        result[~small] = 3.0 * occupation - 3.0 * debye / x
    return result


def _debye_thermal_term(
    temperature: ArrayLike,
    theta: ArrayLike,
) -> np.ndarray:
    """Return :math:`T D_3(\theta/T)` in kelvin."""
    temperatures, theta_values = np.broadcast_arrays(
        np.asarray(temperature, dtype=np.float64),
        np.asarray(theta, dtype=np.float64),
    )
    if not np.all(np.isfinite(temperatures)) or np.any(temperatures < 0.0):
        raise ValueError("temperature values must be finite and non-negative")
    if not np.all(np.isfinite(theta_values)) or np.any(theta_values <= 0.0):
        raise ValueError("Debye temperatures must be finite and positive")
    result = np.zeros_like(temperatures)
    positive = temperatures > 0.0
    if np.any(positive):
        argument = theta_values[positive] / temperatures[positive]
        result[positive] = temperatures[positive] * debye_function_3(argument)
    return result


def _debye_thermal_term_temperature_derivative(
    temperature: ArrayLike,
    theta: ArrayLike,
) -> np.ndarray:
    """Return the derivative of :math:`T D_3(\theta/T)` with respect to T."""
    temperatures, theta_values = np.broadcast_arrays(
        np.asarray(temperature, dtype=np.float64),
        np.asarray(theta, dtype=np.float64),
    )
    result = np.zeros_like(temperatures)
    positive = temperatures > 0.0
    if np.any(positive):
        argument = theta_values[positive] / temperatures[positive]
        debye = debye_function_3(argument)
        bose_ratio = np.zeros_like(argument)
        moderate = argument < _DEBYE_ASYMPTOTIC_LIMIT
        if np.any(moderate):
            bose_ratio[moderate] = argument[moderate] / np.expm1(argument[moderate])
        result[positive] = 4.0 * debye - 3.0 * bose_ratio
    return result


def _debye_thermal_term_theta_derivative(
    temperature: ArrayLike,
    theta: ArrayLike,
) -> np.ndarray:
    """Return the derivative of :math:`T D_3(\theta/T)` with respect to theta."""
    temperatures, theta_values = np.broadcast_arrays(
        np.asarray(temperature, dtype=np.float64),
        np.asarray(theta, dtype=np.float64),
    )
    result = np.zeros_like(temperatures)
    positive = temperatures > 0.0
    if np.any(positive):
        argument = theta_values[positive] / temperatures[positive]
        result[positive] = _debye_function_3_derivative(argument)
    return result


def parse_thermal_pressure_model(
    model: str | ThermalPressureModel,
) -> ThermalPressureModel:
    """Return a canonical thermal-pressure model specification.

    Parameters
    ----------
    model : str or ThermalPressureModel
        Canonical model tag or documented alias. ``"mgd"`` selects the full
        MGD model. The explicit q-compromise tag is ``"mgd:q-compromise"``.

    Returns
    -------
    ThermalPressureModel
        Canonical passive model specification.

    Raises
    ------
    ValueError
        If the model name is unknown.
    """
    if isinstance(model, ThermalPressureModel):
        return model
    text = str(model).strip().lower()
    normalized = _normalized_token(text)
    if normalized in {
        "einstein",
        "hp",
        "hpeinstein",
        "hollandpowell",
        "hollandpowelleinstein",
    }:
        return ThermalPressureModel(ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN)
    if normalized in {
        "mgd",
        "mgdfull",
        "miegruneisendebye",
        "miegruneisendebyefull",
    }:
        return ThermalPressureModel(
            ThermalPressureFamily.MIE_GRUNEISEN_DEBYE,
            MGDVariant.FULL,
        )
    if normalized in {
        "mgdqcompromise",
        "qcompromise",
        "miegruneisendebyeqcompromise",
    }:
        return ThermalPressureModel(
            ThermalPressureFamily.MIE_GRUNEISEN_DEBYE,
            MGDVariant.Q_COMPROMISE,
        )
    raise ValueError(f"unknown thermal-pressure model: {model!r}")


def thermal_pressure_model_contracts() -> tuple[ThermalPressureModel, ...]:
    """Return all frozen thermal-pressure model contracts.

    Returns
    -------
    tuple of ThermalPressureModel
        Holland--Powell Einstein, full MGD, and q-compromise MGD contracts.

    Notes
    -----
    This registry describes public model contracts. It does not by itself
    assert that every model has an executable evaluator in the current phase.
    """
    return (
        ThermalPressureModel(ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN),
        ThermalPressureModel(
            ThermalPressureFamily.MIE_GRUNEISEN_DEBYE, MGDVariant.FULL
        ),
        ThermalPressureModel(
            ThermalPressureFamily.MIE_GRUNEISEN_DEBYE,
            MGDVariant.Q_COMPROMISE,
        ),
    )


def thermal_pressure_parameter_names(
    model: str | ThermalPressureModel,
) -> tuple[str, ...]:
    """Return the ordered parameters owned by a thermal-pressure model.

    Parameters
    ----------
    model : str or ThermalPressureModel
        Thermal-pressure model contract.

    Returns
    -------
    tuple of str
        Stable parameter names. Reference-isotherm pressure parameters such as
        ``K0`` and ``V0`` remain owned by the pressure EOS and are not repeated.
    """
    specification = parse_thermal_pressure_model(model)
    if specification.family_name is ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN:
        return ("temperature_ref", "alpha_ref", "theta_e")
    if specification.mgd_variant is MGDVariant.Q_COMPROMISE:
        return ("temperature_ref", "theta_d0", "gamma0")
    return ("temperature_ref", "theta_d0", "gamma0", "q")


def _parse_thermal_pressure_family(
    family: str | ThermalPressureFamily,
) -> ThermalPressureFamily:
    if isinstance(family, ThermalPressureFamily):
        return family
    normalized = _normalized_token(str(family))
    if normalized in {"einstein", "hpeinstein", "hollandpowelleinstein"}:
        return ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN
    if normalized in {"mgd", "miegruneisendebye"}:
        return ThermalPressureFamily.MIE_GRUNEISEN_DEBYE
    raise ValueError(f"unknown thermal-pressure family: {family!r}")


def _parse_mgd_variant(variant: str | MGDVariant) -> MGDVariant:
    if isinstance(variant, MGDVariant):
        return variant
    normalized = _normalized_token(str(variant))
    if normalized in {"full", "qrefined", "qfixed"}:
        return MGDVariant.FULL
    if normalized in {"qcompromise", "compromise"}:
        return MGDVariant.Q_COMPROMISE
    raise ValueError(f"unknown MGD variant: {variant!r}")


def _parse_volume_basis(basis: str | MGDVolumeBasis) -> MGDVolumeBasis:
    if isinstance(basis, MGDVolumeBasis):
        return basis
    normalized = _normalized_token(str(basis))
    if normalized in {"cell", "unitcell", "crystallographiccell"}:
        return MGDVolumeBasis.CELL
    if normalized in {
        "molar",
        "molarformulaunit",
        "formulaunitmolar",
        "cm3mol",
    }:
        return MGDVolumeBasis.MOLAR_FORMULA_UNIT
    raise ValueError(f"unknown MGD volume basis: {basis!r}")


def _atoms_in_formula(formula: str) -> float:
    counts = parse_formula(formula)
    atoms = float(sum(float(count) for count in counts.values()))
    return _positive_finite(atoms, "atoms in formula")


def _normalized_token(value: str) -> str:
    normalized = unicodedata.normalize("NFKD", value)
    ascii_text = "".join(
        character for character in normalized if not unicodedata.combining(character)
    )
    return "".join(character for character in ascii_text.lower() if character.isalnum())


def _finite(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_finite(value: float, name: str) -> float:
    result = _finite(value, name)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _nonnegative_finite(value: float, name: str) -> float:
    result = _finite(value, name)
    if result < 0.0:
        raise ValueError(f"{name} cannot be negative")
    return result


__all__ = [
    "MGDNormalization",
    "MGDThermalPressure",
    "MGDParameters",
    "MGDVariant",
    "MGDVolumeBasis",
    "ThermalPressureFamily",
    "ThermalPressureModel",
    "debye_function_3",
    "parse_thermal_pressure_model",
    "thermal_pressure_model_contracts",
    "thermal_pressure_parameter_names",
]
