# -*- coding: utf-8 -*-

"""Thermodynamic helpers for quasi-harmonic approximation workflows.

This module builds harmonic thermodynamic data on the sampled volume grid and
interpolates those quantities to the pressure-temperature equilibrium volumes
obtained by the QHA workflow.  The functions are numerical utilities and do not
emit Quantas events or render text.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import constants as cs

from quantas.core.math.derivative import (
    polynomial_derivative,
    polynomial_derivative_from_coefficients,
    relative_derivative,
)
from quantas.core.math.polynomials import fit_polynomial_result as polynomial_fit
from quantas.core.physics.thermodynamics import (
    entropy,
    free_energy,
    internal_energy,
    isochoric_heat_capacity,
    thermal_energy,
    vibrational_free_energy,
    zero_point_energy,
)
from quantas.core.physics.units import (
    N,
    convert_energy,
    convert_energy_per_temperature,
    convert_frequency,
    convert_pressure,
    convert_temperature,
    convert_volume,
    energy_to_pressure,
)
from quantas.models.thermodynamics import HarmonicThermodynamicResult
from quantas.modules.qha.analysis import pressure_energy_density
from quantas.modules.qha.core.interpolation import (
    fit_polynomial_series,
)
from quantas.modules.qha.core.gruneisen import thermal_gruneisen_from_modes
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult


_PROPERTY_MAP: tuple[tuple[str, str], ...] = (
    ("zero_point_energy", "zero_point_energy"),
    ("thermal_energy", "thermal_energy"),
    ("internal_energy", "internal_energy"),
    ("entropy", "entropy"),
    ("vibrational_free_energy", "vibrational_free_energy"),
    ("free_energy", "free_energy"),
    ("isochoric_heat_capacity", "isochoric_heat_capacity"),
)

THERMAL_EXPANSION_SOURCE_CODES: dict[str, int] = {
    "invalid": 0,
    "mixed_derivative": 1,
    "mode_gruneisen": 2,
    "numerical": 3,
    "numerical_fallback": 4,
}


def calculate_sampled_thermodynamics(
    input_data: QHAInput,
    options: QHAOptions,
) -> HarmonicThermodynamicResult:
    """Calculate harmonic thermodynamic properties at sampled volumes.

    Parameters
    ----------
    input_data : QHAInput
        QHA input containing static energies, phonon frequencies, and q-point
        weights on the sampled volume grid.
    options : QHAOptions
        QHA options defining temperature and unit conventions.

    Returns
    -------
    HarmonicThermodynamicResult
        Harmonic properties sampled on the temperature-volume grid.

    Raises
    ------
    ValueError
        If required input arrays are missing or inconsistent.
    """
    if input_data.volume is None or input_data.energy is None:
        raise ValueError("volume and static energy data are required")
    if input_data.frequencies is None or input_data.weights is None:
        raise ValueError("phonon frequencies and q-point weights are required")

    temperature = np.asarray(options.temperature_grid(), dtype=np.float64)
    temperature_k = np.asarray(
        convert_temperature(temperature, options.temperature_unit, "K"),
        dtype=np.float64,
    )
    frequencies_hz = np.asarray(
        convert_frequency(input_data.frequencies, options.frequency_unit, "Hz"),
        dtype=np.float64,
    )
    weights = input_data.normalized_weights()
    static_energy = np.asarray(input_data.energy, dtype=np.float64)

    uzp = np.asarray(
        convert_energy(
            zero_point_energy(
                np.zeros(1, dtype=np.float64),
                frequencies_hz,
                weights,
            ),
            "kjmol",
            options.energy_unit,
        ),
        dtype=np.float64,
    )
    uth = np.asarray(
        convert_energy(
            thermal_energy(temperature_k, frequencies_hz, weights),
            "kjmol",
            options.energy_unit,
        ),
        dtype=np.float64,
    )
    entropy_values = np.asarray(
        convert_energy_per_temperature(
            entropy(temperature_k, frequencies_hz, weights),
            "J mol^-1 K^-1",
            f"{options.energy_unit} cell^-1 K^-1",
        ),
        dtype=np.float64,
    )
    cv = np.asarray(
        convert_energy_per_temperature(
            isochoric_heat_capacity(temperature_k, frequencies_hz, weights),
            "J mol^-1 K^-1",
            f"{options.energy_unit} cell^-1 K^-1",
        ),
        dtype=np.float64,
    )
    fvib = np.asarray(
        convert_energy(
            vibrational_free_energy(temperature_k, frequencies_hz, weights),
            "kjmol",
            options.energy_unit,
        ),
        dtype=np.float64,
    )
    result = HarmonicThermodynamicResult(
        jobname=input_data.jobname,
        temperature=temperature,
        volume=np.asarray(input_data.volume, dtype=np.float64).copy(),
        static_energy=static_energy.copy(),
        zero_point_energy=uzp,
        thermal_energy=uth,
        internal_energy=internal_energy(static_energy, uzp[0], uth),
        entropy=entropy_values,
        vibrational_free_energy=fvib,
        free_energy=free_energy(static_energy, fvib),
        isochoric_heat_capacity=cv,
        metadata={
            "module": "qha",
            "method": "sampled_harmonic_thermodynamics",
            "backend": {"thermodynamics": "numpy"},
            "units": {
                "energy": options.energy_unit,
                "entropy": f"{options.energy_unit} cell^-1 K^-1",
                "heat_capacity": f"{options.energy_unit} cell^-1 K^-1",
                "temperature": options.temperature_unit,
                "frequency": options.frequency_unit,
            },
        },
    )
    return result


def free_energy_grid(
    harmonic_result: HarmonicThermodynamicResult,
) -> np.ndarray:
    """Return Helmholtz free energy as a temperature-volume array.

    Parameters
    ----------
    harmonic_result : HarmonicThermodynamicResult
        Harmonic thermodynamic result.

    Returns
    -------
    ndarray
        Two-dimensional Helmholtz free-energy grid.

    Raises
    ------
    ValueError
        If free-energy data are missing or have an unsupported shape.
    """
    if harmonic_result.free_energy is None:
        raise ValueError("Helmholtz free energy is not available")
    array = np.asarray(harmonic_result.free_energy, dtype=np.float64)
    if array.ndim == 1:
        if harmonic_result.temperature is None:
            raise ValueError("temperature grid is required for free energy")
        return np.broadcast_to(
            array,
            (np.asarray(harmonic_result.temperature).size, array.size),
        ).copy()
    if array.ndim == 2:
        return array.copy()
    raise ValueError("free energy must be one- or two-dimensional")


class FrequencyThermodynamicEvaluator:
    """Evaluate harmonic thermodynamics from volume-fitted phonon modes.

    Parameters
    ----------
    input_data : QHAInput
        QHA input containing static energies, mode-resolved frequencies and
        q-point weights.
    options : QHAOptions
        QHA options defining polynomial degrees and unit conventions.

    Raises
    ------
    ValueError
        If required input data are missing or a polynomial fit fails.
    """

    def __init__(self, input_data: QHAInput, options: QHAOptions) -> None:
        """Fit frequency and static-energy models over sampled volumes.

        Parameters
        ----------
        input_data : QHAInput
            Volume-dependent phonon and static-energy input data.
        options : QHAOptions
            Polynomial degrees and unit conventions used by the evaluator.

        Raises
        ------
        ValueError
            If required data are missing or a polynomial fit fails.
        """
        if input_data.frequencies is None or input_data.weights is None:
            raise ValueError(
                "frequency evaluation requires phonons and q-point weights"
            )
        if input_data.volume is None or input_data.energy is None:
            raise ValueError(
                "frequency evaluation requires volume and static energy data"
            )

        self.input_data = input_data
        self.options = options
        self.sampled_volume = np.asarray(input_data.volume, dtype=np.float64)
        frequencies = np.asarray(input_data.frequencies, dtype=np.float64)

        self.frequency_degree = min(
            int(options.frequency_degree),
            max(self.sampled_volume.size - 1, 1),
        )
        self.frequency_fit = fit_polynomial_series(
            self.sampled_volume,
            frequencies,
            self.frequency_degree,
            axis=-1,
            metadata={"quantity": "frequency"},
        )
        if not self.frequency_fit.success or self.frequency_fit.coefficients is None:
            detail = "; ".join(self.frequency_fit.warnings) or "frequency fit failed"
            raise ValueError(detail)

        self.energy_degree = min(
            int(options.energy_degree),
            max(self.sampled_volume.size - 1, 1),
        )
        self.energy_fit = polynomial_fit(
            self.sampled_volume,
            np.asarray(input_data.energy, dtype=np.float64),
            self.energy_degree,
        )
        if not self.energy_fit.success or self.energy_fit.parameters is None:
            raise ValueError(self.energy_fit.message or "static-energy fit failed")

        self.backend_name = "numpy"
        self.weights = input_data.normalized_weights()
        self._frequency_coefficients = np.moveaxis(
            self.frequency_fit.coefficients,
            -1,
            0,
        )

    def frequencies_at(self, volume: np.ndarray | float) -> np.ndarray:
        """Return mode-resolved frequencies at selected volumes.

        Parameters
        ----------
        volume : ndarray or float
            Volumes at which frequencies are evaluated.

        Returns
        -------
        ndarray
            Interpolated frequencies in the input frequency unit.

        Raises
        ------
        ValueError
            If interpolation produces non-finite values.
        """
        values = np.asarray(
            np.polynomial.polynomial.polyval(
                np.asarray(volume, dtype=np.float64),
                self._frequency_coefficients,
            ),
            dtype=np.float64,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("frequency interpolation produced non-finite values")
        return values

    def mode_gruneisen_at(
        self,
        volume: np.ndarray | float,
    ) -> np.ndarray:
        """Return mode Grüneisen parameters from the fitted frequencies.

        The derivative is evaluated from the same polynomial frequency models
        used by the frequency QHA workflow,

        ``gamma = -(V / frequency) * d(frequency) / dV``.

        Parameters
        ----------
        volume : ndarray or float
            Volumes at which mode Grüneisen parameters are evaluated.

        Returns
        -------
        ndarray
            Mode Grüneisen parameters. Non-positive interpolated frequencies
            are represented by ``NaN``.

        Raises
        ------
        ValueError
            If requested volumes are non-positive or non-finite.
        """
        volume_array = np.asarray(volume, dtype=np.float64)
        if not np.all(np.isfinite(volume_array)) or np.any(volume_array <= 0.0):
            raise ValueError("volumes must be finite and positive")
        frequency = self.frequencies_at(volume_array)
        derivative = polynomial_derivative_from_coefficients(
            self._frequency_coefficients,
            volume_array,
            axis=0,
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            gamma = -(volume_array * derivative) / frequency
        return np.where(
            np.isfinite(gamma) & (frequency > 0.0),
            gamma,
            np.nan,
        )

    def static_energy_at(self, volume: np.ndarray | float) -> np.ndarray:
        """Return fitted static energies at selected volumes.

        Parameters
        ----------
        volume : ndarray or float
            Volumes at which static energy is evaluated.

        Returns
        -------
        ndarray
            Static energies in the selected QHA energy unit.
        """
        return np.asarray(
            np.polynomial.polynomial.polyval(
                np.asarray(volume, dtype=np.float64),
                self.energy_fit.parameters,
            ),
            dtype=np.float64,
        )

    def properties_at(
        self,
        volume: np.ndarray | float,
        temperature: float,
    ) -> dict[str, np.ndarray]:
        """Return harmonic properties at selected volume and temperature.

        Parameters
        ----------
        volume : ndarray or float
            Volumes at which the properties are evaluated.
        temperature : float
            Temperature in the QHA temperature unit.

        Returns
        -------
        dict
            Harmonic energy, entropy and heat-capacity arrays.
        """
        volume_array = np.atleast_1d(np.asarray(volume, dtype=np.float64))
        temperature_k = float(
            convert_temperature(
                temperature,
                self.options.temperature_unit,
                "K",
            )
        )
        interpolated_frequencies = self.frequencies_at(volume_array)
        frequency_hz = np.asarray(
            convert_frequency(
                interpolated_frequencies,
                self.options.frequency_unit,
                "Hz",
            ),
            dtype=np.float64,
        )
        local_temperature = np.asarray([temperature_k], dtype=np.float64)
        static_energy = self.static_energy_at(volume_array)

        uzp = np.asarray(
            convert_energy(
                zero_point_energy(local_temperature, frequency_hz, self.weights)[0],
                "kjmol",
                self.options.energy_unit,
            ),
            dtype=np.float64,
        )
        uth = np.asarray(
            convert_energy(
                thermal_energy(local_temperature, frequency_hz, self.weights)[0],
                "kjmol",
                self.options.energy_unit,
            ),
            dtype=np.float64,
        )
        entropy_values = np.asarray(
            convert_energy_per_temperature(
                entropy(local_temperature, frequency_hz, self.weights)[0],
                "J mol^-1 K^-1",
                f"{self.options.energy_unit} cell^-1 K^-1",
            ),
            dtype=np.float64,
        )
        cv = np.asarray(
            convert_energy_per_temperature(
                isochoric_heat_capacity(local_temperature, frequency_hz, self.weights)[
                    0
                ],
                "J mol^-1 K^-1",
                f"{self.options.energy_unit} cell^-1 K^-1",
            ),
            dtype=np.float64,
        )
        fvib = np.asarray(
            convert_energy(
                vibrational_free_energy(local_temperature, frequency_hz, self.weights)[
                    0
                ],
                "kjmol",
                self.options.energy_unit,
            ),
            dtype=np.float64,
        )
        return {
            "zero_point_energy": uzp,
            "thermal_energy": uth,
            "internal_energy": static_energy + uzp + uth,
            "entropy": entropy_values,
            "vibrational_free_energy": fvib,
            "free_energy": static_energy + fvib,
            "isochoric_heat_capacity": cv,
        }

    def free_energy_at(
        self,
        volume: np.ndarray | float,
        temperature: float,
    ) -> np.ndarray:
        """Return Helmholtz free energy at selected states.

        Parameters
        ----------
        volume : ndarray or float
            Volumes at which free energy is evaluated.
        temperature : float
            Temperature in the QHA temperature unit.

        Returns
        -------
        ndarray
            Total Helmholtz free energy.
        """
        return self.properties_at(volume, temperature)["free_energy"]


def calculate_frequency_thermodynamics_at_equilibrium(
    result: QHAResult,
    input_data: QHAInput,
    options: QHAOptions,
    *,
    evaluator: FrequencyThermodynamicEvaluator | None = None,
) -> None:
    """Calculate QHA properties from frequencies evaluated at equilibrium.

    Mode-resolved frequencies are fitted as functions of volume and evaluated
    at each equilibrium volume obtained from the pressure-temperature
    minimization. Harmonic thermodynamic functions are then recalculated from
    those frequencies instead of interpolating precomputed thermodynamic
    properties.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volumes and pressure-temperature
        axes. The object is updated in place.
    input_data : QHAInput
        QHA input containing mode-resolved frequencies, static energies and
        q-point weights.
    options : QHAOptions
        QHA options defining fit degrees and unit conventions.
    evaluator : FrequencyThermodynamicEvaluator or None, optional
        Prepared evaluator. A new one is created when omitted.

    Raises
    ------
    ValueError
        If required phonon data are missing, a polynomial fit fails, or the
        result axes are incompatible.
    """
    if result.equilibrium_volume is None or result.temperature is None:
        return
    prepared = evaluator or FrequencyThermodynamicEvaluator(input_data, options)
    equilibrium_volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    temperature = np.asarray(result.temperature, dtype=np.float64)
    valid_mask = np.ones_like(equilibrium_volume, dtype=bool)
    if result.valid_mask is not None:
        valid_mask = np.asarray(result.valid_mask, dtype=bool)

    shape = equilibrium_volume.shape
    properties = {
        name: np.full(shape, np.nan, dtype=np.float64)
        for name in (
            "zero_point_energy",
            "thermal_energy",
            "internal_energy",
            "entropy",
            "vibrational_free_energy",
            "free_energy",
            "isochoric_heat_capacity",
        )
    }

    for it, temperature_value in enumerate(temperature):
        mask = valid_mask[it]
        if not np.any(mask):
            continue
        property_values = prepared.properties_at(
            equilibrium_volume[it, mask],
            float(temperature_value),
        )
        for name in properties:
            properties[name][it, mask] = property_values[name]

    for name, property_grid in properties.items():
        setattr(result, name, property_grid)

    _finalize_equilibrium_thermodynamics(result, options)
    result.metadata.setdefault("thermodynamics", {})
    result.metadata["thermodynamics"].update(
        {
            "source": "mode-resolved frequencies evaluated at equilibrium volume",
            "scheme": "freq",
            "frequency_degree": prepared.frequency_degree,
            "static_energy_degree": prepared.energy_degree,
            "frequency_fit_status": prepared.frequency_fit.status.value,
            "frequency_fit_failures": len(prepared.frequency_fit.failed_indices),
            "backend": prepared.backend_name,
        }
    )


def interpolate_thermodynamics_at_equilibrium(
    result: QHAResult,
    harmonic_result: HarmonicThermodynamicResult,
    options: QHAOptions,
) -> None:
    """Interpolate sampled thermodynamic properties to equilibrium volumes.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volumes and pressure-temperature axes.
        The object is updated in place.
    harmonic_result : HarmonicThermodynamicResult
        Harmonic thermodynamic properties sampled on the input volume grid.
    options : QHAOptions
        QHA options defining the interpolation polynomial degree and units.

    Raises
    ------
    ValueError
        If required volume or temperature axes are missing or incompatible.
    """
    if (
        result.equilibrium_volume is None
        or result.temperature is None
        or result.pressure is None
    ):
        return
    if harmonic_result.volume is None or harmonic_result.temperature is None:
        raise ValueError(
            "sampled thermodynamic data require volume and temperature axes"
        )

    sampled_volume = np.asarray(harmonic_result.volume, dtype=np.float64)
    temperatures = np.asarray(result.temperature, dtype=np.float64)
    equilibrium_volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    valid_mask = np.ones_like(equilibrium_volume, dtype=bool)
    if result.valid_mask is not None:
        valid_mask = np.asarray(result.valid_mask, dtype=bool)

    for target_attr, source_attr in _PROPERTY_MAP:
        source = getattr(harmonic_result, source_attr)
        if source is None:
            continue
        target = _interpolate_property(
            sampled_volume,
            source,
            equilibrium_volume,
            valid_mask,
            ntemperatures=temperatures.size,
            degree=int(options.energy_degree),
        )
        setattr(result, target_attr, target)

    _finalize_equilibrium_thermodynamics(result, options)
    result.metadata.setdefault("thermodynamics", {})
    result.metadata["thermodynamics"].update(
        {
            "source": "harmonic sampled volume grid",
            "interpolation_degree": int(options.energy_degree),
            "scheme": options.scheme,
        }
    )


def _interpolate_property(
    volume: np.ndarray,
    values: np.ndarray,
    equilibrium_volume: np.ndarray,
    valid_mask: np.ndarray,
    *,
    ntemperatures: int,
    degree: int,
) -> np.ndarray:
    """Return a pressure-temperature property interpolated over volume."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 0:
        return np.full_like(equilibrium_volume, float(array), dtype=np.float64)
    if array.ndim == 1:
        if array.size == volume.size:
            rows = np.broadcast_to(array, (ntemperatures, volume.size))
        elif array.size == ntemperatures:
            out = np.full_like(equilibrium_volume, np.nan, dtype=np.float64)
            for it in range(ntemperatures):
                out[it, valid_mask[it]] = array[it]
            return out
        elif array.size == 1:
            return np.full_like(equilibrium_volume, float(array[0]), dtype=np.float64)
        else:
            raise ValueError(
                "one-dimensional property cannot be mapped to volume or temperature"
            )
    elif array.ndim == 2:
        if array.shape == (ntemperatures, volume.size):
            rows = array
        elif array.shape == (volume.size, ntemperatures):
            rows = array.T
        elif array.shape[0] == 1 and array.shape[1] == volume.size:
            rows = np.broadcast_to(array[0], (ntemperatures, volume.size))
        else:
            raise ValueError(
                "two-dimensional property shape is incompatible with QHA axes"
            )
    else:
        raise ValueError(
            "thermodynamic properties must be scalar, one- or two-dimensional"
        )

    out = np.full_like(equilibrium_volume, np.nan, dtype=np.float64)
    fit_degree = min(int(degree), max(volume.size - 1, 1))
    for it in range(ntemperatures):
        coeffs = np.polynomial.polynomial.polyfit(volume, rows[it], deg=fit_degree)
        mask = valid_mask[it]
        if np.any(mask):
            out[it, mask] = np.polynomial.polynomial.polyval(
                equilibrium_volume[it, mask], coeffs
            )
    return out


def _finalize_equilibrium_thermodynamics(
    result: QHAResult,
    options: QHAOptions,
) -> None:
    """Calculate thermodynamic properties shared by both QHA schemes.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium thermodynamic properties.
    options : QHAOptions
        Unit options used by pressure and heat-capacity transformations.
    """
    _calculate_pressure_shifted_energies(result, options)
    _calculate_thermal_expansion(result, options)
    _calculate_heat_capacity_correction(result, options)
    if options.calculate_gruneisen:
        _calculate_thermodynamic_gruneisen(result, options)


def refresh_thermal_expansion_dependents(
    result: QHAResult,
    options: QHAOptions,
) -> None:
    """Recalculate thermal expansion and all dependent thermodynamic data.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volumes, bulk moduli and heat
        capacities.
    options : QHAOptions
        QHA options selecting the thermal-expansion method and unit system.

    Returns
    -------
    None
        The result object is updated in place.
    """
    _calculate_thermal_expansion(result, options)
    _calculate_heat_capacity_correction(result, options)
    if options.calculate_gruneisen:
        _calculate_thermodynamic_gruneisen(result, options)


def _calculate_pressure_shifted_energies(
    result: QHAResult, options: QHAOptions
) -> None:
    """Calculate enthalpy and Gibbs free energy from interpolated values."""
    if result.equilibrium_volume is None or result.pressure is None:
        return
    pressure = np.asarray(result.pressure, dtype=np.float64)
    volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    pv = np.full_like(volume, np.nan, dtype=np.float64)
    for ip, pressure_value in enumerate(pressure):
        pv[:, ip] = (
            pressure_energy_density(float(pressure_value), options) * volume[:, ip]
        )

    if result.internal_energy is not None:
        result.enthalpy = np.asarray(result.internal_energy, dtype=np.float64) + pv
    if result.free_energy is not None:
        result.gibbs_free_energy = np.asarray(result.free_energy, dtype=np.float64) + pv


def _calculate_thermal_expansion(
    result: QHAResult,
    options: QHAOptions | None = None,
) -> None:
    """Calculate and select the volumetric thermal-expansion coefficient.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volumes and temperatures.
    options : QHAOptions or None, optional
        Unit options used to convert temperature increments to kelvin and to
        select the preferred thermal-expansion method.
    """
    if result.equilibrium_volume is None or result.temperature is None:
        return
    temperature = np.asarray(result.temperature, dtype=np.float64)
    if options is not None:
        temperature = np.asarray(
            convert_temperature(temperature, options.temperature_unit, "K"),
            dtype=np.float64,
        )
    volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    if temperature.size < 2:
        if result.thermal_expansion_numerical is None:
            result.thermal_expansion_numerical = np.full_like(
                volume, np.nan, dtype=np.float64
            )
        _select_thermal_expansion(result, options)
        return
    alpha_numerical = np.full_like(volume, np.nan, dtype=np.float64)
    for ip in range(volume.shape[1]):
        column = volume[:, ip]
        valid = np.isfinite(column) & (column != 0.0)
        if np.count_nonzero(valid) < 2:
            continue
        alpha_numerical[valid, ip] = relative_derivative(
            temperature[valid],
            column[valid],
        )
    result.thermal_expansion_numerical = alpha_numerical
    _select_thermal_expansion(result, options)


def calculate_mixed_thermal_expansion(
    result: QHAResult,
    sampled_volume: np.ndarray,
    sampled_entropy: np.ndarray,
    options: QHAOptions,
) -> None:
    r"""Calculate thermal expansion from the mixed free-energy derivative.

    The mixed derivative is evaluated through the Maxwell relation

    .. math::

        \alpha_V
        = -\frac{1}{K_T}\frac{\partial^2 F}{\partial V\,\partial T}
        = \frac{1}{K_T}\left(\frac{\partial S}{\partial V}\right)_T.

    Entropy is fitted as a function of volume independently at each
    temperature. The fitted polynomial is differentiated analytically and
    evaluated at each pressure-dependent equilibrium volume.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volumes, temperatures and isothermal
        bulk moduli.
    sampled_volume : ndarray
        One-dimensional sampled volume grid.
    sampled_entropy : ndarray
        Entropy sampled on the temperature-volume grid, expressed in the
        native energy unit per cell and kelvin.
    options : QHAOptions
        QHA options defining polynomial degree and unit conventions.

    Returns
    -------
    None
        ``result.thermal_expansion_mixed`` is updated in place.

    Raises
    ------
    ValueError
        If the sampled arrays are incompatible with the QHA result axes.
    """
    required = (
        result.equilibrium_volume,
        result.temperature,
        result.isothermal_bulk_modulus,
    )
    if any(item is None for item in required):
        return

    volume_grid = np.asarray(sampled_volume, dtype=np.float64)
    entropy_grid = np.asarray(sampled_entropy, dtype=np.float64)
    equilibrium_volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    kt = np.asarray(result.isothermal_bulk_modulus, dtype=np.float64)
    temperature = np.asarray(result.temperature, dtype=np.float64)

    if volume_grid.ndim != 1 or volume_grid.size < 2:
        raise ValueError(
            "sampled_volume must be one-dimensional with at least two values"
        )
    if entropy_grid.shape != (temperature.size, volume_grid.size):
        raise ValueError("sampled_entropy must have shape (ntemperatures, nvolumes)")
    if equilibrium_volume.shape != kt.shape:
        raise ValueError("equilibrium volume and bulk modulus shapes must match")
    if equilibrium_volume.shape[0] != temperature.size:
        raise ValueError("equilibrium arrays must match the temperature axis")

    alpha = np.full_like(equilibrium_volume, np.nan, dtype=np.float64)
    fit_degree = min(int(options.energy_degree), int(volume_grid.size) - 1)
    fitted_temperatures = 0

    for it in range(temperature.size):
        entropy_row = entropy_grid[it]
        finite_row = np.isfinite(volume_grid) & np.isfinite(entropy_row)
        if np.count_nonzero(finite_row) <= fit_degree:
            continue
        mask = (
            np.isfinite(equilibrium_volume[it])
            & (equilibrium_volume[it] > 0.0)
            & np.isfinite(kt[it])
            & (kt[it] > 0.0)
        )
        if not np.any(mask):
            continue
        entropy_volume_derivative = polynomial_derivative(
            volume_grid[finite_row],
            entropy_row[finite_row],
            equilibrium_volume[it, mask],
            fit_degree,
        )
        pressure_temperature_derivative = np.asarray(
            energy_to_pressure(
                entropy_volume_derivative,
                options.energy_unit,
                options.volume_unit,
                options.pressure_unit,
            ),
            dtype=np.float64,
        )
        alpha[it, mask] = pressure_temperature_derivative / kt[it, mask]
        fitted_temperatures += 1

    result.thermal_expansion_mixed = alpha
    metadata = result.metadata.setdefault("thermal_expansion", {})
    metadata["mixed_derivative"] = {
        "definition": "-(1 / KT) * d2F/(dV dT) = (1 / KT) * dS/dV",
        "surface_representation": "temperature-wise polynomial entropy-volume fits",
        "polynomial_degree": int(fit_degree),
        "sampled_volumes": int(volume_grid.size),
        "fitted_temperatures": int(fitted_temperatures),
        "pressure_unit": options.pressure_unit,
        "temperature_unit": "K",
        "evaluated_points": int(np.count_nonzero(np.isfinite(alpha))),
        "unresolved_points": int(np.count_nonzero(~np.isfinite(alpha))),
    }


def calculate_mode_thermal_expansion(
    result: QHAResult,
    options: QHAOptions,
) -> None:
    """Calculate thermal expansion from mode-weighted Grüneisen values.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volume, isothermal bulk modulus,
        isochoric heat capacity and heat-capacity-weighted mode Grüneisen
        values.
    options : QHAOptions
        Unit conventions used by the thermodynamic conversion.

    Returns
    -------
    None
        ``result.thermal_expansion_mode`` is updated in place.
    """
    required = (
        result.equilibrium_volume,
        result.isothermal_bulk_modulus,
        result.isochoric_heat_capacity,
        result.mode_weighted_gruneisen,
    )
    if any(item is None for item in required):
        return

    volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    kt = np.asarray(result.isothermal_bulk_modulus, dtype=np.float64)
    cv = np.asarray(result.isochoric_heat_capacity, dtype=np.float64)
    gamma = np.asarray(result.mode_weighted_gruneisen, dtype=np.float64)
    if not (volume.shape == kt.shape == cv.shape == gamma.shape):
        return

    volume_m3 = np.asarray(
        convert_volume(volume, options.volume_unit, "m"), dtype=np.float64
    )
    kt_pa = np.asarray(
        convert_pressure(kt, options.pressure_unit, "Pa"), dtype=np.float64
    )
    cv_j_mol_k = np.asarray(
        convert_energy_per_temperature(
            cv,
            f"{options.energy_unit} cell^-1 K^-1",
            "J mol^-1 K^-1",
        ),
        dtype=np.float64,
    )

    alpha = np.full_like(volume, np.nan, dtype=np.float64)
    valid = (
        np.isfinite(volume_m3)
        & (volume_m3 > 0.0)
        & np.isfinite(kt_pa)
        & (kt_pa > 0.0)
        & np.isfinite(cv_j_mol_k)
        & (cv_j_mol_k >= 0.0)
        & np.isfinite(gamma)
    )
    alpha[valid] = (
        gamma[valid] * cv_j_mol_k[valid] / (kt_pa[valid] * volume_m3[valid] * N)
    )
    if result.temperature is not None:
        temperature_k = np.asarray(
            convert_temperature(result.temperature, options.temperature_unit, "K"),
            dtype=np.float64,
        )
        zero_temperature = np.isclose(temperature_k, 0.0, rtol=0.0, atol=1.0e-12)
        alpha[zero_temperature, :] = 0.0

    result.thermal_expansion_mode = alpha
    metadata = result.metadata.setdefault("thermal_expansion", {})
    metadata["mode_gruneisen"] = {
        "definition": "gamma_mode * Cv / (KT * V)",
        "evaluated_points": int(np.count_nonzero(np.isfinite(alpha))),
        "unresolved_points": int(np.count_nonzero(~np.isfinite(alpha))),
    }


def _select_thermal_expansion(
    result: QHAResult,
    options: QHAOptions | None,
) -> None:
    """Select the requested thermal-expansion estimate with local fallback."""
    numerical = result.thermal_expansion_numerical
    if numerical is None:
        return
    numerical_array = np.asarray(numerical, dtype=np.float64)
    method = "numerical" if options is None else options.thermal_expansion_method

    if method == "mixed_derivative":
        primary = result.thermal_expansion_mixed
        primary_code = THERMAL_EXPANSION_SOURCE_CODES["mixed_derivative"]
    elif method == "mode_gruneisen":
        primary = result.thermal_expansion_mode
        primary_code = THERMAL_EXPANSION_SOURCE_CODES["mode_gruneisen"]
    else:
        primary = numerical_array
        primary_code = THERMAL_EXPANSION_SOURCE_CODES["numerical"]

    selected = np.full_like(numerical_array, np.nan, dtype=np.float64)
    source = np.full(
        numerical_array.shape,
        THERMAL_EXPANSION_SOURCE_CODES["invalid"],
        dtype=np.int8,
    )
    if primary is not None:
        primary_array = np.asarray(primary, dtype=np.float64)
        if primary_array.shape == selected.shape:
            mask = np.isfinite(primary_array)
            selected[mask] = primary_array[mask]
            source[mask] = primary_code

    fallback = ~np.isfinite(selected) & np.isfinite(numerical_array)
    selected[fallback] = numerical_array[fallback]
    source[fallback] = (
        THERMAL_EXPANSION_SOURCE_CODES["numerical"]
        if method == "numerical"
        else THERMAL_EXPANSION_SOURCE_CODES["numerical_fallback"]
    )

    result.thermal_expansion = selected
    result.thermal_expansion_source = source
    source_counts = {
        name: int(np.count_nonzero(source == code))
        for name, code in THERMAL_EXPANSION_SOURCE_CODES.items()
    }
    used_sources = [
        name
        for name in (
            "mixed_derivative",
            "mode_gruneisen",
            "numerical",
            "numerical_fallback",
        )
        if source_counts[name] > 0
    ]
    effective_method = "+".join(used_sources) if used_sources else "unavailable"
    metadata = result.metadata.setdefault("thermal_expansion", {})
    metadata.update(
        {
            "requested_method": method,
            "selected_method": effective_method,
            "fallback_method": "numerical",
            "source_codes": dict(THERMAL_EXPANSION_SOURCE_CODES),
            "source_counts": source_counts,
        }
    )


def _calculate_heat_capacity_correction(result: QHAResult, options: QHAOptions) -> None:
    """Estimate Cp-Cv, Cp and the adiabatic bulk modulus.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volumes, temperatures, isothermal
        bulk moduli, thermal expansion and isochoric heat capacity.
    options : QHAOptions
        Unit options defining the energy, pressure and volume scales.
    """
    required = (
        result.equilibrium_volume,
        result.temperature,
        result.isothermal_bulk_modulus,
        result.thermal_expansion,
        result.isochoric_heat_capacity,
    )
    if any(item is None for item in required):
        return

    volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    temperature = np.asarray(result.temperature, dtype=np.float64)
    kt = np.asarray(result.isothermal_bulk_modulus, dtype=np.float64)
    alpha = np.asarray(result.thermal_expansion, dtype=np.float64)
    cv = np.asarray(result.isochoric_heat_capacity, dtype=np.float64)

    if (
        cv.shape != volume.shape
        or kt.shape != volume.shape
        or alpha.shape != volume.shape
    ):
        return

    temperature_k = np.asarray(
        convert_temperature(temperature, options.temperature_unit, "K"),
        dtype=np.float64,
    )

    kt_pa = np.asarray(
        convert_pressure(kt, options.pressure_unit, "Pa"),
        dtype=np.float64,
    )

    volume_m3 = np.asarray(
        convert_volume(volume, options.volume_unit, "m"),
        dtype=np.float64,
    )

    temperature_grid = temperature_k[:, np.newaxis]

    delta_j_mol_k = alpha**2 * kt_pa * volume_m3 * N * temperature_grid

    delta = np.asarray(
        convert_energy_per_temperature(
            delta_j_mol_k,
            "J mol^-1 K^-1",
            f"{options.energy_unit} cell^-1 K^-1",
        ),
        dtype=np.float64,
    )

    cp = cv + delta

    ks = np.full_like(volume, np.nan, dtype=np.float64)

    zero_temperature = np.isclose(temperature_k, 0.0, rtol=0.0, atol=1.0e-12)
    ##    if np.any(zero_temperature):
    ##        ks[zero_temperature, :] = kt[zero_temperature, :]
    ##
    ##    mask = np.isfinite(cv) & np.isfinite(cp) & (cv != 0.0)
    ##    mask &= ~zero_temperature[:, np.newaxis]
    ##    ks[mask] = kt[mask] * cp[mask] / cv[mask]
    zero_temperature = np.isclose(
        temperature_k,
        0.0,
        rtol=0.0,
        atol=1.0e-12,
    )

    zero_temperature_grid = zero_temperature[:, np.newaxis]

    finite_cv = np.isfinite(cv)
    cv_abs_max = np.nanmax(np.abs(cv[finite_cv])) if np.any(finite_cv) else 0.0
    cv_threshold = max(1.0e-14, 1.0e-10 * cv_abs_max)

    unstable_cv = (~finite_cv) | (cv <= cv_threshold)

    fallback_mask = zero_temperature_grid | unstable_cv
    ks[fallback_mask] = kt[fallback_mask]

    valid_mask = (
        np.isfinite(cv) & np.isfinite(cp) & (cv > cv_threshold) & ~zero_temperature_grid
    )

    ks[valid_mask] = kt[valid_mask] * cp[valid_mask] / cv[valid_mask]

    result.heat_capacity_difference = delta
    result.isobaric_heat_capacity = cp
    result.adiabatic_bulk_modulus = ks


def _calculate_thermodynamic_gruneisen(
    result: QHAResult,
    options: QHAOptions,
) -> None:
    """Calculate the macroscopic thermodynamic Grüneisen parameter.

    Parameters
    ----------
    result : QHAResult
        QHA result containing ``V``, ``K_T``, ``alpha_V`` and ``C_V``.
    options : QHAOptions
        Unit conventions used by the calculation.
    """
    required = (
        result.equilibrium_volume,
        result.isothermal_bulk_modulus,
        result.thermal_expansion,
        result.isochoric_heat_capacity,
    )
    if any(item is None for item in required):
        return
    volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    kt = np.asarray(result.isothermal_bulk_modulus, dtype=np.float64)
    alpha = np.asarray(result.thermal_expansion, dtype=np.float64)
    cv = np.asarray(result.isochoric_heat_capacity, dtype=np.float64)
    if (
        kt.shape != volume.shape
        or alpha.shape != volume.shape
        or cv.shape != volume.shape
    ):
        return

    volume_m3 = np.asarray(
        convert_volume(volume, options.volume_unit, "m"), dtype=np.float64
    )
    kt_pa = np.asarray(
        convert_pressure(kt, options.pressure_unit, "Pa"), dtype=np.float64
    )
    cv_j_mol_k = np.asarray(
        convert_energy_per_temperature(
            cv,
            f"{options.energy_unit} cell^-1 K^-1",
            "J mol^-1 K^-1",
        ),
        dtype=np.float64,
    )
    input_metadata = result.metadata.get("input", {})
    natoms = (
        int(input_metadata.get("natoms", 0)) if isinstance(input_metadata, dict) else 0
    )
    if natoms > 0:
        dulong_petit = 3.0 * float(natoms) * cs.R
    else:
        finite_cv = cv_j_mol_k[np.isfinite(cv_j_mol_k)]
        dulong_petit = float(np.nanmax(np.abs(finite_cv))) if finite_cv.size else 0.0
    threshold = max(
        1.0e-20,
        float(options.gruneisen_min_cv_fraction) * dulong_petit,
    )

    gamma = np.full_like(volume, np.nan, dtype=np.float64)
    valid = (
        np.isfinite(volume_m3)
        & np.isfinite(kt_pa)
        & np.isfinite(alpha)
        & np.isfinite(cv_j_mol_k)
        & (cv_j_mol_k >= threshold)
    )
    gamma[valid] = (
        volume_m3[valid] * kt_pa[valid] * alpha[valid] * N / cv_j_mol_k[valid]
    )
    if result.temperature is not None:
        temperature_k = np.asarray(
            convert_temperature(result.temperature, options.temperature_unit, "K"),
            dtype=np.float64,
        )
        zero_temperature = np.isclose(temperature_k, 0.0, rtol=0.0, atol=1.0e-12)
        gamma[zero_temperature, :] = 0.0
    result.gruneisen = gamma
    result.metadata.setdefault("gruneisen", {})
    result.metadata["gruneisen"].update(
        {
            "thermodynamic_definition": "alphaV * KT * V / Cv",
            "low_heat_capacity_policy": "nan",
            "minimum_cv_fraction_of_dulong_petit": float(
                options.gruneisen_min_cv_fraction
            ),
            "dulong_petit_j_mol_cell_k": float(dulong_petit),
            "cv_threshold_j_mol_k": float(threshold),
            "unresolved_points": int(np.count_nonzero(~np.isfinite(gamma))),
        }
    )


def calculate_mode_gruneisen_at_equilibrium(
    result: QHAResult,
    input_data: QHAInput,
    options: QHAOptions,
    *,
    frequency_evaluator: FrequencyThermodynamicEvaluator | None = None,
) -> None:
    """Calculate mode-resolved and heat-capacity-weighted Grüneisen values.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volumes and temperatures.
    input_data : QHAInput
        Mode-continuous phonon data and q-point weights.
    options : QHAOptions
        Fit degrees, unit conventions and Grüneisen settings.
    frequency_evaluator : FrequencyThermodynamicEvaluator or None, optional
        Prepared frequency evaluator reused from the frequency QHA workflow.

    Raises
    ------
    ValueError
        If mode continuity is unavailable or the mode fits cannot be built.
    """
    if result.equilibrium_volume is None or result.temperature is None:
        return
    if (
        input_data.frequencies is None
        or input_data.volume is None
        or input_data.weights is None
    ):
        raise ValueError(
            "mode Grüneisen analysis requires frequencies, volumes and q-point weights"
        )
    if not input_data.has_mode_continuity():
        raise ValueError(
            f"mode continuity status '{input_data.mode_continuity_status()}' does not permit mode analysis"
        )

    prepared_frequency = frequency_evaluator or FrequencyThermodynamicEvaluator(
        input_data, options
    )
    sampled_gamma = prepared_frequency.mode_gruneisen_at(input_data.volume)
    result.mode_gruneisen = sampled_gamma
    mode_axes = sampled_gamma.shape[:-1]
    flat_sampled = sampled_gamma.reshape((-1, sampled_gamma.shape[-1]))
    usable_mask = np.all(np.isfinite(flat_sampled), axis=1)
    failed_modes = [
        tuple(int(component) for component in np.unravel_index(index, mode_axes))
        for index in np.flatnonzero(~usable_mask)
    ]
    if failed_modes and not options.gruneisen_allow_nonpositive:
        raise ValueError("non-positive frequencies prevent mode Grüneisen evaluation")
    equilibrium_volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    temperature = np.asarray(result.temperature, dtype=np.float64)
    temperature_k = np.asarray(
        convert_temperature(temperature, options.temperature_unit, "K"),
        dtype=np.float64,
    )
    valid_mask = np.ones_like(equilibrium_volume, dtype=bool)
    if result.valid_mask is not None:
        valid_mask = np.asarray(result.valid_mask, dtype=bool)
    weighted = np.full_like(equilibrium_volume, np.nan, dtype=np.float64)
    denominator = np.full_like(equilibrium_volume, np.nan, dtype=np.float64)

    for it, temperature_value in enumerate(temperature_k):
        mask = valid_mask[it] & np.isfinite(equilibrium_volume[it])
        if not np.any(mask):
            continue
        volumes = equilibrium_volume[it, mask]
        gamma_modes = prepared_frequency.mode_gruneisen_at(volumes)
        frequencies = prepared_frequency.frequencies_at(volumes)
        frequencies_hz = np.asarray(
            convert_frequency(frequencies, options.frequency_unit, "Hz"),
            dtype=np.float64,
        )
        thermal = thermal_gruneisen_from_modes(
            gamma_modes,
            frequencies_hz,
            float(temperature_value),
            input_data.normalized_weights(),
        )
        if not thermal.success or thermal.gamma is None or thermal.denominator is None:
            raise ValueError(
                thermal.message or "mode-weighted Grüneisen calculation failed"
            )
        weighted[it, mask] = thermal.gamma
        denominator[it, mask] = thermal.denominator

    result.mode_weighted_gruneisen = weighted
    calculate_mode_thermal_expansion(result, options)
    metadata = result.metadata.setdefault("gruneisen", {})
    metadata.update(
        {
            "mode_continuity": input_data.mode_continuity_status(),
            "mode_degree": int(prepared_frequency.frequency_degree),
            "mode_source": "derivative of fitted frequency-volume polynomials",
            "mode_status": "success" if not failed_modes else "partial",
            "n_modes": int(len(flat_sampled)),
            "n_usable_modes": int(np.count_nonzero(usable_mask)),
            "failed_modes": [list(index) for index in failed_modes],
            "mode_fit_r_squared": _frequency_fit_r_squared_summary(
                prepared_frequency.frequency_fit.fits,
                excluded_modes=set(failed_modes),
            ),
            "mode_weighting": "harmonic_mode_Cv",
            "zero_weight_states": int(np.count_nonzero(denominator <= 0.0)),
        }
    )
    if result.gruneisen is not None:
        macroscopic = np.asarray(result.gruneisen, dtype=np.float64)
        comparison = (
            np.isfinite(macroscopic)
            & np.isfinite(weighted)
            & valid_mask
            & (temperature_k[:, None] > 0.0)
        )
        if np.any(comparison):
            difference = weighted[comparison] - macroscopic[comparison]
            scale = np.maximum(np.abs(macroscopic[comparison]), 1.0e-12)
            metadata["consistency"] = {
                "n_points": int(np.count_nonzero(comparison)),
                "maximum_absolute_difference": float(np.max(np.abs(difference))),
                "rms_difference": float(np.sqrt(np.mean(difference**2))),
                "maximum_relative_difference": float(
                    np.max(np.abs(difference) / scale)
                ),
            }


def _frequency_fit_r_squared_summary(
    fits: dict[tuple[int, ...], Any],
    *,
    excluded_modes: set[tuple[int, ...]] | None = None,
) -> dict[str, float]:
    """Return R-squared statistics for usable frequency fits.

    Parameters
    ----------
    fits : dict
        Individual frequency-fit diagnostics keyed by mode index.
    excluded_modes : set of tuple or None, optional
        Modes excluded from Grüneisen analysis.

    Returns
    -------
    dict
        Minimum, mean, median and maximum R-squared values.
    """
    excluded = excluded_modes or set()
    values = np.asarray(
        [
            fit.r_squared
            for index, fit in fits.items()
            if index not in excluded
            and fit.r_squared is not None
            and np.isfinite(fit.r_squared)
        ],
        dtype=np.float64,
    )
    if not values.size:
        return {}
    return {
        "minimum": float(np.min(values)),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "maximum": float(np.max(values)),
    }
