# -*- coding: utf-8 -*-

"""Data containers used by the quasi-harmonic approximation workflow.

The classes defined here describe normalized input data, workflow options,
calculated properties, uncertainties, and diagnostic records for QHA
calculations. They are passive containers shared by the Python API, command-line
interface, graphical interface, exporters, and tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping

import numpy as np

from quantas.core.math.fitting import FitResult
from quantas.models.phonons import PhononInputData

QHAScheme = Literal["freq", "td"]
QHAMinimization = Literal["poly", "eos"]
QHAFitFailurePolicy = Literal["continue", "stop", "raise"]
QHAFitQualityPolicy = Literal["warn", "stop"]
QHAExtrapolationPolicy = Literal["warn", "fail", "allow"]
QHAUncertaintyMethod = Literal["none", "covariance", "bootstrap", "montecarlo"]
QHAPolynomialDerivativeMethod = Literal["local_grid", "analytic"]
QHAThermalExpansionMethod = Literal[
    "mixed_derivative",
    "mode_gruneisen",
    "numerical",
]
QHAModeContinuity = Literal["verified", "assumed", "unknown", "unreliable"]


@dataclass(slots=True)
class QHAInput(PhononInputData):
    """Input data for a quasi-harmonic approximation calculation.

    Parameters
    ----------
    mode_continuity : {"verified", "assumed", "unknown", "unreliable"}, optional
        Status describing whether phonon-mode ordering is continuous across the
        sampled volume sequence.
    """

    mode_continuity: QHAModeContinuity = "assumed"

    def has_mode_continuity_data(self) -> bool:
        """Return whether the array layout can represent continuous modes.

        Returns
        -------
        bool
            ``True`` when frequencies have shape ``(qpoints, modes, volumes)``
            and the volume axis matches the sampled-volume array.
        """
        if self.frequencies is None or self.volume is None:
            return False
        frequencies = np.asarray(self.frequencies)
        return frequencies.ndim == 3 and frequencies.shape[-1] == self.nvol

    def mode_continuity_status(self) -> str:
        """Return the normalized phonon-mode continuity status.

        Returns
        -------
        str
            One of ``verified``, ``assumed``, ``unknown``, or ``unreliable``.
        """
        return str(self.mode_continuity).strip().lower()

    def has_mode_continuity(self) -> bool:
        """Return whether mode-continuous analysis is permitted.

        Returns
        -------
        bool
            ``True`` when the data layout is suitable and continuity is marked
            either ``verified`` or ``assumed``.
        """
        return self.has_mode_continuity_data() and self.mode_continuity_status() in {
            "verified",
            "assumed",
        }

    def has_verified_mode_continuity(self) -> bool:
        """Return whether mode continuity was explicitly verified.

        Returns
        -------
        bool
            ``True`` only for structurally valid data marked ``verified``.
        """
        return (
            self.has_mode_continuity_data()
            and self.mode_continuity_status() == "verified"
        )

    def validate_shapes(self) -> None:
        """Validate array dimensions required by QHA calculations.

        Raises
        ------
        ValueError
            If structural, energetic, q-point, or phonon arrays are
            inconsistent.
        """
        if self.formula_units <= 0:
            raise ValueError("formula_units must be positive")
        if self.mode_continuity_status() not in {
            "verified",
            "assumed",
            "unknown",
            "unreliable",
        }:
            raise ValueError(
                "mode_continuity must be verified, assumed, unknown, or unreliable"
            )
        if self.volume is None:
            raise ValueError("volume data are required")
        volume = np.asarray(self.volume)
        if volume.ndim != 1:
            raise ValueError("volume must be a one-dimensional array")
        if self.energy is None:
            raise ValueError("static energy data are required")
        energy = np.asarray(self.energy)
        if energy.ndim != 1 or energy.shape[0] != volume.shape[0]:
            raise ValueError("energy must be one-dimensional and match volume")
        if self.frequencies is not None:
            frequencies = np.asarray(self.frequencies)
            if frequencies.ndim != 3:
                raise ValueError("frequencies must have shape (qpoints, nmodes, nvol)")
            if frequencies.shape[-1] != volume.shape[0]:
                raise ValueError("frequency volume axis must match volume")
            if self.qpoints and frequencies.shape[0] != self.qpoints:
                raise ValueError("frequency q-point axis must match qpoints")
        if self.weights is not None:
            weights = np.asarray(self.weights)
            if weights.ndim != 1:
                raise ValueError("weights must be one-dimensional")
            expected = self.qpoints or (
                0 if self.frequencies is None else np.asarray(self.frequencies).shape[0]
            )
            if expected and weights.shape[0] != expected:
                raise ValueError("weights must match the number of q-points")
        if self.qcoords is not None:
            qcoords = np.asarray(self.qcoords)
            if qcoords.ndim != 2 or qcoords.shape[1] != 3:
                raise ValueError("qcoords must have shape (qpoints, 3)")
            expected = self.qpoints or qcoords.shape[0]
            if qcoords.shape[0] != expected:
                raise ValueError("qcoords must match the number of q-points")


@dataclass(slots=True)
class QHAOptions:
    """Options controlling a quasi-harmonic approximation calculation.

    Parameters
    ----------
    temperature_min : float, optional
        Minimum temperature of the calculation range.
    temperature_max : float, optional
        Maximum temperature of the calculation range.
    temperature_step : float, optional
        Temperature increment.
    pressure_min : float, optional
        Minimum pressure of the calculation range.
    pressure_max : float, optional
        Maximum pressure of the calculation range.
    pressure_step : float, optional
        Pressure increment.
    scheme : {"freq", "td"}, optional
        QHA scheme used for volume-dependent thermodynamics.
    minimization : {"poly", "eos"}, optional
        Method used to determine equilibrium volumes at each ``P, T`` point.
    eos : str, optional
        Equation of state used when ``minimization`` is ``"eos"``.
    energy_degree : int, optional
        Polynomial degree used for static energy fits.
    free_energy_degree : int, optional
        Polynomial degree used for Helmholtz free-energy fits.
    frequency_degree : int, optional
        Polynomial degree used for mode-resolved frequency fits.
    structural_degree : int, optional
        Polynomial degree used for the deviatoric logarithmic-strain path when
        structural cells are available. Axial and tensorial thermal expansion
        are calculated automatically from this path.
    polynomial_derivative_method : {"local_grid", "analytic"}, optional
        Method used to calculate polynomial ``K_T`` and ``K'_T`` values.
    polynomial_grid_points : int, optional
        Number of volumes in the local derivative grid.
    polynomial_grid_separation : float, optional
        Adjacent local-grid spacing as a percentage of the equilibrium volume.
    energy_unit : str, optional
        Energy unit used for input static energies and final energy-like
        results.
    volume_unit : str, optional
        Length unit used to express unit-cell volumes.
    frequency_unit : str, optional
        Frequency unit used by the input phonon data.
    temperature_unit : str, optional
        Temperature unit used by the input temperature range.
    pressure_unit : str, optional
        Pressure unit used by the input pressure range.
    debug : bool, optional
        If ``True``, store and expose detailed fit diagnostics.
    store_fit_diagnostics : bool, optional
        If ``True``, retain structured diagnostics from local fits.
    estimate_uncertainties : bool, optional
        If ``True``, request uncertainty estimates where supported.
    uncertainty_method : str, optional
        Method requested for uncertainty estimates.
    uncertainty_relative_step : float, optional
        Relative EOS-parameter step used by linear covariance propagation.
    uncertainty_confidence_level : float, optional
        Confidence probability used for propagated EOS intervals.
    uncertainty_samples : int, optional
        Number of correlated EOS parameter samples used by Monte Carlo
        propagation.
    uncertainty_seed : int or None, optional
        Random seed used by Monte Carlo propagation.
    uncertainty_minimum_success_fraction : float, optional
        Minimum accepted fraction of physical Monte Carlo EOS states.
    extrapolation_policy : {"warn", "fail", "allow"}, optional
        Policy applied when an equilibrium volume is outside the sampled
        volume interval.
    fit_failure_policy : {"continue", "stop", "raise"}, optional
        Policy applied when a local fit fails.
    fit_quality_policy : {"warn", "stop"}, optional
        Policy applied when a local fit is successful but poor.
    max_consecutive_failures : int, optional
        Maximum number of consecutive local fit failures before a workflow
        stops when the failure policy is ``"stop"``.
    calculate_gruneisen : bool, optional
        If ``True``, calculate the thermodynamic Grüneisen parameter from
        ``alpha_V K_T V / C_V``.
    calculate_mode_gruneisen : bool, optional
        If ``True`` and the frequency scheme is selected, calculate
        mode-resolved and heat-capacity-weighted Grüneisen parameters.
    thermal_expansion_method : {"mixed_derivative", "mode_gruneisen", "numerical"}, optional
        Method used to calculate the volumetric thermal-expansion coefficient.
        The mixed derivative of the fitted free-energy surface is the default.
        The mode-Grüneisen method is available only for the frequency scheme,
        while the numerical volume derivative is retained as a fallback.
    gruneisen_allow_nonpositive : bool, optional
        If ``True``, non-positive modes are excluded from mode averages rather
        than invalidating the full calculation.
    gruneisen_min_cv_fraction : float, optional
        Minimum fraction of the Dulong-Petit heat capacity required before the
        macroscopic ``alpha_V K_T V / C_V`` ratio is considered resolved.
    metadata : dict, optional
        Additional caller-defined options.
    """

    temperature_min: float = 298.15
    temperature_max: float = 298.15
    temperature_step: float = 1.0

    pressure_min: float = 0.0
    pressure_max: float = 0.0
    pressure_step: float = 1.0

    scheme: QHAScheme = "freq"
    minimization: QHAMinimization = "poly"
    eos: str = "BM3"
    energy_degree: int = 3
    free_energy_degree: int = 3
    frequency_degree: int = 3
    structural_degree: int = 3
    polynomial_derivative_method: QHAPolynomialDerivativeMethod = "local_grid"
    polynomial_grid_points: int = 5
    polynomial_grid_separation: float = 0.05

    energy_unit: str = "Ha"
    volume_unit: str = "A"
    frequency_unit: str = "cm^-1"
    temperature_unit: str = "K"
    pressure_unit: str = "GPa"

    debug: bool = False
    store_fit_diagnostics: bool = True
    estimate_uncertainties: bool = True
    uncertainty_method: QHAUncertaintyMethod = "covariance"
    uncertainty_relative_step: float = 1.0e-5
    uncertainty_confidence_level: float = 0.95
    uncertainty_samples: int = 10000
    uncertainty_seed: int | None = None
    uncertainty_minimum_success_fraction: float = 0.8
    extrapolation_policy: QHAExtrapolationPolicy = "warn"
    fit_failure_policy: QHAFitFailurePolicy = "stop"
    fit_quality_policy: QHAFitQualityPolicy = "warn"
    max_consecutive_failures: int = 5

    calculate_gruneisen: bool = True
    calculate_mode_gruneisen: bool = True
    thermal_expansion_method: QHAThermalExpansionMethod = "mixed_derivative"
    gruneisen_allow_nonpositive: bool = True
    gruneisen_min_cv_fraction: float = 1.0e-2

    metadata: dict[str, Any] = field(default_factory=dict)

    def temperature_grid(self) -> np.ndarray:
        """Return the temperature grid used by the QHA workflow.

        Returns
        -------
        ndarray
            One-dimensional temperature array generated with the Quantas grid
            convention ``np.arange(min, max + step, step)``.

        Raises
        ------
        ValueError
            If the temperature step is not positive or the maximum temperature
            is smaller than the minimum temperature.
        """
        if self.temperature_step <= 0.0:
            raise ValueError("temperature_step must be positive")
        if self.temperature_max < self.temperature_min:
            raise ValueError(
                "temperature_max must be greater than or equal to temperature_min"
            )
        return np.arange(
            self.temperature_min,
            self.temperature_max + self.temperature_step,
            self.temperature_step,
            dtype=np.float64,
        )

    def pressure_grid(self) -> np.ndarray:
        """Return the pressure grid used by the QHA workflow.

        Returns
        -------
        ndarray
            One-dimensional pressure array generated with the Quantas grid
            convention ``np.arange(min, max + step, step)``.

        Raises
        ------
        ValueError
            If the pressure step is not positive or the maximum pressure is
            smaller than the minimum pressure.
        """
        if self.pressure_step <= 0.0:
            raise ValueError("pressure_step must be positive")
        if self.pressure_max < self.pressure_min:
            raise ValueError(
                "pressure_max must be greater than or equal to pressure_min"
            )
        return np.arange(
            self.pressure_min,
            self.pressure_max + self.pressure_step,
            self.pressure_step,
            dtype=np.float64,
        )

    def requires_mode_continuity(self) -> bool:
        """Return whether the selected options require mode continuity.

        Returns
        -------
        bool
            ``True`` for methods that operate directly on mode-resolved
            frequency changes over volume.
        """
        return self.scheme == "freq"

    def validate(self) -> None:
        """Validate option values before running a QHA workflow.

        Raises
        ------
        ValueError
            If an option has an unsupported value or an invalid numerical
            setting.
        """
        if self.scheme not in {"freq", "td"}:
            raise ValueError("scheme must be 'freq' or 'td'")
        if self.thermal_expansion_method not in {
            "mixed_derivative",
            "mode_gruneisen",
            "numerical",
        }:
            raise ValueError(
                "thermal_expansion_method must be 'mixed_derivative', "
                "'mode_gruneisen', or 'numerical'"
            )
        if self.thermal_expansion_method == "mode_gruneisen" and self.scheme != "freq":
            raise ValueError(
                "thermal_expansion_method='mode_gruneisen' requires scheme='freq'"
            )
        if self.minimization not in {"poly", "eos"}:
            raise ValueError("minimization must be 'poly' or 'eos'")
        if self.energy_degree < 1:
            raise ValueError("energy_degree must be at least 1")
        if self.free_energy_degree < 1:
            raise ValueError("free_energy_degree must be at least 1")
        if self.frequency_degree < 1:
            raise ValueError("frequency_degree must be at least 1")
        if self.structural_degree < 1:
            raise ValueError("structural_degree must be at least 1")
        if not 0.0 <= self.gruneisen_min_cv_fraction < 1.0:
            raise ValueError("gruneisen_min_cv_fraction must be in [0, 1)")
        if self.polynomial_derivative_method not in {"local_grid", "analytic"}:
            raise ValueError(
                "polynomial_derivative_method must be 'local_grid' or 'analytic'"
            )
        if self.polynomial_grid_points < 3 or self.polynomial_grid_points % 2 == 0:
            raise ValueError(
                "polynomial_grid_points must be an odd integer of at least 3"
            )
        if self.polynomial_grid_separation <= 0.0:
            raise ValueError("polynomial_grid_separation must be positive")
        if self.max_consecutive_failures < 1:
            raise ValueError("max_consecutive_failures must be at least 1")
        if self.uncertainty_method not in {
            "none",
            "covariance",
            "bootstrap",
            "montecarlo",
        }:
            raise ValueError("unsupported uncertainty_method")
        if self.uncertainty_relative_step <= 0.0:
            raise ValueError("uncertainty_relative_step must be positive")
        if not 0.0 < self.uncertainty_confidence_level < 1.0:
            raise ValueError(
                "uncertainty_confidence_level must be between zero and one"
            )
        if self.uncertainty_samples < 2:
            raise ValueError("uncertainty_samples must be at least 2")
        if not 0.0 < self.uncertainty_minimum_success_fraction <= 1.0:
            raise ValueError("uncertainty_minimum_success_fraction must be in (0, 1]")
        self.temperature_grid()
        self.pressure_grid()


@dataclass(slots=True)
class QHAFitRecord:
    """Diagnostic record for a local fit in a QHA workflow.

    Parameters
    ----------
    quantity : str
        Name of the fitted quantity.
    method : str
        Fitting or minimization method.
    temperature : float or None, optional
        Temperature associated with the fit.
    pressure : float or None, optional
        Pressure associated with the fit.
    fit : FitResult or None, optional
        Structured fit result.
    success : bool, optional
        Whether the local operation produced a usable result.
    message : str, optional
        Human-readable diagnostic message.
    metadata : dict, optional
        Additional diagnostic information.
    """

    quantity: str
    method: str
    temperature: float | None = None
    pressure: float | None = None
    fit: FitResult | None = None
    success: bool = True
    message: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        """Return the diagnostic record as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the diagnostic record.
        """
        return {
            "quantity": self.quantity,
            "method": self.method,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "fit": None if self.fit is None else self.fit.as_dict(),
            "success": self.success,
            "message": self.message,
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class QHAFailedPoint:
    """Description of a failed point in the pressure-temperature grid.

    Parameters
    ----------
    temperature : float
        Temperature at which the local workflow failed.
    pressure : float
        Pressure at which the local workflow failed.
    stage : str
        Workflow stage where the failure occurred.
    message : str
        Human-readable explanation of the failure.
    diagnostics : dict, optional
        Additional diagnostic details.
    """

    temperature: float
    pressure: float
    stage: str
    message: str
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        """Return the failed point as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the failed point.
        """
        return {
            "temperature": self.temperature,
            "pressure": self.pressure,
            "stage": self.stage,
            "message": self.message,
            "diagnostics": dict(self.diagnostics),
        }


@dataclass(slots=True)
class QHAResult:
    """Results of a quasi-harmonic approximation calculation.

    Parameters
    ----------
    jobname : str, optional
        Name or short description of the calculation.
    temperature : ndarray or None, optional
        Temperature grid used in the calculation.
    pressure : ndarray or None, optional
        Pressure grid used in the calculation.
    volume : ndarray or None, optional
        Sampled input volumes.
    static_energy : ndarray or None, optional
        Static energies associated with the sampled volumes.
    equilibrium_volume : ndarray or None, optional
        Equilibrium volume at each pressure-temperature point, historically
        stored as ``VT``.
    zero_point_energy : ndarray or None, optional
        Zero-point vibrational energy, historically stored as ``Uzp``.
    thermal_energy : ndarray or None, optional
        Thermal contribution to internal energy, historically stored as
        ``Uth``.
    internal_energy : ndarray or None, optional
        Total internal energy, historically stored as ``Utot``.
    entropy : ndarray or None, optional
        Entropy, historically stored as ``S``.
    vibrational_free_energy : ndarray or None, optional
        Vibrational Helmholtz free energy, historically stored as ``Fvib``.
    free_energy : ndarray or None, optional
        Total Helmholtz free energy, historically stored as ``F``.
    isochoric_heat_capacity : ndarray or None, optional
        Isochoric heat capacity, historically stored as ``Cv``.
    isobaric_heat_capacity : ndarray or None, optional
        Isobaric heat capacity, historically stored as ``Cp``.
    heat_capacity_difference : ndarray or None, optional
        Difference between isobaric and isochoric heat capacities, historically
        stored as ``Cp-Cv``.
    isothermal_bulk_modulus : ndarray or None, optional
        Isothermal bulk modulus, historically stored as ``KT``.
    adiabatic_bulk_modulus : ndarray or None, optional
        Adiabatic bulk modulus, historically stored as ``KS``.
    bulk_modulus_derivative : ndarray or None, optional
        First pressure derivative of the bulk modulus, historically stored as
        ``Kp``.
    thermal_expansion : ndarray or None, optional
        Volumetric thermal-expansion coefficient, historically stored as
        ``alphaV``.
    thermal_expansion_mixed : ndarray or None, optional
        Thermal expansion from the mixed derivative of the fitted free-energy
        surface.
    thermal_expansion_mode : ndarray or None, optional
        Thermal expansion reconstructed from heat-capacity-weighted mode
        Grüneisen parameters.
    thermal_expansion_numerical : ndarray or None, optional
        Thermal expansion from the numerical derivative of equilibrium volume.
    thermal_expansion_source : ndarray or None, optional
        Integer source code used for each selected thermal-expansion value.
    equilibrium_lattice : ndarray or None, optional
        Crystallographic direct-lattice matrices at every pressure-temperature
        point, with shape ``(nT, nP, 3, 3)`` and vectors stored by rows.
    lattice_parameters : ndarray or None, optional
        ``a, b, c, alpha, beta, gamma`` at every pressure-temperature point.
    lattice_parameter_derivatives : ndarray or None, optional
        Temperature derivatives of lattice lengths and angles.
    axial_thermal_expansion : ndarray or None, optional
        Linear thermal-expansion coefficients of crystallographic axes
        ``a``, ``b``, and ``c``.
    thermal_expansion_tensor : ndarray or None, optional
        Symmetric Cartesian thermal-expansion tensor.
    structural_extrapolation_mask : ndarray or None, optional
        Mask marking equilibrium volumes outside the sampled structural path.
    enthalpy : ndarray or None, optional
        Enthalpy, historically stored as ``H``.
    gibbs_free_energy : ndarray or None, optional
        Gibbs free energy, historically stored as ``G``.
    gruneisen : ndarray or None, optional
        Thermodynamic Grüneisen parameter calculated from
        ``alpha_V K_T V / C_V``.
    mode_weighted_gruneisen : ndarray or None, optional
        Heat-capacity-weighted mode Grüneisen parameter on the pressure-
        temperature grid.
    mode_gruneisen : ndarray or None, optional
        Mode-resolved Grüneisen parameters on the sampled volume grid.
    uncertainties : dict, optional
        Property uncertainties keyed by result name. Suggested keys are
        ``sigma_VT``, ``sigma_KT``, and analogous names for other properties.
    fit_records : list of QHAFitRecord, optional
        Structured diagnostics from local fits.
    failed_points : list of QHAFailedPoint, optional
        Failed pressure-temperature points.
    valid_mask : ndarray or None, optional
        Boolean mask marking valid pressure-temperature points.
    completed : bool, optional
        Whether the workflow reached the end of the requested grid.
    metadata : dict, optional
        Additional module-specific metadata.
    """

    jobname: str = "Unknown"
    temperature: np.ndarray | None = None
    pressure: np.ndarray | None = None
    volume: np.ndarray | None = None
    static_energy: np.ndarray | None = None

    equilibrium_volume: np.ndarray | None = None
    zero_point_energy: np.ndarray | None = None
    thermal_energy: np.ndarray | None = None
    internal_energy: np.ndarray | None = None
    entropy: np.ndarray | None = None
    vibrational_free_energy: np.ndarray | None = None
    free_energy: np.ndarray | None = None
    isochoric_heat_capacity: np.ndarray | None = None
    isobaric_heat_capacity: np.ndarray | None = None
    heat_capacity_difference: np.ndarray | None = None
    isothermal_bulk_modulus: np.ndarray | None = None
    adiabatic_bulk_modulus: np.ndarray | None = None
    bulk_modulus_derivative: np.ndarray | None = None
    thermal_expansion: np.ndarray | None = None
    thermal_expansion_mixed: np.ndarray | None = None
    thermal_expansion_mode: np.ndarray | None = None
    thermal_expansion_numerical: np.ndarray | None = None
    thermal_expansion_source: np.ndarray | None = None
    equilibrium_lattice: np.ndarray | None = None
    lattice_parameters: np.ndarray | None = None
    lattice_parameter_derivatives: np.ndarray | None = None
    axial_thermal_expansion: np.ndarray | None = None
    thermal_expansion_tensor: np.ndarray | None = None
    structural_extrapolation_mask: np.ndarray | None = None
    enthalpy: np.ndarray | None = None
    gibbs_free_energy: np.ndarray | None = None
    gruneisen: np.ndarray | None = None
    mode_weighted_gruneisen: np.ndarray | None = None
    mode_gruneisen: np.ndarray | None = None

    uncertainties: dict[str, np.ndarray] = field(default_factory=dict)
    fit_records: list[QHAFitRecord] = field(default_factory=list)
    failed_points: list[QHAFailedPoint] = field(default_factory=list)
    valid_mask: np.ndarray | None = None
    completed: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)

    def as_property_dict(self) -> dict[str, np.ndarray]:
        """Return calculated properties using the standard QHA property symbols.

        Returns
        -------
        dict
            Dictionary containing available properties with keys such as
            ``Uzp``, ``Uth``, ``Utot``, ``S``, ``Fvib``, ``F``, ``Cv``, ``VT``,
            ``Cp``, ``KT``, ``KS``, ``Kp``, ``alphaV``, ``H``, ``G``, and
            ``gruneisen``.
        """
        mapping = {
            "Uzp": self.zero_point_energy,
            "Uth": self.thermal_energy,
            "Utot": self.internal_energy,
            "S": self.entropy,
            "Fvib": self.vibrational_free_energy,
            "F": self.free_energy,
            "Cv": self.isochoric_heat_capacity,
            "VT": self.equilibrium_volume,
            "Cp": self.isobaric_heat_capacity,
            "Cp-Cv": self.heat_capacity_difference,
            "KT": self.isothermal_bulk_modulus,
            "KS": self.adiabatic_bulk_modulus,
            "Kp": self.bulk_modulus_derivative,
            "alphaV": self.thermal_expansion,
            "alphaV_mixed": self.thermal_expansion_mixed,
            "alphaV_mode": self.thermal_expansion_mode,
            "alphaV_numerical": self.thermal_expansion_numerical,
            "lattice": self.equilibrium_lattice,
            "cell_parameters": self.lattice_parameters,
            "cell_parameter_dT": self.lattice_parameter_derivatives,
            "alphaABC": self.axial_thermal_expansion,
            "alpha_tensor": self.thermal_expansion_tensor,
            "H": self.enthalpy,
            "G": self.gibbs_free_energy,
            "gruneisen": self.gruneisen,
            "mode_weighted_gruneisen": self.mode_weighted_gruneisen,
        }
        return {key: value for key, value in mapping.items() if value is not None}

    def uncertainty_dict(self) -> dict[str, np.ndarray]:
        """Return property uncertainties stored in the result.

        Returns
        -------
        dict
            Dictionary of uncertainty arrays keyed by property name.
        """
        return dict(self.uncertainties)

    def diagnostics_as_dict(self) -> list[dict[str, Any]]:
        """Return fit diagnostics as serializable dictionaries.

        Returns
        -------
        list of dict
            Serialized diagnostic records.
        """
        return [record.as_dict() for record in self.fit_records]

    def failed_points_as_dict(self) -> list[dict[str, Any]]:
        """Return failed pressure-temperature points as dictionaries.

        Returns
        -------
        list of dict
            Serialized failed-point records.
        """
        return [point.as_dict() for point in self.failed_points]

    def add_fit_record(self, record: QHAFitRecord) -> None:
        """Append a local fit diagnostic record.

        Parameters
        ----------
        record : QHAFitRecord
            Diagnostic record to store.
        """
        self.fit_records.append(record)

    def add_failed_point(self, point: QHAFailedPoint) -> None:
        """Append a failed pressure-temperature point.

        Parameters
        ----------
        point : QHAFailedPoint
            Failed point to store.
        """
        self.failed_points.append(point)

    def has_thermodynamic_data(self) -> bool:
        """Return whether any QHA property has been calculated.

        Returns
        -------
        bool
            ``True`` if at least one property array is available, otherwise
            ``False``.
        """
        return bool(self.as_property_dict())

    def partial_copy(self, valid_mask: np.ndarray) -> QHAResult:
        """Return a result containing only valid pressure-temperature points.

        Parameters
        ----------
        valid_mask : ndarray
            Boolean mask with the same shape as the pressure-temperature result
            arrays.

        Returns
        -------
        QHAResult
            Result object where arrays matching the mask shape are filtered.

        Raises
        ------
        ValueError
            If ``valid_mask`` is not a boolean array.
        """
        mask = np.asarray(valid_mask)
        if mask.dtype != np.bool_:
            raise ValueError("valid_mask must be a boolean array")

        def filter_value(value: np.ndarray | None) -> np.ndarray | None:
            if value is None:
                return None
            array = np.asarray(value)
            if array.shape[: mask.ndim] == mask.shape:
                return array[mask]
            return value

        uncertainties: dict[str, np.ndarray] = {}
        for key, value in self.uncertainties.items():
            filtered = filter_value(value)
            if filtered is not None:
                uncertainties[key] = filtered
        return QHAResult(
            jobname=self.jobname,
            temperature=self.temperature,
            pressure=self.pressure,
            volume=self.volume,
            static_energy=self.static_energy,
            equilibrium_volume=filter_value(self.equilibrium_volume),
            zero_point_energy=filter_value(self.zero_point_energy),
            thermal_energy=filter_value(self.thermal_energy),
            internal_energy=filter_value(self.internal_energy),
            entropy=filter_value(self.entropy),
            vibrational_free_energy=filter_value(self.vibrational_free_energy),
            free_energy=filter_value(self.free_energy),
            isochoric_heat_capacity=filter_value(self.isochoric_heat_capacity),
            isobaric_heat_capacity=filter_value(self.isobaric_heat_capacity),
            heat_capacity_difference=filter_value(self.heat_capacity_difference),
            isothermal_bulk_modulus=filter_value(self.isothermal_bulk_modulus),
            adiabatic_bulk_modulus=filter_value(self.adiabatic_bulk_modulus),
            bulk_modulus_derivative=filter_value(self.bulk_modulus_derivative),
            thermal_expansion=filter_value(self.thermal_expansion),
            thermal_expansion_mixed=filter_value(self.thermal_expansion_mixed),
            thermal_expansion_mode=filter_value(self.thermal_expansion_mode),
            thermal_expansion_numerical=filter_value(self.thermal_expansion_numerical),
            thermal_expansion_source=filter_value(self.thermal_expansion_source),
            equilibrium_lattice=filter_value(self.equilibrium_lattice),
            lattice_parameters=filter_value(self.lattice_parameters),
            lattice_parameter_derivatives=filter_value(
                self.lattice_parameter_derivatives
            ),
            axial_thermal_expansion=filter_value(self.axial_thermal_expansion),
            thermal_expansion_tensor=filter_value(self.thermal_expansion_tensor),
            structural_extrapolation_mask=filter_value(
                self.structural_extrapolation_mask
            ),
            enthalpy=filter_value(self.enthalpy),
            gibbs_free_energy=filter_value(self.gibbs_free_energy),
            gruneisen=filter_value(self.gruneisen),
            mode_weighted_gruneisen=filter_value(self.mode_weighted_gruneisen),
            mode_gruneisen=self.mode_gruneisen,
            uncertainties=uncertainties,
            fit_records=list(self.fit_records),
            failed_points=list(self.failed_points),
            valid_mask=mask.copy(),
            completed=self.completed,
            metadata=dict(self.metadata),
        )


def diagnostics_to_dict(
    diagnostics: list[QHAFitRecord] | tuple[QHAFitRecord, ...],
) -> list[dict[str, Any]]:
    """Convert diagnostic records to dictionaries.

    Parameters
    ----------
    diagnostics : sequence of QHAFitRecord
        Diagnostic records to serialize.

    Returns
    -------
    list of dict
        Serialized diagnostic records.
    """
    return [record.as_dict() for record in diagnostics]


def result_metadata_from_options(options: QHAOptions) -> dict[str, Any]:
    """Return workflow metadata derived from QHA options.

    Parameters
    ----------
    options : QHAOptions
        Options controlling the QHA workflow.

    Returns
    -------
    dict
        Metadata describing the selected QHA scheme, minimization method, units,
        and diagnostic policies.
    """
    return {
        "scheme": options.scheme,
        "minimization": options.minimization,
        "eos": options.eos,
        "energy_degree": options.energy_degree,
        "free_energy_degree": options.free_energy_degree,
        "frequency_degree": options.frequency_degree,
        "structural_degree": options.structural_degree,
        "calculate_gruneisen": options.calculate_gruneisen,
        "calculate_mode_gruneisen": options.calculate_mode_gruneisen,
        "thermal_expansion_method": options.thermal_expansion_method,
        "gruneisen_allow_nonpositive": options.gruneisen_allow_nonpositive,
        "gruneisen_min_cv_fraction": options.gruneisen_min_cv_fraction,
        "units": {
            "energy": options.energy_unit,
            "volume": options.volume_unit,
            "frequency": options.frequency_unit,
            "temperature": options.temperature_unit,
            "pressure": options.pressure_unit,
            "entropy": f"{options.energy_unit} cell^-1 K^-1",
            "heat_capacity": f"{options.energy_unit} cell^-1 K^-1",
        },
        "diagnostics": {
            "debug": options.debug,
            "store_fit_diagnostics": options.store_fit_diagnostics,
            "estimate_uncertainties": options.estimate_uncertainties,
            "uncertainty_method": options.uncertainty_method,
            "uncertainty_relative_step": options.uncertainty_relative_step,
            "uncertainty_confidence_level": options.uncertainty_confidence_level,
            "uncertainty_samples": options.uncertainty_samples,
            "uncertainty_seed": options.uncertainty_seed,
            "uncertainty_minimum_success_fraction": (
                options.uncertainty_minimum_success_fraction
            ),
            "extrapolation_policy": options.extrapolation_policy,
            "fit_failure_policy": options.fit_failure_policy,
            "fit_quality_policy": options.fit_quality_policy,
            "max_consecutive_failures": options.max_consecutive_failures,
        },
    }


def merge_uncertainties(
    *items: Mapping[str, np.ndarray | float | int | list[float] | tuple[float, ...]],
) -> dict[str, np.ndarray]:
    """Merge uncertainty mappings into NumPy arrays.

    Parameters
    ----------
    *items : mapping
        One or more dictionaries containing uncertainty values.

    Returns
    -------
    dict
        Merged dictionary where every value is converted to ``ndarray``.
    """
    merged: dict[str, np.ndarray] = {}
    for item in items:
        for key, value in item.items():
            merged[key] = np.asarray(value, dtype=np.float64)
    return merged
