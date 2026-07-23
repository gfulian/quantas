# -*- coding: utf-8 -*-

"""User-configurable scientific policies for thermoelastic workflows."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Literal
import numpy as np
from quantas.core.math.fitting import FitMethod
from .types import (
    ThermoelasticAdiabaticMode,
    ThermoelasticExtrapolationPolicy,
    ThermoelasticFitFailurePolicy,
    ThermoelasticQualityPolicy,
    ThermoelasticReportLevel,
    ThermoelasticStabilityPolicy,
)


@dataclass(slots=True)
class ThermoelasticOptions:
    """Control quasi-static fitting and pressure-temperature reconstruction.

    Parameters
    ----------
    reference_eos : str, optional
        Integrated static energy EOS used to obtain fixed ``V0``, ``K0`` and
        ``Kprime``. ``BM3`` is the production default.
    finite_strain_order : {2, 3}, optional
        Eulerian finite-strain order used for ``Cij(V)``.
    fit_method : FitMethod, optional
        Regression method for independent elastic components. The current
        implementation accepts OLS because CRYSTAL SOEC outputs do not include
        observation uncertainties.
    max_iterations : int or None, optional
        Maximum model evaluations for each component fit.
    zero_tolerance : float, optional
        Components whose maximum absolute symmetry-averaged value does not
        exceed this tolerance, in GPa, are retained as exact zeros.
    relative_error_floor : float, optional
        Positive GPa floor used in relative-residual diagnostics.
    extrapolation_policy : {"fail", "warn", "allow"}, optional
        Policy for QHA volumes outside the elastic-volume interval.
    fit_failure_policy : {"stop", "continue", "raise"}, optional
        Response to failed independent-component fits. ``stop`` returns a
        diagnostic result without reconstructing tensors; ``continue`` uses
        NaN values for failed components; ``raise`` aborts the calculation.
    report_level : {"standard", "extended", "debug"}, optional
        Default frontend-neutral reporting detail.
    propagate_reference_eos_covariance : bool, optional
        Include the fitted static-EOS covariance in predicted uncertainties.
    propagate_volume_uncertainty : bool, optional
        Include QHA equilibrium-volume uncertainties when available.
    solver_debug : bool, optional
        Retain bounded numerical model-evaluation traces in every fit record.
    quality_policy : {"fail", "warn", "allow"}, optional
        Response to fit-support diagnostics. ``fail`` rejects unsupported
        component fits, ``warn`` records warnings without changing numerical
        results, and ``allow`` records diagnostics silently.
    minimum_fit_points : int, optional
        Minimum preferred number of observations for each two-parameter
        component fit.
    minimum_eulerian_strain_span : float, optional
        Minimum preferred sampled Eulerian finite-strain span.
    maximum_design_condition_number : float, optional
        Maximum preferred condition number of the column-normalized design
        matrix.
    maximum_relative_symmetry_spread : float, optional
        Maximum preferred disagreement among symmetry-equivalent observations,
        relative to the local component magnitude.
    maximum_leave_one_out_parameter_change : float, optional
        Maximum preferred scale-aware parameter change when omitting one
        elastic-volume point.
    maximum_order_parameter_change : float, optional
        Maximum preferred scale-aware parameter change between second- and
        third-order finite-strain descriptions.
    stability_policy : {"fail", "warn", "allow"}, optional
        Response when reconstructed Wallace stiffness matrices are unstable or
        indeterminate.
    stability_tolerance : float, optional
        Minimum accepted stiffness eigenvalue in GPa for the generic positive-
        definiteness criterion.
    adiabatic_mode : {"auto", "off", "require"}, optional
        Control conversion of reconstructed isothermal tensors. ``auto``
        converts when complete QHA ``C_V`` and anisotropic expansion data are
        available, ``off`` disables conversion, and ``require`` raises when
        those data are unavailable or invalid.
    propagate_adiabatic_uncertainty : bool, optional
        Propagate available isothermal-tensor, volume, heat-capacity, and
        thermal-expansion uncertainties to the adiabatic tensor using a
        first-order delta method.
    """

    reference_eos: str = "BM3"
    finite_strain_order: Literal[2, 3] = 3
    fit_method: FitMethod = FitMethod.OLS
    max_iterations: int | None = None
    zero_tolerance: float = 1.0e-10
    relative_error_floor: float = 1.0e-8
    extrapolation_policy: ThermoelasticExtrapolationPolicy = "warn"
    fit_failure_policy: ThermoelasticFitFailurePolicy = "stop"
    report_level: ThermoelasticReportLevel = "standard"
    propagate_reference_eos_covariance: bool = True
    propagate_volume_uncertainty: bool = True
    solver_debug: bool = False
    quality_policy: ThermoelasticQualityPolicy = "warn"
    minimum_fit_points: int = 4
    minimum_eulerian_strain_span: float = 5.0e-3
    maximum_design_condition_number: float = 1.0e6
    maximum_relative_symmetry_spread: float = 1.0e-2
    maximum_leave_one_out_parameter_change: float = 5.0e-1
    maximum_order_parameter_change: float = 5.0e-1
    stability_policy: ThermoelasticStabilityPolicy = "warn"
    stability_tolerance: float = 0.0
    adiabatic_mode: ThermoelasticAdiabaticMode = "auto"
    propagate_adiabatic_uncertainty: bool = True

    def __post_init__(self) -> None:
        """Normalize enums and validate numerical controls."""
        self.fit_method = FitMethod(self.fit_method)
        if self.fit_method is not FitMethod.OLS:
            raise ValueError("thermoelastic component fitting currently supports OLS")
        if self.finite_strain_order not in (2, 3):
            raise ValueError("finite_strain_order must be 2 or 3")
        if self.max_iterations is not None and self.max_iterations <= 0:
            raise ValueError("max_iterations must be positive when provided")
        if not np.isfinite(self.zero_tolerance) or self.zero_tolerance < 0.0:
            raise ValueError("zero_tolerance must be finite and non-negative")
        if (
            not np.isfinite(self.relative_error_floor)
            or self.relative_error_floor <= 0.0
        ):
            raise ValueError("relative_error_floor must be finite and positive")
        if self.extrapolation_policy not in {"fail", "warn", "allow"}:
            raise ValueError("invalid extrapolation_policy")
        if self.fit_failure_policy not in {"stop", "continue", "raise"}:
            raise ValueError("invalid fit_failure_policy")
        if self.report_level not in {"standard", "extended", "debug"}:
            raise ValueError("invalid report_level")
        if self.quality_policy not in {"fail", "warn", "allow"}:
            raise ValueError("invalid quality_policy")
        if self.minimum_fit_points < 3:
            raise ValueError("minimum_fit_points must be at least 3")
        for name in (
            "minimum_eulerian_strain_span",
            "maximum_design_condition_number",
            "maximum_relative_symmetry_spread",
            "maximum_leave_one_out_parameter_change",
            "maximum_order_parameter_change",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
            setattr(self, name, value)
        if self.maximum_design_condition_number <= 0.0:
            raise ValueError("maximum_design_condition_number must be positive")
        if self.stability_policy not in {"fail", "warn", "allow"}:
            raise ValueError("invalid stability_policy")
        if self.adiabatic_mode not in {"auto", "off", "require"}:
            raise ValueError("invalid adiabatic_mode")
        if not np.isfinite(self.stability_tolerance):
            raise ValueError("stability_tolerance must be finite")
        self.stability_tolerance = float(self.stability_tolerance)
        if self.report_level == "debug":
            self.solver_debug = True


__all__ = ["ThermoelasticOptions"]
