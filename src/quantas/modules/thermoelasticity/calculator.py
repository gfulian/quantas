# -*- coding: utf-8 -*-

"""Calibration calculator for quasi-static thermoelasticity.

``run`` has one meaning in every frontend: fit the authoritative static EOS
and the independent volume-dependent elastic components.  Pressure,
temperature, and depth evaluation belongs exclusively to the analysis engine.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from quantas.core.events import EventLevel, Observer
from quantas.core.physics.units import convert_volume
from quantas.references import module_citation_keys
from quantas.models import BasicCalculator, InputData, ResultData, ResultMetadata
from quantas.modules.thermoelasticity.adiabatic import (
    StandardAdiabaticInputs,
    standard_adiabatic_inputs_from_qha,
)
from quantas.modules.thermoelasticity.analysis import (
    density_from_volume,
    independent_component_labels,
    normalized_cell_mass,
    qha_grid_in_standard_units,
)
from quantas.modules.thermoelasticity.fitting import (
    fit_elastic_components,
    fit_reference_static_eos,
)
from quantas.modules.thermoelasticity.models import (
    ElasticComponentFit,
    ReferenceEOSFit,
    ThermoelasticContext,
    ThermoelasticOptions,
    ThermoelasticResult,
)


@dataclass(slots=True)
class _CalibrationSourceFields:
    """Normalized QHA fields retained by the calibration archive."""

    temperature: np.ndarray
    pressure: np.ndarray
    volume: np.ndarray
    density: np.ndarray
    sigma_volume: np.ndarray | None
    outside_elastic_support: np.ndarray
    outside_count: int
    mass_kg: float
    mass_relative_spread: float
    adiabatic: StandardAdiabaticInputs | None


class ThermoelasticityCalculator(BasicCalculator):
    """Calibrate a reusable quasi-static thermoelastic model.

    Parameters
    ----------
    context : ThermoelasticContext
        Validated pairing of a CRYSTAL SOEC volume series and QHA result.
    options : ThermoelasticOptions or None, optional
        Fit, diagnostics, and uncertainty controls.
    observer : Observer or None, optional
        Quantas event observer.

    Notes
    -----
    The calculator never evaluates ``Cij(P,T)``.  It stores the source QHA
    fields required by :class:`ThermoelasticAnalysisEngine`, so CLI, Python,
    and future GUI frontends all share one analysis path.
    """

    module = "thermoelasticity"
    method = "quasistatic"

    def __init__(
        self,
        context: ThermoelasticContext,
        options: ThermoelasticOptions | None = None,
        observer: Observer | None = None,
    ) -> None:
        """Initialize the calibration calculator."""
        self.context = context
        self.thermoelastic_options = (
            ThermoelasticOptions() if options is None else options
        )
        input_data = InputData(
            source=context.input_data.source,
            data={
                "jobname": context.input_data.jobname,
                "method": context.input_data.method,
                "elastic_symmetry": context.input_data.elastic_series.elastic_symmetry,
                "elastic_points": len(context.input_data.elastic_series.points),
                "elastic_volume_bounds_A3": list(
                    context.input_data.elastic_series.volume_bounds
                ),
                "qha_source": (
                    None
                    if context.qha_result_data.input_data is None
                    else str(context.qha_result_data.input_data.source)
                ),
                "context": dict(context.metadata),
            },
        )
        super().__init__(
            input_data=input_data,
            options=asdict(self.thermoelastic_options),
            observer=observer,
        )

    def prepare(self) -> None:
        """Validate the static and QHA fields required by calibration."""
        super().prepare()
        if not self.context.has_complete_quasistatic_inputs:
            missing = ", ".join(self.context.missing_qha_fields)
            raise ValueError(
                f"QHA data are incomplete for quasi-static thermoelasticity: {missing}"
            )
        qha = self.context.qha
        if qha.temperature is None or qha.pressure is None:
            raise ValueError("QHA temperature and pressure grids are required")
        if qha.equilibrium_volume is None:
            raise ValueError("QHA equilibrium volumes are required")
        if qha.volume is None or qha.static_energy is None:
            raise ValueError("QHA sampled static volume-energy data are required")

    def run(self) -> ResultData:
        """Fit the reference EOS and independent elastic components.

        Returns
        -------
        ResultData
            Calibration archive containing fit diagnostics and source QHA
            fields, but no reconstructed stiffness grid.
        """
        qha_options = dict(self.context.qha_result_data.options)
        reference_eos = self._fit_reference_eos(qha_options)
        component_fits, failed_labels = self._fit_components(reference_eos)
        source = self._prepare_source_fields(qha_options)
        labels = independent_component_labels(
            self.context.input_data.elastic_series.elastic_symmetry
        )
        result = self._build_calibration_result(
            reference_eos=reference_eos,
            component_fits=component_fits,
            failed_labels=failed_labels,
            labels=labels,
            source=source,
        )
        return ResultData(
            metadata=ResultMetadata(module=self.module, method=self.method),
            input_data=self.input_data,
            options=dict(self.options),
            results={"thermoelasticity": result},
            warnings=list(self.warnings),
        )

    def _fit_reference_eos(self, qha_options: dict[str, Any]) -> ReferenceEOSFit:
        """Fit the authoritative QHA static energy-volume field."""
        qha = self.context.qha
        assert qha.volume is not None
        assert qha.static_energy is not None
        volume_unit = str(qha_options.get("volume_unit", "A"))
        energy_unit = str(qha_options.get("energy_unit", "Ha"))
        sampled_volume = np.asarray(
            convert_volume(qha.volume, volume_unit, "A"), dtype=np.float64
        )
        static_energy = np.asarray(qha.static_energy, dtype=np.float64)
        if sampled_volume.ndim != 1 or static_energy.shape != sampled_volume.shape:
            raise ValueError("QHA sampled volume and static energy must be aligned")
        self.emit("Fitting the authoritative static energy-volume EOS")
        reference = fit_reference_static_eos(
            sampled_volume,
            static_energy,
            eos=self.thermoelastic_options.reference_eos,
            energy_unit=energy_unit,
            volume_unit="A",
            max_iterations=self.thermoelastic_options.max_iterations,
        )
        self.emit(
            "Static reference EOS fitted",
            level=EventLevel.RESULT,
            data={
                "eos": reference.eos,
                "V0_A3": reference.reference_volume,
                "K0_GPa": reference.bulk_modulus,
                "Kprime": reference.bulk_modulus_derivative,
                "fit": reference.fit.as_dict(),
            },
        )
        return reference

    def _fit_components(
        self,
        reference_eos: ReferenceEOSFit,
    ) -> tuple[dict[str, ElasticComponentFit], tuple[str, ...]]:
        """Fit independent components and apply support/failure policies."""
        options = self.thermoelastic_options
        self.emit("Fitting independent elastic components")
        records = fit_elastic_components(
            self.context.input_data.elastic_series,
            reference_eos,
            options,
        )
        self._apply_fit_quality_policy(records)
        active = [record for record in records.values() if record.active]
        for index, record in enumerate(active, start=1):
            self.emit(
                f"Elastic component {record.label} fit processed",
                level=EventLevel.PROGRESS,
                progress=index / max(len(active), 1),
                data={"component": record.label},
            )
            if record.fit is not None:
                self.emit(
                    f"Elastic component {record.label}: {record.fit.status.value}",
                    level=EventLevel.RESULT,
                    data={
                        "kind": "thermoelastic_component_fit",
                        "component": record.label,
                        "fit": record.fit.as_dict(),
                    },
                )
        failed = tuple(
            label
            for label, record in records.items()
            if record.active and (record.fit is None or not record.fit.success)
        )
        if failed and options.fit_failure_policy == "raise":
            raise RuntimeError("elastic-component fits failed: " + ", ".join(failed))
        if failed and options.fit_failure_policy == "stop":
            self.add_warning(
                "calibration contains failed component fits: " + ", ".join(failed)
            )
        return records, failed

    def _prepare_source_fields(
        self,
        qha_options: dict[str, Any],
    ) -> _CalibrationSourceFields:
        """Normalize and validate source QHA fields archived for analysis."""
        qha = self.context.qha
        assert qha.temperature is not None
        assert qha.pressure is not None
        assert qha.equilibrium_volume is not None
        temperature_unit = str(qha_options.get("temperature_unit", "K"))
        pressure_unit = str(qha_options.get("pressure_unit", "GPa"))
        volume_unit = str(qha_options.get("volume_unit", "A"))
        temperature, pressure, volume = qha_grid_in_standard_units(
            qha.temperature,
            qha.pressure,
            qha.equilibrium_volume,
            temperature_unit=temperature_unit,
            pressure_unit=pressure_unit,
            volume_unit=volume_unit,
        )
        sigma_volume = self._qha_volume_uncertainty(volume_unit)
        series = self.context.input_data.elastic_series
        mass_kg, mass_relative_spread = normalized_cell_mass(series)
        if mass_relative_spread > 1.0e-5:
            self.add_warning(
                "normalized cell mass inferred from CRYSTAL density-volume "
                f"pairs varies by {100.0 * mass_relative_spread:.4g}%"
            )
        lower, upper = series.volume_bounds
        outside = np.asarray(
            (~np.isfinite(volume)) | (volume < lower) | (volume > upper),
            dtype=np.bool_,
        )
        outside_count = int(np.count_nonzero(outside))
        if outside_count:
            self.add_warning(
                f"{outside_count} archived QHA states lie outside the sampled "
                "elastic-volume interval; analysis policies will control their "
                "later evaluation"
            )
        adiabatic = standard_adiabatic_inputs_from_qha(
            qha,
            qha_options,
            volume.shape,
        )
        if self.thermoelastic_options.adiabatic_mode == "require" and adiabatic is None:
            raise ValueError(
                "adiabatic_mode='require' needs QHA isochoric heat capacity "
                "and thermal-expansion tensor fields"
            )
        return _CalibrationSourceFields(
            temperature=temperature,
            pressure=pressure,
            volume=volume,
            density=density_from_volume(volume, mass_kg),
            sigma_volume=sigma_volume,
            outside_elastic_support=outside,
            outside_count=outside_count,
            mass_kg=mass_kg,
            mass_relative_spread=mass_relative_spread,
            adiabatic=adiabatic,
        )

    def _build_calibration_result(
        self,
        *,
        reference_eos: ReferenceEOSFit,
        component_fits: dict[str, ElasticComponentFit],
        failed_labels: tuple[str, ...],
        labels: tuple[str, ...],
        source: _CalibrationSourceFields,
    ) -> ThermoelasticResult:
        """Assemble the passive calibration payload."""
        series = self.context.input_data.elastic_series
        lower, upper = series.volume_bounds
        adiabatic = source.adiabatic
        return ThermoelasticResult(
            jobname=self.context.input_data.jobname,
            reference_eos=reference_eos,
            component_fits=component_fits,
            independent_labels=labels,
            temperature=source.temperature,
            pressure=source.pressure,
            equilibrium_volume=source.volume,
            density=source.density,
            independent_stiffness=None,
            sigma_independent_stiffness=None,
            independent_stiffness_covariance=None,
            stiffness_isothermal=None,
            sigma_stiffness_isothermal=None,
            extrapolation_mask=source.outside_elastic_support,
            sigma_equilibrium_volume=source.sigma_volume,
            qha_extrapolation_mask=np.zeros_like(
                source.outside_elastic_support, dtype=np.bool_
            ),
            profiles={},
            isochoric_heat_capacity_cell=(
                None if adiabatic is None else adiabatic.heat_capacity
            ),
            sigma_isochoric_heat_capacity_cell=(
                None if adiabatic is None else adiabatic.sigma_heat_capacity
            ),
            thermal_expansion_tensor=(
                None if adiabatic is None else adiabatic.thermal_expansion_tensor
            ),
            sigma_thermal_expansion_tensor=(
                None if adiabatic is None else adiabatic.sigma_thermal_expansion_tensor
            ),
            stiffness_adiabatic=None,
            sigma_stiffness_adiabatic=None,
            adiabatic_correction=None,
            adiabatic_thermal_stress=None,
            adiabatic_valid_mask=None,
            stability=None,
            completed=not failed_labels,
            metadata=self._calibration_metadata(
                component_fits=component_fits,
                failed_labels=failed_labels,
                source=source,
                lower=lower,
                upper=upper,
                elastic_symmetry=series.elastic_symmetry,
                labels=labels,
                frame_normalization=dict(
                    series.metadata.get("frame_normalization", {})
                ),
            ),
        )

    def _calibration_metadata(
        self,
        *,
        component_fits: dict[str, ElasticComponentFit],
        failed_labels: tuple[str, ...],
        source: _CalibrationSourceFields,
        lower: float,
        upper: float,
        elastic_symmetry: str,
        labels: tuple[str, ...],
        frame_normalization: dict[str, Any],
    ) -> dict[str, Any]:
        """Return calibration provenance and scientific-policy metadata."""
        options = self.thermoelastic_options
        return {
            "analysis_stage": "calibration",
            "approximation": "quasi-static",
            "tensor_thermodynamic_condition": "not_evaluated",
            "adiabatic_mode": options.adiabatic_mode,
            "adiabatic_inputs_available": bool(source.adiabatic is not None),
            "reference_eos_source": (
                "QHA sampled static DFT energy-volume data; no thermal or "
                "zero-point contribution"
            ),
            "reference_eos_state": "static 0 K, P=0 reference",
            "wallace_convention": (
                "CRYSTAL PRESSURE-corrected Wallace stress-strain coefficients; "
                "no second pressure correction is applied"
            ),
            "citation_keys": list(module_citation_keys("thermoelasticity")),
            "elastic_symmetry": elastic_symmetry,
            "independent_component_labels": list(labels),
            "failed_component_fits": list(failed_labels),
            "fit_failure_policy": options.fit_failure_policy,
            "elastic_volume_min_A3": lower,
            "elastic_volume_max_A3": upper,
            "source_qha_states_outside_elastic_support": source.outside_count,
            "qha_grid_reconstructed": False,
            "fit_quality_policy": options.quality_policy,
            "fit_quality_counts": self._fit_quality_counts(component_fits),
            "normalized_cell_mass_kg": source.mass_kg,
            "normalized_cell_mass_relative_spread": source.mass_relative_spread,
            "frame_normalization": frame_normalization,
            "uncertainty_propagation": {
                "component_parameter_covariance": True,
                "shared_reference_eos_covariance": (
                    options.propagate_reference_eos_covariance
                ),
                "two_stage_eos_parameter_sensitivity": (
                    "symmetric finite differences of exact conditional component refits"
                ),
                "qha_volume_uncertainty_archived": bool(
                    source.sigma_volume is not None
                ),
            },
        }

    def _qha_volume_uncertainty(self, volume_unit: str) -> np.ndarray | None:
        """Return QHA volume uncertainty converted to angstrom cubed."""
        if not self.thermoelastic_options.propagate_volume_uncertainty:
            return None
        uncertainties = getattr(self.context.qha, "uncertainties", None)
        if not isinstance(uncertainties, dict):
            return None
        value = next(
            (
                uncertainties[key]
                for key in ("sigma_VT", "sigma_equilibrium_volume", "VT")
                if key in uncertainties
            ),
            None,
        )
        if value is None:
            return None
        sigma = np.asarray(convert_volume(value, volume_unit, "A"), dtype=np.float64)
        assert self.context.qha.equilibrium_volume is not None
        if sigma.shape != np.asarray(self.context.qha.equilibrium_volume).shape:
            self.add_warning(
                "QHA equilibrium-volume uncertainty has an incompatible shape "
                "and was excluded from propagation"
            )
            return None
        if np.any(~np.isfinite(sigma)) or np.any(sigma < 0.0):
            self.add_warning(
                "QHA equilibrium-volume uncertainty contains invalid values "
                "and was excluded from propagation"
            )
            return None
        return sigma

    def _apply_fit_quality_policy(self, component_fits: dict[str, Any]) -> None:
        """Apply the configured response to fit-support diagnostics."""
        unsupported = [
            label
            for label, record in component_fits.items()
            if record.quality is not None and record.quality.level == "unsupported"
        ]
        caution = [
            label
            for label, record in component_fits.items()
            if record.quality is not None and record.quality.level == "caution"
        ]
        if unsupported and self.thermoelastic_options.quality_policy == "fail":
            raise ValueError(
                "scientifically unsupported elastic-component fits: "
                + ", ".join(unsupported)
            )
        if self.thermoelastic_options.quality_policy == "warn":
            if unsupported:
                self.add_warning(
                    "scientifically unsupported elastic-component fits were "
                    "retained without numerical modification: " + ", ".join(unsupported)
                )
            if caution:
                self.add_warning(
                    "elastic-component fits with scientific-support cautions: "
                    + ", ".join(caution)
                )

    @staticmethod
    def _fit_quality_counts(component_fits: dict[str, Any]) -> dict[str, int]:
        """Return counts by stable fit-quality classification."""
        counts = {
            "supported": 0,
            "caution": 0,
            "unsupported": 0,
            "not_applicable": 0,
        }
        for record in component_fits.values():
            if record.quality is not None:
                counts[record.quality.level] += 1
        return counts


__all__ = ["ThermoelasticityCalculator"]
