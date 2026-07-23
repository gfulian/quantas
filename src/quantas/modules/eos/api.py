# -*- coding: utf-8 -*-

"""Frontend-neutral Python API for EOS fitting workflows."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from quantas.core.math.fitting import (
    BaseFitModel,
    EffectiveVarianceFitter,
    FitMethod,
    FitObservations,
    FitOptions,
    FitResult,
    FitSolver,
    LeastSquaresFitter,
    OrthogonalDistanceFitter,
    ParameterDefinition,
    ParameterMap,
)
from quantas.core.physics.eos import (
    EOSModel,
    PVTCouplingFamily,
    PVTModel,
    TemperatureEOSModel,
)

from .domains.pv import (
    AxialEOSFitModel,
    PressureEOSFitModel,
    axial_to_volume_parameters,
    build_axial_parameter_map,
    build_pressure_parameter_map,
)
from .models import (
    EOSDataset,
    EOSFitDomain,
    EOSFitRequest,
    EOSFitResult,
    EOSSeries,
)
from .domains.vt import (
    TemperatureEOSFitModel,
    build_temperature_parameter_map,
)
from .domains.pvt import (
    PVTEOSFitModel,
    build_pvt_parameter_map,
)


@dataclass(slots=True)
class _PreparedPVFit:
    """Validated numerical inputs for one volumetric or axial P-V fit."""

    series: EOSSeries
    observations: FitObservations
    selected: FitObservations
    parameter_map: ParameterMap
    model: BaseFitModel
    options: FitOptions
    axial: bool = False


@dataclass(slots=True)
class _PreparedPVTFit:
    """Validated numerical inputs for one global P--V--T fit."""

    observations: FitObservations
    selected: FitObservations
    parameter_map: ParameterMap
    model: PVTEOSFitModel
    options: FitOptions


@dataclass(slots=True)
class _PreparedVTFit:
    """Validated numerical inputs for one volumetric or axial V-T fit."""

    series: EOSSeries
    observations: FitObservations
    selected: FitObservations
    parameter_map: ParameterMap
    model: TemperatureEOSFitModel
    options: FitOptions
    axial: bool = False


class EOSFitter:
    """Execute EOS fit requests through shared numerical services.

    The class is the scientific authority used by Python callers, future batch
    workflows, the interactive CLI, and a future GUI. It contains no Click,
    Rich, HDF5, or rendering dependencies.
    """

    def __init__(self) -> None:
        self._least_squares = LeastSquaresFitter()
        self._effective_variance = EffectiveVarianceFitter(self._least_squares)
        self._orthogonal_distance = OrthogonalDistanceFitter()

    def validate_request(self, dataset: EOSDataset, request: EOSFitRequest) -> None:
        """Validate and prepare one request without running a numerical solver.

        Parameters
        ----------
        dataset : EOSDataset
            Normalized input data.
        request : EOSFitRequest
            Candidate scientific and numerical request.

        Raises
        ------
        ValueError, TypeError, NotImplementedError
            If the dataset, model, parameters, uncertainties, or solver
            contract cannot support the request.

        Notes
        -----
        This method is intended for dry-run validation by declarative
        frontends. It constructs the same observations and parameter map used
        by :meth:`fit`, but does not invoke SciPy or an ODR backend.
        """
        request = _request_with_dataset_default_selection(dataset, request)
        _validate_supported_request(request)
        if request.domain is EOSFitDomain.PRESSURE_VOLUME:
            prepared: Any = _prepare_pv_fit(dataset, request)
        elif request.domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            prepared = _prepare_pvt_fit(dataset, request)
        else:
            prepared = _prepare_vt_fit(dataset, request)
        selected = prepared.selected
        method = prepared.options.method
        if method is FitMethod.WLS:
            selected.require_positive_sigma("y")
        elif method is FitMethod.EFFECTIVE_VARIANCE:
            if selected.sigma_x is None and selected.sigma_y is None:
                raise ValueError(
                    "effective variance requires independent- or dependent-variable uncertainties"
                )
        elif method is FitMethod.ODR:
            selected.require_positive_sigma("x")
            selected.require_positive_sigma("y")

    def parameter_definitions(
        self,
        dataset: EOSDataset,
        request: EOSFitRequest,
    ) -> tuple[ParameterDefinition, ...]:
        """Return the resolved parameter contract for one fit request.

        Parameters
        ----------
        dataset : EOSDataset
            Normalized input data used to estimate unspecified parameters.
        request : EOSFitRequest
            Candidate model, target, constraints, and numerical options.

        Returns
        -------
        tuple of ParameterDefinition
            Complete reporting-order definitions, including free, fixed,
            implied, and derived parameters.

        Raises
        ------
        ValueError, TypeError, NotImplementedError
            If the dataset and request cannot define a valid fitting problem.

        Notes
        -----
        This method does not invoke a numerical solver. It is intended for
        session controllers and other frontends that need model-aware initial
        values and bounds without duplicating EOS parameter logic.
        """
        _validate_supported_request(request)
        if request.domain is EOSFitDomain.PRESSURE_VOLUME:
            prepared: Any = _prepare_pv_fit(dataset, request)
        elif request.domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            prepared = _prepare_pvt_fit(dataset, request)
        else:
            prepared = _prepare_vt_fit(dataset, request)
        return tuple(prepared.parameter_map.definitions)

    def fit(self, dataset: EOSDataset, request: EOSFitRequest) -> EOSFitResult:
        """Execute one frontend-neutral EOS fit.

        Parameters
        ----------
        dataset : EOSDataset
            Normalized input data.
        request : EOSFitRequest
            Model, target, constraints, mask, and statistical strategy.

        Returns
        -------
        EOSFitResult
            Structured fit, predictions, warnings, and provenance.

        Raises
        ------
        NotImplementedError
            If the requested scientific domain or fitting strategy is outside
            the implemented pressure-volume and volume-temperature milestones.
        ValueError
            If the selected dataset cannot support the requested fit.
        """
        request = _request_with_dataset_default_selection(dataset, request)
        _validate_supported_request(request)
        if request.domain is EOSFitDomain.PRESSURE_VOLUME:
            prepared_pv = _prepare_pv_fit(dataset, request)
            solver = self._solver_for(prepared_pv.options.method)
            fit = solver.fit_problem(
                prepared_pv.model,
                prepared_pv.observations,
                prepared_pv.parameter_map,
                prepared_pv.options,
            )
            fit = _augment_axial_fit(prepared_pv, fit)
            return _build_pv_result(dataset, request, prepared_pv, fit)

        if request.domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            prepared_pvt = _prepare_pvt_fit(dataset, request)
            solver = self._solver_for(prepared_pvt.options.method)
            fit = solver.fit_problem(
                prepared_pvt.model,
                prepared_pvt.observations,
                prepared_pvt.parameter_map,
                prepared_pvt.options,
            )
            return _build_pvt_result(dataset, request, prepared_pvt, fit)

        prepared_vt = _prepare_vt_fit(dataset, request)
        solver = self._solver_for(prepared_vt.options.method)
        fit = solver.fit_problem(
            prepared_vt.model,
            prepared_vt.observations,
            prepared_vt.parameter_map,
            prepared_vt.options,
        )
        fit = _augment_axial_vt_fit(prepared_vt, fit)
        return _build_vt_result(dataset, request, prepared_vt, fit)

    def _solver_for(self, method: FitMethod) -> FitSolver:
        """Return the general numerical service for one statistical method."""
        if method in {FitMethod.OLS, FitMethod.WLS}:
            return self._least_squares
        if method is FitMethod.EFFECTIVE_VARIANCE:
            return self._effective_variance
        if method is FitMethod.ODR:
            return self._orthogonal_distance
        raise NotImplementedError(f"unsupported EOS fitting method: {method.value}")


def _request_with_dataset_default_selection(
    dataset: EOSDataset, request: EOSFitRequest
) -> EOSFitRequest:
    """Apply non-destructive dataset exclusions to an unmasked API request."""
    if request.mask is not None or dataset.excluded_npoints == 0:
        return request
    metadata = dict(request.metadata)
    metadata.setdefault(
        "selection",
        {
            "base": "default",
            "total": dataset.npoints,
            "selected": dataset.selected_npoints,
            "excluded": dataset.excluded_npoints,
            "groups": list(dataset.group_summary()),
        },
    )
    return replace(
        request,
        mask=dataset.selection_mask(),
        metadata=metadata,
    )


def _validate_supported_request(request: EOSFitRequest) -> None:
    """Reject requests outside the implemented P-V and V-T milestones."""
    if request.domain not in {
        EOSFitDomain.PRESSURE_VOLUME,
        EOSFitDomain.VOLUME_TEMPERATURE,
        EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE,
    }:
        raise NotImplementedError(
            "the current EOS milestone implements P-V, V-T, and P-V-T fitting"
        )
    if request.target not in {"volume", "a", "b", "c"}:
        raise NotImplementedError(
            "EOS fitting supports volume and the linear cell parameters a, b, c"
        )
    if request.options.method not in {
        FitMethod.OLS,
        FitMethod.WLS,
        FitMethod.EFFECTIVE_VARIANCE,
        FitMethod.ODR,
    }:
        raise NotImplementedError(
            "the current EOS milestone implements OLS, WLS, effective variance, and ODR"
        )


def _prepare_pv_fit(
    dataset: EOSDataset,
    request: EOSFitRequest,
) -> _PreparedPVFit:
    """Build observations, parameters, model, and options for one P-V fit."""
    dataset.require_variable_coordinate(
        "pressure",
        purpose="a pressure-volume equation of state",
        mask=request.mask,
    )
    dataset.require_variable_coordinate(
        request.target,
        purpose="a pressure-volume equation of state",
        mask=request.mask,
    )
    if not isinstance(request.model, EOSModel):
        raise TypeError("P-V requests require a pressure EOS model")
    series = dataset.series(
        request.target,
        independent="pressure",
        mask=request.mask,
    )
    if request.target == "volume":
        observations = series.pressure_observations()
        selected = observations.selected()
        parameter_map = build_pressure_parameter_map(
            request.model,
            selected.x,
            selected.y,
            request.constraints,
            pressure_unit=series.units.get("pressure", "GPa"),
            volume_unit=series.units.get("volume", "angstrom^3"),
        )
        model: BaseFitModel = PressureEOSFitModel(request.model)
        axial = False
    else:
        observations = series.axial_pressure_observations()
        selected = observations.selected()
        selected_lengths = np.cbrt(selected.x)
        parameter_map = build_axial_parameter_map(
            request.model,
            selected_lengths,
            selected.y,
            request.constraints,
            pressure_unit=series.units.get("pressure", "GPa"),
            length_unit=series.units.get(request.target, "angstrom"),
        )
        model = AxialEOSFitModel(request.model, request.target)
        axial = True
    if selected.size <= parameter_map.n_free:
        raise ValueError(
            "P-V fit requires more selected observations than free parameters"
        )
    return _PreparedPVFit(
        series=series,
        observations=observations,
        selected=selected,
        parameter_map=parameter_map,
        model=model,
        options=_fit_options(dataset, request),
        axial=axial,
    )


def _prepare_vt_fit(
    dataset: EOSDataset,
    request: EOSFitRequest,
) -> _PreparedVTFit:
    """Build observations, parameters, model, and options for one V-T fit."""
    dataset.require_variable_coordinate(
        "temperature",
        purpose="a volume-temperature equation",
        mask=request.mask,
    )
    dataset.require_variable_coordinate(
        request.target,
        purpose="a volume-temperature equation",
        mask=request.mask,
    )
    if not isinstance(request.model, TemperatureEOSModel):
        raise TypeError("V-T requests require a temperature EOS model")
    series = dataset.series(
        request.target,
        independent="temperature",
        mask=request.mask,
    )
    axial = request.target in {"a", "b", "c"}
    observations = (
        series.axial_temperature_observations()
        if axial
        else series.temperature_observations()
    )
    selected = observations.selected()
    model = TemperatureEOSFitModel(
        request.model,
        request.target,
        axial=axial,
    )
    parameter_map = build_temperature_parameter_map(
        request.model,
        selected.x,
        selected.y,
        request.constraints,
        value_unit=(
            series.units.get(series.target, "1")
            if axial
            else observations.y_unit or "1"
        ),
        axial=axial,
    )
    if selected.size <= parameter_map.n_free:
        raise ValueError(
            "V-T fit requires more selected observations than free parameters"
        )
    return _PreparedVTFit(
        series=series,
        observations=observations,
        selected=selected,
        parameter_map=parameter_map,
        model=model,
        options=_fit_options(dataset, request),
        axial=axial,
    )


def _prepare_pvt_fit(
    dataset: EOSDataset,
    request: EOSFitRequest,
) -> _PreparedPVTFit:
    """Build observations, parameters, model, and options for P--V--T."""
    for coordinate in ("pressure", "volume", "temperature"):
        dataset.require_variable_coordinate(
            coordinate,
            purpose="a pressure-volume-temperature equation",
            mask=request.mask,
        )
    if not isinstance(request.model, PVTModel):
        raise TypeError("P-V-T requests require a compositional PVTModel")
    observations = dataset.pvt_observations(mask=request.mask)
    selected = observations.selected()
    volume = selected.x[0]
    temperature = selected.x[1]
    parameter_map = build_pvt_parameter_map(
        request.model,
        volume,
        temperature,
        selected.y,
        request.constraints,
        pressure_unit=observations.y_unit or "GPa",
        volume_unit=dataset.units.get("volume", "angstrom^3"),
    )
    if selected.size <= parameter_map.n_free:
        raise ValueError(
            "P-V-T fit requires more selected observations than free parameters"
        )
    return _PreparedPVTFit(
        observations=observations,
        selected=selected,
        parameter_map=parameter_map,
        model=PVTEOSFitModel(request.model),
        options=_fit_options(dataset, request),
    )


def _fit_options(dataset: EOSDataset, request: EOSFitRequest) -> FitOptions:
    """Attach workflow provenance to generic mathematical options."""
    return request.options.to_fit_options(
        {
            "domain": request.domain.value,
            "target": request.target,
            "request_id": request.request_id,
            "dataset_jobname": dataset.jobname,
            "dataset_source": (None if dataset.source is None else str(dataset.source)),
        }
    )


def _augment_axial_fit(prepared: _PreparedPVFit, fit: FitResult) -> FitResult:
    """Add length-space diagnostics while preserving solver coordinates."""
    if not prepared.axial or not fit.success or fit.parameters is None:
        return fit
    if not isinstance(prepared.model, AxialEOSFitModel):
        raise TypeError("axial fit preparation requires AxialEOSFitModel")

    selected_length = np.cbrt(prepared.selected.x)
    diagnostics = fit.diagnostics
    if diagnostics is not None:
        metadata = dict(diagnostics.metadata)
        metadata.update(
            {
                "solver_coordinate": f"{prepared.series.target}^3",
                "solver_coordinate_unit": prepared.observations.x_unit,
                "reported_length": prepared.series.target,
                "reported_length_unit": prepared.series.units.get(
                    prepared.series.target
                ),
                "selected_length": selected_length.tolist(),
                "derivative_length": prepared.model.derivative_length(
                    selected_length, fit.parameters
                ).tolist(),
                "linear_modulus": prepared.model.linear_modulus(
                    selected_length, fit.parameters
                ).tolist(),
            }
        )
        if diagnostics.x_corrections is not None:
            adjusted_q = np.asarray(
                fit.metadata.get(
                    "adjusted_x",
                    prepared.selected.x + diagnostics.x_corrections,
                ),
                dtype=np.float64,
            )
            adjusted_length = np.cbrt(adjusted_q)
            metadata["adjusted_length"] = adjusted_length.tolist()
            metadata["length_corrections"] = (
                adjusted_length - selected_length
            ).tolist()
            metadata["x_corrections_are_in"] = "cubed-length coordinates"
        diagnostics = replace(diagnostics, metadata=metadata)

    auxiliary = axial_to_volume_parameters(fit.parameters)
    return replace(
        fit,
        diagnostics=diagnostics,
        metadata={
            **fit.metadata,
            "linear_eos": True,
            "linear_target": prepared.series.target,
            "coordinate_transform": "q=x^3",
            "auxiliary_volume_parameters": auxiliary,
        },
    )


def _augment_axial_vt_fit(
    prepared: _PreparedVTFit,
    fit: FitResult,
) -> FitResult:
    """Add original-length diagnostics to an axial V-T fit."""
    if not prepared.axial or not fit.success or fit.parameters is None:
        return fit
    diagnostics = fit.diagnostics
    selected_length = np.cbrt(prepared.selected.y)
    if diagnostics is not None:
        metadata = dict(diagnostics.metadata)
        metadata.update(
            {
                "solver_response": f"{prepared.series.target}^3",
                "solver_response_unit": prepared.observations.y_unit,
                "reported_length": prepared.series.target,
                "reported_length_unit": prepared.series.units.get(
                    prepared.series.target
                ),
                "selected_length": selected_length.tolist(),
                "auxiliary_expansion_coefficient": prepared.model.expansion_coefficient(
                    prepared.selected.x, fit.parameters
                ).tolist(),
                "linear_expansion_coefficient": prepared.model.linear_expansion_coefficient(
                    prepared.selected.x, fit.parameters
                ).tolist(),
            }
        )
        if diagnostics.y_corrections is not None:
            adjusted_q = prepared.selected.y + diagnostics.y_corrections
            adjusted_length = np.cbrt(adjusted_q)
            metadata["adjusted_length"] = adjusted_length.tolist()
            metadata["length_corrections"] = (
                adjusted_length - selected_length
            ).tolist()
            metadata["y_corrections_are_in"] = "cubed-length coordinates"
        diagnostics = replace(diagnostics, metadata=metadata)
    return replace(
        fit,
        diagnostics=diagnostics,
        metadata={
            **fit.metadata,
            "linear_eos": True,
            "linear_target": prepared.series.target,
            "coordinate_transform": "q=x^3",
            "coefficient_space": "auxiliary_cubed_length",
        },
    )


def _build_pv_result(
    dataset: EOSDataset,
    request: EOSFitRequest,
    prepared: _PreparedPVFit,
    fit: FitResult,
) -> EOSFitResult:
    """Wrap a general fit in the EOS-domain result contract."""
    predictions = _pressure_predictions(prepared, fit)
    warnings = _fit_warnings(fit)
    metadata = _result_metadata(dataset, prepared)
    derived = _derived_parameters(prepared, fit)
    return EOSFitResult(
        request=request,
        fit=fit,
        predictions=predictions,
        derived=derived,
        warnings=warnings,
        metadata=metadata,
    )


def _pressure_predictions(
    prepared: _PreparedPVFit,
    fit: FitResult,
) -> dict[str, np.ndarray]:
    """Calculate pressure over every original observation after success."""
    if not fit.success or fit.parameters is None:
        return {}
    coordinate = np.asarray(prepared.series.y, dtype=np.float64)
    if prepared.axial:
        coordinate = coordinate**3
    return {"pressure": prepared.model.evaluate(coordinate, fit.parameters)}


def _derived_parameters(
    prepared: _PreparedPVFit,
    fit: FitResult,
) -> dict[str, float]:
    """Return auxiliary volumetric values for a successful linear fit."""
    if not prepared.axial or not fit.success or fit.parameters is None:
        return {}
    auxiliary = axial_to_volume_parameters(fit.parameters)
    return {f"auxiliary_{name}": value for name, value in auxiliary.items()}


def _build_vt_result(
    dataset: EOSDataset,
    request: EOSFitRequest,
    prepared: _PreparedVTFit,
    fit: FitResult,
) -> EOSFitResult:
    """Wrap one general V-T fit in the EOS-domain result contract."""
    predictions = _temperature_predictions(prepared, fit)
    warnings = _fit_warnings(fit)
    metadata = _vt_result_metadata(dataset, prepared)
    derived = _vt_derived_parameters(prepared, fit)
    return EOSFitResult(
        request=request,
        fit=fit,
        predictions=predictions,
        derived=derived,
        warnings=warnings,
        metadata=metadata,
    )


def _temperature_predictions(
    prepared: _PreparedVTFit,
    fit: FitResult,
) -> dict[str, np.ndarray]:
    """Calculate fitted structure and expansion over all original temperatures."""
    if not fit.success or fit.parameters is None:
        return {}
    temperature = np.asarray(prepared.series.x, dtype=np.float64)
    solver_value = prepared.model.evaluate(temperature, fit.parameters)
    auxiliary_alpha = prepared.model.expansion_coefficient(temperature, fit.parameters)
    if not prepared.axial:
        return {
            prepared.series.target: solver_value,
            "expansion_coefficient": auxiliary_alpha,
            "temperature_derivative": prepared.model.derivative_x(
                temperature, fit.parameters
            ),
        }
    length = np.cbrt(solver_value)
    return {
        prepared.series.target: length,
        f"{prepared.series.target}^3": solver_value,
        "auxiliary_expansion_coefficient": auxiliary_alpha,
        "linear_expansion_coefficient": auxiliary_alpha / 3.0,
        "temperature_derivative_auxiliary": prepared.model.derivative_x(
            temperature, fit.parameters
        ),
    }


def _vt_derived_parameters(
    prepared: _PreparedVTFit,
    fit: FitResult,
) -> dict[str, float]:
    """Return exact expansion values at the fixed reference condition."""
    if not fit.success or fit.parameters is None:
        return {}
    mapping = {
        name: float(value)
        for name, value in zip(fit.parameter_names, fit.parameters, strict=True)
    }
    tref = float(mapping["temperature_ref"])
    alpha = float(
        prepared.model.expansion_coefficient(
            np.asarray([tref], dtype=np.float64), fit.parameters
        )[0]
    )
    derived = {
        "expansion_coefficient_at_reference": alpha,
        "reference_temperature": tref,
    }
    if prepared.axial:
        derived["length_ref"] = float(mapping["L0"])
        derived["linear_expansion_coefficient_at_reference"] = alpha / 3.0
    return derived


def _vt_result_metadata(
    dataset: EOSDataset,
    prepared: _PreparedVTFit,
) -> dict[str, Any]:
    """Return scientific metadata for one V-T result."""
    series = prepared.series
    selected = prepared.selected
    mask = np.asarray(series.mask, dtype=np.bool_)
    classification = dataset.classify(mask=series.mask)
    metadata: dict[str, Any] = {
        "relationship": f"{series.target}(temperature)",
        "residual_definition": (
            f"observed_{series.target}^3-calculated_{series.target}^3"
            if prepared.axial
            else f"observed_{series.target}-calculated_{series.target}"
        ),
        "solver_response": (f"{series.target}^3" if prepared.axial else series.target),
        "temperature_unit": series.units.get("temperature", "K"),
        "target_unit": series.units.get(series.target, "1"),
        "selected_mask": mask.tolist(),
        "sampled_temperature_range": [
            float(np.min(selected.x)),
            float(np.max(selected.x)),
        ],
        "parameter_map": prepared.parameter_map.as_dict(),
        "dataset_classification": classification.as_dict(),
        "reference_pressure": classification.reference_pressure,
        "is_isobaric": classification.is_isobaric,
        "input_target_scale": series.metadata.get("target_scale", "absolute"),
        "input_metadata": dict(dataset.metadata),
        "linear_eos": prepared.axial,
    }
    if not prepared.axial:
        metadata["sampled_target_range"] = [
            float(np.min(selected.y)),
            float(np.max(selected.y)),
        ]
        return metadata
    selected_length = np.cbrt(selected.y)
    metadata.update(
        {
            "linear_target": series.target,
            "length_unit": series.units.get(series.target, "1"),
            "coordinate_transform": "q=x^3",
            "uncertainty_transform": "sigma_q=3*x^2*sigma_x",
            "coefficient_space": "auxiliary_cubed_length",
            "sampled_length_range": [
                float(np.min(selected_length)),
                float(np.max(selected_length)),
            ],
            "sampled_cubed_length_range": [
                float(np.min(selected.y)),
                float(np.max(selected.y)),
            ],
        }
    )
    return metadata


def _build_pvt_result(
    dataset: EOSDataset,
    request: EOSFitRequest,
    prepared: _PreparedPVTFit,
    fit: FitResult,
) -> EOSFitResult:
    """Wrap a global P--V--T fit in the EOS-domain result contract."""
    predictions: dict[str, np.ndarray] = {}
    derived: dict[str, float] = {}
    if fit.success and fit.parameters is not None:
        predictions["pressure"] = prepared.model.evaluate(
            prepared.observations.x, fit.parameters
        )
        pressure_parameters, thermal_parameters, coupling_parameters = (
            prepared.model.split_parameters(fit.parameters)
        )
        core = prepared.model.core
        volume = prepared.observations.x[0]
        temperature = prepared.observations.x[1]
        predictions["bulk_modulus"] = core.bulk_modulus(
            prepared.model.pvt_model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            volume,
            temperature,
        )
        predictions["reference_volume"] = core.reference_volume(
            prepared.model.pvt_model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            temperature,
        )
        predictions["zero_pressure_bulk_modulus"] = core.zero_pressure_bulk_modulus(
            prepared.model.pvt_model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            temperature,
        )
        predictions["zero_pressure_expansion_coefficient"] = (
            core.expansion_coefficient_zero_pressure(
                prepared.model.pvt_model,
                pressure_parameters,
                thermal_parameters,
                coupling_parameters,
                temperature,
            )
        )
        predictions["zero_pressure_dK0_dT"] = (
            core.zero_pressure_bulk_modulus_temperature_derivative(
                prepared.model.pvt_model,
                pressure_parameters,
                thermal_parameters,
                coupling_parameters,
                temperature,
            )
        )
        if (
            prepared.model.pvt_model.coupling_family
            is PVTCouplingFamily.THERMAL_PRESSURE
        ):
            predictions["thermal_pressure"] = core.thermal_pressure_contribution(
                prepared.model.pvt_model,
                pressure_parameters,
                coupling_parameters,
                volume,
                temperature,
            )
        tref = float(
            dict(zip(fit.parameter_names, fit.parameters, strict=True))[
                "temperature_ref"
            ]
        )
        derived = {
            "reference_temperature": tref,
            "zero_pressure_expansion_at_reference": float(
                core.expansion_coefficient_zero_pressure(
                    prepared.model.pvt_model,
                    pressure_parameters,
                    thermal_parameters,
                    coupling_parameters,
                    np.asarray([tref], dtype=np.float64),
                )[0]
            ),
            "dK0_dT_at_reference": float(
                core.zero_pressure_bulk_modulus_temperature_derivative(
                    prepared.model.pvt_model,
                    pressure_parameters,
                    thermal_parameters,
                    coupling_parameters,
                    np.asarray([tref], dtype=np.float64),
                )[0]
            ),
        }
    selected = prepared.selected
    classification = dataset.classify(mask=request.mask)
    controlled_parameters = [
        definition.name
        for definition in prepared.parameter_map.definitions
        if definition.state.value == "fixed"
        and definition.metadata.get("source") == "user_fixed"
        and definition.metadata.get("default_state") == "free"
    ]
    metadata = {
        "relationship": "pressure(volume,temperature)",
        "residual_definition": "observed_pressure-calculated_pressure",
        "coordinate_order": ["volume", "temperature"],
        "pressure_unit": dataset.units.get("pressure", "GPa"),
        "volume_unit": dataset.units.get("volume", "angstrom^3"),
        "temperature_unit": dataset.units.get("temperature", "K"),
        "sampled_pressure_range": [
            float(np.min(selected.y)),
            float(np.max(selected.y)),
        ],
        "sampled_volume_range": [
            float(np.min(selected.x[0])),
            float(np.max(selected.x[0])),
        ],
        "sampled_temperature_range": [
            float(np.min(selected.x[1])),
            float(np.max(selected.x[1])),
        ],
        "parameter_map": prepared.parameter_map.as_dict(),
        "pvt_model": prepared.model.pvt_model.as_dict(),
        "dataset_classification": classification.as_dict(),
        "fit_mode": "controlled" if controlled_parameters else "global",
        "controlled_fit": bool(controlled_parameters),
        "controlled_parameters": controlled_parameters,
        "input_metadata": dict(dataset.metadata),
    }
    return EOSFitResult(
        request=request,
        fit=fit,
        predictions=predictions,
        derived=derived,
        warnings=_fit_warnings(fit),
        metadata=metadata,
    )


def _fit_warnings(fit: FitResult) -> list[str]:
    """Return backend warnings plus stable EOS-specific bound warnings."""
    messages = list(fit.warnings)
    diagnostics = fit.diagnostics
    if diagnostics is None or not diagnostics.parameters_at_bounds:
        return messages
    at_bounds = [
        name
        for name, flag in zip(
            fit.parameter_names,
            diagnostics.parameters_at_bounds,
            strict=True,
        )
        if flag
    ]
    if at_bounds:
        messages.append("parameters at bounds: " + ", ".join(at_bounds))
    return messages


def _result_metadata(
    dataset: EOSDataset,
    prepared: _PreparedPVFit,
) -> dict[str, Any]:
    """Return serialization-ready scientific and sampling metadata."""
    series = prepared.series
    selected = prepared.selected
    mask = np.asarray(series.mask, dtype=np.bool_)
    common: dict[str, Any] = {
        "relationship": f"pressure({series.target})",
        "residual_definition": "observed_pressure-calculated_pressure",
        "pressure_unit": series.units.get("pressure", "GPa"),
        "selected_mask": mask.tolist(),
        "sampled_pressure_range": [
            float(np.min(selected.y)),
            float(np.max(selected.y)),
        ],
        "parameter_map": prepared.parameter_map.as_dict(),
        "dataset_classification": dataset.classify(mask=series.mask).as_dict(),
        "input_metadata": dict(dataset.metadata),
    }
    if not prepared.axial:
        common.update(
            {
                "volume_unit": series.units.get("volume", "angstrom^3"),
                "sampled_volume_range": [
                    float(np.min(selected.x)),
                    float(np.max(selected.x)),
                ],
                "linear_eos": False,
            }
        )
        return common

    selected_length = np.cbrt(selected.x)
    common.update(
        {
            "linear_eos": True,
            "linear_target": series.target,
            "length_unit": series.units.get(series.target, "angstrom"),
            "coordinate_transform": "q=x^3",
            "uncertainty_transform": "sigma_q=3*x^2*sigma_x",
            "sampled_length_range": [
                float(np.min(selected_length)),
                float(np.max(selected_length)),
            ],
            "sampled_cubed_length_range": [
                float(np.min(selected.x)),
                float(np.max(selected.x)),
            ],
            "solver_coordinate_unit": prepared.observations.x_unit,
            "linear_parameter_relation": {
                "M0": "3*K0_auxiliary",
                "MP": "3*KP_auxiliary",
                "MPP": "3*KPP_auxiliary",
                "L0": "cuberoot(V0_auxiliary)",
            },
        }
    )
    return common


__all__ = ["EOSFitter"]
