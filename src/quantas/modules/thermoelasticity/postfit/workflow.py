# -*- coding: utf-8 -*-

"""Post-fit pressure-temperature reconstruction from thermoelastic archives.

This module evaluates the already fitted quasi-static elastic model without
repeating either the static EOS fit or the independent-component fits.  The
QHA equilibrium-volume surface archived by the fitting run supplies
``V(P, T)`` through rectilinear interpolation or controlled extrapolation.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Sequence, cast

import numpy as np
from numpy.typing import ArrayLike

from quantas.models import ResultData
from quantas.modules.thermoelasticity.models import (
    ThermoelasticDepthProfile,
    ThermoelasticExtrapolationPolicy,
    ThermoelasticOptions,
    ThermoelasticResult,
)


from .evaluation import evaluate_thermoelastic_grid, evaluate_thermoelastic_profile
from .options import thermoelastic_options_from_mapping
from .policies import _append_stability_warning


def analyze_thermoelastic_result(
    result_data: ResultData,
    *,
    temperature: ArrayLike | None = None,
    pressure: ArrayLike | None = None,
    profiles: Sequence[ThermoelasticDepthProfile] = (),
    extrapolation_policy: str | None = None,
) -> ResultData:
    """Create a post-fit archive for a grid and optional depth profiles."""
    source = result_data.results.get("thermoelasticity")
    if not isinstance(source, ThermoelasticResult):
        raise ValueError("result does not contain a thermoelasticity payload")
    if not source.completed:
        raise ValueError("thermoelastic component fitting did not complete")
    options = thermoelastic_options_from_mapping(result_data.options)
    if extrapolation_policy is not None:
        options = replace(
            options,
            extrapolation_policy=cast(
                ThermoelasticExtrapolationPolicy,
                extrapolation_policy,
            ),
        )
    target_temperature = source.temperature if temperature is None else temperature
    target_pressure = source.pressure if pressure is None else pressure
    payload = evaluate_thermoelastic_grid(
        source,
        target_temperature,
        target_pressure,
        options=options,
    )
    profile_results = {
        profile.name: evaluate_thermoelastic_profile(source, profile, options=options)
        for profile in profiles
    }
    payload.profiles = profile_results
    payload.metadata["depth_profiles"] = list(profile_results)
    warnings = list(result_data.warnings)
    qha_count = int(
        np.count_nonzero(np.asarray(payload.qha_extrapolation_mask, dtype=np.bool_))
    )
    elastic_count = int(np.count_nonzero(payload.extrapolation_mask))
    if options.extrapolation_policy == "warn":
        if qha_count:
            warnings.append(
                f"{qha_count} requested pressure-temperature states lie outside "
                "the archived QHA coordinate grid"
            )
        if elastic_count:
            warnings.append(
                f"{elastic_count} requested pressure-temperature states lie outside "
                "the sampled elastic-volume interval"
            )
        for name, profile in profile_results.items():
            profile_count = int(
                np.count_nonzero(
                    profile.qha_extrapolation_mask | profile.elastic_extrapolation_mask
                )
            )
            if profile_count:
                warnings.append(
                    f"{profile_count} states in depth profile '{name}' require "
                    "QHA-coordinate or elastic-volume extrapolation"
                )
    if options.stability_policy == "warn":
        _append_stability_warning(
            warnings,
            payload.stability,
            "requested pressure-temperature grid",
        )
        for name, profile in profile_results.items():
            _append_stability_warning(
                warnings,
                profile.stability,
                f"depth profile '{name}'",
            )
    return ResultData(
        metadata=result_data.metadata,
        input_data=result_data.input_data,
        options=dict(result_data.options),
        results={"thermoelasticity": payload},
        warnings=warnings,
        events=list(result_data.events),
    )


def analyze_thermoelastic_profiles(
    result_data: ResultData,
    profiles: Sequence[ThermoelasticDepthProfile],
    *,
    options: ThermoelasticOptions | None = None,
    extrapolation_policy: str | None = None,
) -> ResultData:
    """Evaluate only geological profiles without reconstructing a P-T grid.

    Parameters
    ----------
    result_data : ResultData
        Completed thermoelastic fit archive.
    profiles : sequence of ThermoelasticDepthProfile
        One or more named depth-pressure-temperature paths.
    extrapolation_policy : {"fail", "warn", "allow"} or None, optional
        Override the policy stored in the fit archive.

    Returns
    -------
    ResultData
        Copy of the fit result retaining its archived QHA support fields and
        containing evaluated profile results, while leaving grid stiffness
        arrays absent.

    Raises
    ------
    ValueError
        If fitting did not complete or no profile was supplied.
    """
    import copy

    source = result_data.results.get("thermoelasticity")
    if not isinstance(source, ThermoelasticResult):
        raise ValueError("result does not contain a thermoelasticity payload")
    if not source.completed:
        raise ValueError("thermoelastic component fitting did not complete")
    if not profiles:
        raise ValueError("at least one depth profile is required")
    options = (
        thermoelastic_options_from_mapping(result_data.options)
        if options is None
        else options
    )
    if extrapolation_policy is not None:
        options = replace(
            options,
            extrapolation_policy=cast(
                ThermoelasticExtrapolationPolicy,
                extrapolation_policy,
            ),
        )
    payload = copy.deepcopy(source)
    payload.independent_stiffness = None
    payload.sigma_independent_stiffness = None
    payload.independent_stiffness_covariance = None
    payload.stiffness_isothermal = None
    payload.sigma_stiffness_isothermal = None
    payload.stiffness_adiabatic = None
    payload.sigma_stiffness_adiabatic = None
    payload.adiabatic_correction = None
    payload.adiabatic_thermal_stress = None
    payload.adiabatic_valid_mask = None
    payload.stability = None
    profile_results = {
        profile.name: evaluate_thermoelastic_profile(source, profile, options=options)
        for profile in profiles
    }
    payload.profiles = profile_results
    payload.metadata["depth_profiles"] = list(profile_results)
    payload.metadata["analysis_stage"] = "profiles_only"
    payload.metadata["qha_grid_reconstructed"] = False
    warnings = list(result_data.warnings)
    if options.extrapolation_policy == "warn":
        for name, profile in profile_results.items():
            count = int(
                np.count_nonzero(
                    profile.qha_extrapolation_mask | profile.elastic_extrapolation_mask
                )
            )
            if count:
                warnings.append(
                    f"{count} states in depth profile '{name}' require "
                    "QHA-coordinate or elastic-volume extrapolation"
                )
    if options.stability_policy == "warn":
        for name, profile in profile_results.items():
            _append_stability_warning(
                warnings,
                profile.stability,
                f"depth profile '{name}'",
            )
    return ResultData(
        metadata=result_data.metadata,
        input_data=result_data.input_data,
        options=dict(result_data.options),
        results={"thermoelasticity": payload},
        warnings=warnings,
        events=list(result_data.events),
    )


__all__ = ["analyze_thermoelastic_profiles", "analyze_thermoelastic_result"]
