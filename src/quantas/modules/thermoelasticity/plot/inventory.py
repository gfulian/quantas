# -*- coding: utf-8 -*-

"""Result-aware public plot inventory for thermoelastic workflows.

The thermoelastic workflow accumulates scientific capabilities. A calibration
archive can contain elastic-volume fits only, a point or grid analysis adds
reconstructed tensors, and profile analyses add one or more depth-dependent
paths. Discovery therefore reports the plot families that can actually be
built from the supplied result instead of assuming one fixed workflow stage.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from numpy.typing import NDArray

from quantas.models import (
    PlotContextDescriptor,
    PlotInventory,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
    ResultData,
)
from quantas.modules.thermoelasticity.models import ThermoelasticResult
from quantas.modules.thermoelasticity.plot.common import extract_thermoelastic_result
from quantas.modules.thermoelasticity.plot.components import resolve_components

_COMPONENT_GROUPS = (
    "independent",
    "normal",
    "shear",
    "coupling",
    "offdiagonal",
    "all",
)


def describe_thermoelastic_plots(
    result: ResultData | ThermoelasticResult,
) -> PlotInventory:
    """Describe plot families buildable from one thermoelastic result.

    Parameters
    ----------
    result : ResultData or ThermoelasticResult
        Complete public result envelope or direct thermoelastic payload. The
        isothermal--adiabatic comparison family requires a complete envelope
        because its builder re-evaluates the calibrated model through the
        archived public analysis context.

    Returns
    -------
    PlotInventory
        Result-aware stiffness properties, workflow stages, available plot
        families, exact stored grids, components, profiles, and scientific
        selection contexts.
    """
    payload = extract_thermoelastic_result(result)
    has_envelope = isinstance(result, ResultData)

    successful_fit_components = _successful_fit_components(payload)
    tensor_components = _tensor_components(payload)
    profile_components = _profile_components(payload)
    all_components = _ordered_union(
        successful_fit_components,
        tensor_components,
        profile_components,
    )

    has_temperature_curve = payload.temperature.size >= 2
    has_pressure_curve = payload.pressure.size >= 2
    has_two_dimensional_grid = has_temperature_curve and has_pressure_curve
    has_isothermal_grid = (
        has_two_dimensional_grid and payload.stiffness_isothermal is not None
    )
    has_domain = has_two_dimensional_grid
    has_profiles = bool(payload.profiles)

    grid_conditions = _grid_tensor_conditions(payload)
    profile_conditions = _profile_tensor_conditions(payload)
    compare_pressures, compare_temperatures = _compare_coordinates(payload)
    compare_axes: list[str] = []
    if has_envelope and compare_pressures and has_temperature_curve:
        compare_axes.append("temperature")
    if has_envelope and compare_temperatures and has_pressure_curve:
        compare_axes.append("pressure")
    has_compare = bool(compare_axes)

    representations: list[PlotRepresentationDescriptor] = []
    if successful_fit_components:
        representations.append(
            PlotRepresentationDescriptor(
                key="fit",
                name="Elastic-volume calibration fits",
                plot_kind="panel",
                description=(
                    "Observed C_IJ(V) values, cold finite-strain fit, propagated "
                    "confidence band, and optional residual panel."
                ),
                property_keys=("elastic_stiffness",),
                supported_contexts=("fit_component", "component_group"),
                constraints=(
                    "Only independent components with a successful archived fit are available.",
                    "Disabling residuals may return line specifications instead of panels.",
                ),
            )
        )
    if has_isothermal_grid:
        pt_property_keys = ["elastic_stiffness"]
        if _grid_uncertainty_available(payload):
            pt_property_keys.extend(
                ("stiffness_uncertainty", "relative_stiffness_uncertainty")
            )
        representations.append(
            PlotRepresentationDescriptor(
                key="pt",
                name="Pressure-temperature stiffness maps",
                plot_kind="contour",
                description=(
                    "Contour maps of selected tensor components on the stored "
                    "pressure-temperature grid."
                ),
                property_keys=tuple(pt_property_keys),
                supported_contexts=(
                    "stiffness_component",
                    "component_group",
                    "pt_tensor_condition",
                    "pt_quantity",
                    "temperature_grid",
                    "pressure_grid",
                ),
                constraints=(
                    "At least two temperatures and two pressures are required.",
                    "Faceted layouts may return a panel specification containing contour maps.",
                    "Uncertainty quantities require the matching tensor uncertainty field.",
                ),
            )
        )
    if has_profiles:
        profile_property_keys = ["elastic_stiffness", "relative_stiffness_change"]
        if _profile_uncertainty_available(payload):
            profile_property_keys.append("stiffness_uncertainty")
        representations.append(
            PlotRepresentationDescriptor(
                key="profile",
                name="Depth-dependent thermoelastic profiles",
                plot_kind="line",
                description=(
                    "Absolute stiffness or relative stiffness change along an "
                    "archived geothermobarometric depth path."
                ),
                property_keys=tuple(profile_property_keys),
                supported_contexts=(
                    "stiffness_component",
                    "component_group",
                    "profile_name",
                    "profile_tensor_condition",
                    "profile_mode",
                ),
                constraints=(
                    "A profile name is required when several profiles are archived.",
                    "Relative mode uses the first depth point unless an explicit reference depth is supplied.",
                    "Faceted layouts may return a panel specification containing line plots.",
                ),
            )
        )
    if has_compare:
        representations.append(
            PlotRepresentationDescriptor(
                key="compare",
                name="Isothermal-adiabatic stiffness comparison",
                plot_kind="line",
                description=(
                    "Direct comparison of C^T_IJ and C^S_IJ while pressure or "
                    "temperature is held fixed."
                ),
                property_keys=("elastic_stiffness",),
                supported_contexts=(
                    "stiffness_component",
                    "component_group",
                    "compare_axis",
                    "compare_fixed_pressure",
                    "compare_fixed_temperature",
                    "temperature_grid",
                    "pressure_grid",
                ),
                constraints=(
                    "Both isothermal and adiabatic stiffness fields are required.",
                    "Exactly one pressure or temperature coordinate is fixed.",
                    "The calibrated model is re-evaluated at the requested coordinate; nearest-grid snapping is not used.",
                    "Faceted layouts may return a panel specification containing line plots.",
                ),
            )
        )
    if has_domain:
        representations.append(
            PlotRepresentationDescriptor(
                key="domain",
                name="Calibration and analysis domain",
                plot_kind="contour",
                description=(
                    "Equilibrium-volume field on the stored P-T grid with "
                    "extrapolation masks and optional archived profile paths."
                ),
                property_keys=("equilibrium_volume",),
                supported_contexts=(
                    "temperature_grid",
                    "pressure_grid",
                    "profile_name",
                ) if has_profiles else ("temperature_grid", "pressure_grid"),
                constraints=(
                    "At least two temperatures and two pressures are required.",
                    "Profile overlays are available only for archived profile analyses.",
                ),
            )
        )

    representation_keys = {item.key for item in representations}
    properties: list[PlotPropertyDescriptor] = []
    stiffness_representations = tuple(
        key for key in ("fit", "pt", "profile", "compare") if key in representation_keys
    )
    if stiffness_representations:
        properties.append(
            PlotPropertyDescriptor(
                key="elastic_stiffness",
                name="Elastic stiffness component",
                symbol_math=r"C_{IJ}",
                symbol_plain="C_IJ",
                unit="GPa",
                description=(
                    "Symmetric Voigt stiffness component under the selected "
                    "thermodynamic condition."
                ),
                category="elastic_tensor",
                components=tuple(label.lower() for label in all_components),
                representations=stiffness_representations,
            )
        )
    uncertainty_representations = tuple(
        key
        for key in ("pt", "profile")
        if key in representation_keys
        and (
            (key == "pt" and _grid_uncertainty_available(payload))
            or (key == "profile" and _profile_uncertainty_available(payload))
        )
    )
    if uncertainty_representations:
        properties.append(
            PlotPropertyDescriptor(
                key="stiffness_uncertainty",
                name="One-sigma elastic stiffness uncertainty",
                symbol_math=r"\sigma(C_{IJ})",
                symbol_plain="σ(C_IJ)",
                unit="GPa",
                description=(
                    "One-standard-deviation uncertainty propagated for a "
                    "selected stiffness component."
                ),
                category="uncertainty",
                components=tuple(label.lower() for label in all_components),
                representations=uncertainty_representations,
            )
        )
    if "pt" in representation_keys and _grid_uncertainty_available(payload):
        properties.append(
            PlotPropertyDescriptor(
                key="relative_stiffness_uncertainty",
                name="Relative elastic stiffness uncertainty",
                symbol_math=r"100\,\sigma(C_{IJ})/|C_{IJ}|",
                symbol_plain="100 σ(C_IJ)/|C_IJ|",
                unit="%",
                description=(
                    "One-sigma uncertainty relative to the absolute component "
                    "magnitude; undefined at near-zero values."
                ),
                category="uncertainty",
                components=tuple(label.lower() for label in tensor_components),
                representations=("pt",),
            )
        )
    if "profile" in representation_keys:
        properties.append(
            PlotPropertyDescriptor(
                key="relative_stiffness_change",
                name="Relative elastic stiffness change",
                symbol_math=r"\Delta C_{IJ}/C_{IJ,ref}",
                symbol_plain="ΔC_IJ/C_IJ,ref",
                unit="%",
                description=(
                    "Percentage stiffness change relative to a selected depth "
                    "state on an archived profile."
                ),
                category="relative_property",
                components=tuple(label.lower() for label in profile_components),
                representations=("profile",),
            )
        )
    if "domain" in representation_keys:
        properties.append(
            PlotPropertyDescriptor(
                key="equilibrium_volume",
                name="Equilibrium volume",
                symbol_math="V",
                symbol_plain="V",
                unit="angstrom^3",
                description=(
                    "QHA equilibrium volume used to evaluate the calibrated "
                    "elastic model over the P-T domain."
                ),
                category="thermodynamic_state",
                representations=("domain",),
            )
        )

    contexts: list[PlotContextDescriptor] = []
    stages = _workflow_stages(payload)
    if stages:
        contexts.append(
            PlotContextDescriptor(
                key="workflow_stage",
                name="Available workflow stage",
                description=(
                    "Cumulative scientific stages represented by the archived "
                    "thermoelastic result."
                ),
                values=stages,
                selectable=False,
            )
        )
    if successful_fit_components:
        contexts.append(
            PlotContextDescriptor(
                key="fit_component",
                name="Fitted elastic stiffness component",
                description=(
                    "Independent Voigt component with a successful archived "
                    "elastic-volume fit."
                ),
                values=successful_fit_components,
                default=successful_fit_components[0],
            )
        )
    if all_components:
        default_component = _default_component(payload, all_components)
        contexts.append(
            PlotContextDescriptor(
                key="stiffness_component",
                name="Elastic stiffness component",
                description="Canonical symmetric Voigt component selected for plotting.",
                values=all_components,
                default=default_component,
            )
        )
        groups = _available_component_groups(payload)
        contexts.append(
            PlotContextDescriptor(
                key="component_group",
                name="Elastic component group",
                description=(
                    "Semantic component selection used when explicit components "
                    "are not supplied."
                ),
                values=groups,
                default="independent" if "independent" in groups else groups[0],
            )
        )
    if payload.temperature.size:
        contexts.append(
            PlotContextDescriptor(
                key="temperature_grid",
                name="Stored temperature grid",
                description="Exact temperature coordinates archived in the result.",
                values=tuple(float(value) for value in payload.temperature),
                unit="K",
                selectable=False,
            )
        )
    if payload.pressure.size:
        contexts.append(
            PlotContextDescriptor(
                key="pressure_grid",
                name="Stored pressure grid",
                description="Exact pressure coordinates archived in the result.",
                values=tuple(float(value) for value in payload.pressure),
                unit="GPa",
                selectable=False,
            )
        )
    if "pt" in representation_keys:
        contexts.extend(
            (
                PlotContextDescriptor(
                    key="pt_tensor_condition",
                    name="P-T tensor condition",
                    description="Thermodynamic stiffness field mapped on the P-T grid.",
                    values=grid_conditions,
                    default="isothermal",
                    required=True,
                ),
                PlotContextDescriptor(
                    key="pt_quantity",
                    name="P-T mapped quantity",
                    description="Value or propagated uncertainty represented by color.",
                    values=_pt_quantities(payload),
                    default="value",
                    required=True,
                ),
            )
        )
    if has_profiles:
        contexts.extend(
            (
                PlotContextDescriptor(
                    key="profile_name",
                    name="Archived profile",
                    description="Named geothermobarometric path stored in the result.",
                    values=tuple(sorted(payload.profiles)),
                    default=(
                        next(iter(sorted(payload.profiles)))
                        if len(payload.profiles) == 1
                        else None
                    ),
                    required=len(payload.profiles) > 1,
                ),
                PlotContextDescriptor(
                    key="profile_tensor_condition",
                    name="Profile tensor condition",
                    description="Thermodynamic stiffness field followed along depth.",
                    values=profile_conditions,
                    default="isothermal",
                    required=True,
                ),
                PlotContextDescriptor(
                    key="profile_mode",
                    name="Profile quantity mode",
                    description=(
                        "Absolute stiffness or percentage change relative to a "
                        "reference depth."
                    ),
                    values=("absolute", "relative"),
                    default="absolute",
                    required=True,
                ),
            )
        )
    if has_compare:
        contexts.append(
            PlotContextDescriptor(
                key="compare_axis",
                name="Comparison varying coordinate",
                description=(
                    "Coordinate varied while the complementary thermodynamic "
                    "coordinate is held fixed."
                ),
                values=tuple(compare_axes),
                default=compare_axes[0],
                required=True,
            )
        )
        if compare_pressures:
            contexts.append(
                PlotContextDescriptor(
                    key="compare_fixed_pressure",
                    name="Validated fixed pressures",
                    description=(
                        "Stored pressure coordinates whose complete temperature "
                        "sections are valid for isothermal-adiabatic comparison. "
                        "The public builder also accepts other exact coordinates "
                        "and re-evaluates the calibrated model."
                    ),
                    values=compare_pressures,
                    unit="GPa",
                    default=compare_pressures[0],
                )
            )
        if compare_temperatures:
            contexts.append(
                PlotContextDescriptor(
                    key="compare_fixed_temperature",
                    name="Validated fixed temperatures",
                    description=(
                        "Stored temperature coordinates whose complete pressure "
                        "sections are valid for isothermal-adiabatic comparison. "
                        "The public builder also accepts other exact coordinates "
                        "and re-evaluates the calibrated model."
                    ),
                    values=compare_temperatures,
                    unit="K",
                    default=compare_temperatures[0],
                )
            )

    warnings = _inventory_warnings(
        payload,
        successful_fit_components=successful_fit_components,
        has_two_dimensional_grid=has_two_dimensional_grid,
        has_isothermal_grid=has_isothermal_grid,
        has_compare=has_compare,
        has_envelope=has_envelope,
    )
    return PlotInventory(
        module="thermoelasticity",
        properties=tuple(properties),
        representations=tuple(representations),
        contexts=tuple(contexts),
        warnings=warnings,
    )


def _successful_fit_components(result: ThermoelasticResult) -> tuple[str, ...]:
    return tuple(
        label
        for label in result.independent_labels
        if label in result.component_fits
        and result.component_fits[label].fit is not None
        and result.component_fits[label].parameters is not None
    )


def _tensor_components(result: ThermoelasticResult) -> tuple[str, ...]:
    if result.stiffness_isothermal is None:
        return ()
    try:
        return resolve_components(result, group="all")
    except ValueError:
        return ()


def _profile_components(result: ThermoelasticResult) -> tuple[str, ...]:
    if not result.profiles:
        return ()
    # Profile tensors use the same symmetry reconstruction and component labels
    # as the parent result. The parent tensor inventory remains the authoritative
    # component catalogue when available.
    tensor = _tensor_components(result)
    return tensor if tensor else tuple(result.independent_labels)


def _ordered_union(*groups: Iterable[str]) -> tuple[str, ...]:
    ordered: list[str] = []
    seen: set[str] = set()
    for group in groups:
        for value in group:
            if value not in seen:
                ordered.append(value)
                seen.add(value)
    return tuple(ordered)


def _available_component_groups(result: ThermoelasticResult) -> tuple[str, ...]:
    groups: list[str] = []
    for group in _COMPONENT_GROUPS:
        try:
            values = resolve_components(result, group=group)  # type: ignore[arg-type]
        except ValueError:
            continue
        if values:
            groups.append(group)
    return tuple(groups)


def _default_component(
    result: ThermoelasticResult,
    available: tuple[str, ...],
) -> str:
    for label in result.independent_labels:
        if label in available:
            return label
    return available[0]


def _grid_tensor_conditions(result: ThermoelasticResult) -> tuple[str, ...]:
    conditions: list[str] = []
    if result.stiffness_isothermal is not None:
        conditions.append("isothermal")
    if result.stiffness_adiabatic is not None:
        conditions.append("adiabatic")
    return tuple(conditions)


def _profile_tensor_conditions(result: ThermoelasticResult) -> tuple[str, ...]:
    conditions = ["isothermal"]
    if result.profiles and all(
        profile.stiffness_adiabatic is not None for profile in result.profiles.values()
    ):
        conditions.append("adiabatic")
    return tuple(conditions)


def _grid_uncertainty_available(result: ThermoelasticResult) -> bool:
    arrays = (
        result.sigma_stiffness_isothermal,
        result.sigma_stiffness_adiabatic,
    )
    return any(
        value is not None and np.any(np.isfinite(np.asarray(value, dtype=np.float64)))
        for value in arrays
    )


def _profile_uncertainty_available(result: ThermoelasticResult) -> bool:
    return any(
        np.any(np.isfinite(profile.sigma_stiffness_isothermal))
        or (
            profile.sigma_stiffness_adiabatic is not None
            and np.any(np.isfinite(profile.sigma_stiffness_adiabatic))
        )
        for profile in result.profiles.values()
    )


def _pt_quantities(result: ThermoelasticResult) -> tuple[str, ...]:
    quantities = ["value"]
    if _grid_uncertainty_available(result):
        quantities.extend(("uncertainty", "relative-uncertainty"))
    return tuple(quantities)


def _workflow_stages(result: ThermoelasticResult) -> tuple[str, ...]:
    stages: list[str] = []
    if result.component_fits:
        stages.append("calibration")
    if result.stiffness_isothermal is not None:
        if result.temperature.size == 1 and result.pressure.size == 1:
            stages.append("point")
        elif result.temperature.size >= 1 and result.pressure.size >= 1:
            stages.append("grid")
    if result.profiles:
        stages.append("profile")
    return tuple(stages)


def _compare_coordinates(
    result: ThermoelasticResult,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if result.stiffness_isothermal is None or result.stiffness_adiabatic is None:
        return (), ()
    valid = result.adiabatic_valid_mask
    mask: NDArray[np.bool_]
    if valid is None:
        mask = np.ones(
            (result.temperature.size, result.pressure.size),
            dtype=np.bool_,
        )
    else:
        mask = np.asarray(valid, dtype=np.bool_)
    pressures = tuple(
        float(result.pressure[index])
        for index in range(result.pressure.size)
        if bool(np.all(mask[:, index]))
    )
    temperatures = tuple(
        float(result.temperature[index])
        for index in range(result.temperature.size)
        if bool(np.all(mask[index, :]))
    )
    return pressures, temperatures


def _inventory_warnings(
    result: ThermoelasticResult,
    *,
    successful_fit_components: tuple[str, ...],
    has_two_dimensional_grid: bool,
    has_isothermal_grid: bool,
    has_compare: bool,
    has_envelope: bool,
) -> tuple[str, ...]:
    warnings: list[str] = []
    failed = tuple(
        label for label in result.independent_labels if label not in successful_fit_components
    )
    if result.component_fits and failed:
        warnings.append(
            "Independent fit plots exclude components without a successful "
            f"numerical fit: {', '.join(failed)}."
        )
    if result.stiffness_isothermal is not None and not has_two_dimensional_grid:
        warnings.append(
            "The reconstructed result is a point or one-dimensional section; "
            "P-T contour and domain plots require at least two temperatures and "
            "two pressures."
        )
    if has_two_dimensional_grid and not has_isothermal_grid:
        warnings.append(
            "The result contains a P-T volume domain but no reconstructed "
            "isothermal stiffness field; stiffness maps are unavailable."
        )
    if result.stiffness_isothermal is not None and result.stiffness_adiabatic is None:
        warnings.append(
            "Adiabatic stiffness is unavailable; isothermal-adiabatic comparison "
            "plots are not advertised."
        )
    elif result.stiffness_adiabatic is not None and not has_compare:
        warnings.append(
            "No complete valid stored pressure or temperature section is "
            "available for isothermal-adiabatic comparison."
        )
    if result.stiffness_adiabatic is not None and not has_envelope:
        warnings.append(
            "Comparison discovery requires the complete ResultData envelope "
            "because the public builder re-evaluates the calibrated model."
        )
    return tuple(warnings)


__all__ = ["describe_thermoelastic_plots"]
