# -*- coding: utf-8 -*-

"""Result-aware public plot inventory for seismic-wave workflows."""

from __future__ import annotations

from dataclasses import dataclass

from quantas.core.physics.seismic import SamplingLevel, WaveMode
from quantas.models import (
    PlotContextDescriptor,
    PlotInventory,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
)
from quantas.modules.seismic.models import SeismicResult


@dataclass(frozen=True, slots=True)
class _SeismicPropertyDefinition:
    """Static metadata and availability level for one seismic scalar field."""

    key: str
    name: str
    symbol_math: str
    symbol_plain: str
    unit: str | None
    description: str
    category: str
    components: tuple[str, ...]
    minimum_level: SamplingLevel
    summary: bool = False


_PHASE_PROPERTIES = (
    _SeismicPropertyDefinition(
        "phase_v_p",
        "P-wave phase velocity",
        "V_P",
        "Vₚ",
        "km s^-1",
        "Phase velocity of the quasi-longitudinal acoustic mode.",
        "phase_velocity",
        ("v_p",),
        SamplingLevel.PHASE,
        True,
    ),
    _SeismicPropertyDefinition(
        "phase_v_s1",
        "Fast shear-wave phase velocity",
        "V_{S1}",
        "Vₛ₁",
        "km s^-1",
        "Phase velocity of the faster quasi-shear acoustic mode.",
        "phase_velocity",
        ("v_s1",),
        SamplingLevel.PHASE,
        True,
    ),
    _SeismicPropertyDefinition(
        "phase_v_s2",
        "Slow shear-wave phase velocity",
        "V_{S2}",
        "Vₛ₂",
        "km s^-1",
        "Phase velocity of the slower quasi-shear acoustic mode.",
        "phase_velocity",
        ("v_s2",),
        SamplingLevel.PHASE,
        True,
    ),
    _SeismicPropertyDefinition(
        "shear_anisotropy",
        "Directional shear-wave anisotropy",
        "A_S",
        "Aₛ",
        "%",
        "Normalized directional separation of the fast and slow shear modes.",
        "shear_diagnostic",
        ("v_s1", "v_s2"),
        SamplingLevel.PHASE,
        True,
    ),
    _SeismicPropertyDefinition(
        "shear_splitting",
        "Shear-wave velocity splitting",
        r"\Delta V_S",
        "ΔVₛ",
        "km s^-1",
        "Absolute directional phase-velocity difference V_S1 minus V_S2.",
        "shear_diagnostic",
        ("v_s1", "v_s2"),
        SamplingLevel.PHASE,
    ),
    _SeismicPropertyDefinition(
        "phase_v_p_over_v_s1",
        "P-to-fast-shear phase-velocity ratio",
        "V_P/V_{S1}",
        "Vₚ/Vₛ₁",
        None,
        "Directional ratio of P-wave and fast shear-wave phase velocities.",
        "velocity_ratio",
        ("v_p", "v_s1"),
        SamplingLevel.PHASE,
        True,
    ),
    _SeismicPropertyDefinition(
        "phase_v_p_over_v_s2",
        "P-to-slow-shear phase-velocity ratio",
        "V_P/V_{S2}",
        "Vₚ/Vₛ₂",
        None,
        "Directional ratio of P-wave and slow shear-wave phase velocities.",
        "velocity_ratio",
        ("v_p", "v_s2"),
        SamplingLevel.PHASE,
        True,
    ),
)


def _mode_definitions(
    *,
    prefix: str,
    name: str,
    symbol_template: str,
    plain_template: str,
    unit: str | None,
    description: str,
    category: str,
    minimum_level: SamplingLevel,
) -> tuple[_SeismicPropertyDefinition, ...]:
    """Create the three mode-resolved definitions for one field family."""
    definitions = []
    labels = {
        WaveMode.V_P: ("P", "P", "P", "ₚ", "v_p"),
        WaveMode.V_S1: ("S1", "{S1}", "S1", "ₛ₁", "v_s1"),
        WaveMode.V_S2: ("S2", "{S2}", "S2", "ₛ₂", "v_s2"),
    }
    for mode in (WaveMode.V_P, WaveMode.V_S1, WaveMode.V_S2):
        math_label, math_subscript, name_label, plain_suffix, component = labels[mode]
        definitions.append(
            _SeismicPropertyDefinition(
                key=f"{prefix}_{mode.value}",
                name=name.format(mode=name_label),
                symbol_math=symbol_template.format(
                    mode=math_label,
                    subscript=math_subscript,
                ),
                symbol_plain=plain_template.format(mode=plain_suffix),
                unit=unit,
                description=description.format(mode=name_label),
                category=category,
                components=(component,),
                minimum_level=minimum_level,
            )
        )
    return tuple(definitions)


_GROUP_PROPERTIES = _mode_definitions(
    prefix="group",
    name="{mode}-wave group velocity",
    symbol_template="V_{{g,{mode}}}",
    plain_template="Vg,{mode}",
    unit="km s^-1",
    description="Magnitude of the {mode}-wave group-velocity vector.",
    category="group_velocity",
    minimum_level=SamplingLevel.GROUP,
)
_POWER_FLOW_PROPERTIES = _mode_definitions(
    prefix="power_flow",
    name="{mode}-wave power-flow angle",
    symbol_template=r"\psi_{subscript}",
    plain_template="ψ{mode}",
    unit="degree",
    description="Angle between phase and group directions for the {mode} mode.",
    category="power_flow",
    minimum_level=SamplingLevel.GROUP,
)
_ENHANCEMENT_PROPERTIES = _mode_definitions(
    prefix="log10_enhancement",
    name="{mode}-wave logarithmic enhancement",
    symbol_template=r"\log_{{10}}(A_{subscript})",
    plain_template="log₁₀(A{mode})",
    unit=None,
    description="Base-10 logarithm of the {mode}-wave focusing enhancement.",
    category="enhancement",
    minimum_level=SamplingLevel.ENHANCEMENT,
)
_PROPERTY_DEFINITIONS = (
    *_PHASE_PROPERTIES,
    *_GROUP_PROPERTIES,
    *_POWER_FLOW_PROPERTIES,
    *_ENHANCEMENT_PROPERTIES,
)
_LEVEL_ORDER = {
    SamplingLevel.PHASE: 0,
    SamplingLevel.GROUP: 1,
    SamplingLevel.ENHANCEMENT: 2,
}


def available_seismic_plot_properties(
    result: SeismicResult,
) -> tuple[PlotPropertyDescriptor, ...]:
    """Return scalar properties actually available in a seismic result."""
    definitions = tuple(
        item
        for item in _PROPERTY_DEFINITIONS
        if _definition_is_available(result, item)
    )
    return tuple(
        PlotPropertyDescriptor(
            key=item.key,
            name=item.name,
            symbol_math=item.symbol_math,
            symbol_plain=item.symbol_plain,
            unit=item.unit,
            description=item.description,
            category=item.category,
            components=item.components,
            representations=tuple(
                representation
                for representation in (
                    "spherical_map",
                    "spherical_summary" if item.summary else None,
                    "property_surface_3d",
                )
                if representation is not None
            ),
        )
        for item in definitions
    )


def describe_seismic_plots(result: SeismicResult) -> PlotInventory:
    """Describe seismic maps, summaries, and surfaces for one result."""
    properties = available_seismic_plot_properties(result)
    property_keys = tuple(item.key for item in properties)
    summary_keys = tuple(
        item.key
        for item in properties
        if "spherical_summary" in item.representations
    )
    tracking_values = (False, True) if result.field.tracking is not None else (False,)
    surface_types = ["phase", "slowness"]
    if result.field.group is not None:
        surface_types.append("group")

    contexts = (
        PlotContextDescriptor(
            key="projection",
            name="Spherical projection",
            description="Projection used for two-dimensional directional maps.",
            values=("equal_area", "stereographic"),
            default="equal_area",
            required=True,
        ),
        PlotContextDescriptor(
            key="sampled_hemisphere",
            name="Sampled hemisphere",
            description="Angular domain stored in the seismic result.",
            values=(result.grid.hemisphere.value,),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="sampling_level",
            name="Sampling level",
            description="Highest acoustic field available in the result.",
            values=(result.field.level.value,),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="extrema_markers",
            name="Sampled extrema markers",
            description="Attach sampled minimum and maximum directions.",
            values=(False, True),
            default=True,
        ),
        PlotContextDescriptor(
            key="polarization_overlay",
            name="Polarization-axis overlay",
            description="Attach decimated tracked polarization axes when available.",
            values=tracking_values,
            default=False,
        ),
        PlotContextDescriptor(
            key="surface_geometry",
            name="Surface geometry",
            description=(
                "Physical geometry uses the natural wavefront radius or vector; "
                "unit-sphere geometry carries only the scalar field."
            ),
            values=("physical", "unit_sphere"),
            default="unit_sphere",
            required=True,
        ),
        PlotContextDescriptor(
            key="antipodal_completion",
            name="Antipodal completion",
            description="Mirror a hemispherical result to show the complete surface.",
            values=(False, True),
            default=True,
        ),
        PlotContextDescriptor(
            key="surface_type",
            name="Acoustic surface type",
            description="Physical acoustic quantity represented by an explicit surface.",
            values=tuple(surface_types),
            required=True,
        ),
        PlotContextDescriptor(
            key="wave_mode",
            name="Acoustic mode",
            description="Quasi-longitudinal or split quasi-shear acoustic branch.",
            values=("v_p", "v_s1", "v_s2"),
            required=True,
        ),
    )
    representations = (
        PlotRepresentationDescriptor(
            key="spherical_map",
            name="Directional spherical map",
            plot_kind="spherical_map",
            description="Two-dimensional map of one sampled directional scalar field.",
            property_keys=property_keys,
            supported_contexts=(
                "projection",
                "sampled_hemisphere",
                "sampling_level",
                "extrema_markers",
                "polarization_overlay",
            ),
            constraints=(
                "Polarization overlays require tracked axes and a mode-resolved field.",
            ),
        ),
        PlotRepresentationDescriptor(
            key="spherical_summary",
            name="Seismic six-panel summary",
            plot_kind="spherical_summary",
            description=(
                "Fixed summary of phase velocities, directional shear anisotropy, "
                "and P-to-S velocity ratios."
            ),
            property_keys=summary_keys,
            supported_contexts=(
                "projection",
                "sampled_hemisphere",
                "extrema_markers",
                "polarization_overlay",
            ),
        ),
        PlotRepresentationDescriptor(
            key="property_surface_3d",
            name="Scalar-property three-dimensional surface",
            plot_kind="surface",
            description=(
                "Any available sampled scalar property on a unit sphere or its "
                "natural physical carrier when one exists."
            ),
            property_keys=property_keys,
            supported_contexts=(
                "surface_geometry",
                "sampled_hemisphere",
                "antipodal_completion",
                "polarization_overlay",
            ),
            constraints=(
                "Polarization overlays are supported only on unit-sphere geometry.",
            ),
        ),
        PlotRepresentationDescriptor(
            key="acoustic_surface_3d",
            name="Mode-resolved acoustic surface",
            plot_kind="surface",
            description=(
                "Phase-velocity, slowness, or group-wavefront surface selected "
                "by acoustic mode."
            ),
            supported_contexts=(
                "surface_type",
                "wave_mode",
                "surface_geometry",
                "sampled_hemisphere",
                "antipodal_completion",
                "polarization_overlay",
            ),
            constraints=(
                "Group surfaces are advertised only when group fields are stored.",
                "Polarization overlays are supported only on unit-sphere geometry.",
            ),
        ),
    )
    warnings = ()
    if result.field.tracking is None:
        warnings = (
            "Tracked polarization axes are unavailable; polarization-overlay "
            "contexts are restricted to False.",
        )
    return PlotInventory(
        module="seismic",
        properties=properties,
        representations=representations,
        contexts=contexts,
        warnings=warnings,
    )


def _definition_is_available(
    result: SeismicResult,
    definition: _SeismicPropertyDefinition,
) -> bool:
    """Return whether the result contains the required acoustic field level."""
    if _LEVEL_ORDER[result.field.level] < _LEVEL_ORDER[definition.minimum_level]:
        return False
    if definition.minimum_level is SamplingLevel.GROUP:
        return result.field.group is not None
    if definition.minimum_level is SamplingLevel.ENHANCEMENT:
        return result.field.enhancement is not None
    return True


__all__ = [
    "available_seismic_plot_properties",
    "describe_seismic_plots",
]
