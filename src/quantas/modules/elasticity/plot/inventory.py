# -*- coding: utf-8 -*-

"""Result-aware public plot inventory for elasticity workflows."""

from __future__ import annotations

from dataclasses import dataclass

from quantas.models import (
    PlotContextDescriptor,
    PlotInventory,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
)
from quantas.modules.elasticity.models import ElasticityResult


_PLANES = ("xy", "xz", "yz")


@dataclass(frozen=True, slots=True)
class _ElasticityPropertyDefinition:
    """Static scientific metadata for one directional property family."""

    key: str
    stored_name: str
    name: str
    symbol_math: str
    symbol_plain: str
    unit: str | None
    description: str
    components: tuple[str, ...]


_PROPERTY_DEFINITIONS = (
    _ElasticityPropertyDefinition(
        key="young",
        stored_name="young_modulus",
        name="Young's modulus",
        symbol_math="E",
        symbol_plain="E",
        unit="GPa",
        description="Directional tensile stiffness along the sampled direction.",
        components=("value",),
    ),
    _ElasticityPropertyDefinition(
        key="compressibility",
        stored_name="linear_compressibility",
        name="Linear compressibility",
        symbol_math=r"\beta",
        symbol_plain="β",
        unit="TPa^-1",
        description="Signed longitudinal strain response to hydrostatic pressure.",
        components=("positive", "negative"),
    ),
    _ElasticityPropertyDefinition(
        key="shear",
        stored_name="shear_modulus",
        name="Shear modulus",
        symbol_math="G",
        symbol_plain="G",
        unit="GPa",
        description=(
            "Minimum and maximum shear stiffness over transverse directions for "
            "each longitudinal direction."
        ),
        components=("minimum", "maximum"),
    ),
    _ElasticityPropertyDefinition(
        key="poisson",
        stored_name="poisson_ratio",
        name="Poisson's ratio",
        symbol_math=r"\nu",
        symbol_plain="ν",
        unit=None,
        description=(
            "Extremal transverse-to-longitudinal strain ratio for each "
            "longitudinal direction."
        ),
        components=("negative", "minimum_positive", "maximum"),
    ),
)


def describe_elasticity_plots(result: ElasticityResult) -> PlotInventory:
    """Describe plots that can be built from one elasticity result.

    Two-dimensional polar plots are advertised only when all three principal
    planes contain the required stored property. Three-dimensional surfaces are
    advertised when a stiffness matrix is available and the result is not
    explicitly marked mechanically unstable. Such surfaces may be reused from
    persistence or calculated transiently without mutating the result.
    """
    polar_keys = tuple(
        item.key for item in _PROPERTY_DEFINITIONS if _has_2d_property(result, item)
    )
    surface_available = result.stiffness is not None and (
        result.stability is None or result.stability.is_stable
    )
    surface_keys = (
        tuple(item.key for item in _PROPERTY_DEFINITIONS)
        if surface_available
        else ()
    )

    representations: list[PlotRepresentationDescriptor] = []
    contexts: list[PlotContextDescriptor] = []
    if polar_keys:
        contexts.append(
            PlotContextDescriptor(
                key="principal_plane",
                name="Principal plane",
                description=(
                    "Cartesian principal planes included as panels in each "
                    "two-dimensional polar specification."
                ),
                values=_PLANES,
                selectable=False,
            )
        )
        representations.append(
            PlotRepresentationDescriptor(
                key="polar_2d",
                name="Principal-plane polar sections",
                plot_kind="polar",
                description=(
                    "Directional property sections on the xy, xz, and yz "
                    "planes."
                ),
                property_keys=polar_keys,
                supported_contexts=("principal_plane",),
                constraints=(
                    "All three principal planes must be present for a property.",
                ),
            )
        )
    if surface_keys:
        contexts.append(
            PlotContextDescriptor(
                key="surface_geometry",
                name="Surface geometry",
                description=(
                    "Physical-radius surfaces encode the property in radius; "
                    "unit-sphere surfaces encode it only as a scalar field."
                ),
                values=("physical", "unit_sphere"),
                default="physical",
                required=True,
            )
        )
        representations.append(
            PlotRepresentationDescriptor(
                key="surface_3d",
                name="Directional three-dimensional surface",
                plot_kind="surface",
                description=(
                    "Three-dimensional directional surfaces constructed from "
                    "the stored stiffness tensor."
                ),
                property_keys=surface_keys,
                supported_contexts=("surface_geometry",),
                constraints=(
                    "A mechanically stable stiffness matrix is required.",
                    "Supplying sampling options may trigger a transient in-memory calculation.",
                ),
            )
        )

    representation_by_property = {
        item.key: tuple(
            representation.key
            for representation in representations
            if item.key in representation.property_keys
        )
        for item in _PROPERTY_DEFINITIONS
    }
    properties = tuple(
        PlotPropertyDescriptor(
            key=item.key,
            name=item.name,
            symbol_math=item.symbol_math,
            symbol_plain=item.symbol_plain,
            unit=item.unit,
            description=item.description,
            category="directional_elasticity",
            components=item.components,
            representations=representation_by_property[item.key],
        )
        for item in _PROPERTY_DEFINITIONS
        if representation_by_property[item.key]
    )

    warnings: list[str] = []
    if result.properties_2d and not polar_keys:
        warnings.append(
            "Stored two-dimensional elasticity data are incomplete and cannot "
            "produce a principal-plane polar specification."
        )
    if result.stiffness is None:
        warnings.append(
            "Three-dimensional surfaces are unavailable because the result "
            "does not contain a stiffness matrix."
        )
    elif result.stability is not None and not result.stability.is_stable:
        warnings.append(
            "Three-dimensional surfaces are unavailable for a mechanically "
            "unstable stiffness matrix."
        )

    return PlotInventory(
        module="elasticity",
        properties=properties,
        representations=tuple(representations),
        contexts=tuple(contexts),
        warnings=tuple(warnings),
    )


def _has_2d_property(
    result: ElasticityResult,
    definition: _ElasticityPropertyDefinition,
) -> bool:
    """Return whether one property is present on all principal planes."""
    return all(
        plane in result.properties_2d
        and definition.stored_name in result.properties_2d[plane]
        for plane in _PLANES
    )


__all__ = ["describe_elasticity_plots"]
