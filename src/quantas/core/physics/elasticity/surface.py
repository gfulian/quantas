# -*- coding: utf-8 -*-

"""Three-dimensional sampling of directional elastic properties.

The functions in this module convert a neutral batched directional field into
radial surface data.  Numerical sampling is delegated to
:mod:`quantas.core.physics.elasticity.sampling`; this module only separates
signed branches and constructs Cartesian surface coordinates.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass, field
from typing import Literal, TypeAlias

import numpy as np

from quantas.core.geometry import SphericalGrid, create_spherical_grid

from .sampling import (
    SUPPORTED_DIRECTIONAL_PROPERTIES,
    ElasticDirectionalField,
    ElasticDirectionalProperty,
    exact_transverse_extrema,
    sample_elastic_directional_field,
)
from .tensor import ElasticTensor


ElasticSurfaceProperty: TypeAlias = ElasticDirectionalProperty
ProgressCallback = Callable[[int, int], None]

SUPPORTED_SURFACE_PROPERTIES: tuple[ElasticSurfaceProperty, ...] = (
    SUPPORTED_DIRECTIONAL_PROPERTIES
)


@dataclass(frozen=True, slots=True)
class TransverseExtrema:
    """Minimum and maximum over the transverse measurement angle.

    Parameters
    ----------
    minimum, maximum : float
        Extremal property values.
    minimum_angle, maximum_angle : float
        Transverse angles in radians at which the extrema occur.
    """

    minimum: float
    maximum: float
    minimum_angle: float
    maximum_angle: float


@dataclass(frozen=True, slots=True)
class DirectionalSurface:
    """One three-dimensional elastic-property surface.

    Parameters
    ----------
    key : str
        Stable branch identifier, for example ``"shear_minimum"``.
    property_name : str
        Parent property name.
    branch : str
        Human-readable branch name such as ``"minimum"`` or ``"negative"``.
    unit : str
        Physical unit of ``values`` and radial coordinates.
    theta, phi : ndarray
        Two-dimensional angular grids in radians.
    radius : ndarray
        Non-negative radial distance used to construct the surface geometry.
    values : ndarray
        Original signed physical values used for color mapping and tooltips.
    x, y, z : ndarray
        Cartesian surface coordinates.
    metadata : dict, optional
        Additional frontend-neutral information.
    """

    key: str
    property_name: str
    branch: str
    unit: str
    theta: np.ndarray
    phi: np.ndarray
    radius: np.ndarray
    values: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(slots=True)
class ElasticitySurfaceCollection:
    """Collection of elasticity surfaces and non-fatal warnings.

    Parameters
    ----------
    surfaces : dict, optional
        Surfaces indexed by stable branch key.
    warnings : list of str, optional
        Messages for branches that could not or did not need to be generated.
    metadata : dict, optional
        Numerical sampling metadata and diagnostics.
    """

    surfaces: dict[str, DirectionalSurface] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, object] = field(default_factory=dict)


def sample_elasticity_surfaces(
    tensor: ElasticTensor,
    *,
    ntheta: int = 61,
    nphi: int = 121,
    properties: Iterable[ElasticSurfaceProperty] = SUPPORTED_SURFACE_PROPERTIES,
    batch_size: int = 65536,
    progress_callback: ProgressCallback | None = None,
) -> ElasticitySurfaceCollection:
    """Sample selected three-dimensional elastic-property surfaces.

    Parameters
    ----------
    tensor : ElasticTensor
        Stable three-dimensional elastic tensor.
    ntheta, nphi : int, optional
        Polar and azimuthal grid sizes.
    properties : iterable of str, optional
        Any subset of ``young``, ``compressibility``, ``shear`` and ``poisson``.
    batch_size : int, optional
        Maximum number of directions evaluated in one NumPy batch.
    progress_callback : callable or None, optional
        Callback receiving ``(current, total)`` after batched numerical work.

    Returns
    -------
    ElasticitySurfaceCollection
        Directional surfaces. Empty sign branches are omitted and
        reported through ``warnings``.

    Raises
    ------
    ValueError
        If a property or numerical sampling option is invalid.
    """
    grid = create_spherical_grid(ntheta, nphi)
    selected = _normalize_properties(properties)
    field = sample_elastic_directional_field(
        tensor,
        grid.theta_grid,
        grid.phi_grid,
        properties=selected,
        batch_size=batch_size,
        progress_callback=progress_callback,
    )
    collection = ElasticitySurfaceCollection(
        metadata={
            "sampling": "exact_batched",
            "ntheta": grid.shape[0],
            "nphi": grid.shape[1],
            "batch_size": field.batch_size,
            "diagnostics": asdict(field.diagnostics),
        }
    )
    if "young" in selected:
        _append_young_surface(collection, grid, field)
    if "compressibility" in selected:
        _append_compressibility_surfaces(collection, grid, field)
    if "shear" in selected:
        _append_shear_surfaces(collection, grid, field)
    if "poisson" in selected:
        _append_poisson_surfaces(collection, grid, field)
    return collection


def _append_young_surface(
    collection: ElasticitySurfaceCollection,
    grid: SphericalGrid,
    field: ElasticDirectionalField,
) -> None:
    """Append the Young-modulus radial surface."""
    assert field.young_modulus is not None
    collection.surfaces["young"] = _make_surface(
        grid,
        key="young",
        property_name="young",
        branch="value",
        unit="GPa",
        radius=field.young_modulus,
        values=field.young_modulus,
    )


def _append_compressibility_surfaces(
    collection: ElasticitySurfaceCollection,
    grid: SphericalGrid,
    field: ElasticDirectionalField,
) -> None:
    """Append signed linear-compressibility surfaces."""
    assert field.linear_compressibility is not None
    _append_signed_surfaces(
        collection,
        grid,
        property_name="compressibility",
        unit="TPa^-1",
        values=field.linear_compressibility,
    )


def _append_shear_surfaces(
    collection: ElasticitySurfaceCollection,
    grid: SphericalGrid,
    field: ElasticDirectionalField,
) -> None:
    """Append minimum and maximum shear-modulus surfaces."""
    assert field.shear_minimum is not None
    assert field.shear_maximum is not None
    for key, branch, values in (
        ("shear_minimum", "minimum", field.shear_minimum),
        ("shear_maximum", "maximum", field.shear_maximum),
    ):
        collection.surfaces[key] = _make_surface(
            grid,
            key=key,
            property_name="shear",
            branch=branch,
            unit="GPa",
            radius=values,
            values=values,
        )


def _append_poisson_surfaces(
    collection: ElasticitySurfaceCollection,
    grid: SphericalGrid,
    field: ElasticDirectionalField,
) -> None:
    """Append signed minimum and maximum Poisson-ratio surfaces."""
    assert field.poisson_minimum is not None
    assert field.poisson_maximum is not None
    negative = np.minimum(field.poisson_minimum, 0.0)
    positive_minimum = np.maximum(field.poisson_minimum, 0.0)
    branches = (
        (
            "poisson_negative",
            "negative",
            np.abs(negative),
            np.where(negative < 0.0, negative, np.nan),
            "No negative Poisson-ratio values were found.",
        ),
        (
            "poisson_minimum",
            "minimum",
            positive_minimum,
            np.where(positive_minimum > 0.0, positive_minimum, np.nan),
            "No positive minimum Poisson-ratio branch was found.",
        ),
        (
            "poisson_maximum",
            "maximum",
            np.abs(field.poisson_maximum),
            field.poisson_maximum,
            "No non-zero maximum Poisson-ratio branch was found.",
        ),
    )
    for key, branch, radius, values, warning in branches:
        _append_branch(
            collection,
            grid,
            key=key,
            property_name="poisson",
            branch=branch,
            unit="",
            radius=radius,
            values=values,
            empty_warning=warning,
        )


def transverse_extrema(
    tensor: ElasticTensor,
    theta: float,
    phi: float,
    *,
    kind: Literal["shear", "poisson"],
    coarse_points: int = 37,
) -> TransverseExtrema:
    """Determine exact transverse extrema for one longitudinal direction.

    The retained ``coarse_points`` argument is accepted for source
    compatibility but is no longer used.  Transverse extrema are solved
    algebraically from a projected symmetric ``2 x 2`` quadratic form.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor.
    theta, phi : float
        Longitudinal direction in radians.
    kind : {"shear", "poisson"}
        Property optimized over the transverse angle.
    coarse_points : int, optional
        Deprecated compatibility argument. Values below eight remain invalid
        to preserve the historical validation contract.

    Returns
    -------
    TransverseExtrema
        Exact extremal values and corresponding transverse angles.

    Raises
    ------
    ValueError
        If ``kind`` is unsupported or ``coarse_points`` is below eight.
    """
    if coarse_points < 8:
        raise ValueError("coarse_points must be at least 8.")
    exact = exact_transverse_extrema(
        tensor,
        theta,
        phi,
        kind=kind,
    )
    return TransverseExtrema(
        minimum=exact.minimum,
        maximum=exact.maximum,
        minimum_angle=exact.minimum_angle,
        maximum_angle=exact.maximum_angle,
    )


def _normalize_properties(
    properties: Iterable[ElasticSurfaceProperty],
) -> tuple[ElasticSurfaceProperty, ...]:
    """Return unique supported properties while preserving order."""
    selected: list[ElasticSurfaceProperty] = []
    for property_name in properties:
        if property_name not in SUPPORTED_SURFACE_PROPERTIES:
            raise ValueError(f"unsupported 3D elasticity property: {property_name}")
        if property_name not in selected:
            selected.append(property_name)
    return tuple(selected)


def _append_signed_surfaces(
    collection: ElasticitySurfaceCollection,
    grid: SphericalGrid,
    *,
    property_name: str,
    unit: str,
    values: np.ndarray,
) -> None:
    """Append separate positive and negative branches."""
    positive = np.maximum(values, 0.0)
    negative = np.minimum(values, 0.0)
    _append_branch(
        collection,
        grid,
        key=f"{property_name}_positive",
        property_name=property_name,
        branch="positive",
        unit=unit,
        radius=positive,
        values=np.where(positive > 0.0, positive, np.nan),
        empty_warning=f"No positive {property_name} values were found.",
    )
    _append_branch(
        collection,
        grid,
        key=f"{property_name}_negative",
        property_name=property_name,
        branch="negative",
        unit=unit,
        radius=np.abs(negative),
        values=np.where(negative < 0.0, negative, np.nan),
        empty_warning=f"No negative {property_name} values were found.",
    )


def _append_branch(
    collection: ElasticitySurfaceCollection,
    grid: SphericalGrid,
    *,
    key: str,
    property_name: str,
    branch: str,
    unit: str,
    radius: np.ndarray,
    values: np.ndarray,
    empty_warning: str,
) -> None:
    """Append a non-empty branch or record a warning."""
    finite_radius = np.asarray(radius, dtype=float)
    if not np.any(np.isfinite(finite_radius) & (finite_radius > 1.0e-14)):
        collection.warnings.append(empty_warning)
        return
    collection.surfaces[key] = _make_surface(
        grid,
        key=key,
        property_name=property_name,
        branch=branch,
        unit=unit,
        radius=finite_radius,
        values=values,
    )


def _make_surface(
    grid: SphericalGrid,
    *,
    key: str,
    property_name: str,
    branch: str,
    unit: str,
    radius: np.ndarray,
    values: np.ndarray,
) -> DirectionalSurface:
    """Construct Cartesian coordinates for one radial surface."""
    radius_array = np.asarray(radius, dtype=float)
    value_array = np.asarray(values, dtype=float)
    if radius_array.shape != grid.shape or value_array.shape != grid.shape:
        raise ValueError("surface arrays must match the spherical grid shape.")
    x = radius_array * grid.directions[..., 0]
    y = radius_array * grid.directions[..., 1]
    z = radius_array * grid.directions[..., 2]
    return DirectionalSurface(
        key=key,
        property_name=property_name,
        branch=branch,
        unit=unit,
        theta=grid.theta_grid.copy(),
        phi=grid.phi_grid.copy(),
        radius=radius_array.copy(),
        values=value_array.copy(),
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        z=np.asarray(z, dtype=float),
        metadata={
            "ntheta": grid.shape[0],
            "nphi": grid.shape[1],
            "sampling": "exact_batched",
        },
    )
