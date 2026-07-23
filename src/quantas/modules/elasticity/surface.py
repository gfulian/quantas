# -*- coding: utf-8 -*-

"""High-level construction and selection of elasticity surfaces."""

from __future__ import annotations

from collections.abc import Callable, Iterable

from quantas.core.physics.elasticity import (
    ElasticTensor,
    ElasticitySurfaceCollection,
    ElasticSurfaceProperty,
    sample_elasticity_surfaces,
    specialize_elastic_tensor,
)
from quantas.models import ResultData
from quantas.modules.elasticity.models import (
    ElasticityResult,
    ElasticitySurfaceOptions,
)


ProgressCallback = Callable[[int, int], None]


def calculate_elasticity_surfaces(
    source: ElasticTensor | ElasticityResult | ResultData,
    options: ElasticitySurfaceOptions | None = None,
    progress_callback: ProgressCallback | None = None,
) -> ElasticitySurfaceCollection:
    """Calculate three-dimensional elasticity surfaces in memory.

    Parameters
    ----------
    source : ElasticTensor, ElasticityResult or ResultData
        Elastic tensor or a completed elasticity result containing stiffness
        data.
    options : ElasticitySurfaceOptions or None, optional
        Spherical sampling options.
    progress_callback : callable or None, optional
        Callback receiving ``(current, total)``.

    Returns
    -------
    ElasticitySurfaceCollection
        Newly calculated numerical surfaces. The returned collection is not
        added to ``source`` automatically.

    Raises
    ------
    ValueError
        If the source does not contain a stable elasticity result.
    """
    opts = options or ElasticitySurfaceOptions()
    tensor = _resolve_tensor(source)
    return sample_elasticity_surfaces(
        tensor,
        ntheta=opts.ntheta,
        nphi=opts.nphi,
        properties=opts.properties,
        batch_size=opts.batch_size,
        progress_callback=progress_callback,
    )


def resolve_elasticity_surfaces(
    source: ElasticTensor | ElasticityResult | ResultData,
    options: ElasticitySurfaceOptions | None = None,
    *,
    properties: Iterable[ElasticSurfaceProperty] | None = None,
    progress_callback: ProgressCallback | None = None,
) -> ElasticitySurfaceCollection:
    """Return persisted surfaces when possible or calculate them transiently.

    Parameters
    ----------
    source : ElasticTensor, ElasticityResult or ResultData
        Elastic tensor or completed elasticity result.
    options : ElasticitySurfaceOptions or None, optional
        Explicit sampling options. Supplying options always requests a fresh
        transient calculation, allowing the plot workflow to override the
        stored grid.
    properties : iterable of str or None, optional
        Optional property subset used when persisted surfaces are selected.
    progress_callback : callable or None, optional
        Callback receiving ``(current, total)`` during transient sampling.

    Returns
    -------
    ElasticitySurfaceCollection
        Persisted or transient numerical surfaces.
    """
    result = _resolve_result(source)
    selected = tuple(properties) if properties is not None else None
    if options is None and result is not None and result.properties_3d is not None:
        available = {
            surface.property_name for surface in result.properties_3d.surfaces.values()
        }
        if selected is None or set(selected).issubset(available):
            return select_elasticity_surfaces(result.properties_3d, selected)

    if options is None:
        options = ElasticitySurfaceOptions(
            properties=(
                ElasticitySurfaceOptions().properties if selected is None else selected
            )
        )
    elif properties is not None:
        options = ElasticitySurfaceOptions(
            ntheta=options.ntheta,
            nphi=options.nphi,
            properties=tuple(properties),
            batch_size=options.batch_size,
        )
    return calculate_elasticity_surfaces(
        source,
        options=options,
        progress_callback=progress_callback,
    )


def select_elasticity_surfaces(
    collection: ElasticitySurfaceCollection,
    properties: Iterable[ElasticSurfaceProperty] | None,
) -> ElasticitySurfaceCollection:
    """Return a property-filtered view of a surface collection.

    Parameters
    ----------
    collection : ElasticitySurfaceCollection
        Source collection.
    properties : iterable of str or None
        Parent property names to retain. ``None`` retains every surface.

    Returns
    -------
    ElasticitySurfaceCollection
        New collection sharing immutable surface objects with the source.
    """
    if properties is None:
        selected = None
    else:
        selected = tuple(dict.fromkeys(properties))
    surfaces = {
        key: surface
        for key, surface in collection.surfaces.items()
        if selected is None or surface.property_name in selected
    }
    return ElasticitySurfaceCollection(
        surfaces=surfaces,
        warnings=list(collection.warnings),
        metadata=dict(collection.metadata),
    )


def _resolve_result(
    source: ElasticTensor | ElasticityResult | ResultData,
) -> ElasticityResult | None:
    """Return an elasticity result payload when one is available."""
    if isinstance(source, ElasticTensor):
        return None
    if isinstance(source, ResultData):
        payload = source.results.get("elasticity")
        if source.metadata.module != "elasticity" or not isinstance(
            payload, ElasticityResult
        ):
            raise ValueError("ResultData does not contain a valid elasticity result.")
        return payload
    return source


def _resolve_tensor(
    source: ElasticTensor | ElasticityResult | ResultData,
) -> ElasticTensor:
    """Return a specialized elastic tensor from a supported source object."""
    if isinstance(source, ElasticTensor):
        return source
    result = _resolve_result(source)
    assert result is not None
    if result.stiffness is None:
        raise ValueError("The elasticity result does not contain a stiffness matrix.")
    if result.stability is not None and not result.stability.is_stable:
        raise ValueError(
            "Three-dimensional surfaces require a mechanically stable stiffness matrix."
        )
    tensor = ElasticTensor(result.stiffness)
    return specialize_elastic_tensor(tensor, result.crystal_system or "triclinic")
