# -*- coding: utf-8 -*-

"""Read and write the elasticity payload in native Quantas HDF5 files.

Generic metadata, normalized inputs, options, diagnostics, warnings, and events
remain the responsibility of :mod:`quantas.io.hdf5`.
"""

from __future__ import annotations

import h5py
import numpy as np

from quantas.core.physics.elasticity import (
    DirectionalExtrema,
    DirectionalSurface,
    ElasticAverages,
    ElasticitySurfaceCollection,
    IsotropicElasticProperties,
    StabilityResult,
)
from quantas.io.hdf5 import (
    decode_text,
    read_array_dataset,
    read_group_mapping,
    write_array_dataset,
    write_numeric_attribute,
    write_mapping,
)
from quantas.modules.elasticity.models import ElasticityResult


def write_elasticity_payload(h5: h5py.File, result: ElasticityResult) -> h5py.Group:
    """Write the elasticity-specific result payload.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file.
    result : ElasticityResult
        Passive elasticity result object.

    Returns
    -------
    h5py.Group
        Created ``results`` group.
    """
    group = h5.create_group("results")
    group.attrs["jobname"] = result.jobname
    group.attrs["crystal_system"] = result.crystal_system or "unknown"
    write_array_dataset(group, "stiffness", result.stiffness)
    write_array_dataset(group, "compliance", result.compliance)

    if result.averages is not None:
        averages = group.create_group("averages")
        for name in ("voigt", "reuss", "hill"):
            write_array_dataset(
                averages, name, getattr(result.averages, name).as_array()
            )

    if result.stability is not None:
        stability = group.create_group("stability")
        stability.attrs["is_stable"] = result.stability.is_stable
        write_numeric_attribute(stability, "tolerance", result.stability.tolerance)
        write_array_dataset(stability, "eigenvalues", result.stability.eigenvalues)

    variations = group.create_group("variations")
    for name, variation in result.variations.items():
        item = variations.create_group(name)
        write_numeric_attribute(item, "minimum", variation.minimum)
        write_numeric_attribute(item, "maximum", variation.maximum)
        write_numeric_attribute(item, "anisotropy", variation.anisotropy)
        write_array_dataset(item, "minimum_axis", variation.minimum_axis)
        write_array_dataset(item, "maximum_axis", variation.maximum_axis)
        if variation.minimum_measurement_axis is not None:
            write_array_dataset(
                item,
                "minimum_measurement_axis",
                variation.minimum_measurement_axis,
            )
        if variation.maximum_measurement_axis is not None:
            write_array_dataset(
                item,
                "maximum_measurement_axis",
                variation.maximum_measurement_axis,
            )

    properties_2d = group.create_group("properties_2d")
    for plane, plane_data in result.properties_2d.items():
        plane_group = properties_2d.create_group(plane)
        for key, value in plane_data.items():
            write_array_dataset(plane_group, key, value, compression=True)

    if result.properties_3d is not None:
        _write_surface_collection(group, result.properties_3d)

    metadata = group.create_group("metadata")
    write_mapping(metadata, result.metadata)
    return group


def read_elasticity_payload(group: h5py.Group) -> ElasticityResult:
    """Reconstruct the elasticity-specific payload from ``/results``.

    Parameters
    ----------
    group : h5py.Group
        HDF5 ``results`` group.

    Returns
    -------
    ElasticityResult
        Reconstructed passive elasticity result model.
    """
    result = ElasticityResult(
        jobname=decode_text(group.attrs.get("jobname", "Unknown")),
        crystal_system=decode_text(group.attrs.get("crystal_system", "unknown")),
        stiffness=_read_array(group, "stiffness"),
        compliance=_read_array(group, "compliance"),
        averages=_read_averages(group),
        stability=_read_stability(group),
    )

    if "variations" in group:
        for name, item in group["variations"].items():
            result.add_variation(
                name,
                DirectionalExtrema(
                    minimum=float(item.attrs["minimum"]),
                    maximum=float(item.attrs["maximum"]),
                    anisotropy=float(item.attrs["anisotropy"]),
                    minimum_axis=_required_axis(item, "minimum_axis"),
                    maximum_axis=_required_axis(item, "maximum_axis"),
                    minimum_measurement_axis=_optional_axis(
                        _read_array(item, "minimum_measurement_axis")
                    ),
                    maximum_measurement_axis=_optional_axis(
                        _read_array(item, "maximum_measurement_axis")
                    ),
                ),
            )

    data_group = group.get("properties_2d", group.get("polar"))
    if data_group is not None:
        for plane, plane_group in data_group.items():
            for key, value in plane_group.items():
                result.add_2d_data(plane, key, value[()])

    if "properties_3d" in group:
        result.properties_3d = _read_surface_collection(group["properties_3d"])

    if "metadata" in group:
        result.metadata.update(read_group_mapping(group["metadata"]))
    return result


def _write_surface_collection(
    group: h5py.Group,
    collection: ElasticitySurfaceCollection,
) -> None:
    """Write one compressed three-dimensional surface collection."""
    item = group.create_group("properties_3d")
    item.attrs["schema_version"] = "1.0"
    metadata = item.create_group("metadata")
    write_mapping(metadata, collection.metadata)
    write_array_dataset(
        item,
        "warnings",
        np.asarray(collection.warnings, dtype=object),
    )
    if not collection.surfaces:
        item.create_group("surfaces")
        return

    first = next(iter(collection.surfaces.values()))
    write_array_dataset(item, "theta", first.theta, unit="radian", compression=True)
    write_array_dataset(item, "phi", first.phi, unit="radian", compression=True)
    surfaces = item.create_group("surfaces")
    for key, surface in collection.surfaces.items():
        if (
            surface.theta.shape != first.theta.shape
            or surface.phi.shape != first.phi.shape
            or not np.allclose(surface.theta, first.theta, rtol=0.0, atol=0.0)
            or not np.allclose(surface.phi, first.phi, rtol=0.0, atol=0.0)
        ):
            raise ValueError("All persisted elasticity surfaces must share one grid.")
        branch = surfaces.create_group(key)
        branch.attrs["property_name"] = surface.property_name
        branch.attrs["branch"] = surface.branch
        branch.attrs["unit"] = surface.unit
        write_array_dataset(branch, "radius", surface.radius, compression=True)
        write_array_dataset(branch, "values", surface.values, compression=True)
        branch_metadata = branch.create_group("metadata")
        write_mapping(branch_metadata, surface.metadata)


def _read_surface_collection(group: h5py.Group) -> ElasticitySurfaceCollection:
    """Read a compressed three-dimensional surface collection."""
    metadata = read_group_mapping(group["metadata"]) if "metadata" in group else {}
    warnings: list[str] = []
    if "warnings" in group:
        warnings = [decode_text(value) for value in group["warnings"][()]]
    collection = ElasticitySurfaceCollection(
        warnings=warnings,
        metadata=metadata,
    )
    surfaces_group = group.get("surfaces")
    if surfaces_group is None or len(surfaces_group) == 0:
        return collection
    theta = np.asarray(group["theta"][()], dtype=float)
    phi = np.asarray(group["phi"][()], dtype=float)
    if theta.shape != phi.shape or theta.ndim != 2:
        raise ValueError("Stored elasticity 3D angular grids are inconsistent.")
    directions = np.stack(
        (
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ),
        axis=-1,
    )
    for key, branch in surfaces_group.items():
        radius = np.asarray(branch["radius"][()], dtype=float)
        values = np.asarray(branch["values"][()], dtype=float)
        if radius.shape != theta.shape or values.shape != theta.shape:
            raise ValueError(
                f"Stored elasticity 3D surface '{key}' has an invalid shape."
            )
        branch_metadata = (
            read_group_mapping(branch["metadata"]) if "metadata" in branch else {}
        )
        collection.surfaces[key] = DirectionalSurface(
            key=key,
            property_name=decode_text(branch.attrs["property_name"]),
            branch=decode_text(branch.attrs["branch"]),
            unit=decode_text(branch.attrs.get("unit", "")),
            theta=theta.copy(),
            phi=phi.copy(),
            radius=radius,
            values=values,
            x=radius * directions[..., 0],
            y=radius * directions[..., 1],
            z=radius * directions[..., 2],
            metadata=branch_metadata,
        )
    return collection


def _read_averages(group: h5py.Group) -> ElasticAverages | None:
    """Read structured Voigt-Reuss-Hill averages."""
    if "averages" not in group:
        return None
    averages = group["averages"]
    return ElasticAverages(
        voigt=_properties_from_array(averages["voigt"][()]),
        reuss=_properties_from_array(averages["reuss"][()]),
        hill=_properties_from_array(averages["hill"][()]),
    )


def _properties_from_array(values: np.ndarray) -> IsotropicElasticProperties:
    """Build isotropic elastic properties from ``K, E, G, nu`` values."""
    return IsotropicElasticProperties(
        bulk_modulus=float(values[0]),
        young_modulus=float(values[1]),
        shear_modulus=float(values[2]),
        poisson_ratio=float(values[3]),
    )


def _read_stability(group: h5py.Group) -> StabilityResult | None:
    """Read a positive-definiteness result."""
    if "stability" not in group:
        return None
    stability = group["stability"]
    return StabilityResult(
        is_stable=bool(stability.attrs["is_stable"]),
        eigenvalues=np.asarray(stability["eigenvalues"][()], dtype=float),
        tolerance=float(stability.attrs.get("tolerance", 0.0)),
    )


def _read_array(group: h5py.Group, key: str) -> np.ndarray | None:
    """Read an optional HDF5 numerical dataset."""
    value = read_array_dataset(group, key, required=False)
    return None if value is None else np.asarray(value)


def _required_axis(group: h5py.Group, key: str) -> list[float]:
    """Read a required Cartesian direction."""
    value = _read_array(group, key)
    if value is None:
        raise ValueError(f"Missing required direction dataset: {key}")
    return value.tolist()


def _optional_axis(value: np.ndarray | None) -> list[float] | None:
    """Convert an optional direction array to a list."""
    return None if value is None else value.tolist()
