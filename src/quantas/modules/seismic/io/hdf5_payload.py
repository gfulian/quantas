# -*- coding: utf-8 -*-

"""HDF5 payload readers and writers for seismic workflow results.

This module owns only the SEISMIC-specific scientific payload stored under
``/results`` plus module-local diagnostics metadata. Generic Quantas metadata,
input data, options, warnings and events are handled by :mod:`quantas.io.hdf5`.
"""

from __future__ import annotations

from typing import Any

import h5py
import numpy as np
from numpy.typing import NDArray

from quantas.core.geometry import Hemisphere, SphericalGrid, spherical_direction
from quantas.core.physics.elasticity import (
    ElasticAverages,
    IsotropicElasticProperties,
    StabilityResult,
)
from quantas.core.physics.seismic import IsotropicSeismicVelocities
from quantas.core.physics.seismic import (
    EnhancementFieldResult,
    GroupFieldResult,
    PhaseFieldResult,
    PolarizationTrackingResult,
    SamplingLevel,
    SeismicFieldResult,
)
from quantas.io.hdf5 import (
    decode_text,
    read_array_dataset,
    read_mapping,
    require_attr,
    require_group,
    require_text_attr,
    require_unit,
    write_array_dataset,
    write_numeric_attribute,
    write_mapping,
)
from quantas.modules.seismic.models import SeismicResult


_MODE_ORDER = "v_s2,v_s1,v_p"
_PAIR_ORDER = "v_s2-v_s1,v_s1-v_p"
_BRANCH_ORDER = "shear_a,shear_b,p"


def write_seismic_payload(h5: h5py.File, result: SeismicResult) -> h5py.Group:
    """Write the module-specific result tree."""
    group = h5.create_group("results")
    group.attrs["jobname"] = result.jobname
    write_numeric_attribute(group, "density", result.density)
    group.attrs["density_unit"] = "kg m^-3"
    group.attrs["level"] = result.field.level.value
    group.attrs["batch_size"] = result.field.batch_size
    _write_array(group, "stiffness", result.stiffness, unit="GPa")
    _write_stability(group, result)
    _write_averages(group, result)
    _write_isotropic_reference(group, result)
    _write_grid(group, result)
    _write_fields(group, result)
    return group


def write_seismic_diagnostics(
    diagnostics: h5py.Group, result: SeismicResult
) -> h5py.Group:
    """Write SEISMIC-specific diagnostics metadata.

    Parameters
    ----------
    diagnostics : h5py.Group
        Shared diagnostics group created by the HDF5 envelope.
    result : SeismicResult
        SEISMIC result containing module-level metadata.

    Returns
    -------
    h5py.Group
        Created ``diagnostics/seismic`` group.
    """
    seismic_diagnostics = diagnostics.create_group("seismic")
    write_mapping(seismic_diagnostics, result.metadata)
    return seismic_diagnostics


def _write_stability(group: h5py.Group, result: SeismicResult) -> None:
    """Write the stiffness positive-definiteness diagnostic."""
    item = group.create_group("stability")
    item.attrs["is_stable"] = result.stability.is_stable
    write_numeric_attribute(item, "tolerance", result.stability.tolerance)
    _write_array(
        item,
        "eigenvalues",
        result.stability.eigenvalues,
        unit="GPa",
    )


def _write_averages(group: h5py.Group, result: SeismicResult) -> None:
    """Write Voigt, Reuss, and Hill elastic estimates."""
    item = group.create_group("averages")
    item.attrs["columns"] = "bulk_modulus,young_modulus,shear_modulus,poisson_ratio"
    item.attrs["modulus_unit"] = "GPa"
    item.attrs["poisson_ratio_unit"] = "dimensionless"
    for name in ("voigt", "reuss", "hill"):
        _write_array(item, name, getattr(result.averages, name).as_array())


def _write_isotropic_reference(
    group: h5py.Group,
    result: SeismicResult,
) -> None:
    """Write Hill-average isotropic reference velocities."""
    item = group.create_group("isotropic_reference")
    write_numeric_attribute(item, "v_s", result.isotropic_velocities.shear)
    write_numeric_attribute(item, "v_p", result.isotropic_velocities.compressional)
    item.attrs["unit"] = "km s^-1"


def _write_grid(group: h5py.Group, result: SeismicResult) -> None:
    """Write angular coordinates and Cartesian directions."""
    grid = result.grid
    item = group.create_group("grid")
    item.attrs["hemisphere"] = grid.hemisphere.value
    item.attrs["ntheta"] = grid.shape[0]
    item.attrs["nphi"] = grid.shape[1]
    item.attrs["azimuth_endpoint_included"] = False
    _write_array(item, "theta", grid.theta, unit="rad")
    _write_array(item, "phi", grid.phi, unit="rad")
    _write_array(item, "directions", grid.directions, unit="dimensionless")


def _write_fields(group: h5py.Group, result: SeismicResult) -> None:
    """Write phase, group, enhancement, and tracking arrays."""
    fields = group.create_group("fields")
    _write_phase(fields, result)
    if result.field.group is not None:
        _write_group(fields, result)
    if result.field.enhancement is not None:
        _write_enhancement(fields, result)
    if result.field.tracking is not None:
        _write_tracking(fields, result)


def _write_phase(fields: h5py.Group, result: SeismicResult) -> None:
    """Write sampled phase-wave quantities."""
    phase = result.field.phase
    item = fields.create_group("phase")
    item.attrs["mode_order"] = _MODE_ORDER
    item.attrs["mode_symbols"] = "V_S2,V_S1,V_P"
    item.attrs["pair_order"] = _PAIR_ORDER
    _write_array(item, "eigenvalues", phase.eigenvalues, unit="km^2 s^-2")
    _write_array(item, "phase_speeds", phase.phase_speeds, unit="km s^-1")
    _write_array(
        item,
        "polarizations",
        phase.polarizations,
        unit="dimensionless",
    )
    _write_array(
        item,
        "mode_eigenvalue_gaps",
        phase.mode_eigenvalue_gaps,
        unit="km^2 s^-2",
    )
    _write_array(
        item,
        "mode_relative_eigenvalue_gaps",
        phase.mode_relative_eigenvalue_gaps,
        unit="dimensionless",
    )
    _write_array(
        item,
        "pair_eigenvalue_gaps",
        phase.pair_eigenvalue_gaps,
        unit="km^2 s^-2",
    )
    _write_array(
        item,
        "pair_relative_eigenvalue_gaps",
        phase.pair_relative_eigenvalue_gaps,
        unit="dimensionless",
    )
    _write_array(item, "valid_mask", phase.valid_mask)
    _write_array(item, "clamped_mask", phase.clamped_mask)
    _write_array(item, "degeneracy_mask", phase.degeneracy_mask)
    _write_array(
        item,
        "pair_degeneracy_mask",
        phase.pair_degeneracy_mask,
    )
    _write_array(
        item,
        "eigenvalue_thresholds",
        phase.eigenvalue_thresholds,
        unit="km^2 s^-2",
    )
    _write_array(
        item,
        "degeneracy_thresholds",
        phase.degeneracy_thresholds,
        unit="km^2 s^-2",
    )


def _write_group(fields: h5py.Group, result: SeismicResult) -> None:
    """Write sampled group-velocity quantities."""
    assert result.field.group is not None
    group = result.field.group
    item = fields.create_group("group")
    item.attrs["mode_order"] = _MODE_ORDER
    item.attrs["mode_symbols"] = "V_S2,V_S1,V_P"
    _write_array(
        item,
        "eigenvalue_gradients",
        group.eigenvalue_gradients,
        unit="km^2 s^-2",
    )
    _write_array(
        item,
        "group_velocities",
        group.group_velocities,
        unit="km s^-1",
    )
    _write_array(item, "group_speeds", group.group_speeds, unit="km s^-1")
    _write_array(
        item,
        "ray_directions",
        group.ray_directions,
        unit="dimensionless",
    )
    _write_array(
        item,
        "power_flow_angles",
        group.power_flow_angles,
        unit="rad",
    )
    _write_array(item, "valid_mask", group.valid_mask)
    _write_array(item, "resolved_mask", group.resolved_mask)


def _write_enhancement(fields: h5py.Group, result: SeismicResult) -> None:
    """Write sampled curvature and logarithmic enhancement quantities."""
    assert result.field.enhancement is not None
    enhancement = result.field.enhancement
    item = fields.create_group("enhancement")
    item.attrs["mode_order"] = _MODE_ORDER
    item.attrs["mode_symbols"] = "V_S2,V_S1,V_P"
    _write_array(
        item,
        "eigenvalue_hessians",
        enhancement.eigenvalue_hessians,
        unit="km^2 s^-2",
    )
    _write_array(
        item,
        "ray_direction_gradients",
        enhancement.ray_direction_gradients,
        unit="dimensionless",
    )
    _write_array(
        item,
        "area_factors",
        enhancement.area_factors,
        unit="dimensionless",
    )
    _write_array(
        item,
        "caustic_thresholds",
        enhancement.caustic_thresholds,
        unit="dimensionless",
    )
    dataset = _write_array(
        item,
        "log10_enhancement",
        enhancement.log10_enhancement,
        unit="dimensionless",
    )
    dataset.attrs["representation"] = "log10(A)"
    _write_array(item, "valid_mask", enhancement.valid_mask)
    _write_array(item, "resolved_mask", enhancement.resolved_mask)
    _write_array(item, "finite_mask", enhancement.finite_mask)
    _write_array(
        item,
        "caustic_candidate_mask",
        enhancement.caustic_candidate_mask,
    )


def _write_tracking(fields: h5py.Group, result: SeismicResult) -> None:
    """Write tracked axial-polarization branches."""
    assert result.field.tracking is not None
    tracking = result.field.tracking
    item = fields.create_group("tracking")
    item.attrs["branch_order"] = _BRANCH_ORDER
    _write_array(
        item,
        "polarizations",
        tracking.polarizations,
        unit="dimensionless",
    )
    _write_array(
        item,
        "branch_mode_indices",
        tracking.branch_mode_indices,
    )
    _write_array(item, "sign_flip_mask", tracking.sign_flip_mask)
    _write_array(
        item,
        "continuity_scores",
        tracking.continuity_scores,
        unit="dimensionless",
    )
    _write_array(item, "resolved_mask", tracking.resolved_mask)
    _write_array(
        item,
        "segment_start_mask",
        tracking.segment_start_mask,
    )
    _write_array(item, "shear_swap_mask", tracking.shear_swap_mask)
    _write_array(
        item,
        "shear_permutation_ambiguous_mask",
        tracking.shear_permutation_ambiguous_mask,
    )
    _write_array(
        item,
        "shear_subspace_rotation_mask",
        tracking.shear_subspace_rotation_mask,
    )


def _write_array(
    group: h5py.Group,
    key: str,
    value: Any,
    *,
    unit: str | None = None,
) -> h5py.Dataset:
    """Write a seismic payload dataset using shared HDF5 mechanics."""
    return write_array_dataset(group, key, value, unit=unit, compression=True)


def read_seismic_payload(h5: h5py.File) -> SeismicResult:
    """Reconstruct the module-specific seismic result payload."""
    group = h5["results"]
    require_attr(group, "jobname")
    require_attr(group, "density")
    require_attr(group, "level")
    require_attr(group, "batch_size")
    require_text_attr(group, "density_unit", "kg m^-3")

    density = float(group.attrs["density"])
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("Stored seismic density must be finite and positive.")
    stiffness = _read_array(group, "stiffness", dtype=float, shape=(6, 6))
    if not np.all(np.isfinite(stiffness)):
        raise ValueError("Stored seismic stiffness must contain finite values.")
    require_unit(group["stiffness"], "GPa")
    stability = _read_stability(group)
    averages = _read_averages(group)
    isotropic = _read_isotropic_reference(group)
    grid = _read_grid(group)
    field = _read_field(group, grid)

    metadata: dict[str, Any] = {}
    diagnostics = h5.get("diagnostics")
    if diagnostics is not None and "seismic" in diagnostics:
        metadata = read_mapping(diagnostics["seismic"])

    return SeismicResult(
        jobname=decode_text(group.attrs["jobname"]),
        density=density,
        stiffness=stiffness,
        stability=stability,
        averages=averages,
        isotropic_velocities=isotropic,
        grid=grid,
        field=field,
        metadata=metadata,
    )


def _read_stability(group: h5py.Group) -> StabilityResult:
    """Read the positive-definiteness diagnostic."""
    item = require_group(group, "stability")
    values = _read_array(item, "eigenvalues", dtype=float, shape=(6,))
    require_unit(item["eigenvalues"], "GPa")
    return StabilityResult(
        is_stable=bool(require_attr(item, "is_stable")),
        eigenvalues=values,
        tolerance=float(require_attr(item, "tolerance")),
    )


def _read_averages(group: h5py.Group) -> ElasticAverages:
    """Read Voigt, Reuss, and Hill isotropic estimates."""
    item = require_group(group, "averages")
    return ElasticAverages(
        voigt=_properties_from_array(_read_array(item, "voigt", float, (4,))),
        reuss=_properties_from_array(_read_array(item, "reuss", float, (4,))),
        hill=_properties_from_array(_read_array(item, "hill", float, (4,))),
    )


def _properties_from_array(values: NDArray[np.float64]) -> IsotropicElasticProperties:
    """Construct isotropic properties from ``K, E, G, nu`` values."""
    return IsotropicElasticProperties(
        bulk_modulus=float(values[0]),
        young_modulus=float(values[1]),
        shear_modulus=float(values[2]),
        poisson_ratio=float(values[3]),
    )


def _read_isotropic_reference(group: h5py.Group) -> IsotropicSeismicVelocities:
    """Read isotropic shear and compressional reference velocities."""
    item = require_group(group, "isotropic_reference")
    require_text_attr(item, "unit", "km s^-1")
    return IsotropicSeismicVelocities(
        shear=float(require_attr(item, "v_s")),
        compressional=float(require_attr(item, "v_p")),
    )


def _read_grid(group: h5py.Group) -> SphericalGrid:
    """Read the regular spherical sampling grid."""
    item = require_group(group, "grid")
    domain = Hemisphere(decode_text(require_attr(item, "hemisphere")))
    theta = _read_array(item, "theta", float)
    phi = _read_array(item, "phi", float)
    directions = _read_array(
        item,
        "directions",
        float,
        shape=(theta.size, phi.size, 3),
    )
    require_unit(item["theta"], "rad")
    require_unit(item["phi"], "rad")
    require_unit(item["directions"], "dimensionless")
    if int(require_attr(item, "ntheta")) != theta.size:
        raise ValueError("Stored ntheta does not match the theta dataset.")
    if int(require_attr(item, "nphi")) != phi.size:
        raise ValueError("Stored nphi does not match the phi dataset.")
    if bool(require_attr(item, "azimuth_endpoint_included")):
        raise ValueError("The scientific azimuth grid must not duplicate 2π.")
    if theta.size < 2 or phi.size < 3:
        raise ValueError("Stored spherical grid dimensions are too small.")
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(phi)):
        raise ValueError("Stored spherical angles must be finite.")
    if not np.all(np.diff(theta) > 0.0) or not np.all(np.diff(phi) > 0.0):
        raise ValueError("Stored spherical coordinates must be strictly increasing.")
    if phi[0] < -1.0e-14 or phi[-1] >= 2.0 * np.pi:
        raise ValueError("Stored azimuths must lie in [0, 2π).")
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    expected = spherical_direction(theta_grid, phi_grid)
    if not np.allclose(directions, expected, rtol=0.0, atol=2.0e-12):
        raise ValueError("Stored Cartesian directions do not match theta and phi.")
    return SphericalGrid(
        theta=_readonly(theta),
        phi=_readonly(phi),
        theta_grid=_readonly(theta_grid),
        phi_grid=_readonly(phi_grid),
        directions=_readonly(directions),
        hemisphere=domain,
    )


def _read_field(group: h5py.Group, grid: SphericalGrid) -> SeismicFieldResult:
    """Read phase, group, enhancement, and tracking fields."""
    fields = require_group(group, "fields")
    level = SamplingLevel(decode_text(require_attr(group, "level")))
    batch_size = int(require_attr(group, "batch_size"))
    if batch_size < 1:
        raise ValueError("Stored batch_size must be a positive integer.")
    phase = _read_phase_field(require_group(fields, "phase"), grid)
    group_field = (
        _read_group_field(require_group(fields, "group"), phase)
        if level in (SamplingLevel.GROUP, SamplingLevel.ENHANCEMENT)
        else None
    )
    enhancement = (
        _read_enhancement_field(
            require_group(fields, "enhancement"),
            group_field,
        )
        if level is SamplingLevel.ENHANCEMENT and group_field is not None
        else None
    )
    tracking = (
        _read_tracking(require_group(fields, "tracking"), phase)
        if "tracking" in fields
        else None
    )
    if level is SamplingLevel.PHASE and "group" in fields:
        raise ValueError("Phase-only results must not contain a group field.")
    if level is not SamplingLevel.ENHANCEMENT and "enhancement" in fields:
        raise ValueError("The stored sampling level does not include enhancement.")
    return SeismicFieldResult(
        grid=grid,
        level=level,
        phase=phase,
        group=group_field,
        enhancement=enhancement,
        tracking=tracking,
        batch_size=batch_size,
    )


def _read_phase_field(group: h5py.Group, grid: SphericalGrid) -> PhaseFieldResult:
    """Read sampled phase quantities."""
    require_text_attr(group, "mode_order", _MODE_ORDER)
    require_text_attr(group, "mode_symbols", "V_S2,V_S1,V_P")
    require_text_attr(group, "pair_order", _PAIR_ORDER)
    n = grid.size
    result = PhaseFieldResult(
        directions=grid.flat_directions,
        eigenvalues=_read_array(group, "eigenvalues", float, (n, 3)),
        phase_speeds=_read_array(group, "phase_speeds", float, (n, 3)),
        polarizations=_read_array(group, "polarizations", float, (n, 3, 3)),
        mode_eigenvalue_gaps=_read_array(group, "mode_eigenvalue_gaps", float, (n, 3)),
        mode_relative_eigenvalue_gaps=_read_array(
            group, "mode_relative_eigenvalue_gaps", float, (n, 3)
        ),
        pair_eigenvalue_gaps=_read_array(group, "pair_eigenvalue_gaps", float, (n, 2)),
        pair_relative_eigenvalue_gaps=_read_array(
            group, "pair_relative_eigenvalue_gaps", float, (n, 2)
        ),
        valid_mask=_read_array(group, "valid_mask", bool, (n, 3)),
        clamped_mask=_read_array(group, "clamped_mask", bool, (n, 3)),
        degeneracy_mask=_read_array(group, "degeneracy_mask", bool, (n, 3)),
        pair_degeneracy_mask=_read_array(group, "pair_degeneracy_mask", bool, (n, 2)),
        eigenvalue_thresholds=_read_array(group, "eigenvalue_thresholds", float, (n,)),
        degeneracy_thresholds=_read_array(group, "degeneracy_thresholds", float, (n,)),
    )
    require_unit(group["eigenvalues"], "km^2 s^-2")
    require_unit(group["phase_speeds"], "km s^-1")
    require_unit(group["polarizations"], "dimensionless")
    for key in (
        "mode_eigenvalue_gaps",
        "pair_eigenvalue_gaps",
        "eigenvalue_thresholds",
        "degeneracy_thresholds",
    ):
        require_unit(group[key], "km^2 s^-2")
    for key in ("mode_relative_eigenvalue_gaps", "pair_relative_eigenvalue_gaps"):
        require_unit(group[key], "dimensionless")
    return result


def _read_group_field(
    group: h5py.Group,
    phase: PhaseFieldResult,
) -> GroupFieldResult:
    """Read sampled group-velocity quantities."""
    require_text_attr(group, "mode_order", _MODE_ORDER)
    require_text_attr(group, "mode_symbols", "V_S2,V_S1,V_P")
    n = phase.n_points
    result = GroupFieldResult(
        phase=phase,
        eigenvalue_gradients=_read_array(
            group, "eigenvalue_gradients", float, (n, 3, 3)
        ),
        group_velocities=_read_array(group, "group_velocities", float, (n, 3, 3)),
        group_speeds=_read_array(group, "group_speeds", float, (n, 3)),
        ray_directions=_read_array(group, "ray_directions", float, (n, 3, 3)),
        power_flow_angles=_read_array(group, "power_flow_angles", float, (n, 3)),
        valid_mask=_read_array(group, "valid_mask", bool, (n, 3)),
        resolved_mask=_read_array(group, "resolved_mask", bool, (n, 3)),
    )
    require_unit(group["group_velocities"], "km s^-1")
    require_unit(group["group_speeds"], "km s^-1")
    require_unit(group["power_flow_angles"], "rad")
    require_unit(group["eigenvalue_gradients"], "km^2 s^-2")
    require_unit(group["ray_directions"], "dimensionless")
    return result


def _read_enhancement_field(
    group: h5py.Group,
    group_field: GroupFieldResult,
) -> EnhancementFieldResult:
    """Read curvature and logarithmic enhancement quantities."""
    require_text_attr(group, "mode_order", _MODE_ORDER)
    require_text_attr(group, "mode_symbols", "V_S2,V_S1,V_P")
    n = group_field.n_points
    log10_enhancement = _read_array(group, "log10_enhancement", float, (n, 3))
    require_text_attr(group["log10_enhancement"], "representation", "log10(A)")
    for key in (
        "ray_direction_gradients",
        "area_factors",
        "caustic_thresholds",
        "log10_enhancement",
    ):
        require_unit(group[key], "dimensionless")
    require_unit(group["eigenvalue_hessians"], "km^2 s^-2")
    with np.errstate(over="ignore", invalid="ignore"):
        enhancement = np.power(10.0, log10_enhancement)
    enhancement = _readonly(enhancement)
    return EnhancementFieldResult(
        group=group_field,
        eigenvalue_hessians=_read_array(
            group, "eigenvalue_hessians", float, (n, 3, 3, 3)
        ),
        ray_direction_gradients=_read_array(
            group, "ray_direction_gradients", float, (n, 3, 3, 3)
        ),
        area_factors=_read_array(group, "area_factors", float, (n, 3)),
        caustic_thresholds=_read_array(group, "caustic_thresholds", float, (n, 3)),
        enhancement=enhancement,
        log10_enhancement=log10_enhancement,
        valid_mask=_read_array(group, "valid_mask", bool, (n, 3)),
        resolved_mask=_read_array(group, "resolved_mask", bool, (n, 3)),
        finite_mask=_read_array(group, "finite_mask", bool, (n, 3)),
        caustic_candidate_mask=_read_array(
            group, "caustic_candidate_mask", bool, (n, 3)
        ),
    )


def _read_tracking(
    group: h5py.Group,
    phase: PhaseFieldResult,
) -> PolarizationTrackingResult:
    """Read tracked axial-polarization branches."""
    require_text_attr(group, "branch_order", _BRANCH_ORDER)
    require_unit(group["polarizations"], "dimensionless")
    require_unit(group["continuity_scores"], "dimensionless")
    n = phase.n_points
    return PolarizationTrackingResult(
        directions=phase.directions,
        polarizations=_read_array(group, "polarizations", float, (n, 3, 3)),
        branch_mode_indices=_read_array(group, "branch_mode_indices", np.int64, (n, 3)),
        sign_flip_mask=_read_array(group, "sign_flip_mask", bool, (n, 3)),
        continuity_scores=_read_array(group, "continuity_scores", float, (n, 3)),
        resolved_mask=_read_array(group, "resolved_mask", bool, (n, 3)),
        segment_start_mask=_read_array(group, "segment_start_mask", bool, (n, 3)),
        shear_swap_mask=_read_array(group, "shear_swap_mask", bool, (n,)),
        shear_permutation_ambiguous_mask=_read_array(
            group, "shear_permutation_ambiguous_mask", bool, (n,)
        ),
        shear_subspace_rotation_mask=_read_array(
            group, "shear_subspace_rotation_mask", bool, (n,)
        ),
    )


def _read_array(
    group: h5py.Group,
    key: str,
    dtype: Any,
    shape: tuple[int, ...] | None = None,
) -> NDArray[Any]:
    """Read a required numerical array with shared HDF5 validation."""
    array = read_array_dataset(
        group,
        key,
        dtype=dtype,
        shape=shape,
        readonly=True,
    )
    assert array is not None
    return array


def _readonly(array: NDArray[Any]) -> NDArray[Any]:
    """Return a read-only NumPy copy."""
    result = np.array(array, copy=True)
    result.setflags(write=False)
    return result
