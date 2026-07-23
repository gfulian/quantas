# -*- coding: utf-8 -*-

"""HDF5 payload for quasi-static thermoelastic results."""

from __future__ import annotations

from typing import Any, cast

import h5py
import numpy as np

from quantas.core.math.fitting import fit_result_from_dict
from quantas.core.physics.elasticity import StabilityFieldResult
from quantas.io.hdf5 import (
    decode_text,
    read_group_mapping,
    write_array_dataset,
    write_mapping,
)
from quantas.modules.thermoelasticity.models import (
    ElasticComponentFit,
    ReferenceEOSFit,
    ThermoelasticFitQuality,
    ThermoelasticProfileResult,
    ThermoelasticQualityLevel,
    ThermoelasticResult,
)

THERMOELASTIC_HDF5_SCHEMA = "1.0"


def write_thermoelastic_payload(
    h5: h5py.File,
    result: ThermoelasticResult,
) -> h5py.Group:
    """Write the complete thermoelastic scientific payload and diagnostics.

    Parameters
    ----------
    h5 : h5py.File
        Open Quantas HDF5 envelope.
    result : ThermoelasticResult
        Calibration or analysis result to serialize.

    Returns
    -------
    h5py.Group
        Created ``/results`` group.
    """
    group = h5.create_group("results")
    group.attrs["thermoelastic_schema_version"] = THERMOELASTIC_HDF5_SCHEMA
    group.attrs["jobname"] = result.jobname
    group.attrs["completed"] = bool(result.completed)
    labels = np.asarray(
        result.independent_labels,
        dtype=h5py.string_dtype(encoding="utf-8"),
    )
    group.create_dataset("independent_labels", data=labels)
    _write_result_fields(group, result)
    _write_reference_eos(group, result.reference_eos)
    _write_component_fits(group, result.component_fits)
    _write_profiles(group, result.profiles)
    metadata = group.create_group("metadata")
    write_mapping(metadata, result.metadata)
    return group


def read_thermoelastic_payload(group: h5py.Group) -> ThermoelasticResult:
    """Read a complete thermoelastic result from an HDF5 ``/results`` group.

    Parameters
    ----------
    group : h5py.Group
        Native thermoelastic ``/results`` group.

    Returns
    -------
    ThermoelasticResult
        Reconstructed passive scientific result.

    Raises
    ------
    ValueError
        If the thermoelastic payload schema is unsupported.
    """
    schema = decode_text(
        group.attrs.get("thermoelastic_schema_version", THERMOELASTIC_HDF5_SCHEMA)
    )
    if schema != THERMOELASTIC_HDF5_SCHEMA:
        raise ValueError(f"unsupported thermoelastic HDF5 payload schema: {schema}")
    reference = _read_reference_eos(group)
    components = _read_component_fits(group)
    profiles = _read_profiles(group)
    labels = tuple(
        decode_text(value) for value in np.atleast_1d(group["independent_labels"][()])
    )
    return _read_result_fields(
        group,
        reference=reference,
        components=components,
        profiles=profiles,
        labels=labels,
    )


def _write_result_fields(group: h5py.Group, result: ThermoelasticResult) -> None:
    """Write grid-independent and optional tensor fields."""
    write_array_dataset(group, "temperature", result.temperature, unit="K")
    write_array_dataset(group, "pressure", result.pressure, unit="GPa")
    write_array_dataset(
        group, "equilibrium_volume", result.equilibrium_volume, unit="angstrom^3"
    )
    write_array_dataset(group, "density", result.density, unit="kg m^-3")
    write_array_dataset(group, "elastic_extrapolation", result.extrapolation_mask)
    write_array_dataset(
        group,
        "qha_extrapolation",
        np.asarray(result.qha_extrapolation_mask, dtype=np.bool_),
    )
    optional = {
        "sigma_equilibrium_volume": result.sigma_equilibrium_volume,
        "independent_stiffness": result.independent_stiffness,
        "sigma_independent_stiffness": result.sigma_independent_stiffness,
        "independent_stiffness_covariance": result.independent_stiffness_covariance,
        "stiffness_isothermal": result.stiffness_isothermal,
        "sigma_stiffness_isothermal": result.sigma_stiffness_isothermal,
        "isochoric_heat_capacity_cell": result.isochoric_heat_capacity_cell,
        "sigma_isochoric_heat_capacity_cell": (
            result.sigma_isochoric_heat_capacity_cell
        ),
        "thermal_expansion_tensor": result.thermal_expansion_tensor,
        "sigma_thermal_expansion_tensor": result.sigma_thermal_expansion_tensor,
        "stiffness_adiabatic": result.stiffness_adiabatic,
        "sigma_stiffness_adiabatic": result.sigma_stiffness_adiabatic,
        "adiabatic_correction": result.adiabatic_correction,
        "adiabatic_thermal_stress": result.adiabatic_thermal_stress,
        "adiabatic_valid_mask": (
            None
            if result.adiabatic_valid_mask is None
            else np.asarray(result.adiabatic_valid_mask, dtype=np.bool_)
        ),
    }
    for name, value in optional.items():
        _write_optional_array(group, name, value)
    _write_stability(group, result.stability)


def _write_reference_eos(group: h5py.Group, reference: ReferenceEOSFit) -> None:
    """Write the static reference equation-of-state fit."""
    node = group.create_group("reference_eos")
    node.attrs["eos"] = reference.eos
    node.attrs["reference_volume_A3"] = reference.reference_volume
    node.attrs["bulk_modulus_GPa"] = reference.bulk_modulus
    node.attrs["bulk_modulus_derivative"] = reference.bulk_modulus_derivative
    node.attrs["bulk_modulus_second_derivative_GPa-1"] = (
        reference.bulk_modulus_second_derivative
    )
    _write_optional_array(node, "covariance_V0_K0_Kprime", reference.covariance)
    fit_group = node.create_group("fit")
    write_mapping(fit_group, reference.fit.as_dict())
    metadata = node.create_group("metadata")
    write_mapping(metadata, reference.metadata)


def _write_component_fits(
    group: h5py.Group,
    components: dict[str, ElasticComponentFit],
) -> None:
    """Write independent elastic-component fits and diagnostics."""
    component_group = group.create_group("component_fits")
    for label, component in components.items():
        item = component_group.create_group(label)
        item.attrs["active"] = bool(component.active)
        item.attrs["zero_by_tolerance"] = bool(component.zero_by_tolerance)
        item.attrs["wallace_delta"] = component.wallace_delta
        write_array_dataset(
            item, "entries", np.asarray(component.entries, dtype=np.float64)
        )
        for name in (
            "volumes",
            "pressures",
            "observed",
            "fitted",
            "residuals",
            "relative_residuals",
            "symmetry_spread",
        ):
            write_array_dataset(item, name, getattr(component, name))
        if component.fit is not None:
            component_fit = item.create_group("fit")
            write_mapping(component_fit, component.fit.as_dict())
        metadata = item.create_group("metadata")
        write_mapping(metadata, component.metadata)
        if component.quality is not None:
            quality = item.create_group("quality")
            write_mapping(quality, component.quality.as_dict())


def _write_profiles(
    group: h5py.Group,
    profiles: dict[str, ThermoelasticProfileResult],
) -> None:
    """Write all named geological profile fields."""
    profiles_group = group.create_group("profiles")
    for name, profile in profiles.items():
        item = profiles_group.create_group(_validated_hdf5_name(name))
        for field_name in (
            "depth",
            "pressure",
            "temperature",
            "volume",
            "density",
            "independent_stiffness",
            "sigma_independent_stiffness",
            "independent_stiffness_covariance",
            "stiffness_isothermal",
            "sigma_stiffness_isothermal",
            "qha_extrapolation_mask",
            "elastic_extrapolation_mask",
        ):
            write_array_dataset(
                item, field_name, getattr(profile, field_name), compression=True
            )
        for field_name in (
            "stiffness_adiabatic",
            "sigma_stiffness_adiabatic",
            "adiabatic_correction",
            "adiabatic_thermal_stress",
            "adiabatic_valid_mask",
        ):
            _write_optional_array(item, field_name, getattr(profile, field_name))
        metadata = item.create_group("metadata")
        write_mapping(metadata, profile.metadata)
        _write_stability(item, profile.stability)


def _read_reference_eos(group: h5py.Group) -> ReferenceEOSFit:
    """Read the static reference equation-of-state fit."""
    node = group["reference_eos"]
    assert isinstance(node, h5py.Group)
    fit_group = node["fit"]
    assert isinstance(fit_group, h5py.Group)
    return ReferenceEOSFit(
        eos=decode_text(node.attrs["eos"]),
        reference_volume=float(node.attrs["reference_volume_A3"]),
        bulk_modulus=float(node.attrs["bulk_modulus_GPa"]),
        bulk_modulus_derivative=float(node.attrs["bulk_modulus_derivative"]),
        bulk_modulus_second_derivative=float(
            node.attrs["bulk_modulus_second_derivative_GPa-1"]
        ),
        covariance=_read_optional_array(node, "covariance_V0_K0_Kprime"),
        fit=fit_result_from_dict(read_group_mapping(fit_group)),
        metadata=_read_metadata(node),
    )


def _read_component_fits(group: h5py.Group) -> dict[str, ElasticComponentFit]:
    """Read independent elastic-component fits and diagnostics."""
    components: dict[str, ElasticComponentFit] = {}
    component_group = group["component_fits"]
    assert isinstance(component_group, h5py.Group)
    for label, node in component_group.items():
        assert isinstance(node, h5py.Group)
        entries_array = np.asarray(node["entries"], dtype=np.float64)
        entries = tuple(
            (int(row[0]), int(row[1]), float(row[2])) for row in entries_array
        )
        component_fit = None
        if "fit" in node:
            fit_node = node["fit"]
            assert isinstance(fit_node, h5py.Group)
            component_fit = fit_result_from_dict(read_group_mapping(fit_node))
        components[label] = ElasticComponentFit(
            label=label,
            entries=entries,
            wallace_delta=float(node.attrs["wallace_delta"]),
            volumes=np.asarray(node["volumes"], dtype=np.float64),
            pressures=np.asarray(node["pressures"], dtype=np.float64),
            observed=np.asarray(node["observed"], dtype=np.float64),
            fitted=np.asarray(node["fitted"], dtype=np.float64),
            residuals=np.asarray(node["residuals"], dtype=np.float64),
            relative_residuals=np.asarray(node["relative_residuals"], dtype=np.float64),
            symmetry_spread=np.asarray(node["symmetry_spread"], dtype=np.float64),
            fit=component_fit,
            active=bool(node.attrs["active"]),
            zero_by_tolerance=bool(node.attrs["zero_by_tolerance"]),
            metadata=_read_metadata(node),
            quality=_read_fit_quality(node),
        )
    return components


def _read_profiles(group: h5py.Group) -> dict[str, ThermoelasticProfileResult]:
    """Read all named geological profile fields."""
    profiles: dict[str, ThermoelasticProfileResult] = {}
    if "profiles" not in group:
        return profiles
    profiles_group = group["profiles"]
    assert isinstance(profiles_group, h5py.Group)
    for name, node in profiles_group.items():
        assert isinstance(node, h5py.Group)
        profiles[name] = ThermoelasticProfileResult(
            name=name,
            depth=np.asarray(node["depth"], dtype=np.float64),
            pressure=np.asarray(node["pressure"], dtype=np.float64),
            temperature=np.asarray(node["temperature"], dtype=np.float64),
            volume=np.asarray(node["volume"], dtype=np.float64),
            density=np.asarray(node["density"], dtype=np.float64),
            independent_stiffness=np.asarray(
                node["independent_stiffness"], dtype=np.float64
            ),
            sigma_independent_stiffness=np.asarray(
                node["sigma_independent_stiffness"], dtype=np.float64
            ),
            independent_stiffness_covariance=np.asarray(
                node["independent_stiffness_covariance"], dtype=np.float64
            ),
            stiffness_isothermal=np.asarray(
                node["stiffness_isothermal"], dtype=np.float64
            ),
            sigma_stiffness_isothermal=np.asarray(
                node["sigma_stiffness_isothermal"], dtype=np.float64
            ),
            qha_extrapolation_mask=np.asarray(
                node["qha_extrapolation_mask"], dtype=np.bool_
            ),
            elastic_extrapolation_mask=np.asarray(
                node["elastic_extrapolation_mask"], dtype=np.bool_
            ),
            stiffness_adiabatic=_read_optional_array(node, "stiffness_adiabatic"),
            sigma_stiffness_adiabatic=_read_optional_array(
                node, "sigma_stiffness_adiabatic"
            ),
            adiabatic_correction=_read_optional_array(node, "adiabatic_correction"),
            adiabatic_thermal_stress=_read_optional_array(
                node, "adiabatic_thermal_stress"
            ),
            adiabatic_valid_mask=_read_optional_bool_array(
                node, "adiabatic_valid_mask"
            ),
            stability=_read_stability(node),
            metadata=_read_metadata(node),
        )
    return profiles


def _read_result_fields(
    group: h5py.Group,
    *,
    reference: ReferenceEOSFit,
    components: dict[str, ElasticComponentFit],
    profiles: dict[str, ThermoelasticProfileResult],
    labels: tuple[str, ...],
) -> ThermoelasticResult:
    """Read the common calibration or analysis tensor fields."""
    return ThermoelasticResult(
        jobname=decode_text(group.attrs.get("jobname", "Unknown")),
        reference_eos=reference,
        component_fits=components,
        independent_labels=labels,
        temperature=np.asarray(group["temperature"], dtype=np.float64),
        pressure=np.asarray(group["pressure"], dtype=np.float64),
        equilibrium_volume=np.asarray(group["equilibrium_volume"], dtype=np.float64),
        density=np.asarray(group["density"], dtype=np.float64),
        independent_stiffness=_read_optional_array(group, "independent_stiffness"),
        sigma_independent_stiffness=_read_optional_array(
            group, "sigma_independent_stiffness"
        ),
        independent_stiffness_covariance=_read_optional_array(
            group, "independent_stiffness_covariance"
        ),
        stiffness_isothermal=_read_optional_array(group, "stiffness_isothermal"),
        sigma_stiffness_isothermal=_read_optional_array(
            group, "sigma_stiffness_isothermal"
        ),
        extrapolation_mask=np.asarray(group["elastic_extrapolation"], dtype=np.bool_),
        sigma_equilibrium_volume=_read_optional_array(
            group, "sigma_equilibrium_volume"
        ),
        qha_extrapolation_mask=(
            np.asarray(group["qha_extrapolation"], dtype=np.bool_)
            if "qha_extrapolation" in group
            else None
        ),
        profiles=profiles,
        isochoric_heat_capacity_cell=_read_optional_array(
            group, "isochoric_heat_capacity_cell"
        ),
        sigma_isochoric_heat_capacity_cell=_read_optional_array(
            group, "sigma_isochoric_heat_capacity_cell"
        ),
        thermal_expansion_tensor=_read_optional_array(
            group, "thermal_expansion_tensor"
        ),
        sigma_thermal_expansion_tensor=_read_optional_array(
            group, "sigma_thermal_expansion_tensor"
        ),
        stiffness_adiabatic=_read_optional_array(group, "stiffness_adiabatic"),
        sigma_stiffness_adiabatic=_read_optional_array(
            group, "sigma_stiffness_adiabatic"
        ),
        adiabatic_correction=_read_optional_array(group, "adiabatic_correction"),
        adiabatic_thermal_stress=_read_optional_array(
            group, "adiabatic_thermal_stress"
        ),
        adiabatic_valid_mask=_read_optional_bool_array(group, "adiabatic_valid_mask"),
        stability=_read_stability(group),
        completed=bool(group.attrs.get("completed", True)),
        metadata=_read_metadata(group),
    )


def _write_optional_array(
    group: h5py.Group,
    name: str,
    value: np.ndarray | None,
) -> None:
    """Write one optional numerical array with an explicit absence marker."""
    if value is None:
        group.attrs[f"{name}_is_none"] = True
    else:
        write_array_dataset(group, name, value, compression=True)


def _read_optional_array(group: h5py.Group, name: str) -> np.ndarray | None:
    """Read one optional numerical array."""
    if bool(group.attrs.get(f"{name}_is_none", False)) or name not in group:
        return None
    return np.asarray(group[name], dtype=np.float64)


def _read_optional_bool_array(
    group: h5py.Group,
    name: str,
) -> np.ndarray | None:
    """Read one optional Boolean array without a float round trip."""
    if bool(group.attrs.get(f"{name}_is_none", False)) or name not in group:
        return None
    return np.asarray(group[name], dtype=np.bool_)


def _read_metadata(group: h5py.Group) -> dict[str, Any]:
    """Read a conventional child metadata mapping when present."""
    if "metadata" not in group:
        return {}
    node = group["metadata"]
    assert isinstance(node, h5py.Group)
    return read_group_mapping(node)


def _write_stability(
    group: h5py.Group,
    stability: StabilityFieldResult | None,
) -> None:
    """Write one optional mechanical-stability diagnostic field."""
    if stability is None:
        group.attrs["mechanical_stability_is_none"] = True
        return
    node = group.create_group("mechanical_stability")
    node.attrs["criterion"] = stability.criterion
    node.attrs["tolerance_GPa"] = float(stability.tolerance)
    write_array_dataset(node, "eigenvalues", stability.eigenvalues, compression=True)
    write_array_dataset(
        node,
        "minimum_eigenvalue",
        stability.minimum_eigenvalue,
        compression=True,
    )
    write_array_dataset(node, "stable_mask", stability.stable_mask, compression=True)
    write_array_dataset(
        node,
        "indeterminate_mask",
        stability.indeterminate_mask,
        compression=True,
    )
    metadata = node.create_group("metadata")
    write_mapping(metadata, stability.metadata)


def _read_stability(group: h5py.Group) -> StabilityFieldResult | None:
    """Read optional mechanical-stability diagnostics from an archive."""
    if bool(group.attrs.get("mechanical_stability_is_none", False)):
        return None
    if "mechanical_stability" not in group:
        return None
    node = group["mechanical_stability"]
    assert isinstance(node, h5py.Group)
    return StabilityFieldResult(
        eigenvalues=np.asarray(node["eigenvalues"], dtype=np.float64),
        minimum_eigenvalue=np.asarray(node["minimum_eigenvalue"], dtype=np.float64),
        stable_mask=np.asarray(node["stable_mask"], dtype=np.bool_),
        indeterminate_mask=np.asarray(node["indeterminate_mask"], dtype=np.bool_),
        tolerance=float(node.attrs.get("tolerance_GPa", 0.0)),
        criterion=decode_text(
            node.attrs.get(
                "criterion",
                "positive_definite_wallace_stiffness",
            )
        ),
        metadata=_read_metadata(node),
    )


def _read_fit_quality(group: h5py.Group) -> ThermoelasticFitQuality | None:
    """Read optional scientific support diagnostics for one component fit."""
    if "quality" not in group:
        return None
    node = group["quality"]
    assert isinstance(node, h5py.Group)
    data = read_group_mapping(node)
    issues_value = data.get("issues", [])
    issues: tuple[str, ...]
    if isinstance(issues_value, str):
        issues = (issues_value,)
    else:
        issues = tuple(str(value) for value in issues_value)
    return ThermoelasticFitQuality(
        level=cast(
            ThermoelasticQualityLevel,
            str(data.get("level", "unsupported")),
        ),
        issues=issues,
        n_observations=int(data.get("n_observations", 0)),
        n_parameters=int(data.get("n_parameters", 0)),
        degrees_of_freedom=int(data.get("degrees_of_freedom", 0)),
        eulerian_strain_min=float(data.get("eulerian_strain_min", 0.0)),
        eulerian_strain_max=float(data.get("eulerian_strain_max", 0.0)),
        eulerian_strain_span=float(data.get("eulerian_strain_span", 0.0)),
        reference_volume_bracketed=bool(data.get("reference_volume_bracketed", False)),
        reference_volume_distance_fraction=float(
            data.get("reference_volume_distance_fraction", 0.0)
        ),
        design_rank=int(data.get("design_rank", 0)),
        design_condition_number=float(
            data.get("design_condition_number", float("inf"))
        ),
        leverage=np.asarray(data.get("leverage", []), dtype=np.float64),
        maximum_leverage=float(data.get("maximum_leverage", 0.0)),
        maximum_relative_symmetry_spread=float(
            data.get("maximum_relative_symmetry_spread", 0.0)
        ),
        maximum_leave_one_out_parameter_change=_optional_float(
            data.get("maximum_leave_one_out_parameter_change")
        ),
        maximum_order_parameter_change=_optional_float(
            data.get("maximum_order_parameter_change")
        ),
        thresholds=dict(data.get("thresholds", {})),
    )


def _optional_float(value: Any) -> float | None:
    """Normalize an optional numerical mapping value."""
    if value is None:
        return None
    return float(value)


def _validated_hdf5_name(name: str) -> str:
    """Return one safe direct HDF5 child name for a user profile."""
    value = str(name).strip()
    if not value or "/" in value or value in {".", ".."}:
        raise ValueError(
            "thermoelastic profile names must be non-empty direct HDF5 names "
            "without '/'"
        )
    return value


__all__ = [
    "THERMOELASTIC_HDF5_SCHEMA",
    "read_thermoelastic_payload",
    "write_thermoelastic_payload",
]
