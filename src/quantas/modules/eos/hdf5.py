# -*- coding: utf-8 -*-

"""Native HDF5 schema and round-trip helpers for persistent EOS archives."""

from __future__ import annotations

from dataclasses import fields
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping, cast

import h5py
import numpy as np

from quantas import __version__
from quantas.core.math.fitting import (
    CovarianceScaling,
    EffectiveVarianceOptions,
    FitDiagnostics,
    FitMethod,
    FitOptions,
    FitQuality,
    FitResult,
    FitStatus,
    ODRDifferenceScheme,
    OLSOptions,
    OrthogonalDistanceOptions,
    ParameterState,
    SolverOptions,
    WLSOptions,
)
from quantas.core.physics.eos import MGDNormalization, PVTModel
from quantas.io.hdf5 import (
    decode_optional_text,
    decode_text,
    numeric_sort_key,
    parse_datetime,
    read_mapping,
    write_array_dataset,
    write_mapping,
)
from quantas.io.hdf5.envelope import write_precision_metadata

from .history import (
    EOSFitRecord,
    EOSResultSlot,
    EOSSlotState,
    EOSSlotStatus,
    EOSStateEvent,
    EOSStateEventType,
)
from .models import (
    EOSDataset,
    EOSFitDomain,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitResult,
    ParameterConstraint,
)

EOS_ARCHIVE_SCHEMA_VERSION = "1.1"
EOS_SUPPORTED_ARCHIVE_SCHEMA_VERSIONS = frozenset({"1.0", "1.1"})
EOS_ARCHIVE_KIND = "quantas-eos-archive"

_OPTIONS_TYPES: dict[str, type[FitOptions]] = {
    "OLSOptions": OLSOptions,
    "WLSOptions": WLSOptions,
    "EffectiveVarianceOptions": EffectiveVarianceOptions,
    "OrthogonalDistanceOptions": OrthogonalDistanceOptions,
    "ODROptions": OrthogonalDistanceOptions,
}


def initialize_eos_archive(
    h5: h5py.File,
    *,
    creator: str | None = None,
    created_at: datetime | None = None,
) -> None:
    """Initialize an empty native EOS archive.

    Parameters
    ----------
    h5 : h5py.File
        Newly created writable HDF5 file.
    creator : str or None, optional
        Human or workflow identifier creating the archive.
    created_at : datetime or None, optional
        Explicit creation time. UTC now is used by default.

    Raises
    ------
    ValueError
        If the destination is not empty.
    """
    if len(h5) != 0:
        raise ValueError("EOS archive initialization requires an empty HDF5 file")
    timestamp = created_at or datetime.now(timezone.utc)
    metadata = h5.create_group("metadata")
    metadata.attrs["archive_kind"] = EOS_ARCHIVE_KIND
    metadata.attrs["program"] = "quantas"
    metadata.attrs["version"] = __version__
    metadata.attrs["module"] = "eos"
    metadata.attrs["method"] = "persistent EOS archive"
    metadata.attrs["schema_version"] = EOS_ARCHIVE_SCHEMA_VERSION
    metadata.attrs["created_by"] = creator or "unknown"
    metadata.attrs["created_at"] = timestamp.isoformat()
    write_precision_metadata(h5)

    input_group = h5.create_group("input")
    input_group.create_group("datasets")
    session = h5.create_group("session")
    session.attrs["next_dataset_id"] = 1
    session.attrs["next_record_id"] = 1
    session.attrs["next_event_id"] = 1
    session.create_group("records")
    session.create_group("state_events")
    session.create_group("current")
    results = h5.create_group("results")
    results.create_group("accepted")
    events = h5.create_group("events")
    events.attrs["count"] = 0
    h5.flush()


def validate_eos_archive(h5: h5py.File) -> None:
    """Validate the identity and supported schema of an EOS archive.

    Parameters
    ----------
    h5 : h5py.File
        Open source file.

    Raises
    ------
    ValueError
        If required metadata or groups are missing or unsupported.
    """
    required = {
        "metadata",
        "input/datasets",
        "session/records",
        "session/state_events",
        "session/current",
        "results/accepted",
    }
    missing = sorted(path for path in required if path not in h5)
    if missing:
        raise ValueError(f"EOS archive is missing required nodes: {', '.join(missing)}")
    metadata = h5["metadata"].attrs
    kind = decode_text(metadata.get("archive_kind", ""))
    if kind != EOS_ARCHIVE_KIND:
        raise ValueError(f"not a native EOS archive: archive_kind={kind!r}")
    version = decode_text(metadata.get("schema_version", ""))
    if version not in EOS_SUPPORTED_ARCHIVE_SCHEMA_VERSIONS:
        supported = ", ".join(sorted(EOS_SUPPORTED_ARCHIVE_SCHEMA_VERSIONS))
        raise ValueError(
            f"unsupported EOS archive schema {version!r}; "
            f"supported versions are {supported}"
        )


def write_eos_dataset(group: h5py.Group, dataset_id: int, dataset: EOSDataset) -> None:
    """Write one complete input dataset into a newly created group.

    Parameters
    ----------
    group : h5py.Group
        Destination dataset group.
    dataset_id : int
        Positive archive-local identifier.
    dataset : EOSDataset
        Raw and normalized EOS data.
    """
    group.attrs["dataset_id"] = int(dataset_id)
    group.attrs["jobname"] = dataset.jobname
    group.attrs["npoints"] = dataset.npoints
    group.attrs["schema_version"] = "1.1"
    if dataset.source is not None:
        group.attrs["source"] = str(dataset.source)

    source_text = _source_text(dataset.source)
    if source_text is not None:
        text = group.create_dataset("source_text", data=source_text)
        text.attrs["encoding"] = "utf-8"
        group.attrs["source_sha256"] = sha256(source_text.encode("utf-8")).hexdigest()

    raw_group = group.create_group("raw_data")
    raw_columns = dataset.raw_columns or dataset.columns
    raw_units = dataset.raw_units or dataset.units
    for name in dataset.column_names:
        values = raw_columns[name]
        write_array_dataset(
            raw_group,
            name,
            values,
            unit=raw_units.get(name),
            compression=True,
        )

    normalized_group = group.create_group("normalized_data")
    for name in dataset.column_names:
        write_array_dataset(
            normalized_group,
            name,
            dataset.columns[name],
            unit=dataset.units.get(name),
            compression=True,
        )

    uncertainty_group = group.create_group("uncertainties")
    for name in dataset.column_names:
        if name.startswith("sigma_"):
            write_array_dataset(
                uncertainty_group,
                name,
                dataset.columns[name],
                unit=dataset.units.get(name),
                compression=True,
            )

    selection_group = group.create_group("selection")
    write_array_dataset(
        selection_group,
        "default_mask",
        dataset.selection_mask().astype(np.bool_),
        compression=True,
    )
    if dataset.groups is not None:
        write_array_dataset(
            selection_group,
            "groups",
            np.asarray(dataset.groups, dtype=np.int64),
            compression=True,
        )

    write_mapping(group.create_group("units"), dataset.units)
    write_mapping(group.create_group("raw_units"), raw_units)
    write_mapping(group.create_group("provenance"), dataset.provenance)
    write_mapping(group.create_group("metadata"), dataset.metadata)
    write_mapping(group.create_group("classification"), dataset.classify().as_dict())
    write_array_dataset(
        group,
        "column_order",
        np.asarray(dataset.column_names, dtype=object),
    )


def read_eos_dataset(group: h5py.Group) -> EOSDataset:
    """Read one complete input dataset from an EOS archive group."""
    if "normalized_data" not in group:
        raise ValueError(f"missing normalized_data in {group.name}")
    columns = {
        name: np.asarray(node[()], dtype=np.float64)
        for name, node in group["normalized_data"].items()
    }
    raw_columns = {
        name: np.asarray(node[()], dtype=np.float64)
        for name, node in group.get("raw_data", {}).items()
    }
    units = read_mapping(group["units"]) if "units" in group else {}
    raw_units = read_mapping(group["raw_units"]) if "raw_units" in group else {}
    provenance = read_mapping(group["provenance"]) if "provenance" in group else {}
    metadata = (
        _restore_dataset_metadata(read_mapping(group["metadata"]))
        if "metadata" in group
        else {}
    )
    source = decode_optional_text(group.attrs.get("source"))
    selection = group.get("selection")
    default_mask = None
    groups = None
    if selection is not None:
        if "default_mask" in selection:
            default_mask = np.asarray(selection["default_mask"][()], dtype=np.bool_)
        if "groups" in selection:
            groups = np.asarray(selection["groups"][()], dtype=np.int64)
    return EOSDataset(
        jobname=decode_text(group.attrs.get("jobname", "Unknown")),
        columns=columns,
        units={str(key): str(value) for key, value in units.items()},
        raw_columns=raw_columns,
        raw_units={str(key): str(value) for key, value in raw_units.items()},
        provenance={str(key): str(value) for key, value in provenance.items()},
        source=source,
        metadata=dict(metadata),
        default_mask=default_mask,
        groups=groups,
    )


def write_fit_record(group: h5py.Group, record: EOSFitRecord) -> None:
    """Write one immutable fit record into a newly created group."""
    group.attrs["record_id"] = record.record_id
    group.attrs["dataset_id"] = record.dataset_id
    group.attrs["parent_record_id"] = (
        -1 if record.parent_record_id is None else record.parent_record_id
    )
    group.attrs["created_at"] = record.created_at.isoformat()
    group.attrs["slot"] = record.slot.key
    group.attrs["success"] = record.successful
    group.attrs["fit_status"] = record.result.fit.status.value
    group.attrs["immutable"] = True
    if record.note is not None:
        group.attrs["note"] = record.note
    write_mapping(group.create_group("provenance"), record.provenance)
    write_mapping(group.create_group("request"), record.request.as_dict())
    write_mapping(group.create_group("result"), record.result.as_dict())


def read_fit_record(group: h5py.Group) -> EOSFitRecord:
    """Reconstruct one immutable fit record from HDF5."""
    request = eos_fit_request_from_mapping(read_mapping(group["request"]))
    result = eos_fit_result_from_mapping(read_mapping(group["result"]), request=request)
    parent = int(group.attrs.get("parent_record_id", -1))
    return EOSFitRecord(
        record_id=int(group.attrs["record_id"]),
        dataset_id=int(group.attrs["dataset_id"]),
        parent_record_id=None if parent < 0 else parent,
        created_at=parse_datetime(group.attrs.get("created_at")),
        request=request,
        result=result,
        note=decode_optional_text(group.attrs.get("note")),
        provenance=(read_mapping(group["provenance"]) if "provenance" in group else {}),
    )


def write_state_event(group: h5py.Group, event: EOSStateEvent) -> None:
    """Write one immutable append-only EOS state event."""
    group.attrs["event_id"] = event.event_id
    group.attrs["event_type"] = event.event_type.value
    group.attrs["record_id"] = -1 if event.record_id is None else event.record_id
    group.attrs["slot"] = "" if event.slot is None else event.slot.key
    group.attrs["created_at"] = event.created_at.isoformat()
    if event.note is not None:
        group.attrs["note"] = event.note
    group.attrs["immutable"] = True
    write_mapping(group.create_group("metadata"), event.metadata)


def read_state_event(group: h5py.Group) -> EOSStateEvent:
    """Read one append-only EOS state event."""
    record_id = int(group.attrs.get("record_id", -1))
    slot_text = decode_text(group.attrs.get("slot", ""))
    return EOSStateEvent(
        event_id=int(group.attrs["event_id"]),
        event_type=EOSStateEventType(decode_text(group.attrs["event_type"])),
        record_id=None if record_id < 0 else record_id,
        slot=None if not slot_text else EOSResultSlot.parse(slot_text),
        created_at=parse_datetime(group.attrs.get("created_at")),
        note=decode_optional_text(group.attrs.get("note")),
        metadata=read_mapping(group["metadata"]) if "metadata" in group else {},
    )


def write_slot_state(group: h5py.Group, state: EOSSlotState) -> None:
    """Replace the mutable compact state of one result slot."""
    group.attrs["domain"] = state.slot.domain.value
    group.attrs["target"] = state.slot.target
    group.attrs["status"] = state.status.value
    group.attrs["accepted_record_id"] = (
        -1 if state.accepted_record_id is None else state.accepted_record_id
    )
    group.attrs["last_record_id"] = (
        -1 if state.last_record_id is None else state.last_record_id
    )
    if "attempted_record_ids" in group:
        del group["attempted_record_ids"]
    write_array_dataset(
        group,
        "attempted_record_ids",
        np.asarray(state.attempted_record_ids, dtype=np.int64),
    )


def read_slot_state(group: h5py.Group) -> EOSSlotState:
    """Read the mutable compact state of one result slot."""
    accepted = int(group.attrs.get("accepted_record_id", -1))
    last = int(group.attrs.get("last_record_id", -1))
    attempts = (
        tuple(
            int(value)
            for value in np.asarray(group["attempted_record_ids"][()]).ravel()
        )
        if "attempted_record_ids" in group
        else ()
    )
    return EOSSlotState(
        slot=EOSResultSlot(
            EOSFitDomain(decode_text(group.attrs["domain"])),
            decode_text(group.attrs["target"]),
        ),
        status=EOSSlotStatus(decode_text(group.attrs["status"])),
        accepted_record_id=None if accepted < 0 else accepted,
        last_record_id=None if last < 0 else last,
        attempted_record_ids=attempts,
    )


def write_accepted_result(group: h5py.Group, record: EOSFitRecord) -> None:
    """Materialize a compact accepted-result copy with source linkage."""
    group.attrs["record_id"] = record.record_id
    group.attrs["dataset_id"] = record.dataset_id
    group.attrs["slot"] = record.slot.key
    group.attrs["source_record_path"] = f"/session/records/{record.record_id:06d}"
    group.attrs["materialized_at"] = datetime.now(timezone.utc).isoformat()
    write_mapping(group.create_group("request"), record.request.as_dict())
    write_mapping(group.create_group("result"), record.result.as_dict())


def eos_fit_request_from_mapping(values: Mapping[str, Any]) -> EOSFitRequest:
    """Reconstruct :class:`EOSFitRequest` from its serialized mapping."""
    domain = EOSFitDomain(str(values["domain"]))
    model_values = _mapping(values["model"])
    if domain is EOSFitDomain.VOLUME_TEMPERATURE:
        model: Any = str(model_values["tag"])
    elif domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
        pressure = _mapping(model_values["pressure_model"])
        thermal_value = model_values.get("temperature_model")
        thermal = None if thermal_value is None else str(_mapping(thermal_value)["tag"])
        thermal_pressure_value = model_values.get("thermal_pressure_model")
        thermal_pressure = (
            None
            if thermal_pressure_value is None
            else str(_mapping(thermal_pressure_value)["tag"])
        )
        normalization_value = model_values.get("mgd_normalization")
        normalization = None
        if normalization_value is not None:
            payload = _mapping(normalization_value)
            normalization = MGDNormalization(
                volume_basis=str(payload["volume_basis"]),
                atoms_per_unit=float(payload["atoms_per_unit"]),
                formula=_optional_str(payload.get("formula")),
                formula_units_per_cell=_optional_float(
                    payload.get("formula_units_per_cell")
                ),
            )
        model = PVTModel(
            pressure_model=str(pressure["tag"]),
            temperature_model=thermal,
            coupling=str(model_values["coupling"]),
            thermal_pressure_model=thermal_pressure,
            mgd_normalization=normalization,
        )
    else:
        model = str(model_values["tag"])

    constraints: list[ParameterConstraint] = []
    for item in _sequence(values.get("constraints", [])):
        payload = _mapping(item)
        constraints.append(
            ParameterConstraint(
                name=str(payload["name"]),
                state=ParameterState(str(payload["state"])),
                initial_value=_optional_float(payload.get("initial_value")),
                value=_optional_float(payload.get("value")),
                lower_bound=float(payload.get("lower_bound", -np.inf)),
                upper_bound=float(payload.get("upper_bound", np.inf)),
                unit=_optional_str(payload.get("unit")),
                description=str(payload.get("description", "")),
                metadata=dict(_mapping(payload.get("metadata", {}))),
            )
        )
    options_values = _mapping(values["options"])
    solver = fit_options_from_mapping(_mapping(options_values["solver_options"]))
    options = EOSFitOptions(
        solver_options=solver,
        allow_extrapolation=bool(options_values.get("allow_extrapolation", False)),
        metadata=dict(_mapping(options_values.get("metadata", {}))),
    )
    mask_values = values.get("mask")
    mask = None if mask_values is None else np.asarray(mask_values, dtype=np.bool_)
    return EOSFitRequest(
        model=model,
        target=str(values.get("target", "volume")),
        domain=domain,
        constraints=tuple(constraints),
        options=options,
        mask=mask,
        request_id=_optional_str(values.get("request_id")),
        metadata=dict(_mapping(values.get("metadata", {}))),
    )


def fit_options_from_mapping(values: Mapping[str, Any]) -> SolverOptions:
    """Reconstruct one concrete typed solver-options object."""
    type_name = str(values.get("type", ""))
    option_type = _OPTIONS_TYPES.get(type_name)
    if option_type is None:
        raise ValueError(f"unsupported serialized fitting-options type: {type_name!r}")
    allowed = {
        item.name
        for item in fields(option_type)
        if item.init and not item.name.startswith("_")
    }
    payload: dict[str, Any] = {}
    for key, value in values.items():
        if key in allowed and key not in {"absolute_sigma"}:
            payload[key] = value
    if "method" in payload:
        payload["method"] = FitMethod(str(payload["method"]))
    if "covariance_scaling" in payload:
        payload["covariance_scaling"] = CovarianceScaling(
            str(payload["covariance_scaling"])
        )
    if "difference_scheme" in payload:
        payload["difference_scheme"] = ODRDifferenceScheme(
            str(payload["difference_scheme"])
        )
    if "metadata" in payload:
        payload["metadata"] = dict(_mapping(payload["metadata"]))
    for name in (
        "initial_x_corrections",
        "parameter_steps",
        "x_steps",
        "parameter_scales",
        "x_scales",
    ):
        if name in payload and payload[name] is not None:
            payload[name] = np.asarray(payload[name], dtype=np.float64)
    return cast(SolverOptions, option_type(**payload))


def fit_result_from_mapping(values: Mapping[str, Any]) -> FitResult:
    """Reconstruct a complete general fitting result."""
    diagnostics_values = values.get("diagnostics")
    diagnostics = None
    if diagnostics_values is not None:
        payload = _mapping(diagnostics_values)
        diagnostics = FitDiagnostics(
            objective=str(payload.get("objective", "")),
            weighted=bool(payload.get("weighted", False)),
            chi_square=_optional_float(payload.get("chi_square")),
            reduced_chi_square=_optional_float(payload.get("reduced_chi_square")),
            correlation=_optional_array(payload.get("correlation")),
            standardized_residuals=_optional_array(
                payload.get("standardized_residuals")
            ),
            x_corrections=_optional_array(payload.get("x_corrections")),
            y_corrections=_optional_array(payload.get("y_corrections")),
            jacobian_rank=_optional_int(payload.get("jacobian_rank")),
            condition_number=_optional_float(payload.get("condition_number")),
            parameters_at_bounds=tuple(
                bool(item)
                for item in _sequence(payload.get("parameters_at_bounds", []))
            ),
            n_iterations=_optional_int(payload.get("n_iterations")),
            n_evaluations=_optional_int(payload.get("n_evaluations")),
            stop_reason=str(payload.get("stop_reason", "")),
            warnings=[str(item) for item in _sequence(payload.get("warnings", []))],
            metadata=dict(_mapping(payload.get("metadata", {}))),
        )
    return FitResult(
        success=bool(values["success"]),
        status=FitStatus(str(values["status"])),
        quality=FitQuality(str(values["quality"])),
        message=str(values.get("message", "")),
        parameters=_optional_array(values.get("parameters")),
        covariance=_optional_array(values.get("covariance")),
        errors=_optional_array(values.get("errors")),
        fitted=_optional_array(values.get("fitted")),
        residuals=_optional_array(values.get("residuals")),
        rmse=_optional_float(values.get("rmse")),
        mae=_optional_float(values.get("mae")),
        max_abs_error=_optional_float(values.get("max_abs_error")),
        r_squared=_optional_float(values.get("r_squared")),
        n_points=int(values.get("n_points", 0)),
        n_parameters=int(values.get("n_parameters", 0)),
        dof=int(values.get("dof", 0)),
        condition_number=_optional_float(values.get("condition_number")),
        warnings=[str(item) for item in _sequence(values.get("warnings", []))],
        metadata=dict(_mapping(values.get("metadata", {}))),
        method=(
            None if values.get("method") is None else FitMethod(str(values["method"]))
        ),
        parameter_names=tuple(
            str(item) for item in _sequence(values.get("parameter_names", []))
        ),
        parameter_states=tuple(
            ParameterState(str(item))
            for item in _sequence(values.get("parameter_states", []))
        ),
        optimizer_parameters=_optional_array(values.get("optimizer_parameters")),
        optimizer_covariance=_optional_array(values.get("optimizer_covariance")),
        optimizer_errors=_optional_array(values.get("optimizer_errors")),
        diagnostics=diagnostics,
    )


def eos_fit_result_from_mapping(
    values: Mapping[str, Any],
    *,
    request: EOSFitRequest | None = None,
) -> EOSFitResult:
    """Reconstruct a complete EOS-domain fit result."""
    resolved_request = request or eos_fit_request_from_mapping(
        _mapping(values["request"])
    )
    predictions = {
        str(name): np.asarray(data, dtype=np.float64)
        for name, data in _mapping(values.get("predictions", {})).items()
    }
    return EOSFitResult(
        request=resolved_request,
        fit=fit_result_from_mapping(_mapping(values["fit"])),
        predictions=predictions,
        derived={
            str(name): float(value)
            for name, value in _mapping(values.get("derived", {})).items()
        },
        warnings=[str(item) for item in _sequence(values.get("warnings", []))],
        metadata=dict(_mapping(values.get("metadata", {}))),
    )


def sorted_numeric_children(group: h5py.Group) -> tuple[str, ...]:
    """Return HDF5 child names in numeric-aware order."""
    return tuple(sorted(group.keys(), key=numeric_sort_key))


def _restore_dataset_metadata(value: Any) -> Any:
    """Restore reader metadata sequences to immutable Python tuples."""
    if isinstance(value, Mapping):
        return {
            str(key): _restore_dataset_metadata(item) for key, item in value.items()
        }
    if isinstance(value, np.ndarray):
        return tuple(_restore_dataset_metadata(item) for item in value.tolist())
    if isinstance(value, list):
        return tuple(_restore_dataset_metadata(item) for item in value)
    return value


def _source_text(source: str | Path | None) -> str | None:
    if source is None:
        return None
    path = Path(source)
    if not path.is_file():
        return None
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text(encoding="utf-8", errors="replace")


def _mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    raise TypeError(f"expected serialized mapping, received {type(value).__name__}")


def _sequence(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return list(value)
    if isinstance(value, Mapping):
        return [value[key] for key in sorted(value, key=numeric_sort_key)]
    raise TypeError(f"expected serialized sequence, received {type(value).__name__}")


def _optional_float(value: Any) -> float | None:
    return None if value is None else float(value)


def _optional_int(value: Any) -> int | None:
    return None if value is None else int(value)


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_array(value: Any) -> np.ndarray | None:
    return None if value is None else np.asarray(value, dtype=np.float64)


__all__ = [
    "EOS_ARCHIVE_KIND",
    "EOS_ARCHIVE_SCHEMA_VERSION",
    "EOS_SUPPORTED_ARCHIVE_SCHEMA_VERSIONS",
    "eos_fit_request_from_mapping",
    "eos_fit_result_from_mapping",
    "fit_options_from_mapping",
    "fit_result_from_mapping",
    "initialize_eos_archive",
    "read_eos_dataset",
    "read_fit_record",
    "read_slot_state",
    "read_state_event",
    "sorted_numeric_children",
    "validate_eos_archive",
    "write_accepted_result",
    "write_eos_dataset",
    "write_fit_record",
    "write_slot_state",
    "write_state_event",
]
