# -*- coding: utf-8 -*-

"""Session-aware public plot discovery for persistent EOS archives.

The contracts in this module keep archive, dataset, slot, and immutable-record
selection separate from the common :class:`~quantas.models.PlotInventory` used
to describe one selected successful fit.  They intentionally contain only
lightweight metadata and never retain HDF5 handles or numerical fit arrays.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from quantas.models import (
    PlotContextDescriptor,
    PlotInventory,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
)
from quantas.modules.eos.archive import EOSArchive, infer_result_slots
from quantas.modules.eos.diagnostics import EOSDiagnostics
from quantas.modules.eos.history import EOSFitRecord, EOSResultSlot, EOSSlotStatus
from quantas.modules.eos.inspection import EOSRecordDisposition
from quantas.modules.eos.models import EOSDataset, EOSFitDomain

from .spec import EOSPlotter


@dataclass(frozen=True, slots=True)
class EOSDatasetPlotDescriptor:
    """Describe one embedded EOS dataset without exposing numerical arrays.

    Parameters
    ----------
    dataset_id : int
        Archive-local immutable dataset identifier.
    jobname : str
        Human-readable dataset title.
    npoints, selected_npoints, excluded_npoints : int
        Total and input-selection observation counts.
    columns : tuple of str
        Canonical measured and uncertainty columns.
    units : tuple of tuple
        ``(column, unit)`` pairs for available normalized units.
    variable_coordinates, constant_coordinates : tuple of str
        Scientific coordinate classification over selected observations.
    slot_keys : tuple of str
        Result slots that the dataset can scientifically support.
    group_ids : tuple of int
        Optional positive input-group identifiers.
    """

    dataset_id: int
    jobname: str
    npoints: int
    selected_npoints: int
    excluded_npoints: int
    columns: tuple[str, ...]
    units: tuple[tuple[str, str], ...]
    variable_coordinates: tuple[str, ...]
    constant_coordinates: tuple[str, ...]
    slot_keys: tuple[str, ...]
    group_ids: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        """Validate lightweight dataset metadata."""
        if int(self.dataset_id) <= 0:
            raise ValueError("EOS dataset descriptor identifier must be positive")
        if int(self.npoints) <= 0:
            raise ValueError("EOS dataset descriptor must contain observations")
        if int(self.selected_npoints) <= 0:
            raise ValueError("EOS dataset descriptor must select observations")
        if int(self.selected_npoints) + int(self.excluded_npoints) != int(
            self.npoints
        ):
            raise ValueError("EOS dataset descriptor observation counts are inconsistent")
        if not str(self.jobname).strip():
            raise ValueError("EOS dataset descriptor jobname must be non-empty")
        if len(set(self.columns)) != len(self.columns):
            raise ValueError("EOS dataset descriptor columns must be unique")
        if len(set(self.slot_keys)) != len(self.slot_keys):
            raise ValueError("EOS dataset descriptor slots must be unique")

    def unit_for(self, column: str) -> str | None:
        """Return the normalized unit associated with one column."""
        return dict(self.units).get(column)


@dataclass(frozen=True, slots=True)
class EOSRecordPlotDescriptor:
    """Describe plotting availability for one immutable EOS fit record."""

    record_id: int
    dataset_id: int
    slot_key: str
    domain: EOSFitDomain
    target: str
    model_tag: str
    fit_status: str
    successful: bool
    disposition: EOSRecordDisposition
    current_accepted: bool
    representation_keys: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Normalize identifiers, enums, and stable representation keys."""
        if int(self.record_id) <= 0 or int(self.dataset_id) <= 0:
            raise ValueError("EOS record descriptor identifiers must be positive")
        slot = EOSResultSlot.parse(self.slot_key)
        domain = EOSFitDomain(self.domain)
        if slot.domain is not domain or slot.target != str(self.target).lower():
            raise ValueError("EOS record descriptor slot, domain, and target disagree")
        if not str(self.model_tag).strip() or not str(self.fit_status).strip():
            raise ValueError("EOS record descriptor model and status must be non-empty")
        if len(set(self.representation_keys)) != len(self.representation_keys):
            raise ValueError("EOS record representation keys must be unique")
        if not self.successful and self.representation_keys:
            raise ValueError("failed EOS records cannot advertise plot representations")
        object.__setattr__(self, "record_id", int(self.record_id))
        object.__setattr__(self, "dataset_id", int(self.dataset_id))
        object.__setattr__(self, "domain", domain)
        object.__setattr__(self, "target", slot.target)
        object.__setattr__(self, "slot_key", slot.key)
        object.__setattr__(self, "disposition", EOSRecordDisposition(self.disposition))

    @property
    def plottable(self) -> bool:
        """Return whether at least one plot representation is buildable."""
        return bool(self.successful and self.representation_keys)


@dataclass(frozen=True, slots=True)
class EOSSlotPlotDescriptor:
    """Describe one persistent EOS result slot and its record history."""

    key: str
    domain: EOSFitDomain
    target: str
    status: EOSSlotStatus
    accepted_record_id: int | None
    last_record_id: int | None
    record_ids: tuple[int, ...] = ()
    plottable_record_ids: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        """Validate slot identity and record references."""
        slot = EOSResultSlot.parse(self.key)
        domain = EOSFitDomain(self.domain)
        if slot.domain is not domain or slot.target != str(self.target).lower():
            raise ValueError("EOS slot descriptor key, domain, and target disagree")
        if len(set(self.record_ids)) != len(self.record_ids):
            raise ValueError("EOS slot record identifiers must be unique")
        if not set(self.plottable_record_ids).issubset(self.record_ids):
            raise ValueError("EOS plottable records must belong to the slot history")
        for record_identifier in (
            *self.record_ids,
            *self.plottable_record_ids,
        ):
            if int(record_identifier) <= 0:
                raise ValueError("EOS slot record identifiers must be positive")
        for state_identifier in (self.accepted_record_id, self.last_record_id):
            if (
                state_identifier is not None
                and int(state_identifier) not in self.record_ids
            ):
                raise ValueError("EOS slot state must reference its record history")
        object.__setattr__(self, "key", slot.key)
        object.__setattr__(self, "domain", domain)
        object.__setattr__(self, "target", slot.target)
        object.__setattr__(self, "status", EOSSlotStatus(self.status))


@dataclass(frozen=True, slots=True)
class EOSArchivePlotInventory:
    """Lightweight session-aware plot inventory for one EOS archive.

    ``selected_plots`` is populated only for the explicit record, accepted slot,
    or unique accepted result selected by :func:`describe_eos_plots`.  This
    preserves lazy archive browsing while exposing all dataset, slot, record,
    and representation choices needed by a CLI, GUI, notebook, or script.
    """

    path: Path
    schema_version: str
    datasets: tuple[EOSDatasetPlotDescriptor, ...]
    slots: tuple[EOSSlotPlotDescriptor, ...]
    records: tuple[EOSRecordPlotDescriptor, ...]
    event_count: int
    selected_record_id: int | None = None
    selected_plots: PlotInventory | None = None
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate archive-wide uniqueness and selected-record consistency."""
        dataset_ids = tuple(item.dataset_id for item in self.datasets)
        record_ids = tuple(item.record_id for item in self.records)
        slot_keys = tuple(item.key for item in self.slots)
        if len(set(dataset_ids)) != len(dataset_ids):
            raise ValueError("EOS archive inventory contains duplicate datasets")
        if len(set(record_ids)) != len(record_ids):
            raise ValueError("EOS archive inventory contains duplicate records")
        if len(set(slot_keys)) != len(slot_keys):
            raise ValueError("EOS archive inventory contains duplicate slots")
        if int(self.event_count) < 0:
            raise ValueError("EOS archive inventory event count cannot be negative")
        if self.selected_record_id is None and self.selected_plots is not None:
            raise ValueError("selected EOS plots require a selected record")
        if self.selected_record_id is not None:
            if int(self.selected_record_id) not in record_ids:
                raise ValueError("selected EOS record is absent from archive inventory")
            if self.selected_plots is None:
                raise ValueError("selected EOS record requires a plot inventory")
            if self.selected_plots.module != "eos":
                raise ValueError("selected EOS plot inventory has the wrong module")
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(self, "event_count", int(self.event_count))

    def dataset_by_id(self, dataset_id: int) -> EOSDatasetPlotDescriptor:
        """Return one embedded dataset descriptor."""
        for item in self.datasets:
            if item.dataset_id == int(dataset_id):
                return item
        raise KeyError(f"unknown EOS dataset {dataset_id}")

    def slot_by_key(self, key: str) -> EOSSlotPlotDescriptor:
        """Return one result-slot descriptor."""
        canonical = EOSResultSlot.parse(key).key
        for item in self.slots:
            if item.key == canonical:
                return item
        raise KeyError(f"unknown EOS result slot {canonical!r}")

    def record_by_id(self, record_id: int) -> EOSRecordPlotDescriptor:
        """Return one immutable record descriptor."""
        for item in self.records:
            if item.record_id == int(record_id):
                return item
        raise KeyError(f"unknown EOS fit record {record_id}")


_REPRESENTATION_DESCRIPTIONS = {
    "fit": "Observed data and the fitted P-V or V-T relation.",
    "residuals": "Physical fit residuals against the natural control coordinate.",
    "standardized_residuals": (
        "Dimensionless residuals normalized by the effective observation uncertainty."
    ),
    "normalized_pressure": (
        "Finite-strain normalized-pressure diagnostic for compatible P-V models."
    ),
    "coverage": "Pressure-temperature sampling coverage of a P-V-T dataset.",
    "isotherms": "Calculated pressure-volume curves at selected temperatures.",
    "isobars": "Calculated volume-temperature curves at selected pressures.",
}


def describe_eos_plots(
    archive_path: str | Path,
    *,
    slot: str | EOSResultSlot | None = None,
    record_id: int | None = None,
) -> EOSArchivePlotInventory:
    """Describe archive, history, and plots available from one EOS archive.

    Parameters
    ----------
    archive_path : str or Path
        Native EOS HDF5 archive.
    slot : str, EOSResultSlot, or None, optional
        Accepted ``domain/target`` slot to select for detailed plot discovery.
    record_id : int or None, optional
        Explicit immutable record to select. When both selectors are omitted,
        a unique accepted result is selected automatically; zero or multiple
        accepted results leave the archive inventory unselected and add a
        non-fatal warning.

    Returns
    -------
    EOSArchivePlotInventory
        Lightweight datasets, slots, record history, and a detailed common
        :class:`PlotInventory` for the selected successful record.
    """
    warnings: list[str] = []
    with EOSArchive(archive_path) as archive:
        inspection = archive.inspect()
        datasets = tuple(
            _dataset_descriptor(dataset_id, archive.dataset(dataset_id))
            for dataset_id in archive.dataset_ids
        )
        selected_record, selection_warning = _select_record(
            archive,
            slot=slot,
            record_id=record_id,
        )
        if selection_warning is not None:
            warnings.append(selection_warning)

        record_descriptors: list[EOSRecordPlotDescriptor] = []
        representation_by_record: dict[int, tuple[str, ...]] = {}
        for slot_inspection in inspection.slots:
            for record_inspection in slot_inspection.records:
                record = record_inspection.record
                representations: tuple[str, ...] = ()
                if record.successful:
                    try:
                        plotter = EOSPlotter(record, archive.dataset(record.dataset_id))
                        representations = tuple(
                            _inventory_key(item)
                            for item in plotter.available_plot_types()
                        )
                    except (ValueError, RuntimeError, NotImplementedError) as exc:
                        warnings.append(
                            f"EOS record #{record.record_id} plot discovery: {exc}"
                        )
                representation_by_record[record.record_id] = representations
                record_descriptors.append(
                    EOSRecordPlotDescriptor(
                        record_id=record.record_id,
                        dataset_id=record.dataset_id,
                        slot_key=record.slot.key,
                        domain=record.request.domain,
                        target=record.request.target,
                        model_tag=record_inspection.model_tag,
                        fit_status=record.result.fit.status.value,
                        successful=record.successful,
                        disposition=record_inspection.disposition,
                        current_accepted=record_inspection.is_current_accepted,
                        representation_keys=representations,
                    )
                )

        records_by_slot = {
            slot_inspection.state.slot.key: tuple(
                item.record_id
                for item in record_descriptors
                if item.slot_key == slot_inspection.state.slot.key
            )
            for slot_inspection in inspection.slots
        }
        slots = tuple(
            EOSSlotPlotDescriptor(
                key=slot_inspection.state.slot.key,
                domain=slot_inspection.state.slot.domain,
                target=slot_inspection.state.slot.target,
                status=slot_inspection.state.status,
                accepted_record_id=slot_inspection.state.accepted_record_id,
                last_record_id=slot_inspection.state.last_record_id,
                record_ids=records_by_slot[slot_inspection.state.slot.key],
                plottable_record_ids=tuple(
                    identifier
                    for identifier in records_by_slot[slot_inspection.state.slot.key]
                    if representation_by_record.get(identifier)
                ),
            )
            for slot_inspection in inspection.slots
        )

        selected_plots = None
        selected_record_id = None
        if selected_record is not None:
            selected_record_id = selected_record.record_id
            selected_plots = _record_plot_inventory(
                selected_record,
                archive.dataset(selected_record.dataset_id),
                disposition=next(
                    item.disposition
                    for item in record_descriptors
                    if item.record_id == selected_record.record_id
                ),
            )

        return EOSArchivePlotInventory(
            path=Path(archive_path),
            schema_version=inspection.schema_version,
            datasets=datasets,
            slots=slots,
            records=tuple(record_descriptors),
            event_count=inspection.event_count,
            selected_record_id=selected_record_id,
            selected_plots=selected_plots,
            warnings=tuple(warnings),
        )


def _dataset_descriptor(
    dataset_id: int,
    dataset: EOSDataset,
) -> EOSDatasetPlotDescriptor:
    """Build one lightweight dataset descriptor."""
    classification = dataset.classify()
    return EOSDatasetPlotDescriptor(
        dataset_id=dataset_id,
        jobname=dataset.jobname,
        npoints=dataset.npoints,
        selected_npoints=dataset.selected_npoints,
        excluded_npoints=dataset.excluded_npoints,
        columns=dataset.column_names,
        units=tuple(
            (name, dataset.units[name])
            for name in dataset.column_names
            if name in dataset.units
        ),
        variable_coordinates=classification.variable_coordinates,
        constant_coordinates=classification.constant_coordinates,
        slot_keys=tuple(slot.key for slot in infer_result_slots(dataset)),
        group_ids=dataset.group_ids,
    )


def _select_record(
    archive: EOSArchive,
    *,
    slot: str | EOSResultSlot | None,
    record_id: int | None,
) -> tuple[EOSFitRecord | None, str | None]:
    """Resolve an explicit, accepted-slot, or unique accepted record."""
    if record_id is not None:
        record = archive.record(int(record_id))
        if slot is not None and record.slot != EOSResultSlot.parse(slot):
            raise ValueError("record_id does not belong to the requested EOS slot")
        if not record.successful:
            raise ValueError("EOS plot discovery requires a successful fit record")
        return record, None
    if slot is not None:
        resolved = EOSResultSlot.parse(slot)
        accepted_record = archive.accepted(resolved)
        if accepted_record is None:
            return None, f"EOS result slot {resolved.key!r} has no accepted record"
        return accepted_record, None
    accepted = tuple(
        archive.record(state.accepted_record_id)
        for state in archive.slots()
        if state.accepted_record_id is not None
    )
    if len(accepted) == 1:
        return accepted[0], None
    if not accepted:
        return None, "EOS archive has no accepted record selected for plotting"
    keys = ", ".join(record.slot.key for record in accepted)
    return (
        None,
        "EOS archive has multiple accepted records; select a slot or record_id "
        f"for detailed plot discovery: {keys}",
    )


def _record_plot_inventory(
    record: EOSFitRecord,
    dataset: EOSDataset,
    *,
    disposition: EOSRecordDisposition,
) -> PlotInventory:
    """Build the common detailed plot inventory for one successful record."""
    plotter = EOSPlotter(record, dataset)
    available = tuple(_inventory_key(item) for item in plotter.available_plot_types())
    diagnostic = EOSDiagnostics(record, dataset).build(include_normalized_pressure=True)
    representations = tuple(
        _representation_descriptor(key, record)
        for key in available
    )
    property_keys = {
        key
        for representation in representations
        for key in representation.property_keys
    }
    properties = tuple(
        _property_descriptor(
            key,
            record,
            dataset,
            diagnostic.units,
            diagnostic.metadata,
            available,
        )
        for key in _ordered_property_keys(property_keys)
    )
    contexts = [
        PlotContextDescriptor(
            key="record_id",
            name="Immutable fit record",
            description="Archive-global immutable fit-record identifier.",
            values=(record.record_id,),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="dataset_id",
            name="Embedded dataset",
            description="Dataset used by the selected immutable fit record.",
            values=(record.dataset_id,),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="result_slot",
            name="Result slot",
            description="Persistent scientific domain/target result slot.",
            values=(record.slot.key,),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="fit_domain",
            name="Fit domain",
            description="Natural scientific coordinate domain of the fit.",
            values=(record.request.domain.value,),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="fit_target",
            name="Fit target",
            description="Dependent coordinate addressed by the result slot.",
            values=(record.request.target,),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="model_tag",
            name="EOS model",
            description="Normalized EOS formulation used by the selected record.",
            values=(str(getattr(record.request.model, "tag", record.request.model)),),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="fit_status",
            name="Fit status",
            description="Numerical solver status stored in the immutable record.",
            values=(record.result.fit.status.value,),
            selectable=False,
        ),
        PlotContextDescriptor(
            key="record_disposition",
            name="Record disposition",
            description="Current scientific decision derived from archive history.",
            values=(disposition.value,),
            selectable=False,
        ),
    ]
    if "isotherms" in available:
        contexts.append(
            PlotContextDescriptor(
                key="isotherm_temperature",
                name="Isotherm temperature",
                description=(
                    "Finite temperature supplied to the P-V-T calculator. An empty "
                    "selection requests representative temperatures from the sampled "
                    "range; interpolation of stored observations is not performed."
                ),
                unit=dataset.units.get("temperature"),
            )
        )
    if "isobars" in available:
        contexts.append(
            PlotContextDescriptor(
                key="isobar_pressure",
                name="Isobar pressure",
                description=(
                    "Finite pressure supplied to the P-V-T calculator. An empty "
                    "selection requests representative pressures from the sampled "
                    "range; interpolation of stored observations is not performed."
                ),
                unit=dataset.units.get("pressure"),
            )
        )
    if dataset.group_ids:
        contexts.append(
            PlotContextDescriptor(
                key="data_group",
                name="Input data group",
                description="Positive group identifiers stored with the observations.",
                values=dataset.group_ids,
                selectable=False,
            )
        )
    return PlotInventory(
        module="eos",
        properties=properties,
        representations=representations,
        contexts=tuple(contexts),
        warnings=tuple(str(item) for item in diagnostic.warnings),
    )


def _representation_descriptor(
    key: str,
    record: EOSFitRecord,
) -> PlotRepresentationDescriptor:
    """Describe one canonical EOS plot representation."""
    property_keys: tuple[str, ...]
    supported_contexts: tuple[str, ...] = (
        "record_id",
        "dataset_id",
        "result_slot",
        "fit_domain",
        "fit_target",
        "model_tag",
        "fit_status",
        "record_disposition",
    )
    constraints: tuple[str, ...] = ()
    if key == "fit":
        property_keys = (
            "pressure"
            if record.request.domain is EOSFitDomain.PRESSURE_VOLUME
            else _target_property_key(record.request.target),
        )
    elif key == "residuals":
        property_keys = ("residual",)
    elif key == "standardized_residuals":
        property_keys = ("standardized_residual",)
    elif key == "normalized_pressure":
        property_keys = ("normalized_pressure",)
        constraints = (
            "Available only for compatible finite-strain P-V formulations.",
        )
    elif key == "coverage":
        property_keys = ("sampling_coverage",)
    elif key == "isotherms":
        property_keys = ("pressure",)
        supported_contexts = (*supported_contexts, "isotherm_temperature")
        constraints = (
            "Curves are calculated from the fitted P-V-T model, not interpolated from observations.",
        )
    elif key == "isobars":
        property_keys = ("volume",)
        supported_contexts = (*supported_contexts, "isobar_pressure")
        constraints = (
            "Curves are calculated from the fitted P-V-T model, not interpolated from observations.",
        )
    else:  # Defensive guard against an unregistered builder kind.
        raise ValueError(f"unknown EOS inventory representation {key!r}")
    return PlotRepresentationDescriptor(
        key=key,
        name=key.replace("_", " ").title(),
        plot_kind="line",
        description=_REPRESENTATION_DESCRIPTIONS[key],
        property_keys=property_keys,
        supported_contexts=supported_contexts,
        constraints=constraints,
    )


def _property_descriptor(
    key: str,
    record: EOSFitRecord,
    dataset: EOSDataset,
    diagnostic_units: dict[str, str],
    diagnostic_metadata: dict[str, object],
    representations: tuple[str, ...],
) -> PlotPropertyDescriptor:
    """Describe one scalar or sampling quantity used by EOS plots."""
    compatible = tuple(
        representation
        for representation in representations
        if key in _representation_property_keys(representation, record)
    )
    if key == "pressure":
        return PlotPropertyDescriptor(
            key=key,
            name="Pressure",
            symbol_math="P",
            symbol_plain="P",
            unit=dataset.units.get("pressure"),
            description="Observed or model-calculated pressure.",
            category="state_variable",
            representations=compatible,
        )
    if key == "volume":
        return PlotPropertyDescriptor(
            key=key,
            name="Volume",
            symbol_math="V",
            symbol_plain="V",
            unit=dataset.units.get("volume"),
            description="Observed or model-calculated volume.",
            category="state_variable",
            representations=compatible,
        )
    if key.startswith("axis_"):
        axis = key.removeprefix("axis_")
        return PlotPropertyDescriptor(
            key=key,
            name=f"Cell parameter {axis}",
            symbol_math=axis,
            symbol_plain=axis,
            unit=dataset.units.get(axis),
            description=f"Observed or model-calculated lattice parameter {axis}.",
            category="structural_coordinate",
            representations=compatible,
        )
    if key == "residual":
        return PlotPropertyDescriptor(
            key=key,
            name="Fit residual",
            symbol_math=r"\Delta y",
            symbol_plain="Δy",
            unit=diagnostic_units.get("residual"),
            description="Observed minus calculated response in physical units.",
            category="diagnostic",
            representations=compatible,
        )
    if key == "standardized_residual":
        return PlotPropertyDescriptor(
            key=key,
            name="Standardized residual",
            symbol_math="r_{std}",
            symbol_plain="r_std",
            unit=None,
            description="Residual normalized by the effective observation uncertainty.",
            category="diagnostic",
            representations=compatible,
        )
    if key == "normalized_pressure":
        normalized = diagnostic_metadata.get("normalized_pressure", {})
        metadata = normalized if isinstance(normalized, dict) else {}
        symbol = str(metadata.get("normalized_pressure_symbol", "F"))
        return PlotPropertyDescriptor(
            key=key,
            name="Normalized pressure",
            symbol_math=symbol,
            symbol_plain=symbol,
            unit=diagnostic_units.get("normalized_pressure"),
            description="Finite-strain normalized pressure used for EOS-order diagnosis.",
            category="diagnostic",
            representations=compatible,
        )
    if key == "sampling_coverage":
        return PlotPropertyDescriptor(
            key=key,
            name="Pressure-temperature sampling coverage",
            symbol_math="(P,T)",
            symbol_plain="P–T",
            unit=None,
            description="Locations of included and excluded observations in P-T space.",
            category="sampling",
            representations=compatible,
        )
    raise ValueError(f"unknown EOS inventory property {key!r}")


def _representation_property_keys(
    key: str,
    record: EOSFitRecord,
) -> tuple[str, ...]:
    """Return properties associated with one representation key."""
    if key == "fit":
        return (
            "pressure"
            if record.request.domain is EOSFitDomain.PRESSURE_VOLUME
            else _target_property_key(record.request.target),
        )
    return {
        "residuals": ("residual",),
        "standardized_residuals": ("standardized_residual",),
        "normalized_pressure": ("normalized_pressure",),
        "coverage": ("sampling_coverage",),
        "isotherms": ("pressure",),
        "isobars": ("volume",),
    }[key]


def _target_property_key(target: str) -> str:
    """Return a stable plot-property key for one EOS target."""
    normalized = str(target).strip().lower()
    if normalized == "volume":
        return "volume"
    if normalized in {"a", "b", "c"}:
        return f"axis_{normalized}"
    raise ValueError(f"unsupported EOS plot target {target!r}")


def _ordered_property_keys(keys: set[str]) -> tuple[str, ...]:
    """Return deterministic scientific property ordering."""
    order = (
        "pressure",
        "volume",
        "axis_a",
        "axis_b",
        "axis_c",
        "sampling_coverage",
        "residual",
        "standardized_residual",
        "normalized_pressure",
    )
    return tuple(key for key in order if key in keys)


def _inventory_key(plot_type: str) -> str:
    """Convert legacy hyphenated EOS plot identifiers to descriptor keys."""
    return str(plot_type).strip().lower().replace("-", "_")


__all__ = [
    "EOSArchivePlotInventory",
    "EOSDatasetPlotDescriptor",
    "EOSRecordPlotDescriptor",
    "EOSSlotPlotDescriptor",
    "describe_eos_plots",
]
