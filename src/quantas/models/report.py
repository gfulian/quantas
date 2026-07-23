# -*- coding: utf-8 -*-
"""Frontend-neutral report models shared by Quantas workflows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class ReportTable:
    """Describe one frontend-neutral table of scientific results.

    Parameters
    ----------
    title : str
        Human-readable table title.
    columns : list of str
        Ordered column labels.
    rows : list of list
        Ordered table rows. Cells retain raw numerical or textual values until
        a frontend renderer formats them.
    metadata : dict, optional
        Frontend-neutral formatting and provenance metadata.
    """

    title: str
    columns: list[str]
    rows: list[list[Any]]
    metadata: dict[str, Any] = field(default_factory=dict)


def mapping_table(
    title: str,
    values: Mapping[str, Any],
    *,
    key_column: str = "Property",
    value_column: str = "Value",
) -> ReportTable:
    """Build a compact two-column table from a mapping.

    Nested arrays and mappings are summarized rather than expanded into their
    full Python representation. This keeps reports readable while the complete
    values remain available in the result envelope and native HDF5 output.

    Parameters
    ----------
    title : str
        Table title.
    values : mapping
        Values to render in insertion order.
    key_column, value_column : str, optional
        Header labels.

    Returns
    -------
    ReportTable
        Neutral table containing one row per mapping item.
    """
    rows = [[_humanize_key(key), _report_value(value)] for key, value in values.items()]
    return ReportTable(title=title, columns=[key_column, value_column], rows=rows)


def input_data_table(
    title: str,
    values: Mapping[str, Any],
    *,
    source: str | Path | None = None,
) -> ReportTable:
    """Build a readable table for one normalized Quantas input envelope.

    Structural volume-series data are represented by a concise scientific
    summary. The complete structure remains stored in :class:`InputData` and
    in native result files.

    Parameters
    ----------
    title : str
        Table title.
    values : mapping
        Normalized input mapping.
    source : str, Path, or None, optional
        Original input source.

    Returns
    -------
    ReportTable
        Input summary suitable for terminal and persistent reports.
    """
    rows: list[list[Any]] = []
    for key, value in values.items():
        if key == "structure":
            continue
        rows.append([_humanize_key(key), _report_value(value)])

    structure = values.get("structure")
    if isinstance(structure, Mapping):
        rows.extend(_structure_summary_rows(structure))
    elif values.get("has_structure") and structure is not None:
        rows.append(["Structure", _report_value(structure)])

    if source is not None and not any(row[0] == "Source" for row in rows):
        rows.append(["Source", str(source)])

    return ReportTable(title=title, columns=["Property", "Value"], rows=rows)


def _structure_summary_rows(structure: Mapping[str, Any]) -> list[list[Any]]:
    """Return concise report rows for a structural volume series."""
    rows: list[list[Any]] = [["Structure", "available"]]

    representation = structure.get("representation")
    if representation:
        rows.append(["Structure representation", str(representation)])
    orientation = structure.get("orientation")
    if orientation:
        rows.append(["Structure orientation", str(orientation)])
    reference_index = structure.get("reference_index")
    if reference_index is not None:
        rows.append(["Structure reference", f"volume index {reference_index}"])

    atomic_numbers = structure.get("atomic_numbers")
    atom_count = _sequence_length(atomic_numbers)
    if atom_count is not None:
        rows.append(["Structure atoms", atom_count])

    volume_series = structure.get("volume_series")
    if isinstance(volume_series, Mapping):
        sample_count = _sequence_length(volume_series.get("volume"))
        if sample_count is not None:
            rows.append(["Structure volume samples", sample_count])

    normalization = structure.get("normalization")
    if isinstance(normalization, Mapping):
        basis = normalization.get("basis", "unknown")
        source_basis = normalization.get("source_basis", "unknown")
        repetitions = normalization.get("repetitions")
        text = f"{basis} from {source_basis}"
        if repetitions is not None:
            text += f" ({repetitions} repetitions)"
        rows.append(["Cell normalization", text])

    symmetry = structure.get("symmetry")
    if isinstance(symmetry, Mapping):
        symbol = str(symmetry.get("international_symbol", "")).strip()
        number = symmetry.get("space_group_number")
        if symbol and number:
            rows.append(["Space group", f"{symbol} ({number})"])
        elif symbol or number:
            rows.append(["Space group", symbol or number])

    reconstruction = structure.get("reconstruction")
    if isinstance(reconstruction, Sequence) and not isinstance(
        reconstruction, (str, bytes, bytearray)
    ):
        statuses = [
            str(item.get("status", "unknown"))
            for item in reconstruction
            if isinstance(item, Mapping)
        ]
        if statuses:
            unique = sorted(set(statuses))
            rows.append(["Structure reconstruction", ", ".join(unique)])

    return rows


def _humanize_key(value: str) -> str:
    """Convert an internal mapping key to a readable label."""
    return str(value).replace("_", " ").strip().capitalize()


def _report_value(value: Any) -> Any:
    """Return a compact deterministic representation for a report cell."""
    if value is None:
        return "none"
    shape = getattr(value, "shape", None)
    if shape is not None:
        return f"array{tuple(shape)}"
    if isinstance(value, Mapping):
        if not value:
            return "empty mapping"
        if all(_is_scalar(item) for item in value.values()):
            text = ", ".join(f"{key}={item}" for key, item in value.items())
            if "\n" not in text and len(text) <= 120:
                return text
        keys = ", ".join(str(key) for key in list(value)[:6])
        suffix = ", ..." if len(value) > 6 else ""
        return f"mapping[{len(value)}] ({keys}{suffix})"
    if isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    ):
        if len(value) <= 8 and all(_is_scalar(item) for item in value):
            return ", ".join(str(item) for item in value)
        return f"sequence[{len(value)}]"
    return value


def _sequence_length(value: Any) -> int | None:
    """Return the first-axis length of a sequence-like value when available."""
    shape = getattr(value, "shape", None)
    if shape is not None and len(shape) > 0:
        return int(shape[0])
    if isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    ):
        return len(value)
    return None


def _is_scalar(value: Any) -> bool:
    """Return whether a value is suitable for compact inline rendering."""
    return value is None or isinstance(value, (str, bytes, bool, int, float, complex))


__all__ = ["ReportTable", "input_data_table", "mapping_table"]
