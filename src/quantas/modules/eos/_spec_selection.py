# -*- coding: utf-8 -*-

"""Dataset target and observation-selection resolution for EOS specs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .models import EOSDataset, EOSFitDomain
from .spec import EOSSpecError, _Entry

_TARGETS = ("volume", "a", "b", "c")


def _resolve_targets(
    entry: _Entry,
    domain: EOSFitDomain,
    dataset: EOSDataset,
    source: Path | None,
    section: str,
) -> tuple[str, ...]:
    parts = tuple(
        dict.fromkeys(
            part.strip().lower() for part in entry.value.split(",") if part.strip()
        )
    )
    if not parts:
        raise EOSSpecError(
            "targets is empty", source=source, line=entry.line, section=section
        )
    unknown = [part for part in parts if part not in {*_TARGETS, "all"}]
    if unknown:
        raise EOSSpecError(
            "unknown target(s): " + ", ".join(unknown),
            source=source,
            line=entry.line,
            section=section,
        )
    if "all" in parts and len(parts) != 1:
        raise EOSSpecError(
            "targets = all cannot be combined with explicit targets",
            source=source,
            line=entry.line,
            section=section,
        )
    if parts == ("all",):
        resolved = (
            ("volume",)
            if domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE
            else tuple(name for name in _TARGETS if dataset.has(name))
        )
    else:
        resolved = parts
    if domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE and resolved != ("volume",):
        raise EOSSpecError(
            "P-V-T fitting currently supports only targets = volume",
            source=source,
            line=entry.line,
            section=section,
        )
    missing = tuple(name for name in resolved if not dataset.has(name))
    if missing:
        raise EOSSpecError(
            "dataset does not contain target(s): " + ", ".join(missing),
            source=source,
            line=entry.line,
            section=section,
        )
    if not resolved:
        raise EOSSpecError(
            "no supported target is available in the dataset",
            source=source,
            line=entry.line,
            section=section,
        )
    return tuple(resolved)


def _resolve_selection(
    entries: dict[str, _Entry],
    dataset: EOSDataset,
    *,
    source: Path | None,
    section: str,
) -> tuple[Any, dict[str, Any]]:
    """Resolve dataset defaults, groups, and row selections for one job."""

    base_entry = entries.get("selection_base")
    base = "default" if base_entry is None else base_entry.value.strip().lower()
    if base not in {"default", "all"}:
        assert base_entry is not None
        raise EOSSpecError(
            "selection_base must be 'default' or 'all'",
            source=source,
            line=base_entry.line,
            section=section,
        )
    mask = (
        dataset.selection_mask()
        if base == "default"
        else np.ones(dataset.npoints, dtype=np.bool_)
    )

    groups_entry = entries.get("groups")
    include_groups_entry = entries.get("include_groups")
    if groups_entry is not None and include_groups_entry is not None:
        raise EOSSpecError(
            "groups and include_groups are aliases and cannot both be declared",
            source=source,
            line=include_groups_entry.line,
            section=section,
        )
    selected_groups_entry = groups_entry or include_groups_entry
    included_groups: tuple[int, ...] | None = None
    excluded_groups: tuple[int, ...] = ()
    if selected_groups_entry is not None:
        if selected_groups_entry.value.strip().lower() != "all":
            included_groups = _parse_integer_list(
                selected_groups_entry, source, section, label="groups"
            )
            _require_dataset_groups(dataset, selected_groups_entry, source, section)
            _validate_group_ids(
                included_groups, dataset, selected_groups_entry, source, section
            )
            assert dataset.groups is not None
            mask &= np.isin(dataset.groups, included_groups)
    exclude_groups_entry = entries.get("exclude_groups")
    if exclude_groups_entry is not None:
        excluded_groups = _parse_integer_list(
            exclude_groups_entry, source, section, label="exclude_groups"
        )
        _require_dataset_groups(dataset, exclude_groups_entry, source, section)
        _validate_group_ids(
            excluded_groups, dataset, exclude_groups_entry, source, section
        )
        assert dataset.groups is not None
        mask &= ~np.isin(dataset.groups, excluded_groups)

    included_rows: tuple[int, ...] | None = None
    excluded_rows: tuple[int, ...] = ()
    include_rows_entry = entries.get("include_rows")
    if include_rows_entry is not None:
        included_rows = _parse_integer_list(
            include_rows_entry, source, section, label="include_rows"
        )
        row_mask = _row_mask(
            included_rows, dataset.npoints, include_rows_entry, source, section
        )
        mask &= row_mask
    exclude_rows_entry = entries.get("exclude_rows")
    if exclude_rows_entry is not None:
        excluded_rows = _parse_integer_list(
            exclude_rows_entry, source, section, label="exclude_rows"
        )
        row_mask = _row_mask(
            excluded_rows, dataset.npoints, exclude_rows_entry, source, section
        )
        mask &= ~row_mask

    if not np.any(mask):
        location = (
            exclude_rows_entry
            or include_rows_entry
            or exclude_groups_entry
            or selected_groups_entry
            or base_entry
        )
        raise EOSSpecError(
            "data selection excludes every observation",
            source=source,
            line=None if location is None else location.line,
            section=section,
        )
    selected = int(np.count_nonzero(mask))
    metadata: dict[str, Any] = {
        "base": base,
        "total": dataset.npoints,
        "selected": selected,
        "excluded": dataset.npoints - selected,
        "include_groups": None if included_groups is None else list(included_groups),
        "exclude_groups": list(excluded_groups),
        "include_rows": None if included_rows is None else list(included_rows),
        "exclude_rows": list(excluded_rows),
    }
    if dataset.groups is not None:
        metadata["groups"] = list(dataset.group_summary(mask))
    return mask, metadata


def _parse_integer_list(
    entry: _Entry,
    source: Path | None,
    section: str,
    *,
    label: str,
) -> tuple[int, ...]:
    """Parse one comma-separated positive integer list."""
    tokens = [token.strip() for token in entry.value.split(",") if token.strip()]
    if not tokens:
        raise EOSSpecError(
            f"{label} must contain at least one positive integer",
            source=source,
            line=entry.line,
            section=section,
        )
    values: list[int] = []
    for token in tokens:
        try:
            value = int(token)
        except ValueError as exc:
            raise EOSSpecError(
                f"{label} must contain positive integers",
                source=source,
                line=entry.line,
                section=section,
            ) from exc
        if value < 1 or str(value) != token.lstrip("+"):
            raise EOSSpecError(
                f"{label} must contain positive integers",
                source=source,
                line=entry.line,
                section=section,
            )
        if value not in values:
            values.append(value)
    return tuple(values)


def _require_dataset_groups(
    dataset: EOSDataset, entry: _Entry, source: Path | None, section: str
) -> None:
    if dataset.groups is None:
        raise EOSSpecError(
            "group selection requires a GROUP column in the EOS dataset",
            source=source,
            line=entry.line,
            section=section,
        )


def _validate_group_ids(
    values: tuple[int, ...],
    dataset: EOSDataset,
    entry: _Entry,
    source: Path | None,
    section: str,
) -> None:
    unknown = tuple(value for value in values if value not in dataset.group_ids)
    if unknown:
        raise EOSSpecError(
            "unknown data group(s): " + ", ".join(str(value) for value in unknown),
            source=source,
            line=entry.line,
            section=section,
        )


def _row_mask(
    rows: tuple[int, ...],
    npoints: int,
    entry: _Entry,
    source: Path | None,
    section: str,
) -> Any:
    """Return a one-based row selection mask."""

    invalid = tuple(value for value in rows if value > npoints)
    if invalid:
        raise EOSSpecError(
            "row selection exceeds the dataset size: "
            + ", ".join(str(value) for value in invalid),
            source=source,
            line=entry.line,
            section=section,
        )
    mask = np.zeros(npoints, dtype=np.bool_)
    mask[np.asarray(rows, dtype=np.int64) - 1] = True
    return mask
