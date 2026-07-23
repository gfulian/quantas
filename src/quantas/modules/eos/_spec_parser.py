# -*- coding: utf-8 -*-

"""Lexical and structural parser for EOS specification text."""

from __future__ import annotations

import hashlib
from pathlib import Path
import re

from .spec import (
    EOS_SPEC_SIGNATURE,
    EOS_SPEC_VERSION,
    EOSSpecDocument,
    EOSSpecError,
    EOSSpecInputOptions,
    _Entry,
    _Section,
)

_ALLOWED_SECTIONS = {
    "metadata",
    "input",
    "batch",
    "defaults",
    "defaults.pv",
    "defaults.vt",
    "defaults.pvt",
    "presentation",
}
_INPUT_KEYS = {"pressure_unit", "length_unit", "temperature_unit"}
_BATCH_KEYS = {"failure_policy"}
_PRESENTATION_KEYS = {"detail", "show_uncertainties", "max_data_rows"}
_SELECTION_KEYS = {
    "selection_base",
    "groups",
    "include_groups",
    "exclude_groups",
    "include_rows",
    "exclude_rows",
}
_COMMON_FIT_KEYS = {
    "solver",
    "covariance_scaling",
    "max_iterations",
    "ftol",
    "xtol",
    "gtol",
    "inner_max_iterations",
    "odr_difference",
    "odr_ndigit",
    "accept",
    "replace_accepted",
} | _SELECTION_KEYS
_MGD_MODEL_KEYS = {
    "thermal_pressure_model",
    "volume_basis",
    "atoms_per_cell",
    "atoms_per_formula_unit",
    "formula",
    "formula_units_per_cell",
}
_JOB_KEYS = (
    _COMMON_FIT_KEYS
    | {
        "domain",
        "targets",
        "model",
        "pv_model",
        "vt_model",
        "coupling",
        "note",
    }
    | _MGD_MODEL_KEYS
)
_DOMAIN_KEYS = {
    "defaults.pv": _COMMON_FIT_KEYS | {"model"},
    "defaults.vt": _COMMON_FIT_KEYS | {"model"},
    "defaults.pvt": (
        _COMMON_FIT_KEYS | {"pv_model", "vt_model", "coupling"} | _MGD_MODEL_KEYS
    ),
}
_PARAMETER_PREFIXES = ("fix.", "initial.", "bound.")


def read_eos_spec_document(path: str | Path) -> EOSSpecDocument:
    """Read and parse an EOS specification from any regular text file.

    Parameters
    ----------
    path : str or Path
        Specification path.  The suffix has no semantic meaning.

    Returns
    -------
    EOSSpecDocument
        Parsed document with input-unit declarations available immediately.

    Raises
    ------
    EOSSpecError
        If the signature, sections, keys, or scalar values are invalid.
    OSError
        If the file cannot be read.
    """
    source = Path(path)
    text = source.read_text(encoding="utf-8-sig")
    return parse_eos_spec_document(text, source=source)


def parse_eos_spec_document(
    text: str,
    *,
    source: str | Path | None = None,
) -> EOSSpecDocument:
    """Parse EOS specification text without requiring a dataset.

    Parameters
    ----------
    text : str
        Complete specification text.
    source : str or Path or None, optional
        Diagnostic source label.

    Returns
    -------
    EOSSpecDocument
        Strictly parsed document.
    """
    source_path = None if source is None else Path(source)
    sections = _parse_sections(str(text), source_path)
    by_name = {
        section.name: section
        for section in sections
        if not section.name.startswith("job ")
    }
    metadata = {
        entry.key: entry.value
        for entry in by_name.get(
            "metadata", _Section("metadata", "metadata", 0, ())
        ).entries
    }
    input_options = _parse_input_options(by_name.get("input"), source_path)
    return EOSSpecDocument(
        source=source_path,
        source_text=str(text),
        version=EOS_SPEC_VERSION,
        metadata=metadata,
        input_options=input_options,
        source_sha256=hashlib.sha256(str(text).encode("utf-8")).hexdigest(),
        _sections=sections,
    )


def _parse_sections(text: str, source: Path | None) -> tuple[_Section, ...]:
    lines = text.splitlines()
    signature_index = next(
        (index for index, line in enumerate(lines) if line.strip()), None
    )
    if signature_index is None:
        raise EOSSpecError("empty EOS specification", source=source)
    if lines[signature_index].strip() != EOS_SPEC_SIGNATURE:
        raise EOSSpecError(
            f"first non-empty line must be exactly {EOS_SPEC_SIGNATURE!r}",
            source=source,
            line=signature_index + 1,
        )

    sections: list[_Section] = []
    current_name: str | None = None
    current_display: str | None = None
    current_line = 0
    current_entries: list[_Entry] = []
    seen_sections: set[str] = set()

    def finish() -> None:
        nonlocal current_name, current_display, current_line, current_entries
        if current_name is not None and current_display is not None:
            sections.append(
                _Section(
                    current_name, current_display, current_line, tuple(current_entries)
                )
            )
        current_name = None
        current_display = None
        current_line = 0
        current_entries = []

    for number, raw in enumerate(
        lines[signature_index + 1 :], start=signature_index + 2
    ):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        section_match = re.fullmatch(r"\[([^\]]+)\]", stripped)
        if section_match is not None:
            finish()
            display = " ".join(section_match.group(1).strip().split())
            normalized = display.lower()
            if normalized.startswith("job "):
                job_name = display[4:].strip()
                if not job_name:
                    raise EOSSpecError(
                        "job section requires a name", source=source, line=number
                    )
                normalized = f"job {job_name.lower()}"
            elif normalized not in _ALLOWED_SECTIONS:
                raise EOSSpecError(
                    f"unknown section [{display}]",
                    source=source,
                    line=number,
                )
            if normalized in seen_sections:
                raise EOSSpecError(
                    f"duplicate section [{display}]",
                    source=source,
                    line=number,
                )
            seen_sections.add(normalized)
            current_name = normalized
            current_display = display
            current_line = number
            continue
        if current_name is None:
            raise EOSSpecError(
                "key-value entry appears before any section",
                source=source,
                line=number,
            )
        if "=" not in raw:
            raise EOSSpecError(
                "entries must use KEY = VALUE",
                source=source,
                line=number,
                section=current_display,
            )
        raw_key, raw_value = raw.split("=", 1)
        key = _normalize_key(raw_key)
        value = raw_value.strip()
        if not key:
            raise EOSSpecError(
                "empty key", source=source, line=number, section=current_display
            )
        if not value:
            raise EOSSpecError(
                f"key {raw_key.strip()!r} has an empty value",
                source=source,
                line=number,
                section=current_display,
            )
        if any(entry.key == key for entry in current_entries):
            raise EOSSpecError(
                f"duplicate key {raw_key.strip()!r}",
                source=source,
                line=number,
                section=current_display,
            )
        assert current_display is not None
        _validate_key(current_name, key, source, number, current_display)
        current_entries.append(_Entry(key, value, number))
    finish()
    return tuple(sections)


def _normalize_key(value: str) -> str:
    text = value.strip()
    lower = text.lower()
    for prefix in _PARAMETER_PREFIXES:
        if lower.startswith(prefix):
            return prefix + text[len(prefix) :].strip()
    return lower


def _validate_key(
    section: str,
    key: str,
    source: Path | None,
    line: int,
    display: str,
) -> None:
    if section == "metadata":
        return
    if section == "input":
        allowed = _INPUT_KEYS
    elif section == "batch":
        allowed = _BATCH_KEYS
    elif section == "presentation":
        allowed = _PRESENTATION_KEYS
    elif section == "defaults":
        allowed = _COMMON_FIT_KEYS
    elif section in _DOMAIN_KEYS:
        allowed = _DOMAIN_KEYS[section]
    else:
        allowed = _JOB_KEYS
    if key in allowed:
        return
    if key.lower().startswith(_PARAMETER_PREFIXES):
        if section in {
            "defaults",
            "defaults.pv",
            "defaults.vt",
            "defaults.pvt",
        } or section.startswith("job "):
            suffix = key.split(".", 1)[1].strip()
            if suffix:
                return
    raise EOSSpecError(
        f"unknown key {key!r}",
        source=source,
        line=line,
        section=display,
    )


def _parse_input_options(
    section: _Section | None, source: Path | None
) -> EOSSpecInputOptions:
    if section is None:
        return EOSSpecInputOptions()
    values = section.mapping()
    temperature = values.get("temperature_unit")
    if temperature is not None and temperature.value.strip().upper() not in {
        "K",
        "C",
        "F",
    }:
        raise EOSSpecError(
            "temperature_unit must be K, C, or F",
            source=source,
            line=temperature.line,
            section=section.display_name,
        )
    return EOSSpecInputOptions(
        pressure_unit=None
        if values.get("pressure_unit") is None
        else values["pressure_unit"].value,
        length_unit=None
        if values.get("length_unit") is None
        else values["length_unit"].value,
        temperature_unit=None if temperature is None else temperature.value.upper(),
    )
