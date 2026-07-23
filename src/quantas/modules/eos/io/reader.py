# -*- coding: utf-8 -*-

"""Keyword-driven text reader for EOS datasets."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from quantas.core.physics.units import convert_length, convert_pressure, convert_volume

from quantas.models import BasicReader
from quantas.modules.eos.models import (
    EOSDataset,
    parse_eos_crystal_system,
)


_DEFAULT_UNITS = {
    "pressure": "GPa",
    "sigma_pressure": "GPa",
    "temperature": "K",
    "sigma_temperature": "K",
    "volume": "angstrom^3",
    "sigma_volume": "angstrom^3",
    "a": "angstrom",
    "sigma_a": "angstrom",
    "b": "angstrom",
    "sigma_b": "angstrom",
    "c": "angstrom",
    "sigma_c": "angstrom",
    "alpha": "degree",
    "sigma_alpha": "degree",
    "beta": "degree",
    "sigma_beta": "degree",
    "gamma": "degree",
    "sigma_gamma": "degree",
}

_COLUMN_ALIASES = {
    "p": "pressure",
    "pressure": "pressure",
    "sigmap": "sigma_pressure",
    "sigp": "sigma_pressure",
    "sigmapressure": "sigma_pressure",
    "t": "temperature",
    "temp": "temperature",
    "temperature": "temperature",
    "sigmat": "sigma_temperature",
    "sigt": "sigma_temperature",
    "sigmatemperature": "sigma_temperature",
    "v": "volume",
    "vol": "volume",
    "volume": "volume",
    "sigmav": "sigma_volume",
    "sigv": "sigma_volume",
    "sigmavolume": "sigma_volume",
    "e": "energy",
    "energy": "energy",
    "sigmae": "sigma_energy",
    "sige": "sigma_energy",
    "sigmaenergy": "sigma_energy",
    "a": "a",
    "sigmaa": "sigma_a",
    "siga": "sigma_a",
    "b": "b",
    "sigmab": "sigma_b",
    "sigb": "sigma_b",
    "c": "c",
    "sigmac": "sigma_c",
    "sigc": "sigma_c",
    "alpha": "alpha",
    "sigmaalpha": "sigma_alpha",
    "sigalpha": "sigma_alpha",
    "beta": "beta",
    "sigmabeta": "sigma_beta",
    "sigbeta": "sigma_beta",
    "gamma": "gamma",
    "sigmagamma": "sigma_gamma",
    "siggamma": "sigma_gamma",
    "group": "group",
    "grp": "group",
    "datasetgroup": "group",
    "use": "use",
    "include": "use",
    "included": "use",
}

_SELECTION_COLUMNS = frozenset({"group", "use"})

_METADATA_KEYWORDS = frozenset(
    {
        "JOB",
        "TITLE",
        "COMMENT",
        "SYSTEM",
        "TSCALE",
        "VSCALE",
        "LSCALE",
        "UNITS",
        "PROVENANCE",
        "FORMAT",
        "DATA",
    }
)

_DIMENSIONLESS_ALIASES = frozenset(
    {"1", "dimensionless", "unitless", "ratio", "relative"}
)


@dataclass(frozen=True, slots=True)
class _ParsedDataRow:
    """One mixed numeric/selection EOS input row."""

    values: tuple[float, ...]
    included: bool
    group: int | None
    excluded_by_marker: bool


class EOSInputFileReader(BasicReader[EOSDataset]):
    """Read a keyword-directed Quantas- or EosFit-style EOS text file.

    The reader accepts both the historical Quantas layout, in which ``FORMAT``
    and ``DATA`` occur on separate lines, and EosFit-style layouts in which
    column names follow ``FORMAT`` and the numeric table starts immediately.
    Keyword matching is case-insensitive and an optional trailing colon is
    accepted, for example ``FORMAT: T, V``.

    Full-line comments may start with ``#`` or with the ``COMMENT`` keyword.
    Inline text following ``#`` is ignored during numeric parsing. Comment text
    is retained in dataset metadata for provenance and future HDF5 persistence.

    Parameters
    ----------
    eos_input : str, Path or None, optional
        File loaded during construction.
    """

    def __init__(
        self,
        eos_input: str | Path | None = None,
        *,
        pressure_unit: str | None = None,
        length_unit: str | None = None,
        temperature_unit: str | None = None,
    ) -> None:
        super().__init__()
        self.dataset: EOSDataset | None = None
        self.pressure_unit = None if pressure_unit is None else str(pressure_unit)
        self.length_unit = None if length_unit is None else str(length_unit)
        self.temperature_unit = (
            None
            if temperature_unit is None
            else _normalize_temperature_scale(temperature_unit)
        )
        if eos_input is not None:
            self.load(eos_input)

    def load(self, filename: str | Path) -> EOSDataset:
        """Parse and normalize one EOS text dataset.

        Parameters
        ----------
        filename : str or Path
            Input file path.

        Returns
        -------
        EOSDataset
            Normalized table with canonical column names and ``float64`` data.

        Raises
        ------
        OSError
            If the file cannot be read.
        ValueError
            If keywords, column definitions, scales, or table values are
            invalid.
        """
        path = Path(filename)
        self.completed = False
        self.error = None
        self.dataset = None
        try:
            lines = path.read_text(encoding="utf-8").splitlines()
            dataset = self._parse(lines, source=path)
        except Exception as exc:
            self.error = str(exc)
            raise
        self.dataset = dataset
        self.completed = True
        return dataset

    def _parse(self, lines: list[str], *, source: Path) -> EOSDataset:
        """Parse text lines into an :class:`EOSDataset`."""
        jobname = "Unknown"
        system: str | None = None
        provenance: str | None = None
        temperature_scale = "K"
        volume_scale = "absolute"
        linear_scale = "absolute"
        format_columns: list[str] | None = None
        column_units: dict[str, str] = {}
        comments: list[str] = []
        rows: list[_ParsedDataRow] = []
        data_started = False
        index = 0

        while index < len(lines):
            cleaned, line_comment = _line_content(lines[index])
            if line_comment:
                comments.append(line_comment)
            if not cleaned:
                index += 1
                continue

            if data_started:
                data_parts = cleaned.split(maxsplit=1)
                if _normalize_keyword(data_parts[0]) == "COMMENT":
                    value = data_parts[1].strip() if len(data_parts) == 2 else ""
                    if not value:
                        raise ValueError(
                            f"EOS COMMENT keyword has no value on line {index + 1}."
                        )
                    comments.append(value)
                    index += 1
                    continue
                if format_columns is None:
                    raise ValueError("EOS DATA block requires a FORMAT definition.")
                rows.append(_parse_data_row(cleaned, format_columns, index + 1))
                index += 1
                continue

            keyword_parts = cleaned.split(maxsplit=1)
            keyword = _normalize_keyword(keyword_parts[0])
            remainder_text = keyword_parts[1] if len(keyword_parts) == 2 else ""
            if keyword in {
                "JOB",
                "TITLE",
                "COMMENT",
                "SYSTEM",
                "PROVENANCE",
                "TSCALE",
                "VSCALE",
                "LSCALE",
            }:
                remainder = [remainder_text] if remainder_text else []
            else:
                remainder = _split_tokens(remainder_text)

            if keyword in {"JOB", "TITLE"}:
                value, index = _keyword_value(lines, index, remainder, keyword)
                jobname = value
            elif keyword == "COMMENT":
                value, index = _keyword_value(lines, index, remainder, keyword)
                comments.append(value)
            elif keyword == "SYSTEM":
                value, index = _keyword_value(lines, index, remainder, keyword)
                system = parse_eos_crystal_system(value).value
            elif keyword == "PROVENANCE":
                value, index = _keyword_value(lines, index, remainder, keyword)
                provenance = value
            elif keyword == "TSCALE":
                value, index = _keyword_value(lines, index, remainder, keyword)
                temperature_scale = _normalize_temperature_scale(value)
            elif keyword == "VSCALE":
                value, index = _keyword_value(lines, index, remainder, keyword)
                volume_scale, unit = _parse_quantity_scale(
                    value,
                    quantity="volume",
                )
                if unit is not None:
                    column_units.update(
                        {name: unit for name in ("volume", "sigma_volume")}
                    )
            elif keyword == "LSCALE":
                value, index = _keyword_value(lines, index, remainder, keyword)
                linear_scale, unit = _parse_quantity_scale(
                    value,
                    quantity="length",
                )
                if unit is not None:
                    column_units.update(
                        {
                            name: unit
                            for name in (
                                "a",
                                "sigma_a",
                                "b",
                                "sigma_b",
                                "c",
                                "sigma_c",
                            )
                        }
                    )
            elif keyword == "UNITS":
                values, index = _keyword_tokens(lines, index, remainder, keyword)
                column_units.update(_parse_units(values))
            elif keyword == "FORMAT":
                values, index = _keyword_tokens(lines, index, remainder, keyword)
                format_columns = _parse_format(values)
            elif keyword == "DATA":
                if remainder:
                    raise ValueError("EOS DATA keyword does not accept arguments.")
                if format_columns is None:
                    raise ValueError("EOS FORMAT must be declared before DATA.")
                data_started = True
                index += 1
            elif format_columns is not None:
                data_started = True
                rows.append(_parse_data_row(cleaned, format_columns, index + 1))
                index += 1
            else:
                raise ValueError(
                    f"Unknown EOS input keyword '{keyword_parts[0]}' on line "
                    f"{index + 1}."
                )

        if format_columns is None:
            raise ValueError("EOS input does not define a FORMAT line.")
        if not rows:
            raise ValueError("EOS input does not contain any numeric data rows.")
        if self.temperature_unit is not None:
            temperature_scale = self.temperature_unit
        if self.pressure_unit is not None:
            for name in ("pressure", "sigma_pressure"):
                if name in format_columns:
                    column_units[name] = self.pressure_unit
        if self.length_unit is not None:
            for name in ("a", "sigma_a", "b", "sigma_b", "c", "sigma_c"):
                if name in format_columns:
                    column_units[name] = self.length_unit
            for name in ("volume", "sigma_volume"):
                if name in format_columns:
                    column_units[name] = f"{self.length_unit}^3"

        return _build_dataset(
            source=source,
            jobname=jobname,
            system=system,
            provenance=provenance,
            temperature_scale=temperature_scale,
            volume_scale=volume_scale,
            linear_scale=linear_scale,
            format_columns=format_columns,
            column_units=column_units,
            comments=comments,
            rows=rows,
        )


def _build_dataset(
    *,
    source: Path,
    jobname: str,
    system: str | None,
    provenance: str | None,
    temperature_scale: str,
    volume_scale: str,
    linear_scale: str,
    format_columns: list[str],
    column_units: dict[str, str],
    comments: list[str],
    rows: list[_ParsedDataRow],
) -> EOSDataset:
    """Build normalized dataset arrays and metadata from parsed input state."""
    data_columns = [name for name in format_columns if name not in _SELECTION_COLUMNS]
    matrix = np.asarray([row.values for row in rows], dtype=np.float64)
    raw_columns = {
        name: matrix[:, column_index].copy()
        for column_index, name in enumerate(data_columns)
    }
    default_mask = np.asarray([row.included for row in rows], dtype=np.bool_)
    groups = (
        None
        if "group" not in format_columns
        else np.asarray([row.group for row in rows], dtype=np.int64)
    )
    columns = {name: values.copy() for name, values in raw_columns.items()}
    raw_units = _resolve_units(
        raw_columns,
        column_units,
        temperature_scale=temperature_scale,
        normalized_temperature=False,
        volume_scale=volume_scale,
        linear_scale=linear_scale,
    )
    _normalize_physical_units(
        columns,
        raw_units,
        temperature_scale=temperature_scale,
        volume_scale=volume_scale,
        linear_scale=linear_scale,
    )
    units = _resolve_units(
        columns,
        column_units,
        temperature_scale=temperature_scale,
        normalized_temperature=True,
        volume_scale=volume_scale,
        linear_scale=linear_scale,
    )
    for name in ("pressure", "sigma_pressure"):
        if name in columns:
            units[name] = "GPa"
    for name in ("a", "sigma_a", "b", "sigma_b", "c", "sigma_c"):
        if name in columns and linear_scale == "absolute":
            units[name] = "angstrom"
    for name in ("volume", "sigma_volume"):
        if name in columns and volume_scale == "absolute":
            raw_unit = raw_units.get(name, "angstrom^3")
            units[name] = (
                "cm^3/mol" if _is_molar_volume_unit(raw_unit) else "angstrom^3"
            )
    column_provenance = (
        {name: provenance for name in columns} if provenance is not None else {}
    )
    metadata: dict[str, object] = {
        "source_format": "quantas-eos-table-v1",
        "format_columns": tuple(format_columns),
        "input_temperature_scale": temperature_scale,
        "input_volume_scale": volume_scale,
        "input_linear_scale": linear_scale,
        "column_scales": _column_scales(
            columns, volume_scale=volume_scale, linear_scale=linear_scale
        ),
        "comments": tuple(comments),
        "selection": {
            "default_selected": int(np.count_nonzero(default_mask)),
            "default_excluded": int(default_mask.size - np.count_nonzero(default_mask)),
            "excluded_by_marker": int(sum(row.excluded_by_marker for row in rows)),
            "use_column": "use" in format_columns,
            "group_column": "group" in format_columns,
        },
        "unit_overrides": {
            "pressure": raw_units.get("pressure"),
            "length": raw_units.get("a") or raw_units.get("b") or raw_units.get("c"),
            "temperature": temperature_scale if "temperature" in columns else None,
        },
    }
    if system is not None:
        crystal_system = parse_eos_crystal_system(system)
        metadata["crystal_system"] = crystal_system.value
        metadata["independent_cell_axes"] = crystal_system.independent_axes
    return EOSDataset(
        jobname=jobname,
        columns=columns,
        units=units,
        raw_columns=raw_columns,
        raw_units=raw_units,
        provenance=column_provenance,
        source=source,
        metadata=metadata,
        default_mask=default_mask,
        groups=groups,
    )


def read_eos_input(
    filename: str | Path,
    *,
    pressure_unit: str | None = None,
    length_unit: str | None = None,
    temperature_unit: str | None = None,
) -> EOSDataset:
    """Read and normalize one keyword-directed EOS input file.

    Parameters
    ----------
    filename : str or Path
        Input file path. File extensions are not used to select the parser.
    pressure_unit, length_unit, temperature_unit : str or None, optional
        Explicit input-unit overrides. When omitted, declarations in the file
        are used, followed by the EOS defaults GPa, Angstrom, and kelvin.
        Normalized in-memory values always use GPa, Angstrom/Angstrom^3, and K.

    Returns
    -------
    EOSDataset
        Raw input columns together with normalized ``float64`` data.

    Raises
    ------
    OSError
        If the file cannot be opened.
    ValueError
        If the input or unit declarations are invalid.
    """
    return EOSInputFileReader(
        pressure_unit=pressure_unit,
        length_unit=length_unit,
        temperature_unit=temperature_unit,
    ).load(filename)


def _normalize_token(token: str) -> str:
    """Normalize a format or unit field name for alias lookup."""
    return re.sub(r"[^a-z0-9]", "", token.lower())


def _normalize_keyword(token: str) -> str:
    """Normalize one case-insensitive keyword and optional trailing colon."""
    return token.strip().rstrip(":").upper()


def _canonical_column(token: str) -> str:
    """Return the canonical EOS column name for one alias."""
    normalized = _normalize_token(token)
    try:
        return _COLUMN_ALIASES[normalized]
    except KeyError as exc:
        raise ValueError(f"Unsupported EOS FORMAT column: {token}") from exc


def _line_content(line: str) -> tuple[str, str | None]:
    """Return parseable content and optional ``#`` comment text."""
    content, separator, comment = line.partition("#")
    comment_text = comment.strip() if separator and comment.strip() else None
    return content.strip(), comment_text


def _clean_line(line: str) -> str:
    """Strip ``#`` comments and surrounding whitespace from one line."""
    return _line_content(line)[0]


def _split_tokens(line: str) -> list[str]:
    """Split comma- or whitespace-delimited input text."""
    return line.replace(",", " ").split()


def _next_content_line(lines: list[str], index: int, keyword: str) -> tuple[str, int]:
    """Return the next non-comment line and the following index."""
    candidate = index + 1
    while candidate < len(lines):
        cleaned = _clean_line(lines[candidate])
        if cleaned:
            first = _normalize_keyword(_split_tokens(cleaned)[0])
            if first in _METADATA_KEYWORDS:
                raise ValueError(f"EOS {keyword} keyword has no value.")
            return cleaned, candidate + 1
        candidate += 1
    raise ValueError(f"EOS {keyword} keyword has no value.")


def _keyword_value(
    lines: list[str],
    index: int,
    remainder: list[str],
    keyword: str,
) -> tuple[str, int]:
    """Read a free-text keyword value from the same or next content line."""
    if remainder:
        return " ".join(remainder), index + 1
    value, next_index = _next_content_line(lines, index, keyword)
    return value, next_index


def _keyword_tokens(
    lines: list[str],
    index: int,
    remainder: list[str],
    keyword: str,
) -> tuple[list[str], int]:
    """Read tokenized keyword arguments from the same or next content line."""
    if remainder:
        return remainder, index + 1
    value, next_index = _next_content_line(lines, index, keyword)
    return _split_tokens(value), next_index


def _parse_format(tokens: list[str]) -> list[str]:
    """Normalize a FORMAT declaration and reject duplicate quantities."""
    local_tokens = list(tokens)
    if local_tokens and local_tokens[0].isdigit():
        lines_per_record = int(local_tokens.pop(0))
        if lines_per_record != 1:
            raise ValueError("EOS currently supports one numeric line per data record.")
    if not local_tokens:
        raise ValueError("EOS FORMAT must define at least one data column.")
    columns = [_canonical_column(token) for token in local_tokens]
    duplicates = sorted({name for name in columns if columns.count(name) > 1})
    if not any(name not in _SELECTION_COLUMNS for name in columns):
        raise ValueError("EOS FORMAT must define at least one scientific data column.")
    if duplicates:
        duplicate_text = ", ".join(duplicates)
        raise ValueError(f"EOS FORMAT defines duplicate columns: {duplicate_text}")
    return columns


def _parse_units(tokens: list[str]) -> dict[str, str]:
    """Parse ``COLUMN=UNIT`` declarations."""
    if not tokens:
        raise ValueError("EOS UNITS must contain COLUMN=UNIT declarations.")
    units: dict[str, str] = {}
    for token in tokens:
        if "=" not in token:
            raise ValueError("EOS UNITS entries must use COLUMN=UNIT syntax.")
        column_token, unit = token.split("=", maxsplit=1)
        if not unit:
            raise ValueError(f"EOS UNITS entry has no unit: {token}")
        column = _canonical_column(column_token)
        if column in _SELECTION_COLUMNS:
            raise ValueError(f"EOS UNITS cannot assign a unit to {column!r}.")
        if column in units:
            raise ValueError(f"EOS UNITS defines '{column}' more than once.")
        units[column] = unit
    return units


def _parse_data_row(
    line: str,
    format_columns: list[str],
    line_number: int,
) -> _ParsedDataRow:
    """Parse one mixed EOS table row and non-destructive selection markers."""
    tokens = _split_tokens(line)
    excluded_by_marker = bool(tokens and tokens[-1] == "*")
    if excluded_by_marker:
        tokens = tokens[:-1]
    width = len(format_columns)
    if len(tokens) != width:
        suffix = " plus an optional trailing '*'"
        raise ValueError(
            f"EOS data row {line_number} has {len(tokens)} values; expected "
            f"{width}{suffix}."
        )
    values: list[float] = []
    included = True
    group: int | None = None
    for column, token in zip(format_columns, tokens, strict=True):
        if column == "use":
            included = _parse_use_token(token, line_number)
            continue
        if column == "group":
            group = _parse_group_token(token, line_number)
            continue
        try:
            values.append(float(token))
        except ValueError as exc:
            raise ValueError(
                f"EOS data row {line_number} contains a non-numeric value "
                f"for column {column!r}."
            ) from exc
    if "group" in format_columns and group is None:  # pragma: no cover
        raise RuntimeError("EOS group column was not parsed.")
    return _ParsedDataRow(
        values=tuple(values),
        included=included and not excluded_by_marker,
        group=group,
        excluded_by_marker=excluded_by_marker,
    )


def _parse_use_token(token: str, line_number: int) -> bool:
    """Return the boolean meaning of one USE-column token."""
    normalized = token.strip().lower()
    if normalized in {"1", "yes", "true", "on", "include", "included"}:
        return True
    if normalized in {"0", "no", "false", "off", "exclude", "excluded"}:
        return False
    raise ValueError(f"EOS data row {line_number} USE value must be 1/0 or yes/no.")


def _parse_group_token(token: str, line_number: int) -> int:
    """Return one positive integer group identifier."""
    try:
        numeric = float(token)
    except ValueError as exc:
        raise ValueError(
            f"EOS data row {line_number} GROUP value must be a positive integer."
        ) from exc
    integer = int(numeric)
    if numeric != integer or integer < 1:
        raise ValueError(
            f"EOS data row {line_number} GROUP value must be a positive integer."
        )
    return integer


def _normalize_temperature_scale(value: str) -> str:
    """Return canonical K, C, or F temperature-scale notation."""
    normalized = value.strip().lower().replace("°", "")
    aliases = {
        "k": "K",
        "kelvin": "K",
        "c": "C",
        "celsius": "C",
        "f": "F",
        "fahrenheit": "F",
    }
    try:
        return aliases[normalized]
    except KeyError as exc:
        raise ValueError(f"Unsupported EOS temperature scale: {value}") from exc


def _parse_quantity_scale(
    value: str,
    *,
    quantity: str,
) -> tuple[str, str | None]:
    """Normalize a scale declaration and optional legacy unit shorthand.

    ``VSCALE`` and ``LSCALE`` historically accepted either a scale descriptor
    or the physical unit of an absolute quantity.  Quantas preserves that
    useful input convention while storing scale and unit as separate metadata.

    Parameters
    ----------
    value : str
        Keyword value supplied by the input file.
    quantity : {"volume", "length"}
        Quantity governed by the scale declaration.

    Returns
    -------
    tuple of (str, str or None)
        Canonical scale and an optional unit declaration to apply to the
        corresponding value and uncertainty columns.

    Raises
    ------
    ValueError
        If ``quantity`` is unsupported or a normalized scale is malformed.
    """
    stripped = value.strip()
    normalized = stripped.lower().replace(" ", "")
    absolute_aliases = {"absolute", "none", "raw"}
    if quantity == "volume":
        normalized_aliases = {"v/v0", "v/vref", "normalized", "relative"}
        canonical = "V/V0"
    elif quantity == "length":
        normalized_aliases = {"l/l0", "l/lref", "normalized", "relative"}
        canonical = "L/L0"
    else:  # pragma: no cover - private call contract
        raise ValueError(f"Unsupported EOS scale quantity: {quantity}")
    if normalized in absolute_aliases:
        return "absolute", None
    if normalized in normalized_aliases:
        return canonical, None
    if "/" in normalized and normalized.startswith(("v/", "l/")):
        keyword = "VSCALE" if quantity == "volume" else "LSCALE"
        raise ValueError(f"Unsupported EOS {keyword} value: {value}")
    return "absolute", stripped


def _resolve_units(
    columns: dict[str, np.ndarray],
    declared_units: dict[str, str],
    *,
    temperature_scale: str,
    normalized_temperature: bool,
    volume_scale: str,
    linear_scale: str,
) -> dict[str, str]:
    """Resolve raw or normalized units while enforcing scale consistency."""
    units = {
        name: declared_units.get(name, _DEFAULT_UNITS[name])
        for name in columns
        if name in declared_units or name in _DEFAULT_UNITS
    }
    if "temperature" in columns:
        units["temperature"] = "K" if normalized_temperature else temperature_scale
    if "sigma_temperature" in columns:
        units["sigma_temperature"] = (
            "K" if normalized_temperature else temperature_scale
        )
    if volume_scale != "absolute":
        _set_dimensionless_units(
            units,
            declared_units,
            columns,
            ("volume", "sigma_volume"),
            keyword="VSCALE",
        )
    if linear_scale != "absolute":
        _set_dimensionless_units(
            units,
            declared_units,
            columns,
            ("a", "sigma_a", "b", "sigma_b", "c", "sigma_c"),
            keyword="LSCALE",
        )
    return units


def _set_dimensionless_units(
    units: dict[str, str],
    declared_units: dict[str, str],
    columns: dict[str, np.ndarray],
    names: tuple[str, ...],
    *,
    keyword: str,
) -> None:
    """Assign dimensionless units and reject contradictory declarations."""
    for name in names:
        if name not in columns:
            continue
        declared = declared_units.get(name)
        if (
            declared is not None
            and declared.strip().lower() not in _DIMENSIONLESS_ALIASES
        ):
            raise ValueError(
                f"EOS {keyword} declares normalized data but UNITS assigns "
                f"'{declared}' to '{name}'."
            )
        units[name] = "dimensionless"


def _column_scales(
    columns: dict[str, np.ndarray],
    *,
    volume_scale: str,
    linear_scale: str,
) -> dict[str, str]:
    """Return the declared scale associated with each stored data column."""
    scales: dict[str, str] = {}
    for name in ("volume", "sigma_volume"):
        if name in columns:
            scales[name] = volume_scale
    for name in ("a", "sigma_a", "b", "sigma_b", "c", "sigma_c"):
        if name in columns:
            scales[name] = linear_scale
    return scales


def _normalize_physical_units(
    columns: dict[str, np.ndarray],
    raw_units: dict[str, str],
    *,
    temperature_scale: str,
    volume_scale: str,
    linear_scale: str,
) -> None:
    """Convert supported physical input columns to EOS internal units."""
    _convert_temperature(columns, temperature_scale)
    for name in ("pressure", "sigma_pressure"):
        if name in columns:
            source_unit = raw_units.get(name, "GPa")
            if not _unit_is(source_unit, {"gpa", "gigapascal"}):
                columns[name] = np.asarray(
                    convert_pressure(columns[name], source_unit, "GPa"),
                    dtype=np.float64,
                )
    if linear_scale == "absolute":
        for name in ("a", "sigma_a", "b", "sigma_b", "c", "sigma_c"):
            if name in columns:
                source_unit = _length_unit_base(raw_units.get(name, "angstrom"))
                if not _unit_is(source_unit, {"angstrom", "ang", "a", "å"}):
                    columns[name] = np.asarray(
                        convert_length(columns[name], source_unit, "angstrom"),
                        dtype=np.float64,
                    )
    if volume_scale == "absolute":
        for name in ("volume", "sigma_volume"):
            if name in columns:
                raw_unit = raw_units.get(name, "angstrom^3")
                if _is_molar_volume_unit(raw_unit):
                    continue
                source_unit = _length_unit_base(raw_unit)
                if not _unit_is(source_unit, {"angstrom", "ang", "a", "å"}):
                    columns[name] = np.asarray(
                        convert_volume(columns[name], source_unit, "angstrom"),
                        dtype=np.float64,
                    )



def _is_molar_volume_unit(unit: str) -> bool:
    """Return whether a unit denotes cubic centimetres per mole."""
    normalized = (
        str(unit)
        .strip()
        .lower()
        .replace(" ", "")
        .replace("³", "^3")
        .replace("mol⁻¹", "/mol")
        .replace("mol^-1", "/mol")
        .replace("mol-1", "/mol")
    )
    return normalized in {
        "cm^3/mol",
        "cm3/mol",
        "cc/mol",
        "cm^3permol",
        "cm3permol",
    }

def _unit_is(unit: str, aliases: set[str]) -> bool:
    """Return whether a unit label is already one of the target aliases."""
    return str(unit).strip().lower().replace(" ", "") in aliases


def _length_unit_base(unit: str) -> str:
    """Return the length unit defining a length or cubic-volume label."""
    normalized = str(unit).strip()
    lowered = normalized.lower().replace(" ", "")
    for suffix in ("^3", "**3", "³", "3"):
        if lowered.endswith(suffix):
            normalized = normalized[: -len(suffix)]
            break
    return normalized.rstrip("*/ ")


def _convert_temperature(columns: dict[str, np.ndarray], scale: str) -> None:
    """Convert temperature values and uncertainties to Kelvin in place."""
    if "temperature" in columns:
        if scale == "C":
            columns["temperature"] = columns["temperature"] + 273.15
        elif scale == "F":
            columns["temperature"] = (columns["temperature"] - 32.0) * (
                5.0 / 9.0
            ) + 273.15
    if "sigma_temperature" in columns and scale == "F":
        columns["sigma_temperature"] = columns["sigma_temperature"] * (5.0 / 9.0)


__all__ = ["EOSInputFileReader", "read_eos_input"]
