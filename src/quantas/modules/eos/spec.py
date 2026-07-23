# -*- coding: utf-8 -*-

"""Public contracts for strict declarative EOS batch specifications.

The EOS specification format is a small, INI-like language designed for
hand-written scientific batch plans. It is independent of filename suffixes,
Click, terminal rendering, and HDF5. Syntax parsing and dataset-dependent
resolution are implemented in private companion modules so this module remains
a compact public facade.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .batch import EOSBatchPlan
from .models import EOSDataset
from .report import EOSReportOptions

EOS_SPEC_SIGNATURE = "# QUANTAS EOS SPEC 1"
EOS_SPEC_VERSION = 1


class EOSSpecError(ValueError):
    """Error raised for invalid EOS specification syntax or semantics.

    Parameters
    ----------
    message : str
        Human-readable explanation.
    source : path-like or None, optional
        Specification source shown in the diagnostic.
    line : int or None, optional
        One-based source line associated with the error.
    section : str or None, optional
        Section active when the error was detected.
    """

    def __init__(
        self,
        message: str,
        *,
        source: str | Path | None = None,
        line: int | None = None,
        section: str | None = None,
    ) -> None:
        self.message = str(message)
        self.source = None if source is None else Path(source)
        self.line = line
        self.section = section
        location = ""
        if self.source is not None:
            location = str(self.source)
        if line is not None:
            location = f"{location}:{line}" if location else f"line {line}"
        if section is not None:
            location = f"{location} [{section}]" if location else f"[{section}]"
        super().__init__(f"{location}: {self.message}" if location else self.message)


@dataclass(frozen=True, slots=True)
class EOSSpecInputOptions:
    """Input-unit overrides declared by an EOS specification.

    Parameters
    ----------
    pressure_unit, length_unit, temperature_unit : str or None
        Units passed to the EOS input reader before plan resolution.
    """

    pressure_unit: str | None = None
    length_unit: str | None = None
    temperature_unit: str | None = None

    def as_dict(self) -> dict[str, str | None]:
        """Return a serialization-ready mapping."""
        return {
            "pressure_unit": self.pressure_unit,
            "length_unit": self.length_unit,
            "temperature_unit": self.temperature_unit,
        }


@dataclass(frozen=True, slots=True)
class _Entry:
    """One normalized key-value entry retained with its source line."""

    key: str
    value: str
    line: int


@dataclass(frozen=True, slots=True)
class _Section:
    """One parsed specification section retained in source order."""

    name: str
    display_name: str
    line: int
    entries: tuple[_Entry, ...]

    def mapping(self) -> dict[str, _Entry]:
        """Return entries indexed by normalized key."""
        return {entry.key: entry for entry in self.entries}


@dataclass(frozen=True, slots=True)
class EOSSpecDocument:
    """Parsed but dataset-independent EOS specification document.

    Parameters
    ----------
    source : Path or None
        Origin of the specification.
    source_text : str
        Complete original text.
    version : int
        EOS specification schema version.
    metadata : dict
        Free-form provenance from ``[metadata]``.
    input_options : EOSSpecInputOptions
        Unit overrides available before reading the data file.
    source_sha256 : str
        SHA-256 digest of ``source_text`` encoded as UTF-8.
    """

    source: Path | None
    source_text: str
    version: int
    metadata: dict[str, str]
    input_options: EOSSpecInputOptions
    source_sha256: str
    _sections: tuple[_Section, ...] = field(repr=False)

    def as_dict(self) -> dict[str, Any]:
        """Return document provenance without unresolved job internals."""
        return {
            "source": None if self.source is None else str(self.source),
            "version": self.version,
            "metadata": dict(self.metadata),
            "input": self.input_options.as_dict(),
            "source_sha256": self.source_sha256,
        }


@dataclass(frozen=True, slots=True)
class EOSResolvedSpec:
    """Dataset-resolved EOS specification.

    Parameters
    ----------
    document : EOSSpecDocument
        Parsed source document.
    plan : EOSBatchPlan
        Fully typed and ordered batch plan.
    report_options : EOSReportOptions
        Frontend-neutral presentation preferences.
    """

    document: EOSSpecDocument
    plan: EOSBatchPlan
    report_options: EOSReportOptions


def read_eos_spec(path: str | Path) -> EOSSpecDocument:
    """Read and parse an EOS specification from any regular text file.

    Parameters
    ----------
    path : str or Path
        Specification path. The suffix has no semantic meaning.

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
    from ._spec_parser import read_eos_spec_document

    return read_eos_spec_document(path)


def parse_eos_spec(
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
    from ._spec_parser import parse_eos_spec_document

    return parse_eos_spec_document(text, source=source)


def resolve_eos_spec(
    document: EOSSpecDocument,
    dataset: EOSDataset,
) -> EOSResolvedSpec:
    """Resolve defaults, targets, models, constraints, and solver contracts.

    Parameters
    ----------
    document : EOSSpecDocument
        Parsed EOS specification.
    dataset : EOSDataset
        Normalized dataset used to expand ``targets = all`` and validate target
        availability.

    Returns
    -------
    EOSResolvedSpec
        Typed batch plan and report options.

    Raises
    ------
    EOSSpecError
        If inherited settings are inconsistent or unavailable for ``dataset``.
    """
    from ._spec_resolver import resolve_eos_spec_document

    return resolve_eos_spec_document(document, dataset)


def canonical_eos_parameter_name(value: str) -> str:
    """Return the documented canonical EOS parameter name for one alias.

    Parameters
    ----------
    value : str
        User-facing parameter name or supported alias.

    Returns
    -------
    str
        Canonical Quantas parameter name when an alias is known.
    """
    from ._spec_resolver import canonical_parameter_name

    return canonical_parameter_name(value)


__all__ = [
    "EOS_SPEC_SIGNATURE",
    "EOS_SPEC_VERSION",
    "EOSResolvedSpec",
    "EOSSpecDocument",
    "EOSSpecError",
    "EOSSpecInputOptions",
    "canonical_eos_parameter_name",
    "parse_eos_spec",
    "read_eos_spec",
    "resolve_eos_spec",
]
