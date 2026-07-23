# -*- coding: utf-8 -*-

"""YAML reader for Quantas quasi-harmonic approximation input files."""

from __future__ import annotations

from pathlib import Path
from typing import cast

from quantas.io.phonons import PhononInputFileReader
from quantas.models.phonons import PhononInputData
from quantas.modules.qha.models import QHAInput, QHAModeContinuity


def _mode_continuity(value: object) -> QHAModeContinuity:
    """Return one validated mode-continuity contract value."""
    normalized = str(value).strip().lower()
    allowed = {"verified", "assumed", "unknown", "unreliable"}
    if normalized not in allowed:
        raise ValueError(
            "mode_continuity must be verified, assumed, unknown, or unreliable"
        )
    return cast(QHAModeContinuity, normalized)


class QHAInputFileReader(PhononInputFileReader):
    """Reader for Quantas QHA phonon YAML input files."""

    def to_qha_input(self, source: str | Path | None = None) -> QHAInput:
        """Convert parsed YAML data to a normalized QHA input container.

        Parameters
        ----------
        source : str or Path or None, optional
            Source path stored in the returned input object.

        Returns
        -------
        QHAInput
            Normalized quasi-harmonic approximation input data.

        Raises
        ------
        ValueError
            If the reader has not successfully loaded a valid input file.
        """
        data = super().to_input(source=source)
        raw = self.data or {}
        continuity = _mode_continuity(
            raw.get(
                "mode_continuity",
                raw.get("phonon_mode_continuity", "assumed"),
            )
        )
        metadata = dict(data.metadata)
        metadata["mode_continuity"] = continuity
        metadata["source_format"] = "quantas-phonon-yaml"
        return QHAInput(
            jobname=data.jobname,
            natoms=data.natoms,
            formula_units=data.formula_units,
            supercell=data.supercell,
            qpoints=data.qpoints,
            volume=data.volume,
            energy=data.energy,
            frequencies=data.frequencies,
            weights=data.weights,
            qcoords=data.qcoords,
            structure=data.structure,
            source=data.source,
            metadata=metadata,
            mode_continuity=continuity,
        )


def phonon_to_qha_input(input_data: PhononInputData) -> QHAInput:
    """Convert normalized phonon input data to a QHA input container.

    Parameters
    ----------
    input_data : PhononInputData
        Normalized volume-dependent phonon data.

    Returns
    -------
    QHAInput
        Quasi-harmonic input data using the same arrays and metadata.
    """
    metadata = dict(input_data.metadata)
    continuity = _mode_continuity(metadata.get("mode_continuity", "assumed"))
    metadata.setdefault("format", "quantas-phonon-yaml")
    metadata["source_format"] = "quantas-phonon-yaml"
    return QHAInput(
        jobname=input_data.jobname,
        natoms=int(input_data.natoms),
        formula_units=int(input_data.formula_units),
        supercell=input_data.supercell,
        qpoints=int(input_data.qpoints),
        volume=input_data.volume,
        energy=input_data.energy,
        frequencies=input_data.frequencies,
        weights=input_data.weights,
        qcoords=input_data.qcoords,
        structure=input_data.structure,
        source=input_data.source,
        metadata=metadata,
        mode_continuity=continuity,
    )


def read_qha_input(filename: str | Path) -> QHAInput:
    """Read a Quantas QHA phonon YAML input file.

    Parameters
    ----------
    filename : str or Path
        Path to the YAML input file.

    Returns
    -------
    QHAInput
        Normalized QHA input data.

    Raises
    ------
    ValueError
        If the input file cannot be parsed or validated.
    """
    reader = QHAInputFileReader(filename)
    if not reader.completed:
        raise ValueError(reader.error or "Unable to read QHA input file")
    return reader.to_qha_input(source=Path(filename))
