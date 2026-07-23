# -*- coding: utf-8 -*-

"""Input and HDF5 readers for the harmonic-approximation module."""

from __future__ import annotations

from pathlib import Path


from quantas.io.hdf5 import (
    read_events,
    read_input_data,
    read_options,
    read_result_metadata,
    read_warnings,
)
from quantas.io.phonons import PhononInputFileReader
from quantas.models import (
    RESULT_SCHEMA_VERSION,
    ResultData,
    ResultMetadata,
    validate_result_schema,
)
from quantas.models.phonons import PhononInputData
from quantas.modules.ha.io.hdf5_payload import (
    migrate_schema_1_0_thermodynamic_units,
    read_current_ha_payload,
    read_historical_ha_payload,
)
from quantas.modules.ha.models import HAInput


class HAInputFileReader(PhononInputFileReader):
    """Reader for Quantas harmonic phonon YAML input files."""

    def to_input(self, source: str | Path | None = None) -> HAInput:
        """Convert parsed YAML data to a normalized HA input container.

        Parameters
        ----------
        source : str or Path or None, optional
            Source path stored in the returned input object.

        Returns
        -------
        HAInput
            Normalized harmonic-approximation input data.

        Raises
        ------
        ValueError
            If the reader has not successfully loaded a valid input file.
        """
        data = super().to_input(source=source)
        return HAInput(
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
            metadata=dict(data.metadata),
        )


def phonon_to_ha_input(input_data: PhononInputData) -> HAInput:
    """Convert normalized phonon data to a harmonic input container.

    Parameters
    ----------
    input_data : PhononInputData
        Normalized volume-dependent phonon data.

    Returns
    -------
    HAInput
        Harmonic-approximation input data using the same arrays and metadata.
    """
    return HAInput(
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
        metadata=dict(input_data.metadata),
    )


def read_ha_input(filename: str | Path) -> HAInput:
    """Read a Quantas harmonic phonon YAML input file.

    Parameters
    ----------
    filename : str or Path
        Path to the YAML input file.

    Returns
    -------
    HAInput
        Normalized harmonic-approximation input data.

    Raises
    ------
    ValueError
        If the input file cannot be parsed or validated.
    """
    reader = HAInputFileReader(filename)
    if not reader.completed:
        raise ValueError(reader.error or "Unable to read HA input file")
    return reader.to_input(source=Path(filename))


def read_ha_hdf5(filename: str | Path) -> ResultData:
    """Read a complete Quantas harmonic HDF5 result.

    Parameters
    ----------
    filename : str or Path
        Native Quantas HA HDF5 file. Historical HA files are accepted through
        an explicit migration path.

    Returns
    -------
    ResultData
        Generic result containing a :class:`HAResult` under the ``"ha"`` key.

    Raises
    ------
    ValueError
        If the file does not contain recognizable harmonic results.
    """
    import h5py

    path = Path(filename)
    with h5py.File(path, "r") as h5:
        if "metadata" in h5:
            metadata = read_result_metadata(h5)
            if metadata.schema_version == RESULT_SCHEMA_VERSION:
                validate_result_schema(h5, expected_module="ha")
                payload = read_current_ha_payload(h5["results"])
                return ResultData(
                    metadata=metadata,
                    input_data=read_input_data(h5, path),
                    options=read_options(h5),
                    results={"ha": payload},
                    warnings=read_warnings(h5),
                    events=read_events(h5),
                )

        payload, schema_version = read_historical_ha_payload(h5, path)
        warnings: list[str] = []
        if migrate_schema_1_0_thermodynamic_units(payload, schema_version):
            warnings.append(
                "Converted historical HA entropy and heat-capacity arrays to "
                "true native energy units per cell and kelvin."
            )
        return ResultData(
            metadata=ResultMetadata(
                module="ha",
                method="harmonic",
                schema_version=str(schema_version),
            ),
            input_data=None,
            options={},
            results={"ha": payload},
            warnings=warnings,
            events=[],
        )
