# -*- coding: utf-8 -*-

"""YAML reader for quasi-static thermoelastic input data."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import yaml

from quantas.models.structures import CrystalStructure, SymmetryMetadata
from quantas.modules.thermoelasticity.models import (
    ElasticVolumePoint,
    ElasticVolumeSeries,
    ThermoelasticInput,
)


class ThermoelasticInputReader:
    """Read and validate Quantas thermoelastic YAML input data.

    The reader converts schema-versioned YAML mappings into passive
    :class:`~quantas.modules.thermoelasticity.models.ThermoelasticInput`
    contracts and applies the same numerical validation used by the Python
    input generator.
    """

    def load(self, filename: str | Path) -> ThermoelasticInput:
        """Read a quasi-static thermoelastic YAML input file.

        Parameters
        ----------
        filename : str or Path
            YAML input path.

        Returns
        -------
        ThermoelasticInput
            Validated input contract.

        Raises
        ------
        OSError
            If the input file cannot be read.
        ValueError
            If required keys or arrays are invalid.
        """
        path = Path(filename)
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ValueError("thermoelastic input must contain a YAML mapping")
        return _mapping_to_input(raw, source=path)


def read_thermoelastic_input(filename: str | Path) -> ThermoelasticInput:
    """Read a Quantas quasi-static thermoelastic YAML input file.

    Parameters
    ----------
    filename : str or Path
        YAML input path.

    Returns
    -------
    ThermoelasticInput
        Validated input contract.
    """
    return ThermoelasticInputReader().load(filename)


def _mapping_to_input(
    data: dict[str, Any],
    *,
    source: Path,
) -> ThermoelasticInput:
    """Convert a parsed YAML mapping to passive input contracts."""
    schema = _required_mapping(data, "schema")
    if schema.get("name") != "quantas-thermoelastic-input":
        raise ValueError("unsupported thermoelastic input schema name")
    schema_version = str(schema.get("version"))
    if schema_version != "1.0":
        raise ValueError(
            "unsupported thermoelastic input schema version; regenerate the "
            "input with the current thermoelastic inpgen command"
        )
    if str(data.get("method", "")).strip().lower() != "quasistatic":
        raise ValueError("only method: quasistatic is supported")

    reference = _required_mapping(data, "reference")
    structure_data = _required_mapping(reference, "structure")
    structure = CrystalStructure(
        lattice=np.asarray(_required(structure_data, "lattice"), dtype=np.float64),
        fractional_positions=np.asarray(
            _required(structure_data, "fractional_positions"),
            dtype=np.float64,
        ),
        atomic_numbers=np.asarray(
            _required(structure_data, "atomic_numbers"),
            dtype=np.int64,
        ),
        label="thermoelastic reference structure",
        metadata={"source": str(reference.get("source", ""))},
    )
    declared_natom = int(structure_data.get("natom", structure.natoms))
    if declared_natom != structure.natoms:
        raise ValueError("reference structure natom does not match coordinates")

    symmetry_data = _required_mapping(reference, "symmetry")
    symmetry = SymmetryMetadata(
        space_group_number=int(symmetry_data.get("space_group_number", 0)),
        international_symbol=str(symmetry_data.get("international_symbol", "")),
        hall_number=int(symmetry_data.get("hall_number", 0)),
        hall_symbol=str(symmetry_data.get("hall_symbol", "")),
        choice=str(symmetry_data.get("choice", "")),
        point_group=str(symmetry_data.get("point_group", "")),
        symprec=float(symmetry_data.get("symprec", 1.0e-5)),
        angle_tolerance=float(symmetry_data.get("angle_tolerance", -1.0)),
    )

    point_data = data.get("elastic_data")
    if not isinstance(point_data, list) or not point_data:
        raise ValueError("elastic_data must be a non-empty YAML sequence")
    points: list[ElasticVolumePoint] = []
    for index, mapping in enumerate(point_data):
        if not isinstance(mapping, dict):
            raise ValueError(f"elastic_data[{index}] must be a mapping")
        frame = mapping.get("frame")
        if not isinstance(frame, dict):
            raise ValueError(
                f"elastic_data[{index}].frame is required by schema 1.0; "
                "regenerate the input with thermoelastic inpgen"
            )
        required_frame_keys = (
            "status",
            "method",
            "rotation_to_reference",
            "principal_logarithmic_strain",
            "source_lattice",
        )
        missing_frame = [key for key in required_frame_keys if key not in frame]
        if missing_frame:
            missing_text = ", ".join(missing_frame)
            raise ValueError(
                f"elastic_data[{index}].frame is incomplete ({missing_text}); "
                "regenerate the input with thermoelastic inpgen"
            )
        normalized_frame = dict(frame)
        for key in (
            "rotation_to_reference",
            "principal_logarithmic_strain",
            "source_lattice",
        ):
            normalized_frame[key] = np.asarray(normalized_frame[key], dtype=np.float64)
        if np.asarray(normalized_frame["rotation_to_reference"]).shape != (3, 3):
            raise ValueError(
                f"elastic_data[{index}].frame.rotation_to_reference must be 3x3"
            )
        if np.asarray(normalized_frame["source_lattice"]).shape != (3, 3):
            raise ValueError(f"elastic_data[{index}].frame.source_lattice must be 3x3")
        if np.asarray(normalized_frame["principal_logarithmic_strain"]).shape != (3,):
            raise ValueError(
                f"elastic_data[{index}].frame.principal_logarithmic_strain must have length 3"
            )
        point_metadata: dict[str, Any] = {"frame_normalization": normalized_frame}
        points.append(
            ElasticVolumePoint(
                source=str(_required(mapping, "source")),
                pressure=float(_required(mapping, "pressure")),
                stress_pressure=float(
                    mapping.get("stress_pressure", _required(mapping, "pressure"))
                ),
                volume=float(_required(mapping, "volume")),
                density=float(_required(mapping, "density")),
                energy=float(_required(mapping, "energy")),
                prestress_applied=True,
                lattice=np.asarray(_required(mapping, "lattice"), dtype=np.float64),
                stiffness=np.asarray(
                    _required(mapping, "stiffness"),
                    dtype=np.float64,
                ),
                metadata=point_metadata,
            )
        )

    reference_index = int(reference.get("index", 0))
    frame_records = [point.metadata["frame_normalization"] for point in points]
    frame_metadata = {
        "status": "normalized",
        "method": "right_polar_decomposition_corotation",
        "reference_index": reference_index,
        "maximum_removed_rotation_degrees": float(
            max(
                (
                    float(frame.get("removed_rotation_degrees", 0.0))
                    for frame in frame_records
                ),
                default=0.0,
            )
        ),
        "maximum_ordered_atom_displacement_A": float(
            max(
                (
                    float(frame.get("maximum_ordered_atom_displacement_A", 0.0))
                    for frame in frame_records
                ),
                default=0.0,
            )
        ),
        "source": "schema_1.0_frame_metadata",
    }
    series = ElasticVolumeSeries(
        points=tuple(points),
        reference_structure=structure,
        symmetry=symmetry,
        elastic_symmetry=str(symmetry_data.get("elastic_system", "unknown")),
        reference_index=reference_index,
        orientation=str(
            _required_mapping(data, "conventions").get(
                "tensor_orientation",
                "crystal",
            )
        ),
        metadata={
            "interface": str(data.get("interface", "unknown")),
            "schema_version": schema_version,
            "frame_normalization": frame_metadata,
        },
    )
    return ThermoelasticInput(
        jobname=str(data.get("job", "Unknown")),
        elastic_series=series,
        method="quasistatic",
        source=source,
        metadata={
            "schema_name": str(schema.get("name")),
            "schema_version": schema_version,
        },
    )


def _required(mapping: dict[str, Any], key: str) -> Any:
    """Return a required key or raise a readable validation error."""
    if key not in mapping:
        raise ValueError(f"missing required thermoelastic input key: {key}")
    return mapping[key]


def _required_mapping(mapping: dict[str, Any], key: str) -> dict[str, Any]:
    """Return a required child mapping."""
    value = _required(mapping, key)
    if not isinstance(value, dict):
        raise ValueError(f"thermoelastic input key '{key}' must be a mapping")
    return value


__all__ = ["ThermoelasticInputReader", "read_thermoelastic_input"]
