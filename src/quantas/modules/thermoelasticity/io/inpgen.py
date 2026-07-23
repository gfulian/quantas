# -*- coding: utf-8 -*-

"""Generate readable quasi-static thermoelastic YAML input files."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np

from quantas.core.physics.elasticity import detect_elastic_symmetry
from quantas.interfaces.crystal.elasticity import CrystalElasticityReader
from quantas.references import method_citation_keys, render_citation_inline
from quantas.models.structures import CrystalStructure, SymmetryMetadata
from quantas.modules.thermoelasticity.models import (
    ElasticVolumePoint,
    ElasticVolumeSeries,
    ThermoelasticInput,
)
from quantas.modules.thermoelasticity.frames import (
    maximum_ordered_fractional_displacement,
    normalize_elastic_frame,
)


THERMOELASTIC_INPUT_SCHEMA = "1.0"


class ThermoelasticInputCreator:
    """Create a quasi-static thermoelastic input from elastic output files.

    Parameters
    ----------
    interface : str, optional
        Electronic-structure interface.  The first implementation supports
        only ``"crystal"``.
    symprec : float, optional
        Cartesian tolerance in angstrom used for structural symmetry analysis.
    angle_tolerance : float, optional
        Angular tolerance in degrees used for structural symmetry analysis.
    elastic_tolerance : float, optional
        Absolute tolerance in GPa used to detect the elastic crystal system.
    pressure_tolerance : float, optional
        Maximum permitted difference in GPa between CRYSTAL's ``PRESSURE``
        keyword value and the pressure printed for the corrected elastic
        properties.
    structure_correspondence_tolerance : float, optional
        Maximum ordered-atom displacement in angstrom, after wrapping
        fractional differences into the nearest periodic image.  This guards
        against fitting unrelated structural paths while allowing internal
        relaxation along the same phase branch.
    """

    def __init__(
        self,
        *,
        interface: str = "crystal",
        symprec: float = 1.0e-5,
        angle_tolerance: float = -1.0,
        elastic_tolerance: float = 1.0e-3,
        pressure_tolerance: float = 5.0e-2,
        structure_correspondence_tolerance: float = 5.0e-1,
    ) -> None:
        if interface.strip().lower() != "crystal":
            raise ValueError("only the CRYSTAL thermoelastic interface is available")
        self.interface = "crystal"
        self.symprec = float(symprec)
        self.angle_tolerance = float(angle_tolerance)
        self.elastic_tolerance = float(elastic_tolerance)
        self.pressure_tolerance = float(pressure_tolerance)
        self.structure_correspondence_tolerance = float(
            structure_correspondence_tolerance
        )
        if (
            not np.isfinite(self.structure_correspondence_tolerance)
            or self.structure_correspondence_tolerance < 0.0
        ):
            raise ValueError(
                "structure_correspondence_tolerance must be finite and non-negative"
            )

    def create(
        self,
        sources: str | Path | Sequence[str | Path],
        *,
        jobname: str = "Quantas quasi-static thermoelastic input",
        is_list: bool = False,
        reference: int | None = None,
    ) -> ThermoelasticInput:
        """Read CRYSTAL outputs and return normalized thermoelastic input.

        Parameters
        ----------
        sources : str, Path, or sequence
            Output path, list-file path, or explicit sequence of output paths.
        jobname : str, optional
            Human-readable description stored in the YAML file.
        is_list : bool, optional
            Interpret a scalar ``sources`` value as a text file containing one
            output path per line.
        reference : int or None, optional
            Reference index after sorting by increasing volume.  If omitted,
            the point with the smallest absolute pressure is selected.

        Returns
        -------
        ThermoelasticInput
            Validated input contract.

        Raises
        ------
        ValueError
            If a file is invalid, lacks CRYSTAL pre-stress corrections, or is
            inconsistent with the rest of the series.
        """
        files = self._normalize_sources(sources, is_list=is_list)
        if not files:
            raise ValueError("no CRYSTAL elastic output files were provided")

        records: list[
            tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]
        ] = []
        for path in files:
            reader = CrystalElasticityReader(
                path,
                symprec=self.symprec,
                angle_tolerance=self.angle_tolerance,
            )
            if not reader.completed:
                raise ValueError(reader.error or f"unable to read {path}")
            if not reader.prestress_applied:
                raise ValueError(
                    f"{path}: CRYSTAL's PRESSURE keyword is required so that "
                    "hydrostatic pre-stress terms are included"
                )
            if not np.isclose(
                reader.pressure_keyword_value,
                reader.pressure,
                rtol=0.0,
                atol=self.pressure_tolerance,
            ):
                raise ValueError(
                    f"{path}: PRESSURE keyword value {reader.pressure_keyword_value:.6g} "
                    f"GPa differs from elastic pressure {reader.pressure:.6g} GPa"
                )
            structure = reader.structure
            symmetry = reader.symmetry
            if structure is None or symmetry is None:
                raise ValueError(f"{path}: final structure or symmetry is unavailable")
            elastic_symmetry = detect_elastic_symmetry(
                reader.stiffness,
                tolerance=self.elastic_tolerance,
            )
            point = ElasticVolumePoint(
                source=path.name,
                pressure=reader.pressure,
                stress_pressure=reader.stress_pressure,
                volume=reader.volume,
                density=reader.density,
                energy=reader.energy,
                stiffness=reader.stiffness,
                lattice=structure.lattice,
                prestress_applied=reader.prestress_applied,
            )
            records.append((point, structure, symmetry, elastic_symmetry))

        records.sort(key=lambda item: item[0].volume)
        self._validate_records(records)
        points = tuple(item[0] for item in records)
        if reference is None:
            reference_index = int(
                np.argmin(np.abs([point.pressure for point in points]))
            )
        else:
            reference_index = int(reference)
            if reference_index < 0 or reference_index >= len(points):
                raise ValueError("reference index is outside the sorted elastic series")

        reference_structure = records[reference_index][1]
        reference_symmetry = records[reference_index][2]
        elastic_symmetry = records[reference_index][3]
        records, frame_metadata = self._normalize_frames(
            records,
            reference_index=reference_index,
        )
        points = tuple(item[0] for item in records)
        series = ElasticVolumeSeries(
            points=points,
            reference_structure=reference_structure,
            symmetry=reference_symmetry,
            elastic_symmetry=elastic_symmetry,
            reference_index=reference_index,
            orientation="crystal",
            metadata={
                "interface": self.interface,
                "ordering": "volume-ascending",
                "pressure_keyword_required": True,
                "symprec": self.symprec,
                "angle_tolerance": self.angle_tolerance,
                "elastic_tolerance_gpa": self.elastic_tolerance,
                "frame_normalization": frame_metadata,
                "structure_correspondence_tolerance_A": (
                    self.structure_correspondence_tolerance
                ),
            },
        )
        return ThermoelasticInput(
            jobname=str(jobname),
            elastic_series=series,
            method="quasistatic",
            metadata={
                "schema_name": "quantas-thermoelastic-input",
                "schema_version": THERMOELASTIC_INPUT_SCHEMA,
            },
        )

    def write(
        self,
        input_data: ThermoelasticInput,
        outfile: str | Path,
    ) -> Path:
        """Write a readable thermoelastic YAML file.

        Parameters
        ----------
        input_data : ThermoelasticInput
            Input contract to serialize.
        outfile : str or Path
            Destination YAML path.

        Returns
        -------
        Path
            Written path.
        """
        path = Path(outfile)
        path.write_text(format_thermoelastic_yaml(input_data), encoding="utf-8")
        return path

    @staticmethod
    def _normalize_sources(
        sources: str | Path | Sequence[str | Path],
        *,
        is_list: bool,
    ) -> list[Path]:
        """Return explicit source paths, resolving list entries locally."""
        if is_list:
            if isinstance(sources, Sequence) and not isinstance(sources, (str, Path)):
                raise ValueError("is_list requires one list-file path")
            list_path = Path(sources)
            files: list[Path] = []
            for line in list_path.read_text(encoding="utf-8").splitlines():
                value = line.strip()
                if not value or value.startswith("#"):
                    continue
                item = Path(value)
                if not item.is_absolute():
                    item = list_path.parent / item
                files.append(item)
        elif isinstance(sources, (str, Path)):
            files = [Path(sources)]
        else:
            files = [Path(value) for value in sources]

        resolved: list[Path] = []
        seen: set[Path] = set()
        for path in files:
            normalized = path.expanduser().resolve()
            if normalized in seen:
                continue
            if not normalized.is_file():
                raise ValueError(f"elastic output file does not exist: {path}")
            seen.add(normalized)
            resolved.append(normalized)
        return resolved

    @staticmethod
    def _validate_records(
        records: list[
            tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]
        ],
    ) -> None:
        """Validate structural, crystallographic, and elastic consistency."""
        reference_point, reference_structure, reference_symmetry, elastic_symmetry = (
            records[0]
        )
        del reference_point
        for point, structure, symmetry, current_elastic_symmetry in records[1:]:
            if not np.array_equal(
                structure.atomic_numbers,
                reference_structure.atomic_numbers,
            ):
                raise ValueError(
                    f"{point.source}: primitive atomic species/order differs across inputs"
                )
            if symmetry.space_group_number != reference_symmetry.space_group_number:
                raise ValueError(
                    f"{point.source}: space group {symmetry.space_group_number} differs "
                    f"from reference {reference_symmetry.space_group_number}"
                )
            if (
                reference_symmetry.hall_number > 0
                and symmetry.hall_number > 0
                and symmetry.hall_number != reference_symmetry.hall_number
            ):
                raise ValueError(
                    f"{point.source}: Hall number {symmetry.hall_number} differs "
                    f"from reference {reference_symmetry.hall_number}"
                )
            if (
                reference_symmetry.choice
                and symmetry.choice
                and symmetry.choice != reference_symmetry.choice
            ):
                raise ValueError(
                    f"{point.source}: crystallographic setting choice "
                    f"{symmetry.choice!r} differs from reference "
                    f"{reference_symmetry.choice!r}"
                )
            if current_elastic_symmetry != elastic_symmetry:
                raise ValueError(
                    f"{point.source}: elastic symmetry {current_elastic_symmetry} differs "
                    f"from reference {elastic_symmetry}"
                )
        volumes = np.asarray([item[0].volume for item in records], dtype=np.float64)
        if np.any(np.diff(volumes) <= 1.0e-10):
            raise ValueError("elastic outputs contain duplicate or unresolved volumes")

    def _normalize_frames(
        self,
        records: list[
            tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]
        ],
        *,
        reference_index: int,
    ) -> tuple[
        list[tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]],
        dict[str, Any],
    ]:
        """Co-rotate every tensor into the selected reference Cartesian frame.

        Parameters
        ----------
        records : list
            Parsed and volume-sorted elastic records.
        reference_index : int
            Index defining the common Cartesian frame.

        Returns
        -------
        tuple
            Normalized records and series-level diagnostics.

        Raises
        ------
        ValueError
            If atom correspondence is inconsistent or a frame transformation
            is improper.
        """
        reference_structure = records[reference_index][1]
        normalized: list[
            tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]
        ] = []
        angles: list[float] = []
        displacements: list[float] = []
        for point, structure, symmetry, elastic_symmetry in records:
            displacement = maximum_ordered_fractional_displacement(
                reference_structure.fractional_positions,
                structure.fractional_positions,
                reference_structure.lattice,
            )
            if displacement > self.structure_correspondence_tolerance:
                raise ValueError(
                    f"{point.source}: maximum ordered-atom displacement "
                    f"{displacement:.6g} A exceeds the configured structural "
                    f"path tolerance {self.structure_correspondence_tolerance:.6g} A"
                )
            frame = normalize_elastic_frame(
                point.lattice,
                point.stiffness,
                reference_structure.lattice,
            )
            metadata = {
                **point.metadata,
                "frame_normalization": {
                    **frame.metadata,
                    "status": "normalized",
                    "rotation_to_reference": frame.rotation_to_reference.copy(),
                    "removed_rotation_degrees": frame.removed_rotation_degrees,
                    "principal_logarithmic_strain": (
                        frame.principal_logarithmic_strain.copy()
                    ),
                    "source_lattice": point.lattice.copy(),
                    "maximum_ordered_atom_displacement_A": displacement,
                },
            }
            normalized_point = ElasticVolumePoint(
                source=point.source,
                pressure=point.pressure,
                stress_pressure=point.stress_pressure,
                volume=point.volume,
                density=point.density,
                energy=point.energy,
                stiffness=frame.stiffness,
                lattice=frame.lattice,
                prestress_applied=point.prestress_applied,
                metadata=metadata,
            )
            normalized.append((normalized_point, structure, symmetry, elastic_symmetry))
            angles.append(frame.removed_rotation_degrees)
            displacements.append(displacement)
        return normalized, {
            "status": "normalized",
            "method": "right_polar_decomposition_corotation",
            "reference_index": int(reference_index),
            "maximum_removed_rotation_degrees": float(max(angles, default=0.0)),
            "maximum_ordered_atom_displacement_A": float(
                max(displacements, default=0.0)
            ),
            "reference": "; ".join(
                render_citation_inline(key)
                for key in method_citation_keys("wallace_stress_strain")
            ),
            "citation_keys": list(method_citation_keys("wallace_stress_strain")),
        }


def create_thermoelastic_input(
    sources: str | Path | Sequence[str | Path],
    outfile: str | Path,
    *,
    interface: str = "crystal",
    is_list: bool = False,
    jobname: str = "Quantas quasi-static thermoelastic input",
    reference: int | None = None,
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
    elastic_tolerance: float = 1.0e-3,
    pressure_tolerance: float = 5.0e-2,
    structure_correspondence_tolerance: float = 5.0e-1,
) -> Path:
    """Create a quasi-static thermoelastic YAML input through the Python API.

    Parameters
    ----------
    sources : str, Path, or sequence
        CRYSTAL output path, list-file path, or explicit output sequence.
    outfile : str or Path
        Destination YAML file.
    interface : str, optional
        Electronic-structure interface.  Currently ``"crystal"`` only.
    is_list : bool, optional
        Interpret a scalar source as a list file.
    jobname : str, optional
        Input description.
    reference : int or None, optional
        Reference index after volume sorting.
    symprec : float, optional
        spglib Cartesian tolerance in angstrom.
    angle_tolerance : float, optional
        spglib angular tolerance in degrees.
    elastic_tolerance : float, optional
        Elastic symmetry tolerance in GPa.
    pressure_tolerance : float, optional
        Maximum permitted difference, in GPa, between the CRYSTAL ``PRESSURE``
        keyword and the pressure reported for the corrected elastic tensor.
    structure_correspondence_tolerance : float, optional
        Maximum ordered-atom displacement in angstrom along the structural
        path.

    Returns
    -------
    Path
        Written YAML path.
    """
    creator = ThermoelasticInputCreator(
        interface=interface,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        elastic_tolerance=elastic_tolerance,
        pressure_tolerance=pressure_tolerance,
        structure_correspondence_tolerance=structure_correspondence_tolerance,
    )
    input_data = creator.create(
        sources,
        jobname=jobname,
        is_list=is_list,
        reference=reference,
    )
    return creator.write(input_data, outfile)


def format_thermoelastic_yaml(input_data: ThermoelasticInput) -> str:
    """Serialize thermoelastic input with compact row-wise vectors and matrices.

    Parameters
    ----------
    input_data : ThermoelasticInput
        Validated thermoelastic input contract.

    Returns
    -------
    str
        Readable YAML text.
    """
    series = input_data.elastic_series
    structure = series.reference_structure
    symmetry = series.symmetry
    lines: list[str] = [
        "schema:",
        "  name: quantas-thermoelastic-input",
        f"  version: '{THERMOELASTIC_INPUT_SCHEMA}'",
        f"job: {_yaml_string(input_data.jobname)}",
        f"method: {input_data.method}",
        "interface: crystal",
        "conventions:",
        "  strain: eulerian-finite-strain",
        "  stiffness: wallace-hydrostatic-stress-strain",
        "  prestress: applied-by-crystal-pressure-keyword",
        "  tensor_orientation: crystal",
        "  voigt_order: [ 11, 22, 33, 23, 13, 12 ]",
        "units:",
        "  pressure: GPa",
        "  stiffness: GPa",
        "  volume: angstrom^3",
        "  density: kg/m^3",
        "  energy: hartree",
        "reference:",
        f"  index: {series.reference_index}",
        f"  source: {_yaml_string(str(series.points[series.reference_index].source))}",
        "  structure:",
        f"    natom: {structure.natoms}",
        "    lattice:",
    ]
    lines.extend(_matrix_lines(structure.lattice, indent="    ", precision=12))
    atomic_numbers = ", ".join(str(int(value)) for value in structure.atomic_numbers)
    lines.append(f"    atomic_numbers: [ {atomic_numbers} ]")
    lines.append("    fractional_positions:")
    lines.extend(
        _matrix_lines(
            structure.fractional_positions,
            indent="    ",
            precision=12,
        )
    )
    lines.extend(
        [
            "  symmetry:",
            f"    space_group_number: {symmetry.space_group_number}",
            f"    international_symbol: {_yaml_string(symmetry.international_symbol)}",
            f"    hall_number: {symmetry.hall_number}",
            f"    hall_symbol: {_yaml_string(symmetry.hall_symbol)}",
            f"    choice: {_yaml_string(symmetry.choice)}",
            f"    point_group: {_yaml_string(symmetry.point_group)}",
            f"    elastic_system: {series.elastic_symmetry}",
            f"    symprec: {symmetry.symprec:.6E}",
            f"    angle_tolerance: {symmetry.angle_tolerance:.6f}",
            "elastic_data:",
        ]
    )
    for point in series.points:
        lines.extend(
            [
                f"- source: {_yaml_string(str(point.source))}",
                f"  pressure: {point.pressure: .8f}",
                f"  stress_pressure: {point.stress_pressure: .8f}",
                f"  volume: {point.volume: .12f}",
                f"  density: {point.density: .8f}",
                f"  energy: {point.energy: .12E}",
                "  lattice:",
            ]
        )
        lines.extend(_matrix_lines(point.lattice, indent="  ", precision=12))
        lines.append("  stiffness:")
        lines.extend(_matrix_lines(point.stiffness, indent="  ", precision=8))
        frame = point.metadata.get("frame_normalization")
        if isinstance(frame, dict):
            lines.extend(
                [
                    "  frame:",
                    f"    status: {_yaml_string(str(frame.get('status', 'unknown')))}",
                    f"    method: {_yaml_string(str(frame.get('method', 'unknown')))}",
                    "    removed_rotation_degrees: "
                    f"{float(frame.get('removed_rotation_degrees', 0.0)): .12E}",
                    "    maximum_ordered_atom_displacement_A: "
                    f"{float(frame.get('maximum_ordered_atom_displacement_A', 0.0)): .12E}",
                    "    rotation_to_reference:",
                ]
            )
            lines.extend(
                _matrix_lines(
                    frame.get("rotation_to_reference", np.eye(3)),
                    indent="    ",
                    precision=12,
                )
            )
            principal = np.asarray(
                frame.get("principal_logarithmic_strain", np.zeros(3)),
                dtype=np.float64,
            )
            principal_text = ", ".join(f"{value: .12E}" for value in principal)
            lines.append(f"    principal_logarithmic_strain: [ {principal_text} ]")
            lines.append("    source_lattice:")
            lines.extend(
                _matrix_lines(
                    frame.get("source_lattice", point.lattice),
                    indent="    ",
                    precision=12,
                )
            )
    return "\n".join(lines) + "\n"


def _matrix_lines(
    matrix: Any,
    *,
    indent: str,
    precision: int,
) -> list[str]:
    """Format a two-dimensional array with one compact row per YAML line."""
    array = np.asarray(matrix, dtype=np.float64)
    rows: list[str] = []
    threshold = 0.5 * 10.0 ** (-precision)
    for row in array:
        cleaned = [
            0.0 if abs(float(value)) < threshold else float(value) for value in row
        ]
        values = ", ".join(f"{value: .{precision}f}" for value in cleaned)
        rows.append(f"{indent}- [ {values} ]")
    return rows


def _yaml_string(value: str) -> str:
    """Return a single-quoted YAML scalar with escaped apostrophes."""
    return "'" + str(value).replace("'", "''") + "'"


__all__ = [
    "THERMOELASTIC_INPUT_SCHEMA",
    "ThermoelasticInputCreator",
    "create_thermoelastic_input",
    "format_thermoelastic_yaml",
]
