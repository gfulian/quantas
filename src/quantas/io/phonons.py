# -*- coding: utf-8 -*-

"""Reader for normalized Quantas volume-dependent phonon YAML input files.

The reader produces :class:`~quantas.models.phonons.PhononInputData`, a neutral
container that can be consumed by harmonic, quasi-harmonic, phonon-analysis,
and future thermoelastic workflows.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import yaml

from quantas.models.reader import BasicReader
from quantas.models.phonons import PhononInputData
from quantas.models.structures import (
    CellNormalization,
    CrystalStructure,
    StructureReconstructionDiagnostics,
    StructureVolumeSeries,
    SymmetryMetadata,
)


class PhononInputFileReader(BasicReader):
    """
    Reader for Quantas volume-dependent phonon YAML input files.

    Parameters
    ----------
    filename : str or Path or None, optional
        Path to the YAML input file. If provided, the file is read during
        initialization.

    Attributes
    ----------
    completed : bool
        Flag set to ``True`` when the input file has been read successfully.
    error : str or None
        Error message generated during parsing or validation, if any.
    """

    def __init__(self, filename: str | Path | None = None) -> None:
        """Initialize the reader and optionally load an input file.

        Parameters
        ----------
        filename : str or Path or None, optional
            Path to a YAML input file read during initialization.
        """
        super().__init__()
        self._data: dict[str, Any] | None = None
        self._source: Path | None = None

        if filename is not None:
            self.load(filename)

    def load(self, filename: str | Path) -> None:
        """
        Read a Quantas phonon YAML input file.

        Parameters
        ----------
        filename : str or Path
            Path to the YAML input file.
        """
        filename = Path(filename)
        self.completed = False
        self.error = None
        self._data = None
        self._source = filename

        try:
            with filename.open("r", encoding="utf-8") as stream:
                self._data = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            self.error = self._format_yaml_error(exc)
            return
        except UnicodeDecodeError:
            self.error = "The provided input is not a YAML file."
            return
        except OSError as exc:
            self.error = f"Unable to read input file: {exc}"
            return

        self._check()

    def _check(self) -> None:
        """
        Check that the minimum phonon input fields are present.
        """
        if self._data is None:
            return
        if not isinstance(self._data, dict):
            self.error = "The YAML input should contain a mapping."
            return
        if "volume" not in self._data:
            self.error = "No volume values in input file"
            return
        if "energy" not in self._data:
            self.error = "No energy values in input file"
            return
        if "phonon" not in self._data:
            self.error = "No phonon data in input file"
            return

        self.completed = True

    @staticmethod
    def _format_yaml_error(exc: yaml.YAMLError) -> str:
        """
        Format a YAML parser error for user-facing diagnostics.

        Parameters
        ----------
        exc : yaml.YAMLError
            Error raised by PyYAML.

        Returns
        -------
        str
            Human-readable error message.
        """
        if hasattr(exc, "problem_mark"):
            error = "Error in parsing YAML file: \n"
            error += "  - " + str(exc.problem_mark) + "\n"
            error += "  - " + str(exc.problem) + "\n"
            error += "\nPlease, correct data and retry."
            return error

        return (
            "Error in parsing YAML file: unknown problem."
            "\n\nPlease, correct data and retry."
        )

    @staticmethod
    def _get_float(item: Any) -> np.float64 | np.ndarray:
        """
        Convert a YAML scalar or sequence to float data.

        Parameters
        ----------
        item : object
            YAML item to convert. Historical Quantas inputs often store
            multiple values as a one-item list containing a whitespace-separated
            string.

        Returns
        -------
        np.float64 or ndarray
            Converted scalar or one-dimensional array.

        Raises
        ------
        ValueError
            If the value cannot be converted to float.
        """
        try:
            return np.float64(item)
        except (TypeError, ValueError):
            pass

        if isinstance(item, (list, tuple)):
            if len(item) == 1 and isinstance(item[0], str):
                return np.array(item[0].split(), dtype=np.float64)
            return np.array(item, dtype=np.float64)

        if isinstance(item, str):
            values = item.split()
            if len(values) == 1:
                return np.float64(values[0])
            return np.array(values, dtype=np.float64)

        return np.array(item, dtype=np.float64)

    @staticmethod
    def _get_int(item: Any) -> np.int64 | np.ndarray:
        """
        Convert a YAML scalar or sequence to integer data.

        Parameters
        ----------
        item : object
            YAML item to convert.

        Returns
        -------
        np.int64 or ndarray
            Converted scalar or one-dimensional array.

        Raises
        ------
        ValueError
            If the value cannot be converted to integer.
        """
        try:
            return np.int64(item)
        except (TypeError, ValueError):
            pass

        if isinstance(item, (list, tuple)):
            if len(item) == 1 and isinstance(item[0], str):
                return np.array(item[0].split(), dtype=np.int64)
            return np.array(item, dtype=np.int64)

        if isinstance(item, str):
            values = item.split()
            if len(values) == 1:
                return np.int64(values[0])
            return np.array(values, dtype=np.int64)

        return np.array(item, dtype=np.int64)

    @staticmethod
    def _get_string(item: Any) -> str:
        """
        Convert a YAML item to a string.

        Parameters
        ----------
        item : object
            YAML item to convert.

        Returns
        -------
        str
            Converted string with trailing whitespace removed.
        """
        return str(item).rstrip()

    @property
    def data(self) -> dict[str, Any] | None:
        """
        Return the raw YAML data.

        Returns
        -------
        dict or None
            Parsed YAML mapping, or ``None`` if no data were loaded.
        """
        return self._data

    def _require_data(self) -> dict[str, Any]:
        """Return parsed data or raise when no file has been loaded."""
        if self._data is None:
            raise RuntimeError("No phonon input data have been loaded.")
        return self._data

    @property
    def jobname(self) -> str:
        """
        Return the job name.

        Returns
        -------
        str
            Job name stored in the input file, or ``"Unknown"`` if absent.
        """
        if self._data is not None and "job" in self._data:
            return self._get_string(self._data["job"])
        return "Unknown"

    @property
    def natoms(self) -> np.int64:
        """
        Return the number of atoms in the thermodynamic normalization cell.

        Returns
        -------
        np.int64
            Number of atoms.
        """
        return np.int64(self._get_int(self._require_data()["natom"]))

    @property
    def formula_units(self) -> np.int64:
        """
        Return the number of formula units in the thermodynamic normalization cell.

        Returns
        -------
        np.int64
            Number of formula units. Inputs without this field default to one.

        Raises
        ------
        ValueError
            If the stored value is not positive.
        """
        value: Any = 1
        if self._data is not None:
            for key in ("formula_units", "z", "Z"):
                if key in self._data:
                    value = self._data[key]
                    break
        formula_units = np.int64(self._get_int(value))
        if formula_units <= 0:
            raise ValueError("formula_units must be positive")
        return formula_units

    @property
    def kpoints(self) -> int:
        """
        Return the number of sampled k-points from the supercell matrix.

        Returns
        -------
        int
            Rounded determinant of the supercell matrix.
        """
        return int(np.around(np.linalg.det(self._require_data()["supercell"]), 0))

    @property
    def qpoints(self) -> np.int64:
        """
        Return the number of q-points.

        Returns
        -------
        np.int64
            Number of q-points.
        """
        return np.int64(self._get_int(self._require_data()["qpoints"]))

    @property
    def nvol(self) -> int:
        """
        Return the number of volumes in the input file.

        Returns
        -------
        int
            Number of volumes.
        """
        volume = self.volume
        if isinstance(volume, np.ndarray):
            return int(volume.shape[0])
        return 1

    @property
    def volume(self) -> np.float64 | np.ndarray:
        """
        Return unit-cell volumes.

        Returns
        -------
        np.float64 or ndarray
            Volume value or values as stored in the YAML input.
        """
        return self._get_float(self._require_data()["volume"])

    @property
    def energy(self) -> np.float64 | np.ndarray:
        """
        Return static energies.

        Returns
        -------
        np.float64 or ndarray
            Static energy value or values as stored in the YAML input.
        """
        return self._get_float(self._require_data()["energy"])

    @property
    def frequencies(self) -> np.ndarray:
        """
        Return phonon frequencies.

        Returns
        -------
        ndarray
            Phonon frequency matrix with shape ``(qpoints, natoms * 3, nvol)``.
        """
        matrix = np.zeros(
            (int(self.qpoints), int(self.natoms) * 3, self.nvol),
            dtype=np.float64,
        )
        for i, qpoint in enumerate(self._require_data()["phonon"]):
            for j, frequency in enumerate(qpoint["band"]):
                matrix[i, j] = np.array(
                    self._get_float(frequency["frequency"]),
                    dtype=np.float64,
                )
        return matrix

    @property
    def weights(self) -> np.ndarray:
        """
        Return q-point weights.

        Returns
        -------
        ndarray
            Q-point weights with shape ``(qpoints,)``.
        """
        weights = np.zeros(int(self.qpoints), dtype=np.float64)
        for i, qpoint in enumerate(self._require_data()["phonon"]):
            weights[i] = self._get_float(qpoint["weight"])
        return weights

    @property
    def qcoords(self) -> np.ndarray | None:
        """Return fractional q-point coordinates when available.

        Returns
        -------
        ndarray or None
            Fractional coordinates with shape ``(qpoints, 3)``, or ``None``
            when the input does not provide q-point positions.
        """
        if self._data is None:
            return None
        phonons = self._data.get("phonon", [])
        if not phonons or any(
            "q-position" not in qpoint or qpoint.get("q-position") is None
            for qpoint in phonons
        ):
            return None
        return np.asarray(
            [qpoint["q-position"] for qpoint in phonons],
            dtype=np.float64,
        )

    @property
    def structure(self) -> StructureVolumeSeries | None:
        """Return the optional compact primitive structural path.

        Returns
        -------
        StructureVolumeSeries or None
            Structural lattices, coordinates, symmetry, and normalization, or
            ``None`` for historical inputs without the new structure block.
        """
        if self._data is None or "structure" not in self._data:
            return None
        return _structure_series_from_mapping(self._data["structure"])

    @property
    def total_q_points(self) -> np.float64:
        """
        Return the sum of q-point weights.

        Returns
        -------
        np.float64
            Sum of the q-point weights.
        """
        return np.float64(self.weights.sum())

    def to_input(self, source: str | Path | None = None) -> PhononInputData:
        """
        Convert parsed YAML data to a normalized phonon input container.

        Parameters
        ----------
        source : str or Path or None, optional
            Source path to store in the returned input object. If omitted, the
            path passed to :meth:`load` is used when available.

        Returns
        -------
        PhononInputData
            Normalized volume-dependent phonon input data.

        Raises
        ------
        ValueError
            If the reader has not successfully loaded a valid input file.
        """
        if not self.completed:
            raise ValueError(self.error or "Phonon input data have not been loaded")

        supercell = None
        if self._data is not None and "supercell" in self._data:
            supercell = np.asarray(self._require_data()["supercell"], dtype=np.float64)

        return PhononInputData(
            jobname=self.jobname,
            natoms=int(self.natoms),
            formula_units=int(self.formula_units),
            supercell=supercell,
            qpoints=int(self.qpoints),
            volume=np.atleast_1d(np.asarray(self.volume, dtype=np.float64)),
            energy=np.atleast_1d(np.asarray(self.energy, dtype=np.float64)),
            frequencies=np.asarray(self.frequencies, dtype=np.float64),
            weights=np.asarray(self.weights, dtype=np.float64),
            qcoords=(
                None
                if self.qcoords is None
                else np.asarray(self.qcoords, dtype=np.float64)
            ),
            structure=self.structure,
            source=source if source is not None else self._source,
            metadata={
                "format": "quantas-phonon-yaml",
                "formula_units_per_cell": int(self.formula_units),
                "natoms_per_formula_unit": (
                    float(self.natoms) / float(self.formula_units)
                ),
                "q_position_source": (
                    None if self._data is None else self._data.get("q_position_source")
                ),
                "q_position_convention": (
                    None
                    if self._data is None
                    else self._data.get("q_position_convention")
                ),
            },
        )


def _structure_series_from_mapping(mapping: Any) -> StructureVolumeSeries:
    """Convert a YAML structure mapping to passive Quantas models.

    Parameters
    ----------
    mapping : object
        Parsed YAML structure block.

    Returns
    -------
    StructureVolumeSeries
        Normalized structural volume series.

    Raises
    ------
    ValueError
        If required structural fields are absent or inconsistent.
    """
    if not isinstance(mapping, dict):
        raise ValueError("structure must be a YAML mapping")
    reference_raw = mapping.get("reference")
    series_raw = mapping.get("volume_series")
    normalization_raw = mapping.get("normalization")
    if not isinstance(reference_raw, dict):
        raise ValueError("structure.reference must be a mapping")
    if not isinstance(series_raw, dict):
        raise ValueError("structure.volume_series must be a mapping")
    if not isinstance(normalization_raw, dict):
        raise ValueError("structure.normalization must be a mapping")

    atomic_numbers = np.asarray(
        mapping.get("atomic_numbers", reference_raw.get("atomic_numbers", [])),
        dtype=np.int64,
    )
    reference = CrystalStructure(
        lattice=np.asarray(reference_raw["lattice"], dtype=np.float64),
        fractional_positions=np.asarray(
            reference_raw["fractional_positions"],
            dtype=np.float64,
        ),
        atomic_numbers=atomic_numbers,
        label=str(reference_raw.get("label", "reference primitive")),
        metadata=dict(reference_raw.get("metadata", {})),
    )
    expansion = np.asarray(
        normalization_raw.get("expansion_matrix", np.eye(3)),
        dtype=np.int64,
    )
    normalization = CellNormalization(
        basis=str(normalization_raw.get("basis", "primitive")),
        source_basis=str(normalization_raw.get("source_basis", "unknown")),
        expansion_matrix=expansion,
        repetitions=int(
            normalization_raw.get(
                "repetitions",
                round(abs(float(np.linalg.det(expansion)))),
            )
        ),
        source_atoms=int(normalization_raw.get("source_atoms", reference.natoms)),
        normalized_atoms=int(
            normalization_raw.get("normalized_atoms", reference.natoms)
        ),
    )

    symmetry = None
    symmetry_raw = mapping.get("symmetry")
    if isinstance(symmetry_raw, dict):
        symmetry = SymmetryMetadata(
            space_group_number=int(symmetry_raw.get("space_group_number", 0)),
            international_symbol=str(symmetry_raw.get("international_symbol", "")),
            hall_number=int(symmetry_raw.get("hall_number", 0)),
            hall_symbol=str(symmetry_raw.get("hall_symbol", "")),
            choice=str(symmetry_raw.get("choice", "")),
            point_group=str(symmetry_raw.get("point_group", "")),
            symprec=float(symmetry_raw.get("symprec", 1.0e-5)),
            angle_tolerance=float(symmetry_raw.get("angle_tolerance", -1.0)),
            equivalent_atoms=(
                None
                if "equivalent_atoms" not in symmetry_raw
                else np.asarray(symmetry_raw["equivalent_atoms"], dtype=np.int64)
            ),
            transformation_matrix=(
                None
                if "transformation_matrix" not in symmetry_raw
                else np.asarray(
                    symmetry_raw["transformation_matrix"],
                    dtype=np.float64,
                )
            ),
            origin_shift=(
                None
                if "origin_shift" not in symmetry_raw
                else np.asarray(symmetry_raw["origin_shift"], dtype=np.float64)
            ),
        )

    diagnostics: list[StructureReconstructionDiagnostics] = []
    for item in mapping.get("reconstruction", []):
        diagnostics.append(
            StructureReconstructionDiagnostics(
                status=str(item.get("status", "unknown")),
                source_atoms=int(item.get("source_atoms", 0)),
                reconstructed_atoms=int(item.get("reconstructed_atoms", 0)),
                expected_repetitions=int(item.get("expected_repetitions", 1)),
                minimum_replica_count=int(item.get("minimum_replica_count", 0)),
                maximum_replica_count=int(item.get("maximum_replica_count", 0)),
                maximum_translation_residual=float(
                    item.get("maximum_translation_residual_angstrom", 0.0)
                ),
                rms_translation_residual=float(
                    item.get("rms_translation_residual_angstrom", 0.0)
                ),
                message=str(item.get("message", "")),
            )
        )

    transformations = mapping.get("transformations", {})
    primitive_to_crystallographic = None
    if (
        isinstance(transformations, dict)
        and "primitive_to_crystallographic" in transformations
    ):
        primitive_to_crystallographic = np.asarray(
            transformations["primitive_to_crystallographic"],
            dtype=np.float64,
        )

    return StructureVolumeSeries(
        reference=reference,
        lattices=np.asarray(series_raw["lattice"], dtype=np.float64),
        fractional_positions=np.asarray(
            series_raw["fractional_positions"],
            dtype=np.float64,
        ),
        volumes=np.asarray(series_raw["volume"], dtype=np.float64),
        normalization=normalization,
        symmetry=symmetry,
        primitive_to_crystallographic=primitive_to_crystallographic,
        diagnostics=tuple(diagnostics),
        orientation=str(mapping.get("orientation", "crystal")),
        reference_index=int(mapping.get("reference_index", 0)),
        metadata=dict(mapping.get("metadata", {})),
    )
