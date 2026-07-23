# -*- coding: utf-8 -*-

"""Passive structural data contracts shared by Quantas workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(slots=True)
class CrystalStructure:
    """Compact periodic crystal structure.

    Parameters
    ----------
    lattice : array_like
        Direct lattice vectors in Cartesian coordinates, stored by rows with
        shape ``(3, 3)`` and expressed in angstrom.
    fractional_positions : array_like
        Fractional atomic coordinates with shape ``(natoms, 3)``.
    atomic_numbers : array_like
        Atomic numbers with shape ``(natoms,)``.
    label : str, optional
        Human-readable description of the structure.
    metadata : dict, optional
        Additional parser or provenance metadata.

    Raises
    ------
    ValueError
        If array shapes are inconsistent or the lattice volume is zero.
    """

    lattice: FloatArray
    fractional_positions: FloatArray
    atomic_numbers: IntArray
    label: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays to Quantas structural dtypes and validate shapes."""
        self.lattice = np.asarray(self.lattice, dtype=np.float64)
        self.fractional_positions = np.asarray(
            self.fractional_positions,
            dtype=np.float64,
        )
        self.atomic_numbers = np.asarray(self.atomic_numbers, dtype=np.int64)
        if self.lattice.shape != (3, 3):
            raise ValueError("lattice must have shape (3, 3)")
        if (
            self.fractional_positions.ndim != 2
            or self.fractional_positions.shape[1] != 3
        ):
            raise ValueError("fractional_positions must have shape (natoms, 3)")
        if self.atomic_numbers.ndim != 1:
            raise ValueError("atomic_numbers must be one-dimensional")
        if self.fractional_positions.shape[0] != self.atomic_numbers.shape[0]:
            raise ValueError(
                "atomic numbers and positions contain different atom counts"
            )
        if abs(float(np.linalg.det(self.lattice))) <= np.finfo(np.float64).eps:
            raise ValueError("lattice volume must be non-zero")

    @property
    def natoms(self) -> int:
        """Return the number of atoms in the structure.

        Returns
        -------
        int
            Number of stored atoms.
        """
        return int(self.atomic_numbers.shape[0])

    @property
    def volume(self) -> float:
        """Return the positive cell volume in cubic angstrom.

        Returns
        -------
        float
            Absolute determinant of the direct lattice matrix.
        """
        return abs(float(np.linalg.det(self.lattice)))

    def spglib_cell(
        self,
    ) -> tuple[list[list[float]], list[list[float]], list[int]]:
        """Return the tuple representation consumed by spglib.

        Returns
        -------
        tuple
            ``(lattice, fractional_positions, atomic_numbers)`` with lattice
            vectors stored by rows.
        """
        return (
            self.lattice.tolist(),
            self.fractional_positions.tolist(),
            self.atomic_numbers.tolist(),
        )

    def as_dict(self) -> dict[str, Any]:
        """Return a recursively serializable mapping.

        Returns
        -------
        dict
            Structural arrays and metadata suitable for YAML or HDF5 output.
        """
        return {
            "label": self.label,
            "lattice": self.lattice.copy(),
            "fractional_positions": self.fractional_positions.copy(),
            "atomic_numbers": self.atomic_numbers.copy(),
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class CellNormalization:
    """Describe the cell basis used to normalize phonon thermodynamics.

    Parameters
    ----------
    basis : str
        Basis to which energy, volume, atom count, and phonon modes refer.
        Typical values are ``"primitive"`` and ``"phonon_supercell"``.
    source_basis : str
        Cell representation printed by the source electronic-structure code.
    expansion_matrix : array_like
        Integer matrix that expands the primitive cell into the source
        supercell.
    repetitions : int
        Number of primitive-cell repetitions, normally the absolute
        determinant of ``expansion_matrix``.
    source_atoms : int
        Atom count in the source cell.
    normalized_atoms : int
        Atom count in the thermodynamic normalization cell.

    Raises
    ------
    ValueError
        If the expansion matrix or repetition count is invalid.
    """

    basis: str
    source_basis: str
    expansion_matrix: IntArray
    repetitions: int
    source_atoms: int
    normalized_atoms: int

    def __post_init__(self) -> None:
        """Normalize the expansion matrix and validate normalization counts."""
        self.expansion_matrix = np.asarray(self.expansion_matrix, dtype=np.int64)
        if self.expansion_matrix.shape != (3, 3):
            raise ValueError("expansion_matrix must have shape (3, 3)")
        if self.repetitions <= 0:
            raise ValueError("repetitions must be positive")
        if self.source_atoms < 0 or self.normalized_atoms < 0:
            raise ValueError("atom counts cannot be negative")

    def as_dict(self) -> dict[str, Any]:
        """Return a recursively serializable normalization mapping.

        Returns
        -------
        dict
            Cell-normalization metadata.
        """
        return {
            "basis": self.basis,
            "source_basis": self.source_basis,
            "expansion_matrix": self.expansion_matrix.copy(),
            "repetitions": int(self.repetitions),
            "source_atoms": int(self.source_atoms),
            "normalized_atoms": int(self.normalized_atoms),
        }


@dataclass(slots=True)
class SymmetryMetadata:
    """Symmetry information determined for a reference crystal structure.

    Parameters
    ----------
    space_group_number : int
        International space-group number.
    international_symbol : str
        International short Hermann--Mauguin symbol.
    hall_number : int
        spglib Hall-number identifier.
    hall_symbol : str
        Hall symbol.
    choice : str
        Setting or origin choice reported by spglib.
    point_group : str
        Crystallographic point-group symbol.
    symprec : float
        Cartesian tolerance in angstrom used for symmetry detection.
    angle_tolerance : float
        Angular tolerance in degrees used for symmetry detection. A negative
        value denotes spglib's internal default.
    equivalent_atoms : array_like or None, optional
        Mapping from each atom to its crystallographic equivalence class.
    transformation_matrix : array_like or None, optional
        Transformation matrix reported by spglib for its standard setting.
    origin_shift : array_like or None, optional
        Origin shift reported by spglib.
    """

    space_group_number: int = 0
    international_symbol: str = ""
    hall_number: int = 0
    hall_symbol: str = ""
    choice: str = ""
    point_group: str = ""
    symprec: float = 1.0e-5
    angle_tolerance: float = -1.0
    equivalent_atoms: IntArray | None = None
    transformation_matrix: FloatArray | None = None
    origin_shift: FloatArray | None = None

    def __post_init__(self) -> None:
        """Normalize optional symmetry arrays."""
        if self.equivalent_atoms is not None:
            self.equivalent_atoms = np.asarray(self.equivalent_atoms, dtype=np.int64)
        if self.transformation_matrix is not None:
            self.transformation_matrix = np.asarray(
                self.transformation_matrix,
                dtype=np.float64,
            )
        if self.origin_shift is not None:
            self.origin_shift = np.asarray(self.origin_shift, dtype=np.float64)

    def as_dict(self) -> dict[str, Any]:
        """Return a recursively serializable symmetry mapping.

        Returns
        -------
        dict
            Symmetry identifiers, tolerances, and transformations.
        """
        data: dict[str, Any] = {
            "space_group_number": int(self.space_group_number),
            "international_symbol": self.international_symbol,
            "hall_number": int(self.hall_number),
            "hall_symbol": self.hall_symbol,
            "choice": self.choice,
            "point_group": self.point_group,
            "symprec": float(self.symprec),
            "angle_tolerance": float(self.angle_tolerance),
        }
        if self.equivalent_atoms is not None:
            data["equivalent_atoms"] = self.equivalent_atoms.copy()
        if self.transformation_matrix is not None:
            data["transformation_matrix"] = self.transformation_matrix.copy()
        if self.origin_shift is not None:
            data["origin_shift"] = self.origin_shift.copy()
        return data


@dataclass(slots=True)
class StructureReconstructionDiagnostics:
    """Diagnostics for reduction of a source supercell to a primitive cell.

    Parameters
    ----------
    status : str
        Reconstruction status, normally ``"exact"``, ``"averaged"``, or
        ``"inconsistent"``.
    source_atoms : int
        Number of atoms in the source supercell.
    reconstructed_atoms : int
        Number of atoms in the reconstructed primitive cell.
    expected_repetitions : int
        Expected number of translational copies of each primitive atom.
    minimum_replica_count : int
        Smallest number of source atoms assigned to one primitive atom.
    maximum_replica_count : int
        Largest number of source atoms assigned to one primitive atom.
    maximum_translation_residual : float
        Largest Cartesian displacement from the translationally averaged
        position, in angstrom.
    rms_translation_residual : float
        Root-mean-square displacement from translationally averaged positions,
        in angstrom.
    message : str, optional
        Additional diagnostic explanation.
    """

    status: str
    source_atoms: int
    reconstructed_atoms: int
    expected_repetitions: int
    minimum_replica_count: int
    maximum_replica_count: int
    maximum_translation_residual: float
    rms_translation_residual: float
    message: str = ""

    def as_dict(self) -> dict[str, Any]:
        """Return a recursively serializable diagnostics mapping.

        Returns
        -------
        dict
            Reconstruction counts and residuals.
        """
        return {
            "status": self.status,
            "source_atoms": int(self.source_atoms),
            "reconstructed_atoms": int(self.reconstructed_atoms),
            "expected_repetitions": int(self.expected_repetitions),
            "minimum_replica_count": int(self.minimum_replica_count),
            "maximum_replica_count": int(self.maximum_replica_count),
            "maximum_translation_residual_angstrom": float(
                self.maximum_translation_residual
            ),
            "rms_translation_residual_angstrom": float(self.rms_translation_residual),
            "message": self.message,
        }


@dataclass(slots=True)
class StructureVolumeSeries:
    """Primitive structural path sampled as a function of volume.

    Parameters
    ----------
    reference : CrystalStructure
        Reference primitive structure used to maintain atom correspondence and
        cell orientation across the sampled volumes.
    lattices : array_like
        Primitive direct lattice matrices with shape ``(nvol, 3, 3)``.
    fractional_positions : array_like
        Primitive fractional coordinates with shape
        ``(nvol, natoms, 3)``.
    volumes : array_like
        Primitive-cell volumes with shape ``(nvol,)``.
    normalization : CellNormalization
        Description of the thermodynamic normalization cell.
    symmetry : SymmetryMetadata or None, optional
        Symmetry metadata for the reference structure.
    primitive_to_crystallographic : array_like or None, optional
        CRYSTAL transformation from the primitive to crystallographic cell.
    diagnostics : sequence of StructureReconstructionDiagnostics, optional
        Reconstruction diagnostics for each sampled volume.
    orientation : str, optional
        Description of the Cartesian/basis convention. ``"crystal"`` denotes
        continuity with the CRYSTAL source orientation.
    reference_index : int, optional
        Index of the reference structure within the sampled series.
    metadata : dict, optional
        Additional provenance metadata.
    source_lattices : array_like or None, optional
        Full source-cell lattices retained in memory for optional provenance
        storage.
    source_fractional_positions : sequence or None, optional
        Full source-cell fractional positions retained in memory for optional
        provenance storage. Atom counts may differ among entries, so a tuple is
        used instead of a rectangular array.

    Raises
    ------
    ValueError
        If array shapes or atom counts are inconsistent.
    """

    reference: CrystalStructure
    lattices: FloatArray
    fractional_positions: FloatArray
    volumes: FloatArray
    normalization: CellNormalization
    symmetry: SymmetryMetadata | None = None
    primitive_to_crystallographic: FloatArray | None = None
    diagnostics: tuple[StructureReconstructionDiagnostics, ...] = ()
    orientation: str = "crystal"
    reference_index: int = 0
    metadata: dict[str, Any] = field(default_factory=dict)
    source_lattices: FloatArray | None = None
    source_fractional_positions: tuple[FloatArray, ...] | None = None

    def __post_init__(self) -> None:
        """Normalize structural arrays and validate the sampled path."""
        self.lattices = np.asarray(self.lattices, dtype=np.float64)
        self.fractional_positions = np.asarray(
            self.fractional_positions,
            dtype=np.float64,
        )
        self.volumes = np.asarray(self.volumes, dtype=np.float64)
        if self.lattices.ndim != 3 or self.lattices.shape[1:] != (3, 3):
            raise ValueError("lattices must have shape (nvol, 3, 3)")
        if self.fractional_positions.ndim != 3 or self.fractional_positions.shape[
            1:
        ] != (self.reference.natoms, 3):
            raise ValueError("fractional_positions must have shape (nvol, natoms, 3)")
        if self.volumes.ndim != 1:
            raise ValueError("volumes must be one-dimensional")
        nvol = self.lattices.shape[0]
        if self.fractional_positions.shape[0] != nvol or self.volumes.shape[0] != nvol:
            raise ValueError(
                "structural volume-series arrays have inconsistent lengths"
            )
        if self.reference_index < 0 or self.reference_index >= nvol:
            raise ValueError("reference_index is outside the structural series")
        if self.primitive_to_crystallographic is not None:
            self.primitive_to_crystallographic = np.asarray(
                self.primitive_to_crystallographic,
                dtype=np.float64,
            )
            if self.primitive_to_crystallographic.shape != (3, 3):
                raise ValueError("primitive_to_crystallographic must have shape (3, 3)")
        if self.source_lattices is not None:
            self.source_lattices = np.asarray(self.source_lattices, dtype=np.float64)
            if self.source_lattices.shape != self.lattices.shape:
                raise ValueError("source_lattices must have shape (nvol, 3, 3)")
        if self.source_fractional_positions is not None:
            self.source_fractional_positions = tuple(
                np.asarray(item, dtype=np.float64)
                for item in self.source_fractional_positions
            )
            if len(self.source_fractional_positions) != nvol:
                raise ValueError(
                    "source_fractional_positions must contain one entry per volume"
                )

    @property
    def nvol(self) -> int:
        """Return the number of sampled structures.

        Returns
        -------
        int
            Number of volume points.
        """
        return int(self.volumes.shape[0])

    @property
    def natoms(self) -> int:
        """Return the primitive atom count.

        Returns
        -------
        int
            Number of atoms in each compact primitive structure.
        """
        return self.reference.natoms

    def has_source_supercells(self) -> bool:
        """Return whether full source-cell structures are retained.

        Returns
        -------
        bool
            ``True`` when both source lattices and positions are available.
        """
        return (
            self.source_lattices is not None
            and self.source_fractional_positions is not None
        )

    def as_dict(self, *, include_source: bool = False) -> dict[str, Any]:
        """Return a recursively serializable structural-series mapping.

        Parameters
        ----------
        include_source : bool, optional
            Include full source supercells when they are retained. The compact
            primitive path is always included.

        Returns
        -------
        dict
            Structure, symmetry, normalization, and reconstruction data.
        """
        data: dict[str, Any] = {
            "representation": "primitive",
            "orientation": self.orientation,
            "reference_index": int(self.reference_index),
            "atomic_numbers": self.reference.atomic_numbers.copy(),
            "reference": self.reference.as_dict(),
            "volume_series": {
                "volume": self.volumes.copy(),
                "lattice": self.lattices.copy(),
                "fractional_positions": self.fractional_positions.copy(),
            },
            "normalization": self.normalization.as_dict(),
            "metadata": dict(self.metadata),
        }
        if self.symmetry is not None:
            data["symmetry"] = self.symmetry.as_dict()
        if self.primitive_to_crystallographic is not None:
            data["transformations"] = {
                "primitive_to_crystallographic": (
                    self.primitive_to_crystallographic.copy()
                )
            }
        if self.diagnostics:
            data["reconstruction"] = [item.as_dict() for item in self.diagnostics]
        if include_source and self.has_source_supercells():
            source_lattices = self.source_lattices
            source_positions = self.source_fractional_positions
            if source_lattices is None or source_positions is None:
                raise RuntimeError("source-supercell state is internally inconsistent")
            data["source_series"] = {
                "lattice": source_lattices.copy(),
                "fractional_positions": [item.copy() for item in source_positions],
            }
        return data
