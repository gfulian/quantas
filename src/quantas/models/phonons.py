# -*- coding: utf-8 -*-

"""Passive data containers for volume-dependent phonon calculations."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from quantas.models.structures import StructureVolumeSeries


def default_phonon_input_units() -> dict[str, str]:
    """Return the historical Quantas units for HA/QHA input datasets.

    Returns
    -------
    dict
        Energy, volume, frequency, and structural length unit labels.
    """
    return {
        "energy": "Ha",
        "volume": "angstrom^3",
        "frequency": "cm^-1",
        "length": "angstrom",
    }


@dataclass(slots=True)
class PhononModeData:
    """Phonon frequencies and normalized displacement eigenvectors.

    Parameters
    ----------
    frequencies : array_like
        Phonon frequencies with shape ``(qpoints, nmodes)``.
    eigenvectors : array_like
        Complex displacement eigenvectors with shape
        ``(qpoints, nmodes, natoms, 3)``.
    atom_symbols : tuple of str
        Chemical symbols in the atom ordering used by ``eigenvectors``.
    frequency_unit : str
        Unit of ``frequencies``.
    eigenvector_normalization : str
        Explicit convention used for the stored eigenvectors.
    metadata : dict, optional
        Source and parser provenance.

    Raises
    ------
    ValueError
        If array shapes, atom counts, units, or eigenvector norms are invalid.
    """

    frequencies: np.ndarray
    eigenvectors: np.ndarray
    atom_symbols: tuple[str, ...]
    frequency_unit: str = "cm^-1"
    eigenvector_normalization: str = "mass-weighted-unit"
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays to Quantas dtypes and validate mode dimensions."""
        self.frequencies = np.asarray(self.frequencies, dtype=np.float64)
        self.eigenvectors = np.asarray(self.eigenvectors, dtype=np.complex128)
        self.atom_symbols = tuple(str(symbol).strip() for symbol in self.atom_symbols)
        self.frequency_unit = str(self.frequency_unit).strip()
        self.eigenvector_normalization = str(self.eigenvector_normalization).strip()
        if self.frequencies.ndim != 2:
            raise ValueError("phonon frequencies must have shape (qpoints, nmodes)")
        if self.eigenvectors.ndim != 4 or self.eigenvectors.shape[-1] != 3:
            raise ValueError(
                "phonon eigenvectors must have shape (qpoints, nmodes, natoms, 3)"
            )
        if self.eigenvectors.shape[:2] != self.frequencies.shape:
            raise ValueError(
                "phonon frequencies and eigenvectors have incompatible shapes"
            )
        if self.eigenvectors.shape[2] != len(self.atom_symbols):
            raise ValueError(
                "atom symbols and phonon eigenvectors use different atom counts"
            )
        if not self.frequency_unit:
            raise ValueError("phonon frequency unit cannot be empty")
        if not self.eigenvector_normalization:
            raise ValueError("phonon eigenvector normalization cannot be empty")
        norms = np.linalg.norm(
            self.eigenvectors.reshape(
                self.eigenvectors.shape[0],
                self.eigenvectors.shape[1],
                -1,
            ),
            axis=2,
        )
        if not np.all(np.isfinite(norms)) or not np.allclose(
            norms,
            1.0,
            rtol=1.0e-10,
            atol=1.0e-12,
        ):
            raise ValueError("phonon eigenvectors must be normalized to unit norm")

    @property
    def qpoints(self) -> int:
        """Return the number of stored phonon q-points."""
        return int(self.frequencies.shape[0])

    @property
    def nmodes(self) -> int:
        """Return the number of modes stored at each q-point."""
        return int(self.frequencies.shape[1])

    @property
    def natoms(self) -> int:
        """Return the number of atoms represented by each eigenvector."""
        return int(self.eigenvectors.shape[2])


@dataclass(slots=True)
class PhononInputData:
    """Normalized structural, energetic, and phonon input data.

    Parameters
    ----------
    jobname : str, optional
        Name or short description of the calculation.
    natoms : int, optional
        Number of atoms in the thermodynamic normalization cell.
    formula_units : int, optional
        Number of chemical formula units in the thermodynamic normalization cell.
    supercell : ndarray or None, optional
        Supercell matrix used to sample phonons, with shape ``(3, 3)``.
    qpoints : int, optional
        Number of phonon q-points stored in the input data.
    volume : ndarray or None, optional
        Unit-cell volumes with shape ``(nvol,)``.
    energy : ndarray or None, optional
        Static energies with shape ``(nvol,)``.
    frequencies : ndarray or None, optional
        Phonon frequencies with shape ``(qpoints, nmodes, nvol)``.
    weights : ndarray or None, optional
        Q-point weights with shape ``(qpoints,)``.
    qcoords : ndarray or None, optional
        Q-point fractional coordinates with shape ``(qpoints, 3)``.
    structure : StructureVolumeSeries or None, optional
        Compact primitive structural path, symmetry, and normalization data.
    units : dict, optional
        Explicit physical units carried by the input dataset. Common keys are
        ``"energy"``, ``"volume"``, ``"frequency"``, and ``"length"``.
    source : str or Path or None, optional
        Source file from which the data were read.
    metadata : dict, optional
        Additional parser-specific or workflow-specific metadata.
    """

    jobname: str = "Unknown"
    natoms: int = 0
    formula_units: int = 1
    supercell: np.ndarray | None = None
    qpoints: int = 0

    volume: np.ndarray | None = None
    energy: np.ndarray | None = None
    frequencies: np.ndarray | None = None
    weights: np.ndarray | None = None
    qcoords: np.ndarray | None = None
    structure: StructureVolumeSeries | None = None
    units: dict[str, str] = field(default_factory=default_phonon_input_units)

    source: str | Path | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def natoms_per_formula_unit(self) -> float:
        """Return the number of atoms per chemical formula unit.

        Returns
        -------
        float
            Number of atoms in one formula unit.

        Raises
        ------
        ValueError
            If ``formula_units`` is not positive.
        """
        if self.formula_units <= 0:
            raise ValueError("formula_units must be positive")
        return float(self.natoms) / float(self.formula_units)

    @property
    def nvol(self) -> int:
        """Return the number of sampled volumes.

        Returns
        -------
        int
            Number of volume points, or ``0`` when no volume array is stored.
        """
        if self.volume is None:
            return 0
        return int(np.asarray(self.volume).shape[0])

    @property
    def nmodes(self) -> int:
        """Return the number of phonon modes per q-point.

        Returns
        -------
        int
            Number of modes, or ``0`` when no frequency array is stored.
        """
        if self.frequencies is None:
            return 0
        return int(np.asarray(self.frequencies).shape[1])

    @property
    def kpoints(self) -> int:
        """Return the number of sampled k-points from the supercell matrix.

        Returns
        -------
        int
            Rounded determinant of the supercell matrix, or ``0`` when absent.
        """
        if self.supercell is None:
            return 0
        return int(np.around(np.linalg.det(np.asarray(self.supercell)), 0))

    @property
    def total_q_points(self) -> float:
        """Return the sum of q-point weights.

        Returns
        -------
        float
            Sum of the stored weights, or ``0.0`` when no weights are stored.
        """
        if self.weights is None:
            return 0.0
        return float(np.asarray(self.weights, dtype=np.float64).sum())

    def normalized_weights(self) -> np.ndarray:
        """Return q-point weights normalized by their sum.

        Returns
        -------
        ndarray
            Normalized q-point weights.

        Raises
        ------
        ValueError
            If weights are unavailable or their sum is not positive.
        """
        if self.weights is None:
            raise ValueError("q-point weights are not available")
        weights = np.asarray(self.weights, dtype=np.float64)
        total = float(weights.sum())
        if total <= 0.0:
            raise ValueError("the sum of q-point weights must be positive")
        return weights / total

    def has_structure(self) -> bool:
        """Return whether a compact structural volume series is available.

        Returns
        -------
        bool
            ``True`` when structural lattices and atomic coordinates are
            stored alongside phonon data.
        """
        return self.structure is not None

    def has_phonons(self) -> bool:
        """Return whether phonon frequencies are available.

        Returns
        -------
        bool
            ``True`` when a frequency array is stored.
        """
        return self.frequencies is not None
