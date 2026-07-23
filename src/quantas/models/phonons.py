# -*- coding: utf-8 -*-

"""Passive data containers for volume-dependent phonon calculations."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from quantas.models.structures import StructureVolumeSeries


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
