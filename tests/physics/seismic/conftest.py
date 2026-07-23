"""Shared fixtures for SEISMIC formula-reference tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tests.reference.seismic_reference import SeismicFormulaReference


DATA = Path(__file__).parent / "data" / "hydroxylapatite.dat"


def read_elasticity_input(path: Path) -> tuple[np.ndarray, float]:
    """Read a triangular Quantas elasticity fixture.

    Parameters
    ----------
    path : pathlib.Path
        Input data path.

    Returns
    -------
    tuple
        Symmetric stiffness matrix in GPa and density in kg m-3.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    matrix = np.zeros((6, 6), dtype=float)
    for row, line in enumerate(lines[1:7]):
        values = np.asarray([float(value) for value in line.split()], dtype=float)
        matrix[row, row:] = values
    matrix += np.triu(matrix, 1).T
    return matrix, float(lines[7])


@pytest.fixture(scope="module")
def hydroxylapatite_data() -> tuple[np.ndarray, float]:
    """Return the characterized hydroxylapatite stiffness and density."""
    return read_elasticity_input(DATA)


@pytest.fixture(scope="module")
def hydroxylapatite_reference(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> SeismicFormulaReference:
    """Return the frozen formula reference for hydroxylapatite."""
    stiffness, density = hydroxylapatite_data
    return SeismicFormulaReference(stiffness, density)


@pytest.fixture(scope="module")
def generic_directions() -> np.ndarray:
    """Return non-symmetry directions used by directional tests."""
    directions = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [2.0, -1.0, 3.0],
            [-2.0, 3.0, 1.0],
            [1.0, -3.0, 2.0],
        ],
        dtype=float,
    )
    return directions / np.linalg.norm(directions, axis=1)[:, None]
