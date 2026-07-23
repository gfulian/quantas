"""Shared fixtures for elasticity-core tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor


DATA = Path(__file__).parent / "data" / "hydroxylapatite.dat"


def _read_reference_matrix() -> np.ndarray:
    lines = DATA.read_text(encoding="utf-8").splitlines()
    rows = lines[1:7]
    matrix = np.zeros((6, 6), dtype=float)
    for i, line in enumerate(rows):
        values = [float(value) for value in line.split()]
        matrix[i, i:] = values
    matrix += np.triu(matrix, 1).T
    return matrix


@pytest.fixture
def hydroxylapatite_tensor() -> ElasticTensor:
    """Return the characterized hydroxylapatite elastic tensor."""
    return ElasticTensor(_read_reference_matrix())
