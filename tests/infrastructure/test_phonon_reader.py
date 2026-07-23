"""Tests for the workflow-independent phonon YAML reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.io import PhononInputFileReader
from quantas.models import PhononInputData

DATA = (
    Path(__file__).resolve().parents[1]
    / "modules"
    / "ha"
    / "data"
    / "mgo_b3lyp_qha.yaml"
)


def test_shared_reader_returns_neutral_phonon_data():
    reader = PhononInputFileReader(DATA)
    result = reader.to_input()

    assert reader.completed is True
    assert isinstance(result, PhononInputData)
    assert result.natoms == 2
    assert result.qpoints == 32
    assert result.nvol == 11
    assert result.frequencies.shape == (32, 6, 11)
    assert result.qcoords is not None
    assert result.qcoords.shape == (32, 3)
    np.testing.assert_allclose(result.qcoords[0], [0.0, 0.0, 0.0])
    assert result.metadata["format"] == "quantas-phonon-yaml"


def test_shared_reader_requires_a_completed_parse(tmp_path):
    reader = PhononInputFileReader()

    with pytest.raises(ValueError, match="Phonon input data"):
        reader.to_input()

    broken = tmp_path / "broken.yaml"
    broken.write_text("volume: 10.0\n", encoding="utf-8")
    reader.load(broken)

    assert reader.completed is False
    with pytest.raises(ValueError, match="No energy values"):
        reader.to_input()
