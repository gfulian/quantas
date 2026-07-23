from pathlib import Path

import numpy as np
import pytest

from quantas.core.events import ListObserver
from quantas.models import PhononInputData, ResultData
from quantas.modules.qha.api import normalize_qha_input, read_qha_input, run_qha
from quantas.modules.qha.io.reader import QHAInputFileReader, phonon_to_qha_input
from quantas.modules.qha.models import QHAInput, QHAOptions


def _qha_input() -> QHAInput:
    return QHAInput(
        jobname="api test",
        natoms=1,
        qpoints=1,
        volume=np.array([9.8, 10.0, 10.2, 10.4], dtype=float),
        energy=np.array([0.040, 0.000, 0.035, 0.140], dtype=float),
        frequencies=np.ones((1, 3, 4), dtype=float),
        weights=np.array([1.0], dtype=float),
    )


def _phonon_input() -> PhononInputData:
    return PhononInputData(
        jobname="ha source",
        natoms=1,
        qpoints=1,
        volume=np.array([9.8, 10.0, 10.2, 10.4], dtype=float),
        energy=np.array([0.040, 0.000, 0.035, 0.140], dtype=float),
        frequencies=np.ones((1, 3, 4), dtype=float),
        weights=np.array([1.0], dtype=float),
        metadata={"format": "test"},
    )


def _options() -> QHAOptions:
    return QHAOptions(
        temperature_min=300.0,
        temperature_max=300.0,
        pressure_min=0.0,
        pressure_max=0.0,
        scheme="td",
        minimization="poly",
        free_energy_degree=3,
    )


def test_phonon_input_can_be_converted_to_qha_input():
    phonon_input = _phonon_input()

    qha_input = phonon_to_qha_input(phonon_input)

    assert qha_input.jobname == "ha source"
    assert qha_input.nvol == 4
    assert qha_input.qpoints == 1
    assert qha_input.nmodes == 3
    np.testing.assert_allclose(qha_input.volume, phonon_input.volume)
    np.testing.assert_allclose(qha_input.energy, phonon_input.energy)
    np.testing.assert_allclose(qha_input.frequencies, phonon_input.frequencies)
    assert qha_input.metadata["source_format"] == "quantas-phonon-yaml"


def test_normalize_qha_input_accepts_qha_input():
    qha_input = _qha_input()

    normalized = normalize_qha_input(qha_input)

    assert normalized is qha_input


def test_normalize_qha_input_accepts_phonon_input():
    normalized = normalize_qha_input(_phonon_input())

    assert isinstance(normalized, QHAInput)
    assert normalized.jobname == "ha source"


def test_normalize_qha_input_rejects_unknown_objects():
    with pytest.raises(TypeError):
        normalize_qha_input(object())


def test_run_qha_accepts_qha_input_and_returns_result_data():
    observer = ListObserver()

    result = run_qha(_qha_input(), _options(), observer=observer)

    assert isinstance(result, ResultData)
    assert result.metadata.module == "qha"
    assert result.results["qha"].completed is True
    assert result.results["qha"].equilibrium_volume.shape == (1, 1)
    assert observer.events


def test_qha_reader_reads_existing_yaml_input():
    filename = Path(__file__).parent / "data/mgo_b3lyp_qha.yaml"
    if not filename.exists():
        pytest.skip("shared MgO QHA input is not available")

    qha_input = read_qha_input(filename)

    assert qha_input.nvol == 11
    assert qha_input.qpoints == 32
    assert qha_input.has_mode_continuity()


def test_qha_input_file_reader_exposes_qha_input_conversion():
    filename = Path(__file__).parent / "data/mgo_b3lyp_qha.yaml"
    if not filename.exists():
        pytest.skip("shared MgO QHA input is not available")

    reader = QHAInputFileReader(filename)
    qha_input = reader.to_qha_input()

    assert reader.completed
    assert qha_input.nvol == 11
    assert qha_input.has_phonons()
