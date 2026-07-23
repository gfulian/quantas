"""Tests for the public HA Python API."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.events import EventLevel, ListObserver
from quantas.io import PhononInputFileReader
from quantas.models import ResultData
from quantas.renderers.tables.text import render_tables
from quantas.modules.ha.api import (
    build_ha_report,
    normalize_ha_input,
    read_ha_input,
    run_ha,
)
from quantas.modules.ha.models import HAInput, HAOptions


DATA = Path(__file__).parent / "data/mgo_b3lyp_qha.yaml"


def test_read_ha_input_returns_normalized_input() -> None:
    """read_ha_input should convert a YAML file to HAInput."""
    ha_input = read_ha_input(DATA)

    assert isinstance(ha_input, HAInput)
    assert ha_input.jobname == "Quasi-Harmonic analysis of periclase (MgO)"
    assert ha_input.natoms == 2
    assert ha_input.qpoints == 32
    assert ha_input.nvol == 11
    assert ha_input.frequencies.shape == (32, 6, 11)
    assert ha_input.source == DATA


def test_read_ha_input_raises_for_missing_file(tmp_path: Path) -> None:
    """read_ha_input should expose reader failures as ValueError."""
    missing = tmp_path / "missing.yaml"

    with pytest.raises(ValueError, match="Unable to read input file"):
        read_ha_input(missing)


def test_run_ha_accepts_input_object() -> None:
    """run_ha should accept an already normalized HAInput object."""
    ha_input = read_ha_input(DATA)
    options = HAOptions(
        temperature_min=100.0,
        temperature_max=300.0,
        temperature_step=100.0,
    )

    result = run_ha(ha_input, options=options)

    assert isinstance(result, ResultData)
    assert result.metadata.module == "ha"
    assert result.metadata.method == "harmonic"
    assert result.results["ha"].temperature.tolist() == [100.0, 200.0, 300.0]
    assert result.results["ha"].jobname == "Quasi-Harmonic analysis of periclase (MgO)"


def test_run_ha_accepts_file_path() -> None:
    """run_ha should read the YAML file when a path is passed."""
    result = run_ha(
        DATA,
        options=HAOptions(
            temperature_min=298.15,
            temperature_max=298.15,
            temperature_step=1.0,
        ),
    )

    assert isinstance(result, ResultData)
    assert result.results["ha"].volume.shape == (11,)
    assert result.results["ha"].zero_point_energy.shape == (1, 11)
    assert result.results["ha"].thermal_energy.shape == (1, 11)
    assert result.results["ha"].free_energy.shape == (1, 11)


def test_run_ha_emits_events_when_observer_is_provided() -> None:
    """run_ha should forward observers to HACalculator."""
    observer = ListObserver()

    run_ha(
        DATA,
        options=HAOptions(
            temperature_min=100.0,
            temperature_max=100.0,
            temperature_step=1.0,
        ),
        observer=observer,
    )

    assert any(event.level is EventLevel.INFO for event in observer.events)
    assert any(event.level is EventLevel.PROGRESS for event in observer.events)
    assert any(event.level is EventLevel.RESULT for event in observer.events)


def test_run_ha_is_silent_by_default(capsys: pytest.CaptureFixture[str]) -> None:
    """The API workflow should not print when no observer is supplied."""
    run_ha(
        DATA,
        options=HAOptions(
            temperature_min=298.15,
            temperature_max=298.15,
            temperature_step=1.0,
        ),
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_run_ha_rejects_unsupported_input_type() -> None:
    """run_ha should fail explicitly for unsupported input types."""
    with pytest.raises(TypeError, match="PhononInputData instance"):
        run_ha(object())  # type: ignore[arg-type]


def test_normalize_ha_input_accepts_neutral_phonon_data() -> None:
    """HA accepts the shared phonon input container directly."""
    phonon_data = PhononInputFileReader(DATA).to_input()

    normalized = normalize_ha_input(phonon_data)

    assert normalized.jobname == phonon_data.jobname
    assert normalized.qpoints == phonon_data.qpoints
    assert normalized.qcoords is not None
    np.testing.assert_allclose(normalized.qcoords, phonon_data.qcoords)


def test_public_report_supports_multivolume_results() -> None:
    result = run_ha(
        DATA,
        options=HAOptions(
            temperature_min=0.0,
            temperature_max=300.0,
            temperature_step=100.0,
        ),
    )

    tables = build_ha_report(result)
    text = render_tables(tables)

    assert any(table.title == "Static and zero-point energies" for table in tables)
    assert any(table.title == "Helmholtz free energy" for table in tables)
    assert "V=15.495427" in text
    assert "V=21.164485" in text
    assert "property temperature dimension is incompatible" not in text
