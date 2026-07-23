from __future__ import annotations

from datetime import datetime, timezone

import h5py
import numpy as np
import pytest

from quantas.models import InputData, ResultData, ResultMetadata
from quantas.modules.qha.io.export import QHAHDF5Export, QHATableExport
from quantas.modules.qha.models import QHAResult


def make_result_data() -> ResultData:
    qha = QHAResult(
        jobname="MgO",
        temperature=np.array([300.0, 400.0]),
        pressure=np.array([0.0, 5.0]),
        volume=np.array([18.0, 19.0, 20.0]),
        static_energy=np.array([-10.0, -10.1, -10.0]),
        equilibrium_volume=np.array([[18.9, 18.5], [19.1, 18.7]]),
        isothermal_bulk_modulus=np.array([[160.0, 170.0], [150.0, 160.0]]),
        gibbs_free_energy=np.array([[-9.0, -8.7], [-8.5, -8.2]]),
        uncertainties={"sigma_VT": np.full((2, 2), 0.02)},
        valid_mask=np.array([[True, True], [True, False]]),
        completed=False,
        metadata={
            "scheme": "freq",
            "minimization": "eos",
            "units": {
                "energy": "Ha",
                "volume": "A^3",
                "pressure": "GPa",
                "temperature": "K",
            },
        },
    )
    return ResultData(
        metadata=ResultMetadata(
            module="qha",
            method="quasi-harmonic approximation",
            created_at=datetime(2026, 1, 1, tzinfo=timezone.utc),
        ),
        input_data=InputData(source="input.yaml", data={"jobname": "MgO"}),
        options={"scheme": "freq", "minimization": "eos"},
        results={"qha": qha},
        warnings=["example warning"],
    )


def test_qha_hdf5_export_writes_metadata_and_datasets(tmp_path) -> None:
    result = make_result_data()
    filename = tmp_path / "qha_result"

    exporter = QHAHDF5Export()
    exporter.export(result, filename, report_text="report")

    output = tmp_path / "qha_result.hdf5"
    assert exporter.completed is True
    assert output.exists()

    with h5py.File(output, "r") as h5:
        assert h5["metadata"].attrs["module"] == "qha"
        assert h5["results"].attrs["jobname"] == "MgO"
        assert h5["results"].attrs["completed"] == np.False_
        np.testing.assert_allclose(
            h5["results/equilibrium_volume"][...], [[18.9, 18.5], [19.1, 18.7]]
        )
        np.testing.assert_allclose(
            h5["results/uncertainties/sigma_VT"][...], np.full((2, 2), 0.02)
        )
        assert "VT" not in h5
        assert h5["diagnostics/report_text"].asstr()[()] == "report"
        assert h5["metadata"].attrs["schema_version"] == "2.0"


def test_qha_table_export_writes_property_with_uncertainty(tmp_path) -> None:
    result = make_result_data()
    filename = tmp_path / "volume_table"

    exporter = QHATableExport()
    exporter.export(result, filename, property_name="VT")

    output = tmp_path / "volume_table.dat"
    assert exporter.completed is True
    text = output.read_text(encoding="utf8")
    assert "Quantas QHA human-readable table export" in text
    assert "sigma_VT" in text
    assert "300.00" in text
    assert "5.00" in text
    assert "18.500000" in text


def test_qha_table_export_rejects_missing_property(tmp_path) -> None:
    result = make_result_data()
    with pytest.raises(ValueError):
        QHATableExport().export(result, tmp_path / "entropy", property_name="S")


def test_qha_export_requires_qha_result(tmp_path) -> None:
    result = make_result_data()
    result.results["qha"] = object()
    with pytest.raises(TypeError):
        QHAHDF5Export().export(result, tmp_path / "bad.hdf5")
