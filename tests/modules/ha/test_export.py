from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from quantas.modules.ha.api import read_ha_input, run_ha
from quantas.modules.ha.io.export import (
    HAHDF5Export,
    HATableExport,
    convert_property_values,
)
from quantas.modules.ha.io.reader import read_ha_hdf5
from quantas.modules.ha.models import HAOptions

DATA = Path(__file__).parent / "data/mgo_b3lyp_qha.yaml"


def _result():
    input_data = read_ha_input(DATA)
    options = HAOptions(
        temperature_min=0.0, temperature_max=100.0, temperature_step=100.0
    )
    return run_ha(input_data, options=options)


def test_hdf5_export_writes_quantas2_layout(tmp_path):
    result = _result()
    outfile = tmp_path / "mgo_ha.hdf5"

    exporter = HAHDF5Export()
    exporter.export(result, outfile, report_text="example report")

    assert exporter.completed is True
    with h5py.File(outfile, "r") as h5:
        assert h5["metadata"].attrs["module"] == "ha"
        assert h5["metadata"].attrs["method"] == "harmonic"
        assert h5["metadata"].attrs["schema_version"] == "2.0"
        assert "input" in h5
        assert "options" in h5
        assert "results" in h5
        assert h5["results"].attrs["jobname"] == result.results["ha"].jobname
        assert h5["diagnostics/report_text"].asstr()[()] == "example report"
        assert h5["input/data"].attrs["formula_units"] == 1
        assert h5["results/metadata/units"].attrs["entropy"] == "Ha cell^-1 K^-1"
        assert h5["results/metadata/units"].attrs["heat_capacity"] == "Ha cell^-1 K^-1"
        assert (
            h5["results/metadata"].attrs["thermodynamic_unit_convention"]
            == "native_energy_per_cell_per_kelvin"
        )

        for key in [
            "volume",
            "temperature",
            "static_energy",
            "zero_point_energy",
            "thermal_energy",
            "internal_energy",
            "entropy",
            "vibrational_free_energy",
            "free_energy",
            "isochoric_heat_capacity",
        ]:
            assert key in h5["results"]

        np.testing.assert_allclose(
            h5["results/free_energy"][:], result.results["ha"].free_energy
        )


def test_read_ha_hdf5_roundtrip(tmp_path):
    result = _result()
    outfile = tmp_path / "roundtrip.hdf5"
    HAHDF5Export().export(result, outfile)

    loaded_data = read_ha_hdf5(outfile)
    loaded = loaded_data.results["ha"]

    assert loaded.jobname == result.results["ha"].jobname
    np.testing.assert_allclose(loaded.temperature, result.results["ha"].temperature)
    np.testing.assert_allclose(loaded.volume, result.results["ha"].volume)
    np.testing.assert_allclose(loaded.free_energy, result.results["ha"].free_energy)
    np.testing.assert_allclose(
        loaded.isochoric_heat_capacity, result.results["ha"].isochoric_heat_capacity
    )


def test_table_export_writes_selected_property(tmp_path):
    result = _result()
    outfile = tmp_path / "free_energy.dat"

    exporter = HATableExport()
    exporter.export(result, outfile, property_name="F")

    text = outfile.read_text(encoding="utf8")
    assert exporter.completed is True
    assert "# Quantas HA table export" in text
    assert "Property: F - Helmholtz free energy" in text
    assert "T / K" in text


def test_heat_capacity_export_uses_formula_unit_normalization():
    native = np.array([2.0e-5])
    converted, unit = convert_property_values(
        value=native,
        attr="isochoric_heat_capacity",
        units={"heat_capacity": "Ha cell^-1 K^-1"},
        target_unit="J/mol/K",
        normalization={"formula_units_per_cell": 2},
    )

    expected = native * 2625499.6394798254 / 2.0
    np.testing.assert_allclose(converted, expected, rtol=1.0e-12)
    assert unit == "J mol^-1 K^-1"


def test_read_ha_hdf5_migrates_schema_1_0_implicit_milli_units(tmp_path):
    result = _result()
    ha = result.results["ha"]
    reference_entropy = np.asarray(ha.entropy).copy()
    reference_cv = np.asarray(ha.isochoric_heat_capacity).copy()
    ha.entropy = reference_entropy * 1.0e3
    ha.isochoric_heat_capacity = reference_cv * 1.0e3
    ha.metadata.pop("thermodynamic_unit_convention", None)
    ha.metadata.pop("normalization", None)
    ha.metadata["units"]["entropy"] = "Ha K^-1"
    ha.metadata["units"]["heat_capacity"] = "Ha K^-1"
    result.metadata.schema_version = "1.0"

    outfile = tmp_path / "ha_schema_1_0.hdf5"
    HAHDF5Export().export(result, outfile)
    loaded_data = read_ha_hdf5(outfile)
    loaded = loaded_data.results["ha"]

    np.testing.assert_allclose(loaded.entropy, reference_entropy)
    np.testing.assert_allclose(loaded.isochoric_heat_capacity, reference_cv)
    assert loaded.metadata["units"]["entropy"] == "Ha cell^-1 K^-1"
    assert loaded.metadata["unit_migration"]["scale_factor"] == 1.0e-3
