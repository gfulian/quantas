"""Tests for the Quantas EOS keyword-driven input reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.modules.eos import EOSDataset, EOSInputFileReader, read_eos_input


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "eos.dat"
    path.write_text(text, encoding="utf-8")
    return path


def test_original_quantas_layout_is_normalized(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        """
        # Quartz dataset
        JOB
        Quartz experimental compression
        FORMAT
        P V sigmaP sigmaV
        DATA
        0.001 112.981 0.000001 0.002
        0.430 111.725 0.008000 0.014
        """,
    )

    dataset = read_eos_input(path)

    assert dataset.jobname == "Quartz experimental compression"
    assert dataset.column_names == (
        "pressure",
        "volume",
        "sigma_pressure",
        "sigma_volume",
    )
    assert dataset.npoints == 2
    assert dataset.source == path
    assert all(values.dtype == np.float64 for values in dataset.columns.values())
    np.testing.assert_allclose(dataset.column("pressure"), [0.001, 0.430])
    np.testing.assert_allclose(dataset.column("volume"), [112.981, 111.725])
    assert dataset.units["pressure"] == "GPa"
    assert dataset.units["volume"] == "angstrom^3"


def test_eosfit_style_inline_format_and_comma_table_are_supported(
    tmp_path: Path,
) -> None:
    path = _write(
        tmp_path,
        """
        TITLE Kalsilite high-P data
        SYSTEM hexagonal
        TSCALE C
        FORMAT 1 T P SIGP A SIGA C SIGC
        25,0.000,0.004,5.16026,0.00021,8.71661,0.00021
        50,0.328,0.004,5.14963,0.00020,8.70565,0.00019
        """,
    )

    dataset = read_eos_input(path)

    assert dataset.jobname == "Kalsilite high-P data"
    assert dataset.metadata["crystal_system"] == "hexagonal"
    assert dataset.metadata["input_temperature_scale"] == "C"
    np.testing.assert_allclose(dataset.column("temperature"), [298.15, 323.15])
    np.testing.assert_allclose(dataset.raw_columns["temperature"], [25.0, 50.0])
    assert dataset.raw_units["temperature"] == "C"
    np.testing.assert_allclose(dataset.column("a"), [5.16026, 5.14963])
    np.testing.assert_allclose(dataset.column("c"), [8.71661, 8.70565])
    assert dataset.units["temperature"] == "K"


def test_format_order_controls_table_mapping(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        """
        FORMAT V SIGV P SIGP B SIGB
        112.0 0.02 0.0 0.01 7.0 0.001
        110.0 0.03 2.0 0.02 6.9 0.002
        """,
    )

    dataset = read_eos_input(path)

    np.testing.assert_allclose(dataset.column("pressure"), [0.0, 2.0])
    np.testing.assert_allclose(dataset.column("sigma_pressure"), [0.01, 0.02])
    np.testing.assert_allclose(dataset.column("b"), [7.0, 6.9])


def test_long_and_short_aliases_are_canonicalized(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        """
        FORMAT pressure volume sigp sigv temperature sigt
        DATA
        0.0 100.0 0.01 0.02 300.0 0.5
        """,
    )

    dataset = read_eos_input(path)

    assert set(dataset.columns) == {
        "pressure",
        "volume",
        "sigma_pressure",
        "sigma_volume",
        "temperature",
        "sigma_temperature",
    }


def test_fahrenheit_values_and_uncertainties_are_converted_to_kelvin(
    tmp_path: Path,
) -> None:
    path = _write(
        tmp_path,
        """
        TSCALE F
        FORMAT T SIGT V
        32.0 1.8 100.0
        212.0 3.6 101.0
        """,
    )

    dataset = read_eos_input(path)

    np.testing.assert_allclose(dataset.column("temperature"), [273.15, 373.15])
    np.testing.assert_allclose(dataset.column("sigma_temperature"), [1.0, 2.0])
    np.testing.assert_allclose(dataset.raw_columns["temperature"], [32.0, 212.0])
    np.testing.assert_allclose(dataset.raw_columns["sigma_temperature"], [1.8, 3.6])
    assert dataset.raw_units["temperature"] == "F"


def test_units_and_provenance_metadata_are_retained(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        """
        PROVENANCE Diamond-anvil-cell run 4
        UNITS P=kbar V=nm^3 SIGP=kbar SIGV=nm^3
        FORMAT P V SIGP SIGV
        0.0 0.10 0.1 0.001
        """,
    )

    dataset = read_eos_input(path)

    assert dataset.units == {
        "pressure": "GPa",
        "volume": "angstrom^3",
        "sigma_pressure": "GPa",
        "sigma_volume": "angstrom^3",
    }
    assert dataset.raw_units == {
        "pressure": "kbar",
        "volume": "nm^3",
        "sigma_pressure": "kbar",
        "sigma_volume": "nm^3",
    }
    np.testing.assert_allclose(dataset.column("pressure"), [0.0])
    np.testing.assert_allclose(dataset.column("sigma_pressure"), [0.01])
    np.testing.assert_allclose(dataset.column("volume"), [100.0])
    np.testing.assert_allclose(dataset.column("sigma_volume"), [1.0])
    assert set(dataset.provenance.values()) == {"Diamond-anvil-cell run 4"}


def test_dataset_selects_pressure_volume_series_with_uncertainties(
    tmp_path: Path,
) -> None:
    path = _write(
        tmp_path,
        """
        FORMAT P V SIGP SIGV A SIGA
        0.0 100.0 0.01 0.02 5.0 0.001
        1.0 99.0 0.02 0.03 4.9 0.002
        """,
    )
    dataset = read_eos_input(path)

    series = dataset.series("volume", mask=np.array([True, False]))

    assert series.independent == "pressure"
    assert series.target == "volume"
    assert series.size == 2
    assert series.selected_size == 1
    np.testing.assert_allclose(series.sigma_x, [0.01, 0.02])
    np.testing.assert_allclose(series.sigma_y, [0.02, 0.03])


def test_reader_status_is_updated_on_success_and_failure(tmp_path: Path) -> None:
    valid = _write(tmp_path, "FORMAT P V\n0.0 10.0\n")
    reader = EOSInputFileReader()

    dataset = reader.load(valid)

    assert reader.completed is True
    assert reader.error is None
    assert reader.dataset is dataset

    invalid = _write(tmp_path, "FORMAT P V\nnot-a-number 10.0\n")
    with pytest.raises(ValueError, match="non-numeric"):
        reader.load(invalid)
    assert reader.completed is False
    assert reader.error is not None
    assert reader.dataset is None


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("DATA\n0 1\n", "FORMAT must be declared"),
        ("FORMAT\nP P\n0 1\n", "duplicate columns"),
        ("FORMAT P UNKNOWN\n0 1\n", "Unsupported EOS FORMAT column"),
        ("FORMAT 2 P V\n0 1\n", "one numeric line per data record"),
        ("FORMAT P V\n0 1 2\n", "expected 2"),
        ("FORMAT P V\n", "does not contain any numeric data"),
        ("UNITS P GPa\nFORMAT P V\n0 1\n", "COLUMN=UNIT"),
        ("TSCALE rankine\nFORMAT T V\n1 1\n", "temperature scale"),
    ],
)
def test_invalid_input_is_rejected(
    tmp_path: Path,
    text: str,
    message: str,
) -> None:
    path = _write(tmp_path, text)

    with pytest.raises(ValueError, match=message):
        read_eos_input(path)


def test_negative_uncertainty_is_rejected(tmp_path: Path) -> None:
    path = _write(tmp_path, "FORMAT P V SIGP\n0.0 10.0 -0.1\n")

    with pytest.raises(ValueError, match="cannot be negative"):
        read_eos_input(path)


def test_uncertainty_without_measured_quantity_is_rejected(tmp_path: Path) -> None:
    path = _write(tmp_path, "FORMAT P SIGV\n0.0 0.1\n")

    with pytest.raises(ValueError, match="requires 'volume'"):
        read_eos_input(path)


def test_dataset_rejects_inconsistent_direct_construction() -> None:
    with pytest.raises(ValueError, match="equal length"):
        EOSDataset(
            columns={
                "pressure": np.array([0.0, 1.0]),
                "volume": np.array([10.0]),
            }
        )


def test_real_eosfit7_quartz_file_is_read_without_preprocessing() -> None:
    path = Path(__file__).with_name("data") / "PV_quartz.dat"

    dataset = read_eos_input(path)

    assert dataset.jobname == "PV_Quartz"
    assert dataset.npoints == 22
    assert dataset.column_names == (
        "pressure",
        "volume",
        "sigma_pressure",
        "sigma_volume",
    )
    np.testing.assert_allclose(dataset.column("pressure")[[0, -1]], [0.0, 8.905])
    np.testing.assert_allclose(dataset.column("volume")[[0, -1]], [112.981, 96.989])
    assert dataset.metadata["format_columns"] == dataset.column_names


def test_original_quantas_topaz_file_retains_volume_and_axes() -> None:
    path = Path(__file__).with_name("data") / "PV_topaz.dat"

    dataset = read_eos_input(path)

    assert dataset.jobname == "Topaz (Gatta et al., 2014)"
    assert dataset.npoints == 25
    assert dataset.column_names == (
        "pressure",
        "sigma_pressure",
        "volume",
        "sigma_volume",
        "a",
        "sigma_a",
        "b",
        "sigma_b",
        "c",
        "sigma_c",
    )
    np.testing.assert_allclose(dataset.column("volume")[[0, -1]], [345.46, 281.70])
    np.testing.assert_allclose(dataset.column("a")[[0, -1]], [4.6627, 4.360])
    np.testing.assert_allclose(dataset.column("b")[[0, -1]], [8.8343, 8.348])
    np.testing.assert_allclose(dataset.column("c")[[0, -1]], [8.3867, 7.739])
