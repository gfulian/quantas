"""V-T input compatibility and dataset-classification tests."""

from pathlib import Path

import numpy as np
import pytest

from quantas.modules.eos import (
    EOSCoordinateVariation,
    EOSFitRequest,
    EOSFitter,
    read_eos_input,
)


_DATA = Path(__file__).with_name("data")


def test_rutile_reader_accepts_comments_colon_format_and_normalized_scales() -> None:
    dataset = read_eos_input(_DATA / "rutile.dat")

    assert dataset.jobname == "Rutile Unit cell data as f(T)"
    assert dataset.npoints == 59
    assert dataset.column_names == (
        "temperature",
        "sigma_temperature",
        "pressure",
        "sigma_pressure",
        "volume",
        "sigma_volume",
        "a",
        "sigma_a",
        "c",
        "sigma_c",
    )
    assert dataset.units["temperature"] == "K"
    assert dataset.units["pressure"] == "GPa"
    for name in ("volume", "sigma_volume", "a", "sigma_a", "c", "sigma_c"):
        assert dataset.units[name] == "dimensionless"
    assert dataset.metadata["input_volume_scale"] == "V/V0"
    assert dataset.metadata["input_linear_scale"] == "L/L0"
    assert dataset.metadata["column_scales"]["volume"] == "V/V0"
    assert dataset.metadata["column_scales"]["a"] == "L/L0"
    comments = dataset.metadata["comments"]
    assert "Data selected by Zaffiro et al (2019)" in comments[0]
    assert any("Hummer et al" in comment for comment in comments)
    np.testing.assert_allclose(dataset.column("pressure"), 1.0e-4)


def test_rutile_is_classified_as_isobaric_vt_data() -> None:
    dataset = read_eos_input(_DATA / "rutile.dat")

    classification = dataset.classify()

    pressure = classification.profile("pressure")
    assert pressure.variation is EOSCoordinateVariation.CONSTANT
    assert pressure.reference_value == pytest.approx(1.0e-4)
    assert classification.is_isobaric is True
    assert classification.reference_pressure == pytest.approx(1.0e-4)
    assert classification.is_isothermal is False
    assert classification.reference_temperature is None
    assert classification.variable_coordinates == ("temperature", "volume", "a", "c")
    assert classification.constant_coordinates == ("pressure",)


def test_rutile_provides_temperature_volume_and_linear_series() -> None:
    dataset = read_eos_input(_DATA / "rutile.dat")

    volume = dataset.series("volume", independent="temperature")
    axis_a = dataset.series("a", independent="temperature")

    assert volume.x.size == 59
    assert volume.units["temperature"] == "K"
    assert volume.units["volume"] == "dimensionless"
    assert volume.metadata["target_scale"] == "V/V0"
    assert axis_a.units["a"] == "dimensionless"
    assert axis_a.metadata["target_scale"] == "L/L0"
    assert volume.sigma_x is not None
    assert volume.sigma_y is not None


def test_constant_pressure_is_rejected_before_pv_solver_dispatch() -> None:
    dataset = read_eos_input(_DATA / "rutile.dat")

    with pytest.raises(ValueError, match="pressure column is constant"):
        EOSFitter().fit(dataset, EOSFitRequest(model="BM3"))


def test_triclinic_reader_accepts_simple_vt_table_without_pressure() -> None:
    dataset = read_eos_input(_DATA / "T_triclinic.dat")

    assert dataset.npoints == 186
    assert dataset.column_names == ("temperature", "volume")
    assert dataset.units == {"temperature": "K", "volume": "angstrom^3"}
    classification = dataset.classify()
    assert classification.variable_coordinates == ("temperature", "volume")
    assert classification.constant_coordinates == ()
    assert classification.is_isobaric is False
    assert classification.reference_pressure is None
    vt = dataset.series("volume", independent="temperature")
    assert vt.x[0] == pytest.approx(91.8534)
    assert vt.y[-1] == pytest.approx(678.902)


def test_coordinate_classification_respects_mask() -> None:
    dataset = read_eos_input(_DATA / "PV_quartz.dat")
    mask = np.zeros(dataset.npoints, dtype=bool)
    mask[:2] = True
    pressure = dataset.column("pressure")
    pressure[1] = pressure[0]

    profile = dataset.coordinate_profile("pressure", mask=mask)

    assert profile.is_constant
    assert profile.npoints == 2


def test_normalized_scale_rejects_dimensional_unit_declaration(tmp_path: Path) -> None:
    source = tmp_path / "bad_scale.dat"
    source.write_text(
        "TITLE bad scale\nVSCALE V/V0\nUNITS V=angstrom^3\nFORMAT T V\n300 1.0\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="VSCALE declares normalized data"):
        read_eos_input(source)
