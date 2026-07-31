"""End-to-end public lifecycle contracts outside the CLI layer."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from quantas.api import elasticity, ha, qha, seismic
from quantas.api.common import ResultData, ResultMetadata
from quantas.cli.elasticity import elasticity as elasticity_cli
from quantas.cli.ha import ha as ha_cli
from quantas.cli.qha import qha as qha_cli
from quantas.cli.seismic import seismic as seismic_cli


CRYSTAL_ELASTICITY_OUTPUT = (
    Path(__file__).parents[1]
    / "interfaces"
    / "data"
    / "calcite_crystal_elastcon_excerpt.out"
)


def _write_vasp_elasticity_output(
    path: Path,
    *,
    include_density: bool = True,
) -> Path:
    """Write a compact VASP elasticity output accepted by the public reader."""
    matrix = np.diag([2000.0, 2100.0, 2200.0, 700.0, 800.0, 900.0])
    lines: list[str] = []
    if include_density:
        lines.extend(
            [
                "POMASS = 24.305; ZVAL = 2.000",
                "ions per type = 2",
                "volume of cell : 80.000000",
            ]
        )
    lines.extend(
        [
            "SYMMETRIZED ELASTIC MODULI (kBar)",
            "separator",
            "separator",
        ]
    )
    lines.extend(
        f"{row + 1} " + " ".join(str(value) for value in matrix[row])
        for row in range(6)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_elasticity_public_input_generator_is_runnable(tmp_path: Path) -> None:
    """The CLI input generator has a complete public Python equivalent."""
    output = elasticity.create_input(
        CRYSTAL_ELASTICITY_OUTPUT,
        tmp_path / "calcite-input",
        interface="crystal",
        jobname="Calcite public API",
    )

    assert output == tmp_path / "calcite-input.dat"
    parsed = elasticity.read_input(output)
    assert parsed.jobname == "Calcite public API"
    assert parsed.stiffness is not None
    assert parsed.stiffness[0, 4] == 20.67


def test_elasticity_cli_input_matches_public_api(tmp_path: Path) -> None:
    """Elasticity CLI and API generate identical normalized input text."""
    api_path = elasticity.create_input(
        CRYSTAL_ELASTICITY_OUTPUT,
        tmp_path / "api-input",
        interface="crystal",
        jobname="Public lifecycle",
    )
    cli_path = tmp_path / "cli-input.dat"
    response = CliRunner().invoke(
        elasticity_cli,
        [
            "inpgen",
            str(CRYSTAL_ELASTICITY_OUTPUT),
            "--output",
            str(cli_path),
            "--interface",
            "crystal",
        ],
        input="Public lifecycle\n",
    )
    assert response.exit_code == 0, response.output
    assert cli_path.read_bytes() == api_path.read_bytes()


def test_seismic_cli_input_matches_public_api(tmp_path: Path) -> None:
    """SEISMIC CLI and API generate identical normalized input text."""
    api_path = seismic.create_input(
        CRYSTAL_ELASTICITY_OUTPUT,
        tmp_path / "api-seismic",
        interface="crystal",
        jobname="SEISMIC lifecycle",
    )
    cli_path = tmp_path / "cli-seismic.dat"
    response = CliRunner().invoke(
        seismic_cli,
        [
            "inpgen",
            str(CRYSTAL_ELASTICITY_OUTPUT),
            "--output",
            str(cli_path),
            "--interface",
            "crystal",
        ],
        input="SEISMIC lifecycle\n",
    )

    assert response.exit_code == 0, response.output
    assert cli_path.read_bytes() == api_path.read_bytes()


def test_seismic_crystal_input_generator_preserves_density(tmp_path: Path) -> None:
    """SEISMIC owns a public generator for the shared elastic input format."""
    output = seismic.create_input(
        CRYSTAL_ELASTICITY_OUTPUT,
        tmp_path / "calcite-seismic",
        interface="crystal",
        jobname="Calcite seismic API",
    )

    parsed = seismic.read_input(output)
    assert parsed.jobname == "Calcite seismic API"
    assert parsed.stiffness is not None
    assert parsed.stiffness[0, 4] == 20.67
    assert parsed.density == 2680.0


def test_seismic_and_elasticity_generators_share_input_text(tmp_path: Path) -> None:
    """Both public namespaces expose byte-identical shared input generation."""
    elastic_path = elasticity.create_input(
        CRYSTAL_ELASTICITY_OUTPUT,
        tmp_path / "elasticity-input",
        interface="crystal",
        jobname="Shared input",
    )
    seismic_path = seismic.create_input(
        CRYSTAL_ELASTICITY_OUTPUT,
        tmp_path / "seismic-input",
        interface="crystal",
        jobname="Shared input",
    )

    assert seismic_path.read_bytes() == elastic_path.read_bytes()


def test_seismic_vasp_generator_preserves_density(
    tmp_path: Path,
) -> None:
    """VASP mass, population, and final volume metadata produce kg m^-3."""
    source = _write_vasp_elasticity_output(tmp_path / "OUTCAR")

    output = seismic.create_input(
        source,
        tmp_path / "vasp-seismic",
        interface="vasp",
        jobname="VASP seismic API",
    )
    parsed = seismic.read_input(output)

    expected_density = 2.0 * 24.305 / 80.0 * 1660.53906660
    assert parsed.density == pytest.approx(expected_density, rel=5.0e-7)
    assert parsed.stiffness is not None
    assert parsed.stiffness[0, 0] == 200.0


def test_seismic_cli_vasp_input(tmp_path: Path) -> None:
    """The CLI exposes the public VASP generator with density preservation."""
    source = _write_vasp_elasticity_output(tmp_path / "OUTCAR")
    output = tmp_path / "cli-vasp.dat"
    response = CliRunner().invoke(
        seismic_cli,
        [
            "inpgen",
            str(source),
            "--output",
            str(output),
            "--interface",
            "vasp",
        ],
        input="VASP SEISMIC CLI\n",
    )

    assert response.exit_code == 0, response.output
    parsed = seismic.read_input(output)
    expected_density = 2.0 * 24.305 / 80.0 * 1660.53906660
    assert parsed.jobname == "VASP SEISMIC CLI"
    assert parsed.density == pytest.approx(expected_density, rel=5.0e-7)


def test_seismic_vasp_input_generator_rejects_missing_density(
    tmp_path: Path,
) -> None:
    """A stiffness-only VASP OUTCAR cannot become a runnable SEISMIC input."""
    source = _write_vasp_elasticity_output(
        tmp_path / "OUTCAR-no-density",
        include_density=False,
    )

    with pytest.raises(ValueError, match="finite positive density"):
        seismic.create_input(source, tmp_path / "invalid", interface="vasp")


def test_qha_input_generator_uses_shared_phonon_contract(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """QHA exposes the shared generator from its own scientific namespace."""
    captured: dict[str, object] = {}

    def fake_create(source, destination, **kwargs):
        captured.update(source=source, destination=destination, **kwargs)
        return Path(destination)

    monkeypatch.setattr(qha, "_create_phonon_input", fake_create)
    output = qha.create_input(
        "phonons.out",
        tmp_path / "qha.yaml",
        interface="crystal-qha",
        jobname="Public QHA",
        formula_units=4,
    )

    assert output == tmp_path / "qha.yaml"
    assert captured == {
        "source": "phonons.out",
        "destination": tmp_path / "qha.yaml",
        "interface": "crystal-qha",
        "is_list": False,
        "reference": 0,
        "jobname": "Public QHA",
        "formula_units": 4,
    }


def test_public_table_writers_cover_elasticity_ha_and_qha(tmp_path: Path) -> None:
    """Public table writers reproduce the three historical CLI exporters."""
    angles = np.linspace(0.0, 2.0 * np.pi, 5)
    elasticity_result = ResultData(
        metadata=ResultMetadata(module="elasticity", method="second_order"),
        results={
            "elasticity": elasticity.Result(
                jobname="Synthetic",
                crystal_system="cubic",
                properties_2d={
                    "xy": {
                        "theta": np.full(5, np.pi / 2.0),
                        "phi": angles,
                        "young_modulus": np.linspace(100.0, 110.0, 5),
                    }
                },
            )
        },
    )
    elasticity_output = elasticity.write_table(
        elasticity_result,
        tmp_path / "elasticity-table",
    )
    assert elasticity_output.suffix == ".dat"
    assert "theta_deg" in elasticity_output.read_text(encoding="utf-8")

    ha_result = ResultData(
        metadata=ResultMetadata(module="ha", method="harmonic"),
        results={
            "ha": ha.Result(
                jobname="Synthetic HA",
                temperature=np.array([0.0, 100.0]),
                volume=np.array([10.0]),
                free_energy=np.array([[1.0], [2.0]]),
                metadata={
                    "units": {
                        "temperature": "K",
                        "volume": "A^3",
                        "energy": "Ha",
                    }
                },
            )
        },
    )
    ha_output = ha.write_table(ha_result, tmp_path / "ha-table", property_name="F")
    assert ha_output.suffix == ".dat"
    assert "Helmholtz free energy" in ha_output.read_text(encoding="utf-8")

    qha_result = ResultData(
        metadata=ResultMetadata(module="qha", method="quasi-harmonic approximation"),
        results={
            "qha": qha.Result(
                jobname="Synthetic QHA",
                temperature=np.array([300.0, 600.0]),
                pressure=np.array([0.0, 5.0]),
                equilibrium_volume=np.array([[10.0, 9.8], [10.2, 10.0]]),
                uncertainties={"sigma_VT": np.full((2, 2), 0.01)},
                metadata={
                    "units": {
                        "temperature": "K",
                        "pressure": "GPa",
                        "volume": "A^3",
                    }
                },
            )
        },
    )
    qha_output = qha.write_table(
        qha_result,
        tmp_path / "qha-table",
        property_name="VT",
        file_format="csv",
    )
    assert qha_output.suffix == ".csv"
    text = qha_output.read_text(encoding="utf-8")
    assert "sigma_VT" in text
    assert "," in text


def test_cli_table_exports_match_public_api(tmp_path: Path) -> None:
    """HA and QHA CLI exports are byte-identical to their public writers."""
    ha_result = ResultData(
        metadata=ResultMetadata(module="ha", method="harmonic"),
        results={
            "ha": ha.Result(
                jobname="Synthetic HA",
                temperature=np.array([0.0, 100.0]),
                volume=np.array([10.0]),
                free_energy=np.array([[1.0], [2.0]]),
                metadata={
                    "units": {
                        "temperature": "K",
                        "volume": "A^3",
                        "energy": "Ha",
                    }
                },
            )
        },
    )
    ha_hdf5 = ha.write_result(ha_result, tmp_path / "ha-result")
    ha_api = ha.write_table(ha_result, tmp_path / "ha-api", property_name="F")
    ha_cli_path = tmp_path / "ha-cli.dat"
    response = CliRunner().invoke(
        ha_cli,
        [
            "export",
            str(ha_hdf5),
            "--output",
            str(ha_cli_path),
            "--property",
            "F",
        ],
    )
    assert response.exit_code == 0, response.output
    assert ha_cli_path.read_bytes() == ha_api.read_bytes()

    qha_result = ResultData(
        metadata=ResultMetadata(
            module="qha",
            method="quasi-harmonic approximation",
        ),
        results={
            "qha": qha.Result(
                jobname="Synthetic QHA",
                temperature=np.array([300.0, 600.0]),
                pressure=np.array([0.0, 5.0]),
                equilibrium_volume=np.array([[10.0, 9.8], [10.2, 10.0]]),
                uncertainties={"sigma_VT": np.full((2, 2), 0.01)},
                metadata={
                    "units": {
                        "temperature": "K",
                        "pressure": "GPa",
                        "volume": "A^3",
                    }
                },
            )
        },
    )
    qha_hdf5 = qha.write_result(qha_result, tmp_path / "qha-result")
    qha_api = qha.write_table(
        qha_result,
        tmp_path / "qha-api.csv",
        property_name="VT",
        file_format="csv",
    )
    qha_cli_path = tmp_path / "qha-cli.csv"
    response = CliRunner().invoke(
        qha_cli,
        [
            "export",
            str(qha_hdf5),
            "--output",
            str(qha_cli_path),
            "--property",
            "VT",
            "--format",
            "csv",
        ],
    )
    assert response.exit_code == 0, response.output
    assert qha_cli_path.read_bytes() == qha_api.read_bytes()
