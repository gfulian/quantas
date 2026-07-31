"""Tests for external-code elasticity interfaces."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.interfaces.crystal.elasticity import CrystalElasticityReader
from quantas.interfaces.vasp.elasticity import VASPElasticityReader


DATA = Path(__file__).parent / "data"


def test_crystal_reader_collects_tensor_and_density(
    tmp_path,
) -> None:
    """CRYSTAL ELASTCON output is normalized to a symmetric GPa matrix."""
    rows = [
        "| 100 10 20 0 0 0 |",
        "| 110 30 0 0 0 |",
        "| 120 0 0 0 |",
        "| 40 0 0 |",
        "| 50 0 |",
        "| 60 |",
    ]
    text = "\n".join(
        [
            "ELASTCON OPTION",
            "GEOMETRY NOW FULLY CONSISTENT WITH THE GROUP",
            "PRIMITIVE CELL - TEST 3.178 g/cm3",
            "FINAL RESULTS START",
            "SYMMETRIZED ELASTIC CONSTANTS",
            "header",
            *rows,
        ]
    )
    filename = tmp_path / "crystal.out"
    filename.write_text(text + "\n", encoding="utf-8")

    reader = CrystalElasticityReader(filename)

    assert reader.completed is True
    assert reader.error is None
    assert reader.density == 3178.0
    assert reader.stiffness[0, 0] == 100.0
    assert reader.stiffness[0, 2] == reader.stiffness[2, 0] == 20.0
    assert reader.stiffness[5, 5] == 60.0


def test_vasp_reader_prefers_relaxed_moduli(
    tmp_path,
) -> None:
    """VASP relaxed-ion elastic moduli are selected and converted to GPa."""
    clamped = np.diag([1000.0, 1100.0, 1200.0, 400.0, 500.0, 600.0])
    relaxed = np.diag([2000.0, 2100.0, 2200.0, 700.0, 800.0, 900.0])

    def block(header: str, matrix: np.ndarray) -> list[str]:
        lines = [header, "separator", "separator"]
        lines.extend(
            f"{row + 1} " + " ".join(str(value) for value in matrix[row])
            for row in range(6)
        )
        return lines

    filename = tmp_path / "OUTCAR"
    filename.write_text(
        "\n".join(
            [
                "POMASS = 24.305; ZVAL = 2.000",
                "ions per type = 2",
                "volume of cell : 80.000000",
                *block("SYMMETRIZED ELASTIC MODULI (kBar)", clamped),
                *block("TOTAL ELASTIC MODULI (kBar)", relaxed),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    reader = VASPElasticityReader(filename)

    assert reader.completed is True
    assert reader.error is None
    expected = relaxed / 10.0
    expected[[3, 5]] = expected[[5, 3]]
    expected[:, [3, 5]] = expected[:, [5, 3]]
    np.testing.assert_allclose(reader.stiffness, expected)
    expected_density = 2.0 * 24.305 / 80.0 * 1660.53906660
    assert reader.density == pytest.approx(expected_density)


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_crystal_calcite_reader_preserves_reported_source_components() -> None:
    """The CRYSTAL interface parses the source-frame matrix without rotation."""
    reader = CrystalElasticityReader(DATA / "calcite_crystal_elastcon_excerpt.out")

    assert reader.completed is True
    assert reader.error is None
    assert reader.density == 2680.0
    assert reader.stiffness[0, 3] == pytest.approx(0.0)
    assert reader.stiffness[0, 4] == pytest.approx(20.670)
    assert reader.stiffness[3, 5] == pytest.approx(-20.670)
    assert reader.stiffness[4, 5] == pytest.approx(0.0)
