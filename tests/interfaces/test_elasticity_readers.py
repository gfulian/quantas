"""Tests for external-code elasticity interfaces."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.interfaces.crystal import markers
from quantas.interfaces.crystal.elasticity import CrystalElasticityReader
from quantas.interfaces.vasp.elasticity import VASPElasticityReader


DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize(
    "option_marker",
    markers.ELASTICITY_OPTION_MARKERS,
    ids=("elastcon", "elapiezo"),
)
def test_crystal_reader_collects_tensor_and_density(
    tmp_path,
    option_marker: str,
) -> None:
    """CRYSTAL elastic outputs are normalized to a symmetric GPa matrix."""
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
            option_marker,
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


def test_crystal_reader_uses_unstrained_elastic_reference(tmp_path: Path) -> None:
    """Geometry, energy, and stress come from the state before elastic strains."""

    def geometry_block(marker: str, a: float, density: float) -> list[str]:
        volume = a**3
        return [
            marker,
            "LATTICE PARAMETERS (ANGSTROMS AND DEGREES)",
            (
                "PRIMITIVE CELL - CENTRING CODE 1/0 "
                f"VOLUME= {volume:.10f} - DENSITY {density:.6f} g/cm^3"
            ),
            "        A              B              C           ALPHA      BETA       GAMMA",
            f" {a:.12f} {a:.12f} {a:.12f} 90.000000 90.000000 90.000000",
            "ATOMS IN THE ASYMMETRIC UNIT    1 - ATOMS IN THE UNIT CELL:    1",
            "     ATOM                 X/A                 Y/B                 Z/C",
            "      1 T   8 O     0.000000000000E+00 0.000000000000E+00 0.000000000000E+00",
            "DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)",
            "          X                    Y                    Z",
            f" {a:.12f} 0.000000000000 0.000000000000",
            f" 0.000000000000 {a:.12f} 0.000000000000",
            f" 0.000000000000 0.000000000000 {a:.12f}",
        ]

    lines = ["ELAPIEZO OPTION", "COORPRT"]
    lines.extend(geometry_block("GEOMETRY FOR WAVE FUNCTION - TEST", 4.2, 3.0))
    lines.extend(
        [
            "PRESSURE IN GIGAPASCAL: 4.00000000",
            "TOTAL ENERGY(DFT)(AU)(  3) -19.500000000000",
        ]
    )
    lines.extend(geometry_block("FINAL OPTIMIZED GEOMETRY - TEST", 4.0, 3.2))
    lines.extend(
        [
            "TOTAL ENERGY(DFT)(AU)(  2) -20.000000000000",
            "PRESSURE IN GIGAPASCAL: 1.50000000",
            "VOLUME OF THE CELL: 64.0000000000",
            "DENSITY OF THE CRYSTAL = 3.20000000 g/cm^3",
            "STRAIN MATRIX   1 :",
        ]
    )
    lines.extend(geometry_block("FINAL OPTIMIZED GEOMETRY - TEST", 5.0, 4.0))
    lines.extend(
        [
            "TOTAL ENERGY(DFT)(AU)(  2) -19.000000000000",
            "PRESSURE IN GIGAPASCAL: 9.00000000",
            "DENSITY OF THE CRYSTAL = 4.00000000 g/cm^3",
            "FINAL RESULTS START",
            "SYMMETRIZED ELASTIC CONSTANTS",
            "header",
            "200 80 70 0 0 0",
            "190 65 0 0 0",
            "180 0 0 0",
            "60 0 0",
            "55 0",
            "50",
        ]
    )
    output = tmp_path / "elastic-coorprt.out"
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")

    reader = CrystalElasticityReader(output)

    assert reader.completed is True
    assert reader.error is None
    assert reader.structure is not None
    assert reader.structure.volume == pytest.approx(64.0)
    assert reader.volume == pytest.approx(64.0)
    assert reader.density == pytest.approx(3200.0)
    assert reader.energy == pytest.approx(-20.0)
    assert reader.stress_pressure == pytest.approx(1.5)
