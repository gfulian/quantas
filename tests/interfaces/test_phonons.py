from __future__ import annotations

import numpy as np
import pytest

from quantas.interfaces.crystal.phonons import CrystalPhononReader


QPOINT_BLOCK = """
 *******************************************************************************
 * LIST OF THE K POINTS USED FOF PHONON DISPERSION.                            *
 *  K       WEIGHT       COORD                                                 *
 *      1         1.        0  0  0                                            *
 *      2         2.        0  0  1                                            *
 *      3         3.        1  2  2                                            *
 * WITH SHRINKING FACTORS: IS1 =     2 IS2 =     3 IS3 =     4                 *
 *******************************************************************************
  DISPERSION K POINT NUMBER     1 COORD:  R(  0  0  0 )    WEIGHT:    1.
  DISPERSION K POINT NUMBER     2 COORD:  C(  0  0  1 )    WEIGHT:    2.
  DISPERSION K POINT NUMBER     3 COORD:  C(  1  2  2 )    WEIGHT:    3.
"""


def test_crystal_qpoint_table_preserves_order_weights_and_fractional_positions(
    tmp_path,
) -> None:
    filename = tmp_path / "qmesh.out"
    filename.write_text(QPOINT_BLOCK, encoding="utf-8")
    reader = CrystalPhononReader()

    qpoints, coordinates, weights, shrinkf = reader.set_q_mesh(filename)
    reader.qpoints = qpoints
    reader.qcoords = coordinates
    reader.weights = weights
    reader.shrinkf = shrinkf

    assert qpoints == 3
    np.testing.assert_allclose(weights, [1.0, 2.0, 3.0])
    np.testing.assert_allclose(shrinkf, [2.0, 3.0, 4.0])
    np.testing.assert_allclose(
        reader.qcoords_fractional,
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.25], [0.5, 2.0 / 3.0, 0.5]],
    )
    assert reader.q_position_source == "crystal-dispersion-table"


def test_crystal_qpoint_parser_rejects_table_dispersion_mismatch(tmp_path) -> None:
    filename = tmp_path / "bad_qmesh.out"
    filename.write_text(
        QPOINT_BLOCK.replace(
            "NUMBER     3 COORD:  C(  1  2  2 )",
            "NUMBER     3 COORD:  C(  1  1  2 )",
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="different coordinates or ordering"):
        CrystalPhononReader().set_q_mesh(filename)
