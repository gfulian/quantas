from pathlib import Path

import numpy as np

from quantas.modules.elasticity.io.reader import ElasticityInputFileReader


DATA = Path(__file__).parent / "data" / "hydroxylapatite.dat"


def test_elasticity_input_reader_reads_hydroxylapatite():
    reader = ElasticityInputFileReader()
    reader.load(DATA)

    assert reader.completed is True
    assert reader.error is None
    assert reader.jobname == "Hydroxylapatite"
    assert not hasattr(reader, "density")
    assert reader.stiffness.shape == (6, 6)


def test_elasticity_input_reader_reconstructs_upper_triangular_matrix():
    reader = ElasticityInputFileReader()
    reader.load(DATA)

    stiffness = reader.stiffness

    assert np.allclose(stiffness, stiffness.T)
    np.testing.assert_allclose(stiffness[0, 0], 187.208)
    np.testing.assert_allclose(stiffness[0, 1], 65.193)
    np.testing.assert_allclose(stiffness[2, 2], 222.658)
    np.testing.assert_allclose(stiffness[5, 5], 61.007)


def test_elasticity_input_reader_reconstructs_lower_triangular_matrix(tmp_path):
    """Conventional lower-triangular input is assembled without reordering."""
    expected = np.array(
        [
            [10.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            [1.0, 20.0, 6.0, 7.0, 8.0, 9.0],
            [2.0, 6.0, 30.0, 10.0, 11.0, 12.0],
            [3.0, 7.0, 10.0, 40.0, 13.0, 14.0],
            [4.0, 8.0, 11.0, 13.0, 50.0, 15.0],
            [5.0, 9.0, 12.0, 14.0, 15.0, 60.0],
        ]
    )
    rows = [
        " ".join(str(value) for value in expected[index, : index + 1])
        for index in range(6)
    ]
    filename = tmp_path / "lower.dat"
    filename.write_text(
        "Lower triangular\n" + "\n".join(rows) + "\n3000\n", encoding="utf-8"
    )

    reader = ElasticityInputFileReader(filename)

    assert reader.completed is True
    np.testing.assert_allclose(reader.stiffness, expected)
    assert not hasattr(reader, "density")
