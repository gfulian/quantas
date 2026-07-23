from __future__ import annotations

import numpy as np

from quantas.io.phonons import PhononInputFileReader
from quantas.modules.ha.io.reader import HAInputFileReader


YAML_INPUT = """
job: HA YAML test
natom: 2
supercell:
  - [2, 0, 0]
  - [0, 2, 0]
  - [0, 0, 1]
qpoints: 2
volume:
  - "10.0 11.0"
energy:
  - "-20.0 -19.5"
phonon:
  - weight: 1
    band:
      - frequency:
          - "1.0 1.1"
      - frequency:
          - "2.0 2.1"
      - frequency:
          - "3.0 3.1"
      - frequency:
          - "4.0 4.1"
      - frequency:
          - "5.0 5.1"
      - frequency:
          - "6.0 6.1"
  - weight: 3
    band:
      - frequency:
          - "7.0 7.1"
      - frequency:
          - "8.0 8.1"
      - frequency:
          - "9.0 9.1"
      - frequency:
          - "10.0 10.1"
      - frequency:
          - "11.0 11.1"
      - frequency:
          - "12.0 12.1"
"""


def test_yaml_reader_preserves_shared_phonon_contract(tmp_path):
    input_file = tmp_path / "ha_input.yaml"
    input_file.write_text(YAML_INPUT, encoding="utf-8")

    reader = HAInputFileReader(input_file)

    assert reader.completed is True
    assert reader.error is None
    assert reader.jobname == "HA YAML test"
    assert reader.natoms == 2
    assert reader.kpoints == 4
    assert reader.qpoints == 2
    assert reader.nvol == 2
    assert reader.total_q_points == 4.0

    np.testing.assert_allclose(reader.volume, np.array([10.0, 11.0]))
    np.testing.assert_allclose(reader.energy, np.array([-20.0, -19.5]))
    np.testing.assert_allclose(reader.weights, np.array([1.0, 3.0]))

    frequencies = reader.frequencies
    assert frequencies.shape == (2, 6, 2)
    np.testing.assert_allclose(frequencies[0, 0], np.array([1.0, 1.1]))
    np.testing.assert_allclose(frequencies[1, 5], np.array([12.0, 12.1]))


def test_ha_reader_uses_shared_phonon_reader(tmp_path):
    input_file = tmp_path / "ha_input.yaml"
    input_file.write_text(YAML_INPUT, encoding="utf-8")

    reader = HAInputFileReader(input_file)

    assert isinstance(reader, PhononInputFileReader)
    assert reader.completed is True
    assert reader.qpoints == 2


def test_reader_reports_missing_required_fields(tmp_path):
    input_file = tmp_path / "broken.yaml"
    input_file.write_text("job: broken\nvolume: 10.0\n", encoding="utf-8")

    reader = HAInputFileReader(input_file)

    assert reader.completed is False
    assert reader.error == "No energy values in input file"


def test_reader_defaults_formula_units_to_one(tmp_path):
    input_file = tmp_path / "ha_input.yaml"
    input_file.write_text(YAML_INPUT, encoding="utf-8")

    reader = HAInputFileReader(input_file)
    input_data = reader.to_input()

    assert reader.formula_units == 1
    assert input_data.formula_units == 1
    assert input_data.natoms_per_formula_unit == 2.0


def test_reader_accepts_explicit_formula_units(tmp_path):
    input_file = tmp_path / "ha_input.yaml"
    input_file.write_text(
        YAML_INPUT.replace("natom: 2", "natom: 2\nformula_units: 2"),
        encoding="utf-8",
    )

    input_data = HAInputFileReader(input_file).to_input()

    assert input_data.formula_units == 2
    assert input_data.natoms_per_formula_unit == 1.0
