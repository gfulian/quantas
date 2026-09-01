from __future__ import annotations

from pathlib import Path

import numpy as np

from quantas.models.phonons import PhononModeData
from quantas.modules.ha.io.inpgen import HAInputCreator, create_ha_input
from quantas.modules.ha.io.reader import HAInputFileReader


class FakePhononReader:
    """Small completed reader exposing the historical phonon attributes."""

    completed = True
    error = None
    natom = 2
    dim = np.identity(3, dtype=int)
    qpoints = 2
    qcoords = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    weights = np.array([1, 3], dtype=int)
    shrinkf = np.array([2, 2, 2], dtype=int)
    nphonon = 6
    units = {
        "energy": "Ha",
        "volume": "angstrom^3",
        "frequency": "cm^-1",
        "length": "angstrom",
    }

    def __init__(self, filename: str | Path) -> None:
        path = Path(filename)
        value = float(path.read_text(encoding="utf-8").strip())
        self.volume = 10.0 + value
        self.energy = -100.0 - value
        self.phonons = np.zeros((self.qpoints, self.nphonon), dtype=float)
        for iq in range(self.qpoints):
            for imode in range(self.nphonon):
                self.phonons[iq, imode] = 100.0 + 10.0 * iq + imode + value

    def phonons_array(self) -> np.ndarray:
        """Return frequencies using the interface-reader array contract."""
        return np.asarray(self.phonons, dtype=np.float64)


class TrackingFakePhononReader(FakePhononReader):
    """Fake reader exposing deliberately permuted phonon eigenvectors."""

    qpoints = 1
    qcoords = np.array([[0, 0, 0]], dtype=float)
    weights = np.array([1], dtype=int)
    nphonon = 3
    natom = 1

    def __init__(self, filename: str | Path) -> None:
        path = Path(filename)
        value = int(float(path.read_text(encoding="utf-8").strip()))
        self.volume = 10.0 + value
        self.energy = -100.0 - value
        basis = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
        if value == 0:
            self.phonons = np.array([[100.0, 200.0, 300.0]])
            vectors = basis
        else:
            self.phonons = np.array([[200.0 + value, 300.0 + value, 100.0 + value]])
            vectors = basis[[1, 2, 0]]
        self.mode_data = PhononModeData(
            frequencies=self.phonons.copy(),
            eigenvectors=vectors[np.newaxis, ...],
            atom_symbols=("H",),
        )


class FakeCrystalQHAReader:
    """Small completed reader exposing the historical CRYSTAL-QHA attributes."""

    completed = True
    error = None
    points = 3
    natom = 1
    dim = np.diag([2, 2, 2])
    qpoints = 1
    qcoords = np.array([[0, 0, 0]], dtype=float)
    weights = np.array([8], dtype=int)
    shrinkf = np.ones(3, dtype=int)
    nphonon = 3
    units = {
        "energy": "Ha",
        "volume": "angstrom^3",
        "frequency": "cm^-1",
        "length": "angstrom",
    }
    volume = np.array([9.0, 10.0, 11.0], dtype=float)
    energy = np.array([-10.0, -11.0, -10.5], dtype=float)
    phonons = np.array(
        [[[100.0, 101.0, 102.0], [200.0, 201.0, 202.0], [300.0, 301.0, 302.0]]]
    )
    mode_continuity = "verified"
    mode_continuity_metadata = {"method": "crystal-qha", "source": "crystal"}

    def __init__(self, filename: str | Path) -> None:
        self.filename = Path(filename)


def test_create_ha_input_from_crystal_qha_single_file(tmp_path):
    source = tmp_path / "qha.out"
    source.write_text("fake", encoding="utf-8")
    outfile = tmp_path / "qha.yaml"

    create_ha_input(
        source,
        outfile,
        interface="crystal-qha",
        jobname="Single QHA",
        interface_filter={"crystal-qha": FakeCrystalQHAReader},
    )

    reader = HAInputFileReader(outfile)
    assert reader.completed
    ha_input = reader.to_input()
    assert ha_input.jobname == "Single QHA"
    assert ha_input.natoms == 1
    assert ha_input.qpoints == 1
    np.testing.assert_allclose(ha_input.volume, [9.0, 10.0, 11.0])
    assert ha_input.frequencies.shape == (1, 3, 3)
    np.testing.assert_allclose(ha_input.frequencies[0, 1], [200.0, 201.0, 202.0])


def test_creator_rejects_crystal_qha_file_list(tmp_path):
    file_list = tmp_path / "files.lst"
    file_list.write_text("one.out\n", encoding="utf-8")
    creator = HAInputCreator(
        interface="crystal-qha",
        interface_filter={"crystal-qha": FakeCrystalQHAReader},
    )

    completed, error = creator.read(file_list, is_list=True)

    assert completed is False
    assert "Only one CRYSTAL output file" in error


def test_creator_rejects_invalid_reference(tmp_path):
    f0 = tmp_path / "v0.out"
    f0.write_text("0", encoding="utf-8")
    creator = HAInputCreator(
        interface="crystal",
        interface_filter={"crystal": FakePhononReader},
    )

    completed, error = creator.read(f0, is_list=False, reference=2)

    assert completed is False
    assert error == "Invalid reference provided"


def test_creator_reports_missing_file(tmp_path):
    missing = tmp_path / "missing.out"
    creator = HAInputCreator(
        interface="crystal",
        interface_filter={"crystal": FakePhononReader},
    )

    completed, error = creator.read(missing)

    assert completed is False
    assert "does not exists" in error


def test_input_generator_writes_formula_units_explicitly(tmp_path):
    source = tmp_path / "qha.out"
    source.write_text("fake", encoding="utf-8")
    outfile = tmp_path / "qha.yaml"

    create_ha_input(
        source,
        outfile,
        interface="crystal-qha",
        jobname="Z test",
        formula_units=4,
        interface_filter={"crystal-qha": FakeCrystalQHAReader},
    )

    text = outfile.read_text(encoding="utf-8")
    input_data = HAInputFileReader(outfile).to_input()

    assert "formula_units: 4\n" in text
    assert "mode_continuity: verified\n" in text
    assert input_data.formula_units == 4


def test_q_positions_are_fractional(tmp_path):
    first = tmp_path / "v0.out"
    second = tmp_path / "v1.out"
    first.write_text("0", encoding="utf-8")
    second.write_text("1", encoding="utf-8")
    file_list = tmp_path / "phonons.lst"
    file_list.write_text("v0.out\nv1.out\n", encoding="utf-8")
    outfile = tmp_path / "phonons.yaml"

    create_ha_input(
        file_list,
        outfile,
        interface="crystal",
        is_list=True,
        jobname="fractional q points",
        interface_filter={"crystal": FakePhononReader},
    )

    reader = HAInputFileReader(outfile)
    assert reader.completed, reader.error
    input_data = reader.to_input()
    assert input_data.qcoords is not None
    np.testing.assert_allclose(
        input_data.qcoords,
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        rtol=0.0,
        atol=0.0,
    )
    assert input_data.metadata["q_position_source"] == "reader-provided"


def test_file_list_entries_are_resolved_relative_to_the_list_file(
    tmp_path, monkeypatch
):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "v0.out").write_text("0", encoding="utf-8")
    file_list = data_dir / "phonons.lst"
    file_list.write_text("v0.out\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    creator = HAInputCreator(
        interface="crystal",
        interface_filter={"crystal": FakePhononReader},
    )
    completed, error = creator.read(file_list, is_list=True)

    assert completed is True, error
    assert creator.files == [data_dir / "v0.out"]


def test_multiple_file_generator_rejects_q_point_order_changes(tmp_path):
    first = tmp_path / "v0.out"
    second = tmp_path / "v1.out"
    first.write_text("0", encoding="utf-8")
    second.write_text("1", encoding="utf-8")
    creator = HAInputCreator(
        interface="crystal",
        interface_filter={"crystal": FakePhononReader},
    )
    completed, error = creator.read(first)
    assert completed is True, error
    creator.phondata.append(FakePhononReader(second))
    creator.phondata[1].qcoords = np.array(
        [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=float,
    )

    import pytest

    with pytest.raises(ValueError, match="different q-point coordinates"):
        creator.to_dict(jobname="bad q mesh")


def test_input_generator_writes_units_and_source_provenance(tmp_path):
    source = tmp_path / "qha.out"
    source.write_text("fake", encoding="utf-8")
    outfile = tmp_path / "qha.yaml"

    create_ha_input(
        source,
        outfile,
        interface="crystal-qha",
        jobname="unit metadata",
        interface_filter={"crystal-qha": FakeCrystalQHAReader},
    )

    reader = HAInputFileReader(outfile)
    input_data = reader.to_input()

    assert reader.units == FakeCrystalQHAReader.units
    assert input_data.units == FakeCrystalQHAReader.units
    assert input_data.metadata["provenance"]["interface"] == "crystal-qha"
    assert input_data.metadata["provenance"]["reference_index"] == 0
    assert input_data.metadata["provenance"]["sources"] == [str(source)]


def test_multiple_file_generator_marks_missing_eigenvectors_unknown(tmp_path):
    first = tmp_path / "v0.out"
    second = tmp_path / "v1.out"
    first.write_text("0", encoding="utf-8")
    second.write_text("1", encoding="utf-8")
    file_list = tmp_path / "phonons.lst"
    file_list.write_text("v0.out\nv1.out\n", encoding="utf-8")

    creator = HAInputCreator(
        interface="crystal",
        interface_filter={"crystal": FakePhononReader},
    )
    completed, error = creator.read(file_list, is_list=True)
    assert completed is True, error

    data = creator.to_dict(jobname="unknown continuity")

    assert data["mode_continuity"] == "unknown"
    assert data["mode_continuity_metadata"]["reason"] == (
        "phonon_eigenvectors_unavailable"
    )


def test_multiple_file_generator_reorders_verified_mode_branches(tmp_path):
    first = tmp_path / "v0.out"
    second = tmp_path / "v1.out"
    first.write_text("0", encoding="utf-8")
    second.write_text("1", encoding="utf-8")
    file_list = tmp_path / "phonons.lst"
    file_list.write_text("v0.out\nv1.out\n", encoding="utf-8")

    creator = HAInputCreator(
        interface="crystal",
        interface_filter={"crystal": TrackingFakePhononReader},
    )
    completed, error = creator.read(file_list, is_list=True)
    assert completed is True, error

    data = creator.to_dict(jobname="tracked continuity", reference=0)

    assert data["mode_continuity"] == "verified"
    assert data["mode_continuity_metadata"]["method"] == "eigenvector_overlap"
    frequencies = [band["frequency"] for band in data["phonon"][0]["band"]]
    np.testing.assert_allclose(
        frequencies,
        [[100.0, 101.0], [200.0, 201.0], [300.0, 301.0]],
    )
