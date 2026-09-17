"""Characterization tests for primitive-cell VASP Gamma phonons."""

from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from click.testing import CliRunner

from quantas.interfaces.vasp.phonons import VaspPhononReader
from quantas.models.structures import PrimitiveCellReduction, SymmetryMetadata
from quantas.cli.ha import ha
from quantas.modules.ha.io.inpgen import HAInputCreator


DATA = Path(__file__).parent / "data"
EOS_XML = DATA / "vasp_mgo_eos_00_v544.xml"

_EIGENVECTORS = (
    (0.0, -0.775525289, 0.0, 0.0, 0.631316502, 0.0),
    (0.775525289, 0.0, 0.0, -0.631316502, 0.0, 0.0),
    (0.0, 0.0, 0.775525289, 0.0, 0.0, -0.631316502),
    (-0.128832336, 0.435976292, -0.438049572, -0.158260927, 0.535564394, -0.538111264),
    (-0.609249644, -0.0144263302, 0.164824992, -0.748417798, -0.0177216718, 0.202475223),
    (-0.103815346, -0.456373838, -0.423681272, -0.127529418, -0.560621259, -0.520460878),
)


def _identity_reduction(structure, **_kwargs):
    return PrimitiveCellReduction(
        structure=structure,
        repetitions=1,
        source_to_primitive=np.eye(3),
        source_atoms=structure.natoms,
        source_volume=structure.volume,
    )


def _symmetry(_structure, **_kwargs):
    return SymmetryMetadata(
        space_group_number=225,
        international_symbol="Fm-3m",
        point_group="m-3m",
    )


def _write_gamma_run(tmp_path: Path) -> Path:
    run = tmp_path / "mgo-gamma"
    run.mkdir()
    tree = ET.parse(EOS_XML)
    root = tree.getroot()
    calculation = root.find("calculation")
    assert calculation is not None
    dynmat = ET.SubElement(calculation, "dynmat")
    eigenvalues = ET.SubElement(dynmat, "v", {"name": "eigenvalues"})
    eigenvalues.text = (
        " -5.09619014E-01 -5.09619014E-01 -5.09619014E-01"
        "  1.34026769E-06  1.34026769E-06  1.34026769E-06"
    )
    varray = ET.SubElement(dynmat, "varray", {"name": "eigenvectors"})
    for values in _EIGENVECTORS:
        row = ET.SubElement(varray, "v")
        row.text = " ".join(f"{value:.12E}" for value in values)
    tree.write(run / "vasprun.xml", encoding="utf-8", xml_declaration=True)

    (run / "OUTCAR").write_text(
        """
 Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------

   1 f  =   11.160240 THz    70.121855 2PiTHz  372.265520 cm-1    46.155059 meV
   2 f  =   11.160240 THz    70.121855 2PiTHz  372.265520 cm-1    46.155059 meV
   3 f  =   11.160240 THz    70.121855 2PiTHz  372.265520 cm-1    46.155059 meV
   4 f/i=    0.018099 THz     0.113717 2PiTHz    0.603706 cm-1     0.074850 meV
   5 f/i=    0.018099 THz     0.113717 2PiTHz    0.603706 cm-1     0.074850 meV
   6 f/i=    0.018099 THz     0.113717 2PiTHz    0.603706 cm-1     0.074850 meV

 General timing and accounting informations for this job:
""".lstrip(),
        encoding="utf-8",
    )
    return run


@pytest.mark.interfaces
def test_vasp_gamma_reader_identifies_rigid_translations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """MgO VASP-5 Gamma data become one normalized backend-neutral q-point."""
    import quantas.interfaces.vasp.phonons as module

    monkeypatch.setattr(module, "reduce_to_primitive_cell", _identity_reduction)
    monkeypatch.setattr(module, "analyze_symmetry", _symmetry)
    reader = VaspPhononReader(_write_gamma_run(tmp_path))

    assert reader.completed is True
    assert reader.error is None
    assert reader.qpoints == 1
    assert reader.natom == 2
    assert reader.nphonon == 6
    np.testing.assert_array_equal(reader.dim, np.eye(3, dtype=np.int64))
    np.testing.assert_allclose(reader.qcoords, 0.0)
    np.testing.assert_allclose(reader.weights, 1.0)
    np.testing.assert_allclose(
        reader.phonons_array()[0],
        [372.265520, 372.265520, 372.265520, 0.0, 0.0, 0.0],
    )
    assert reader.translation_indices == (3, 4, 5)
    assert np.min(reader.translation_projection_scores[3:]) > 0.99
    assert np.max(reader.translation_projection_scores[:3]) < 0.25

    modes = reader.mode_data
    assert modes.eigenvector_normalization == "mass-weighted-unit"
    assert modes.metadata["dispersion_support"] == "not_implemented"
    assert modes.metadata["raw_frequencies_cm^-1"][3:] == pytest.approx(
        [-0.603706, -0.603706, -0.603706]
    )
    norms = np.linalg.norm(modes.eigenvectors.reshape(1, 6, -1), axis=2)
    np.testing.assert_allclose(norms, 1.0, atol=1.0e-12)


@pytest.mark.interfaces
def test_vasp_gamma_reader_rejects_nonprimitive_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Folded supercell modes are not silently exposed as primitive Gamma data."""
    import quantas.interfaces.vasp.phonons as module

    def repeated(structure, **_kwargs):
        return PrimitiveCellReduction(
            structure=structure,
            repetitions=2,
            source_to_primitive=np.eye(3),
            source_atoms=structure.natoms,
            source_volume=structure.volume,
        )

    monkeypatch.setattr(module, "reduce_to_primitive_cell", repeated)
    reader = VaspPhononReader(_write_gamma_run(tmp_path))

    assert reader.completed is False
    assert reader.error is not None
    assert "dispersion" in reader.error
    assert "not implemented" in reader.error


@pytest.mark.interfaces
def test_shared_phonon_input_generator_accepts_vasp_run_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The shared HA/QHA generator consumes a VASP calculation directory."""
    import quantas.interfaces.vasp.phonons as module

    monkeypatch.setattr(module, "reduce_to_primitive_cell", _identity_reduction)
    monkeypatch.setattr(module, "analyze_symmetry", _symmetry)
    run = _write_gamma_run(tmp_path)
    creator = HAInputCreator(interface="vasp")

    ok, error = creator.read(run)
    assert ok is True, error
    data = creator.to_dict("MgO VASP Gamma", formula_units=1)

    assert data["qpoints"] == 1
    assert data["supercell"] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert data["units"]["energy"] == "eV"
    assert data["units"]["frequency"] == "cm^-1"
    assert data["phonon"][0]["q-position"] == [0.0, 0.0, 0.0]
    assert data["phonon"][0]["band"][3]["frequency"] == [0.0]
    assert data["provenance"]["interface"] == "vasp"
    assert data["provenance"]["energy"]["unit"] == "eV"


@pytest.mark.interfaces
def test_ha_inpgen_cli_accepts_vasp_run_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The shared HA/QHA Click boundary accepts a VASP calculation directory."""
    source = tmp_path / "vasp-run"
    source.mkdir()
    destination = tmp_path / "phonons.yaml"

    def fake_create_input(filename, outfile, **kwargs):
        assert Path(filename) == source
        assert kwargs["interface"] == "vasp"
        Path(outfile).write_text("job: test\n", encoding="utf-8")
        return Path(outfile)

    monkeypatch.setattr(
        "quantas.cli.phonon_input.create_ha_input",
        fake_create_input,
    )
    monkeypatch.setattr(
        "quantas.cli.phonon_input._read_generated_structure_summary",
        lambda _filename: None,
    )

    result = CliRunner().invoke(
        ha,
        [
            "inpgen",
            str(source),
            "--interface",
            "vasp",
            "--output",
            str(destination),
            "--jobname",
            "test",
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    assert destination.exists()
