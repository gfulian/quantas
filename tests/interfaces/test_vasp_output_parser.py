"""Characterization tests for the generic VASP run-output interface."""

from __future__ import annotations

from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from quantas.interfaces.vasp import (
    VaspOutputParser,
    VaspRunDocument,
    resolve_vasp_run_source,
)
from quantas.models.computation import RunTerminationStatus


DATA = Path(__file__).resolve().parent / "data"
OPT_XML = DATA / "vasp_mgo_opt_step01_v544.xml"
OPT_OUTCAR = DATA / "vasp_mgo_opt_step01_v544.OUTCAR"
EOS_XML = DATA / "vasp_mgo_eos_00_v544.xml"
EOS_OUTCAR = DATA / "vasp_mgo_eos_00_v544.OUTCAR"


def _run_directory(
    tmp_path: Path,
    xml_source: Path = OPT_XML,
    outcar_source: Path | None = OPT_OUTCAR,
) -> Path:
    """Create one compact VASP run directory from characterized fixtures."""
    run = tmp_path / "vasp_run"
    run.mkdir()
    shutil.copyfile(xml_source, run / "vasprun.xml")
    if outcar_source is not None:
        shutil.copyfile(outcar_source, run / "OUTCAR")
    return run


@pytest.mark.interfaces
def test_vasp_run_source_resolves_directory_and_primary_files(tmp_path: Path) -> None:
    """A calculation directory should be the canonical VASP source unit."""
    run = _run_directory(tmp_path)

    source = resolve_vasp_run_source(run)

    assert source.directory == run
    assert source.vasprun_xml == run / "vasprun.xml"
    assert source.outcar == run / "OUTCAR"
    assert resolve_vasp_run_source(run / "vasprun.xml") == source
    assert resolve_vasp_run_source(run / "OUTCAR") == source


@pytest.mark.interfaces
def test_vasp_run_source_requires_vasprun_xml(tmp_path: Path) -> None:
    """OUTCAR-only ingestion is not silently treated as a complete run document."""
    run = tmp_path / "vasp_run"
    run.mkdir()
    (run / "OUTCAR").write_text("VASP output\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="vasprun.xml"):
        resolve_vasp_run_source(run)


@pytest.mark.interfaces
def test_vasp_document_characterizes_real_v544_metadata(tmp_path: Path) -> None:
    """Generator, INCAR, and atom order should survive structured parsing."""
    run = _run_directory(tmp_path)
    document = VaspRunDocument(run)

    assert document.generator()["program"] == "vasp"
    assert document.generator()["version"] == "5.4.4.18Apr17-6-g9f103f2a35"
    assert document.parameter("IBRION") == 2
    assert document.parameter("ISIF") == 3
    assert document.parameter("ISMEAR") == 1
    assert document.parameter("NSW") == 40
    assert document.parameter("NELM") == 60
    assert document.atom_symbols() == ("O", "Mg")
    np.testing.assert_array_equal(document.atomic_numbers(), [8, 12])
    assert len(document.ionic_step_nodes()) == 6


@pytest.mark.interfaces
def test_vasp_parser_reconstructs_real_mgo_structures(tmp_path: Path) -> None:
    """VASP fractional coordinates must remain float64 and preserve atom order."""
    run = _run_directory(tmp_path)
    parser = VaspOutputParser(run)

    initial = parser.initial_structure()
    final = parser.final_structure()

    assert initial.lattice.dtype == np.float64
    assert initial.fractional_positions.dtype == np.float64
    assert initial.atomic_numbers.dtype == np.int64
    np.testing.assert_array_equal(initial.atomic_numbers, [8, 12])
    np.testing.assert_allclose(
        initial.fractional_positions,
        [[0.5, 0.5, 0.5], [0.0, 0.0, 0.0]],
    )
    assert initial.metadata["reported_volume_angstrom3"] == pytest.approx(
        19.29438437
    )
    assert final.metadata["reported_volume_angstrom3"] == pytest.approx(19.26345154)
    assert initial.volume == pytest.approx(19.2943843915)
    assert final.volume == pytest.approx(19.2634515263)


@pytest.mark.interfaces
def test_vasp_parser_resolves_v544_outer_energy_bug(tmp_path: Path) -> None:
    """VASP 5.4.4 outer energy tags must not corrupt E or E0."""
    run = _run_directory(tmp_path)
    parser = VaspOutputParser(run)

    steps = parser.ionic_steps()
    last = steps[-1]

    assert len(steps) == 6
    assert last.electronic_steps == 4
    assert last.energies.free_energy.value == pytest.approx(-11.87757624)
    assert last.energies.energy_without_entropy.value == pytest.approx(-11.87945457)
    assert last.energies.sigma_zero_energy.value == pytest.approx(-11.87820235)
    assert last.energies.entropy_term == pytest.approx(0.00187834)
    assert last.energies.metadata["resolution"] == (
        "last_scstep_plus_outer_free_energy_shift"
    )
    assert last.energies.metadata["outer_energy_tags_consistent"] is False
    assert last.energies.metadata["known_vasp5_outer_energy_tag_bug"] is True
    assert last.metadata["outcar_energy_crosscheck"] == "matched"


@pytest.mark.interfaces
def test_vasp_parser_characterizes_forces_stress_and_convergence(tmp_path: Path) -> None:
    """Generic ionic states should retain mechanical data without unit conversion."""
    run = _run_directory(tmp_path)
    parser = VaspOutputParser(run)

    steps = parser.ionic_steps()

    assert parser.termination().status is RunTerminationStatus.NORMAL
    assert parser.ionic_converged() is True
    assert all(step.metadata["electronic_converged"] is True for step in steps)
    assert all(step.metadata["ionic_converged"] is True for step in steps)
    assert steps[-1].forces is not None
    assert steps[-1].forces.shape == (2, 3)
    assert steps[-1].stress is not None
    assert steps[-1].stress.shape == (3, 3)
    assert steps[-1].forces.dtype == np.float64
    assert steps[-1].stress.dtype == np.float64


@pytest.mark.interfaces
def test_vasp_parser_characterizes_fixed_cell_mgo_state(tmp_path: Path) -> None:
    """The MgO EOS reference run should parse as one converged ionic state."""
    run = _run_directory(tmp_path, EOS_XML, EOS_OUTCAR)
    parser = VaspOutputParser(run)

    steps = parser.ionic_steps()

    assert parser.document.parameter("ISIF") == 2
    assert parser.document.parameter("ISMEAR") == -5
    assert len(steps) == 1
    assert steps[0].structure.volume == pytest.approx(19.2819239986)
    assert steps[0].energies.free_energy.value == pytest.approx(-11.87213470)
    assert steps[0].energies.energy_without_entropy.value == pytest.approx(
        -11.87213470
    )
    assert steps[0].energies.sigma_zero_energy.value == pytest.approx(-11.87213470)
    assert steps[0].metadata["outcar_energy_crosscheck"] == "matched"


@pytest.mark.interfaces
def test_vasp_parser_can_use_xml_without_outcar_crosscheck(tmp_path: Path) -> None:
    """OUTCAR is complementary rather than mandatory for structured ingestion."""
    run = _run_directory(tmp_path, outcar_source=None)
    parser = VaspOutputParser(run)

    assert parser.source.outcar is None
    assert parser.termination().status is RunTerminationStatus.NORMAL
    assert parser.ionic_converged() is None
    assert len(parser.ionic_steps()) == 6


@pytest.mark.interfaces
def test_vasp_parser_rejects_disagreement_between_xml_and_outcar(tmp_path: Path) -> None:
    """Paired primary outputs must not disagree silently on energy."""
    run = _run_directory(tmp_path)
    outcar = (run / "OUTCAR").read_text(encoding="utf-8")
    outcar = outcar.replace("-11.87757624", "-10.87757624", 1)
    (run / "OUTCAR").write_text(outcar, encoding="utf-8")

    with pytest.raises(ValueError, match="disagrees between vasprun.xml"):
        VaspOutputParser(run).ionic_steps()


@pytest.mark.interfaces
def test_vasp_parser_accepts_flat_ionic_layout(tmp_path: Path) -> None:
    """The parser should also understand the flat ionic layout documented by VASP."""
    tree = ET.parse(EOS_XML)
    root = tree.getroot()
    calculation = root.find("calculation")
    assert calculation is not None
    scsteps = calculation.findall("scstep")
    assert scsteps
    last_energy = scsteps[-1].find("energy")
    assert last_energy is not None

    insertion = list(root).index(calculation)
    root.remove(calculation)
    flat_nodes = [
        calculation.find("structure"),
        calculation.find("varray[@name='forces']"),
        calculation.find("varray[@name='stress']"),
        last_energy,
    ]
    for offset, node in enumerate(flat_nodes):
        assert node is not None
        root.insert(insertion + offset, node)

    run = tmp_path / "vasp_run"
    run.mkdir()
    tree.write(run / "vasprun.xml", encoding="utf-8", xml_declaration=True)

    step = VaspOutputParser(run).ionic_steps()[0]

    assert step.metadata["source_layout"] == "flat"
    assert step.electronic_steps == 0
    assert step.energies.metadata["resolution"] == "outer_ionic_energy"
    assert step.energies.sigma_zero_energy.value == pytest.approx(-11.87213470)


@pytest.mark.interfaces
def test_vasp_energy_resolution_transfers_outer_shift(
    tmp_path: Path,
) -> None:
    """An outer additive correction should shift F, E, and E0 equally."""
    tree = ET.parse(EOS_XML)
    calculation = tree.getroot().find("calculation")
    assert calculation is not None
    outer = calculation.find("energy/i[@name='e_fr_energy']")
    assert outer is not None and outer.text is not None
    outer.text = f"{float(outer.text) - 0.125:.8f}"

    run = tmp_path / "vasp_run"
    run.mkdir()
    tree.write(run / "vasprun.xml", encoding="utf-8", xml_declaration=True)

    energies = VaspOutputParser(run).ionic_steps()[0].energies

    assert energies.free_energy.value == pytest.approx(-11.99713470)
    assert energies.energy_without_entropy.value == pytest.approx(-11.99713470)
    assert energies.sigma_zero_energy.value == pytest.approx(-11.99713470)
    assert energies.metadata["outer_free_energy_shift_eV"] == pytest.approx(-0.125)


@pytest.mark.interfaces
def test_vasp_document_rejects_malformed_xml(tmp_path: Path) -> None:
    """A truncated or malformed vasprun.xml must fail explicitly."""
    run = tmp_path / "vasp_run"
    run.mkdir()
    (run / "vasprun.xml").write_text("<modeling><generator>", encoding="utf-8")

    with pytest.raises(ValueError, match="malformed VASP XML document"):
        VaspRunDocument(run)
