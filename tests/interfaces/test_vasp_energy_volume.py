"""VASP Energy EOS adapter characterization and normalization tests."""

from __future__ import annotations

from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

from click.testing import CliRunner
import numpy as np
import pytest

from quantas.api import eos
from quantas.cli.main import main
from quantas.core.physics.units import convert_energy
from quantas.interfaces.vasp.energy_volume import (
    VaspEnergyVolumeReader,
    read_vasp_energy_volume,
)
from quantas.models.structures import (
    PrimitiveCellReduction,
    SymmetryMetadata,
)
from quantas.modules.eos.io.inpgen import EOSEnergyInputCreator


DATA = Path(__file__).resolve().parent / "data"
EOS_XML = DATA / "vasp_mgo_eos_00_v544.xml"
EOS_OUTCAR = DATA / "vasp_mgo_eos_00_v544.OUTCAR"
OPT_XML = DATA / "vasp_mgo_opt_step01_v544.xml"
OPT_OUTCAR = DATA / "vasp_mgo_opt_step01_v544.OUTCAR"


def _run_directory(
    tmp_path: Path,
    *,
    name: str,
    xml_source: Path = EOS_XML,
    outcar_source: Path | None = EOS_OUTCAR,
) -> Path:
    """Create one compact characterized VASP run directory."""
    run = tmp_path / name
    run.mkdir()
    shutil.copyfile(xml_source, run / "vasprun.xml")
    if outcar_source is not None:
        shutil.copyfile(outcar_source, run / "OUTCAR")
    return run


def _identity_reduction(structure, **kwargs):
    """Return a deterministic primitive identity reduction for parser tests."""
    del kwargs
    return PrimitiveCellReduction(
        structure=structure,
        repetitions=1,
        source_to_primitive=np.eye(3),
        source_atoms=structure.natoms,
        source_volume=structure.volume,
    )


def _mgo_symmetry(structure, **kwargs):
    """Return MgO Fm-3m metadata without requiring spglib in this test layer."""
    del structure, kwargs
    return SymmetryMetadata(
        space_group_number=225,
        international_symbol="Fm-3m",
        transformation_matrix=np.diag([0.5, 0.5, 1.0]),
    )


def _patch_symmetry(monkeypatch: pytest.MonkeyPatch) -> None:
    """Patch only the spglib boundary while retaining real VASP fixtures."""
    monkeypatch.setattr(
        "quantas.interfaces.vasp.energy_volume.reduce_to_primitive_cell",
        _identity_reduction,
    )
    monkeypatch.setattr(
        "quantas.interfaces.vasp.energy_volume.analyze_symmetry",
        _mgo_symmetry,
    )
    monkeypatch.setattr(
        "quantas.modules.eos.io.inpgen.analyze_symmetry",
        _mgo_symmetry,
    )


def _modified_single_state_run(
    tmp_path: Path,
    *,
    name: str,
    volume: float,
    energy: float,
    ismear: int = -5,
    kmesh: tuple[int, int, int] = (12, 12, 12),
) -> Path:
    """Create a second compact single-state run from the real MgO fixture."""
    tree = ET.parse(EOS_XML)
    root = tree.getroot()
    reference_volume = 19.2819239986
    scale = (float(volume) / reference_volume) ** (1.0 / 3.0)
    for structure in root.findall("structure") + root.findall("calculation/structure"):
        basis = structure.find("crystal/varray[@name='basis']")
        assert basis is not None
        for vector in basis.findall("v"):
            values = [float(item) * scale for item in (vector.text or "").split()]
            vector.text = " " + " ".join(f"{value:.12f}" for value in values) + " "
        volume_node = structure.find("crystal/i[@name='volume']")
        assert volume_node is not None
        volume_node.text = f" {volume:.12f} "
    for item in root.findall(".//i[@name='ISMEAR']"):
        item.text = f" {ismear} "
    divisions = root.find("kpoints/generation/v[@name='divisions']")
    assert divisions is not None
    divisions.text = " " + " ".join(str(value) for value in kmesh) + " "
    for item in root.findall(".//energy/i"):
        if item.get("name") in {"e_fr_energy", "e_wo_entrp", "e_0_energy"}:
            item.text = f" {energy:.12f} "
        elif item.get("name") == "eentropy":
            item.text = " 0.000000000000 "

    run = tmp_path / name
    run.mkdir()
    tree.write(run / "vasprun.xml", encoding="utf-8", xml_declaration=True)
    return run


@pytest.mark.interfaces
def test_vasp_eos_selects_sigma_zero_energy_and_primitive_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The EOS adapter should select E0 while preserving explicit semantics."""
    _patch_symmetry(monkeypatch)
    tree = ET.parse(EOS_XML)
    root = tree.getroot()
    for item in root.findall(".//i[@name='ISMEAR']"):
        item.text = " 0 "
    for item in root.findall(".//i[@name='SIGMA']"):
        item.text = " 0.1 "
    for energy_node in root.findall(".//energy"):
        values = {
            "e_fr_energy": -11.70,
            "e_wo_entrp": -11.90,
            "e_0_energy": -11.80,
            "eentropy": 0.20,
        }
        for name, value in values.items():
            item = energy_node.find(f"i[@name='{name}']")
            if item is not None:
                item.text = f" {value:.8f} "
    run = tmp_path / "eos00"
    run.mkdir()
    tree.write(run / "vasprun.xml", encoding="utf-8", xml_declaration=True)

    result = read_vasp_energy_volume(run)
    point = result.series.points[0]

    assert result.run_kind == "single_state"
    assert result.source_repetitions == 1
    np.testing.assert_allclose(result.source_to_primitive, np.eye(3))
    assert point.energy.value == pytest.approx(-11.80)
    assert point.energy.unit == "eV"
    assert point.energy.metadata["energy_selection"] == "sigma_zero_energy"
    assert point.energy.metadata["vasp_tag"] == "e_0_energy"
    assert point.volume == pytest.approx(19.2819239986)
    assert "quantity=e_0_energy" in result.energy_signature
    assert "VASP_VERSION=5.4.4.18Apr17-6-g9f103f2a35" in result.energy_signature
    assert "ISMEAR=0" in result.energy_signature
    assert "SIGMA=0.1" in result.energy_signature
    assert "KPOINTS:mode=Gamma" in result.energy_signature
    assert "KPOINTS:divisions=12,12,12" in result.energy_signature
    assert any(item.startswith("POTENTIALS=PAW_PBE O") for item in result.energy_signature)


@pytest.mark.interfaces
def test_vasp_energy_volume_rejects_optimization_history(
    tmp_path: Path,
) -> None:
    """A multi-step optimization must not be silently flattened into one EOS point."""
    run = _run_directory(
        tmp_path,
        name="optimization",
        xml_source=OPT_XML,
        outcar_source=OPT_OUTCAR,
    )

    with pytest.raises(ValueError, match="exactly one ionic state"):
        VaspEnergyVolumeReader(run).read()


@pytest.mark.interfaces
def test_vasp_eos_rejects_fixed_occupancy_energy_surface(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """ISMEAR=-2 constrained occupations are outside the ground-state EOS policy."""
    _patch_symmetry(monkeypatch)
    run = _modified_single_state_run(
        tmp_path, name="fixed_occ", volume=19.0, energy=-11.9, ismear=-2
    )

    with pytest.raises(ValueError, match="ISMEAR=-2"):
        VaspEnergyVolumeReader(run).read()


@pytest.mark.interfaces
def test_vasp_eos_collection_rejects_mixed_smearing_semantics(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One E(V) dataset cannot mix independently defined smearing surfaces."""
    _patch_symmetry(monkeypatch)
    run_a = _modified_single_state_run(
        tmp_path, name="tetra", volume=19.0, energy=-11.9, ismear=-5
    )
    run_b = _modified_single_state_run(
        tmp_path, name="tetra_alt", volume=20.0, energy=-11.8, ismear=-4
    )

    with pytest.raises(ValueError, match="incompatible energy semantics/settings"):
        EOSEnergyInputCreator(interface="vasp").read([run_a, run_b])


@pytest.mark.interfaces
def test_vasp_eos_collection_rejects_mixed_kpoint_sampling(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One E(V) dataset should not silently mix Brillouin-zone samplings."""
    _patch_symmetry(monkeypatch)
    run_a = _modified_single_state_run(
        tmp_path, name="mesh12", volume=19.0, energy=-11.9, kmesh=(12, 12, 12)
    )
    run_b = _modified_single_state_run(
        tmp_path, name="mesh10", volume=20.0, energy=-11.8, kmesh=(10, 10, 10)
    )

    with pytest.raises(ValueError, match="incompatible energy semantics/settings"):
        EOSEnergyInputCreator(interface="vasp").read([run_a, run_b])


@pytest.mark.interfaces
def test_vasp_eos_cli_accepts_directory_list(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The Click adapter should expose the same VASP directory-list workflow."""
    _patch_symmetry(monkeypatch)
    _modified_single_state_run(
        tmp_path, name="state01", volume=18.0, energy=-11.80
    )
    _modified_single_state_run(
        tmp_path, name="state02", volume=20.0, energy=-11.90
    )
    listing = tmp_path / "vasp-runs.txt"
    listing.write_text("state01\nstate02\n", encoding="utf-8")
    output = tmp_path / "mgo-vasp-eos.dat"

    result = CliRunner().invoke(
        main,
        [
            "eos",
            "inpgen",
            str(listing),
            "--interface",
            "vasp",
            "--list",
            "-o",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output.is_file()
    assert "Interface" in result.output
    assert "vasp" in result.output
    assert "States" in result.output


@pytest.mark.interfaces
def test_vasp_eos_input_generation_accepts_directory_list(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The shared EOS generator should consume a list of VASP run directories."""
    _patch_symmetry(monkeypatch)
    run_a = _modified_single_state_run(
        tmp_path, name="state01", volume=18.0, energy=-11.80
    )
    run_b = _modified_single_state_run(
        tmp_path, name="state02", volume=20.0, energy=-11.90
    )
    listing = tmp_path / "vasp-runs.txt"
    listing.write_text("state01\nstate02\n", encoding="utf-8")

    written = eos.create_input(
        listing,
        tmp_path / "mgo-vasp-eos",
        interface="vasp",
        is_list=True,
    )
    dataset = eos.read_input(written)

    assert dataset.npoints == 2
    assert dataset.units["energy"] == "Ha"
    assert dataset.column("volume").tolist() == pytest.approx([18.0, 20.0])
    np.testing.assert_allclose(
        dataset.column("energy"),
        convert_energy(np.array([-11.80, -11.90]), "eV", "Ha"),
    )
    assert dataset.metadata["crystal_reference"] == "primitive"
    assert dataset.metadata["space_group_number"] == 225
