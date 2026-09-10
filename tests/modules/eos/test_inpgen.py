"""Energy EOS input generation and mixed-source collection."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from quantas.api import eos
from quantas.cli.main import main


def _cell(volume: float, species: tuple[int, ...]) -> str:
    a = volume ** (1.0 / 3.0)
    rows = "\n".join(
        f" {i} T {z} X {i/10:.6f} 0.0 0.0"
        for i, z in enumerate(species, start=1)
    )
    return f"""PRIMITIVE CELL - CENTRING CODE 1/0 VOLUME= {volume:.6f}
 A B C ALPHA BETA GAMMA
 {a:.8f} {a:.8f} {a:.8f} 90.0 90.0 90.0
 ATOMS IN THE UNIT CELL: {len(species)}
{rows}
"""


def _static(volume: float, energy: float, species: tuple[int, ...]) -> str:
    return f"""GEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY 3
{_cell(volume, species)}
 CYC 3 ETOT(AU) {energy + 0.01:.14E} DETOT -1.0E-8
 == SCF ENDED - CONVERGENCE ON ENERGY E(AU) {energy + 0.01:.14E} CYCLES 3
 TOTAL ENERGY + DISP (AU) {energy:.14E}
 EEEEEEEEEE TERMINATION
"""


def _native(states: tuple[tuple[float, float], ...], species: tuple[int, ...]) -> str:
    blocks: list[str] = []
    for index, (volume, energy) in enumerate(states, start=1):
        blocks.append(
            f"""ATOM ONLY OPTIMIZATION - POINT {index}
 * OPT END - CONVERGED * E(AU): {energy - 0.02:.14E} POINTS {index} *
 FINAL OPTIMIZED GEOMETRY - DIMENSIONALITY OF THE SYSTEM 3
 {_cell(volume, species)}
 CYC 5 ETOT(AU) {energy + 0.01:.14E} DETOT -1.0E-9
 == SCF ENDED - CONVERGENCE ON ENERGY E(AU) {energy + 0.01:.14E} CYCLES 5
 TOTAL ENERGY + DISP (AU) {energy:.14E}
"""
        )
    summary = "\n".join(
        f" {volume:.6f} {energy:.12E}" for volume, energy in sorted(states)
    )
    return "\n".join(blocks) + f"""
 SORTING VOLUMES/ENERGIES
 VOLUME (A^3) ENERGY (a.u.)
{summary}
 +++++++ FITTING USING ALL POINTS +++++++
 EEEEEEEEEE TERMINATION
"""


def test_list_flattens_single_and_native_eos(tmp_path: Path) -> None:
    """One list may combine a single state with one native CRYSTAL EOS series."""
    native = tmp_path / "native.out"
    extra = tmp_path / "extra.out"
    native.write_text(_native(((10.0, -10.0), (12.0, -10.2)), (12, 8)), encoding="utf-8")
    extra.write_text(_static(14.0, -9.8, (12, 8)), encoding="utf-8")
    listing = tmp_path / "files.txt"
    listing.write_text("native.out\nextra.out\n", encoding="utf-8")

    written = eos.create_input(listing, tmp_path / "combined", is_list=True)
    dataset = eos.read_input(written)

    assert written.suffix == ".dat"
    assert dataset.npoints == 3
    assert dataset.column("volume").tolist() == pytest.approx([10.0, 12.0, 14.0])
    assert dataset.column("energy").tolist() == pytest.approx([-10.0, -10.2, -9.8])
    for name in ("a", "b", "c", "alpha", "beta", "gamma"):
        assert name in dataset.columns


def test_list_rejects_mixed_composition(tmp_path: Path) -> None:
    native = tmp_path / "native.out"
    other = tmp_path / "other.out"
    native.write_text(_native(((10.0, -10.0), (12.0, -10.2)), (12, 8)), encoding="utf-8")
    other.write_text(_static(14.0, -20.0, (6, 8, 8)), encoding="utf-8")
    listing = tmp_path / "files.txt"
    listing.write_text("native.out\nother.out\n", encoding="utf-8")

    with pytest.raises(ValueError, match="atom counts|chemical compositions"):
        eos.create_input(listing, tmp_path / "mixed.dat", is_list=True)


def test_list_rejects_mixed_corrections(tmp_path: Path) -> None:
    native = tmp_path / "native.out"
    plain = tmp_path / "plain.out"
    native.write_text(_native(((10.0, -10.0), (12.0, -10.2)), (12, 8)), encoding="utf-8")
    plain_text = _static(14.0, -9.8, (12, 8)).replace(
        " TOTAL ENERGY + DISP (AU) -9.80000000000000E+00\n",
        "",
    ).replace(
        "-9.79000000000000E+00",
        "-9.80000000000000E+00",
    )
    plain.write_text(plain_text, encoding="utf-8")
    listing = tmp_path / "files.txt"
    listing.write_text("native.out\nplain.out\n", encoding="utf-8")

    with pytest.raises(ValueError, match="incompatible energy corrections"):
        eos.create_input(listing, tmp_path / "mixed.dat", is_list=True)


def test_cli_inpgen_writes_energy_dataset(tmp_path: Path) -> None:
    source = tmp_path / "single.out"
    source.write_text(_static(18.8, -275.17, (12, 8)), encoding="utf-8")
    target = tmp_path / "energy.dat"

    result = CliRunner().invoke(
        main,
        ["eos", "inpgen", str(source), "-o", str(target)],
    )

    assert result.exit_code == 0, result.output
    assert target.is_file()
    assert "Energy EOS input summary" in result.output
    dataset = eos.read_input(target)
    assert dataset.npoints == 1
    assert dataset.column("volume")[0] == pytest.approx(18.8)
