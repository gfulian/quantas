"""CRYSTAL structure--energy parsing for Energy EOS input generation."""

from __future__ import annotations

from pathlib import Path

import pytest

from quantas.interfaces.crystal.energy_volume import read_crystal_energy_volume


def _cell_block(volume: float, a: float, species: tuple[int, ...]) -> str:
    atom_rows = "\n".join(
        f" {index:3d} T {number:3d} X  {index/10:.6f} 0.000000 0.000000"
        for index, number in enumerate(species, start=1)
    )
    return f"""PRIMITIVE CELL - CENTRING CODE 1/0 VOLUME= {volume:.6f}
 A B C ALPHA BETA GAMMA
 {a:.8f} {a:.8f} {a:.8f} 90.000000 90.000000 90.000000
 ATOMS IN THE UNIT CELL: {len(species)}
{atom_rows}
"""


def _static_output(
    volume: float,
    energy: float,
    *,
    corrections: str = "",
    species: tuple[int, ...] = (12, 8),
) -> str:
    a = volume ** (1.0 / 3.0)
    corrected = (
        f"TOTAL ENERGY + {corrections} (AU) {energy:.14E}\n" if corrections else ""
    )
    scf_energy = energy + 0.01 if corrections else energy
    return f"""GEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY 3
{_cell_block(volume, a, species)}
 CYC   4 ETOT(AU) {scf_energy:.14E} DETOT -1.0E-8
 == SCF ENDED - CONVERGENCE ON ENERGY E(AU) {scf_energy:.14E} CYCLES 4
{corrected} EEEEEEEEEE TERMINATION
"""


def _native_eos_output(
    states: tuple[tuple[float, float], ...],
    *,
    corrections: str = "DISP",
    species: tuple[int, ...] = (12, 8),
) -> str:
    blocks: list[str] = []
    for index, (volume, energy) in enumerate(states, start=1):
        a = volume ** (1.0 / 3.0)
        blocks.append(
            f"""ATOM ONLY OPTIMIZATION - POINT {index}
 * OPT END - CONVERGED * E(AU): {energy - 0.02:.14E} POINTS {index} *
 FINAL OPTIMIZED GEOMETRY - DIMENSIONALITY OF THE SYSTEM 3
 {_cell_block(volume, a, species)}
 CYC   5 ETOT(AU) {energy + 0.01:.14E} DETOT -1.0E-9
 == SCF ENDED - CONVERGENCE ON ENERGY E(AU) {energy + 0.01:.14E} CYCLES 5
 TOTAL ENERGY + {corrections} (AU) {energy:.14E}
"""
        )
    summary = "\n".join(f" {volume:.6f} {energy:.12E}" for volume, energy in sorted(states))
    return "\n".join(blocks) + f"""
 SORTING VOLUMES/ENERGIES
 VOLUME (A^3) ENERGY (a.u.)
{summary}
 +++++++ FITTING USING ALL POINTS +++++++
 EEEEEEEEEE TERMINATION
"""


def test_static_output_returns_one_state(tmp_path: Path) -> None:
    source = tmp_path / "static.out"
    source.write_text(_static_output(18.8, -275.17), encoding="utf-8")

    result = read_crystal_energy_volume(source)

    assert result.run_kind == "static"
    assert result.corrections == ()
    assert result.series.npoints == 1
    assert result.series.energies[0] == pytest.approx(-275.17)


def test_optimization_uses_post_geometry_total(tmp_path: Path) -> None:
    volume = 11.0
    energy = -10.25
    a = volume ** (1.0 / 3.0)
    text = f"""ATOM ONLY OPTIMIZATION - POINT 1
 * OPT END - CONVERGED * E(AU): -10.00 POINTS 1 *
 FINAL OPTIMIZED GEOMETRY - DIMENSIONALITY OF THE SYSTEM 3
 {_cell_block(volume, a, (12, 8))}
 CYC 4 ETOT(AU) -10.20 DETOT -1.0E-8
 == SCF ENDED - CONVERGENCE ON ENERGY E(AU) -10.20 CYCLES 4
 TOTAL ENERGY + DISP (AU) {energy:.14E}
 EEEEEEEEEE TERMINATION
"""
    source = tmp_path / "opt.out"
    source.write_text(text, encoding="utf-8")

    result = read_crystal_energy_volume(source)

    assert result.run_kind == "optimization"
    assert result.corrections == ("DISP",)
    assert result.series.npoints == 1
    assert result.series.energies[0] == pytest.approx(energy)


def test_native_eos_matches_sorted_states(tmp_path: Path) -> None:
    source = tmp_path / "eos.out"
    source.write_text(
        _native_eos_output(((12.0, -10.2), (10.0, -10.0))),
        encoding="utf-8",
    )

    result = read_crystal_energy_volume(source)

    assert result.run_kind == "native_eos"
    assert result.corrections == ("DISP",)
    assert result.series.npoints == 2
    assert [point.metadata["reported_volume_angstrom3"] for point in result.series.points] == pytest.approx([10.0, 12.0])
    assert result.series.energies.tolist() == pytest.approx([-10.0, -10.2])


def test_native_eos_rejects_summary_mismatch(tmp_path: Path) -> None:
    text = _native_eos_output(((10.0, -10.0), (12.0, -10.2)))
    text = text.replace("10.000000 -1.000000000000E+01", "10.000000 -9.900000000000E+00")
    source = tmp_path / "bad.out"
    source.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="summary energy disagrees"):
        read_crystal_energy_volume(source)
