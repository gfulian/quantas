"""Characterization tests for CRYSTAL phonon displacement eigenvectors."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from quantas.core.chemistry import atomic_mass
from quantas.core.numerics.phonon_tracking import track_phonon_modes
from quantas.interfaces.crystal.phonon_modes import CrystalPhononModeParser
from quantas.interfaces.crystal.qha import CrystalQHAReader


ROOT = Path(__file__).resolve().parents[2]
PHONON_OUTPUT = (
    ROOT / "examples/qha/crystal-phonons/dol_pbe0_crystal_pdisp_p01.out"
)


def _constant_mode_block(frequencies: list[float], natoms: int) -> list[str]:
    """Build one compact CRYSTAL-like eigenvector block for testing."""
    values = " ".join(f"{0.01 * (index + 1):9.4f}" for index in range(len(frequencies)))
    lines = [
        " FREQ(CM**-1) " + " ".join(f"{value:9.2f}" for value in frequencies),
        "",
    ]
    for atom in range(1, natoms + 1):
        lines.append(f" AT. {atom:3d} H  X {values}")
        lines.append(f"            Y {values}")
        lines.append(f"            Z {values}")
    return lines


def test_crystal_mode_parser_handles_partial_final_block() -> None:
    """The final printed block may contain fewer than six phonon modes."""
    lines = [" NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES (IN BOHR)"]
    lines.extend(_constant_mode_block([1, 2, 3, 4, 5, 6], natoms=3))
    lines.extend(_constant_mode_block([7, 8, 9], natoms=3))

    data = CrystalPhononModeParser(lines).parse(nmodes=9)

    assert data is not None
    assert data.frequencies.shape == (1, 9)
    assert data.eigenvectors.shape == (1, 9, 3, 3)
    np.testing.assert_allclose(data.frequencies[0], np.arange(1.0, 10.0))
    np.testing.assert_allclose(
        np.linalg.norm(data.eigenvectors.reshape(1, 9, -1), axis=2),
        np.ones((1, 9)),
    )


def test_crystal_mode_parser_reads_real_dispersion_vectors() -> None:
    """Real and complex CRYSTAL modes should share one neutral representation."""
    data = CrystalPhononModeParser(PHONON_OUTPUT).parse(nmodes=30)

    assert data is not None
    assert data.frequencies.shape == (27, 30)
    assert data.eigenvectors.shape == (27, 30, 10, 3)
    assert data.atom_symbols[:4] == ("Ca", "Mg", "C", "C")
    assert data.frequency_unit == "cm^-1"
    assert data.eigenvector_normalization == "mass-weighted-unit"
    np.testing.assert_allclose(
        data.frequencies[0, :6],
        [0, 0, 0, 146.58, 157.03, 157.03],
    )
    np.testing.assert_allclose(
        np.linalg.norm(data.eigenvectors.reshape(27, 30, -1), axis=2),
        np.ones((27, 30)),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    assert np.allclose(data.eigenvectors[0].imag, 0.0)
    assert np.max(np.abs(data.eigenvectors[1].imag)) > 0.0

QHA_OUTPUT = ROOT / "examples/qha/crystal-qha/mgo-b3lyp-crystal-qha.out"


def test_crystal_mode_parser_reads_native_qha_mode_sets() -> None:
    """A monolithic CRYSTAL QHA output may contain many Gamma mode sets."""
    data = CrystalPhononModeParser(QHA_OUTPUT).parse(nmodes=192)

    assert data is not None
    assert data.frequencies.shape == (11, 192)
    assert data.eigenvectors.shape == (11, 192, 64, 3)
    np.testing.assert_allclose(
        np.linalg.norm(data.eigenvectors.reshape(11, 192, -1), axis=2),
        np.ones((11, 192)),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    translation_ratio = abs(
        data.eigenvectors[0, 0, 0, 0] / data.eigenvectors[0, 0, 1, 0]
    )
    expected_ratio = np.sqrt(atomic_mass("Mg") / atomic_mass("O"))
    np.testing.assert_allclose(translation_ratio, expected_ratio, rtol=1.0e-12)


def test_native_crystal_qha_reports_verified_mode_continuity() -> None:
    """The native QHA output should expose CRYSTAL's continuity result."""
    reader = CrystalQHAReader()

    assert reader.has_verified_mode_continuity(QHA_OUTPUT)


def test_native_qha_mode_sets_are_trackable_without_ambiguity() -> None:
    """Quantas overlaps should reproduce a safe native-QHA tracking case."""
    data = CrystalPhononModeParser(QHA_OUTPUT).parse(nmodes=192)
    assert data is not None
    frequencies = data.frequencies.reshape(11, 1, 192)
    eigenvectors = data.eigenvectors.reshape(11, 1, 192, 64, 3)

    result = track_phonon_modes(
        frequencies,
        eigenvectors,
        np.arange(11, dtype=float),
        reference_index=5,
    )

    assert result.verified
    assert result.ambiguous_assignments == 0
    assert result.low_overlap_assignments == 0
    assert np.isnan(result.minimum_overlap)
    assert result.minimum_subspace_singular_value > 0.9
