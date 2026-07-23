from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.geometry import (
    analyze_symmetry,
    find_primitive_structure,
    reconstruct_primitive_structure,
)
from quantas.interfaces.crystal.geometry import CrystalGeometryParser
from quantas.models import CrystalStructure


pytestmark = [pytest.mark.interfaces, pytest.mark.fast]

DATA = Path(__file__).parent / "data" / "dolomite_crystal23_geometry_no_coorprt.out"


def _diagonal_supercell(structure: CrystalStructure, size: int) -> CrystalStructure:
    """Build an exact diagonal supercell without relying on spglib helpers."""
    expansion = np.diag([size, size, size])
    positions = []
    numbers = []
    for shift in np.ndindex(size, size, size):
        shift_array = np.asarray(shift, dtype=np.float64)
        for position, atomic_number in zip(
            structure.fractional_positions,
            structure.atomic_numbers,
            strict=True,
        ):
            positions.append((position + shift_array) / float(size))
            numbers.append(int(atomic_number))
    return CrystalStructure(
        lattice=expansion @ structure.lattice,
        fractional_positions=np.asarray(positions, dtype=np.float64),
        atomic_numbers=np.asarray(numbers, dtype=np.int64),
        label=f"{size}x{size}x{size} exact supercell",
    )


def test_geometry_parser_reads_primitive_without_coorprt():
    parser = CrystalGeometryParser(DATA)

    structure = parser.initial_primitive_cell()
    symmetry = analyze_symmetry(structure)

    assert parser.has_coorprt is False
    assert structure.natoms == 10
    assert structure.volume == pytest.approx(106.790803, abs=1.0e-5)
    np.testing.assert_array_equal(parser.expansion_matrix(), np.diag([3, 3, 3]))
    np.testing.assert_allclose(
        parser.primitive_to_crystallographic(),
        [[1.0, 0.0, 1.0], [-1.0, 1.0, 1.0], [0.0, -1.0, 1.0]],
    )
    assert parser.space_group_number() == 148
    assert symmetry.space_group_number == 148
    assert symmetry.international_symbol == "R-3"


@pytest.mark.parametrize("size", [2, 3, 4, 5])
def test_spglib_and_explicit_reduction_are_independent_of_supercell_parity(size):
    reference = CrystalGeometryParser(DATA).initial_primitive_cell()
    source = _diagonal_supercell(reference, size)
    expansion = np.diag([size, size, size])

    primitive_spglib = find_primitive_structure(source, symprec=1.0e-7)
    primitive_explicit, diagnostics = reconstruct_primitive_structure(
        source,
        expansion,
        reference,
        exact_tolerance=1.0e-7,
    )

    assert primitive_spglib.natoms == reference.natoms == 10
    assert primitive_explicit.natoms == reference.natoms
    assert diagnostics.status == "exact"
    assert diagnostics.minimum_replica_count == size**3
    assert diagnostics.maximum_replica_count == size**3
    assert diagnostics.maximum_translation_residual < 1.0e-10


def test_translationally_perturbed_supercell_is_averaged_not_hidden_by_symprec():
    reference = CrystalGeometryParser(DATA).initial_primitive_cell()
    size = 2
    source = _diagonal_supercell(reference, size)
    expansion = np.diag([size, size, size])

    rng = np.random.default_rng(1729)
    primitive_lattice = reference.lattice
    cartesian_noise = rng.normal(scale=1.5e-3, size=source.fractional_positions.shape)
    primitive_noise = cartesian_noise @ np.linalg.inv(primitive_lattice)
    source.fractional_positions = source.fractional_positions + primitive_noise / float(
        size
    )

    primitive, diagnostics = reconstruct_primitive_structure(
        source,
        expansion,
        reference,
    )

    assert primitive.natoms == 10
    assert diagnostics.status == "averaged"
    assert 0.0 < diagnostics.maximum_translation_residual < 0.05
    assert diagnostics.minimum_replica_count == diagnostics.maximum_replica_count == 8
