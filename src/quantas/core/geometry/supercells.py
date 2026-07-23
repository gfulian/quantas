# -*- coding: utf-8 -*-

"""Supercell reduction and translational-reconstruction utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.geometry.cells import (
    minimum_image_delta,
    wrap_fractional,
)
from quantas.models.structures import (
    CrystalStructure,
    StructureReconstructionDiagnostics,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def supercell_repetitions(expansion_matrix: ArrayLike) -> int:
    """Return the number of primitive repetitions in a supercell.

    Parameters
    ----------
    expansion_matrix : array_like
        Integer expansion matrix with shape ``(3, 3)``.

    Returns
    -------
    int
        Absolute rounded determinant of the expansion matrix.

    Raises
    ------
    ValueError
        If the matrix is invalid or singular.
    """
    matrix = np.asarray(expansion_matrix, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise ValueError("expansion_matrix must have shape (3, 3)")
    determinant = float(np.linalg.det(matrix))
    repetitions = int(round(abs(determinant)))
    if repetitions <= 0 or not np.isclose(abs(determinant), repetitions, atol=1.0e-8):
        raise ValueError("expansion_matrix must have a non-zero integer determinant")
    return repetitions


def primitive_lattice_from_supercell(
    supercell_lattice: ArrayLike,
    expansion_matrix: ArrayLike,
) -> FloatArray:
    """Recover the primitive lattice in the source Cartesian orientation.

    CRYSTAL defines supercell vectors by ``B = D A`` when direct vectors are
    stored by rows, where ``D`` is the expansion matrix, ``A`` the primitive
    lattice, and ``B`` the supercell lattice.

    Parameters
    ----------
    supercell_lattice : array_like
        Supercell direct lattice vectors with shape ``(3, 3)``.
    expansion_matrix : array_like
        Integer primitive-to-supercell expansion matrix.

    Returns
    -------
    ndarray
        Primitive direct lattice matrix with shape ``(3, 3)``.

    Raises
    ------
    ValueError
        If matrices do not have shape ``(3, 3)`` or are singular.
    """
    supercell = np.asarray(supercell_lattice, dtype=np.float64)
    expansion = np.asarray(expansion_matrix, dtype=np.float64)
    if supercell.shape != (3, 3) or expansion.shape != (3, 3):
        raise ValueError("lattice and expansion_matrix must have shape (3, 3)")
    supercell_repetitions(expansion)
    return np.asarray(np.linalg.solve(expansion, supercell), dtype=np.float64)


def fold_supercell_positions(
    fractional_positions: ArrayLike,
    expansion_matrix: ArrayLike,
) -> FloatArray:
    """Fold supercell fractional coordinates into the primitive cell.

    Parameters
    ----------
    fractional_positions : array_like
        Supercell fractional coordinates with shape ``(natoms, 3)``.
    expansion_matrix : array_like
        Primitive-to-supercell expansion matrix.

    Returns
    -------
    ndarray
        Primitive fractional coordinates wrapped into ``[0, 1)``.
    """
    positions = np.asarray(fractional_positions, dtype=np.float64)
    expansion = np.asarray(expansion_matrix, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("fractional_positions must have shape (natoms, 3)")
    if expansion.shape != (3, 3):
        raise ValueError("expansion_matrix must have shape (3, 3)")
    return wrap_fractional(positions @ expansion)


def reconstruct_primitive_structure(
    source: CrystalStructure,
    expansion_matrix: ArrayLike,
    reference: CrystalStructure,
    *,
    exact_tolerance: float = 1.0e-6,
    acceptance_tolerance: float = 5.0e-2,
) -> tuple[CrystalStructure, StructureReconstructionDiagnostics]:
    """Reconstruct a compact primitive structure from a source supercell.

    Source atoms are folded through the known expansion matrix, assigned to
    reference primitive atoms of the same species using minimum-image
    Cartesian distances, and translational copies are averaged around the
    reference positions. This deliberately avoids increasing a free symmetry
    tolerance until a desired primitive atom count is obtained.

    Parameters
    ----------
    source : CrystalStructure
        Source supercell printed by the electronic-structure code.
    expansion_matrix : array_like
        Primitive-to-supercell expansion matrix.
    reference : CrystalStructure
        Reference primitive structure that fixes atom correspondence and basis
        orientation.
    exact_tolerance : float, optional
        Maximum residual in angstrom below which reconstruction is labelled
        ``"exact"``.
    acceptance_tolerance : float, optional
        Maximum residual in angstrom accepted as a translationally averaged
        structure.

    Returns
    -------
    CrystalStructure
        Reconstructed primitive structure in the source Cartesian orientation.
    StructureReconstructionDiagnostics
        Replica counts and residual displacements.

    Raises
    ------
    ValueError
        If species or replica counts are inconsistent, or residuals exceed the
        acceptance tolerance.
    """
    expansion = np.asarray(expansion_matrix, dtype=np.float64)
    repetitions = supercell_repetitions(expansion)
    primitive_lattice = primitive_lattice_from_supercell(source.lattice, expansion)
    folded = fold_supercell_positions(source.fractional_positions, expansion)
    reference_positions = wrap_fractional(reference.fractional_positions)

    groups: list[list[int]] = [[] for _ in range(reference.natoms)]
    for source_index, (position, atomic_number) in enumerate(
        zip(folded, source.atomic_numbers, strict=True)
    ):
        candidates = np.flatnonzero(reference.atomic_numbers == atomic_number)
        if candidates.size == 0:
            raise ValueError(
                f"source atom {source_index} has atomic number {atomic_number} "
                "not present in the reference primitive cell"
            )
        deltas = minimum_image_delta(position - reference_positions[candidates])
        distances = np.linalg.norm(deltas @ primitive_lattice, axis=1)
        groups[int(candidates[int(np.argmin(distances))])].append(source_index)

    counts = np.asarray([len(group) for group in groups], dtype=np.int64)
    if np.any(counts != repetitions):
        raise ValueError(
            "translational reconstruction produced inconsistent replica counts: "
            f"expected {repetitions}, observed {counts.tolist()}"
        )

    averaged = np.zeros_like(reference_positions)
    residuals: list[float] = []
    for reference_index, group in enumerate(groups):
        group_positions = folded[np.asarray(group, dtype=np.int64)]
        deltas = minimum_image_delta(
            group_positions - reference_positions[reference_index]
        )
        mean_delta = np.mean(deltas, axis=0)
        averaged[reference_index] = wrap_fractional(
            reference_positions[reference_index] + mean_delta
        )
        group_residuals = minimum_image_delta(deltas - mean_delta)
        residuals.extend(
            np.linalg.norm(group_residuals @ primitive_lattice, axis=1).tolist()
        )

    residual_array = np.asarray(residuals, dtype=np.float64)
    maximum = float(np.max(residual_array)) if residual_array.size else 0.0
    rms = (
        float(np.sqrt(np.mean(residual_array * residual_array)))
        if residual_array.size
        else 0.0
    )
    if maximum <= exact_tolerance:
        status = "exact"
        message = "Translational copies are equivalent within numerical tolerance."
    elif maximum <= acceptance_tolerance:
        status = "averaged"
        message = (
            "Translational copies were averaged because the optimized source "
            "supercell weakly breaks primitive translations."
        )
    else:
        status = "inconsistent"
        message = (
            "Translational residual exceeds the accepted reconstruction "
            f"tolerance of {acceptance_tolerance:.6g} angstrom."
        )
        raise ValueError(message)

    reconstructed = CrystalStructure(
        lattice=primitive_lattice,
        fractional_positions=averaged,
        atomic_numbers=reference.atomic_numbers.copy(),
        label="reconstructed primitive",
        metadata={
            "source_label": source.label,
            "reference_label": reference.label,
            "reconstruction_status": status,
        },
    )
    diagnostics = StructureReconstructionDiagnostics(
        status=status,
        source_atoms=source.natoms,
        reconstructed_atoms=reconstructed.natoms,
        expected_repetitions=repetitions,
        minimum_replica_count=int(counts.min()) if counts.size else 0,
        maximum_replica_count=int(counts.max()) if counts.size else 0,
        maximum_translation_residual=maximum,
        rms_translation_residual=rms,
        message=message,
    )
    return reconstructed, diagnostics
