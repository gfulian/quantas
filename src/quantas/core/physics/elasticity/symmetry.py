# -*- coding: utf-8 -*-

"""Detect elastic symmetry and specialize second-order stiffness tensors.

The routines compare a Voigt stiffness matrix with the independent-component
patterns allowed by each crystallographic system.
"""

from __future__ import annotations

import itertools as it
from dataclasses import dataclass
from typing import TYPE_CHECKING, Mapping, Sequence

import numpy as np

if TYPE_CHECKING:
    from .tensor import ElasticTensor


ELASTIC_SYMMETRY_PATTERNS = {
    "cubic": np.array(
        [
            [1, 7, 7, 0, 0, 0],
            [7, 1, 7, 0, 0, 0],
            [7, 7, 1, 0, 0, 0],
            [0, 0, 0, 4, 0, 0],
            [0, 0, 0, 0, 4, 0],
            [0, 0, 0, 0, 0, 4],
        ]
    ),
    "hexagonal": np.array(
        [
            [1, 7, 8, 9, 10, 0],
            [7, 1, 8, 0, -9, 0],
            [8, 8, 3, 0, 0, 0],
            [9, -9, 0, 4, 0, 0],
            [10, 0, 0, 0, 4, 0],
            [0, 0, 0, 0, 0, 6],
        ]
    ),
    "trigonal_high": np.array(
        [
            [1, 7, 8, 9, 10, 0],
            [7, 1, 8, 0, -9, 0],
            [8, 8, 3, 0, 0, 0],
            [9, -9, 0, 4, 0, 0],
            [10, 0, 0, 0, 4, 0],
            [0, 0, 0, 0, 0, 6],
        ]
    ),
    "trigonal_low": np.array(
        [
            [1, 7, 8, 9, 10, 0],
            [7, 1, 8, -9, -10, 0],
            [8, 8, 3, 0, 0, 0],
            [9, -9, 0, 4, 0, -10],
            [10, -10, 0, 0, 4, 9],
            [0, 0, 0, -10, 9, 6],
        ]
    ),
    "tetragonal_high": np.array(
        [
            [1, 7, 8, 0, 0, 0],
            [7, 1, 8, 0, 0, 0],
            [8, 8, 3, 0, 0, 0],
            [0, 0, 0, 4, 0, 0],
            [0, 0, 0, 0, 4, 0],
            [0, 0, 0, 0, 0, 6],
        ]
    ),
    "tetragonal_low": np.array(
        [
            [1, 7, 8, 0, 0, 11],
            [7, 1, 8, 0, 0, -11],
            [8, 8, 3, 0, 0, 0],
            [0, 0, 0, 4, 0, 0],
            [0, 0, 0, 0, 4, 0],
            [11, -11, 0, 0, 0, 6],
        ]
    ),
    "orthorhombic": np.array(
        [
            [1, 7, 8, 0, 0, 0],
            [7, 2, 12, 0, 0, 0],
            [8, 12, 3, 0, 0, 0],
            [0, 0, 0, 4, 0, 0],
            [0, 0, 0, 0, 5, 0],
            [0, 0, 0, 0, 0, 6],
        ]
    ),
    "monoclinic": np.array(
        [
            [1, 7, 8, 0, 10, 0],
            [7, 2, 12, 0, 14, 0],
            [8, 12, 3, 0, 17, 0],
            [0, 0, 0, 4, 0, 20],
            [10, 14, 17, 0, 5, 0],
            [0, 0, 0, 20, 0, 6],
        ]
    ),
    "triclinic": np.array(
        [
            [1, 7, 8, 9, 10, 11],
            [7, 2, 12, 13, 14, 15],
            [8, 12, 3, 16, 17, 18],
            [9, 13, 16, 4, 19, 20],
            [10, 14, 17, 19, 5, 21],
            [11, 15, 18, 20, 21, 6],
        ]
    ),
}


@dataclass(frozen=True, slots=True)
class ElasticComponentDefinition:
    """Describe one independent or derived Voigt stiffness component.

    Parameters
    ----------
    label : str
        Canonical component label such as ``"C11"``.
    entries : tuple of tuple
        Matrix entries represented by this component. Each item is
        ``(row, column, multiplier)`` and satisfies
        ``C[row, column] = multiplier * value``.
    derived_from : tuple of tuple, optional
        Linear combination ``(label, coefficient)`` used for components that
        are constrained by symmetry rather than fitted independently.
    """

    label: str
    entries: tuple[tuple[int, int, float], ...]
    derived_from: tuple[tuple[str, float], ...] = ()

    @property
    def derived(self) -> bool:
        """Return whether the component is obtained from other components."""
        return bool(self.derived_from)


def elastic_component_definitions(
    symmetry: str,
) -> tuple[ElasticComponentDefinition, ...]:
    """Return canonical independent and derived components for a symmetry.

    The definitions follow Quantas' Voigt order ``11, 22, 33, 23, 13, 12``.
    Symmetry-equivalent entries, including sign changes, are grouped together.
    For hexagonal and trigonal systems, ``C66`` is represented by the exact
    relation ``(C11-C12)/2`` and is therefore not fitted independently.

    Parameters
    ----------
    symmetry : str
        Elastic-symmetry identifier.

    Returns
    -------
    tuple of ElasticComponentDefinition
        Independent components followed by derived components.

    Raises
    ------
    KeyError
        If the symmetry is unknown.
    """
    if symmetry not in ELASTIC_SYMMETRY_PATTERNS:
        raise KeyError(f"Unknown crystal system: {symmetry}")
    pattern = ELASTIC_SYMMETRY_PATTERNS[symmetry]
    ids: list[int] = []
    for i in range(6):
        for j in range(i, 6):
            component_id = abs(int(pattern[i, j]))
            if component_id and component_id not in ids:
                ids.append(component_id)

    derived_ids = (
        {6} if symmetry in {"hexagonal", "trigonal_high", "trigonal_low"} else set()
    )
    definitions: list[ElasticComponentDefinition] = []
    for component_id in ids:
        if component_id in derived_ids:
            continue
        entries: list[tuple[int, int, float]] = []
        representative: tuple[int, int] | None = None
        for i in range(6):
            for j in range(i, 6):
                value = int(pattern[i, j])
                if abs(value) != component_id:
                    continue
                multiplier = 1.0 if value > 0 else -1.0
                entries.append((i, j, multiplier))
                if representative is None and multiplier > 0.0:
                    representative = (i, j)
        if representative is None:
            representative = (entries[0][0], entries[0][1])
        label = f"C{representative[0] + 1}{representative[1] + 1}"
        definitions.append(
            ElasticComponentDefinition(label=label, entries=tuple(entries))
        )

    if derived_ids:
        definitions.append(
            ElasticComponentDefinition(
                label="C66",
                entries=((5, 5, 1.0),),
                derived_from=(("C11", 0.5), ("C12", -0.5)),
            )
        )
    return tuple(definitions)


def extract_independent_stiffness_components(
    stiffness: np.ndarray,
    symmetry: str,
) -> tuple[dict[str, float], dict[str, float]]:
    """Average symmetry-equivalent entries of one stiffness matrix.

    Parameters
    ----------
    stiffness : ndarray
        Symmetric ``(6, 6)`` stiffness matrix.
    symmetry : str
        Elastic-symmetry identifier.

    Returns
    -------
    tuple of dict
        Independent component values and maximum absolute symmetry spreads.
        Derived components are omitted from the first mapping.

    Raises
    ------
    ValueError
        If the stiffness matrix is invalid.
    """
    matrix = np.asarray(stiffness, dtype=np.float64)
    if matrix.shape != (6, 6) or not np.all(np.isfinite(matrix)):
        raise ValueError("stiffness must be a finite matrix with shape (6, 6)")
    if not np.allclose(matrix, matrix.T, rtol=0.0, atol=1.0e-10):
        raise ValueError("stiffness matrix must be symmetric")
    values: dict[str, float] = {}
    spreads: dict[str, float] = {}
    for definition in elastic_component_definitions(symmetry):
        if definition.derived:
            continue
        samples = np.asarray(
            [multiplier * matrix[i, j] for i, j, multiplier in definition.entries],
            dtype=np.float64,
        )
        mean = float(np.mean(samples))
        values[definition.label] = mean
        spreads[definition.label] = float(np.max(np.abs(samples - mean)))
    return values, spreads


def reconstruct_stiffness_from_components(
    components: Mapping[str, float],
    symmetry: str,
) -> np.ndarray:
    """Reconstruct a symmetry-constrained stiffness matrix.

    Parameters
    ----------
    components : mapping
        Independent component values keyed by canonical labels.
    symmetry : str
        Elastic-symmetry identifier.

    Returns
    -------
    ndarray
        Symmetric ``(6, 6)`` stiffness matrix.

    Raises
    ------
    KeyError
        If an independent component is missing.
    """
    matrix = np.zeros((6, 6), dtype=np.float64)
    resolved: dict[str, float] = {
        str(key): float(value) for key, value in components.items()
    }
    for definition in elastic_component_definitions(symmetry):
        if definition.derived:
            value = sum(
                resolved[label] * coefficient
                for label, coefficient in definition.derived_from
            )
            resolved[definition.label] = float(value)
        else:
            if definition.label not in resolved:
                raise KeyError(
                    f"missing independent elastic component {definition.label}"
                )
            value = resolved[definition.label]
        for i, j, multiplier in definition.entries:
            matrix[i, j] = multiplier * value
            matrix[j, i] = multiplier * value
    return matrix


def stiffness_component_linear_coefficients(
    symmetry: str,
    labels: Sequence[str],
) -> np.ndarray:
    """Return the linear map from independent components to a Voigt matrix.

    Parameters
    ----------
    symmetry : str
        Elastic-symmetry identifier.
    labels : sequence of str
        Independent component order.

    Returns
    -------
    ndarray
        Array with shape ``(6, 6, ncomponents)``.
    """
    order = tuple(str(label) for label in labels)
    index = {label: position for position, label in enumerate(order)}
    coefficients = np.zeros((6, 6, len(order)), dtype=np.float64)
    for definition in elastic_component_definitions(symmetry):
        vector = np.zeros(len(order), dtype=np.float64)
        if definition.derived:
            for label, coefficient in definition.derived_from:
                if label not in index:
                    raise KeyError(f"missing independent elastic component {label}")
                vector[index[label]] += coefficient
        else:
            if definition.label not in index:
                raise KeyError(
                    f"missing independent elastic component {definition.label}"
                )
            vector[index[definition.label]] = 1.0
        for i, j, multiplier in definition.entries:
            coefficients[i, j] = multiplier * vector
            coefficients[j, i] = multiplier * vector
    return coefficients


def detect_elastic_symmetry(matrix: np.ndarray, tolerance: float = 1.0e-3) -> str:
    """
    Return the crystal system compatible with an elastic matrix.

    Parameters
    ----------
    matrix : ndarray
        Elastic stiffness matrix in Voigt notation. The expected shape is
        ``(6, 6)``.
    tolerance : float, optional
        Numerical tolerance used to compare tensor components.

    Returns
    -------
    str
        Name of the detected elastic crystal system.

    Raises
    ------
    ValueError
        If the matrix does not have shape ``(6, 6)`` or if no compatible
        crystal system is found.
    """
    matrix = np.asarray(matrix, dtype=float)

    if matrix.shape != (6, 6):
        raise ValueError("The elastic stiffness matrix must have shape (6, 6).")

    for symmetry in ELASTIC_SYMMETRY_PATTERNS:
        if is_elastic_symmetry(matrix, symmetry, tolerance=tolerance):
            return symmetry

    raise ValueError("No compatible crystal system was found.")


def is_elastic_symmetry(
    matrix: np.ndarray,
    symmetry: str,
    tolerance: float = 1.0e-3,
) -> bool:
    """
    Check whether an elastic matrix is compatible with a crystal system.

    Parameters
    ----------
    matrix : ndarray
        Elastic stiffness matrix in Voigt notation. The expected shape is
        ``(6, 6)``.
    symmetry : str
        Name of the crystal system to test.
    tolerance : float, optional
        Numerical tolerance used to compare tensor components.

    Returns
    -------
    bool
        ``True`` if the matrix is compatible with the requested symmetry,
        otherwise ``False``.

    Raises
    ------
    ValueError
        If the matrix does not have shape ``(6, 6)``.
    KeyError
        If the requested symmetry is not available.
    """
    matrix = np.asarray(matrix, dtype=float)

    if matrix.shape != (6, 6):
        raise ValueError("The elastic stiffness matrix must have shape (6, 6).")

    if symmetry not in ELASTIC_SYMMETRY_PATTERNS:
        raise KeyError(f"Unknown crystal system: {symmetry}")

    component_map = _build_component_map(symmetry)

    for component_id in range(22):
        indexes = [
            index for index, value in component_map.items() if value == component_id
        ]

        if not indexes:
            continue

        if component_id == 0:
            for index in indexes:
                if abs(matrix[index]) > tolerance:
                    return False
        else:
            for index_a, index_b in it.combinations(indexes, 2):
                if abs(matrix[index_a] - matrix[index_b]) > tolerance:
                    return False

    if symmetry == "hexagonal":
        c66 = matrix[5, 5]
        expected_c66 = 0.5 * (matrix[0, 0] - matrix[0, 1])
        if abs(c66 - expected_c66) > tolerance:
            return False

    return True


def _build_component_map(symmetry: str) -> dict[tuple[int, int], int]:
    """
    Build the full component map for a crystal system.

    Parameters
    ----------
    symmetry : str
        Name of the crystal system.

    Returns
    -------
    dict
        Dictionary mapping matrix indexes to component identifiers.
    """
    symmetry_map = ELASTIC_SYMMETRY_PATTERNS[symmetry]
    component_map = {}

    for i in range(6):
        for j in range(i, 6):
            component_map[(i, j)] = symmetry_map[i, j]

    for i, j in list(component_map):
        component_map[(j, i)] = component_map[(i, j)]

    return component_map


def specialize_elastic_tensor(
    tensor: "ElasticTensor",
    crystal_system: str,
) -> "ElasticTensor":
    """Return a symmetry-specialized elastic tensor when available.

    Parameters
    ----------
    tensor : ElasticTensor
        General elastic tensor.
    crystal_system : str
        Detected crystal-system name.

    Returns
    -------
    ElasticTensor
        Orthorhombic specialization for orthorhombic tensors, otherwise the
        original tensor.
    """
    from .tensor import OrthorhombicElasticTensor

    if crystal_system == "orthorhombic":
        return OrthorhombicElasticTensor(tensor)
    return tensor


__all__ = [
    "ELASTIC_SYMMETRY_PATTERNS",
    "ElasticComponentDefinition",
    "detect_elastic_symmetry",
    "elastic_component_definitions",
    "extract_independent_stiffness_components",
    "is_elastic_symmetry",
    "reconstruct_stiffness_from_components",
    "specialize_elastic_tensor",
    "stiffness_component_linear_coefficients",
]
