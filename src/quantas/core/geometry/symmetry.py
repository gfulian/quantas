# -*- coding: utf-8 -*-

"""spglib-backed crystal-symmetry analysis and primitive-cell utilities."""

from __future__ import annotations

from typing import Any

import numpy as np

from quantas.models.structures import CrystalStructure, SymmetryMetadata


def _dataset_value(dataset: Any, name: str, default: Any = None) -> Any:
    """Return a value from old mapping-style or new attribute-style datasets."""
    if dataset is None:
        return default
    if hasattr(dataset, name):
        return getattr(dataset, name)
    try:
        return dataset[name]
    except (KeyError, TypeError):
        return default


def analyze_symmetry(
    structure: CrystalStructure,
    *,
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
) -> SymmetryMetadata:
    """Determine crystallographic symmetry using spglib.

    Parameters
    ----------
    structure : CrystalStructure
        Structure expressed by row-wise lattice vectors, fractional positions,
        and atomic numbers.
    symprec : float, optional
        Cartesian symmetry tolerance in angstrom.
    angle_tolerance : float, optional
        Angular tolerance in degrees. Negative values request spglib's default.

    Returns
    -------
    SymmetryMetadata
        Space-group identifiers, equivalence mapping, and standard-setting
        transformations.

    Raises
    ------
    ImportError
        If spglib is unavailable.
    ValueError
        If spglib cannot determine a symmetry dataset.
    """
    try:
        import spglib
    except ImportError as exc:  # pragma: no cover - dependency contract
        raise ImportError(
            "spglib is required for Quantas structural symmetry analysis"
        ) from exc

    if hasattr(spglib, "OLD_ERROR_HANDLING"):
        spglib.OLD_ERROR_HANDLING = False
    if hasattr(spglib, "error") and hasattr(spglib.error, "OLD_ERROR_HANDLING"):
        spglib.error.OLD_ERROR_HANDLING = False
    try:
        dataset = spglib.get_symmetry_dataset(
            structure.spglib_cell(),
            symprec=float(symprec),
            angle_tolerance=float(angle_tolerance),
        )
    except spglib.SpglibError as exc:
        raise ValueError("spglib could not determine a symmetry dataset") from exc
    if dataset is None:
        raise ValueError("spglib could not determine a symmetry dataset")
    hall_number = int(_dataset_value(dataset, "hall_number", 0))
    hall_symbol = ""
    if hall_number > 0:
        spacegroup_type = spglib.get_spacegroup_type(hall_number)
        if spacegroup_type is not None:
            hall_symbol = str(_dataset_value(spacegroup_type, "hall_symbol", ""))
    return SymmetryMetadata(
        space_group_number=int(_dataset_value(dataset, "number", 0)),
        international_symbol=str(_dataset_value(dataset, "international", "")),
        hall_number=hall_number,
        hall_symbol=hall_symbol,
        choice=str(_dataset_value(dataset, "choice", "")),
        point_group=str(_dataset_value(dataset, "pointgroup", "")),
        symprec=float(symprec),
        angle_tolerance=float(angle_tolerance),
        equivalent_atoms=np.asarray(
            _dataset_value(dataset, "equivalent_atoms", []),
            dtype=np.int64,
        ),
        transformation_matrix=np.asarray(
            _dataset_value(dataset, "transformation_matrix", np.eye(3)),
            dtype=np.float64,
        ),
        origin_shift=np.asarray(
            _dataset_value(dataset, "origin_shift", np.zeros(3)),
            dtype=np.float64,
        ),
    )


def find_primitive_structure(
    structure: CrystalStructure,
    *,
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
    no_idealize: bool = True,
) -> CrystalStructure:
    """Return a primitive structure found by spglib.

    This function is intended for validation and fallback use. Workflows that
    possess an explicit source-code expansion matrix should use that matrix to
    preserve the source orientation rather than independently standardizing
    every sampled volume.

    Parameters
    ----------
    structure : CrystalStructure
        Input periodic structure.
    symprec : float, optional
        Cartesian symmetry tolerance in angstrom.
    angle_tolerance : float, optional
        Angular tolerance in degrees.
    no_idealize : bool, optional
        Preserve non-idealized coordinates and lattice as far as spglib allows.

    Returns
    -------
    CrystalStructure
        Primitive structure returned by spglib.

    Raises
    ------
    ImportError
        If spglib is unavailable.
    ValueError
        If primitive standardization fails.
    """
    try:
        import spglib
    except ImportError as exc:  # pragma: no cover - dependency contract
        raise ImportError(
            "spglib is required for Quantas structural symmetry analysis"
        ) from exc
    if hasattr(spglib, "OLD_ERROR_HANDLING"):
        spglib.OLD_ERROR_HANDLING = False
    if hasattr(spglib, "error") and hasattr(spglib.error, "OLD_ERROR_HANDLING"):
        spglib.error.OLD_ERROR_HANDLING = False
    try:
        primitive = spglib.standardize_cell(
            structure.spglib_cell(),
            to_primitive=True,
            no_idealize=bool(no_idealize),
            symprec=float(symprec),
            angle_tolerance=float(angle_tolerance),
        )
    except spglib.SpglibError as exc:
        raise ValueError("spglib could not construct a primitive structure") from exc
    if primitive is None:
        raise ValueError("spglib could not construct a primitive structure")
    lattice, positions, numbers = primitive
    return CrystalStructure(
        lattice=np.asarray(lattice, dtype=np.float64),
        fractional_positions=np.asarray(positions, dtype=np.float64),
        atomic_numbers=np.asarray(numbers, dtype=np.int64),
        label=f"spglib primitive of {structure.label}".strip(),
        metadata={"symprec": float(symprec), "no_idealize": bool(no_idealize)},
    )
