# -*- coding: utf-8 -*-

"""Neutral report-table builders for elasticity results."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

import numpy as np

from quantas.core.geometry import TensorRotation

from quantas.models.report import ReportTable as _ReportTable

if TYPE_CHECKING:
    from quantas.models.report import ReportTable
from quantas.modules.elasticity.models import (
    ElasticityInput,
    ElasticityOptions,
    ElasticityResult,
)


_SINGLE_DIRECTION_PROPERTIES = frozenset(
    {"young_modulus", "linear_compressibility"}
)
_PAIRED_DIRECTION_PROPERTIES = frozenset({"shear_modulus", "poisson_ratio"})


def input_table(
    input_data: ElasticityInput,
    result: ElasticityResult | None = None,
) -> ReportTable:
    """Build a neutral table describing elasticity input data."""
    rows: list[list[object]] = [
        ["Job name", input_data.jobname],
        [
            "Source",
            str(input_data.source) if input_data.source is not None else "Unknown",
        ],
    ]
    if result is not None:
        rows.insert(1, ["Crystal system", result.crystal_system or "Unknown"])
    return _ReportTable(
        title="Elasticity input", columns=["Property", "Value"], rows=rows
    )


def options_table(options: ElasticityOptions) -> ReportTable:
    """Build a neutral table describing scientific elasticity options."""
    return _ReportTable(
        title="Elasticity options",
        columns=["Option", "Value"],
        rows=[
            ["Stiffness unit", options.pressure_unit],
            ["Calculate 2D properties", options.calculate_2d],
            ["Angular points (ntheta)", options.ntheta],
            [
                "Tensor rotation",
                "none" if options.rotation is None else options.rotation.kind.value,
            ],
            ["Numerical precision", "float64"],
        ],
    )


def stiffness_table(
    result: ElasticityResult,
    *,
    title: str = "Stiffness matrix (GPa)",
) -> ReportTable:
    """Build a neutral stiffness-matrix table.

    Parameters
    ----------
    result : ElasticityResult
        Elasticity result containing the analysis-frame matrix.
    title : str, optional
        Table title.

    Returns
    -------
    ReportTable
        Matrix table.
    """
    return _matrix_table(title, result.stiffness, "C")


def source_stiffness_table(matrix: np.ndarray) -> ReportTable:
    """Build a source-frame stiffness-matrix table.

    Parameters
    ----------
    matrix : ndarray
        Source stiffness matrix in GPa.

    Returns
    -------
    ReportTable
        Source-frame matrix table.
    """
    return _matrix_table("Stiffness matrix before rotation (GPa)", matrix, "C")


def tensor_rotation_table(rotation: TensorRotation) -> ReportTable:
    """Build a table documenting a user-defined tensor transformation.

    Parameters
    ----------
    rotation : TensorRotation
        Transformation applied before analysis.

    Returns
    -------
    ReportTable
        Rotation provenance and matrix components.
    """
    return tensor_rotation_metadata_table(rotation.as_mapping())


def tensor_rotation_metadata_table(
    metadata: Mapping[str, Any],
) -> ReportTable:
    """Build a transformation table from serialized workflow metadata.

    Parameters
    ----------
    metadata : mapping
        Tensor-frame metadata containing ``component_transform``.

    Returns
    -------
    ReportTable
        Rotation provenance and matrix components.

    Raises
    ------
    ValueError
        If the stored transformation matrix is missing or malformed.
    """
    matrix = np.asarray(metadata.get("component_transform"), dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError("Tensor rotation metadata must contain a 3 x 3 matrix.")
    rows: list[list[object]] = [
        ["Kind", metadata.get("kind", "matrix"), "", ""],
        [
            "Convention",
            metadata.get(
                "convention",
                "C'_ijkl = R_ia R_jb R_kc R_ld C_abcd",
            ),
            "",
            "",
        ],
    ]
    angles = metadata.get("angles")
    if angles is not None:
        values = np.asarray(angles, dtype=float)
        rows.append(
            [
                "Fixed-axis XYZ angles",
                f"{values[0]:.10g}",
                f"{values[1]:.10g}",
                f"{values[2]:.10g} {metadata.get('angle_unit', '')}".strip(),
            ]
        )
    for index, row in enumerate(matrix, start=1):
        rows.append(
            [
                f"R{index}",
                f"{row[0]:.12f}",
                f"{row[1]:.12f}",
                f"{row[2]:.12f}",
            ]
        )
    return _ReportTable(
        title="Tensor component transformation",
        columns=["Entry", "1", "2", "3"],
        rows=rows,
        metadata={
            "source_frame": metadata.get("source_frame", "source"),
            "analysis_frame": metadata.get("analysis_frame", "rotated"),
            "transformed": True,
        },
    )


def compliance_table(result: ElasticityResult) -> ReportTable:
    """Build a neutral compliance-matrix table."""
    return _matrix_table("Compliance matrix (GPa^-1)", result.compliance, "S")


def averages_table(result: ElasticityResult) -> ReportTable:
    """Build a Voigt-Reuss-Hill average-property table."""
    rows: list[list[Any]] = []
    if result.averages is not None:
        for label, properties in (
            ("Voigt", result.averages.voigt),
            ("Reuss", result.averages.reuss),
            ("Hill", result.averages.hill),
        ):
            rows.append(
                [
                    label,
                    f"{properties.bulk_modulus:.4f}",
                    f"{properties.young_modulus:.4f}",
                    f"{properties.shear_modulus:.4f}",
                    f"{properties.poisson_ratio:.6f}",
                ]
            )
    return _ReportTable(
        title="Voigt-Reuss-Hill average properties",
        columns=["Scheme", "K (GPa)", "E (GPa)", "G (GPa)", "nu"],
        rows=rows,
    )


def stability_table(result: ElasticityResult) -> ReportTable:
    """Build a positive-definiteness and eigenvalue table."""
    rows: list[list[Any]] = []
    metadata: dict[str, object] = {"positive_definite": False}
    if result.stability is not None:
        metadata["positive_definite"] = result.stability.is_stable
        metadata["tolerance"] = result.stability.tolerance
        for index, value in enumerate(result.stability.eigenvalues, start=1):
            rows.append([index, f"{value:.6f}"])
    return _ReportTable(
        title="Mechanical stability",
        columns=["Eigenvalue", "Value (GPa)"],
        rows=rows,
        metadata=metadata,
    )


def variations_table(
    result: ElasticityResult,
    *,
    paired_directions: bool | None = None,
) -> ReportTable:
    """Build a normalized table of directional elastic extrema.

    Parameters
    ----------
    result : ElasticityResult
        Elasticity result containing directional extrema.
    paired_directions : bool or None, optional
        Select properties that require only the primary direction ``a``
        (``False``) or the orthogonal pair ``(a, b)`` (``True``). If ``None``,
        include all properties and expose the optional ``b`` direction.

    Returns
    -------
    ReportTable
        Frontend-neutral long-form extrema table. Each minimum and maximum is
        represented by one self-contained row.
    """
    items = [
        (name, variation)
        for name, variation in result.variations.items()
        if paired_directions is None
        or _uses_paired_directions(name, variation) is paired_directions
    ]
    include_measurement_axis = paired_directions is not False
    rows: list[list[Any]] = []
    for name, variation in items:
        label = _directional_property_label(name)
        extrema = (
            (
                "Minimum",
                variation.minimum,
                variation.minimum_axis,
                variation.minimum_measurement_axis,
            ),
            (
                "Maximum",
                variation.maximum,
                variation.maximum_axis,
                variation.maximum_measurement_axis,
            ),
        )
        for extremum, value, primary_axis, measurement_axis in extrema:
            row: list[Any] = [
                label,
                extremum,
                value,
                _format_axis(primary_axis),
            ]
            if include_measurement_axis:
                row.append(_format_axis(measurement_axis, missing=""))
            row.append(variation.anisotropy)
            rows.append(row)

    if paired_directions is False:
        title = "Directional extrema — single-direction properties"
        notes = [
            "a is the primary direction; its Cartesian components refer to "
            "the analysis frame."
        ]
    elif paired_directions is True:
        title = "Directional extrema — paired-direction properties"
        notes = [
            "a is the primary direction; b is the orthogonal transverse "
            "measurement direction (a · b = 0)."
        ]
    else:
        title = "Directional elastic extrema"
        notes = [
            "a is the primary direction; where present, b is the orthogonal "
            "transverse measurement direction (a · b = 0)."
        ]

    columns = ["Property", "Extremum", "Value", "a: primary direction"]
    column_formats: list[str | None] = [None, None, ".4f", None]
    column_alignments = ["left", "left", "right", "left"]
    if include_measurement_axis:
        columns.append("b: transverse direction")
        column_formats.append(None)
        column_alignments.append("left")
    columns.append("Ratio max/min")
    column_formats.append(".4f")
    column_alignments.append("right")

    direction_roles = {"a": "primary"}
    if include_measurement_axis:
        direction_roles["b"] = "orthogonal transverse measurement"

    return _ReportTable(
        title=title,
        columns=columns,
        rows=rows,
        metadata={
            "column_formats": column_formats,
            "column_alignments": column_alignments,
            "notes": notes,
            "direction_roles": direction_roles,
        },
    )


def variations_tables(result: ElasticityResult) -> list[ReportTable]:
    """Build separate extrema tables for one- and two-direction properties.

    Parameters
    ----------
    result : ElasticityResult
        Elasticity result containing directional extrema.

    Returns
    -------
    list of ReportTable
        Ordered non-empty tables. Properties depending only on ``a`` precede
        properties depending on the orthogonal pair ``(a, b)``.
    """
    tables = [
        variations_table(result, paired_directions=False),
        variations_table(result, paired_directions=True),
    ]
    return [table for table in tables if table.rows]


def build_elasticity_report(
    input_data: ElasticityInput,
    options: ElasticityOptions,
    result: ElasticityResult,
) -> list[ReportTable]:
    """Build the complete neutral elasticity report."""
    tables = [
        input_table(input_data, result),
        options_table(options),
    ]
    if options.rotation is not None and input_data.stiffness is not None:
        tables.extend(
            [
                source_stiffness_table(np.asarray(input_data.stiffness, dtype=float)),
                tensor_rotation_table(options.rotation),
                stiffness_table(
                    result,
                    title="Stiffness matrix after rotation (GPa)",
                ),
            ]
        )
    else:
        tables.append(stiffness_table(result))
    tables.extend(
        [
            compliance_table(result),
            averages_table(result),
            stability_table(result),
        ]
    )
    if result.variations:
        tables.extend(variations_tables(result))
    return tables


def _matrix_table(
    title: str,
    matrix: np.ndarray | None,
    row_prefix: str,
) -> ReportTable:
    """Build a neutral table from a ``6 x 6`` matrix."""
    columns = [""] + [f"{row_prefix}{index}" for index in range(1, 7)]
    rows: list[list[Any]] = []
    if matrix is not None:
        for index, row in enumerate(np.asarray(matrix), start=1):
            rows.append([f"{row_prefix}{index}"] + [f"{value:.6f}" for value in row])
    return _ReportTable(title=title, columns=columns, rows=rows)


def _format_axis(axis: list[float] | None, *, missing: str = "None") -> str:
    """Format an optional Cartesian direction with aligned components."""
    if axis is None:
        return missing
    values = [0.0 if abs(value) < 5.0e-4 else value for value in axis]
    return "[{: >6.3f}, {: >6.3f}, {: >6.3f}]".format(*values)


def _uses_paired_directions(name: str, variation: Any) -> bool:
    """Return whether a property depends on an orthogonal direction pair."""
    if name in _SINGLE_DIRECTION_PROPERTIES:
        return False
    if name in _PAIRED_DIRECTION_PROPERTIES:
        return True
    return (
        variation.minimum_measurement_axis is not None
        or variation.maximum_measurement_axis is not None
    )


def _directional_property_label(name: str) -> str:
    """Return a readable directional-property label with its physical unit."""
    labels = {
        "young_modulus": "Young's modulus (GPa)",
        "linear_compressibility": "Linear compressibility (TPa^-1)",
        "shear_modulus": "Shear modulus (GPa)",
        "poisson_ratio": "Poisson's ratio",
    }
    return labels.get(name, name.replace("_", " ").strip().capitalize())
