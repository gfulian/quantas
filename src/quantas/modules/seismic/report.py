# -*- coding: utf-8 -*-

"""Neutral report-table builders for sampled seismic-wave results."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

from quantas.core.geometry import cartesian_to_spherical
from quantas.core.physics.elasticity import detect_elastic_symmetry
from quantas.core.physics.seismic import (
    MODE_INDEX,
    MODE_ORDER,
    MODE_PAIR_INDEX,
    MODE_PAIR_ORDER,
    MODE_PAIR_SYMBOLS,
    MODE_SYMBOLS,
    ModePair,
    WaveMode,
)
from quantas.models import ReportTable
from quantas.modules.seismic.models import SeismicResult


@dataclass(frozen=True, slots=True)
class SampledExtremum:
    """One extremum identified on the sampled spherical grid.

    Parameters
    ----------
    value : float
        Extremal scalar value.
    direction : ndarray or None
        Representative tensor-frame Cartesian direction. ``None`` indicates
        that the quantity is constant over all eligible sampled directions.
    ray_direction : ndarray or None
        Representative energy-flow direction when available.
    theta_degree, phi_degree : float or None
        Polar and azimuthal angles of the representative direction in degrees.
    multiplicity : int
        Number of unique antipodal sampled axes sharing the extremal value.
    all_directions : bool
        Whether every eligible sampled direction has the same value.
    """

    value: float
    direction: NDArray[np.float64] | None
    ray_direction: NDArray[np.float64] | None
    theta_degree: float | None
    phi_degree: float | None
    multiplicity: int
    all_directions: bool


@dataclass(frozen=True, slots=True)
class SampledVariation:
    """Minimum, maximum, and anisotropy of a sampled scalar field.

    Parameters
    ----------
    minimum, maximum : SampledExtremum
        Sampled extrema.
    anisotropy_percent : float
        Symmetric percentage anisotropy,
        ``200 * (maximum - minimum) / (maximum + minimum)``.
    maximum_to_minimum : float
        Ratio between maximum and minimum values.
    """

    minimum: SampledExtremum
    maximum: SampledExtremum
    anisotropy_percent: float
    maximum_to_minimum: float


def build_seismic_report_tables(
    result: SeismicResult,
    *,
    level: Literal["standard", "extended", "debug"] = "extended",
) -> list[ReportTable]:
    """Build a detailed frontend-neutral report for a seismic result.

    Parameters
    ----------
    result : SeismicResult
        Complete sampled seismic-wave result.
    level : {"standard", "extended", "debug"}, optional
        Scientific report detail. ``standard`` retains the principal input,
        stability, isotropic-reference, extrema, and diagnostic summaries.
        ``extended`` adds conventions, power-flow details, and candidate axes.
        ``debug`` additionally includes the full stiffness matrix and available
        tensor-frame provenance.

    Returns
    -------
    list of ReportTable
        Ordered report tables describing input, extrema, anisotropy, wave
        splitting, power flow, enhancement, and numerical diagnostics.
    """
    if level not in {"standard", "extended", "debug"}:
        raise ValueError("seismic report level must be standard, extended, or debug")

    tables = [
        calculation_summary_table(result),
        elastic_reference_table(result),
        stability_table(result),
        isotropic_reference_table(result),
        phase_velocity_extrema_table(result),
        phase_derived_properties_table(result),
    ]
    if result.field.group is not None:
        tables.append(group_velocity_extrema_table(result))
    if result.field.enhancement is not None:
        tables.append(enhancement_table(result))
    tables.append(diagnostics_table(result))
    if level in {"extended", "debug"}:
        tables.insert(4, acoustic_mode_conventions_table())
        if result.field.group is not None:
            tables.append(power_flow_table(result))
        candidates = acoustic_axis_candidates_table(result)
        if candidates.rows:
            tables.append(candidates)
    return tables


def stiffness_matrix_table(
    matrix: NDArray[np.float64],
    *,
    title: str,
) -> ReportTable:
    """Build a neutral seismic stiffness-matrix table.

    Parameters
    ----------
    matrix : ndarray
        Stiffness matrix in GPa.
    title : str
        Table title identifying the tensor frame.

    Returns
    -------
    ReportTable
        Matrix table.
    """
    values = np.asarray(matrix, dtype=float)
    rows = [
        [f"C{index}"] + [f"{value:.6f}" for value in row]
        for index, row in enumerate(values, start=1)
    ]
    return ReportTable(
        title=title,
        columns=[""] + [f"C{index}" for index in range(1, 7)],
        rows=rows,
    )


def tensor_rotation_metadata_table(
    metadata: Mapping[str, Any],
) -> ReportTable:
    """Build a table documenting the tensor transformation used by SEISMIC.

    Parameters
    ----------
    metadata : mapping
        Serialized tensor-frame provenance.

    Returns
    -------
    ReportTable
        Rotation provenance and matrix components.

    Raises
    ------
    ValueError
        If the component transformation is malformed.
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
    return ReportTable(
        title="Tensor component transformation",
        columns=["Entry", "1", "2", "3"],
        rows=rows,
        metadata={
            "source_frame": metadata.get("source_frame", "source"),
            "analysis_frame": metadata.get("analysis_frame", "rotated"),
            "transformed": True,
        },
    )


def calculation_summary_table(result: SeismicResult) -> ReportTable:
    """Build a table describing the calculation and spherical sampling."""
    try:
        symmetry = detect_elastic_symmetry(result.stiffness)
    except ValueError:
        symmetry = "unknown"
    grid = result.grid
    rows: list[list[object]] = [
        ["Job name", result.jobname],
        ["Elastic symmetry", symmetry],
        ["Density / kg m^-3", f"{result.density:.6f}"],
        [
            "Mechanical stability",
            "stable" if result.stability.is_stable else "unstable",
        ],
        ["Sampling level", result.field.level.value],
        ["Hemisphere", grid.hemisphere.value],
        ["Polar samples", grid.shape[0]],
        ["Azimuthal samples", grid.shape[1]],
        ["Sampled directions", grid.size],
        ["Azimuthal seam duplicated", "no"],
        ["Polarization tracking", "yes" if result.field.tracking is not None else "no"],
    ]
    return ReportTable(
        title="Seismic calculation summary",
        columns=["Property", "Value"],
        rows=rows,
        metadata={
            "direction_frame": "elastic-tensor Cartesian frame",
            "theta_definition": "polar angle from +z",
            "phi_definition": "azimuth from +x toward +y",
        },
    )


def elastic_reference_table(result: SeismicResult) -> ReportTable:
    """Build a table of Voigt, Reuss, and Hill isotropic elastic estimates."""
    rows: list[list[object]] = []
    for label, values in (
        ("Voigt", result.averages.voigt),
        ("Reuss", result.averages.reuss),
        ("Hill", result.averages.hill),
    ):
        rows.append(
            [
                label,
                f"{values.bulk_modulus:.6f}",
                f"{values.shear_modulus:.6f}",
                f"{values.young_modulus:.6f}",
                f"{values.poisson_ratio:.8f}",
            ]
        )
    return ReportTable(
        title="Isotropic elastic reference properties",
        columns=["Scheme", "K / GPa", "G / GPa", "E / GPa", "Poisson ratio"],
        rows=rows,
    )


def stability_table(result: SeismicResult) -> ReportTable:
    """Build a stiffness-eigenvalue table for mechanical stability."""
    rows = [
        [index, f"{value:.8f}"]
        for index, value in enumerate(result.stability.eigenvalues, start=1)
    ]
    return ReportTable(
        title="Mechanical stability",
        columns=["Stiffness eigenvalue", "Value / GPa"],
        rows=rows,
        metadata={
            "positive_definite": result.stability.is_stable,
            "tolerance": result.stability.tolerance,
        },
    )


def isotropic_reference_table(result: SeismicResult) -> ReportTable:
    """Build a table of Hill-average isotropic acoustic velocities."""
    v_s = result.isotropic_velocities.shear
    v_p = result.isotropic_velocities.compressional
    rows = [
        ["V_S", f"{v_s:.8f}"],
        ["V_P", f"{v_p:.8f}"],
        ["V_P / V_S", f"{v_p / v_s:.8f}"],
    ]
    return ReportTable(
        title="Hill-average isotropic velocity reference",
        columns=["Quantity", "Value"],
        rows=rows,
        metadata={"velocity_unit": "km s^-1"},
    )


def acoustic_mode_conventions_table() -> ReportTable:
    """Build a table defining the acoustic-mode symbols used in reports."""
    return ReportTable(
        title="Acoustic-mode conventions",
        columns=["Symbol", "Physical description", "Local phase-speed order"],
        rows=[
            ["V_P", "quasi-longitudinal wave", "fastest"],
            ["V_S1", "fast quasi-shear wave", "intermediate"],
            ["V_S2", "slow quasi-shear wave", "slowest"],
        ],
        metadata={
            "ordering": "V_S2 <= V_S1 <= V_P at each sampled wave normal",
        },
    )


def phase_velocity_extrema_table(result: SeismicResult) -> ReportTable:
    """Build sampled phase-speed extrema and anisotropy for all modes."""
    rows: list[list[object]] = []
    phase = result.field.phase
    for mode in reversed(MODE_ORDER):
        field = phase.for_mode(mode)
        variation = sampled_variation(
            field.phase_speeds,
            field.valid_mask,
            phase.directions,
        )
        rows.append(_variation_row(MODE_SYMBOLS[mode], variation))
    return ReportTable(
        title="Sampled phase-velocity extrema",
        columns=_variation_columns(include_ray=False),
        rows=rows,
        metadata={
            "unit": "km s^-1",
            "anisotropy_definition": "200 * (maximum - minimum) / (maximum + minimum)",
            "direction_frame": "elastic-tensor Cartesian frame",
            "extrema_scope": "sampled grid",
            "notes": [
                "Extrema are selected from the sampled angular grid and are not "
                "continuous optimizations.",
                "Directions are expressed in the Cartesian frame of the input "
                "elastic tensor.",
            ],
        },
    )


def phase_derived_properties_table(result: SeismicResult) -> ReportTable:
    """Build shear splitting and phase-velocity ratio diagnostics."""
    phase = result.field.phase
    i_s2 = MODE_INDEX[WaveMode.V_S2]
    i_s1 = MODE_INDEX[WaveMode.V_S1]
    i_p = MODE_INDEX[WaveMode.V_P]
    speeds = phase.phase_speeds
    valid = phase.valid_mask

    derived: list[tuple[str, NDArray[np.float64], NDArray[np.bool_], str, bool]] = []
    shear_mask = valid[:, i_s1] & valid[:, i_s2]
    delta_s = speeds[:, i_s1] - speeds[:, i_s2]
    shear_anisotropy = _safe_divide(
        200.0 * delta_s,
        speeds[:, i_s1] + speeds[:, i_s2],
    )
    derived.extend(
        [
            ("V_S1 - V_S2", delta_s, shear_mask, "km s^-1", False),
            (
                "Directional S-wave anisotropy",
                shear_anisotropy,
                shear_mask,
                "%",
                False,
            ),
        ]
    )
    for label, denominator_index in (
        ("V_P / V_S1", i_s1),
        ("V_P / V_S2", i_s2),
    ):
        mask = valid[:, i_p] & valid[:, denominator_index]
        ratio = _safe_divide(speeds[:, i_p], speeds[:, denominator_index])
        derived.append((label, ratio, mask, "dimensionless", True))

    rows: list[list[object]] = []
    for label, values, mask, unit, report_range_anisotropy in derived:
        variation = sampled_variation(values, mask, phase.directions)
        eligible = values[mask & np.isfinite(values)]
        mean = float(np.mean(eligible)) if eligible.size else float("nan")
        percentile = (
            float(np.percentile(eligible, 95.0)) if eligible.size else float("nan")
        )
        range_anisotropy = (
            _format_number(variation.anisotropy_percent)
            if report_range_anisotropy
            else "not applicable"
        )
        rows.append(
            [
                label,
                unit,
                _format_number(variation.minimum.value),
                _format_direction(
                    variation.minimum.direction,
                    all_directions=variation.minimum.all_directions,
                ),
                _format_optional_number(variation.minimum.theta_degree),
                _format_optional_number(variation.minimum.phi_degree),
                _format_number(variation.maximum.value),
                _format_direction(
                    variation.maximum.direction,
                    all_directions=variation.maximum.all_directions,
                ),
                _format_optional_number(variation.maximum.theta_degree),
                _format_optional_number(variation.maximum.phi_degree),
                _format_number(mean),
                _format_number(percentile),
                range_anisotropy,
            ]
        )
    return ReportTable(
        title="Derived phase-velocity diagnostics",
        columns=[
            "Quantity",
            "Unit",
            "Minimum",
            "Minimum wave normal",
            "Theta min / degree",
            "Phi min / degree",
            "Maximum",
            "Maximum wave normal",
            "Theta max / degree",
            "Phi max / degree",
            "Mean",
            "95th percentile",
            "Range anisotropy / %",
        ],
        rows=rows,
        metadata={
            "directional_s_anisotropy_definition": (
                "200 * (V_S1 - V_S2) / (V_S1 + V_S2)"
            ),
            "range_anisotropy_definition": (
                "200 * (maximum - minimum) / (maximum + minimum)"
            ),
            "extrema_scope": "sampled grid",
            "notes": [
                "Range anisotropy is reported for velocity ratios only; it is "
                "not applied to splitting or to a quantity already expressed "
                "as an anisotropy percentage."
            ],
        },
    )


def group_velocity_extrema_table(result: SeismicResult) -> ReportTable:
    """Build sampled group-speed extrema and associated ray directions."""
    assert result.field.group is not None
    rows: list[list[object]] = []
    group = result.field.group
    for mode in reversed(MODE_ORDER):
        field = group.for_mode(mode)
        mask = field.valid_mask & field.resolved_mask
        variation = sampled_variation(
            field.group_speeds,
            mask,
            group.directions,
            ray_directions=field.ray_directions,
        )
        rows.append(_variation_row(MODE_SYMBOLS[mode], variation, include_ray=True))
    return ReportTable(
        title="Sampled group-velocity extrema",
        columns=_variation_columns(include_ray=True),
        rows=rows,
        metadata={
            "unit": "km s^-1",
            "anisotropy_definition": "200 * (maximum - minimum) / (maximum + minimum)",
            "wave_normal_frame": "elastic-tensor Cartesian frame",
            "ray_direction_frame": "elastic-tensor Cartesian frame",
            "extrema_scope": "resolved sampled modes",
            "notes": [
                "Group-speed extrema report both the wave normal and the energy-flow "
                "direction.",
                "Extrema are selected from the sampled grid.",
            ],
        },
    )


def power_flow_table(result: SeismicResult) -> ReportTable:
    """Build power-flow angle extrema and distribution statistics."""
    assert result.field.group is not None
    rows: list[list[object]] = []
    group = result.field.group
    for mode in reversed(MODE_ORDER):
        field = group.for_mode(mode)
        mask = field.valid_mask & field.resolved_mask
        values_degree = np.degrees(field.power_flow_angles)
        variation = sampled_variation(
            values_degree,
            mask,
            group.directions,
            ray_directions=field.ray_directions,
        )
        eligible = values_degree[mask & np.isfinite(values_degree)]
        if eligible.size:
            mean = float(np.mean(eligible))
            rms = float(np.sqrt(np.mean(np.square(eligible))))
            percentile = float(np.percentile(eligible, 95.0))
        else:
            mean = rms = percentile = float("nan")
        rows.append(
            [
                MODE_SYMBOLS[mode],
                _format_number(variation.minimum.value),
                _format_direction(variation.minimum.direction),
                _format_number(variation.maximum.value),
                _format_direction(variation.maximum.direction),
                _format_direction(variation.maximum.ray_direction),
                _format_number(mean),
                _format_number(rms),
                _format_number(percentile),
            ]
        )
    return ReportTable(
        title="Power-flow angle statistics",
        columns=[
            "Mode",
            "Minimum / degree",
            "Minimum wave normal",
            "Maximum / degree",
            "Maximum wave normal",
            "Ray direction at maximum",
            "Mean / degree",
            "RMS / degree",
            "95th percentile / degree",
        ],
        rows=rows,
        metadata={
            "definition": "angle between wave normal and group-velocity direction",
            "extrema_scope": "resolved sampled modes",
        },
    )


def enhancement_table(result: SeismicResult) -> ReportTable:
    """Build logarithmic enhancement extrema and caustic diagnostics."""
    assert result.field.enhancement is not None
    enhancement = result.field.enhancement
    rows: list[list[object]] = []
    for mode in reversed(MODE_ORDER):
        index = MODE_INDEX[mode]
        field = enhancement.for_mode(mode)
        mask = field.valid_mask & field.resolved_mask & field.finite_mask
        variation = sampled_variation(
            field.log10_enhancement,
            mask,
            enhancement.directions,
            ray_directions=field.group.ray_directions,
        )
        rows.append(
            [
                MODE_SYMBOLS[mode],
                _format_number(variation.minimum.value),
                _format_direction(variation.minimum.direction),
                _format_number(variation.maximum.value),
                _format_direction(variation.maximum.direction),
                _format_direction(variation.maximum.ray_direction),
                int(np.count_nonzero(mask)),
                int(np.count_nonzero(enhancement.caustic_candidate_mask[:, index])),
                int(
                    np.count_nonzero(
                        enhancement.valid_mask[:, index]
                        & ~enhancement.finite_mask[:, index]
                    )
                ),
            ]
        )
    return ReportTable(
        title="Enhancement and focusing diagnostics",
        columns=[
            "Mode",
            "Minimum log10(A)",
            "Minimum wave normal",
            "Maximum log10(A)",
            "Maximum wave normal",
            "Ray direction at maximum",
            "Finite resolved points",
            "Caustic candidates",
            "Non-finite points",
        ],
        rows=rows,
        metadata={
            "enhancement_representation": "log10(A)",
            "extrema_scope": "finite resolved sampled modes",
        },
    )


def diagnostics_table(result: SeismicResult) -> ReportTable:
    """Build a compact table of numerical and physical diagnostics."""
    field = result.field
    phase = field.phase
    rows: list[list[object]] = [
        ["Total sampled directions", field.n_points],
        ["Total phase-mode values", phase.valid_mask.size],
        ["Invalid phase-mode values", int(np.count_nonzero(~phase.valid_mask))],
        ["Clamped phase eigenvalues", int(np.count_nonzero(phase.clamped_mask))],
        ["Degenerate phase-mode values", int(np.count_nonzero(phase.degeneracy_mask))],
        [
            "Shear acoustic-axis candidate points",
            int(np.count_nonzero(phase.shear_axis_candidate_mask)),
        ],
        [
            "V_S1-V_P degeneracy candidate points",
            int(np.count_nonzero(phase.upper_axis_candidate_mask)),
        ],
    ]
    if field.group is not None:
        rows.extend(
            [
                [
                    "Invalid group-mode values",
                    int(np.count_nonzero(~field.group.valid_mask)),
                ],
                [
                    "Unresolved group-mode values",
                    int(
                        np.count_nonzero(
                            field.group.valid_mask & ~field.group.resolved_mask
                        )
                    ),
                ],
            ]
        )
    if field.enhancement is not None:
        rows.extend(
            [
                [
                    "Caustic candidate values",
                    int(np.count_nonzero(field.enhancement.caustic_candidate_mask)),
                ],
                [
                    "Non-finite enhancement values",
                    int(
                        np.count_nonzero(
                            field.enhancement.valid_mask
                            & ~field.enhancement.finite_mask
                        )
                    ),
                ],
            ]
        )
    if field.tracking is not None:
        tracking = field.tracking
        rows.extend(
            [
                [
                    "Polarization sign alignments",
                    int(np.count_nonzero(tracking.sign_flip_mask)),
                ],
                [
                    "Shear-branch exchanges",
                    int(np.count_nonzero(tracking.shear_swap_mask)),
                ],
                [
                    "Ambiguous shear permutations",
                    int(np.count_nonzero(tracking.shear_permutation_ambiguous_mask)),
                ],
                [
                    "Shear-subspace rotations",
                    int(np.count_nonzero(tracking.shear_subspace_rotation_mask)),
                ],
            ]
        )
    return ReportTable(
        title="Sampling and numerical diagnostics",
        columns=["Diagnostic", "Count"],
        rows=rows,
    )


def acoustic_axis_candidates_table(
    result: SeismicResult,
    *,
    maximum_rows: int = 24,
) -> ReportTable:
    """Build a table of the best sampled adjacent-mode degeneracy candidates.

    Parameters
    ----------
    result : SeismicResult
        Complete sampled seismic result.
    maximum_rows : int, optional
        Maximum number of unique antipodal candidate axes displayed.

    Returns
    -------
    ReportTable
        Candidate directions sorted by relative eigenvalue separation.
    """
    phase = result.field.phase
    entries: list[tuple[float, float, ModePair, NDArray[np.float64], int]] = []
    for pair in MODE_PAIR_ORDER:
        pair_index = MODE_PAIR_INDEX[pair]
        mask = phase.pair_degeneracy_mask[:, pair_index]
        indexes = np.flatnonzero(mask)
        grouped: dict[tuple[float, float, float], list[int]] = {}
        for index in indexes:
            axis = _canonical_axis(phase.directions[index])
            key = tuple(np.round(axis, 10))
            grouped.setdefault(key, []).append(int(index))
        for key, members in grouped.items():
            best = min(
                members,
                key=lambda item: phase.pair_relative_eigenvalue_gaps[item, pair_index],
            )
            entries.append(
                (
                    float(phase.pair_relative_eigenvalue_gaps[best, pair_index]),
                    float(phase.pair_eigenvalue_gaps[best, pair_index]),
                    pair,
                    np.asarray(key, dtype=float),
                    len(members),
                )
            )
    entries.sort(key=lambda item: (item[0], item[2].value, tuple(item[3])))
    rows: list[list[object]] = []
    for relative_gap, absolute_gap, pair, direction, multiplicity in entries[
        :maximum_rows
    ]:
        theta, phi = _angles_degree(direction)
        rows.append(
            [
                MODE_PAIR_SYMBOLS[pair],
                _format_direction(direction),
                _format_number(theta),
                _format_number(phi),
                f"{absolute_gap:.8e}",
                f"{relative_gap:.8e}",
                multiplicity,
            ]
        )
    return ReportTable(
        title="Sampled acoustic-axis candidates",
        columns=[
            "Mode pair",
            "Antipodal axis [x, y, z]",
            "Theta / degree",
            "Phi / degree",
            "Eigenvalue gap / km^2 s^-2",
            "Relative gap",
            "Equivalent samples",
        ],
        rows=rows,
        metadata={
            "total_unique_candidates": len(entries),
            "displayed_candidates": min(len(entries), maximum_rows),
            "truncated": len(entries) > maximum_rows,
            "candidate_definition": "sampled adjacent-mode gap below the solver tolerance",
            "direction_frame": "elastic-tensor Cartesian frame",
            "notes": [
                "Candidates are sampled near-degeneracies, not refined acoustic-axis "
                "solutions.",
                "Antipodal directions are reported as one physical axis.",
            ],
        },
    )


def sampled_variation(
    values: NDArray[np.float64],
    mask: NDArray[np.bool_],
    directions: NDArray[np.float64],
    *,
    ray_directions: NDArray[np.float64] | None = None,
) -> SampledVariation:
    """Return sampled extrema and symmetric percentage anisotropy.

    Parameters
    ----------
    values : ndarray
        Scalar values with shape ``(n_points,)``.
    mask : ndarray
        Boolean eligibility mask with shape ``(n_points,)``.
    directions : ndarray
        Wave-normal directions with shape ``(n_points, 3)``.
    ray_directions : ndarray or None, optional
        Energy-flow directions associated with the scalar values.

    Returns
    -------
    SampledVariation
        Extrema, anisotropy percentage, and maximum-to-minimum ratio.

    Raises
    ------
    ValueError
        If array shapes are inconsistent.
    """
    array = np.asarray(values, dtype=float)
    eligible = np.asarray(mask, dtype=bool) & np.isfinite(array)
    if array.ndim != 1 or eligible.shape != array.shape:
        raise ValueError("values and mask must be one-dimensional and equally sized.")
    if np.asarray(directions).shape != (array.size, 3):
        raise ValueError("directions must have shape (n_points, 3).")
    if ray_directions is not None and np.asarray(ray_directions).shape != (
        array.size,
        3,
    ):
        raise ValueError("ray_directions must have shape (n_points, 3).")
    if not np.any(eligible):
        unavailable = SampledExtremum(
            value=float("nan"),
            direction=None,
            ray_direction=None,
            theta_degree=None,
            phi_degree=None,
            multiplicity=0,
            all_directions=False,
        )
        return SampledVariation(
            minimum=unavailable,
            maximum=unavailable,
            anisotropy_percent=float("nan"),
            maximum_to_minimum=float("nan"),
        )

    minimum = _sampled_extremum(
        array,
        eligible,
        directions,
        kind="minimum",
        ray_directions=ray_directions,
    )
    maximum = _sampled_extremum(
        array,
        eligible,
        directions,
        kind="maximum",
        ray_directions=ray_directions,
    )
    denominator = maximum.value + minimum.value
    anisotropy = (
        200.0 * (maximum.value - minimum.value) / denominator
        if denominator != 0.0
        else float("nan")
    )
    ratio = maximum.value / minimum.value if minimum.value != 0.0 else float("inf")
    return SampledVariation(
        minimum=minimum,
        maximum=maximum,
        anisotropy_percent=float(anisotropy),
        maximum_to_minimum=float(ratio),
    )


def _sampled_extremum(
    values: NDArray[np.float64],
    eligible: NDArray[np.bool_],
    directions: NDArray[np.float64],
    *,
    kind: Literal["minimum", "maximum"],
    ray_directions: NDArray[np.float64] | None,
) -> SampledExtremum:
    """Return one deterministic extremum and its sampled multiplicity."""
    selected = values[eligible]
    target = float(np.min(selected) if kind == "minimum" else np.max(selected))
    absolute_tolerance = max(1.0e-12, abs(target) * 1.0e-10)
    tied = eligible & np.isclose(
        values,
        target,
        rtol=1.0e-9,
        atol=absolute_tolerance,
    )
    indexes = np.flatnonzero(tied)
    eligible_indexes = np.flatnonzero(eligible)
    all_directions = indexes.size == eligible_indexes.size
    axes = {_axis_key(directions[index]) for index in indexes}
    multiplicity = len(axes)
    if all_directions:
        return SampledExtremum(
            value=target,
            direction=None,
            ray_direction=None,
            theta_degree=None,
            phi_degree=None,
            multiplicity=multiplicity,
            all_directions=True,
        )
    index = int(indexes[0])
    direction = np.array(directions[index], dtype=float, copy=True)
    direction.setflags(write=False)
    ray = None
    if ray_directions is not None and np.all(np.isfinite(ray_directions[index])):
        ray = np.array(ray_directions[index], dtype=float, copy=True)
        ray.setflags(write=False)
    theta, phi = _angles_degree(direction)
    return SampledExtremum(
        value=target,
        direction=direction,
        ray_direction=ray,
        theta_degree=theta,
        phi_degree=phi,
        multiplicity=multiplicity,
        all_directions=False,
    )


def _variation_columns(*, include_ray: bool) -> list[str]:
    """Return common columns for sampled-variation tables."""
    columns = [
        "Mode",
        "Minimum",
        "Minimum wave normal",
        "Theta min / degree",
        "Phi min / degree",
    ]
    if include_ray:
        columns.append("Ray direction at minimum")
    columns.extend(
        [
            "Maximum",
            "Maximum wave normal",
            "Theta max / degree",
            "Phi max / degree",
        ]
    )
    if include_ray:
        columns.append("Ray direction at maximum")
    columns.extend(
        [
            "Anisotropy / %",
            "Maximum / minimum",
            "Equivalent minima",
            "Equivalent maxima",
        ]
    )
    return columns


def _variation_row(
    label: str,
    variation: SampledVariation,
    *,
    include_ray: bool = False,
) -> list[object]:
    """Build one report row from a sampled variation."""
    return [label, *_variation_row_cells(variation, include_ray=include_ray)]


def _variation_row_cells(
    variation: SampledVariation,
    *,
    include_ray: bool = False,
) -> list[object]:
    """Build report cells for one sampled variation."""
    minimum = variation.minimum
    maximum = variation.maximum
    cells: list[object] = [
        _format_number(minimum.value),
        _format_direction(minimum.direction, all_directions=minimum.all_directions),
        _format_optional_number(minimum.theta_degree),
        _format_optional_number(minimum.phi_degree),
    ]
    if include_ray:
        cells.append(_format_direction(minimum.ray_direction))
    cells.extend(
        [
            _format_number(maximum.value),
            _format_direction(maximum.direction, all_directions=maximum.all_directions),
            _format_optional_number(maximum.theta_degree),
            _format_optional_number(maximum.phi_degree),
        ]
    )
    if include_ray:
        cells.append(_format_direction(maximum.ray_direction))
    cells.extend(
        [
            _format_number(variation.anisotropy_percent),
            _format_number(variation.maximum_to_minimum),
            minimum.multiplicity,
            maximum.multiplicity,
        ]
    )
    return cells


def _safe_divide(
    numerator: NDArray[np.float64],
    denominator: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Return a division with non-zero finite denominators only."""
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    result = np.full(np.broadcast_shapes(numerator.shape, denominator.shape), np.nan)
    np.divide(
        numerator,
        denominator,
        out=result,
        where=np.isfinite(denominator) & (denominator != 0.0),
    )
    return result


def _canonical_axis(direction: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return a deterministic representative of an antipodal Cartesian axis."""
    axis = np.asarray(direction, dtype=float).copy()
    for component in axis:
        if abs(component) > 1.0e-12:
            if component < 0.0:
                axis *= -1.0
            break
    axis[np.abs(axis) < 1.0e-12] = 0.0
    return axis


def _axis_key(direction: NDArray[np.float64]) -> tuple[float, float, float]:
    """Return a rounded hash key for one antipodal axis."""
    axis = np.round(_canonical_axis(direction), 10)
    return (float(axis[0]), float(axis[1]), float(axis[2]))


def _angles_degree(direction: NDArray[np.float64]) -> tuple[float, float]:
    """Return polar and azimuthal angles in degrees for one direction."""
    theta, phi = cartesian_to_spherical(np.asarray(direction, dtype=float))
    return float(np.degrees(theta)), float(np.degrees(phi))


def _format_direction(
    direction: NDArray[np.float64] | None,
    *,
    all_directions: bool = False,
) -> str:
    """Format a Cartesian direction or an unavailable-direction marker."""
    if all_directions:
        return "all sampled directions"
    if direction is None:
        return "unavailable"
    return "[{:.6f}, {:.6f}, {:.6f}]".format(*direction)


def _format_optional_number(value: float | None) -> str:
    """Format an optional floating-point report value."""
    return "unavailable" if value is None else _format_number(value)


def _format_number(value: float) -> str:
    """Format finite and non-finite report values consistently."""
    if np.isnan(value):
        return "unavailable"
    if np.isposinf(value):
        return "inf"
    if np.isneginf(value):
        return "-inf"
    return f"{value:.8f}"
