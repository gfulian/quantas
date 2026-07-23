# -*- coding: utf-8 -*-

"""Batched sampling of directional elastic properties.

This module evaluates directional elastic fields directly from the Cartesian
compliance tensor.  Young's modulus and linear compressibility are contracted
vectorially over arbitrary angular arrays.  The extrema of shear modulus and
Poisson's ratio over the transverse angle are obtained analytically from the
eigenvalues of a symmetric ``2 x 2`` quadratic form in the plane normal to the
longitudinal direction.

No interpolation, fitting, or local numerical optimization is used.  The
batched implementation is therefore deterministic and does not introduce
isolated optimizer artefacts into two- or three-dimensional property fields.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Literal, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .tensor import ElasticTensor


ElasticDirectionalProperty: TypeAlias = Literal[
    "young",
    "compressibility",
    "shear",
    "poisson",
]
ProgressCallback = Callable[[int, int], None]

SUPPORTED_DIRECTIONAL_PROPERTIES: tuple[ElasticDirectionalProperty, ...] = (
    "young",
    "compressibility",
    "shear",
    "poisson",
)


@dataclass(frozen=True, slots=True)
class ElasticSamplingDiagnostics:
    """Numerical diagnostics for a sampled directional field.

    Parameters
    ----------
    minimum_longitudinal_compliance : float or None
        Smallest value of ``S_ijkl n_i n_j n_k n_l`` over the sampled
        directions.  A positive value is required for a mechanically stable
        tensor.
    minimum_transverse_compliance : float or None
        Smallest eigenvalue of the compliance quadratic form projected onto
        the transverse plane.  It controls the maximum sampled shear modulus.
    shear_degeneracy_count : int
        Number of directions for which the two transverse shear eigenvalues
        are numerically degenerate.
    poisson_degeneracy_count : int
        Number of directions for which the two transverse Poisson quadratic
        forms are numerically degenerate.
    """

    minimum_longitudinal_compliance: float | None = None
    minimum_transverse_compliance: float | None = None
    shear_degeneracy_count: int = 0
    poisson_degeneracy_count: int = 0


@dataclass(frozen=True, slots=True)
class ElasticDirectionalField:
    """Batched directional elastic properties on paired angular arrays.

    Parameters
    ----------
    theta, phi : ndarray
        Paired polar and azimuthal angles in radians.  Their shape is retained
        by all sampled property arrays.
    directions : ndarray
        Cartesian unit directions with shape ``theta.shape + (3,)``.
    young_modulus : ndarray or None
        Young's modulus in GPa.
    linear_compressibility : ndarray or None
        Signed linear compressibility in TPa^-1.
    shear_minimum, shear_maximum : ndarray or None
        Exact extrema of shear modulus in GPa over the transverse direction.
    poisson_minimum, poisson_maximum : ndarray or None
        Exact extrema of Poisson's ratio over the transverse direction.
    diagnostics : ElasticSamplingDiagnostics
        Numerical denominator and degeneracy diagnostics.
    batch_size : int
        Maximum number of directions evaluated in one bounded NumPy batch.
    """

    theta: NDArray[np.float64]
    phi: NDArray[np.float64]
    directions: NDArray[np.float64]
    young_modulus: NDArray[np.float64] | None
    linear_compressibility: NDArray[np.float64] | None
    shear_minimum: NDArray[np.float64] | None
    shear_maximum: NDArray[np.float64] | None
    poisson_minimum: NDArray[np.float64] | None
    poisson_maximum: NDArray[np.float64] | None
    diagnostics: ElasticSamplingDiagnostics
    batch_size: int


@dataclass(frozen=True, slots=True)
class _DirectionalBatch:
    """Numerical values and diagnostics for one bounded direction batch."""

    directions: NDArray[np.float64]
    values: dict[str, NDArray[np.float64]]
    minimum_longitudinal_compliance: float | None
    minimum_transverse_compliance: float | None
    shear_degeneracy_count: int
    poisson_degeneracy_count: int


@dataclass(frozen=True, slots=True)
class ExactTransverseExtrema:
    """Exact extrema of one transverse directional property.

    Parameters
    ----------
    minimum, maximum : float
        Extremal property values.
    minimum_angle, maximum_angle : float
        Transverse angles in radians in ``[0, pi)``.  At a numerical
        degeneracy both angles are deterministically set to zero because every
        transverse direction is equivalent.
    """

    minimum: float
    maximum: float
    minimum_angle: float
    maximum_angle: float


def sample_elastic_directional_field(
    tensor: ElasticTensor,
    theta: ArrayLike,
    phi: ArrayLike,
    *,
    properties: Iterable[ElasticDirectionalProperty] = (
        "young",
        "compressibility",
        "shear",
        "poisson",
    ),
    batch_size: int = 65536,
    progress_callback: ProgressCallback | None = None,
) -> ElasticDirectionalField:
    """Sample selected directional elastic properties in bounded batches.

    Parameters
    ----------
    tensor : ElasticTensor
        Mechanically stable elastic tensor.
    theta, phi : array_like
        Paired, broadcast-compatible polar and azimuthal angles in radians.
    properties : iterable of str, optional
        Any subset of ``young``, ``compressibility``, ``shear`` and
        ``poisson``.
    batch_size : int, optional
        Maximum number of directions evaluated in one NumPy batch.
    progress_callback : callable or None, optional
        Callback receiving ``(current, total)`` in units of completed
        direction-property evaluations after each batch.

    Returns
    -------
    ElasticDirectionalField
        Read-only arrays with the broadcast angular shape.

    Raises
    ------
    ValueError
        If angles, properties, or ``batch_size`` are invalid, or if a required
        compliance denominator is non-finite or non-positive.
    """
    selected = _normalize_properties(properties)
    size = _validate_batch_size(batch_size)
    polar, azimuth = _broadcast_angles(theta, phi)
    shape = polar.shape
    flat_polar = polar.reshape(-1)
    flat_azimuth = azimuth.reshape(-1)
    count = int(flat_polar.size)
    total = count * len(selected)
    arrays = _allocate_field_arrays(count, selected)

    minimum_longitudinal: float | None = None
    minimum_transverse: float | None = None
    shear_degeneracy_count = 0
    poisson_degeneracy_count = 0

    for start in range(0, count, size):
        stop = min(start + size, count)
        batch = _sample_direction_batch(
            tensor.compliance_tensor,
            flat_polar[start:stop],
            flat_azimuth[start:stop],
            selected,
        )
        arrays["directions"][start:stop] = batch.directions
        for name, values in batch.values.items():
            arrays[name][start:stop] = values
        minimum_longitudinal = _optional_minimum(
            minimum_longitudinal,
            batch.minimum_longitudinal_compliance,
        )
        minimum_transverse = _optional_minimum(
            minimum_transverse,
            batch.minimum_transverse_compliance,
        )
        shear_degeneracy_count += batch.shear_degeneracy_count
        poisson_degeneracy_count += batch.poisson_degeneracy_count
        _report_progress(stop * len(selected), total, progress_callback)

    diagnostics = ElasticSamplingDiagnostics(
        minimum_longitudinal_compliance=minimum_longitudinal,
        minimum_transverse_compliance=minimum_transverse,
        shear_degeneracy_count=shear_degeneracy_count,
        poisson_degeneracy_count=poisson_degeneracy_count,
    )
    return _build_directional_field(
        polar,
        azimuth,
        arrays,
        shape,
        diagnostics,
        size,
    )


def _sample_direction_batch(
    compliance: NDArray[np.float64],
    theta: NDArray[np.float64],
    phi: NDArray[np.float64],
    selected: tuple[ElasticDirectionalProperty, ...],
) -> _DirectionalBatch:
    """Evaluate one independent batch of directions."""
    directions, first, second = _direction_frame(theta, phi)
    values: dict[str, NDArray[np.float64]] = {}
    longitudinal: NDArray[np.float64] | None = None
    minimum_longitudinal: float | None = None
    minimum_transverse: float | None = None
    shear_degeneracy_count = 0
    poisson_degeneracy_count = 0

    if "young" in selected or "poisson" in selected:
        longitudinal = np.einsum(
            "...i,...j,...k,...l,ijkl->...",
            directions,
            directions,
            directions,
            directions,
            compliance,
            optimize=True,
        )
        _require_positive_finite(longitudinal, name="longitudinal compliance")
        minimum_longitudinal = float(np.min(longitudinal))

    if "young" in selected:
        assert longitudinal is not None
        values["young_modulus"] = np.asarray(1.0 / longitudinal, dtype=float)

    if "compressibility" in selected:
        contracted = np.einsum("ijkk->ij", compliance, optimize=True)
        compressibility = 1000.0 * np.einsum(
            "...i,...j,ij->...",
            directions,
            directions,
            contracted,
            optimize=True,
        )
        _require_finite(compressibility, name="linear compressibility")
        values["linear_compressibility"] = np.asarray(
            compressibility,
            dtype=float,
        )

    if "shear" in selected:
        eigenvalue_minimum, eigenvalue_maximum = _project_transverse_quadratic_form(
            compliance,
            directions,
            first,
            second,
            contraction="shear",
        )
        _require_positive_finite(
            eigenvalue_minimum,
            name="transverse compliance",
        )
        minimum_transverse = float(np.min(eigenvalue_minimum))
        values["shear_minimum"] = np.asarray(
            1.0 / (4.0 * eigenvalue_maximum),
            dtype=float,
        )
        values["shear_maximum"] = np.asarray(
            1.0 / (4.0 * eigenvalue_minimum),
            dtype=float,
        )
        shear_degeneracy_count = _degeneracy_count(
            eigenvalue_minimum,
            eigenvalue_maximum,
        )

    if "poisson" in selected:
        assert longitudinal is not None
        eigenvalue_minimum, eigenvalue_maximum = _project_transverse_quadratic_form(
            compliance,
            directions,
            first,
            second,
            contraction="poisson",
        )
        poisson_minimum = -eigenvalue_maximum / longitudinal
        poisson_maximum = -eigenvalue_minimum / longitudinal
        _require_finite(poisson_minimum, name="minimum Poisson ratio")
        _require_finite(poisson_maximum, name="maximum Poisson ratio")
        values["poisson_minimum"] = np.asarray(poisson_minimum, dtype=float)
        values["poisson_maximum"] = np.asarray(poisson_maximum, dtype=float)
        poisson_degeneracy_count = _degeneracy_count(
            eigenvalue_minimum,
            eigenvalue_maximum,
        )

    return _DirectionalBatch(
        directions=np.asarray(directions, dtype=float),
        values=values,
        minimum_longitudinal_compliance=minimum_longitudinal,
        minimum_transverse_compliance=minimum_transverse,
        shear_degeneracy_count=shear_degeneracy_count,
        poisson_degeneracy_count=poisson_degeneracy_count,
    )


def _allocate_field_arrays(
    count: int,
    selected: tuple[ElasticDirectionalProperty, ...],
) -> dict[str, NDArray[np.float64]]:
    """Allocate flat mutable arrays for a complete directional field."""
    arrays: dict[str, NDArray[np.float64]] = {
        "directions": np.empty((count, 3), dtype=float)
    }
    if "young" in selected:
        arrays["young_modulus"] = np.empty(count, dtype=float)
    if "compressibility" in selected:
        arrays["linear_compressibility"] = np.empty(count, dtype=float)
    if "shear" in selected:
        arrays["shear_minimum"] = np.empty(count, dtype=float)
        arrays["shear_maximum"] = np.empty(count, dtype=float)
    if "poisson" in selected:
        arrays["poisson_minimum"] = np.empty(count, dtype=float)
        arrays["poisson_maximum"] = np.empty(count, dtype=float)
    return arrays


def _build_directional_field(
    theta: NDArray[np.float64],
    phi: NDArray[np.float64],
    arrays: dict[str, NDArray[np.float64]],
    shape: tuple[int, ...],
    diagnostics: ElasticSamplingDiagnostics,
    batch_size: int,
) -> ElasticDirectionalField:
    """Reshape mutable buffers into one read-only directional field."""

    def optional(name: str) -> NDArray[np.float64] | None:
        values = arrays.get(name)
        return None if values is None else _read_only(values.reshape(shape))

    return ElasticDirectionalField(
        theta=_read_only(theta),
        phi=_read_only(phi),
        directions=_read_only(arrays["directions"].reshape(shape + (3,))),
        young_modulus=optional("young_modulus"),
        linear_compressibility=optional("linear_compressibility"),
        shear_minimum=optional("shear_minimum"),
        shear_maximum=optional("shear_maximum"),
        poisson_minimum=optional("poisson_minimum"),
        poisson_maximum=optional("poisson_maximum"),
        diagnostics=diagnostics,
        batch_size=batch_size,
    )


def _optional_minimum(
    current: float | None,
    candidate: float | None,
) -> float | None:
    """Combine optional diagnostic minima."""
    if candidate is None:
        return current
    if current is None:
        return candidate
    return min(current, candidate)


def _validate_batch_size(batch_size: int) -> int:
    """Return a validated positive integer batch size."""
    if isinstance(batch_size, bool) or int(batch_size) != batch_size:
        raise ValueError("batch_size must be a positive integer.")
    size = int(batch_size)
    if size <= 0:
        raise ValueError("batch_size must be a positive integer.")
    return size


def exact_transverse_extrema(
    tensor: ElasticTensor,
    theta: float,
    phi: float,
    *,
    kind: Literal["shear", "poisson"],
) -> ExactTransverseExtrema:
    """Return exact transverse extrema and deterministic extremal angles.

    Parameters
    ----------
    tensor : ElasticTensor
        Mechanically stable elastic tensor.
    theta, phi : float
        Longitudinal polar and azimuthal angles in radians.
    kind : {"shear", "poisson"}
        Property optimized over the transverse plane.

    Returns
    -------
    ExactTransverseExtrema
        Exact extremal values and transverse angles.

    Raises
    ------
    ValueError
        If ``kind`` is unsupported or a compliance denominator is invalid.
    """
    if kind not in {"shear", "poisson"}:
        raise ValueError(f"unsupported transverse property: {kind}")

    polar, azimuth = _broadcast_angles(theta, phi)
    directions, first, second = _direction_frame(polar, azimuth)
    direction = directions.reshape(3)
    basis = np.stack((first.reshape(3), second.reshape(3)), axis=0)
    compliance = tensor.compliance_tensor

    if kind == "shear":
        matrix = np.einsum(
            "ijkl,i,k,aj,bl->ab",
            compliance,
            direction,
            direction,
            basis,
            basis,
            optimize=True,
        )
        eigenvalues, eigenvectors = np.linalg.eigh(_symmetrize_2x2(matrix))
        _require_positive_finite(eigenvalues, name="transverse compliance")
        minimum = float(1.0 / (4.0 * eigenvalues[1]))
        maximum = float(1.0 / (4.0 * eigenvalues[0]))
        minimum_vector = eigenvectors[:, 1]
        maximum_vector = eigenvectors[:, 0]
    else:
        matrix = np.einsum(
            "ijkl,i,j,ak,bl->ab",
            compliance,
            direction,
            direction,
            basis,
            basis,
            optimize=True,
        )
        matrix = _symmetrize_2x2(matrix)
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        denominator = float(
            np.einsum(
                "ijkl,i,j,k,l->",
                compliance,
                direction,
                direction,
                direction,
                direction,
                optimize=True,
            )
        )
        _require_positive_finite(
            np.asarray([denominator]),
            name="longitudinal compliance",
        )
        minimum = float(-eigenvalues[1] / denominator)
        maximum = float(-eigenvalues[0] / denominator)
        minimum_vector = eigenvectors[:, 1]
        maximum_vector = eigenvectors[:, 0]

    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    degenerate = bool(
        abs(float(eigenvalues[1] - eigenvalues[0]))
        <= 128.0 * np.finfo(float).eps * scale
    )
    if degenerate:
        minimum_angle = maximum_angle = 0.0
    else:
        minimum_angle = _transverse_angle(minimum_vector)
        maximum_angle = _transverse_angle(maximum_vector)

    return ExactTransverseExtrema(
        minimum=minimum,
        maximum=maximum,
        minimum_angle=minimum_angle,
        maximum_angle=maximum_angle,
    )


def _normalize_properties(
    properties: Iterable[ElasticDirectionalProperty],
) -> tuple[ElasticDirectionalProperty, ...]:
    """Return unique supported properties while preserving input order."""
    selected: list[ElasticDirectionalProperty] = []
    for property_name in properties:
        if property_name not in SUPPORTED_DIRECTIONAL_PROPERTIES:
            raise ValueError(f"unsupported elasticity property: {property_name}")
        if property_name not in selected:
            selected.append(property_name)
    return tuple(selected)


def _broadcast_angles(
    theta: ArrayLike,
    phi: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Validate and broadcast paired angular arrays."""
    polar = np.asarray(theta, dtype=float)
    azimuth = np.asarray(phi, dtype=float)
    if not np.all(np.isfinite(polar)) or not np.all(np.isfinite(azimuth)):
        raise ValueError("theta and phi must contain finite values.")
    try:
        polar, azimuth = np.broadcast_arrays(polar, azimuth)
    except ValueError as exc:
        raise ValueError("theta and phi must be broadcast-compatible.") from exc
    return np.asarray(polar, dtype=float), np.asarray(azimuth, dtype=float)


def _direction_frame(
    theta: NDArray[np.float64],
    phi: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Construct a longitudinal direction and orthonormal transverse basis."""
    sine_theta = np.sin(theta)
    cosine_theta = np.cos(theta)
    sine_phi = np.sin(phi)
    cosine_phi = np.cos(phi)

    threshold = 8.0 * np.finfo(float).eps
    sine_theta = np.where(np.abs(sine_theta) <= threshold, 0.0, sine_theta)
    cosine_theta = np.where(np.abs(cosine_theta) <= threshold, 0.0, cosine_theta)
    sine_phi = np.where(np.abs(sine_phi) <= threshold, 0.0, sine_phi)
    cosine_phi = np.where(np.abs(cosine_phi) <= threshold, 0.0, cosine_phi)

    direction = np.stack(
        (
            sine_theta * cosine_phi,
            sine_theta * sine_phi,
            cosine_theta,
        ),
        axis=-1,
    )
    first = np.stack(
        (
            cosine_theta * cosine_phi,
            cosine_theta * sine_phi,
            -sine_theta,
        ),
        axis=-1,
    )
    second = np.stack(
        (
            -sine_phi,
            cosine_phi,
            np.zeros_like(phi),
        ),
        axis=-1,
    )
    return direction, first, second


def _project_transverse_quadratic_form(
    compliance: NDArray[np.float64],
    directions: NDArray[np.float64],
    first: NDArray[np.float64],
    second: NDArray[np.float64],
    *,
    contraction: Literal["shear", "poisson"],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return ordered eigenvalues of a transverse ``2 x 2`` quadratic form."""
    if contraction == "shear":
        cartesian = np.einsum(
            "ijkl,...i,...k->...jl",
            compliance,
            directions,
            directions,
            optimize=True,
        )
    else:
        cartesian = np.einsum(
            "ijkl,...i,...j->...kl",
            compliance,
            directions,
            directions,
            optimize=True,
        )

    q00 = np.einsum(
        "...i,...ij,...j->...",
        first,
        cartesian,
        first,
        optimize=True,
    )
    q11 = np.einsum(
        "...i,...ij,...j->...",
        second,
        cartesian,
        second,
        optimize=True,
    )
    q01 = 0.5 * (
        np.einsum(
            "...i,...ij,...j->...",
            first,
            cartesian,
            second,
            optimize=True,
        )
        + np.einsum(
            "...i,...ij,...j->...",
            second,
            cartesian,
            first,
            optimize=True,
        )
    )
    center = 0.5 * (q00 + q11)
    radius = np.hypot(0.5 * (q00 - q11), q01)
    minimum = np.asarray(center - radius, dtype=float)
    maximum = np.asarray(center + radius, dtype=float)
    _require_finite(minimum, name=f"{contraction} transverse eigenvalue")
    _require_finite(maximum, name=f"{contraction} transverse eigenvalue")
    return minimum, maximum


def _require_positive_finite(values: NDArray[np.float64], *, name: str) -> None:
    """Reject invalid or non-positive compliance denominators."""
    _require_finite(values, name=name)
    scale = max(float(np.max(np.abs(values))), 1.0)
    tolerance = 128.0 * np.finfo(float).eps * scale
    if np.any(values <= tolerance):
        minimum = float(np.min(values))
        raise ValueError(
            f"{name} must be positive for directional sampling; "
            f"minimum value is {minimum:.16g}."
        )


def _require_finite(values: NDArray[np.float64], *, name: str) -> None:
    """Reject non-finite numerical results."""
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} contains non-finite values.")


def _degeneracy_count(
    minimum: NDArray[np.float64],
    maximum: NDArray[np.float64],
) -> int:
    """Count numerically degenerate transverse eigenvalue pairs."""
    scale = np.maximum(np.maximum(np.abs(minimum), np.abs(maximum)), 1.0)
    tolerance = 128.0 * np.finfo(float).eps * scale
    return int(np.count_nonzero(np.abs(maximum - minimum) <= tolerance))


def _symmetrize_2x2(matrix: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return the symmetric part of a two-dimensional quadratic form."""
    return np.asarray(0.5 * (matrix + matrix.T), dtype=float)


def _transverse_angle(vector: NDArray[np.float64]) -> float:
    """Return a transverse basis angle modulo the physical period ``pi``."""
    return float(np.arctan2(float(vector[1]), float(vector[0])) % np.pi)


def _read_only(array: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return a copied read-only floating-point array."""
    result = np.array(array, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _report_progress(
    current: int,
    total: int,
    callback: ProgressCallback | None,
) -> None:
    """Report numerical progress through an optional callback."""
    if callback is not None:
        callback(current, total)
