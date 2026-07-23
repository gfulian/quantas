# -*- coding: utf-8 -*-

"""Directional elastic properties and angular sampling utilities."""

from __future__ import annotations

import itertools as it
import math
from collections.abc import Callable, Sequence

import numpy as np

from .sampling import exact_transverse_extrema, sample_elastic_directional_field
from .tensor import ElasticTensor, OrthorhombicElasticTensor


ProgressCallback = Callable[[int, int], None]


def direction_vector(
    theta: float,
    phi: float,
    chi: float | None = None,
) -> np.ndarray:
    """Construct a unit direction from spherical angles.

    Parameters
    ----------
    theta : float
        Polar angle in radians.
    phi : float
        Azimuthal angle in radians.
    chi : float or None, optional
        Rotation defining a transverse measurement direction.

    Returns
    -------
    ndarray
        Cartesian unit vector.
    """
    theta = _as_scalar_angle(theta)
    phi = _as_scalar_angle(phi)
    chi = None if chi is None else _as_scalar_angle(chi)

    if chi is None:
        vector = [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    else:
        vector = [
            np.cos(theta) * np.cos(phi) * np.cos(chi) - np.sin(phi) * np.sin(chi),
            np.cos(theta) * np.sin(phi) * np.cos(chi) + np.cos(phi) * np.sin(chi),
            -np.sin(theta) * np.cos(chi),
        ]
    return np.asarray(vector, dtype=float)


def _as_scalar_angle(value: float | np.ndarray) -> float:
    """Normalize a scalar or one-element optimizer array to ``float``.

    Parameters
    ----------
    value : float or ndarray
        Angular coordinate supplied directly or by SciPy.

    Returns
    -------
    float
        Scalar angular coordinate.

    Raises
    ------
    ValueError
        If the supplied value does not contain exactly one element.
    """
    array = np.asarray(value, dtype=float)
    if array.size != 1:
        raise ValueError("An angular coordinate must contain exactly one value.")
    return float(array.reshape(-1)[0])


def young_modulus(tensor: ElasticTensor, angles: Sequence[float]) -> float:
    """Calculate Young's modulus along a direction.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor.
    angles : sequence of float
        Polar and azimuthal angles ``(theta, phi)``.

    Returns
    -------
    float
        Young's modulus in GPa.
    """
    if isinstance(tensor, OrthorhombicElasticTensor):
        return _orthorhombic_young_modulus(tensor, angles)

    axis = direction_vector(angles[0], angles[1])
    result = 0.0
    compliance = tensor.compliance_tensor
    for i, j, k, ell in it.product((0, 1, 2), repeat=4):
        result += compliance[i, j, k, ell] * axis[i] * axis[j] * axis[k] * axis[ell]
    return float(1.0 / result)


def linear_compressibility(
    tensor: ElasticTensor,
    angles: Sequence[float],
) -> float:
    """Calculate linear compressibility along a direction.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor.
    angles : sequence of float
        Polar and azimuthal angles ``(theta, phi)``.

    Returns
    -------
    float
        Linear compressibility in TPa^-1.
    """
    if isinstance(tensor, OrthorhombicElasticTensor):
        return _orthorhombic_linear_compressibility(tensor, angles)

    axis = direction_vector(angles[0], angles[1])
    result = 0.0
    compliance = tensor.compliance_tensor
    for i, j, k in it.product((0, 1, 2), repeat=3):
        result += axis[i] * axis[j] * compliance[i, j, k, k]
    return float(1000.0 * result)


def shear_modulus(tensor: ElasticTensor, angles: Sequence[float]) -> float:
    """Calculate shear modulus for a direction and transverse axis.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor.
    angles : sequence of float
        Angles ``(theta, phi, chi)``.

    Returns
    -------
    float
        Shear modulus in GPa.
    """
    if isinstance(tensor, OrthorhombicElasticTensor):
        return _orthorhombic_shear_modulus(tensor, angles)

    axis = direction_vector(angles[0], angles[1])
    transverse = direction_vector(angles[0], angles[1], angles[2])
    result = 0.0
    compliance = tensor.compliance_tensor
    for i, j, k, ell in it.product((0, 1, 2), repeat=4):
        result += (
            axis[i]
            * transverse[j]
            * axis[k]
            * transverse[ell]
            * compliance[i, j, k, ell]
        )
    return float(1.0 / (4.0 * result))


def poisson_ratio(tensor: ElasticTensor, angles: Sequence[float]) -> float:
    """Calculate Poisson's ratio for a direction and transverse axis.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor.
    angles : sequence of float
        Angles ``(theta, phi, chi)``.

    Returns
    -------
    float
        Dimensionless Poisson ratio.
    """
    if isinstance(tensor, OrthorhombicElasticTensor):
        return _orthorhombic_poisson_ratio(tensor, angles)

    axis = direction_vector(angles[0], angles[1])
    transverse = direction_vector(angles[0], angles[1], angles[2])
    numerator = 0.0
    denominator = 0.0
    compliance = tensor.compliance_tensor
    for i, j, k, ell in it.product((0, 1, 2), repeat=4):
        numerator += (
            axis[i]
            * axis[j]
            * transverse[k]
            * transverse[ell]
            * compliance[i, j, k, ell]
        )
        denominator += (
            axis[i] * axis[j] * axis[k] * axis[ell] * compliance[i, j, k, ell]
        )
    return float(-numerator / denominator)


def transverse_shear_extrema(
    tensor: ElasticTensor,
    angles: Sequence[float],
) -> tuple[float, float]:
    """Calculate exact transverse shear-modulus extrema.

    The transverse dependence is a two-dimensional quadratic form.  Its
    extrema are obtained algebraically from the projected compliance
    eigenvalues rather than by local numerical optimization.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor.
    angles : sequence of float
        Polar and azimuthal angles ``(theta, phi)``.

    Returns
    -------
    tuple of float
        Minimum and maximum shear modulus in GPa.
    """
    extrema = exact_transverse_extrema(
        tensor,
        float(angles[0]),
        float(angles[1]),
        kind="shear",
    )
    return extrema.minimum, extrema.maximum


def transverse_poisson_extrema(
    tensor: ElasticTensor,
    angles: Sequence[float],
) -> tuple[float, float, float]:
    """Calculate exact transverse Poisson-ratio extrema.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor.
    angles : sequence of float
        Polar and azimuthal angles ``(theta, phi)``.

    Returns
    -------
    tuple of float
        Negative part, positive minimum, and maximum used by polar plots.
    """
    extrema = exact_transverse_extrema(
        tensor,
        float(angles[0]),
        float(angles[1]),
        kind="poisson",
    )
    return (
        min(0.0, extrema.minimum),
        max(0.0, extrema.minimum),
        extrema.maximum,
    )


def sample_young_modulus(
    tensor: ElasticTensor,
    theta: Sequence[float],
    phi: Sequence[float],
    progress_callback: ProgressCallback | None = None,
) -> np.ndarray:
    """Sample Young's modulus on paired angular arrays in one batch."""
    theta_array, phi_array = _paired_angles(theta, phi)
    field = sample_elastic_directional_field(
        tensor,
        theta_array,
        phi_array,
        properties=("young",),
    )
    _replay_progress(theta_array.size, progress_callback)
    assert field.young_modulus is not None
    return np.array(field.young_modulus, dtype=float, copy=True)


def sample_linear_compressibility(
    tensor: ElasticTensor,
    theta: Sequence[float],
    phi: Sequence[float],
    progress_callback: ProgressCallback | None = None,
) -> np.ndarray:
    """Sample signed linear compressibility in one batch."""
    theta_array, phi_array = _paired_angles(theta, phi)
    field = sample_elastic_directional_field(
        tensor,
        theta_array,
        phi_array,
        properties=("compressibility",),
    )
    _replay_progress(theta_array.size, progress_callback)
    assert field.linear_compressibility is not None
    values = field.linear_compressibility
    return np.column_stack((np.maximum(values, 0.0), np.maximum(-values, 0.0)))


def sample_shear_modulus(
    tensor: ElasticTensor,
    theta: Sequence[float],
    phi: Sequence[float],
    progress_callback: ProgressCallback | None = None,
) -> np.ndarray:
    """Sample exact transverse shear-modulus extrema in one batch."""
    theta_array, phi_array = _paired_angles(theta, phi)
    field = sample_elastic_directional_field(
        tensor,
        theta_array,
        phi_array,
        properties=("shear",),
    )
    _replay_progress(theta_array.size, progress_callback)
    assert field.shear_minimum is not None
    assert field.shear_maximum is not None
    return np.column_stack((field.shear_minimum, field.shear_maximum))


def sample_poisson_ratio(
    tensor: ElasticTensor,
    theta: Sequence[float],
    phi: Sequence[float],
    progress_callback: ProgressCallback | None = None,
) -> np.ndarray:
    """Sample exact transverse Poisson-ratio extrema in one batch."""
    theta_array, phi_array = _paired_angles(theta, phi)
    field = sample_elastic_directional_field(
        tensor,
        theta_array,
        phi_array,
        properties=("poisson",),
    )
    _replay_progress(theta_array.size, progress_callback)
    assert field.poisson_minimum is not None
    assert field.poisson_maximum is not None
    minimum = field.poisson_minimum
    return np.column_stack(
        (
            np.minimum(minimum, 0.0),
            np.maximum(minimum, 0.0),
            field.poisson_maximum,
        )
    )


def _replay_progress(
    total: int,
    callback: ProgressCallback | None,
) -> None:
    """Preserve the historical per-direction wrapper callback contract."""
    if callback is None:
        return
    for current in range(1, total + 1):
        callback(current, total)


def _paired_angles(
    theta: Sequence[float],
    phi: Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Validate and return paired angular arrays."""
    theta_array = np.asarray(theta, dtype=float)
    phi_array = np.asarray(phi, dtype=float)
    if theta_array.shape != phi_array.shape:
        raise ValueError("theta and phi arrays must have the same shape.")
    return theta_array, phi_array


def _report_progress(
    current: int,
    total: int,
    progress_callback: ProgressCallback | None,
) -> None:
    """Report numerical progress through an optional callback."""
    if progress_callback is not None:
        progress_callback(current, total)


def _orthorhombic_young_modulus(
    tensor: OrthorhombicElasticTensor,
    angles: Sequence[float],
) -> float:
    """Evaluate the specialized orthorhombic Young-modulus expression."""
    ct2 = math.cos(angles[0]) ** 2
    st2 = 1 - ct2
    cf2 = math.cos(angles[1]) ** 2
    sf2 = 1 - cf2
    s = tensor.compliance_tensor
    s11 = s[0, 0, 0, 0]
    s22 = s[1, 1, 1, 1]
    s33 = s[2, 2, 2, 2]
    s44 = 4 * s[1, 2, 1, 2]
    s55 = 4 * s[0, 2, 0, 2]
    s66 = 4 * s[0, 1, 0, 1]
    s12 = s[0, 0, 1, 1]
    s13 = s[0, 0, 2, 2]
    s23 = s[1, 1, 2, 2]
    return float(
        1
        / (
            ct2**2 * s33
            + 2 * cf2 * ct2 * s13 * st2
            + cf2 * ct2 * s55 * st2
            + 2 * ct2 * s23 * sf2 * st2
            + ct2 * s44 * sf2 * st2
            + cf2**2 * s11 * st2**2
            + 2 * cf2 * s12 * sf2 * st2**2
            + cf2 * s66 * sf2 * st2**2
            + s22 * sf2**2 * st2**2
        )
    )


def _orthorhombic_linear_compressibility(
    tensor: OrthorhombicElasticTensor,
    angles: Sequence[float],
) -> float:
    """Evaluate the specialized orthorhombic compressibility expression."""
    ct2 = math.cos(angles[0]) ** 2
    cf2 = math.cos(angles[1]) ** 2
    s = tensor.compliance_tensor
    s11 = s[0, 0, 0, 0]
    s22 = s[1, 1, 1, 1]
    s33 = s[2, 2, 2, 2]
    s12 = s[0, 0, 1, 1]
    s13 = s[0, 0, 2, 2]
    s23 = s[1, 1, 2, 2]
    return float(
        1000
        * (
            ct2 * (s13 + s23 + s33)
            + (cf2 * (s11 + s12 + s13) + (s12 + s22 + s23) * (1 - cf2)) * (1 - ct2)
        )
    )


def _orthorhombic_shear_modulus(
    tensor: OrthorhombicElasticTensor,
    angles: Sequence[float],
) -> float:
    """Evaluate the specialized orthorhombic shear-modulus expression."""
    ct = math.cos(angles[0])
    ct2 = ct * ct
    st2 = 1 - ct2
    cf = math.cos(angles[1])
    sf = math.sin(angles[1])
    sf2 = sf * sf
    cx = math.cos(angles[2])
    cx2 = cx * cx
    sx = math.sin(angles[2])
    sx2 = 1 - cx2
    s = tensor.compliance_tensor
    s11 = s[0, 0, 0, 0]
    s22 = s[1, 1, 1, 1]
    s33 = s[2, 2, 2, 2]
    s44 = 4 * s[1, 2, 1, 2]
    s55 = 4 * s[0, 2, 0, 2]
    s66 = 4 * s[0, 1, 0, 1]
    s12 = s[0, 0, 1, 1]
    s13 = s[0, 0, 2, 2]
    s23 = s[1, 1, 2, 2]
    result = (
        ct2 * ct2 * cx2 * s44 * sf2
        + cx2 * s44 * sf2 * st2 * st2
        + 4 * cf**3 * ct * cx * (-2 * s11 + 2 * s12 + s66) * sf * st2 * sx
        + 2
        * cf
        * ct
        * cx
        * sf
        * (
            ct2 * (s44 - s55)
            + (
                4 * s13
                - 4 * s23
                - s44
                + s55
                - 4 * s12 * sf2
                + 4 * s22 * sf2
                - 2 * s66 * sf2
            )
            * st2
        )
        * sx
        + s66 * sf2 * sf2 * st2 * sx2
        + cf**4 * st2 * (4 * ct2 * cx2 * s11 + s66 * sx2)
        + ct2
        * (
            2 * cx2 * (2 * s33 + sf2 * (-4 * s23 - s44 + 2 * s22 * sf2)) * st2
            + s55 * sf2 * sx2
        )
        + cf**2
        * (
            ct2 * ct2 * cx2 * s55
            + ct2
            * (-2 * cx2 * (4 * s13 + s55 - 2 * (2 * s12 + s66) * sf2) * st2 + s44 * sx2)
            + st2
            * (cx2 * s55 * st2 + 2 * (2 * s11 - 4 * s12 + 2 * s22 - s66) * sf2 * sx2)
        )
    )
    return float(1 / result)


def _orthorhombic_poisson_ratio(
    tensor: OrthorhombicElasticTensor,
    angles: Sequence[float],
) -> float:
    """Evaluate the specialized orthorhombic Poisson-ratio expression."""
    ct = math.cos(angles[0])
    ct2 = ct * ct
    st2 = 1 - ct2
    cf = math.cos(angles[1])
    sf = math.sin(angles[1])
    cx = math.cos(angles[2])
    sx = math.sin(angles[2])
    s = tensor.compliance_tensor
    s11 = s[0, 0, 0, 0]
    s22 = s[1, 1, 1, 1]
    s33 = s[2, 2, 2, 2]
    s44 = 4 * s[1, 2, 1, 2]
    s55 = 4 * s[0, 2, 0, 2]
    s66 = 4 * s[0, 1, 0, 1]
    s12 = s[0, 0, 1, 1]
    s13 = s[0, 0, 2, 2]
    s23 = s[1, 1, 2, 2]
    numerator = (
        -(ct**2 * cx**2 * s33 * st2)
        - cf**2 * cx**2 * s13 * st2 * st2
        - cx**2 * s23 * sf**2 * st2 * st2
        + ct * cx * s44 * sf * st2 * (ct * cx * sf + cf * sx)
        - ct**2 * s23 * (ct * cx * sf + cf * sx) ** 2
        - cf**2 * s12 * st2 * (ct * cx * sf + cf * sx) ** 2
        - s22 * sf**2 * st2 * (ct * cx * sf + cf * sx) ** 2
        + cf * ct * cx * s55 * st2 * (cf * ct * cx - sf * sx)
        - cf * s66 * sf * st2 * (ct * cx * sf + cf * sx) * (cf * ct * cx - sf * sx)
        - ct**2 * s13 * (cf * ct * cx - sf * sx) ** 2
        - cf**2 * s11 * st2 * (cf * ct * cx - sf * sx) ** 2
        - s12 * sf**2 * st2 * (cf * ct * cx - sf * sx) ** 2
    )
    denominator = (
        ct**4 * s33
        + 2 * cf**2 * ct**2 * s13 * st2
        + cf**2 * ct**2 * s55 * st2
        + 2 * ct**2 * s23 * sf**2 * st2
        + ct**2 * s44 * sf**2 * st2
        + cf**4 * s11 * st2 * st2
        + 2 * cf**2 * s12 * sf**2 * st2 * st2
        + cf**2 * s66 * sf**2 * st2 * st2
        + s22 * sf**4 * st2 * st2
    )
    return float(numerator / denominator)
