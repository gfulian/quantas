# -*- coding: utf-8 -*-

"""Structural-path interpolation and anisotropic thermal expansion.

The routines in this module reconstruct a continuous crystallographic cell
along a one-dimensional volume path.  They are frontend-neutral and do not
know about QHA workflows, command-line interfaces, persistence, or plotting.

The structural path is represented by the deviatoric logarithmic stretch with
respect to one reference cell.  The isotropic logarithmic strain is imposed
analytically from the requested volume.  This construction preserves positive
cell volume and guarantees that the trace of the resulting thermal-expansion
tensor equals the supplied volumetric expansion coefficient, apart from
floating-point roundoff.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import polar

from quantas.core.geometry.cells import lattice_parameters
from quantas.core.math import matrix_exponential_frechet
from quantas.core.math.fitting import FitResult
from quantas.core.math.polynomials import FittedPolynomial, fit_polynomial
from quantas.models.structures import StructureVolumeSeries

FloatArray: TypeAlias = NDArray[np.float64]


@dataclass(slots=True)
class StructuralPathEvaluation:
    """Values reconstructed from a crystallographic structural path.

    Parameters
    ----------
    lattice : ndarray
        Crystallographic direct-lattice matrices with shape
        ``target_shape + (3, 3)`` and vectors stored by rows.
    lattice_parameters : ndarray
        Cell parameters ``a, b, c, alpha, beta, gamma`` with shape
        ``target_shape + (6,)``. Lengths are in angstrom and angles in degrees.
    lattice_parameter_uncertainties : ndarray or None
        Propagated one-standard-deviation uncertainties of the six cell
        parameters with shape ``target_shape + (6,)``.  Length uncertainties
        are in angstrom and angular uncertainties in degrees.
    lattice_parameter_derivatives : ndarray or None
        Temperature derivatives of the six cell parameters with shape
        ``target_shape + (6,)``. Length derivatives are in angstrom per kelvin
        and angular derivatives in degrees per kelvin.
    axial_thermal_expansion : ndarray or None
        Linear expansion coefficients of the crystallographic ``a``, ``b``,
        and ``c`` edges with shape ``target_shape + (3,)`` and units K\\ :sup:`-1`.
    thermal_expansion_tensor : ndarray or None
        Symmetric Cartesian thermal-expansion tensor with shape
        ``target_shape + (3, 3)`` and units K\\ :sup:`-1`.
    trace_residual : ndarray or None
        Difference between the tensor trace and the supplied volumetric
        expansion coefficient at every target point.
    extrapolation_mask : ndarray
        Boolean mask identifying target volumes outside the sampled structural
        interval.
    metadata : dict
        Method, basis, fit, and numerical diagnostics.
    """

    lattice: FloatArray
    lattice_parameters: FloatArray
    lattice_parameter_uncertainties: FloatArray | None
    lattice_parameter_derivatives: FloatArray | None
    axial_thermal_expansion: FloatArray | None
    thermal_expansion_tensor: FloatArray | None
    trace_residual: FloatArray | None
    extrapolation_mask: NDArray[np.bool_]
    metadata: dict[str, Any]


class StructuralPathModel:
    r"""Interpolate crystal shape as a function of volume.

    The model uses a pure-stretch representation. For every sampled lattice
    matrix :math:`A(V)`, the deformation gradient relative to the reference
    lattice is decomposed as :math:`F=RU`. The symmetric positive-definite
    stretch :math:`U` is represented through its logarithm,
    :math:`E=\log U`. Its trace is replaced by the exact volumetric term
    :math:`\log(V/V_r)`, while five independent deviatoric components are
    fitted as functions of :math:`x=\log(V/V_r)`.

    Parameters
    ----------
    series : StructureVolumeSeries
        Primitive structural volume series with a continuous orientation.
    degree : int, optional
        Polynomial degree used for the five deviatoric logarithmic-strain
        components. The effective degree is limited to ``nvol - 1``.

    Raises
    ------
    ValueError
        If the structural series is incomplete, singular, or cannot be fitted.
    """

    _COMPONENTS = ((0, 0), (1, 1), (0, 1), (0, 2), (1, 2))

    def __init__(self, series: StructureVolumeSeries, degree: int = 3) -> None:
        """Build the structural-path interpolation model."""
        if series.nvol < 2:
            raise ValueError("at least two structures are required")
        if degree < 1:
            raise ValueError("structural-path polynomial degree must be positive")
        if not np.all(np.isfinite(series.volumes)) or np.any(series.volumes <= 0.0):
            raise ValueError("structural volumes must be finite and positive")

        self.series = series
        self.degree = min(int(degree), series.nvol - 1)
        self.reference_index = int(series.reference_index)
        self.reference_volume = float(series.volumes[self.reference_index])
        self.basis_matrix, self.basis_source = _crystallographic_basis_matrix(series)
        self.sampled_lattices = np.einsum(
            "ij,vjk->vik",
            self.basis_matrix,
            np.asarray(series.lattices, dtype=np.float64),
        )
        self.reference_lattice = self.sampled_lattices[self.reference_index].copy()
        if abs(float(np.linalg.det(self.reference_lattice))) <= np.finfo(float).eps:
            raise ValueError("reference lattice is singular")

        self._x = np.log(
            np.asarray(series.volumes, dtype=np.float64) / self.reference_volume
        )
        self.is_cubic = _is_cubic_series(series, self.sampled_lattices)
        self._models: list[FittedPolynomial] = []
        self.fit_results: list[FitResult] = []
        if self.is_cubic:
            self.maximum_removed_rotation_degrees = 0.0
        else:
            deviatoric, rotation_angles = self._sampled_deviatoric_log_strain()
            self.maximum_removed_rotation_degrees = float(np.max(rotation_angles))
            for column in range(deviatoric.shape[1]):
                fit, model = fit_polynomial(
                    self._x,
                    deviatoric[:, column],
                    self.degree,
                    scale_coordinate=True,
                )
                self.fit_results.append(fit)
                if not fit.success or model is None:
                    raise ValueError(
                        "structural-path fit failed for deviatoric component "
                        f"{column}: {fit.message}"
                    )
                self._models.append(model)

    def evaluate(
        self,
        volume: FloatArray | float,
        volumetric_expansion: FloatArray | float | None = None,
        volume_uncertainty: FloatArray | float | None = None,
        *,
        include_fit_uncertainty: bool = True,
    ) -> StructuralPathEvaluation:
        """Evaluate lattice parameters and thermal expansion.

        Parameters
        ----------
        volume : array-like or float
            Target primitive-cell volumes in the same volume unit as the
            sampled structural series.
        volumetric_expansion : array-like, float, or None, optional
            Volumetric thermal-expansion coefficient at every target point.
            When omitted, only cell matrices and cell parameters are returned.
        volume_uncertainty : array-like, float, or None, optional
            One-standard-deviation uncertainty of each target volume.  Values
            must be non-negative and broadcastable to ``volume``.
        include_fit_uncertainty : bool, optional
            Propagate the covariance of the independent structural-path fits.
            This source is absent for the exact cubic branch.

        Returns
        -------
        StructuralPathEvaluation
            Reconstructed crystallographic cells and, when requested, axial
            and tensorial thermal-expansion data.

        Raises
        ------
        ValueError
            If target volumes or expansion coefficients are invalid or have
            incompatible shapes.
        """
        target = np.asarray(volume, dtype=np.float64)
        finite_target = np.isfinite(target)
        if np.any(target[finite_target] <= 0.0):
            raise ValueError("finite target volumes must be positive")
        target_shape = target.shape
        flat_volume = target.reshape(-1)
        flat_finite = finite_target.reshape(-1)
        x_values = np.full(flat_volume.shape, np.nan, dtype=np.float64)
        x_values[flat_finite] = np.log(flat_volume[flat_finite] / self.reference_volume)

        alpha_flat: FloatArray | None
        if volumetric_expansion is None:
            alpha_flat = None
        else:
            alpha = np.asarray(volumetric_expansion, dtype=np.float64)
            try:
                alpha = np.broadcast_to(alpha, target_shape)
            except ValueError as exc:
                raise ValueError(
                    "volumetric expansion must be broadcastable to target volumes"
                ) from exc
            alpha_flat = np.asarray(alpha, dtype=np.float64).reshape(-1)

        sigma_volume_flat: FloatArray | None
        if volume_uncertainty is None:
            sigma_volume_flat = None
        else:
            sigma_volume = np.asarray(volume_uncertainty, dtype=np.float64)
            try:
                sigma_volume = np.broadcast_to(sigma_volume, target_shape)
            except ValueError as exc:
                raise ValueError(
                    "volume uncertainty must be broadcastable to target volumes"
                ) from exc
            finite_sigma = np.isfinite(sigma_volume)
            if np.any(sigma_volume[finite_sigma] < 0.0):
                raise ValueError("finite volume uncertainties must be non-negative")
            sigma_volume_flat = np.asarray(sigma_volume, dtype=np.float64).reshape(-1)

        lattice = np.full((flat_volume.size, 3, 3), np.nan, dtype=np.float64)
        parameters = np.full((flat_volume.size, 6), np.nan, dtype=np.float64)
        calculate_uncertainties = sigma_volume_flat is not None or (
            include_fit_uncertainty and not self.is_cubic
        )
        parameter_uncertainties = (
            np.full((flat_volume.size, 6), np.nan, dtype=np.float64)
            if calculate_uncertainties
            else None
        )
        parameter_derivatives = (
            None
            if alpha_flat is None
            else np.full((flat_volume.size, 6), np.nan, dtype=np.float64)
        )
        axial = (
            None
            if alpha_flat is None
            else np.full((flat_volume.size, 3), np.nan, dtype=np.float64)
        )
        tensor = (
            None
            if alpha_flat is None
            else np.full((flat_volume.size, 3, 3), np.nan, dtype=np.float64)
        )
        residual = (
            None
            if alpha_flat is None
            else np.full(flat_volume.size, np.nan, dtype=np.float64)
        )

        for index, x_value in enumerate(x_values):
            if not flat_finite[index]:
                continue
            log_stretch, derivative = self._log_stretch(float(x_value))
            if self.is_cubic:
                scale = float(np.exp(float(x_value) / 3.0))
                stretch = np.eye(3, dtype=np.float64) * scale
                dstretch_dx = np.eye(3, dtype=np.float64) * (scale / 3.0)
            else:
                stretch = _symmetric_exponential(log_stretch)
                dstretch_dx = matrix_exponential_frechet(log_stretch, derivative)
            current_lattice = self.reference_lattice @ stretch
            lattice[index] = current_lattice
            parameters[index] = lattice_parameters(current_lattice)

            if parameter_uncertainties is not None:
                variance: FloatArray = np.zeros(6, dtype=np.float64)
                if sigma_volume_flat is not None:
                    sigma_volume_value = float(sigma_volume_flat[index])
                    if np.isfinite(sigma_volume_value):
                        dlattice_dx = self.reference_lattice @ dstretch_dx
                        dparameters_dx = lattice_parameter_derivatives(
                            current_lattice,
                            dlattice_dx,
                        )
                        sigma_x = sigma_volume_value / flat_volume[index]
                        variance += np.square(dparameters_dx * sigma_x)
                if include_fit_uncertainty and not self.is_cubic:
                    fit_variance = self._fit_parameter_variance(
                        float(x_value),
                        log_stretch,
                        current_lattice,
                    )
                    variance += fit_variance
                parameter_uncertainties[index] = np.sqrt(np.maximum(variance, 0.0))

            if alpha_flat is None:
                continue
            if (
                parameter_derivatives is None
                or axial is None
                or tensor is None
                or residual is None
            ):
                raise RuntimeError(
                    "thermal-expansion output buffers were not initialized"
                )
            alpha_v = float(alpha_flat[index])
            if not np.isfinite(alpha_v):
                parameter_derivatives[index].fill(np.nan)
                axial[index].fill(np.nan)
                tensor[index].fill(np.nan)
                residual[index] = np.nan
                continue

            if self.is_cubic:
                alpha_linear = alpha_v / 3.0
                dlattice_dt = current_lattice * alpha_linear
                alpha_tensor = np.eye(3, dtype=np.float64) * alpha_linear
                axial_values: FloatArray = np.full(3, alpha_linear, dtype=np.float64)
            else:
                dstretch_dt = dstretch_dx * alpha_v
                dlattice_dt = self.reference_lattice @ dstretch_dt
                rate = dstretch_dt @ np.linalg.inv(stretch)
                alpha_tensor = 0.5 * (rate + rate.T)
                axial_values = axial_expansion(current_lattice, dlattice_dt)

            tensor[index] = alpha_tensor
            axial[index] = axial_values
            parameter_derivatives[index] = lattice_parameter_derivatives(
                current_lattice,
                dlattice_dt,
            )
            residual[index] = float(np.trace(alpha_tensor) - alpha_v)

        minimum = float(np.min(self.series.volumes))
        maximum = float(np.max(self.series.volumes))
        extrapolation = finite_target & ((target < minimum) | (target > maximum))
        metadata = {
            "method": "volume_constrained_structural_path",
            "coordinate": "logarithmic_volume",
            "strain_measure": "deviatoric_logarithmic_stretch",
            "polynomial_degree": int(self.degree),
            "reference_index": int(self.reference_index),
            "reference_primitive_volume": float(self.reference_volume),
            "crystallographic_basis_source": self.basis_source,
            "crystallographic_volume_factor": abs(
                float(np.linalg.det(self.basis_matrix))
            ),
            "sampled_primitive_volume_min": minimum,
            "sampled_primitive_volume_max": maximum,
            "maximum_removed_rotation_degrees": self.maximum_removed_rotation_degrees,
            "fit_diagnostics": [fit.as_dict() for fit in self.fit_results],
            "calculation_branch": (
                "cubic_exact" if self.is_cubic else "general_anisotropic_path"
            ),
            "symmetry_shortcut": "alpha_linear=alphaV/3" if self.is_cubic else None,
            "uncertainty_method": (
                "first_order_delta_method"
                if parameter_uncertainties is not None
                else "none"
            ),
            "uncertainty_sources": {
                "equilibrium_volume": bool(
                    sigma_volume_flat is not None
                    and np.any(np.isfinite(sigma_volume_flat))
                ),
                "equilibrium_volume_points": int(
                    0
                    if sigma_volume_flat is None
                    else np.count_nonzero(np.isfinite(sigma_volume_flat))
                ),
                "structural_fit_covariance": bool(
                    include_fit_uncertainty and not self.is_cubic
                ),
                "source_structure_uncertainty": False,
                "cross_covariance_between_deviatoric_fits": False,
            },
        }
        if residual is not None:
            finite = residual[np.isfinite(residual)]
            metadata["maximum_trace_residual"] = (
                float(np.max(np.abs(finite))) if finite.size else np.nan
            )

        return StructuralPathEvaluation(
            lattice=lattice.reshape(target_shape + (3, 3)),
            lattice_parameters=parameters.reshape(target_shape + (6,)),
            lattice_parameter_uncertainties=(
                None
                if parameter_uncertainties is None
                else parameter_uncertainties.reshape(target_shape + (6,))
            ),
            lattice_parameter_derivatives=(
                None
                if parameter_derivatives is None
                else parameter_derivatives.reshape(target_shape + (6,))
            ),
            axial_thermal_expansion=(
                None if axial is None else axial.reshape(target_shape + (3,))
            ),
            thermal_expansion_tensor=(
                None if tensor is None else tensor.reshape(target_shape + (3, 3))
            ),
            trace_residual=(
                None if residual is None else residual.reshape(target_shape)
            ),
            extrapolation_mask=np.asarray(extrapolation, dtype=bool),
            metadata=metadata,
        )

    def _sampled_deviatoric_log_strain(self) -> tuple[FloatArray, FloatArray]:
        """Return independent deviatoric log strains and removed rotations."""
        inverse_reference_transpose = np.linalg.inv(self.reference_lattice.T)
        values = np.empty((self.series.nvol, len(self._COMPONENTS)), dtype=np.float64)
        rotations = np.empty(self.series.nvol, dtype=np.float64)
        for index, sampled in enumerate(self.sampled_lattices):
            deformation = sampled.T @ inverse_reference_transpose
            rotation, stretch = polar(deformation, side="right")
            log_stretch = _symmetric_logarithm(np.asarray(stretch, dtype=np.float64))
            x_value = float(self._x[index])
            deviatoric = log_stretch - np.eye(3, dtype=np.float64) * (x_value / 3.0)
            deviatoric -= np.eye(3, dtype=np.float64) * (np.trace(deviatoric) / 3.0)
            values[index] = _independent_deviatoric_components(deviatoric)
            rotations[index] = _rotation_angle_degrees(
                np.asarray(rotation, dtype=np.float64)
            )
        return values, rotations

    def _log_stretch(self, x_value: float) -> tuple[FloatArray, FloatArray]:
        """Return logarithmic stretch and its derivative with respect to log V."""
        identity = np.eye(3, dtype=np.float64)
        if self.is_cubic:
            return identity * (x_value / 3.0), identity / 3.0
        values = np.asarray(
            [float(model.evaluate(x_value)) for model in self._models],
            dtype=np.float64,
        )
        derivatives = np.asarray(
            [float(model.derivative(x_value)) for model in self._models],
            dtype=np.float64,
        )
        deviatoric = _deviatoric_matrix(values)
        ddeviatoric = _deviatoric_matrix(derivatives)
        return (
            identity * (x_value / 3.0) + deviatoric,
            identity / 3.0 + ddeviatoric,
        )

    def _fit_parameter_variance(
        self,
        x_value: float,
        log_stretch: FloatArray,
        current_lattice: FloatArray,
    ) -> FloatArray:
        """Return parameter variances propagated from structural fit covariance."""
        variance: FloatArray = np.zeros(6, dtype=np.float64)
        for component, (fit, model) in enumerate(zip(self.fit_results, self._models)):
            covariance = fit.covariance
            if covariance is None:
                continue
            covariance_array = np.asarray(covariance, dtype=np.float64)
            nparameters = model.parameters.size
            if covariance_array.shape != (nparameters, nparameters):
                continue
            scaled = float(model.scaled_coordinate(x_value).item())
            jacobian = np.empty((6, nparameters), dtype=np.float64)
            for order in range(nparameters):
                component_direction: FloatArray = np.zeros(5, dtype=np.float64)
                component_direction[component] = scaled**order
                direction = _deviatoric_matrix(component_direction)
                dstretch = matrix_exponential_frechet(log_stretch, direction)
                dlattice = self.reference_lattice @ dstretch
                jacobian[:, order] = lattice_parameter_derivatives(
                    current_lattice,
                    dlattice,
                )
            propagated = jacobian @ covariance_array @ jacobian.T
            variance += np.maximum(np.diag(propagated), 0.0)
        return variance


def axial_expansion(lattice: FloatArray, derivative: FloatArray) -> FloatArray:
    r"""Return linear expansion coefficients of the three lattice edges.

    Parameters
    ----------
    lattice : array-like
        Direct lattice matrix with vectors stored by rows.
    derivative : array-like
        Temperature derivative of the lattice matrix.

    Returns
    -------
    ndarray
        :math:`(1/a_i)(\partial a_i/\partial T)` for ``a``, ``b``, and ``c``.
    """
    matrix = np.asarray(lattice, dtype=np.float64)
    rate = np.asarray(derivative, dtype=np.float64)
    denominator = np.einsum("ij,ij->i", matrix, matrix)
    if np.any(denominator <= np.finfo(float).eps):
        raise ValueError("lattice vectors must have positive length")
    return np.einsum("ij,ij->i", matrix, rate) / denominator


def lattice_parameter_derivatives(
    lattice: FloatArray,
    derivative: FloatArray,
) -> FloatArray:
    """Return temperature derivatives of lengths and interaxial angles.

    Parameters
    ----------
    lattice : array-like
        Direct lattice matrix with vectors stored by rows.
    derivative : array-like
        Temperature derivative of the lattice matrix.

    Returns
    -------
    ndarray
        Derivatives ``da/dT, db/dT, dc/dT, dalpha/dT, dbeta/dT, dgamma/dT``.
        Angular derivatives are expressed in degrees per kelvin.
    """
    matrix = np.asarray(lattice, dtype=np.float64)
    rate = np.asarray(derivative, dtype=np.float64)
    lengths = np.linalg.norm(matrix, axis=1)
    length_rate = axial_expansion(matrix, rate) * lengths
    angle_rate = np.asarray(
        [
            _angle_derivative_degrees(matrix[1], matrix[2], rate[1], rate[2]),
            _angle_derivative_degrees(matrix[0], matrix[2], rate[0], rate[2]),
            _angle_derivative_degrees(matrix[0], matrix[1], rate[0], rate[1]),
        ],
        dtype=np.float64,
    )
    return np.concatenate((length_rate, angle_rate))


def _crystallographic_basis_matrix(
    series: StructureVolumeSeries,
) -> tuple[FloatArray, str]:
    """Return a constant primitive-to-crystallographic row-basis transform."""
    if series.primitive_to_crystallographic is not None:
        return (
            np.asarray(series.primitive_to_crystallographic, dtype=np.float64).T,
            "crystal_output",
        )
    if (
        series.symmetry is not None
        and series.symmetry.transformation_matrix is not None
    ):
        transformation = np.asarray(
            series.symmetry.transformation_matrix,
            dtype=np.float64,
        )
        if abs(float(np.linalg.det(transformation))) > np.finfo(float).eps:
            basis = np.asarray(
                np.linalg.inv(transformation).T,
                dtype=np.float64,
            )
            return basis, "spglib_standard"
    return np.eye(3, dtype=np.float64), "primitive"


def _is_cubic_series(
    series: StructureVolumeSeries,
    lattices: FloatArray,
) -> bool:
    """Return whether symmetry and metrics define a cubic structural path."""
    if series.symmetry is not None:
        number = int(series.symmetry.space_group_number)
        if 195 <= number <= 230:
            return True
    for lattice in np.asarray(lattices, dtype=np.float64):
        parameters = lattice_parameters(lattice)
        lengths = parameters[:3]
        angles = parameters[3:]
        if not np.allclose(lengths, np.mean(lengths), rtol=1.0e-7, atol=1.0e-8):
            return False
        if not np.allclose(angles, 90.0, rtol=0.0, atol=1.0e-6):
            return False
    return True


def _independent_deviatoric_components(matrix: FloatArray) -> FloatArray:
    """Return five independent components of a symmetric traceless matrix."""
    return np.asarray(
        [matrix[0, 0], matrix[1, 1], matrix[0, 1], matrix[0, 2], matrix[1, 2]],
        dtype=np.float64,
    )


def _deviatoric_matrix(values: FloatArray) -> FloatArray:
    """Build a symmetric traceless matrix from five independent components."""
    v = np.asarray(values, dtype=np.float64)
    if v.shape != (5,):
        raise ValueError("five deviatoric components are required")
    return np.asarray(
        [
            [v[0], v[2], v[3]],
            [v[2], v[1], v[4]],
            [v[3], v[4], -v[0] - v[1]],
        ],
        dtype=np.float64,
    )


def _symmetric_logarithm(matrix: FloatArray) -> FloatArray:
    """Return the real logarithm of a symmetric positive-definite matrix."""
    eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (matrix + matrix.T))
    if np.any(eigenvalues <= 0.0):
        raise ValueError("stretch tensor must be positive definite")
    return (eigenvectors * np.log(eigenvalues)) @ eigenvectors.T


def _symmetric_exponential(matrix: FloatArray) -> FloatArray:
    """Return the exponential of a symmetric matrix."""
    eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (matrix + matrix.T))
    return (eigenvectors * np.exp(eigenvalues)) @ eigenvectors.T


def _rotation_angle_degrees(rotation: FloatArray) -> float:
    """Return the proper-rotation angle in degrees for diagnostics."""
    cosine = np.clip((float(np.trace(rotation)) - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(cosine)))


def _angle_degrees(first: FloatArray, second: FloatArray) -> float:
    """Return the angle between two vectors in degrees."""
    denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
    cosine = np.clip(float(np.dot(first, second)) / denominator, -1.0, 1.0)
    return float(np.degrees(np.arccos(cosine)))


def _angle_derivative_degrees(
    first: FloatArray,
    second: FloatArray,
    dfirst: FloatArray,
    dsecond: FloatArray,
) -> float:
    """Return the temperature derivative of the angle between two vectors."""
    norm_first = float(np.linalg.norm(first))
    norm_second = float(np.linalg.norm(second))
    cosine = np.clip(
        float(np.dot(first, second)) / (norm_first * norm_second),
        -1.0,
        1.0,
    )
    sine = float(np.sqrt(max(1.0 - cosine * cosine, 0.0)))
    if sine <= 1.0e-14:
        return np.nan
    dcosine = float(np.dot(dfirst, second) + np.dot(first, dsecond)) / (
        norm_first * norm_second
    ) - cosine * (
        float(np.dot(first, dfirst)) / (norm_first * norm_first)
        + float(np.dot(second, dsecond)) / (norm_second * norm_second)
    )
    return float(-dcosine / sine * 180.0 / np.pi)
