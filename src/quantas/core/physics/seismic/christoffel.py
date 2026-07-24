# -*- coding: utf-8 -*-

"""Christoffel equation for acoustic propagation in elastic crystals."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.physics.elasticity.tensor import ElasticTensor

from .medium import ElasticMedium

from .enhancement import (
    DirectionalEnhancementResult,
    calculate_area_factor,
    calculate_eigenvalue_hessian,
    calculate_ray_direction_gradient,
)
from .group import DirectionalGroupResult
from .modes import DirectionalPhaseResult


_GPA_DENSITY_TO_KM2_S2 = 1000.0


class ChristoffelSolver:
    """Solve the acoustic Christoffel equation for a three-dimensional solid.

    Parameters
    ----------
    medium : ElasticMedium
        Elastic medium containing the stiffness tensor in GPa and the material
        density in kg m^-3.
    eigenvalue_rtol, eigenvalue_atol : float, optional
        Relative and absolute tolerances used to distinguish small negative
        numerical eigenvalues from physically invalid eigenvalues. The
        absolute tolerance is expressed in km^2 s^-2.
    degeneracy_rtol, degeneracy_atol : float, optional
        Relative and absolute tolerances used to identify degenerate
        Christoffel eigenvalues. The absolute tolerance is expressed in
        km^2 s^-2.
    pseudoinverse_rcond : float, optional
        Relative cutoff used by the pseudoinverse in the analytical
        eigenvalue-Hessian calculation.
    caustic_rtol, caustic_atol : float, optional
        Relative and absolute tolerances used to identify area factors close
        to zero. These thresholds are dimensionless.

    Raises
    ------
    TypeError
        If ``medium`` is not an :class:`ElasticMedium` instance.
    ValueError
        If any numerical tolerance is invalid.
    """

    def __init__(
        self,
        medium: ElasticMedium,
        *,
        eigenvalue_rtol: float = 1.0e-10,
        eigenvalue_atol: float = 1.0e-12,
        degeneracy_rtol: float = 1.0e-8,
        degeneracy_atol: float = 1.0e-10,
        pseudoinverse_rcond: float = 1.0e-10,
        caustic_rtol: float = 1.0e-10,
        caustic_atol: float = 1.0e-12,
    ) -> None:
        if not isinstance(medium, ElasticMedium):
            raise TypeError("ChristoffelSolver requires an ElasticMedium instance.")

        self._medium = medium
        self._elastic_tensor = medium.elastic_tensor
        self._density = medium.density
        self._eigenvalue_rtol = _validate_tolerance(eigenvalue_rtol, "eigenvalue_rtol")
        self._eigenvalue_atol = _validate_tolerance(eigenvalue_atol, "eigenvalue_atol")
        self._degeneracy_rtol = _validate_tolerance(degeneracy_rtol, "degeneracy_rtol")
        self._degeneracy_atol = _validate_tolerance(degeneracy_atol, "degeneracy_atol")
        self._pseudoinverse_rcond = _validate_rcond(pseudoinverse_rcond)
        self._caustic_rtol = _validate_tolerance(caustic_rtol, "caustic_rtol")
        self._caustic_atol = _validate_tolerance(caustic_atol, "caustic_atol")

        self._reduced_stiffness = np.asarray(
            self._elastic_tensor.stiffness_tensor
            * (_GPA_DENSITY_TO_KM2_S2 / self._density),
            dtype=float,
        )
        self._christoffel_hessian = _build_christoffel_hessian(self._reduced_stiffness)
        self._gradient_tensor = self._reduced_stiffness + np.transpose(
            self._reduced_stiffness,
            (0, 2, 1, 3),
        )

        for array in (
            self._reduced_stiffness,
            self._christoffel_hessian,
            self._gradient_tensor,
        ):
            array.setflags(write=False)

    @property
    def medium(self) -> ElasticMedium:
        """Return the elastic medium used by the solver."""
        return self._medium

    @property
    def elastic_tensor(self) -> ElasticTensor:
        """Return the elastic tensor used by the solver."""
        return self._elastic_tensor

    @property
    def density(self) -> float:
        """Return the material density in kg m^-3."""
        return self._density

    @property
    def reduced_stiffness(self) -> NDArray[np.float64]:
        """Return ``1000 C_ijkl / density`` in km^2 s^-2."""
        return self._reduced_stiffness

    @property
    def christoffel_hessian(self) -> NDArray[np.float64]:
        """Return the direction-independent Hessian of the Christoffel matrix."""
        return self._christoffel_hessian

    @property
    def eigenvalue_rtol(self) -> float:
        """Return the relative negative-eigenvalue tolerance."""
        return self._eigenvalue_rtol

    @property
    def eigenvalue_atol(self) -> float:
        """Return the absolute negative-eigenvalue tolerance in km^2 s^-2."""
        return self._eigenvalue_atol

    @property
    def degeneracy_rtol(self) -> float:
        """Return the relative eigenvalue-degeneracy tolerance."""
        return self._degeneracy_rtol

    @property
    def degeneracy_atol(self) -> float:
        """Return the absolute eigenvalue-degeneracy tolerance in km^2 s^-2."""
        return self._degeneracy_atol

    @property
    def pseudoinverse_rcond(self) -> float:
        """Return the relative cutoff used by eigenvalue pseudoinverses."""
        return self._pseudoinverse_rcond

    @property
    def caustic_rtol(self) -> float:
        """Return the relative area-factor tolerance for possible caustics."""
        return self._caustic_rtol

    @property
    def caustic_atol(self) -> float:
        """Return the absolute area-factor tolerance for possible caustics."""
        return self._caustic_atol

    def normalize_direction(self, direction: ArrayLike) -> NDArray[np.float64]:
        """Validate and normalize a Cartesian wave-normal direction.

        Parameters
        ----------
        direction : array_like
            Cartesian direction with shape ``(3,)``.

        Returns
        -------
        ndarray
            Read-only unit direction with shape ``(3,)``.

        Raises
        ------
        ValueError
            If the input has the wrong shape, contains non-finite values, or
            has zero norm.
        """
        return _normalize_direction(direction)

    def christoffel_matrix(self, direction: ArrayLike) -> NDArray[np.float64]:
        """Calculate the symmetric Christoffel matrix for one direction.

        Parameters
        ----------
        direction : array_like
            Cartesian wave-normal direction.

        Returns
        -------
        ndarray
            Symmetric matrix in km^2 s^-2 with shape ``(3, 3)``.

        Raises
        ------
        ValueError
            If ``direction`` is invalid.
        """
        q = _normalize_direction(direction)
        matrix = self._christoffel_matrix_from_unit_direction(q)
        matrix.setflags(write=False)
        return matrix

    def _christoffel_matrix_from_unit_direction(
        self,
        direction: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Calculate the Christoffel matrix from a validated unit direction.

        Parameters
        ----------
        direction : ndarray
            Unit wave-normal direction with shape ``(3,)``.

        Returns
        -------
        ndarray
            Symmetric Christoffel matrix in km^2 s^-2.
        """
        matrix = np.einsum(
            "j,ijkl,l->ik",
            direction,
            self._reduced_stiffness,
            direction,
            optimize=True,
        )
        return np.asarray(0.5 * (matrix + matrix.T), dtype=float)

    def christoffel_gradient(
        self,
        direction: ArrayLike,
    ) -> NDArray[np.float64]:
        """Calculate the analytical gradient of the Christoffel matrix.

        Parameters
        ----------
        direction : array_like
            Cartesian wave-normal direction.

        Returns
        -------
        ndarray
            Gradient with shape ``(3, 3, 3)``. The first axis identifies the
            derivative coordinate and the final two axes contain a symmetric
            matrix in km^2 s^-2.

        Raises
        ------
        ValueError
            If ``direction`` is invalid.
        """
        q = _normalize_direction(direction)
        gradient = self._christoffel_gradient_from_unit_direction(q)
        gradient.setflags(write=False)
        return gradient

    def _christoffel_gradient_from_unit_direction(
        self,
        direction: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Calculate the Christoffel gradient from a unit direction.

        Parameters
        ----------
        direction : ndarray
            Unit wave-normal direction with shape ``(3,)``.

        Returns
        -------
        ndarray
            Gradient with shape ``(3, 3, 3)`` in km^2 s^-2.
        """
        gradient = np.einsum(
            "k,iakl->ail",
            direction,
            self._gradient_tensor,
            optimize=True,
        )
        return np.asarray(
            0.5 * (gradient + np.swapaxes(gradient, 1, 2)),
            dtype=float,
        )

    def solve_direction(self, direction: ArrayLike) -> DirectionalPhaseResult:
        """Solve the phase-velocity problem for one wave-normal direction.

        Parameters
        ----------
        direction : array_like
            Cartesian wave-normal direction. The vector is normalized
            internally.

        Returns
        -------
        DirectionalPhaseResult
            Ordered eigenvalues, phase speeds, polarization axes, validity
            masks and degeneracy diagnostics.

        Raises
        ------
        ValueError
            If ``direction`` is invalid.
        """
        q = _normalize_direction(direction)
        christoffel = self._christoffel_matrix_from_unit_direction(q)
        eigenvalues, eigenvectors = np.linalg.eigh(christoffel)
        eigenvalues = np.asarray(eigenvalues, dtype=float)
        polarizations = np.asarray(eigenvectors.T, dtype=float)

        scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).tiny)
        eigenvalue_threshold = self._eigenvalue_atol + self._eigenvalue_rtol * scale
        valid_mask = eigenvalues >= -eigenvalue_threshold
        clamped_mask = valid_mask & (eigenvalues < 0.0)
        values_for_speed = np.where(clamped_mask, 0.0, eigenvalues)
        phase_speeds = np.full(3, np.nan, dtype=float)
        phase_speeds[valid_mask] = np.sqrt(values_for_speed[valid_mask])

        eigenvalue_gaps = np.diff(eigenvalues)
        relative_eigenvalue_gaps = eigenvalue_gaps / scale
        mode_eigenvalue_gaps = np.asarray(
            [
                eigenvalue_gaps[0],
                min(eigenvalue_gaps[0], eigenvalue_gaps[1]),
                eigenvalue_gaps[1],
            ],
            dtype=float,
        )
        mode_relative_eigenvalue_gaps = mode_eigenvalue_gaps / scale
        degeneracy_threshold = self._degeneracy_atol + self._degeneracy_rtol * scale
        degeneracy_mask = mode_eigenvalue_gaps <= degeneracy_threshold

        arrays = (
            q,
            christoffel,
            eigenvalues,
            polarizations,
            phase_speeds,
            eigenvalue_gaps,
            relative_eigenvalue_gaps,
            mode_eigenvalue_gaps,
            mode_relative_eigenvalue_gaps,
            valid_mask,
            clamped_mask,
            degeneracy_mask,
        )
        for array in arrays:
            array.setflags(write=False)

        return DirectionalPhaseResult(
            direction=q,
            christoffel=christoffel,
            eigenvalues=eigenvalues,
            polarizations=polarizations,
            phase_speeds=phase_speeds,
            eigenvalue_gaps=eigenvalue_gaps,
            relative_eigenvalue_gaps=relative_eigenvalue_gaps,
            mode_eigenvalue_gaps=mode_eigenvalue_gaps,
            mode_relative_eigenvalue_gaps=mode_relative_eigenvalue_gaps,
            valid_mask=valid_mask,
            clamped_mask=clamped_mask,
            degeneracy_mask=degeneracy_mask,
            eigenvalue_threshold=float(eigenvalue_threshold),
            degeneracy_threshold=float(degeneracy_threshold),
        )

    def solve_group_direction(
        self,
        direction: ArrayLike,
    ) -> DirectionalGroupResult:
        """Solve phase and group propagation for one wave-normal direction.

        Parameters
        ----------
        direction : array_like
            Cartesian wave-normal direction. The vector is normalized
            internally.

        Returns
        -------
        DirectionalGroupResult
            Phase solution, analytical eigenvalue gradients, group velocities,
            ray directions and power-flow angles.

        Raises
        ------
        ValueError
            If ``direction`` is invalid.
        """
        phase = self.solve_direction(direction)
        gradient = self._christoffel_gradient_from_unit_direction(phase.direction)
        eigenvalue_gradients = np.einsum(
            "mi,aij,mj->ma",
            phase.polarizations,
            gradient,
            phase.polarizations,
            optimize=True,
        )

        valid_mask = (
            phase.valid_mask
            & np.isfinite(phase.phase_speeds)
            & (phase.phase_speeds > 0.0)
        )
        group_velocities = np.full((3, 3), np.nan, dtype=float)
        group_velocities[valid_mask] = eigenvalue_gradients[valid_mask] / (
            2.0 * phase.phase_speeds[valid_mask, np.newaxis]
        )

        group_speeds = np.full(3, np.nan, dtype=float)
        group_speeds[valid_mask] = np.linalg.norm(
            group_velocities[valid_mask],
            axis=1,
        )
        valid_mask = valid_mask & np.isfinite(group_speeds) & (group_speeds > 0.0)

        ray_directions = np.full((3, 3), np.nan, dtype=float)
        ray_directions[valid_mask] = (
            group_velocities[valid_mask] / group_speeds[valid_mask, np.newaxis]
        )

        power_flow_angles = np.full(3, np.nan, dtype=float)
        cosines = np.einsum(
            "mi,i->m",
            ray_directions[valid_mask],
            phase.direction,
            optimize=True,
        )
        cosines = np.clip(cosines, -1.0, 1.0)
        endpoint_tolerance = 8.0 * np.finfo(float).eps
        cosines[np.abs(cosines - 1.0) <= endpoint_tolerance] = 1.0
        cosines[np.abs(cosines + 1.0) <= endpoint_tolerance] = -1.0
        power_flow_angles[valid_mask] = np.arccos(cosines)
        resolved_mask = valid_mask & ~phase.degeneracy_mask

        arrays = (
            eigenvalue_gradients,
            group_velocities,
            group_speeds,
            ray_directions,
            power_flow_angles,
            valid_mask,
            resolved_mask,
        )
        for array in arrays:
            array.setflags(write=False)

        return DirectionalGroupResult(
            phase=phase,
            eigenvalue_gradients=eigenvalue_gradients,
            group_velocities=group_velocities,
            group_speeds=group_speeds,
            ray_directions=ray_directions,
            power_flow_angles=power_flow_angles,
            valid_mask=valid_mask,
            resolved_mask=resolved_mask,
        )

    def solve_enhancement_direction(
        self,
        direction: ArrayLike,
    ) -> DirectionalEnhancementResult:
        """Solve analytical ray curvature and enhancement for one direction.

        Parameters
        ----------
        direction : array_like
            Cartesian wave-normal direction. The vector is normalized
            internally.

        Returns
        -------
        DirectionalEnhancementResult
            Phase and group quantities, eigenvalue Hessians, ray-direction
            gradients, area factors and enhancement diagnostics.

        Raises
        ------
        ValueError
            If ``direction`` is invalid.
        """
        group = self.solve_group_direction(direction)
        phase = group.phase
        gradient = self._christoffel_gradient_from_unit_direction(phase.direction)

        eigenvalue_hessians = np.full((3, 3, 3), np.nan, dtype=float)
        ray_direction_gradients = np.full((3, 3, 3), np.nan, dtype=float)
        area_factors = np.full(3, np.nan, dtype=float)
        caustic_thresholds = np.full(3, np.nan, dtype=float)
        enhancement = np.full(3, np.nan, dtype=float)
        log10_enhancement = np.full(3, np.nan, dtype=float)
        valid_mask = np.zeros(3, dtype=bool)
        resolved_mask = np.asarray(group.resolved_mask, dtype=bool).copy()
        finite_mask = np.zeros(3, dtype=bool)
        caustic_candidate_mask = np.zeros(3, dtype=bool)

        for mode in range(3):
            if not resolved_mask[mode]:
                continue

            hessian = calculate_eigenvalue_hessian(
                phase.christoffel,
                gradient,
                self._christoffel_hessian,
                float(phase.eigenvalues[mode]),
                phase.polarizations[mode],
                pseudoinverse_rcond=self._pseudoinverse_rcond,
            )
            phase_speed = float(phase.phase_speeds[mode])
            group_speed = float(group.group_speeds[mode])
            ray_direction = group.ray_directions[mode]
            denominator = 2.0 * phase_speed * group_speed
            if (
                not np.all(np.isfinite(hessian))
                or not np.isfinite(denominator)
                or denominator <= 0.0
                or not np.all(np.isfinite(ray_direction))
            ):
                continue

            ray_gradient = calculate_ray_direction_gradient(
                hessian,
                ray_direction,
                phase_speed,
                group_speed,
            )
            area_factor, cofactor_scale = calculate_area_factor(
                ray_gradient,
                phase.direction,
            )
            cofactor_scale = max(cofactor_scale, np.finfo(float).tiny)
            caustic_threshold = self._caustic_atol + self._caustic_rtol * cofactor_scale

            if not np.isfinite(area_factor) or area_factor < 0.0:
                continue

            eigenvalue_hessians[mode] = hessian
            ray_direction_gradients[mode] = ray_gradient
            area_factors[mode] = area_factor
            caustic_thresholds[mode] = caustic_threshold
            valid_mask[mode] = True
            caustic_candidate_mask[mode] = area_factor <= caustic_threshold

            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                enhancement[mode] = 1.0 / area_factor
                log10_enhancement[mode] = np.log10(enhancement[mode])
            finite_mask[mode] = bool(
                np.isfinite(enhancement[mode]) and np.isfinite(log10_enhancement[mode])
            )

        arrays = (
            eigenvalue_hessians,
            ray_direction_gradients,
            area_factors,
            caustic_thresholds,
            enhancement,
            log10_enhancement,
            valid_mask,
            resolved_mask,
            finite_mask,
            caustic_candidate_mask,
        )
        for array in arrays:
            array.setflags(write=False)

        return DirectionalEnhancementResult(
            group=group,
            eigenvalue_hessians=eigenvalue_hessians,
            ray_direction_gradients=ray_direction_gradients,
            area_factors=area_factors,
            caustic_thresholds=caustic_thresholds,
            enhancement=enhancement,
            log10_enhancement=log10_enhancement,
            valid_mask=valid_mask,
            resolved_mask=resolved_mask,
            finite_mask=finite_mask,
            caustic_candidate_mask=caustic_candidate_mask,
        )


def _validate_tolerance(value: float, name: str) -> float:
    """Validate one non-negative numerical tolerance.

    Parameters
    ----------
    value : float
        Tolerance value.
    name : str
        Parameter name used in error messages.

    Returns
    -------
    float
        Validated tolerance.

    Raises
    ------
    ValueError
        If the value is non-finite or negative.
    """
    converted = float(value)
    if not np.isfinite(converted) or converted < 0.0:
        raise ValueError(f"{name} must be finite and non-negative.")
    return converted


def _validate_rcond(value: float) -> float:
    """Validate a relative pseudoinverse cutoff.

    Parameters
    ----------
    value : float
        Relative singular-value cutoff.

    Returns
    -------
    float
        Validated cutoff in the half-open interval ``[0, 1)``.

    Raises
    ------
    ValueError
        If the cutoff is non-finite or outside ``[0, 1)``.
    """
    converted = float(value)
    if not np.isfinite(converted) or converted < 0.0 or converted >= 1.0:
        raise ValueError("pseudoinverse_rcond must be finite and in [0, 1).")
    return converted


def _normalize_direction(direction: ArrayLike) -> NDArray[np.float64]:
    """Return a read-only normalized Cartesian direction.

    Parameters
    ----------
    direction : array_like
        Cartesian direction with shape ``(3,)``.

    Returns
    -------
    ndarray
        Unit vector with shape ``(3,)``.

    Raises
    ------
    ValueError
        If the input shape, values or norm are invalid.
    """
    vector = np.asarray(direction, dtype=float)
    if vector.shape != (3,):
        raise ValueError("A wave-normal direction must have shape (3,).")
    if not np.all(np.isfinite(vector)):
        raise ValueError("A wave-normal direction must contain finite values.")
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm == 0.0:
        raise ValueError("A wave-normal direction must be non-zero.")
    normalized = np.asarray(vector / norm, dtype=float)
    normalized.setflags(write=False)
    return normalized


def _build_christoffel_hessian(
    reduced_stiffness: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Build the direction-independent Christoffel Hessian.

    Parameters
    ----------
    reduced_stiffness : ndarray
        Cartesian stiffness tensor scaled to km^2 s^-2.

    Returns
    -------
    ndarray
        Hessian with shape ``(3, 3, 3, 3)``.
    """
    hessian = np.transpose(reduced_stiffness, (1, 2, 0, 3)) + np.transpose(
        reduced_stiffness,
        (2, 1, 0, 3),
    )
    return np.asarray(hessian, dtype=float)
