# -*- coding: utf-8 -*-

r"""Directional acoustic averages and Kieffer cutoff construction."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import NDArray

from .christoffel import ChristoffelSolver
from .sampling import ProgressCallback, solve_phase_directions


@dataclass(frozen=True, slots=True)
class AcousticVelocityAverage:
    """Spherical inverse-cube averages of three acoustic phase velocities.

    Parameters
    ----------
    effective_velocities : ndarray
        Effective velocities in km s^-1, ordered as slow shear, fast shear,
        and quasi-longitudinal.
    inverse_cube_means : ndarray
        Weighted spherical means of ``v^-3`` in s^3 km^-3.
    relative_errors : ndarray
        Relative changes from the coarser quadrature.
    mu_order, phi_order : int
        Final Gauss-Legendre and periodic azimuthal quadrature orders.
    direction_count : int
        Number of directions in the final quadrature.
    degenerate_direction_count : int
        Directions containing at least one degenerate mode.
    clamped_mode_count : int
        Small negative eigenvalues clamped by the Christoffel solver.
    """

    effective_velocities: NDArray[np.float64]
    inverse_cube_means: NDArray[np.float64]
    relative_errors: NDArray[np.float64]
    mu_order: int
    phi_order: int
    direction_count: int
    degenerate_direction_count: int
    clamped_mode_count: int


@dataclass(frozen=True, slots=True)
class KiefferCutoffFrequencies:
    """Kieffer acoustic cutoffs derived from effective velocities.

    Parameters
    ----------
    primitive_volume_angstrom3 : float
        Primitive-cell volume in cubic angstrom.
    brillouin_radius_m_inv : float
        Radius of the equal-volume spherical Brillouin zone in m^-1.
    angular_frequencies_rad_s : ndarray
        Angular cutoff frequencies in rad s^-1.
    frequencies_hz : ndarray
        Ordinary cutoff frequencies in hertz.
    wavenumbers_cm1 : ndarray
        Spectroscopic cutoff wavenumbers in cm^-1.
    """

    primitive_volume_angstrom3: float
    brillouin_radius_m_inv: float
    angular_frequencies_rad_s: NDArray[np.float64]
    frequencies_hz: NDArray[np.float64]
    wavenumbers_cm1: NDArray[np.float64]


def _quadrature(
    order_mu: int, order_phi: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return directions and normalized solid-angle weights."""
    if isinstance(order_mu, bool) or int(order_mu) != order_mu or order_mu < 2:
        raise ValueError("mu_order must be an integer not smaller than 2")
    if isinstance(order_phi, bool) or int(order_phi) != order_phi or order_phi < 4:
        raise ValueError("phi_order must be an integer not smaller than 4")
    mu, mu_weights = np.polynomial.legendre.leggauss(int(order_mu))
    phi = 2.0 * math.pi * np.arange(int(order_phi), dtype=np.float64) / order_phi
    radius = np.sqrt(np.maximum(0.0, 1.0 - mu * mu))
    directions = np.stack(
        (
            (radius[:, None] * np.cos(phi)[None, :]).ravel(),
            (radius[:, None] * np.sin(phi)[None, :]).ravel(),
            np.broadcast_to(mu[:, None], (mu.size, phi.size)).ravel(),
        ),
        axis=1,
    )
    weights = np.repeat(mu_weights / (2.0 * order_phi), int(order_phi))
    return directions, weights


def _average_once(
    solver: ChristoffelSolver,
    mu_order: int,
    phi_order: int,
    batch_size: int,
    progress_callback: ProgressCallback | None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], int, int]:
    """Evaluate one spherical quadrature level."""
    directions, weights = _quadrature(mu_order, phi_order)
    phase = solve_phase_directions(
        solver,
        directions,
        batch_size=batch_size,
        progress_callback=progress_callback,
    )
    speeds = phase.phase_speeds
    if not np.all(phase.valid_mask) or not np.all(np.isfinite(speeds)):
        raise ValueError("invalid Christoffel modes prevent acoustic averaging")
    if np.any(speeds <= 0.0):
        raise ValueError("acoustic phase velocities must be strictly positive")
    inverse_cube = np.sum(weights[:, None] * speeds**-3, axis=0)
    effective = inverse_cube ** (-1.0 / 3.0)
    degeneracies = int(np.count_nonzero(np.any(phase.degeneracy_mask, axis=1)))
    clamped = int(np.count_nonzero(phase.clamped_mask))
    return effective, inverse_cube, degeneracies, clamped


def average_acoustic_phase_velocities(
    solver: ChristoffelSolver,
    *,
    mu_order: int = 12,
    phi_order: int = 24,
    refinement_factor: int = 2,
    batch_size: int = 512,
    progress_callback: ProgressCallback | None = None,
) -> AcousticVelocityAverage:
    r"""Calculate Kieffer's inverse-cube spherical velocity averages.

    The three locally speed-ordered Christoffel modes are averaged as

    .. math::

        u_i = \left[\frac{1}{4\pi}\int_{4\pi}v_i^{-3}\,d\Omega\right]^{-1/3}.

    Gauss-Legendre quadrature is used in ``mu = cos(theta)`` and a periodic
    uniform rule in azimuth.  The returned value uses the refined grid; its
    error estimate is the relative change from the requested coarse grid.

    Parameters
    ----------
    solver : ChristoffelSolver
        Acoustic phase-velocity solver.
    mu_order, phi_order : int, optional
        Coarse quadrature orders.
    refinement_factor : int, optional
        Integer multiplier applied to both final quadrature orders.
    batch_size : int, optional
        Maximum number of directions solved in one batch.
    progress_callback : callable or None, optional
        Callback forwarded during the final refined sampling.

    Returns
    -------
    AcousticVelocityAverage
        Effective velocities and convergence diagnostics.

    Raises
    ------
    TypeError
        If ``solver`` is invalid.
    ValueError
        If quadrature parameters or phase solutions are invalid.
    """
    if not isinstance(solver, ChristoffelSolver):
        raise TypeError("acoustic averaging requires a ChristoffelSolver")
    if (
        isinstance(refinement_factor, bool)
        or int(refinement_factor) != refinement_factor
        or refinement_factor < 1
    ):
        raise ValueError("refinement_factor must be a positive integer")
    coarse, _coarse_inverse, _coarse_degenerate, _coarse_clamped = _average_once(
        solver, mu_order, phi_order, batch_size, None
    )
    final_mu = int(mu_order) * int(refinement_factor)
    final_phi = int(phi_order) * int(refinement_factor)
    effective, inverse_cube, degeneracies, clamped = _average_once(
        solver, final_mu, final_phi, batch_size, progress_callback
    )
    relative_errors = np.abs(effective - coarse) / effective
    for array in (effective, inverse_cube, relative_errors):
        array.setflags(write=False)
    return AcousticVelocityAverage(
        effective_velocities=effective,
        inverse_cube_means=inverse_cube,
        relative_errors=relative_errors,
        mu_order=final_mu,
        phi_order=final_phi,
        direction_count=final_mu * final_phi,
        degenerate_direction_count=degeneracies,
        clamped_mode_count=clamped,
    )


def kieffer_cutoff_frequencies(
    effective_velocities_km_s: NDArray[np.float64],
    primitive_volume_angstrom3: float,
) -> KiefferCutoffFrequencies:
    r"""Convert effective velocities to sine-wave acoustic cutoffs.

    Parameters
    ----------
    effective_velocities_km_s : array_like
        Three positive effective velocities in km s^-1.
    primitive_volume_angstrom3 : float
        Primitive-cell volume in cubic angstrom.

    Returns
    -------
    KiefferCutoffFrequencies
        Brillouin-sphere radius and cutoffs in rad/s, Hz, and cm^-1.

    Raises
    ------
    ValueError
        If velocities or volume are invalid.
    """
    velocities = np.asarray(effective_velocities_km_s, dtype=np.float64)
    volume = float(primitive_volume_angstrom3)
    if (
        velocities.shape != (3,)
        or not np.all(np.isfinite(velocities))
        or np.any(velocities <= 0.0)
    ):
        raise ValueError("effective velocities must be three finite positive values")
    if not np.isfinite(volume) or volume <= 0.0:
        raise ValueError("primitive volume must be finite and positive")
    radius = (6.0 * math.pi**2 / (volume * 1.0e-30)) ** (1.0 / 3.0)
    angular = (2.0 / math.pi) * (velocities * 1.0e3) * radius
    frequencies = angular / (2.0 * math.pi)
    wavenumbers = frequencies / (299_792_458.0 * 100.0)
    for array in (angular, frequencies, wavenumbers):
        array.setflags(write=False)
    return KiefferCutoffFrequencies(
        primitive_volume_angstrom3=volume,
        brillouin_radius_m_inv=radius,
        angular_frequencies_rad_s=angular,
        frequencies_hz=frequencies,
        wavenumbers_cm1=wavenumbers,
    )
