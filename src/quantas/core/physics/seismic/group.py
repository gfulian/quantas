# -*- coding: utf-8 -*-

"""Direction-dependent group velocities and energy-flow directions."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .modes import MODE_INDEX, DirectionalPhaseResult, PhaseModeResult, WaveMode


@dataclass(frozen=True, slots=True)
class GroupModeResult:
    """Group-velocity solution for one acoustic mode.

    Parameters
    ----------
    phase : PhaseModeResult
        Phase solution associated with the selected mode.
    eigenvalue_gradient : ndarray
        Cartesian gradient of the Christoffel eigenvalue in km^2 s^-2.
    group_velocity : ndarray
        Group-velocity vector in km s^-1.
    group_speed : float
        Magnitude of ``group_velocity`` in km s^-1.
    ray_direction : ndarray
        Unit energy-flow direction. Undefined results contain ``nan``.
    power_flow_angle : float
        Angle between the wave-normal and ray directions in radians.
    valid : bool
        Whether the group quantities are numerically defined.
    resolved : bool
        Whether the mode is non-degenerate and uniquely resolved.
    """

    phase: PhaseModeResult
    eigenvalue_gradient: NDArray[np.float64]
    group_velocity: NDArray[np.float64]
    group_speed: float
    ray_direction: NDArray[np.float64]
    power_flow_angle: float
    valid: bool
    resolved: bool


@dataclass(frozen=True, slots=True)
class DirectionalGroupResult:
    """Complete phase and group solution for one wave-normal direction.

    Parameters
    ----------
    phase : DirectionalPhaseResult
        Pointwise Christoffel phase solution.
    eigenvalue_gradients : ndarray
        Eigenvalue gradients with shape ``(3, 3)`` in mode order.
    group_velocities : ndarray
        Group-velocity vectors with shape ``(3, 3)`` in km s^-1.
    group_speeds : ndarray
        Group-speed magnitudes with shape ``(3,)`` in km s^-1.
    ray_directions : ndarray
        Unit energy-flow directions with shape ``(3, 3)``.
    power_flow_angles : ndarray
        Power-flow angles in radians with shape ``(3,)``.
    valid_mask : ndarray
        Boolean mask identifying numerically defined group quantities.
    resolved_mask : ndarray
        Boolean mask identifying non-degenerate group solutions.
    """

    phase: DirectionalPhaseResult
    eigenvalue_gradients: NDArray[np.float64]
    group_velocities: NDArray[np.float64]
    group_speeds: NDArray[np.float64]
    ray_directions: NDArray[np.float64]
    power_flow_angles: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    resolved_mask: NDArray[np.bool_]

    @property
    def direction(self) -> NDArray[np.float64]:
        """Return the normalized wave-normal direction."""
        return self.phase.direction

    @property
    def phase_speeds(self) -> NDArray[np.float64]:
        """Return the phase speeds in km s^-1."""
        return self.phase.phase_speeds

    @property
    def polarizations(self) -> NDArray[np.float64]:
        """Return the row-oriented polarization axes."""
        return self.phase.polarizations

    def for_mode(self, mode: WaveMode | str) -> GroupModeResult:
        """Return the phase and group solution for one acoustic mode.

        Parameters
        ----------
        mode : WaveMode or str
            Mode enum or one of ``"v_s2"``, ``"v_s1"`` and ``"v_p"``.

        Returns
        -------
        GroupModeResult
            Selected mode result.

        Raises
        ------
        ValueError
            If ``mode`` is not a supported acoustic mode.
        """
        resolved = WaveMode(mode)
        index = MODE_INDEX[resolved]
        return GroupModeResult(
            phase=self.phase.for_mode(resolved),
            eigenvalue_gradient=self.eigenvalue_gradients[index],
            group_velocity=self.group_velocities[index],
            group_speed=float(self.group_speeds[index]),
            ray_direction=self.ray_directions[index],
            power_flow_angle=float(self.power_flow_angles[index]),
            valid=bool(self.valid_mask[index]),
            resolved=bool(self.resolved_mask[index]),
        )

    @property
    def is_fully_valid(self) -> bool:
        """Return whether all three group solutions are numerically defined."""
        return bool(np.all(self.valid_mask))

    @property
    def is_fully_resolved(self) -> bool:
        """Return whether all three group solutions are uniquely resolved."""
        return bool(np.all(self.resolved_mask))
