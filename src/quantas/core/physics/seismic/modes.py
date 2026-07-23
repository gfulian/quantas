# -*- coding: utf-8 -*-

"""Acoustic-wave modes and direction-dependent phase results."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import Mapping

import numpy as np
from numpy.typing import NDArray


class WaveMode(str, Enum):
    """Identify the three acoustic modes in ascending phase-speed order.

    The two quasi-shear modes are labelled according to the conventional
    split-shear notation: :math:`V_{S1}` is the faster mode and
    :math:`V_{S2}` is the slower mode.
    """

    V_S2 = "v_s2"
    V_S1 = "v_s1"
    V_P = "v_p"


MODE_ORDER: tuple[WaveMode, WaveMode, WaveMode] = (
    WaveMode.V_S2,
    WaveMode.V_S1,
    WaveMode.V_P,
)

MODE_INDEX: Mapping[WaveMode, int] = MappingProxyType(
    {mode: index for index, mode in enumerate(MODE_ORDER)}
)

MODE_SYMBOLS: Mapping[WaveMode, str] = MappingProxyType(
    {
        WaveMode.V_S2: "V_S2",
        WaveMode.V_S1: "V_S1",
        WaveMode.V_P: "V_P",
    }
)

MODE_DESCRIPTIONS: Mapping[WaveMode, str] = MappingProxyType(
    {
        WaveMode.V_S2: "slow quasi-shear wave",
        WaveMode.V_S1: "fast quasi-shear wave",
        WaveMode.V_P: "quasi-longitudinal wave",
    }
)


class ModePair(str, Enum):
    """Identify adjacent Christoffel eigenvalue pairs."""

    V_S2_V_S1 = "v_s2_v_s1"
    V_S1_V_P = "v_s1_v_p"


MODE_PAIR_ORDER: tuple[ModePair, ModePair] = (
    ModePair.V_S2_V_S1,
    ModePair.V_S1_V_P,
)

MODE_PAIR_INDEX: Mapping[ModePair, int] = MappingProxyType(
    {pair: index for index, pair in enumerate(MODE_PAIR_ORDER)}
)

MODE_PAIR_SYMBOLS: Mapping[ModePair, str] = MappingProxyType(
    {
        ModePair.V_S2_V_S1: "V_S2-V_S1",
        ModePair.V_S1_V_P: "V_S1-V_P",
    }
)


class PolarizationBranch(str, Enum):
    """Identify continuity branches used for polarization tracking."""

    SHEAR_A = "shear_a"
    SHEAR_B = "shear_b"
    P = "p"


BRANCH_ORDER: tuple[
    PolarizationBranch,
    PolarizationBranch,
    PolarizationBranch,
] = (
    PolarizationBranch.SHEAR_A,
    PolarizationBranch.SHEAR_B,
    PolarizationBranch.P,
)

BRANCH_INDEX: Mapping[PolarizationBranch, int] = MappingProxyType(
    {branch: index for index, branch in enumerate(BRANCH_ORDER)}
)


@dataclass(frozen=True, slots=True)
class PhaseModeResult:
    """Phase solution for one acoustic mode and one propagation direction.

    Parameters
    ----------
    mode : WaveMode
        Acoustic-wave mode.
    eigenvalue : float
        Christoffel eigenvalue in km^2 s^-2.
    phase_speed : float
        Phase speed in km s^-1. Invalid modes contain ``nan``.
    polarization : ndarray
        Unit polarization axis with shape ``(3,)``.
    eigenvalue_gap : float
        Smallest absolute eigenvalue separation from an adjacent mode in
        km^2 s^-2.
    relative_eigenvalue_gap : float
        ``eigenvalue_gap`` divided by the largest absolute eigenvalue of the
        directional Christoffel matrix.
    valid : bool
        Whether the eigenvalue is physically admissible under the selected
        tolerance.
    clamped : bool
        Whether a small negative eigenvalue was set to zero before taking the
        square root.
    degenerate : bool
        Whether the mode belongs to a numerically degenerate eigenspace.
    """

    mode: WaveMode
    eigenvalue: float
    phase_speed: float
    polarization: NDArray[np.float64]
    eigenvalue_gap: float
    relative_eigenvalue_gap: float
    valid: bool
    clamped: bool
    degenerate: bool


@dataclass(frozen=True, slots=True)
class DirectionalPhaseResult:
    """Complete phase solution for one wave-normal direction.

    Parameters
    ----------
    direction : ndarray
        Normalized wave-normal direction with shape ``(3,)``.
    christoffel : ndarray
        Symmetric Christoffel matrix in km^2 s^-2.
    eigenvalues : ndarray
        Ascending Christoffel eigenvalues in km^2 s^-2, ordered as
        ``V_S2``, ``V_S1`` and ``V_P``.
    polarizations : ndarray
        Row-oriented polarization axes with shape ``(3, 3)``.
    phase_speeds : ndarray
        Phase speeds in km s^-1.
    eigenvalue_gaps : ndarray
        Adjacent absolute eigenvalue separations with shape ``(2,)``.
    relative_eigenvalue_gaps : ndarray
        Adjacent eigenvalue separations normalized by the largest absolute
        directional eigenvalue.
    mode_eigenvalue_gaps : ndarray
        Smallest adjacent absolute separation for each mode.
    mode_relative_eigenvalue_gaps : ndarray
        Smallest adjacent relative separation for each mode.
    valid_mask : ndarray
        Boolean validity mask in mode order.
    clamped_mask : ndarray
        Boolean mask identifying small negative eigenvalues set to zero.
    degeneracy_mask : ndarray
        Boolean mask identifying modes in a numerically degenerate eigenspace.
    eigenvalue_threshold : float
        Absolute negative-eigenvalue threshold in km^2 s^-2.
    degeneracy_threshold : float
        Absolute eigenvalue-gap threshold in km^2 s^-2.
    """

    direction: NDArray[np.float64]
    christoffel: NDArray[np.float64]
    eigenvalues: NDArray[np.float64]
    polarizations: NDArray[np.float64]
    phase_speeds: NDArray[np.float64]
    eigenvalue_gaps: NDArray[np.float64]
    relative_eigenvalue_gaps: NDArray[np.float64]
    mode_eigenvalue_gaps: NDArray[np.float64]
    mode_relative_eigenvalue_gaps: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    clamped_mask: NDArray[np.bool_]
    degeneracy_mask: NDArray[np.bool_]
    eigenvalue_threshold: float
    degeneracy_threshold: float

    def for_mode(self, mode: WaveMode | str) -> PhaseModeResult:
        """Return the phase solution associated with one acoustic mode.

        Parameters
        ----------
        mode : WaveMode or str
            Mode enum or one of ``"v_s2"``, ``"v_s1"`` and ``"v_p"``.

        Returns
        -------
        PhaseModeResult
            Selected mode result.

        Raises
        ------
        ValueError
            If ``mode`` is not a supported acoustic mode.
        """
        resolved = WaveMode(mode)
        index = MODE_INDEX[resolved]
        return PhaseModeResult(
            mode=resolved,
            eigenvalue=float(self.eigenvalues[index]),
            phase_speed=float(self.phase_speeds[index]),
            polarization=self.polarizations[index],
            eigenvalue_gap=float(self.mode_eigenvalue_gaps[index]),
            relative_eigenvalue_gap=float(self.mode_relative_eigenvalue_gaps[index]),
            valid=bool(self.valid_mask[index]),
            clamped=bool(self.clamped_mask[index]),
            degenerate=bool(self.degeneracy_mask[index]),
        )

    @property
    def invalid_mask(self) -> NDArray[np.bool_]:
        """Return a boolean mask for physically invalid modes."""
        invalid = np.logical_not(self.valid_mask)
        invalid.setflags(write=False)
        return invalid

    @property
    def is_fully_valid(self) -> bool:
        """Return whether all three acoustic modes are physically valid."""
        return bool(np.all(self.valid_mask))

    @property
    def has_degeneracy(self) -> bool:
        """Return whether at least one mode is numerically degenerate."""
        return bool(np.any(self.degeneracy_mask))
