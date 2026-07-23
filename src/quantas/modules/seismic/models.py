# -*- coding: utf-8 -*-

"""Passive input, option and result models for seismic workflows."""

from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field
from pathlib import Path
from typing import Any

import numpy as np

from quantas.core.geometry import Hemisphere, SphericalGrid, TensorRotation
from quantas.core.physics.elasticity import ElasticAverages, StabilityResult
from quantas.core.physics.seismic import (
    IsotropicSeismicVelocities,
    SamplingLevel,
    SeismicFieldResult,
)


@dataclass(slots=True)
class SeismicInput:
    """Input data for acoustic-wave propagation in an elastic crystal.

    Parameters
    ----------
    jobname : str, optional
        Name or short description of the calculation.
    stiffness : ndarray or None, optional
        Elastic stiffness matrix in Voigt notation with shape ``(6, 6)`` and
        values in GPa.
    density : float, optional
        Material density in kg m^-3.
    source : str, Path or None, optional
        Source from which the input was obtained.
    raw : str or None, optional
        Original text input when available.
    """

    jobname: str = "Unknown"
    stiffness: np.ndarray | None = None
    density: float = 0.0
    source: str | Path | None = None
    raw: str | None = None


@dataclass(slots=True)
class SeismicOptions:
    """Scientific options controlling a seismic field calculation.

    Parameters
    ----------
    ntheta, nphi : int, optional
        Polar and azimuthal sampling counts. The azimuthal seam is not
        duplicated.
    hemisphere : Hemisphere, optional
        Polar domain sampled by the workflow.
    level : SamplingLevel, optional
        Highest acoustic quantity to calculate.
    batch_size : int, optional
        Maximum number of directions evaluated in one NumPy batch.
    track_polarization_axes : bool, optional
        Whether to align axial polarizations along a deterministic grid path.
    eigenvalue_rtol, eigenvalue_atol : float, optional
        Tolerances for small negative Christoffel eigenvalues.
    degeneracy_rtol, degeneracy_atol : float, optional
        Tolerances for acoustic-mode degeneracies.
    pseudoinverse_rcond : float, optional
        Relative cutoff used in analytical eigenvalue Hessians.
    caustic_rtol, caustic_atol : float, optional
        Tolerances used to identify possible caustics.
    rotation : TensorRotation or None, optional
        Optional source-to-analysis tensor-component transformation.
    """

    ntheta: int = 91
    nphi: int = 181
    hemisphere: Hemisphere = Hemisphere.UPPER
    level: SamplingLevel = SamplingLevel.ENHANCEMENT
    batch_size: int = 512
    track_polarization_axes: bool = True
    eigenvalue_rtol: float = 1.0e-10
    eigenvalue_atol: float = 1.0e-12
    degeneracy_rtol: float = 1.0e-8
    degeneracy_atol: float = 1.0e-10
    pseudoinverse_rcond: float = 1.0e-10
    caustic_rtol: float = 1.0e-10
    caustic_atol: float = 1.0e-12
    rotation: TensorRotation | None = None

    def __post_init__(self) -> None:
        """Normalize enumerated options to their public enum types."""
        self.hemisphere = Hemisphere(self.hemisphere)
        self.level = SamplingLevel(self.level)


@dataclass(slots=True)
class SeismicResult:
    """Results of a sampled acoustic-wave calculation.

    Parameters
    ----------
    jobname : str
        Name or short description of the calculation.
    density : float
        Material density in kg m^-3.
    stiffness : ndarray
        Elastic stiffness matrix in GPa.
    stability : StabilityResult
        Positive-definiteness result for the stiffness matrix.
    averages : ElasticAverages
        Voigt, Reuss and Hill elastic averages.
    isotropic_velocities : IsotropicSeismicVelocities
        Hill-average shear and compressional reference velocities.
    grid : SphericalGrid
        Sampled angular grid.
    field : SeismicFieldResult
        Sampled phase, group and enhancement fields.
    metadata : dict, optional
        Additional numerical diagnostics.
    """

    jobname: str
    density: float
    stiffness: np.ndarray
    stability: StabilityResult
    averages: ElasticAverages
    isotropic_velocities: IsotropicSeismicVelocities
    grid: SphericalGrid
    field: SeismicFieldResult
    metadata: dict[str, Any] = dataclass_field(default_factory=dict)
