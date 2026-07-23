# -*- coding: utf-8 -*-

"""Passive input, option, and result models for elasticity workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from quantas.core.geometry import TensorRotation
from quantas.core.physics.elasticity import (
    DirectionalExtrema,
    ElasticAverages,
    ElasticitySurfaceCollection,
    ElasticSurfaceProperty,
    StabilityResult,
    SUPPORTED_SURFACE_PROPERTIES,
)


@dataclass(slots=True)
class ElasticityInput:
    """Input data for a second-order elasticity calculation.

    Parameters
    ----------
    jobname : str, optional
        Name or short description of the calculation.
    stiffness : ndarray or None, optional
        Elastic stiffness matrix in Voigt notation, with shape ``(6, 6)`` and
        values in GPa.
    source : str or Path or None, optional
        Source file from which the input was read.
    """

    jobname: str = "Unknown"
    stiffness: np.ndarray | None = None
    source: str | Path | None = None


@dataclass(slots=True)
class ElasticitySurfaceOptions:
    """Options controlling three-dimensional elasticity sampling.

    Parameters
    ----------
    ntheta : int, optional
        Number of polar samples, including both poles.
    nphi : int, optional
        Number of azimuthal samples without duplicating ``2π``.
    properties : tuple of str, optional
        Selected property families. Supported values are ``young``,
        ``compressibility``, ``shear`` and ``poisson``.
    batch_size : int, optional
        Maximum number of directions evaluated in one NumPy batch.

    Raises
    ------
    ValueError
        If a grid, property, or batch option is invalid.
    """

    ntheta: int = 61
    nphi: int = 121
    properties: tuple[ElasticSurfaceProperty, ...] = SUPPORTED_SURFACE_PROPERTIES
    batch_size: int = 65536

    def __post_init__(self) -> None:
        """Validate and normalize surface-sampling options."""
        if self.ntheta < 2:
            raise ValueError("ntheta must be at least 2.")
        if self.nphi < 3:
            raise ValueError("nphi must be at least 3.")
        if self.batch_size < 1:
            raise ValueError("batch_size must be positive.")
        selected: list[ElasticSurfaceProperty] = []
        for property_name in self.properties:
            if property_name not in SUPPORTED_SURFACE_PROPERTIES:
                raise ValueError(f"unsupported 3D elasticity property: {property_name}")
            if property_name not in selected:
                selected.append(property_name)
        self.properties = tuple(selected)


@dataclass(slots=True)
class ElasticityOptions:
    """Scientific options controlling a persisted elasticity calculation.

    Parameters
    ----------
    pressure_unit : str, optional
        Unit used for stiffness and pressure-like quantities. The current
        implementation supports GPa.
    calculate_2d : bool, optional
        Whether to calculate and persist properties on the three principal
        Cartesian planes.
    ntheta : int, optional
        Number of angular points in each closed two-dimensional curve.
    calculate_3d : bool, optional
        Whether to calculate and persist three-dimensional directional data.
    surface_options : ElasticitySurfaceOptions or None, optional
        Grid, property, and batch settings used when ``calculate_3d`` is true.
    rotation : TensorRotation or None, optional
        Optional source-to-analysis tensor-component transformation.
    """

    pressure_unit: str = "GPa"
    calculate_2d: bool = False
    ntheta: int = 361
    calculate_3d: bool = False
    surface_options: ElasticitySurfaceOptions | None = None
    rotation: TensorRotation | None = None

    def __post_init__(self) -> None:
        """Validate persisted elasticity options."""
        if self.ntheta < 2:
            raise ValueError("ntheta must be at least 2.")
        if self.calculate_3d and self.surface_options is None:
            self.surface_options = ElasticitySurfaceOptions()


@dataclass(slots=True)
class ElasticityResult:
    """Results of a second-order elasticity calculation.

    Parameters
    ----------
    jobname : str, optional
        Name or short description of the calculation.
    crystal_system : str or None, optional
        Crystal system inferred from the stiffness matrix.
    stiffness, compliance : ndarray or None, optional
        Elastic matrices in Voigt notation.
    averages : ElasticAverages or None, optional
        Voigt, Reuss, and Hill estimates.
    stability : StabilityResult or None, optional
        Positive-definiteness check and stiffness eigenvalues.
    variations : dict, optional
        Directional extrema for elastic properties.
    properties_2d : dict, optional
        Principal-plane directional data requested by the workflow.
    properties_3d : ElasticitySurfaceCollection or None, optional
        Persisted three-dimensional directional surfaces requested by the
        workflow. Plot-only calculations may remain transient instead.
    metadata : dict, optional
        Additional workflow metadata.
    """

    jobname: str = "Unknown"
    crystal_system: str | None = None
    stiffness: np.ndarray | None = None
    compliance: np.ndarray | None = None
    averages: ElasticAverages | None = None
    stability: StabilityResult | None = None
    variations: dict[str, DirectionalExtrema] = field(default_factory=dict)
    properties_2d: dict[str, Any] = field(default_factory=dict)
    properties_3d: ElasticitySurfaceCollection | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def has_2d_data(self) -> bool:
        """Return whether principal-plane directional data are available."""
        return bool(self.properties_2d)

    def has_3d_data(self) -> bool:
        """Return whether persisted three-dimensional data are available."""
        return self.properties_3d is not None and bool(self.properties_3d.surfaces)

    def add_variation(self, name: str, value: DirectionalExtrema) -> None:
        """Store directional-extrema data under a property name."""
        self.variations[name] = value

    def add_2d_data(self, plane: str, property_name: str, value: Any) -> None:
        """Store one two-dimensional property for a principal plane."""
        self.properties_2d.setdefault(plane, {})[property_name] = value
