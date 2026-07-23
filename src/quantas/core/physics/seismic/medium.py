# -*- coding: utf-8 -*-

"""Elastic media used by acoustic-wave propagation solvers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from quantas.core.physics.elasticity.tensor import ElasticTensor


@dataclass(frozen=True, slots=True)
class ElasticMedium:
    """Associate an elastic stiffness tensor with a material density.

    Parameters
    ----------
    elastic_tensor : ElasticTensor
        Elastic stiffness tensor expressed in GPa.
    density : float
        Material density in kg m^-3.

    Raises
    ------
    TypeError
        If ``elastic_tensor`` is not an :class:`ElasticTensor` instance.
    ValueError
        If ``density`` is non-finite or not strictly positive.
    """

    elastic_tensor: ElasticTensor
    density: float

    def __post_init__(self) -> None:
        """Validate the tensor and normalize the density to ``float``."""
        if not isinstance(self.elastic_tensor, ElasticTensor):
            raise TypeError("ElasticMedium requires an ElasticTensor instance.")
        density = float(self.density)
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError("Density must be finite and positive.")
        object.__setattr__(self, "density", density)
