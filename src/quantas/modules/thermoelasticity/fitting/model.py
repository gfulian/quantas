# -*- coding: utf-8 -*-

"""Fitting services for cold quasi-static elastic components.

This module owns no file or frontend logic. It combines the general Quantas
fitting infrastructure with the Eulerian finite-strain relations implemented in
:mod:`quantas.core.physics.elasticity.quasistatic`.
"""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.math.fitting import (
    BaseFitModel,
)
from quantas.core.physics.elasticity import (
    cold_finite_strain_component,
    cold_finite_strain_component_jacobian,
)


FloatArray: TypeAlias = NDArray[np.float64]


class ColdFiniteStrainComponentModel(BaseFitModel):
    """Two-parameter cold finite-strain model for one Voigt component.

    Parameters
    ----------
    reference_volume : float
        Fixed static-EOS reference volume in angstrom cubed.
    bulk_modulus : float
        Fixed static-EOS bulk modulus in GPa.
    bulk_modulus_derivative : float
        Fixed static-EOS first pressure derivative.
    wallace_delta : float
        Matching Wallace pre-stress coefficient.
    order : {2, 3}
        Finite-strain truncation order.
    label : str
        Component label used in metadata.
    """

    def __init__(
        self,
        *,
        reference_volume: float,
        bulk_modulus: float,
        bulk_modulus_derivative: float,
        wallace_delta: float,
        order: int,
        label: str,
    ) -> None:
        self.reference_volume = float(reference_volume)
        self.bulk_modulus = float(bulk_modulus)
        self.bulk_modulus_derivative = float(bulk_modulus_derivative)
        self.wallace_delta = float(wallace_delta)
        self.order = int(order)
        self.label = str(label)

    @property
    def name(self) -> str:
        """Return the stable model identifier."""
        return "cold_eulerian_finite_strain_component"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return ``C0`` and ``Cprime`` in fitting order."""
        return ("C0", "Cprime")

    def evaluate(self, x: ArrayLike, parameters: ArrayLike) -> FloatArray:
        """Evaluate the component at sampled volumes."""
        values = np.asarray(parameters, dtype=np.float64)
        if values.shape != (2,):
            raise ValueError("component parameters must contain C0 and Cprime")
        return cold_finite_strain_component(
            x,
            reference_volume=self.reference_volume,
            bulk_modulus=self.bulk_modulus,
            bulk_modulus_derivative=self.bulk_modulus_derivative,
            reference_component=float(values[0]),
            component_pressure_derivative=float(values[1]),
            wallace_delta=self.wallace_delta,
            order=self.order,  # type: ignore[arg-type]
        )

    def initial_guess(self, x: ArrayLike, y: ArrayLike) -> FloatArray:
        """Return the exact linear least-squares estimate for model parameters."""
        volumes = np.asarray(x, dtype=np.float64)
        observed = np.asarray(y, dtype=np.float64)
        zero = self.evaluate(volumes, np.asarray([0.0, 0.0], dtype=np.float64))
        column_c0 = (
            self.evaluate(volumes, np.asarray([1.0, 0.0], dtype=np.float64)) - zero
        )
        column_cp = (
            self.evaluate(volumes, np.asarray([0.0, 1.0], dtype=np.float64)) - zero
        )
        design = np.column_stack((column_c0, column_cp))
        parameters, _, rank, _ = np.linalg.lstsq(design, observed - zero, rcond=None)
        if rank < 2 or not np.all(np.isfinite(parameters)):
            reference = float(
                observed[np.argmin(np.abs(volumes - self.reference_volume))]
            )
            return np.asarray([reference, 1.0], dtype=np.float64)
        return np.asarray(parameters, dtype=np.float64)

    def derivative_x(self, x: ArrayLike, parameters: ArrayLike) -> FloatArray:
        """Return the analytical volume derivative of the component."""
        values = np.asarray(parameters, dtype=np.float64)
        jacobian = cold_finite_strain_component_jacobian(
            x,
            reference_volume=self.reference_volume,
            bulk_modulus=self.bulk_modulus,
            bulk_modulus_derivative=self.bulk_modulus_derivative,
            reference_component=float(values[0]),
            component_pressure_derivative=float(values[1]),
            wallace_delta=self.wallace_delta,
            order=self.order,  # type: ignore[arg-type]
        )
        return np.asarray(jacobian[..., 5], dtype=np.float64)

    def metadata(self) -> dict[str, object]:
        """Return fixed-EOS and component metadata."""
        return {
            **super().metadata(),
            "component": self.label,
            "reference_volume": self.reference_volume,
            "bulk_modulus_GPa": self.bulk_modulus,
            "bulk_modulus_derivative": self.bulk_modulus_derivative,
            "wallace_delta": self.wallace_delta,
            "wallace_delta_role": (
                "analytical finite-strain coefficient; no additional pressure "
                "correction of observed CRYSTAL tensors"
            ),
            "finite_strain_order": self.order,
        }


__all__ = ["ColdFiniteStrainComponentModel"]
