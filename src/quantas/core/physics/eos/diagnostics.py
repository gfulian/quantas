# -*- coding: utf-8 -*-

"""Finite-strain and normalized-pressure diagnostics for pressure EOS models.

The transformations in this module are pure numerical operations.  They do not
know about fitted datasets, renderers, archives, or frontends.  Only EOS
families with a scientifically established normalized-pressure representation
are supported: Birch--Murnaghan (Eulerian strain), Natural Strain
(Poirier--Tarantola), and Vinet.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TypeAlias

import numpy as np

from .spec import EOSFamily, EOSModel, parse_eos_model

ArrayLike: TypeAlias = np.ndarray | float


class EOSStrainFamily(str, Enum):
    """Canonical finite-strain representations used for EOS diagnostics."""

    EULERIAN = "eulerian"
    NATURAL = "natural"
    VINET = "vinet"


@dataclass(frozen=True, slots=True)
class EOSStrainTransform:
    """Finite strain and normalized pressure for one EOS family.

    Parameters
    ----------
    family : EOSStrainFamily
        Diagnostic strain representation.
    strain, normalized_pressure : ndarray
        Transformed values.  ``normalized_pressure`` is ``NaN`` where the
        transformation is singular, normally at the exact reference state.
    valid : ndarray
        Boolean mask identifying finite normalized-pressure values.
    dstrain_dcoordinate, dstrain_dreference : ndarray
        First derivatives with respect to the observed and reference
        volume-like coordinates.
    dnormalized_dpressure, dnormalized_dcoordinate, dnormalized_dreference : ndarray
        First derivatives used for uncertainty propagation.
    metadata : dict
        Stable descriptions of the transformation and parameter conventions.
    """

    family: EOSStrainFamily
    strain: np.ndarray
    normalized_pressure: np.ndarray
    valid: np.ndarray
    dstrain_dcoordinate: np.ndarray
    dstrain_dreference: np.ndarray
    dnormalized_dpressure: np.ndarray
    dnormalized_dcoordinate: np.ndarray
    dnormalized_dreference: np.ndarray
    metadata: dict[str, object]

    def __post_init__(self) -> None:
        arrays = (
            self.strain,
            self.normalized_pressure,
            self.valid,
            self.dstrain_dcoordinate,
            self.dstrain_dreference,
            self.dnormalized_dpressure,
            self.dnormalized_dcoordinate,
            self.dnormalized_dreference,
        )
        shape = np.asarray(self.strain).shape
        for value in arrays:
            if np.asarray(value).shape != shape:
                raise ValueError("EOS strain diagnostic arrays must share one shape")
        object.__setattr__(
            self, "strain", np.asarray(self.strain, dtype=np.float64).copy()
        )
        object.__setattr__(
            self,
            "normalized_pressure",
            np.asarray(self.normalized_pressure, dtype=np.float64).copy(),
        )
        object.__setattr__(self, "valid", np.asarray(self.valid, dtype=np.bool_).copy())
        for name in (
            "dstrain_dcoordinate",
            "dstrain_dreference",
            "dnormalized_dpressure",
            "dnormalized_dcoordinate",
            "dnormalized_dreference",
        ):
            object.__setattr__(
                self,
                name,
                np.asarray(getattr(self, name), dtype=np.float64).copy(),
            )
        object.__setattr__(self, "metadata", dict(self.metadata))


class PressureEOSDiagnostics:
    """Evaluate finite-strain and normalized-pressure transformations."""

    def strain_family(self, eos: str | EOSModel) -> EOSStrainFamily:
        """Return the diagnostic strain family supported by ``eos``.

        Raises
        ------
        NotImplementedError
            If the EOS has no established normalized-pressure representation
            in Quantas.
        """
        model = parse_eos_model(eos)
        if model.family is EOSFamily.BIRCH_MURNAGHAN:
            return EOSStrainFamily.EULERIAN
        if model.family is EOSFamily.NATURAL_STRAIN:
            return EOSStrainFamily.NATURAL
        if model.family is EOSFamily.VINET:
            return EOSStrainFamily.VINET
        raise NotImplementedError(
            f"normalized-pressure diagnostics are not defined for {model.tag}"
        )

    def transform(
        self,
        eos: str | EOSModel,
        pressure: ArrayLike,
        coordinate: ArrayLike,
        reference_coordinate: float,
        *,
        singular_tolerance: float | None = None,
    ) -> EOSStrainTransform:
        """Return finite strain, normalized pressure, and first derivatives.

        Parameters
        ----------
        eos : str or EOSModel
            Pressure EOS family and order.  The order does not change the
            diagnostic transformation.
        pressure, coordinate : array-like
            Broadcast-compatible pressure and positive volume-like coordinate.
            Linear EOS data must be supplied as cubed lengths.
        reference_coordinate : float
            Positive reference volume-like coordinate.  For a linear EOS this
            is ``L0**3``.
        singular_tolerance : float or None, optional
            Absolute strain magnitude below which normalized pressure is marked
            undefined.  The default is scale-independent and based on
            ``sqrt(eps)``.

        Returns
        -------
        EOSStrainTransform
            Transformation values and derivatives.
        """
        model = parse_eos_model(eos)
        family = self.strain_family(model)
        p, x = np.broadcast_arrays(
            np.asarray(pressure, dtype=np.float64),
            np.asarray(coordinate, dtype=np.float64),
        )
        if not np.all(np.isfinite(p)) or not np.all(np.isfinite(x)):
            raise ValueError("pressure and EOS coordinates must be finite")
        if np.any(x <= 0.0):
            raise ValueError("EOS diagnostic coordinates must be positive")
        x0 = float(reference_coordinate)
        if not np.isfinite(x0) or x0 <= 0.0:
            raise ValueError("reference EOS coordinate must be finite and positive")
        tolerance = (
            float(np.sqrt(np.finfo(np.float64).eps))
            if singular_tolerance is None
            else float(singular_tolerance)
        )
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("singular_tolerance must be finite and non-negative")

        if family is EOSStrainFamily.EULERIAN:
            result = self._eulerian(p, x, x0, tolerance)
        elif family is EOSStrainFamily.NATURAL:
            result = self._natural(p, x, x0, tolerance)
        else:
            result = self._vinet(p, x, x0, tolerance)
        metadata = {
            **result.metadata,
            "eos_tag": model.tag,
            "eos_family": model.family.value,
            "reference_coordinate": x0,
            "singular_tolerance": tolerance,
        }
        return EOSStrainTransform(
            family=result.family,
            strain=result.strain,
            normalized_pressure=result.normalized_pressure,
            valid=result.valid,
            dstrain_dcoordinate=result.dstrain_dcoordinate,
            dstrain_dreference=result.dstrain_dreference,
            dnormalized_dpressure=result.dnormalized_dpressure,
            dnormalized_dcoordinate=result.dnormalized_dcoordinate,
            dnormalized_dreference=result.dnormalized_dreference,
            metadata=metadata,
        )

    @staticmethod
    def _empty_outputs(shape: tuple[int, ...]) -> tuple[np.ndarray, ...]:
        nan = np.full(shape, np.nan, dtype=np.float64)
        return tuple(nan.copy() for _ in range(4))

    def _eulerian(
        self, p: np.ndarray, x: np.ndarray, x0: float, tolerance: float
    ) -> EOSStrainTransform:
        ratio23 = (x0 / x) ** (2.0 / 3.0)
        strain = 0.5 * (ratio23 - 1.0)
        valid = np.abs(strain) > tolerance
        normalized, d_p, d_x, d_x0 = self._empty_outputs(strain.shape)
        factor = 1.0 + 2.0 * strain
        denominator = 3.0 * strain * factor**2.5
        normalized[valid] = p[valid] / denominator[valid]
        d_p[valid] = 1.0 / denominator[valid]
        dg_df = 3.0 * factor**1.5 * (1.0 + 7.0 * strain)
        d_f_x = -ratio23 / (3.0 * x)
        d_f_x0 = ratio23 / (3.0 * x0)
        d_f = np.full(strain.shape, np.nan, dtype=np.float64)
        d_f[valid] = -p[valid] * dg_df[valid] / denominator[valid] ** 2
        d_x[valid] = d_f[valid] * d_f_x[valid]
        d_x0[valid] = d_f[valid] * d_f_x0[valid]
        return EOSStrainTransform(
            family=EOSStrainFamily.EULERIAN,
            strain=strain,
            normalized_pressure=normalized,
            valid=valid,
            dstrain_dcoordinate=d_f_x,
            dstrain_dreference=d_f_x0,
            dnormalized_dpressure=d_p,
            dnormalized_dcoordinate=d_x,
            dnormalized_dreference=d_x0,
            metadata={
                "strain_symbol": "f_E",
                "normalized_pressure_symbol": "F_E",
                "definition": "Eulerian finite strain and Birch-Murnaghan normalized pressure",
            },
        )

    def _natural(
        self, p: np.ndarray, x: np.ndarray, x0: float, tolerance: float
    ) -> EOSStrainTransform:
        ratio = x0 / x
        strain = np.log(ratio) / 3.0
        valid = np.abs(strain) > tolerance
        normalized, d_p, d_x, d_x0 = self._empty_outputs(strain.shape)
        denominator = 3.0 * ratio * strain
        normalized[valid] = p[valid] / denominator[valid]
        d_p[valid] = 1.0 / denominator[valid]
        d_f_x = -np.ones_like(x) / (3.0 * x)
        d_f_x0 = np.ones_like(x) / (3.0 * x0)
        dg_x = -ratio * (3.0 * strain + 1.0) / x
        dg_x0 = ratio * (3.0 * strain + 1.0) / x0
        d_x[valid] = -p[valid] * dg_x[valid] / denominator[valid] ** 2
        d_x0[valid] = -p[valid] * dg_x0[valid] / denominator[valid] ** 2
        return EOSStrainTransform(
            family=EOSStrainFamily.NATURAL,
            strain=strain,
            normalized_pressure=normalized,
            valid=valid,
            dstrain_dcoordinate=d_f_x,
            dstrain_dreference=d_f_x0,
            dnormalized_dpressure=d_p,
            dnormalized_dcoordinate=d_x,
            dnormalized_dreference=d_x0,
            metadata={
                "strain_symbol": "f_N",
                "normalized_pressure_symbol": "F_N",
                "definition": "natural finite strain and Poirier-Tarantola normalized pressure",
            },
        )

    def _vinet(
        self, p: np.ndarray, x: np.ndarray, x0: float, tolerance: float
    ) -> EOSStrainTransform:
        root = (x / x0) ** (1.0 / 3.0)
        strain = 1.0 - root
        valid = np.abs(strain) > tolerance
        normalized, d_p, d_x, d_x0 = self._empty_outputs(strain.shape)
        scale = (1.0 - strain) ** 2 / (3.0 * strain)
        normalized[valid] = p[valid] * scale[valid]
        d_p[valid] = scale[valid]
        d_f_x = -root / (3.0 * x)
        d_f_x0 = root / (3.0 * x0)
        dscale_df = (strain**2 - 1.0) / (3.0 * strain**2)
        d_x[valid] = p[valid] * dscale_df[valid] * d_f_x[valid]
        d_x0[valid] = p[valid] * dscale_df[valid] * d_f_x0[valid]
        return EOSStrainTransform(
            family=EOSStrainFamily.VINET,
            strain=strain,
            normalized_pressure=normalized,
            valid=valid,
            dstrain_dcoordinate=d_f_x,
            dstrain_dreference=d_f_x0,
            dnormalized_dpressure=d_p,
            dnormalized_dcoordinate=d_x,
            dnormalized_dreference=d_x0,
            metadata={
                "strain_symbol": "f_V",
                "normalized_pressure_symbol": "F_V",
                "definition": "Vinet strain and normalized pressure",
            },
        )


__all__ = ["EOSStrainFamily", "EOSStrainTransform", "PressureEOSDiagnostics"]
