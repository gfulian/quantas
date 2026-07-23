# -*- coding: utf-8 -*-

r"""Isothermal-to-adiabatic conversion of anisotropic elastic tensors.

The conversion implemented here follows the thermoelastic identity

.. math::

   C^S_{IJ} = C^T_{IJ} + \frac{T V}{C_V}\,\lambda_I\lambda_J,
   \qquad
   \lambda_I = \sum_J C^T_{IJ}\alpha_J,

where ``C^T`` and ``C^S`` are isothermal and adiabatic stiffness matrices,
``V`` is the volume associated with ``C_V``, and ``alpha`` is the thermal-
expansion strain vector in engineering-Voigt convention.  The correction is a
positive-semidefinite rank-one update whenever ``T >= 0``, ``V > 0`` and
``C_V > 0``.

References
----------
Canonical citation keys: ``wallace_1972``, ``davies_1974``, and
``waters_bielawski_2016``.

"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.physics.elasticity.validation import validate_stiffness_matrix
from quantas.references import method_citation_keys


FloatArray: TypeAlias = NDArray[np.float64]
BoolArray: TypeAlias = NDArray[np.bool_]


@dataclass(slots=True)
class AdiabaticStiffnessFieldResult:
    r"""Adiabatic stiffness field and conversion diagnostics.

    Parameters
    ----------
    stiffness : ndarray
        Adiabatic stiffness matrices in GPa with shape ``field_shape + (6, 6)``.
    sigma_stiffness : ndarray or None
        First-order one-standard-deviation uncertainties in GPa.
    correction : ndarray
        Rank-one adiabatic correction ``C^S-C^T`` in GPa.
    thermal_stress : ndarray
        Thermal-stress vector ``lambda = C^T alpha`` in GPa K\ :sup:`-1`.
    valid_mask, invalid_mask : ndarray
        Masks over ``field_shape``.  At exactly zero kelvin the conversion is
        valid and returns ``C^S=C^T`` even when heat capacity or expansion data
        are unavailable.
    metadata : dict
        Formula, units, assumptions, and uncertainty provenance.
    """

    stiffness: FloatArray
    sigma_stiffness: FloatArray | None
    correction: FloatArray
    thermal_stress: FloatArray
    valid_mask: BoolArray
    invalid_mask: BoolArray
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays and validate aligned field dimensions."""
        self.stiffness = np.asarray(self.stiffness, dtype=np.float64).copy()
        self.correction = np.asarray(self.correction, dtype=np.float64).copy()
        self.thermal_stress = np.asarray(self.thermal_stress, dtype=np.float64).copy()
        if self.stiffness.ndim < 2 or self.stiffness.shape[-2:] != (6, 6):
            raise ValueError("stiffness must have shape field_shape + (6, 6)")
        if self.correction.shape != self.stiffness.shape:
            raise ValueError("correction must match stiffness")
        field_shape = self.stiffness.shape[:-2]
        if self.thermal_stress.shape != field_shape + (6,):
            raise ValueError("thermal_stress must have shape field_shape + (6,)")
        self.valid_mask = np.asarray(self.valid_mask, dtype=np.bool_).copy()
        self.invalid_mask = np.asarray(self.invalid_mask, dtype=np.bool_).copy()
        if (
            self.valid_mask.shape != field_shape
            or self.invalid_mask.shape != field_shape
        ):
            raise ValueError("adiabatic masks must match the field shape")
        if np.any(self.valid_mask & self.invalid_mask):
            raise ValueError("valid and invalid masks cannot overlap")
        if self.sigma_stiffness is not None:
            sigma = np.asarray(self.sigma_stiffness, dtype=np.float64)
            if sigma.shape != self.stiffness.shape:
                raise ValueError("sigma_stiffness must match stiffness")
            if np.any(sigma < 0.0):
                raise ValueError("sigma_stiffness must be non-negative")
            self.sigma_stiffness = sigma.copy()
        self.metadata = dict(self.metadata)


def thermal_expansion_tensor_to_voigt(tensor: ArrayLike) -> FloatArray:
    r"""Convert a symmetric Cartesian expansion tensor to engineering Voigt form.

    Parameters
    ----------
    tensor : array-like
        Array with shape ``field_shape + (3, 3)`` in K\ :sup:`-1`.

    Returns
    -------
    ndarray
        Expansion strain vectors ordered ``11, 22, 33, 23, 13, 12``.  Shear
        entries contain ``2 alpha_ij`` because Quantas stiffness matrices use
        engineering shear strains.

    Raises
    ------
    ValueError
        If the tensor shape is invalid or finite entries are not symmetric.
    """
    value = np.asarray(tensor, dtype=np.float64)
    if value.ndim < 2 or value.shape[-2:] != (3, 3):
        raise ValueError("thermal-expansion tensor must end with shape (3, 3)")
    transpose = np.swapaxes(value, -1, -2)
    finite = np.isfinite(value) & np.isfinite(transpose)
    if np.any(finite & ~np.isclose(value, transpose, rtol=1.0e-10, atol=1.0e-14)):
        raise ValueError("thermal-expansion tensor must be symmetric")
    return np.stack(
        (
            value[..., 0, 0],
            value[..., 1, 1],
            value[..., 2, 2],
            2.0 * value[..., 1, 2],
            2.0 * value[..., 0, 2],
            2.0 * value[..., 0, 1],
        ),
        axis=-1,
    ).astype(np.float64, copy=False)


def adiabatic_stiffness_field(
    stiffness_isothermal: ArrayLike,
    temperature: ArrayLike,
    volume_m3: ArrayLike,
    heat_capacity_j_per_k: ArrayLike | None,
    thermal_expansion_tensor: ArrayLike | None,
    *,
    sigma_stiffness_isothermal: ArrayLike | None = None,
    sigma_volume_m3: ArrayLike | None = None,
    sigma_heat_capacity_j_per_k: ArrayLike | None = None,
    sigma_thermal_expansion_tensor: ArrayLike | None = None,
    zero_temperature_tolerance: float = 1.0e-12,
) -> AdiabaticStiffnessFieldResult:
    r"""Convert an isothermal stiffness field to adiabatic conditions.

    Parameters
    ----------
    stiffness_isothermal : array-like
        Isothermal stiffness in GPa with shape ``field_shape + (6, 6)``.
    temperature : array-like
        Absolute temperature in K, broadcastable to ``field_shape``.
    volume_m3 : array-like
        Volume represented by each heat-capacity value, in m\ :sup:`3`.
    heat_capacity_j_per_k : array-like or None
        Isochoric heat capacity for the same cell as ``volume_m3``, in J K\ :sup:`-1`.
    thermal_expansion_tensor : array-like or None
        Cartesian thermal-expansion tensor in K\ :sup:`-1` with shape
        ``field_shape + (3, 3)``.
    sigma_stiffness_isothermal, sigma_volume_m3, sigma_heat_capacity_j_per_k,
    sigma_thermal_expansion_tensor : array-like or None, optional
        One-standard-deviation uncertainties.  Cross-covariances are not
        available in the current QHA archive and are therefore assumed zero.
    zero_temperature_tolerance : float, optional
        Temperatures within this absolute tolerance of zero return the exact
        thermodynamic limit ``C^S=C^T``.

    Returns
    -------
    AdiabaticStiffnessFieldResult
        Converted tensors, uncertainty estimates, validity masks, and
        provenance.

    Notes
    -----
    Nonzero-temperature states are invalid when ``C_V`` or the thermal-
    expansion tensor is absent, non-finite, or non-positive where required.
    Invalid states are represented by NaN tensors; Quantas never substitutes
    the isothermal tensor silently.
    """
    stiffness = np.asarray(stiffness_isothermal, dtype=np.float64)
    if stiffness.ndim < 2 or stiffness.shape[-2:] != (6, 6):
        raise ValueError("stiffness_isothermal must end with shape (6, 6)")
    # Validate each finite matrix without forcing invalid field points to raise.
    finite_matrix = np.all(np.isfinite(stiffness), axis=(-2, -1))
    for matrix in stiffness[finite_matrix]:
        validate_stiffness_matrix(matrix, copy=False)
    field_shape = stiffness.shape[:-2]
    t = np.broadcast_to(np.asarray(temperature, dtype=np.float64), field_shape)
    volume = np.broadcast_to(np.asarray(volume_m3, dtype=np.float64), field_shape)
    if not np.isfinite(zero_temperature_tolerance) or zero_temperature_tolerance < 0.0:
        raise ValueError("zero_temperature_tolerance must be finite and non-negative")
    zero_t = np.isfinite(t) & np.isclose(
        t, 0.0, rtol=0.0, atol=float(zero_temperature_tolerance)
    )

    cv = None
    if heat_capacity_j_per_k is not None:
        cv = np.broadcast_to(
            np.asarray(heat_capacity_j_per_k, dtype=np.float64), field_shape
        )
    alpha_voigt = None
    if thermal_expansion_tensor is not None:
        alpha_tensor = np.broadcast_to(
            np.asarray(thermal_expansion_tensor, dtype=np.float64),
            field_shape + (3, 3),
        )
        alpha_voigt = thermal_expansion_tensor_to_voigt(alpha_tensor)

    valid = (
        finite_matrix
        & np.isfinite(t)
        & (t >= 0.0)
        & np.isfinite(volume)
        & (volume > 0.0)
    )
    if cv is None or alpha_voigt is None:
        valid &= zero_t
    else:
        valid_nonzero = (
            np.isfinite(cv) & (cv > 0.0) & np.all(np.isfinite(alpha_voigt), axis=-1)
        )
        valid &= zero_t | valid_nonzero
    invalid = ~valid

    correction_pa = np.full(field_shape + (6, 6), np.nan, dtype=np.float64)
    thermal_stress_gpa = np.full(field_shape + (6,), np.nan, dtype=np.float64)
    adiabatic = np.full_like(stiffness, np.nan, dtype=np.float64)
    correction_pa[zero_t & valid] = 0.0
    thermal_stress_gpa[zero_t & valid] = 0.0
    adiabatic[zero_t & valid] = stiffness[zero_t & valid]

    active = valid & ~zero_t
    if np.any(active):
        assert cv is not None
        assert alpha_voigt is not None
        c_pa = stiffness * 1.0e9
        beta_pa = np.einsum("...ij,...j->...i", c_pa, alpha_voigt)
        factor = np.zeros(field_shape, dtype=np.float64)
        np.divide(t * volume, cv, out=factor, where=active)
        update = factor[..., None, None] * np.einsum(
            "...i,...j->...ij", beta_pa, beta_pa
        )
        correction_pa[active] = update[active]
        thermal_stress_gpa[active] = beta_pa[active] / 1.0e9
        adiabatic[active] = stiffness[active] + update[active] / 1.0e9

    sigma_adiabatic = _propagate_uncertainty(
        stiffness,
        t,
        volume,
        cv,
        alpha_voigt,
        valid,
        zero_t,
        sigma_stiffness_isothermal,
        sigma_volume_m3,
        sigma_heat_capacity_j_per_k,
        sigma_thermal_expansion_tensor,
    )
    return AdiabaticStiffnessFieldResult(
        stiffness=adiabatic,
        sigma_stiffness=sigma_adiabatic,
        correction=correction_pa / 1.0e9,
        thermal_stress=thermal_stress_gpa,
        valid_mask=valid,
        invalid_mask=invalid,
        metadata={
            "formula": "C_S = C_T + (T V / C_V) outer(C_T alpha, C_T alpha)",
            "stiffness_unit": "GPa",
            "temperature_unit": "K",
            "volume_unit": "m^3 per normalized cell",
            "heat_capacity_unit": "J K^-1 per normalized cell",
            "thermal_expansion_unit": "K^-1",
            "voigt_order": "11 22 33 23 13 12",
            "shear_convention": "engineering strain; alpha4=2 alpha23",
            "zero_temperature_limit": "C_S=C_T",
            "invalid_nonzero_temperature_policy": "NaN; no silent isothermal fallback",
            "uncertainty_method": (
                "first-order independent-input delta method"
                if sigma_adiabatic is not None
                else "none"
            ),
            "uncertainty_cross_covariances": False,
            "citation_keys": list(method_citation_keys("adiabatic_elasticity")),
        },
    )


def _propagate_uncertainty(
    stiffness_gpa: FloatArray,
    temperature: FloatArray,
    volume: FloatArray,
    cv: FloatArray | None,
    alpha_voigt: FloatArray | None,
    valid: BoolArray,
    zero_t: BoolArray,
    sigma_stiffness: ArrayLike | None,
    sigma_volume: ArrayLike | None,
    sigma_cv: ArrayLike | None,
    sigma_alpha_tensor: ArrayLike | None,
) -> FloatArray | None:
    """Return first-order uncertainty under independent input errors."""
    provided = any(
        value is not None
        for value in (sigma_stiffness, sigma_volume, sigma_cv, sigma_alpha_tensor)
    )
    if not provided:
        return None
    field_shape = stiffness_gpa.shape[:-2]
    variance = np.zeros(field_shape + (6, 6), dtype=np.float64)
    variance[~valid] = np.nan

    if sigma_stiffness is not None:
        sigma_c = np.broadcast_to(
            np.asarray(sigma_stiffness, dtype=np.float64), stiffness_gpa.shape
        )
        variance[zero_t & valid] = np.square(sigma_c[zero_t & valid])
    else:
        sigma_c = None
        variance[zero_t & valid] = 0.0

    active = valid & ~zero_t
    if not np.any(active):
        return np.sqrt(variance)
    assert cv is not None
    assert alpha_voigt is not None
    c_pa = stiffness_gpa * 1.0e9
    beta = np.einsum("...ij,...j->...i", c_pa, alpha_voigt)
    factor = np.zeros(field_shape, dtype=np.float64)
    np.divide(temperature * volume, cv, out=factor, where=active)

    if sigma_c is not None:
        sigma_c_pa = sigma_c * 1.0e9
        for p in range(6):
            for q in range(p, 6):
                basis: FloatArray = np.zeros((6, 6), dtype=np.float64)
                basis[p, q] = 1.0
                basis[q, p] = 1.0
                if p == q:
                    basis[p, q] = 1.0
                db = np.einsum("ij,...j->...i", basis, alpha_voigt)
                derivative = basis + factor[..., None, None] * (
                    np.einsum("...i,...j->...ij", db, beta)
                    + np.einsum("...i,...j->...ij", beta, db)
                )
                sigma_component = 0.5 * (sigma_c_pa[..., p, q] + sigma_c_pa[..., q, p])
                contribution = np.square(derivative * sigma_component[..., None, None])
                variance[active] += contribution[active]

    if sigma_alpha_tensor is not None:
        sigma_alpha = thermal_expansion_tensor_to_voigt(
            np.broadcast_to(
                np.asarray(sigma_alpha_tensor, dtype=np.float64),
                field_shape + (3, 3),
            )
        )
        for q in range(6):
            column = c_pa[..., :, q]
            derivative = factor[..., None, None] * (
                np.einsum("...i,...j->...ij", column, beta)
                + np.einsum("...i,...j->...ij", beta, column)
            )
            contribution = np.square(derivative * sigma_alpha[..., q, None, None])
            variance[active] += contribution[active]

    update = factor[..., None, None] * np.einsum("...i,...j->...ij", beta, beta)
    if sigma_volume is not None:
        sigma_v = np.broadcast_to(
            np.asarray(sigma_volume, dtype=np.float64), field_shape
        )
        derivative = update / volume[..., None, None]
        contribution = np.square(derivative * sigma_v[..., None, None])
        variance[active] += contribution[active]
    if sigma_cv is not None:
        sigma_c_v = np.broadcast_to(np.asarray(sigma_cv, dtype=np.float64), field_shape)
        derivative = np.zeros_like(update)
        np.divide(
            -update,
            cv[..., None, None],
            out=derivative,
            where=active[..., None, None],
        )
        contribution = np.square(derivative * sigma_c_v[..., None, None])
        variance[active] += contribution[active]

    variance[~active & ~zero_t] = np.nan
    return np.sqrt(np.maximum(variance, 0.0)) / 1.0e9


__all__ = [
    "AdiabaticStiffnessFieldResult",
    "adiabatic_stiffness_field",
    "thermal_expansion_tensor_to_voigt",
]
