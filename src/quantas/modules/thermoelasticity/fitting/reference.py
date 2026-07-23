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

from quantas.core.physics.eos import (
    EnergyEOS,
    resolve_energy_parameters,
    resolved_energy_parameter_covariance,
)
from quantas.core.physics.units import energy_to_pressure
from quantas.references import method_citation_keys
from quantas.modules.thermoelasticity.models import (
    ReferenceEOSFit,
)


FloatArray: TypeAlias = NDArray[np.float64]


def fit_reference_static_eos(
    volume: ArrayLike,
    energy: ArrayLike,
    *,
    eos: str = "BM3",
    energy_unit: str = "Ha",
    volume_unit: str = "A",
    max_iterations: int | None = None,
) -> ReferenceEOSFit:
    """Fit the authoritative static energy-volume EOS used by QSA.

    Parameters
    ----------
    volume, energy : array_like
        QHA sampled static volume and energy data.
    eos : str, optional
        Integrated energy-EOS tag.
    energy_unit : str, optional
        Static-energy unit.
    volume_unit : str, optional
        Length unit whose cube defines the volume unit.
    max_iterations : int or None, optional
        Maximum EOS model evaluations.

    Returns
    -------
    ReferenceEOSFit
        Fixed physical parameters and complete diagnostics.

    Notes
    -----
    ``volume`` and ``energy`` are the sampled static DFT energy--volume
    arrays archived by the QHA workflow.  Consequently ``V0``, ``K0`` and
    ``Kprime`` describe the cold static reference at 0 K and 0 GPa.  Thermal
    vibrational and zero-point energies are not included in this fit.  QHA
    thermal effects enter later only through the equilibrium volume
    :math:`V(P,T)` and, for adiabatic conversion, through heat capacity and
    thermal expansion.

    References
    ----------
    Canonical citation key: ``stixrude_lithgow_bertelloni_2005``.

    Raises
    ------
    ValueError
        If the EOS fit fails or its covariance is incompatible.
    """
    eos_service = EnergyEOS()
    model = eos_service.model(eos)
    volume_values = np.asarray(volume, dtype=np.float64)
    energy_values = np.asarray(energy, dtype=np.float64)
    fit = eos_service.fit(
        model,
        volume_values,
        energy_values,
        maxfev=max_iterations,
    )
    if not fit.success or fit.parameters is None:
        raise ValueError(fit.message or "static reference EOS fit failed")

    raw_parameters = resolve_energy_parameters(model, fit.parameters)
    pressure_factor = float(energy_to_pressure(1.0, energy_unit, volume_unit, "GPa"))
    k0 = float(raw_parameters.K0 * pressure_factor)
    kp = float(raw_parameters.KP)
    kpp = float(raw_parameters.KPP / pressure_factor)
    v0 = float(raw_parameters.V0)

    covariance = None
    resolved_covariance_raw = None
    if fit.covariance is not None:
        resolved_covariance_raw = resolved_energy_parameter_covariance(
            model,
            np.asarray(fit.parameters, dtype=np.float64),
            np.asarray(fit.covariance, dtype=np.float64),
        )
        transform = np.diag(
            np.asarray([1.0, pressure_factor, 1.0, 1.0 / pressure_factor, 1.0])
        )
        converted = transform @ resolved_covariance_raw @ transform.T
        indexes = np.asarray([4, 1, 2], dtype=np.int64)
        covariance = np.asarray(converted[np.ix_(indexes, indexes)], dtype=np.float64)

    return ReferenceEOSFit(
        eos=model.tag,
        reference_volume=v0,
        bulk_modulus=k0,
        bulk_modulus_derivative=kp,
        bulk_modulus_second_derivative=kpp,
        covariance=covariance,
        fit=fit,
        metadata={
            "energy_unit": energy_unit,
            "volume_unit": volume_unit,
            "pressure_unit": "GPa",
            "pressure_conversion_factor": pressure_factor,
            "resolved_parameter_order": ["V0", "K0", "Kprime"],
            "raw_resolved_covariance_available": resolved_covariance_raw is not None,
            "reference_state": "static 0 K, P=0",
            "source_field": "QHA static_energy sampled versus volume",
            "includes_zero_point_energy": False,
            "includes_thermal_vibrational_energy": False,
            "citation_keys": list(method_citation_keys("cold_finite_strain")),
        },
    )


__all__ = ["fit_reference_static_eos"]
