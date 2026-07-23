# -*- coding: utf-8 -*-

"""Passive data containers for harmonic thermodynamic properties."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(slots=True)
class HarmonicThermodynamicResult:
    """Harmonic thermodynamic properties on a temperature-volume grid.

    Parameters
    ----------
    jobname : str, optional
        Name or short description of the calculation.
    temperature : ndarray or None, optional
        Temperature grid.
    volume : ndarray or None, optional
        Unit-cell volumes.
    static_energy : ndarray or None, optional
        Static energies associated with each volume.
    zero_point_energy : ndarray or None, optional
        Zero-point vibrational energy.
    thermal_energy : ndarray or None, optional
        Thermal contribution to internal energy.
    internal_energy : ndarray or None, optional
        Total internal energy.
    entropy : ndarray or None, optional
        Vibrational entropy.
    vibrational_free_energy : ndarray or None, optional
        Vibrational Helmholtz free energy.
    free_energy : ndarray or None, optional
        Total Helmholtz free energy.
    isochoric_heat_capacity : ndarray or None, optional
        Isochoric heat capacity.
    metadata : dict, optional
        Additional calculation metadata.
    """

    jobname: str = "Unknown"
    temperature: np.ndarray | None = None
    volume: np.ndarray | None = None
    static_energy: np.ndarray | None = None

    zero_point_energy: np.ndarray | None = None
    thermal_energy: np.ndarray | None = None
    internal_energy: np.ndarray | None = None
    entropy: np.ndarray | None = None
    vibrational_free_energy: np.ndarray | None = None
    free_energy: np.ndarray | None = None
    isochoric_heat_capacity: np.ndarray | None = None

    metadata: dict[str, Any] = field(default_factory=dict)

    def as_property_dict(self) -> dict[str, np.ndarray]:
        """Return available thermodynamic arrays using compact property keys.

        Returns
        -------
        dict
            Available properties with keys ``Uzp``, ``Uth``, ``Utot``, ``S``,
            ``Fvib``, ``F``, and ``Cv``.
        """
        mapping = {
            "Uzp": self.zero_point_energy,
            "Uth": self.thermal_energy,
            "Utot": self.internal_energy,
            "S": self.entropy,
            "Fvib": self.vibrational_free_energy,
            "F": self.free_energy,
            "Cv": self.isochoric_heat_capacity,
        }
        return {key: value for key, value in mapping.items() if value is not None}

    def has_thermodynamic_data(self) -> bool:
        """Return whether any thermodynamic property has been calculated.

        Returns
        -------
        bool
            ``True`` if at least one thermodynamic result array is available.
        """
        return bool(self.as_property_dict())
