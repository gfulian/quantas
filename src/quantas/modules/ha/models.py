# -*- coding: utf-8 -*-

"""Data containers used by the harmonic-approximation workflow."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from quantas.models.kieffer import KiefferThermodynamicContribution
from quantas.models.phonons import PhononInputData
from quantas.models.thermodynamics import HarmonicThermodynamicResult


HAKiefferContribution = KiefferThermodynamicContribution


@dataclass(slots=True)
class HAInput(PhononInputData):
    """Input data for a harmonic-approximation calculation."""


@dataclass(slots=True)
class HAOptions:
    """Options controlling a harmonic-approximation calculation.

    Parameters
    ----------
    temperature_min : float, optional
        Minimum temperature of the calculation range.
    temperature_max : float, optional
        Maximum temperature of the calculation range.
    temperature_step : float, optional
        Temperature increment.
    energy_unit : str, optional
        Energy unit used for input static energies and final energy-like
        results.
    volume_unit : str, optional
        Length unit used to express unit-cell volumes.
    frequency_unit : str, optional
        Frequency unit used by the input phonon data.
    temperature_unit : str, optional
        Temperature unit used by the input temperature range.
    calculate_thermodynamics : bool, optional
        If ``True``, calculate harmonic thermodynamic properties.
    """

    temperature_min: float = 298.15
    temperature_max: float = 298.15
    temperature_step: float = 1.0

    energy_unit: str = "Ha"
    volume_unit: str = "A"
    frequency_unit: str = "cm^-1"
    temperature_unit: str = "K"

    calculate_thermodynamics: bool = True

    def temperature_grid(self) -> np.ndarray:
        """Return the temperature grid used by the HA workflow.

        Returns
        -------
        ndarray
            One-dimensional temperature array generated with
            ``np.arange(min, max + step, step)``.

        Raises
        ------
        ValueError
            If the temperature range or step is invalid.
        """
        if self.temperature_step <= 0.0:
            raise ValueError("temperature_step must be positive")
        if self.temperature_max < self.temperature_min:
            raise ValueError(
                "temperature_max must be greater than or equal to temperature_min"
            )
        return np.arange(
            self.temperature_min,
            self.temperature_max + self.temperature_step,
            self.temperature_step,
            dtype=np.float64,
        )


@dataclass(slots=True)
class HAResult(HarmonicThermodynamicResult):
    """Results of a harmonic-approximation calculation."""

    kieffer_contribution: HAKiefferContribution | None = None
