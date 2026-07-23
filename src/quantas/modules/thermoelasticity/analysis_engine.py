# -*- coding: utf-8 -*-

"""Reusable analysis engine for calibrated thermoelastic models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

from quantas.models import ResultData
from quantas.modules.thermoelasticity.models import (
    ThermoelasticDepthProfile,
    ThermoelasticOptions,
    ThermoelasticResult,
)
from quantas.modules.thermoelasticity.postfit import (
    analyze_thermoelastic_profiles,
    evaluate_thermoelastic_grid,
    thermoelastic_options_from_mapping,
)


@dataclass(slots=True)
class ThermoelasticAnalysisEngine:
    """Evaluate point, grid, and depth-profile states through one path.

    Parameters
    ----------
    source : ResultData
        Calibration or prior analysis archive containing the thermoelastic
        payload and archived QHA source fields.
    options : ThermoelasticOptions or None, optional
        Analysis policies.  When omitted, options are reconstructed from the
        archive while ignoring obsolete frontend keys.
    """

    source: ResultData
    options: ThermoelasticOptions

    def __init__(
        self,
        source: ResultData,
        options: ThermoelasticOptions | None = None,
    ) -> None:
        payload = source.results.get("thermoelasticity")
        if not isinstance(payload, ThermoelasticResult):
            raise ValueError("result does not contain a thermoelasticity payload")
        self.source = source
        self.options = (
            thermoelastic_options_from_mapping(source.options)
            if options is None
            else options
        )

    @property
    def payload(self) -> ThermoelasticResult:
        """Return the calibrated thermoelastic payload."""
        value = self.source.results.get("thermoelasticity")
        assert isinstance(value, ThermoelasticResult)
        return value

    def evaluate_point(self, pressure: float, temperature: float) -> ResultData:
        """Evaluate one pressure-temperature state."""
        return self.evaluate_grid([pressure], [temperature])

    def evaluate_grid(
        self,
        pressure: ArrayLike,
        temperature: ArrayLike,
    ) -> ResultData:
        """Evaluate a Cartesian pressure-temperature grid."""
        payload = evaluate_thermoelastic_grid(
            self.payload,
            temperature=np.asarray(temperature, dtype=np.float64),
            pressure=np.asarray(pressure, dtype=np.float64),
            options=self.options,
        )
        return ResultData(
            metadata=self.source.metadata,
            input_data=self.source.input_data,
            options=dict(self.source.options),
            results={**self.source.results, "thermoelasticity": payload},
            warnings=list(self.source.warnings),
            events=list(self.source.events),
        )

    def evaluate_profiles(
        self,
        profiles: Sequence[ThermoelasticDepthProfile],
    ) -> ResultData:
        """Evaluate one or more depth-dependent profiles."""
        return analyze_thermoelastic_profiles(
            self.source,
            profiles,
            options=self.options,
        )

    def constant_pressure_section(
        self,
        pressure: float,
    ) -> ThermoelasticResult:
        """Evaluate all archived temperatures at one exact pressure."""
        result = self.evaluate_grid([pressure], self.payload.temperature)
        value = result.results["thermoelasticity"]
        assert isinstance(value, ThermoelasticResult)
        return value

    def constant_temperature_section(
        self,
        temperature: float,
    ) -> ThermoelasticResult:
        """Evaluate all archived pressures at one exact temperature."""
        result = self.evaluate_grid(self.payload.pressure, [temperature])
        value = result.results["thermoelasticity"]
        assert isinstance(value, ThermoelasticResult)
        return value


__all__ = ["ThermoelasticAnalysisEngine"]
