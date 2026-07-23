# -*- coding: utf-8 -*-

"""Frontend-neutral evaluation of fitted EOS records and uncertainties."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from collections.abc import Sequence
from typing import Any, Callable

import numpy as np
from scipy.optimize import brentq

from quantas.core.physics.eos import (
    EOSModel,
    PressureEOS,
    PVTCouplingFamily,
    PVTEOS,
    PVTModel,
    TemperatureEOS,
    TemperatureEOSModel,
)

from .archive import EOSArchive
from .domains.pv import axial_to_volume_parameters
from .history import EOSFitRecord, EOSResultSlot
from .models import EOSDataset, EOSFitDomain
from .domains.pvt import PVTEOSFitModel
from .domains.vt import TemperatureEOSFitModel


@dataclass(frozen=True, slots=True)
class EOSCalculationResult:
    """Tabular properties evaluated from one successful EOS fit record.

    Parameters
    ----------
    record_id : int
        Source immutable fit record.
    slot : EOSResultSlot
        Scientific domain and target represented by the record.
    columns : dict
        Equal-length ``float64`` property columns.
    units : dict
        Unit label for each column when applicable.
    uncertainties : dict, optional
        One-sigma parameter-covariance propagation keyed by property name.
    metadata : dict, optional
        Model, source, sampling, and propagation details.
    warnings : tuple, optional
        Non-fatal extrapolation or propagation warnings.
    """

    record_id: int
    slot: EOSResultSlot
    columns: dict[str, np.ndarray]
    units: dict[str, str]
    uncertainties: dict[str, np.ndarray] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if int(self.record_id) <= 0:
            raise ValueError("EOS calculation record_id must be positive")
        object.__setattr__(self, "record_id", int(self.record_id))
        object.__setattr__(self, "slot", EOSResultSlot.parse(self.slot))
        if not self.columns:
            raise ValueError("EOS calculation result requires property columns")
        normalized: dict[str, np.ndarray] = {}
        size: int | None = None
        for name, values in self.columns.items():
            array = np.asarray(values, dtype=np.float64)
            if array.ndim != 1 or array.size == 0:
                raise ValueError("EOS calculation columns must be non-empty vectors")
            if size is None:
                size = int(array.size)
            elif array.size != size:
                raise ValueError("EOS calculation columns must have equal length")
            if not np.all(np.isfinite(array)):
                raise ValueError(f"EOS calculation column {name!r} is not finite")
            normalized[str(name)] = array.copy()
        uncertainty_values: dict[str, np.ndarray] = {}
        assert size is not None
        for name, values in self.uncertainties.items():
            array = np.asarray(values, dtype=np.float64)
            if array.shape != (size,) or not np.all(np.isfinite(array)):
                raise ValueError(
                    "EOS calculation uncertainties must match the row count"
                )
            if np.any(array < 0.0):
                raise ValueError("EOS calculation uncertainties cannot be negative")
            uncertainty_values[str(name)] = array.copy()
        object.__setattr__(self, "columns", normalized)
        object.__setattr__(
            self,
            "units",
            {str(name): str(value) for name, value in self.units.items()},
        )
        object.__setattr__(self, "uncertainties", uncertainty_values)
        object.__setattr__(self, "metadata", dict(self.metadata))
        object.__setattr__(self, "warnings", tuple(str(item) for item in self.warnings))

    @property
    def nrows(self) -> int:
        """Return the number of calculated states."""
        return int(next(iter(self.columns.values())).size)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready representation."""
        return {
            "record_id": self.record_id,
            "slot": self.slot.as_dict(),
            "columns": {name: value.tolist() for name, value in self.columns.items()},
            "units": dict(self.units),
            "uncertainties": {
                name: value.tolist() for name, value in self.uncertainties.items()
            },
            "metadata": dict(self.metadata),
            "warnings": list(self.warnings),
        }


class EOSCalculator:
    """Evaluate properties from one successful archived or in-memory EOS fit.

    The calculator is independent of Click, Rich, Matplotlib, and HDF5.  The
    :meth:`from_archive` constructor is only a convenience that resolves one
    immutable record and closes the archive before calculations begin.
    """

    def __init__(self, record: EOSFitRecord, dataset: EOSDataset) -> None:
        if not record.result.fit.success or record.result.fit.parameters is None:
            raise ValueError("EOS calculator requires a successful fit record")
        self.record = record
        self.dataset = dataset
        self._pressure = PressureEOS()
        self._temperature = TemperatureEOS()
        self._pvt = PVTEOS()
        self._parameters = np.asarray(record.result.fit.parameters, dtype=np.float64)
        self._parameter_names = tuple(record.result.fit.parameter_names)
        self._covariance = (
            None
            if record.result.fit.covariance is None
            else np.asarray(record.result.fit.covariance, dtype=np.float64)
        )
        if self._parameters.size != len(self._parameter_names):
            raise ValueError("EOS fit record has inconsistent parameter metadata")
        if self._covariance is not None and self._covariance.shape != (
            self._parameters.size,
            self._parameters.size,
        ):
            raise ValueError("EOS fit covariance does not match reported parameters")

    @classmethod
    def from_archive(
        cls,
        path: str | Path,
        *,
        slot: str | EOSResultSlot | None = None,
        record_id: int | None = None,
    ) -> EOSCalculator:
        """Construct a calculator from an explicit or accepted archive record.

        If ``record_id`` is omitted, ``slot`` identifies the current accepted
        result.  If both are omitted, the archive must contain exactly one
        accepted result.
        """
        with EOSArchive(path) as archive:
            if record_id is not None:
                record = archive.record(int(record_id))
                if slot is not None and record.slot != EOSResultSlot.parse(slot):
                    raise ValueError("record_id does not belong to the requested slot")
            else:
                if slot is None:
                    accepted = [
                        state
                        for state in archive.slots()
                        if state.accepted_record_id is not None
                    ]
                    if len(accepted) != 1:
                        keys = ", ".join(state.slot.key for state in accepted) or "none"
                        raise ValueError(
                            "archive must contain exactly one accepted result when "
                            f"--slot is omitted; accepted slots: {keys}"
                        )
                    slot = accepted[0].slot
                accepted_record = archive.accepted(EOSResultSlot.parse(slot))
                if accepted_record is None:
                    raise ValueError(
                        f"EOS result slot {EOSResultSlot.parse(slot).key!r} has no accepted record"
                    )
                record = accepted_record
            dataset = archive.dataset(record.dataset_id)
        return cls(record, dataset)

    @property
    def slot(self) -> EOSResultSlot:
        """Return the scientific result slot evaluated by the calculator."""
        return self.record.slot

    @property
    def parameter_values(self) -> dict[str, float]:
        """Return complete fitted parameters keyed by stable name."""
        return dict(zip(self._parameter_names, self._parameters, strict=True))

    def calculate(
        self,
        *,
        pressure: np.ndarray | Sequence[float] | float | None = None,
        volume: np.ndarray | Sequence[float] | float | None = None,
        temperature: np.ndarray | Sequence[float] | float | None = None,
        propagate_uncertainty: bool = True,
        relative_step: float = 1.0e-5,
    ) -> EOSCalculationResult:
        """Evaluate states for the fitted scientific domain.

        P--V records require exactly one of ``pressure`` or ``volume``.  For an
        axial target, ``volume`` represents the physical fitted length.  V--T
        records require ``temperature`` only.  P--V--T records require
        ``temperature`` and exactly one of ``pressure`` or ``volume``.
        """
        pressure_values = (
            None if pressure is None else np.asarray(pressure, dtype=np.float64)
        )
        volume_values = None if volume is None else np.asarray(volume, dtype=np.float64)
        temperature_values = (
            None if temperature is None else np.asarray(temperature, dtype=np.float64)
        )

        domain = self.record.request.domain
        if domain is EOSFitDomain.PRESSURE_VOLUME:
            evaluator = self._calculate_pv
        elif domain is EOSFitDomain.VOLUME_TEMPERATURE:
            evaluator = self._calculate_vt
        elif domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            evaluator = self._calculate_pvt
        else:
            raise NotImplementedError(
                f"EOS calculator does not support domain {domain.value!r}"
            )
        columns, units, metadata, warnings = evaluator(
            self._parameters,
            pressure=pressure_values,
            volume=volume_values,
            temperature=temperature_values,
        )
        uncertainties: dict[str, np.ndarray] = {}
        if propagate_uncertainty:
            uncertainties, propagation_warnings = self._propagate_uncertainty(
                evaluator,
                columns,
                pressure=pressure_values,
                volume=volume_values,
                temperature=temperature_values,
                relative_step=relative_step,
            )
            warnings.extend(propagation_warnings)
        metadata = {
            **metadata,
            "parameter_order": list(self._parameter_names),
            "uncertainty_method": (
                "parameter-covariance-delta" if uncertainties else "none"
            ),
            "relative_parameter_step": float(relative_step),
            "source_dataset_id": self.record.dataset_id,
            "source_request_id": self.record.request.request_id,
        }
        return EOSCalculationResult(
            record_id=self.record.record_id,
            slot=self.record.slot,
            columns=columns,
            units=units,
            uncertainties=uncertainties,
            metadata=metadata,
            warnings=tuple(dict.fromkeys(warnings)),
        )

    def _calculate_pv(
        self,
        parameters: np.ndarray,
        *,
        pressure: np.ndarray | float | None,
        volume: np.ndarray | float | None,
        temperature: np.ndarray | float | None,
    ) -> tuple[dict[str, np.ndarray], dict[str, str], dict[str, Any], list[str]]:
        if temperature is not None:
            raise ValueError("temperature is not used by a P-V EOS record")
        if (pressure is None) == (volume is None):
            raise ValueError(
                "P-V calculation requires exactly one of pressure or volume"
            )
        request = self.record.request
        model = request.model
        assert isinstance(model, EOSModel)
        axial = request.target in {"a", "b", "c"}
        values = dict(zip(self._parameter_names, parameters, strict=True))
        physical = axial_to_volume_parameters(values) if axial else values
        if pressure is not None:
            p = self._vector(pressure, "pressure")
            q = np.asarray(
                [self._volume_at_pressure(model, physical, float(item)) for item in p],
                dtype=np.float64,
            )
        else:
            coordinate = self._positive_vector(volume, request.target)
            q = coordinate**3 if axial else coordinate
            p = self._pressure.pressure(model, physical, q)
        k = self._pressure.bulk_modulus(model, physical, q)
        kp = self._pressure.bulk_modulus_derivative(model, physical, q)
        kpp = self._pressure.bulk_modulus_second_derivative(model, physical, q)
        coordinate = np.cbrt(q) if axial else q
        modulus = 3.0 * k if axial else k
        derivative = 3.0 * kp if axial else kp
        second = 3.0 * kpp if axial else kpp
        coordinate_name = request.target
        modulus_name = "linear_modulus" if axial else "bulk_modulus"
        derivative_name = (
            "linear_modulus_derivative" if axial else "bulk_modulus_derivative"
        )
        second_name = (
            "linear_modulus_second_derivative"
            if axial
            else "bulk_modulus_second_derivative"
        )
        extrapolated = self._outside_sampled_range(coordinate, coordinate_name)
        columns = {
            "pressure": np.asarray(p, dtype=np.float64),
            coordinate_name: np.asarray(coordinate, dtype=np.float64),
            modulus_name: np.asarray(modulus, dtype=np.float64),
            derivative_name: np.asarray(derivative, dtype=np.float64),
            second_name: np.asarray(second, dtype=np.float64),
            "extrapolated": extrapolated.astype(np.float64),
        }
        units = {
            "pressure": self.dataset.units.get("pressure", "GPa"),
            coordinate_name: self.dataset.units.get(
                coordinate_name, "angstrom" if axial else "angstrom^3"
            ),
            modulus_name: self.dataset.units.get("pressure", "GPa"),
            derivative_name: "1",
            second_name: f"{self.dataset.units.get('pressure', 'GPa')}^-1",
            "extrapolated": "1",
        }
        warnings = self._extrapolation_warnings(extrapolated)
        return (
            columns,
            units,
            {
                "model": model.as_dict(),
                "input_mode": "pressure" if pressure is not None else coordinate_name,
                "linear_eos": axial,
            },
            warnings,
        )

    def _calculate_vt(
        self,
        parameters: np.ndarray,
        *,
        pressure: np.ndarray | float | None,
        volume: np.ndarray | float | None,
        temperature: np.ndarray | float | None,
    ) -> tuple[dict[str, np.ndarray], dict[str, str], dict[str, Any], list[str]]:
        if pressure is not None or volume is not None:
            raise ValueError("V-T calculation requires temperature only")
        if temperature is None:
            raise ValueError("V-T calculation requires temperature")
        temp = self._nonnegative_vector(temperature, "temperature")
        request = self.record.request
        model = request.model
        assert isinstance(model, TemperatureEOSModel)
        axial = request.target in {"a", "b", "c"}
        adapter = TemperatureEOSFitModel(model, request.target, axial=axial)
        auxiliary = adapter.evaluate(temp, parameters)
        alpha_aux = adapter.expansion_coefficient(temp, parameters)
        target = np.cbrt(auxiliary) if axial else auxiliary
        alpha = alpha_aux / 3.0 if axial else alpha_aux
        derivative = alpha * target
        extrapolated = self._outside_temperature_range(temp)
        columns = {
            "temperature": temp,
            request.target: target,
            "expansion_coefficient": alpha,
            "temperature_derivative": derivative,
            "extrapolated": extrapolated.astype(np.float64),
        }
        units = {
            "temperature": self.dataset.units.get("temperature", "K"),
            request.target: self.dataset.units.get(request.target, "1"),
            "expansion_coefficient": "K^-1",
            "temperature_derivative": f"{self.dataset.units.get(request.target, '1')} K^-1",
            "extrapolated": "1",
        }
        return (
            columns,
            units,
            {
                "model": {
                    "family": model.family.value,
                    "variant": None if model.variant is None else model.variant.value,
                    "tag": model.tag,
                },
                "input_mode": "temperature",
                "linear_eos": axial,
            },
            self._extrapolation_warnings(extrapolated),
        )

    def _calculate_pvt(
        self,
        parameters: np.ndarray,
        *,
        pressure: np.ndarray | float | None,
        volume: np.ndarray | float | None,
        temperature: np.ndarray | float | None,
    ) -> tuple[dict[str, np.ndarray], dict[str, str], dict[str, Any], list[str]]:
        if temperature is None:
            raise ValueError("P-V-T calculation requires temperature")
        if (pressure is None) == (volume is None):
            raise ValueError(
                "P-V-T calculation requires exactly one of pressure or volume"
            )
        temp = self._nonnegative_vector(temperature, "temperature")
        model = self.record.request.model
        assert isinstance(model, PVTModel)
        adapter = PVTEOSFitModel(model)
        pressure_parameters, thermal_parameters, coupling_parameters = (
            adapter.split_parameters(parameters)
        )
        if pressure is not None:
            p, temp = np.broadcast_arrays(self._vector(pressure, "pressure"), temp)
            p = p.reshape(-1)
            temp = temp.reshape(-1)
            vol = self._pvt.volume(
                model,
                pressure_parameters,
                thermal_parameters,
                coupling_parameters,
                p,
                temp,
            )
            input_mode = "pressure-temperature"
        else:
            vol, temp = np.broadcast_arrays(
                self._positive_vector(volume, "volume"), temp
            )
            vol = vol.reshape(-1)
            temp = temp.reshape(-1)
            p = self._pvt.pressure(
                model,
                pressure_parameters,
                thermal_parameters,
                coupling_parameters,
                vol,
                temp,
            )
            input_mode = "volume-temperature"
        bulk = self._pvt.bulk_modulus(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            vol,
            temp,
        )
        kp, kpp = self._pvt_bulk_pressure_derivatives(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            vol,
            temp,
        )
        alpha = self._pvt_expansion_coefficient(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            vol,
            temp,
            bulk,
        )
        dkdt = self._pvt_bulk_temperature_derivative_at_pressure(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            p,
            temp,
        )
        reference_volume = self._pvt.reference_volume(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            temp,
        )
        k0_t = self._pvt.zero_pressure_bulk_modulus(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            temp,
        )
        alpha0_t = self._pvt.expansion_coefficient_zero_pressure(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            temp,
        )
        dk0dt = self._pvt.zero_pressure_bulk_modulus_temperature_derivative(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            temp,
        )
        thermal_pressure: np.ndarray | None = None
        if model.coupling_family is PVTCouplingFamily.THERMAL_PRESSURE:
            thermal_pressure = self._pvt.thermal_pressure_contribution(
                model,
                pressure_parameters,
                coupling_parameters,
                vol,
                temp,
            )
        extrapolated = self._outside_pvt_range(p, vol, temp)
        columns = {
            "pressure": np.asarray(p, dtype=np.float64),
            "temperature": np.asarray(temp, dtype=np.float64),
            "volume": np.asarray(vol, dtype=np.float64),
            "bulk_modulus": np.asarray(bulk, dtype=np.float64),
            "bulk_modulus_derivative": np.asarray(kp, dtype=np.float64),
            "bulk_modulus_second_derivative": np.asarray(kpp, dtype=np.float64),
            "expansion_coefficient": np.asarray(alpha, dtype=np.float64),
            "dK_dT_at_pressure": np.asarray(dkdt, dtype=np.float64),
            "reference_volume": np.asarray(reference_volume, dtype=np.float64),
            "zero_pressure_bulk_modulus": np.asarray(k0_t, dtype=np.float64),
            "zero_pressure_expansion_coefficient": np.asarray(
                alpha0_t, dtype=np.float64
            ),
            "zero_pressure_dK0_dT": np.asarray(dk0dt, dtype=np.float64),
            "extrapolated": extrapolated.astype(np.float64),
        }
        if thermal_pressure is not None:
            columns["thermal_pressure"] = np.asarray(
                thermal_pressure,
                dtype=np.float64,
            )
        pressure_unit = self.dataset.units.get("pressure", "GPa")
        units = {
            "pressure": pressure_unit,
            "temperature": self.dataset.units.get("temperature", "K"),
            "volume": self.dataset.units.get("volume", "angstrom^3"),
            "bulk_modulus": pressure_unit,
            "bulk_modulus_derivative": "1",
            "bulk_modulus_second_derivative": f"{pressure_unit}^-1",
            "expansion_coefficient": "K^-1",
            "dK_dT_at_pressure": f"{pressure_unit} K^-1",
            "reference_volume": self.dataset.units.get("volume", "angstrom^3"),
            "zero_pressure_bulk_modulus": pressure_unit,
            "zero_pressure_expansion_coefficient": "K^-1",
            "zero_pressure_dK0_dT": f"{pressure_unit} K^-1",
            "extrapolated": "1",
        }
        if thermal_pressure is not None:
            units["thermal_pressure"] = pressure_unit
        return (
            columns,
            units,
            {
                "model": model.as_dict(),
                "input_mode": input_mode,
            },
            self._extrapolation_warnings(extrapolated),
        )

    def _propagate_uncertainty(
        self,
        evaluator: Callable[
            ..., tuple[dict[str, np.ndarray], dict[str, str], dict[str, Any], list[str]]
        ],
        base_columns: dict[str, np.ndarray],
        *,
        pressure: np.ndarray | float | None,
        volume: np.ndarray | float | None,
        temperature: np.ndarray | float | None,
        relative_step: float,
    ) -> tuple[dict[str, np.ndarray], list[str]]:
        if self._covariance is None:
            return {}, [
                "parameter covariance is unavailable; uncertainties were not propagated"
            ]
        if not np.isfinite(relative_step) or relative_step <= 0.0:
            raise ValueError("relative_step must be finite and positive")
        independent_columns = {"extrapolated"}
        if pressure is not None:
            independent_columns.add("pressure")
        if volume is not None:
            independent_columns.add(self.record.request.target)
        if temperature is not None:
            independent_columns.add("temperature")
        property_names = [
            name for name in base_columns if name not in independent_columns
        ]
        nrows = next(iter(base_columns.values())).size
        jacobians = {
            name: np.zeros((nrows, self._parameters.size), dtype=np.float64)
            for name in property_names
        }
        warnings: list[str] = []
        base = self._parameters
        variances = np.clip(np.diag(self._covariance), 0.0, None)
        for index, variance in enumerate(variances):
            if variance == 0.0:
                continue
            scale = max(abs(float(base[index])), 1.0)
            step = max(relative_step * scale, np.sqrt(np.finfo(np.float64).eps) * scale)
            plus = base.copy()
            minus = base.copy()
            plus[index] += step
            minus[index] -= step
            plus_values = minus_values = None
            try:
                plus_values = evaluator(
                    plus,
                    pressure=pressure,
                    volume=volume,
                    temperature=temperature,
                )[0]
            except (ValueError, FloatingPointError, OverflowError):
                pass
            try:
                minus_values = evaluator(
                    minus,
                    pressure=pressure,
                    volume=volume,
                    temperature=temperature,
                )[0]
            except (ValueError, FloatingPointError, OverflowError):
                pass
            if plus_values is None and minus_values is None:
                warnings.append(
                    f"uncertainty derivative for parameter {self._parameter_names[index]} could not be evaluated"
                )
                continue
            for name in property_names:
                if plus_values is not None and minus_values is not None:
                    derivative = (plus_values[name] - minus_values[name]) / (2.0 * step)
                elif plus_values is not None:
                    derivative = (plus_values[name] - base_columns[name]) / step
                else:
                    assert minus_values is not None
                    derivative = (base_columns[name] - minus_values[name]) / step
                jacobians[name][:, index] = derivative
        uncertainties: dict[str, np.ndarray] = {}
        for name, jacobian in jacobians.items():
            variance = np.einsum("ij,jk,ik->i", jacobian, self._covariance, jacobian)
            uncertainties[name] = np.sqrt(np.clip(variance, 0.0, None))
        return uncertainties, warnings

    def _pvt_bulk_pressure_derivatives(
        self,
        model: PVTModel,
        pressure_parameters: dict[str, float],
        thermal_parameters: dict[str, float] | None,
        coupling_parameters: dict[str, float],
        volume: np.ndarray,
        temperature: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        step = np.cbrt(np.finfo(np.float64).eps) * np.maximum(volume, 1.0)
        lower = np.maximum(volume - step, np.nextafter(0.0, 1.0))
        upper = volume + step
        k_upper = self._pvt.bulk_modulus(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            upper,
            temperature,
        )
        k_lower = self._pvt.bulk_modulus(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            lower,
            temperature,
        )
        p_upper = self._pvt.pressure(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            upper,
            temperature,
        )
        p_lower = self._pvt.pressure(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            lower,
            temperature,
        )
        kp = (k_upper - k_lower) / (p_upper - p_lower)
        kp_upper: list[float] = []
        kp_lower: list[float] = []
        for v, t, h in zip(volume, temperature, step, strict=True):
            for local_v, target in (
                (v + h, kp_upper),
                (max(v - h, np.nextafter(0.0, 1.0)), kp_lower),
            ):
                inner = max(
                    np.cbrt(np.finfo(np.float64).eps) * max(local_v, 1.0), 1.0e-8
                )
                va = max(local_v - inner, np.nextafter(0.0, 1.0))
                vb = local_v + inner
                ka = self._pvt.bulk_modulus(
                    model,
                    pressure_parameters,
                    thermal_parameters,
                    coupling_parameters,
                    [va],
                    [t],
                )[0]
                kb = self._pvt.bulk_modulus(
                    model,
                    pressure_parameters,
                    thermal_parameters,
                    coupling_parameters,
                    [vb],
                    [t],
                )[0]
                pa = self._pvt.pressure(
                    model,
                    pressure_parameters,
                    thermal_parameters,
                    coupling_parameters,
                    [va],
                    [t],
                )[0]
                pb = self._pvt.pressure(
                    model,
                    pressure_parameters,
                    thermal_parameters,
                    coupling_parameters,
                    [vb],
                    [t],
                )[0]
                target.append((kb - ka) / (pb - pa))
        kpp = (np.asarray(kp_upper) - np.asarray(kp_lower)) / (p_upper - p_lower)
        return np.asarray(kp), np.asarray(kpp)

    def _pvt_expansion_coefficient(
        self,
        model: PVTModel,
        pressure_parameters: dict[str, float],
        thermal_parameters: dict[str, float] | None,
        coupling_parameters: dict[str, float],
        volume: np.ndarray,
        temperature: np.ndarray,
        bulk_modulus: np.ndarray,
    ) -> np.ndarray:
        step = np.cbrt(np.finfo(np.float64).eps) * np.maximum(temperature, 1.0)
        lower = np.maximum(temperature - step, 0.0)
        upper = temperature + step
        p_upper = self._pvt.pressure(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            volume,
            upper,
        )
        p_lower = self._pvt.pressure(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            volume,
            lower,
        )
        d_p_d_t = (p_upper - p_lower) / (upper - lower)
        return d_p_d_t / bulk_modulus

    def _pvt_bulk_temperature_derivative_at_pressure(
        self,
        model: PVTModel,
        pressure_parameters: dict[str, float],
        thermal_parameters: dict[str, float] | None,
        coupling_parameters: dict[str, float],
        pressure: np.ndarray,
        temperature: np.ndarray,
    ) -> np.ndarray:
        step = np.cbrt(np.finfo(np.float64).eps) * np.maximum(temperature, 1.0)
        lower = np.maximum(temperature - step, 0.0)
        upper = temperature + step
        v_upper = self._pvt.volume(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            pressure,
            upper,
        )
        v_lower = self._pvt.volume(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            pressure,
            lower,
        )
        k_upper = self._pvt.bulk_modulus(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            v_upper,
            upper,
        )
        k_lower = self._pvt.bulk_modulus(
            model,
            pressure_parameters,
            thermal_parameters,
            coupling_parameters,
            v_lower,
            lower,
        )
        return (k_upper - k_lower) / (upper - lower)

    def _volume_at_pressure(
        self,
        model: EOSModel,
        parameters: dict[str, float],
        pressure: float,
    ) -> float:
        target = float(pressure)
        centre = float(parameters["V0"])

        def residual(value: float) -> float:
            return float(self._pressure.pressure(model, parameters, value)) - target

        f_centre = residual(centre)
        if f_centre == 0.0:
            return centre
        lower = upper = centre
        f_lower = f_upper = f_centre
        for _ in range(96):
            candidate = max(lower / 1.25, np.nextafter(0.0, 1.0))
            try:
                f_candidate = residual(candidate)
            except ValueError:
                f_candidate = f_lower
            if f_candidate * f_lower < 0.0:
                return float(
                    brentq(residual, candidate, lower, xtol=1.0e-12, rtol=1.0e-12)
                )
            lower, f_lower = candidate, f_candidate
            candidate = upper * 1.25
            try:
                f_candidate = residual(candidate)
            except ValueError:
                f_candidate = f_upper
            if f_upper * f_candidate < 0.0:
                return float(
                    brentq(residual, upper, candidate, xtol=1.0e-12, rtol=1.0e-12)
                )
            upper, f_upper = candidate, f_candidate
        raise ValueError(f"could not bracket EOS volume at pressure {pressure:g}")

    def _outside_sampled_range(self, values: np.ndarray, target: str) -> np.ndarray:
        selected = self.record.request.mask
        mask = self.dataset.selection_mask(selected)
        sample = self.dataset.column(target)[mask]
        return np.logical_or(values < np.min(sample), values > np.max(sample))

    def _outside_temperature_range(self, values: np.ndarray) -> np.ndarray:
        selected = self.dataset.selection_mask(self.record.request.mask)
        sample = self.dataset.column("temperature")[selected]
        return np.logical_or(values < np.min(sample), values > np.max(sample))

    def _outside_pvt_range(
        self, pressure: np.ndarray, volume: np.ndarray, temperature: np.ndarray
    ) -> np.ndarray:
        selected = self.dataset.selection_mask(self.record.request.mask)
        outside = np.zeros(pressure.shape, dtype=np.bool_)
        for values, name in (
            (pressure, "pressure"),
            (volume, "volume"),
            (temperature, "temperature"),
        ):
            sample = self.dataset.column(name)[selected]
            outside |= np.logical_or(values < np.min(sample), values > np.max(sample))
        return outside

    @staticmethod
    def _extrapolation_warnings(extrapolated: np.ndarray) -> list[str]:
        count = int(np.count_nonzero(extrapolated))
        return (
            []
            if count == 0
            else [f"{count} calculated state(s) lie outside the sampled fit ranges"]
        )

    @staticmethod
    def _vector(value: np.ndarray | float, name: str) -> np.ndarray:
        array = np.asarray(value, dtype=np.float64)
        if array.ndim == 0:
            array = array.reshape(1)
        if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
            raise ValueError(f"{name} must be a non-empty finite vector")
        return array.copy()

    @classmethod
    def _positive_vector(
        cls, value: np.ndarray | float | None, name: str
    ) -> np.ndarray:
        if value is None:
            raise ValueError(f"{name} is required")
        array = cls._vector(value, name)
        if np.any(array <= 0.0):
            raise ValueError(f"{name} must be positive")
        return array

    @classmethod
    def _nonnegative_vector(cls, value: np.ndarray | float, name: str) -> np.ndarray:
        array = cls._vector(value, name)
        if np.any(array < 0.0):
            raise ValueError(f"{name} cannot be negative")
        return array


__all__ = ["EOSCalculationResult", "EOSCalculator"]
