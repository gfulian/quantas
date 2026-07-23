# -*- coding: utf-8 -*-

"""Frontend-neutral residual and finite-strain diagnostics for EOS records."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from quantas.core.physics.eos import EOSModel, PressureEOSDiagnostics

from .archive import EOSArchive
from .history import EOSFitRecord, EOSResultSlot
from .models import EOSDataset, EOSFitDomain


@dataclass(frozen=True, slots=True)
class EOSDiagnosticResult:
    """Tabular diagnostics derived from one immutable EOS fit record.

    Parameters
    ----------
    record_id : int
        Source fit record identifier.
    slot : EOSResultSlot
        Scientific domain and target.
    columns : dict
        Equal-length diagnostic columns. ``NaN`` is allowed for quantities that
        are undefined, such as normalized pressure at the exact reference
        state or standardized residuals of excluded observations.
    units : dict
        Unit labels associated with columns.
    metadata : dict, optional
        Transformation, residual, and selection semantics.
    warnings : tuple, optional
        Non-fatal diagnostic limitations.
    """

    record_id: int
    slot: EOSResultSlot
    columns: dict[str, np.ndarray]
    units: dict[str, str]
    metadata: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if int(self.record_id) <= 0:
            raise ValueError("EOS diagnostic record_id must be positive")
        object.__setattr__(self, "record_id", int(self.record_id))
        object.__setattr__(self, "slot", EOSResultSlot.parse(self.slot))
        if not self.columns:
            raise ValueError("EOS diagnostic result requires columns")
        normalized: dict[str, np.ndarray] = {}
        size: int | None = None
        for name, values in self.columns.items():
            array = np.asarray(values, dtype=np.float64)
            if array.ndim != 1 or array.size == 0:
                raise ValueError("EOS diagnostic columns must be non-empty vectors")
            if size is None:
                size = int(array.size)
            elif array.size != size:
                raise ValueError("EOS diagnostic columns must have equal length")
            if np.any(np.isinf(array)):
                raise ValueError(f"EOS diagnostic column {name!r} contains infinity")
            normalized[str(name)] = array.copy()
        object.__setattr__(self, "columns", normalized)
        object.__setattr__(
            self,
            "units",
            {str(name): str(unit) for name, unit in self.units.items()},
        )
        object.__setattr__(self, "metadata", dict(self.metadata))
        object.__setattr__(self, "warnings", tuple(str(item) for item in self.warnings))

    @property
    def nrows(self) -> int:
        """Return the number of diagnostic rows."""
        return int(next(iter(self.columns.values())).size)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready representation."""
        return {
            "record_id": self.record_id,
            "slot": self.slot.as_dict(),
            "columns": {name: values.tolist() for name, values in self.columns.items()},
            "units": dict(self.units),
            "metadata": dict(self.metadata),
            "warnings": list(self.warnings),
        }


class EOSDiagnostics:
    """Build scientific diagnostics from successful EOS fit records."""

    def __init__(self, record: EOSFitRecord, dataset: EOSDataset) -> None:
        if not record.result.fit.success or record.result.fit.parameters is None:
            raise ValueError("EOS diagnostics require a successful fit record")
        self.record = record
        self.dataset = dataset
        self._pressure_diagnostics = PressureEOSDiagnostics()

    @classmethod
    def from_archive(
        cls,
        path: str | Path,
        *,
        slot: str | EOSResultSlot | None = None,
        record_id: int | None = None,
    ) -> EOSDiagnostics:
        """Construct diagnostics from an explicit or accepted archive record."""
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

    def build(self, *, include_normalized_pressure: bool = True) -> EOSDiagnosticResult:
        """Return residual diagnostics and optional normalized-pressure data."""
        request = self.record.request
        result = self.record.result
        mask = self.dataset.selection_mask(request.mask)
        npoints = self.dataset.npoints
        group = (
            np.zeros(npoints, dtype=np.float64)
            if self.dataset.groups is None
            else self.dataset.groups.astype(np.float64)
        )
        columns: dict[str, np.ndarray] = {
            "row": np.arange(1, npoints + 1, dtype=np.float64),
            "group": group,
            "included": mask.astype(np.float64),
        }
        units: dict[str, str] = {"row": "1", "group": "1", "included": "1"}
        warnings: list[str] = []

        observed_name, calculated_name, residual_unit = self._residual_columns(
            columns, units
        )
        columns["residual"] = columns[observed_name] - columns[calculated_name]
        units["residual"] = residual_unit
        standardized = np.full(npoints, np.nan, dtype=np.float64)
        fit_diagnostics = result.fit.diagnostics
        if (
            fit_diagnostics is not None
            and fit_diagnostics.standardized_residuals is not None
        ):
            values = np.asarray(
                fit_diagnostics.standardized_residuals, dtype=np.float64
            ).reshape(-1)
            if values.size == int(np.count_nonzero(mask)):
                standardized[mask] = values
            else:
                warnings.append(
                    "standardized residual count does not match the selected observations"
                )
        columns["standardized_residual"] = standardized
        units["standardized_residual"] = "1"

        normalized_metadata: dict[str, Any] = {"available": False}
        if (
            include_normalized_pressure
            and request.domain is EOSFitDomain.PRESSURE_VOLUME
        ):
            try:
                normalized_metadata = self._append_pressure_strain(
                    columns, units, warnings
                )
            except NotImplementedError as exc:
                warnings.append(str(exc))
        metadata = {
            "model": getattr(request.model, "tag", str(request.model)),
            "domain": request.domain.value,
            "target": request.target,
            "residual_definition": result.metadata.get("residual_definition"),
            "standardized_residual_space": result.metadata.get(
                "solver_response",
                result.metadata.get("solver_coordinate", request.target),
            ),
            "selected_observations": int(np.count_nonzero(mask)),
            "excluded_observations": int(np.count_nonzero(~mask)),
            "normalized_pressure": normalized_metadata,
        }
        return EOSDiagnosticResult(
            record_id=self.record.record_id,
            slot=self.record.slot,
            columns=columns,
            units=units,
            metadata=metadata,
            warnings=tuple(dict.fromkeys(warnings)),
        )

    def _residual_columns(
        self,
        columns: dict[str, np.ndarray],
        units: dict[str, str],
    ) -> tuple[str, str, str]:
        request = self.record.request
        predictions = self.record.result.predictions
        if request.domain is EOSFitDomain.PRESSURE_VOLUME:
            observed = self.dataset.column("pressure")
            calculated = predictions["pressure"]
            columns[request.target] = self.dataset.column(request.target)
            columns["observed_pressure"] = observed
            columns["calculated_pressure"] = calculated
            units[request.target] = self.dataset.units.get(request.target, "1")
            units["observed_pressure"] = self.dataset.units.get("pressure", "GPa")
            units["calculated_pressure"] = units["observed_pressure"]
            return (
                "observed_pressure",
                "calculated_pressure",
                units["observed_pressure"],
            )
        if request.domain is EOSFitDomain.VOLUME_TEMPERATURE:
            observed = self.dataset.column(request.target)
            calculated = predictions[request.target]
            columns["temperature"] = self.dataset.column("temperature")
            columns[f"observed_{request.target}"] = observed
            columns[f"calculated_{request.target}"] = calculated
            units["temperature"] = self.dataset.units.get("temperature", "K")
            target_unit = self.dataset.units.get(request.target, "1")
            units[f"observed_{request.target}"] = target_unit
            units[f"calculated_{request.target}"] = target_unit
            return (
                f"observed_{request.target}",
                f"calculated_{request.target}",
                target_unit,
            )
        if request.domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            columns["volume"] = self.dataset.column("volume")
            columns["temperature"] = self.dataset.column("temperature")
            columns["observed_pressure"] = self.dataset.column("pressure")
            columns["calculated_pressure"] = predictions["pressure"]
            units["volume"] = self.dataset.units.get("volume", "angstrom^3")
            units["temperature"] = self.dataset.units.get("temperature", "K")
            pressure_unit = self.dataset.units.get("pressure", "GPa")
            units["observed_pressure"] = pressure_unit
            units["calculated_pressure"] = pressure_unit
            return "observed_pressure", "calculated_pressure", pressure_unit
        raise NotImplementedError(
            f"EOS residual diagnostics do not support {request.domain.value!r}"
        )

    def _append_pressure_strain(
        self,
        columns: dict[str, np.ndarray],
        units: dict[str, str],
        warnings: list[str],
    ) -> dict[str, Any]:
        request = self.record.request
        model = request.model
        assert isinstance(model, EOSModel)
        parameters = self.record.result.parameter_values
        axial = request.target in {"a", "b", "c"}
        coordinate = self.dataset.column(request.target)
        sigma_coordinate = (
            self.dataset.column(f"sigma_{request.target}")
            if self.dataset.has(f"sigma_{request.target}")
            else np.zeros(self.dataset.npoints, dtype=np.float64)
        )
        if axial:
            reference_length = float(parameters["L0"])
            reference_coordinate = reference_length**3
            volume_like = coordinate**3
            sigma_volume_like = 3.0 * coordinate**2 * sigma_coordinate
            sigma_reference = self._parameter_error("L0")
            sigma_reference_coordinate = 3.0 * reference_length**2 * sigma_reference
            modulus_note = (
                "Normalized-pressure values are volume-like; intercepts correspond "
                "to M/3 and slopes to the auxiliary volume EOS."
            )
            warnings.append(modulus_note)
        else:
            reference_coordinate = float(parameters["V0"])
            volume_like = coordinate
            sigma_volume_like = sigma_coordinate
            sigma_reference_coordinate = self._parameter_error("V0")
        pressure = self.dataset.column("pressure")
        calculated_pressure = self.record.result.predictions["pressure"]
        transform = self._pressure_diagnostics.transform(
            model, pressure, volume_like, reference_coordinate
        )
        calculated_transform = self._pressure_diagnostics.transform(
            model, calculated_pressure, volume_like, reference_coordinate
        )
        sigma_pressure = (
            self.dataset.column("sigma_pressure")
            if self.dataset.has("sigma_pressure")
            else np.zeros(self.dataset.npoints, dtype=np.float64)
        )
        sigma_strain = np.sqrt(
            (transform.dstrain_dcoordinate * sigma_volume_like) ** 2
            + (transform.dstrain_dreference * sigma_reference_coordinate) ** 2
        )
        sigma_normalized = np.sqrt(
            (transform.dnormalized_dpressure * sigma_pressure) ** 2
            + (transform.dnormalized_dcoordinate * sigma_volume_like) ** 2
            + (transform.dnormalized_dreference * sigma_reference_coordinate) ** 2
        )
        columns.update(
            {
                "finite_strain": transform.strain,
                "sigma_finite_strain": sigma_strain,
                "normalized_pressure": transform.normalized_pressure,
                "sigma_normalized_pressure": sigma_normalized,
                "calculated_normalized_pressure": calculated_transform.normalized_pressure,
                "normalized_pressure_residual": (
                    transform.normalized_pressure
                    - calculated_transform.normalized_pressure
                ),
                "normalized_pressure_valid": transform.valid.astype(np.float64),
            }
        )
        pressure_unit = self.dataset.units.get("pressure", "GPa")
        units.update(
            {
                "finite_strain": "1",
                "sigma_finite_strain": "1",
                "normalized_pressure": pressure_unit,
                "sigma_normalized_pressure": pressure_unit,
                "calculated_normalized_pressure": pressure_unit,
                "normalized_pressure_residual": pressure_unit,
                "normalized_pressure_valid": "1",
            }
        )
        invalid = int(np.count_nonzero(~transform.valid))
        if invalid:
            warnings.append(
                f"normalized pressure is undefined for {invalid} exact or near-reference observation(s)"
            )
        return {
            "available": True,
            "strain_family": transform.family.value,
            "strain_symbol": transform.metadata["strain_symbol"],
            "normalized_pressure_symbol": transform.metadata[
                "normalized_pressure_symbol"
            ],
            "reference_coordinate": reference_coordinate,
            "reference_coordinate_uncertainty": sigma_reference_coordinate,
            "linear_eos": axial,
            "uncertainty_assumption": (
                "first-order propagation treating measurement errors and the fitted "
                "reference coordinate as independent"
            ),
        }

    def _parameter_error(self, name: str) -> float:
        fit = self.record.result.fit
        if fit.errors is None:
            return 0.0
        try:
            index = fit.parameter_names.index(name)
        except ValueError:
            return 0.0
        value = float(fit.errors[index])
        return value if np.isfinite(value) and value >= 0.0 else 0.0


__all__ = ["EOSDiagnosticResult", "EOSDiagnostics"]
