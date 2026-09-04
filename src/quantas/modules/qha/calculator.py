# -*- coding: utf-8 -*-

"""Calculator for quasi-harmonic approximation workflows.

The :class:`QHACalculator` class connects the QHA numerical workflow to the
common Quantas calculator model. It validates input data, runs the selected
pressure-temperature minimization workflow, converts workflow callbacks into
Quantas events, and returns a generic :class:`~quantas.models.ResultData`
container suitable for the Python API, command-line interface, and graphical
interfaces.
"""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

from quantas.core.events import EventLevel, Observer
from quantas.models import BasicCalculator, InputData, ResultData, ResultMetadata
from quantas.models.kieffer import KiefferVolumeSeries
from quantas.modules.qha.inspect import pressure_volume_preview
from quantas.modules.qha.kieffer import validate_kieffer_qha_applicability
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.qha.structural import calculate_structural_thermal_expansion
from quantas.modules.qha.thermodynamics import (
    FrequencyThermodynamicEvaluator,
    calculate_frequency_thermodynamics_at_equilibrium,
    calculate_mixed_thermal_expansion,
    calculate_mode_gruneisen_at_equilibrium,
    calculate_sampled_thermodynamics,
    free_energy_grid,
    interpolate_thermodynamics_at_equilibrium,
    refresh_thermal_expansion_dependents,
)
from quantas.modules.qha.workflow import (
    QHAWorkflowEvent,
    run_volume_minimization_workflow,
)


class QHACalculator(BasicCalculator):
    """Calculator for quasi-harmonic approximation properties.

    Parameters
    ----------
    qha_input : QHAInput
        Normalized input data used by the QHA workflow.
    options : QHAOptions or None, optional
        Options controlling the calculation. If ``None``, default QHA options
        are used.
    kieffer_cutoffs : KiefferVolumeSeries or None, optional
        Direct acoustic cutoff series for ``scheme="td"``.
    observer : Observer or None, optional
        Observer receiving workflow events. If ``None``, a null observer is
        used by the base calculator.
    """

    module = "qha"
    method = "quasi-harmonic"

    def __init__(
        self,
        qha_input: QHAInput,
        options: QHAOptions | None = None,
        kieffer_cutoffs: KiefferVolumeSeries | None = None,
        observer: Observer | None = None,
    ) -> None:
        """Initialize the QHA calculator."""
        self.qha_input = qha_input
        self.qha_options = options if options is not None else QHAOptions()
        self.kieffer_cutoffs = kieffer_cutoffs

        input_data = InputData(
            source=qha_input.source,
            data={
                "jobname": qha_input.jobname,
                "natoms": qha_input.natoms,
                "formula_units": qha_input.formula_units,
                "natoms_per_formula_unit": qha_input.natoms_per_formula_unit,
                "nvol": qha_input.nvol,
                "qpoints": qha_input.qpoints,
                "nmodes": qha_input.nmodes,
                "total_q_points": qha_input.total_q_points,
                "has_structure": qha_input.has_structure(),
                "structure": (
                    None
                    if qha_input.structure is None
                    else qha_input.structure.as_dict(include_source=True)
                ),
                "has_phonons": qha_input.has_phonons(),
                "has_mode_continuity": qha_input.has_mode_continuity(),
                "mode_continuity": qha_input.mode_continuity_status(),
            },
        )

        super().__init__(
            input_data=input_data,
            options=asdict(self.qha_options),
            observer=observer,
        )

    def prepare(self) -> None:
        """Prepare the quasi-harmonic approximation calculation.

        Raises
        ------
        ValueError
            If the QHA input data or options are incomplete or inconsistent.
        """
        self.emit("Preparing QHA calculation", level=EventLevel.DEBUG)
        self.qha_input.validate_shapes()
        self.qha_options.validate()
        if self.kieffer_cutoffs is not None:
            if self.qha_options.scheme != "td":
                raise ValueError(
                    "Kieffer QHA currently requires scheme='td'; frequency-scheme "
                    "cutoff evaluation at equilibrium volume is not yet available"
                )
            validate_kieffer_qha_applicability(
                self.qha_input,
                self.kieffer_cutoffs,
            )

        if self.qha_options.requires_mode_continuity():
            status = self.qha_input.mode_continuity_status()
            if not self.qha_input.has_mode_continuity_data():
                self.add_warning(
                    "mode-continuity data are not available; the workflow may "
                    "fall back to static-energy minimization"
                )
            elif status in {"unknown", "unreliable"}:
                raise ValueError(
                    f"frequency QHA requires verified or assumed mode continuity; "
                    f"received '{status}'"
                )
            elif status == "assumed":
                self.add_warning(
                    "phonon-mode continuity is assumed but not explicitly verified"
                )

        self.emit(
            "QHA input summary",
            level=EventLevel.RESULT,
            data={
                "kind": "input_summary",
                "jobname": self.qha_input.jobname,
                "natoms": int(self.qha_input.natoms),
                "formula_units": int(self.qha_input.formula_units),
                "natoms_per_formula_unit": float(
                    self.qha_input.natoms_per_formula_unit
                ),
                "nvol": int(self.qha_input.nvol),
                "qpoints": int(self.qha_input.qpoints),
                "nmodes": int(self.qha_input.nmodes),
                "total_q_points": float(self.qha_input.total_q_points),
                "has_phonons": self.qha_input.has_phonons(),
                "has_mode_continuity": self.qha_input.has_mode_continuity(),
                "mode_continuity": self.qha_input.mode_continuity_status(),
                "source": (
                    str(self.qha_input.source)
                    if self.qha_input.source is not None
                    else None
                ),
            },
        )

        self.emit(
            "QHA static data",
            level=EventLevel.RESULT,
            data={
                "kind": "static_data",
                "volume": self.qha_input.volume,
                "energy": self.qha_input.energy,
                "volume_unit": self.qha_options.volume_unit,
                "energy_unit": self.qha_options.energy_unit,
            },
        )

        self.emit(
            "QHA settings selected",
            level=EventLevel.RESULT,
            data={
                "kind": "settings",
                "options": self.qha_options,
            },
        )

        try:
            preview = pressure_volume_preview(self.qha_input, self.qha_options)
        except ValueError as exc:
            self.emit(
                f"Pressure-volume inspection skipped: {exc}",
                level=EventLevel.WARNING,
            )
        else:
            self.emit(
                "Pressure-volume inspection completed",
                level=EventLevel.RESULT,
                data={
                    "kind": "pressure_volume_preview",
                    "preview": preview,
                },
            )
            for warning in preview.warnings:
                self.emit(warning, level=EventLevel.WARNING)

    def run(self) -> ResultData:
        """Run the QHA workflow.

        Returns
        -------
        ResultData
            Generic Quantas result containing the QHA result object and the
            main numerical arrays.
        """
        self.emit("Starting quasi-harmonic analysis", level=EventLevel.INFO)
        if self.qha_options.scheme == "freq":
            self.emit(
                "Fitting phonon frequencies as a function of volume",
                level=EventLevel.INFO,
            )
            self.emit(
                "Frequency-volume fit diagnostics",
                level=EventLevel.RESULT,
                data={
                    "kind": "frequency_fit_report",
                    "input": self.qha_input,
                    "options": self.qha_options,
                },
            )
        else:
            self.emit(
                "Calculating harmonic thermodynamic properties at sampled volumes",
                level=EventLevel.INFO,
            )

        sampled_thermodynamics = None
        free_energy = None
        frequency_evaluator = None
        try:
            sampled_thermodynamics = calculate_sampled_thermodynamics(
                self.qha_input,
                self.qha_options,
                kieffer_cutoffs=self.kieffer_cutoffs,
            )
            free_energy = free_energy_grid(sampled_thermodynamics)
            if self.qha_options.scheme == "freq":
                frequency_evaluator = FrequencyThermodynamicEvaluator(
                    self.qha_input,
                    self.qha_options,
                )
        except ValueError as exc:
            self.add_warning(
                "harmonic thermodynamics could not be evaluated; "
                f"falling back to static-energy minimization ({exc})"
            )

        if self.qha_options.scheme == "td" and sampled_thermodynamics is not None:
            self.emit(
                "Fitting volume-dependent thermodynamic functions",
                level=EventLevel.INFO,
            )
            self.emit(
                "Thermodynamic fit diagnostics",
                level=EventLevel.RESULT,
                data={
                    "kind": "thermodynamic_fit_report",
                    "sampled_thermodynamics": sampled_thermodynamics,
                    "options": self.qha_options,
                },
            )

        self.emit(
            "Minimizing the pressure-shifted Helmholtz free energy",
            level=EventLevel.INFO,
        )
        local_free_energy_evaluator = None
        if frequency_evaluator is not None:

            def local_free_energy_evaluator(
                volume,
                temperature,
                temperature_index,
            ):
                del temperature_index
                return frequency_evaluator.free_energy_at(volume, temperature)

        qha_result = run_volume_minimization_workflow(
            self.qha_input,
            self.qha_options,
            free_energy=free_energy,
            local_free_energy_evaluator=local_free_energy_evaluator,
            callback=self._workflow_event,
        )
        if sampled_thermodynamics is not None:
            qha_result.kieffer_sampled_contribution = (
                sampled_thermodynamics.kieffer_contribution
            )
            if sampled_thermodynamics.kieffer_contribution is not None:
                qha_result.metadata["kieffer"] = dict(
                    sampled_thermodynamics.kieffer_contribution.metadata
                )

        if (
            sampled_thermodynamics is not None
            and sampled_thermodynamics.volume is not None
            and sampled_thermodynamics.entropy is not None
        ):
            self.emit(
                "Calculating thermal expansion from the free-energy surface",
                level=EventLevel.INFO,
            )
            try:
                calculate_mixed_thermal_expansion(
                    qha_result,
                    sampled_thermodynamics.volume,
                    sampled_thermodynamics.entropy,
                    self.qha_options,
                )
            except ValueError as exc:
                self.add_warning(
                    "mixed-derivative thermal expansion is unavailable; "
                    f"the numerical fallback will be used ({exc})"
                )

        if sampled_thermodynamics is not None:
            if self.qha_options.scheme == "freq":
                self.emit(
                    "Recalculating thermodynamic properties from frequencies at equilibrium volumes",
                    level=EventLevel.INFO,
                )
                calculate_frequency_thermodynamics_at_equilibrium(
                    qha_result,
                    self.qha_input,
                    self.qha_options,
                    evaluator=frequency_evaluator,
                )
            else:
                self.emit(
                    "Interpolating thermodynamic properties at equilibrium volumes",
                    level=EventLevel.INFO,
                )
                interpolate_thermodynamics_at_equilibrium(
                    qha_result,
                    sampled_thermodynamics,
                    self.qha_options,
                )

        if (
            self.qha_options.scheme == "freq"
            and (
                self.qha_options.calculate_mode_gruneisen
                or self.qha_options.thermal_expansion_method == "mode_gruneisen"
            )
            and frequency_evaluator is not None
            and qha_result.equilibrium_volume is not None
        ):
            self.emit(
                "Calculating mode-resolved Gruneisen parameters",
                level=EventLevel.INFO,
            )
            try:
                calculate_mode_gruneisen_at_equilibrium(
                    qha_result,
                    self.qha_input,
                    self.qha_options,
                    frequency_evaluator=frequency_evaluator,
                )
            except ValueError as exc:
                self.add_warning(f"mode Gruneisen analysis skipped: {exc}")
            else:
                refresh_thermal_expansion_dependents(
                    qha_result,
                    self.qha_options,
                )
                gruneisen_metadata = qha_result.metadata.get("gruneisen", {})
                failed_modes = gruneisen_metadata.get("failed_modes", [])
                if failed_modes:
                    self.add_warning(
                        f"mode Gruneisen analysis excluded {len(failed_modes)} "
                        "non-positive or failed modes"
                    )
                self.emit(
                    "Mode Gruneisen analysis completed",
                    level=EventLevel.RESULT,
                    data={
                        "kind": "gruneisen_summary",
                        "metadata": gruneisen_metadata,
                    },
                )

        if self.qha_input.structure is not None:
            self.emit(
                "Calculating lattice parameters and anisotropic thermal expansion",
                level=EventLevel.INFO,
            )
            try:
                calculate_structural_thermal_expansion(
                    qha_result,
                    self.qha_input.structure,
                    self.qha_options,
                )
            except ValueError as exc:
                self.add_warning(
                    f"structural thermal-expansion analysis skipped: {exc}"
                )
            else:
                structural_metadata = qha_result.metadata.get(
                    "structural_thermal_expansion",
                    {},
                )
                extrapolated = int(structural_metadata.get("extrapolated_points", 0))
                if extrapolated and self.qha_options.extrapolation_policy == "warn":
                    self.add_warning(
                        "structural-path interpolation extrapolated "
                        f"{extrapolated} pressure-temperature points"
                    )
                self.emit(
                    "Structural thermal-expansion analysis completed",
                    level=EventLevel.RESULT,
                    data={
                        "kind": "structural_thermal_expansion_summary",
                        "metadata": structural_metadata,
                    },
                )

        if qha_result.completed:
            self.emit(
                "QHA workflow completed",
                level=EventLevel.RESULT,
                data={
                    "kind": "workflow_summary",
                    "workflow": qha_result.metadata.get("workflow", {}),
                },
            )
        else:
            self.add_warning(
                qha_result.metadata.get("stop", {}).get(
                    "message", "QHA workflow stopped before completion"
                )
            )

        self.emit(
            "Final QHA pressure-temperature results",
            level=EventLevel.RESULT,
            data={
                "kind": "final_results",
                "result": qha_result,
            },
        )

        return self._build_result_data(qha_result)

    def _workflow_event(self, event: QHAWorkflowEvent) -> None:
        """Convert a workflow event into a Quantas event.

        Parameters
        ----------
        event : QHAWorkflowEvent
            Structured event emitted by the workflow controller.
        """
        level = self._event_level(event)
        data: dict[str, Any] = {
            "kind": f"qha_{event.kind}",
            "temperature": event.temperature,
            "pressure": event.pressure,
            "temperature_unit": self.qha_options.temperature_unit,
            "pressure_unit": self.qha_options.pressure_unit,
            **event.data,
        }

        if event.kind == "warning":
            warning = event.message
            if warning not in self.warnings:
                self.warnings.append(warning)
        elif event.kind in {"point_failed", "stopped"}:
            warning = event.message
            if warning not in self.warnings:
                self.warnings.append(warning)

        self.emit(
            event.message,
            level=level,
            progress=event.progress if event.kind != "fit_record" else None,
            data=data,
        )

    def _event_level(self, event: QHAWorkflowEvent) -> EventLevel:
        """Return the Quantas event level associated with a workflow event.

        Parameters
        ----------
        event : QHAWorkflowEvent
            Structured event emitted by the workflow controller.

        Returns
        -------
        EventLevel
            Event level used by the calculator observer.
        """
        if event.kind in {"started", "point_started"}:
            return EventLevel.DEBUG
        if event.kind == "fit_record":
            return EventLevel.DEBUG
        if event.kind == "warning":
            return EventLevel.WARNING
        if event.kind in {"point_failed", "stopped"}:
            return EventLevel.WARNING
        if event.kind == "point_completed":
            return EventLevel.PROGRESS
        if event.kind == "completed":
            return EventLevel.RESULT
        return EventLevel.INFO

    def _build_result_data(self, qha_result: QHAResult) -> ResultData:
        """Build the generic Quantas result object.

        Parameters
        ----------
        qha_result : QHAResult
            Quasi-harmonic approximation result object.

        Returns
        -------
        ResultData
            Generic Quantas result data container.
        """
        return ResultData(
            metadata=ResultMetadata(module=self.module, method=self.method),
            input_data=self.input_data,
            options=self.options,
            results={"qha": qha_result},
            warnings=self.warnings,
        )
