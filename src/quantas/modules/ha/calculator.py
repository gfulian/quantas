# -*- coding: utf-8 -*-

"""
Calculator for harmonic-approximation workflows.

The :class:`HACalculator` class connects the library-level HA workflow to the
common Quantas calculator model. It validates input data, runs the harmonic
thermodynamic analysis, converts numerical progress callbacks into Quantas
events, and returns a generic :class:`~quantas.models.ResultData` object.

No command-line or graphical frontend code is used in this module.
"""

from __future__ import annotations

from dataclasses import asdict

from quantas.core.events import EventLevel, Observer
from quantas.models import BasicCalculator, InputData, ResultData, ResultMetadata
from quantas.modules.ha.analysis import (
    calculate_thermodynamic_properties,
    validate_input,
)
from quantas.modules.ha.models import HAInput, HAOptions, HAResult


class HACalculator(BasicCalculator):
    """
    Calculator for harmonic-approximation thermodynamic properties.

    Parameters
    ----------
    ha_input : HAInput
        Normalized input data used by the HA workflow.
    options : HAOptions or None, optional
        Options controlling the calculation. If ``None``, default HA options
        are used.
    observer : Observer or None, optional
        Observer receiving workflow events. If ``None``, a null observer is
        used by the base calculator.
    """

    module = "ha"
    method = "harmonic"

    def __init__(
        self,
        ha_input: HAInput,
        options: HAOptions | None = None,
        observer: Observer | None = None,
    ) -> None:
        """
        Initialize the HA calculator.
        """
        self.ha_input = ha_input
        self.ha_options = options if options is not None else HAOptions()
        self._timing_records: list[dict[str, float | str]] = []

        input_data = InputData(
            source=ha_input.source,
            data={
                "jobname": ha_input.jobname,
                "natoms": ha_input.natoms,
                "formula_units": ha_input.formula_units,
                "natoms_per_formula_unit": ha_input.natoms_per_formula_unit,
                "qpoints": ha_input.qpoints,
                "nvol": ha_input.nvol,
                "nmodes": ha_input.nmodes,
                "total_q_points": ha_input.total_q_points,
                "has_structure": ha_input.has_structure(),
                "structure": (
                    None
                    if ha_input.structure is None
                    else ha_input.structure.as_dict(include_source=True)
                ),
            },
        )

        super().__init__(
            input_data=input_data,
            options=asdict(self.ha_options),
            observer=observer,
        )

    def prepare(self) -> None:
        """
        Prepare the harmonic-approximation calculation.

        Raises
        ------
        ValueError
            If the HA input data are incomplete or inconsistent.
        """
        self.emit("Preparing HA calculation", level=EventLevel.DEBUG)
        validate_input(self.ha_input)

        self.emit(
            "HA input summary",
            level=EventLevel.RESULT,
            data={
                "kind": "input_summary",
                "jobname": self.ha_input.jobname,
                "natoms": int(self.ha_input.natoms),
                "formula_units": int(self.ha_input.formula_units),
                "natoms_per_formula_unit": float(self.ha_input.natoms_per_formula_unit),
                "nvol": int(self.ha_input.nvol),
                "qpoints": int(self.ha_input.qpoints),
                "nmodes": int(self.ha_input.nmodes),
                "total_q_points": float(self.ha_input.total_q_points),
                "source": (
                    str(self.ha_input.source)
                    if self.ha_input.source is not None
                    else None
                ),
            },
        )

        self.emit(
            "HA static data",
            level=EventLevel.RESULT,
            data={
                "kind": "static_data",
                "volume": self.ha_input.volume,
                "energy": self.ha_input.energy,
                "volume_unit": self.ha_options.volume_unit,
                "energy_unit": self.ha_options.energy_unit,
            },
        )

        self.emit(
            "HA settings selected",
            level=EventLevel.RESULT,
            data={
                "kind": "settings",
                "options": self.ha_options,
            },
        )

    def run(self) -> ResultData:
        """
        Run the harmonic-approximation workflow.

        Returns
        -------
        ResultData
            Generic Quantas result containing the HA result object and the main
            numerical arrays.
        """
        self.emit("Starting harmonic analysis", level=EventLevel.INFO)

        ha_result = calculate_thermodynamic_properties(
            self.ha_input,
            self.ha_options,
            progress_callback=self._thermodynamic_progress,
            step_callback=self._thermodynamic_step,
            result_callback=self._thermodynamic_result,
            timing_callback=self._thermodynamic_timing,
        )

        self.emit(
            "Harmonic thermodynamic properties calculated",
            level=EventLevel.RESULT,
            data={
                "kind": "thermodynamics",
                "result": ha_result,
            },
        )

        return self._build_result_data(ha_result)

    def _thermodynamic_result(self, label: str, payload: dict) -> None:
        """
        Emit a structured thermodynamic result event.

        Parameters
        ----------
        label : str
            Name of the calculated thermodynamic quantity.
        payload : dict
            Structured numerical result payload produced by ``analysis.py``.
        """
        data = dict(payload)
        data.setdefault("kind", "thermodynamic_property")
        data.setdefault("property", label)

        if data.get("kind") == "thermodynamic_backend" and data.get("warning"):
            self.add_warning(str(data["warning"]))

        self.emit(
            f"Calculated {label.replace('_', ' ')}",
            level=EventLevel.RESULT,
            data=data,
        )

    def _build_result_data(self, ha_result: HAResult) -> ResultData:
        """
        Build the generic Quantas result object.

        Parameters
        ----------
        ha_result : HAResult
            Harmonic-approximation result object.

        Returns
        -------
        ResultData
            Generic Quantas result data container.
        """
        return ResultData(
            metadata=ResultMetadata(module=self.module, method=self.method),
            input_data=self.input_data,
            options=self.options,
            results={"ha": ha_result},
            warnings=self.warnings,
        )

    def _thermodynamic_timing(
        self,
        label: str,
        elapsed_seconds: float,
        backend: str,
    ) -> None:
        """
        Emit timing information for a HA numerical operation.

        Parameters
        ----------
        label : str
            Name of the timed numerical operation.
        elapsed_seconds : float
            Wall-clock elapsed time in seconds.
        backend : str
            Name of the thermodynamic backend used for the operation.
        """
        record: dict[str, float | str] = {
            "label": label,
            "elapsed_seconds": float(elapsed_seconds),
            "backend": backend,
        }
        self._timing_records.append(record)
        self.emit(
            f"Timing {label.replace('_', ' ')}: {elapsed_seconds:.6f} s",
            level=EventLevel.DEBUG,
            data={
                "kind": "timing",
                **record,
            },
        )

    def _thermodynamic_progress(
        self,
        label: str,
        current: int,
        total: int,
    ) -> None:
        """
        Emit progress during harmonic thermodynamic calculations.

        Parameters
        ----------
        label : str
            Label of the completed workflow step.
        current : int
            Number of completed steps.
        total : int
            Total number of workflow steps.
        """
        self.emit(
            label,
            level=EventLevel.PROGRESS,
            progress=current / total,
            data={
                "kind": "thermodynamic_progress",
                "current": current,
                "total": total,
                "label": label,
            },
        )

    def _thermodynamic_step(self, label: str, status: str) -> None:
        """
        Emit an event at the beginning or end of a HA workflow step.

        Parameters
        ----------
        label : str
            Name of the current thermodynamic step.
        status : str
            Step status, for example ``"started"`` or ``"finished"``.
        """
        if status != "started":
            return

        self.emit(
            f"Calculating {label.replace('_', ' ')}",
            level=EventLevel.DEBUG,
            data={
                "kind": "thermodynamic_step",
                "step": label,
                "status": status,
            },
        )
