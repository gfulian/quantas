# -*- coding: utf-8 -*-

"""Calculator for sampled seismic-wave workflows."""

from __future__ import annotations

from dataclasses import asdict

from quantas.core.events import EventLevel, Observer
from quantas.core.geometry import tensor_frame_mapping
from quantas.models import BasicCalculator, InputData, ResultData, ResultMetadata
from quantas.modules.seismic.analysis import run_analysis, validate_input
from quantas.modules.seismic.models import SeismicInput, SeismicOptions, SeismicResult


class SeismicCalculator(BasicCalculator):
    """Run a spherical acoustic-wave propagation analysis.

    Parameters
    ----------
    seismic_input : SeismicInput
        Elastic stiffness matrix and material density.
    options : SeismicOptions or None, optional
        Scientific sampling and numerical options.
    observer : Observer or None, optional
        Observer receiving Quantas events.
    """

    module = "seismic"
    method = "christoffel"

    def __init__(
        self,
        seismic_input: SeismicInput,
        options: SeismicOptions | None = None,
        observer: Observer | None = None,
    ) -> None:
        self.seismic_input = seismic_input
        self.seismic_options = options or SeismicOptions()
        input_data = InputData(
            source=seismic_input.source,
            raw=seismic_input.raw,
            data={
                "jobname": seismic_input.jobname,
                "density": seismic_input.density,
                "stiffness": seismic_input.stiffness,
                "stiffness_unit": "GPa",
                "density_unit": "kg m^-3",
                "tensor_frame": tensor_frame_mapping(self.seismic_options.rotation),
            },
        )
        option_data = asdict(self.seismic_options)
        option_data["hemisphere"] = self.seismic_options.hemisphere.value
        option_data["level"] = self.seismic_options.level.value
        super().__init__(
            input_data=input_data,
            options=option_data,
            observer=observer,
        )

    def prepare(self) -> None:
        """Validate the input and emit the selected scientific settings."""
        self.emit("Preparing seismic calculation", level=EventLevel.DEBUG)
        validate_input(self.seismic_input, self.seismic_options)
        self.emit(
            "Seismic settings selected",
            level=EventLevel.RESULT,
            data={"kind": "settings", "options": self.seismic_options},
        )

    def run(self) -> ResultData:
        """Run the acoustic-wave workflow and return generic results."""
        self.emit(
            "Seismic input data loaded",
            level=EventLevel.RESULT,
            data={"kind": "input", "input": self.seismic_input},
        )
        if self.seismic_options.rotation is not None:
            self.emit("Elastic tensor components transformed to the analysis frame")
        self.emit("Sampling acoustic-wave fields")
        result = run_analysis(
            self.seismic_input,
            self.seismic_options,
            progress_callback=self._emit_progress,
        )
        self.emit(
            "Isotropic reference velocities calculated",
            level=EventLevel.RESULT,
            data={
                "kind": "isotropic_reference",
                "velocities": result.isotropic_velocities,
            },
        )
        self.emit(
            "Seismic fields calculated",
            level=EventLevel.RESULT,
            data={
                "kind": "seismic_field",
                "level": result.field.level.value,
                "n_points": result.field.n_points,
                "diagnostics": result.metadata,
            },
        )
        self._add_diagnostic_warnings(result)
        return ResultData(
            metadata=ResultMetadata(module=self.module, method=self.method),
            input_data=self.input_data,
            options=self.options,
            results={"seismic": result},
            warnings=list(self.warnings),
        )

    def _emit_progress(self, current: int, total: int) -> None:
        """Convert numerical sampling progress into a Quantas event."""
        progress = float(current / total) if total else 1.0
        self.emit(
            "Sampling acoustic-wave fields",
            level=EventLevel.PROGRESS,
            progress=progress,
            data={
                "kind": "sampling_progress",
                "current": current,
                "total": total,
            },
        )

    def _add_diagnostic_warnings(self, result: SeismicResult) -> None:
        """Convert non-fatal field diagnostics into workflow warnings."""
        invalid = int(result.metadata["invalid_phase_points"])
        if invalid:
            self.add_warning(
                f"The sampled field contains {invalid} invalid phase-mode values."
            )
        if result.field.enhancement is not None:
            non_finite = int(result.metadata["non_finite_enhancement_points"])
            if non_finite:
                self.add_warning(
                    "The sampled field contains "
                    f"{non_finite} non-finite enhancement values."
                )
