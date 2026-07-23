# -*- coding: utf-8 -*-

"""Calculator for second-order elasticity workflows."""

from __future__ import annotations

from dataclasses import asdict

from quantas.core.events import EventLevel, Observer
from quantas.core.geometry import tensor_frame_mapping
from quantas.core.physics.elasticity import sample_elasticity_surfaces
from quantas.models import BasicCalculator, InputData, ResultData, ResultMetadata
from quantas.modules.elasticity.analysis import (
    calculate_basic_properties,
    calculate_directional_variations,
    calculate_2d_properties,
    create_elastic_tensor,
    specialize_tensor,
    validate_input,
)
from quantas.modules.elasticity.models import (
    ElasticityInput,
    ElasticityOptions,
    ElasticityResult,
    ElasticitySurfaceOptions,
)


class ElasticityCalculator(BasicCalculator):
    """Run a second-order elastic-tensor analysis.

    Parameters
    ----------
    elasticity_input : ElasticityInput
        Normalized stiffness-matrix input.
    options : ElasticityOptions or None, optional
        Scientific workflow options.
    observer : Observer or None, optional
        Observer receiving Quantas events.
    """

    module = "elasticity"
    method = "second_order"

    def __init__(
        self,
        elasticity_input: ElasticityInput,
        options: ElasticityOptions | None = None,
        observer: Observer | None = None,
    ) -> None:
        self.elasticity_input = elasticity_input
        self.elasticity_options = options or ElasticityOptions()
        input_data = InputData(
            source=elasticity_input.source,
            data={
                "jobname": elasticity_input.jobname,
                "stiffness": elasticity_input.stiffness,
                "stiffness_unit": self.elasticity_options.pressure_unit,
                "tensor_frame": tensor_frame_mapping(self.elasticity_options.rotation),
            },
        )
        super().__init__(
            input_data=input_data,
            options=asdict(self.elasticity_options),
            observer=observer,
        )

    def prepare(self) -> None:
        """Validate input data and emit the selected scientific settings."""
        self.emit("Preparing elasticity calculation", level=EventLevel.DEBUG)
        validate_input(self.elasticity_input, self.elasticity_options)
        self.emit(
            "Elasticity settings selected",
            level=EventLevel.RESULT,
            data={"kind": "settings", "options": self.elasticity_options},
        )

    def run(self) -> ResultData:
        """Run the elasticity workflow and return generic Quantas results."""
        self.emit("Creating elastic tensor", level=EventLevel.DEBUG)
        tensor = create_elastic_tensor(
            self.elasticity_input,
            self.elasticity_options.rotation,
        )
        if self.elasticity_options.rotation is not None:
            self.emit(
                "Elastic tensor components transformed to the analysis frame",
                level=EventLevel.INFO,
            )
            self.emit(
                "Elastic tensor rotation applied",
                level=EventLevel.RESULT,
                data={
                    "kind": "rotation",
                    "source_stiffness": self.elasticity_input.stiffness,
                    "analysis_stiffness": tensor.stiffness,
                    "rotation": self.elasticity_options.rotation,
                },
            )

        self.emit("Calculating basic elastic properties", level=EventLevel.DEBUG)
        elasticity_result = calculate_basic_properties(
            tensor,
            self.elasticity_input,
        )
        elasticity_result.metadata["tensor_frame"] = self.input_data.data[
            "tensor_frame"
        ]

        self.emit(
            "Elasticity input data loaded",
            level=EventLevel.RESULT,
            data={
                "kind": "input",
                "input": self.elasticity_input,
                "result": elasticity_result,
            },
        )
        self.emit(
            "Elastic averages calculated",
            level=EventLevel.RESULT,
            data={"kind": "averages", "result": elasticity_result},
        )
        self.emit(
            "Mechanical stability evaluated",
            level=EventLevel.RESULT,
            data={"kind": "stability", "result": elasticity_result},
        )

        if (
            elasticity_result.stability is None
            or not elasticity_result.stability.is_stable
        ):
            self.add_warning(
                "The elastic stiffness matrix is not positive definite. "
                "The structure is mechanically unstable; directional analysis "
                "skipped."
            )
            return self._build_result_data(elasticity_result)

        tensor = specialize_tensor(tensor, elasticity_result)
        elasticity_result.metadata["tensor_class"] = type(tensor).__name__

        self.emit("Calculating directional extrema")
        calculate_directional_variations(tensor, elasticity_result)
        self.emit(
            "Directional extrema calculated",
            level=EventLevel.RESULT,
            data={"kind": "variations", "result": elasticity_result},
        )

        if self.elasticity_options.calculate_2d:
            self.emit("Calculating two-dimensional elastic properties:")
            calculate_2d_properties(
                tensor,
                elasticity_result,
                self.elasticity_options,
                progress_callback=self._two_dimensional_progress,
                step_callback=self._two_dimensional_step,
            )

        if self.elasticity_options.calculate_3d:
            surface_options = (
                self.elasticity_options.surface_options or ElasticitySurfaceOptions()
            )
            self.emit("Calculating three-dimensional elastic properties")
            elasticity_result.properties_3d = sample_elasticity_surfaces(
                tensor,
                ntheta=surface_options.ntheta,
                nphi=surface_options.nphi,
                properties=surface_options.properties,
                batch_size=surface_options.batch_size,
                progress_callback=self._three_dimensional_progress,
            )
            self.emit(
                "Three-dimensional elastic properties calculated",
                level=EventLevel.RESULT,
                data={
                    "kind": "properties_3d",
                    "surface_options": surface_options,
                    "surface_count": len(elasticity_result.properties_3d.surfaces),
                },
            )

        return self._build_result_data(elasticity_result)

    def _build_result_data(self, result: ElasticityResult) -> ResultData:
        """Build the generic Quantas result container."""
        return ResultData(
            metadata=ResultMetadata(module=self.module, method=self.method),
            input_data=self.input_data,
            options=self.options,
            results={"elasticity": result},
            warnings=self.warnings,
        )

    def _two_dimensional_progress(self, label: str, current: int, total: int) -> None:
        """Convert numerical progress to a Quantas progress event."""
        self.emit(
            label,
            level=EventLevel.PROGRESS,
            progress=current / total,
            data={"current": current, "total": total, "label": label},
        )

    def _three_dimensional_progress(self, current: int, total: int) -> None:
        """Convert 3D batch progress to a Quantas progress event."""
        self.emit(
            "Three-dimensional elasticity sampling",
            level=EventLevel.PROGRESS,
            progress=current / total,
            data={"current": current, "total": total},
        )

    def _two_dimensional_step(self, plane: str, property_name: str) -> None:
        """Emit the beginning of one principal-plane calculation."""
        label = property_name.replace("_", " ")
        self.emit(
            f"  - {label} on ({plane}) plane",
            level=EventLevel.INFO,
            data={
                "kind": "2d_step",
                "plane": plane,
                "property": property_name,
            },
        )
