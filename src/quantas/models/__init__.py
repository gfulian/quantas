# -*- coding: utf-8 -*-

"""Shared data contracts and active-object base classes for Quantas.

The package exposes the frontend-neutral interfaces used by calculators,
readers, exporters, scientific results, reports, plots, and persistence.
"""

from .calculator import BasicCalculator
from .data import InputData, ResultData, ResultMetadata
from .export import BasicExport, BasicHDF5Export
from .elastic_states import (
    ElasticState,
    ElasticStateSeries,
    ElasticTensorKind,
    PressureSource,
    PrestressProvenance,
)
from .kieffer import CutoffVolumeSource, KiefferCutoffState, KiefferVolumeSeries
from .module import ModuleContract
from .phonons import PhononInputData, PhononModeData
from .plot import (
    AxisFieldLayer,
    ColoredPathSeries,
    ColoredPathStyle,
    ContourPlotSpec,
    LinePlotSpec,
    LineStyle,
    PanelPlotSpec,
    PlotAxis,
    PlotBand,
    PlotBandStyle,
    PlotCollection,
    PlotMask,
    PlotSeries,
    PlotSeriesStyle,
    PlotSpan,
    ScalarBackground,
    SecondaryAxis,
    PlotSpec,
    PolarPlotPanel,
    PolarPlotSpec,
    SurfaceLayer,
    SurfacePlotSpec,
    SurfaceStyle,
    SphericalMapSpec,
    SphericalMarker,
    SphericalProjection,
    SphericalSummarySpec,
    VectorFieldLayer,
    VectorFieldStyle,
)
from .plot_inventory import (
    PlotContextDescriptor,
    PlotContextValue,
    PlotInventory,
    PlotKind,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
)
from .reader import BasicReader
from .structures import (
    CellNormalization,
    CrystalStructure,
    StructureReconstructionDiagnostics,
    StructureVolumeSeries,
    SymmetryMetadata,
)
from .report import ReportTable, input_data_table, mapping_table
from .schema import (
    REQUIRED_RESULT_GROUPS,
    RESULT_SCHEMA_VERSION,
    ResultSchema,
    validate_result_schema,
)
from .thermodynamics import HarmonicThermodynamicResult
from .volume_matching import VolumeMatch, VolumeMatchPolicy, match_sampled_volumes

__all__ = [
    "BasicCalculator",
    "BasicExport",
    "BasicHDF5Export",
    "BasicReader",
    "CellNormalization",
    "CrystalStructure",
    "CutoffVolumeSource",
    "ElasticState",
    "ElasticStateSeries",
    "ElasticTensorKind",
    "HarmonicThermodynamicResult",
    "StructureReconstructionDiagnostics",
    "StructureVolumeSeries",
    "SymmetryMetadata",
    "InputData",
    "KiefferCutoffState",
    "KiefferVolumeSeries",
    "ModuleContract",
    "PhononInputData",
    "PhononModeData",
    "PressureSource",
    "PrestressProvenance",
    "AxisFieldLayer",
    "ColoredPathSeries",
    "ColoredPathStyle",
    "ContourPlotSpec",
    "LinePlotSpec",
    "LineStyle",
    "PanelPlotSpec",
    "PlotAxis",
    "PlotBand",
    "PlotBandStyle",
    "PlotCollection",
    "PlotMask",
    "PlotSeries",
    "PlotSeriesStyle",
    "PlotSpan",
    "ScalarBackground",
    "SecondaryAxis",
    "PlotSpec",
    "PlotContextDescriptor",
    "PlotContextValue",
    "PlotInventory",
    "PlotKind",
    "PlotPropertyDescriptor",
    "PlotRepresentationDescriptor",
    "PolarPlotPanel",
    "PolarPlotSpec",
    "SurfaceLayer",
    "SurfacePlotSpec",
    "SurfaceStyle",
    "SphericalMapSpec",
    "SphericalMarker",
    "SphericalProjection",
    "SphericalSummarySpec",
    "VectorFieldLayer",
    "VectorFieldStyle",
    "ResultData",
    "ReportTable",
    "input_data_table",
    "mapping_table",
    "ResultMetadata",
    "REQUIRED_RESULT_GROUPS",
    "RESULT_SCHEMA_VERSION",
    "ResultSchema",
    "validate_result_schema",
    "VolumeMatch",
    "VolumeMatchPolicy",
    "match_sampled_volumes",
]
