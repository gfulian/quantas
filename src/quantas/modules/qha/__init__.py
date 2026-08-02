# -*- coding: utf-8 -*-

"""Public quasi-harmonic approximation workflow API."""

from .api import (
    MODULE_CONTRACT,
    build_qha_inspection_plots,
    build_qha_inspection_report,
    build_qha_plots,
    build_qha_report,
    inspect_qha_input,
    normalize_qha_input,
    read_qha_hdf5,
    read_qha_input,
    write_qha_hdf5,
    run_qha,
)
from .calculator import QHACalculator
from .models import QHAInput, QHAOptions, QHAResult
from .plot import QHAPlotOptions
from .structural import calculate_structural_thermal_expansion
from .validation import (
    PropertyDifference,
    QHAValidationSummary,
    compare_qha_results,
    validate_qha_result,
)

__all__ = [
    "MODULE_CONTRACT",
    "PropertyDifference",
    "QHACalculator",
    "QHAInput",
    "QHAOptions",
    "QHAPlotOptions",
    "QHAResult",
    "QHAValidationSummary",
    "build_qha_plots",
    "build_qha_report",
    "build_qha_inspection_plots",
    "build_qha_inspection_report",
    "calculate_structural_thermal_expansion",
    "compare_qha_results",
    "inspect_qha_input",
    "normalize_qha_input",
    "read_qha_hdf5",
    "read_qha_input",
    "write_qha_hdf5",
    "run_qha",
    "validate_qha_result",
]
