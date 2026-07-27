# -*- coding: utf-8 -*-

"""Frontend-neutral plotting services for EOS archives."""

from .inventory import (
    EOSArchivePlotInventory,
    EOSDatasetPlotDescriptor,
    EOSRecordPlotDescriptor,
    EOSSlotPlotDescriptor,
    describe_eos_plots,
)
from .labels import format_unit, normalized_pressure_labels, property_label
from .spec import EOS_PLOT_TYPES, EOSPlotOptions, EOSPlotter

__all__ = [
    "EOS_PLOT_TYPES",
    "EOSArchivePlotInventory",
    "EOSDatasetPlotDescriptor",
    "EOSPlotOptions",
    "EOSRecordPlotDescriptor",
    "EOSSlotPlotDescriptor",
    "EOSPlotter",
    "describe_eos_plots",
    "format_unit",
    "normalized_pressure_labels",
    "property_label",
]
