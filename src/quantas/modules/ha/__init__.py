# -*- coding: utf-8 -*-

"""Public harmonic-approximation workflow API."""

from .api import (
    MODULE_CONTRACT,
    build_ha_plots,
    build_ha_report,
    create_ha_input,
    normalize_ha_input,
    read_ha_hdf5,
    read_ha_input,
    write_ha_hdf5,
    run_ha,
)
from .calculator import HACalculator
from .models import HAInput, HAOptions, HAResult

__all__ = [
    "HACalculator",
    "HAInput",
    "HAOptions",
    "HAResult",
    "MODULE_CONTRACT",
    "build_ha_plots",
    "build_ha_report",
    "create_ha_input",
    "normalize_ha_input",
    "read_ha_hdf5",
    "read_ha_input",
    "write_ha_hdf5",
    "run_ha",
]
