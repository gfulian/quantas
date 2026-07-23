# -*- coding: utf-8 -*-

"""Uniform public contracts for Quantas scientific modules."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .data import ResultData
from .plot import PlotCollection
from .report import ReportTable


@dataclass(frozen=True, slots=True)
class ModuleContract:
    """Public entry points implemented by a Quantas scientific module.

    Parameters
    ----------
    name : str
        Stable module identifier used by metadata and frontend registries.
    result_key : str
        Key containing the module-specific payload in ``ResultData.results``.
    read_input : callable
        Function reading and normalizing module input.
    run : callable
        Function executing the module workflow.
    read_hdf5 : callable
        Function reconstructing a complete ``ResultData`` from native HDF5.
    write_hdf5 : callable
        Function persisting a complete ``ResultData`` to native HDF5.
    build_report : callable
        Function converting a result to neutral report tables.
    build_plots : callable
        Function converting a result to neutral plot specifications.
    """

    name: str
    result_key: str
    read_input: Callable[..., Any]
    run: Callable[..., ResultData]
    read_hdf5: Callable[..., ResultData]
    write_hdf5: Callable[..., Path]
    build_report: Callable[..., list[ReportTable]]
    build_plots: Callable[..., PlotCollection]

    def validate_result(self, result: ResultData) -> None:
        """Validate that a result belongs to this module.

        Parameters
        ----------
        result : ResultData
            Generic Quantas result.

        Raises
        ------
        ValueError
            If module metadata or the module payload key do not match.
        """
        if result.metadata.module != self.name:
            raise ValueError(
                f"Expected a '{self.name}' result, got '{result.metadata.module}'."
            )
        if self.result_key not in result.results:
            raise ValueError(
                f"Result does not contain the '{self.result_key}' payload."
            )
