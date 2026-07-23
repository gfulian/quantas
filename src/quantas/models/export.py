# -*- coding: utf-8 -*-

"""Base contracts for writing Quantas results to external formats.

Module-specific exporters implement scientific payload details while the shared
classes define common path handling and native-result behavior.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from .data import ResultData


class BasicExport(ABC):
    """
    Basic class for Quantas exporters.

    Attributes
    ----------
    completed : bool
        Flag that is set to ``True`` when the export operation is completed.
    error : str or None
        Error message generated during export, if any.
    """

    completed: bool = False
    error: str | None = None

    @abstractmethod
    def export(self, result: ResultData, filename: str | Path) -> None:
        """
        Export a Quantas result.

        Parameters
        ----------
        result : ResultData
            Result object to be exported.
        filename : str or Path
            Path to the output file.
        """
        raise NotImplementedError


class BasicHDF5Export(BasicExport):
    """
    Basic class for HDF5 exporters.

    This class preserves the historical Quantas idea of using HDF5 files as
    native binary containers for calculation results. Derived classes should
    implement the actual layout of the HDF5 file according to the module that
    generated the result.
    """
