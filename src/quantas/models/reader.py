# -*- coding: utf-8 -*-

"""Base contracts for readers that convert external files into Quantas inputs.

Readers parse and validate source data without running scientific workflows or
depending on command-line and graphical frontends.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Generic, TypeVar


ReaderResult = TypeVar("ReaderResult")


@dataclass
class BasicReader(ABC, Generic[ReaderResult]):
    """
    Basic class for Quantas input readers.

    Attributes
    ----------
    completed : bool
        Flag that is set to ``True`` when the input file has been completely
        read.
    error : str or None
        Error message generated while reading the input file, if any.
    """

    completed: bool = False
    error: str | None = None

    @abstractmethod
    def load(self, filename: str | Path) -> ReaderResult:
        """
        Load an input file.

        Parameters
        ----------
        filename : str or Path
            Path to the input file.

        Returns
        -------
        ReaderResult
            Parsed reader result.
        """
        raise NotImplementedError
