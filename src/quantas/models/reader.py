# -*- coding: utf-8 -*-

"""Base contracts for readers that convert external files into Quantas inputs.

Readers parse and validate source data without running scientific workflows or
depending on command-line and graphical frontends.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, TypeVar


ReaderResult = TypeVar("ReaderResult")


class BasicReader(ABC, Generic[ReaderResult]):
    """
    Basic class for Quantas input readers.

    Parameters
    ----------
    completed : bool, optional
        Initial completion state. The default is ``False``.
    error : str or None, optional
        Initial error message. The default is ``None``.

    Attributes
    ----------
    completed : bool
        Flag that is set to ``True`` when the input file has been completely
        read.
    error : str or None
        Error message generated while reading the input file, if any.

    Notes
    -----
    Readers are active, stateful objects. Their identity is therefore distinct
    from passive data contracts even when two instances currently expose the
    same completion and error state.
    """

    def __init__(
        self,
        completed: bool = False,
        error: str | None = None,
    ) -> None:
        """Initialize the shared reader state.

        Parameters
        ----------
        completed : bool, optional
            Initial completion state.
        error : str or None, optional
            Initial error message.
        """
        self.completed = completed
        self.error = error

    @abstractmethod
    def load(self, filename: str | Path) -> ReaderResult:
        """Load and validate one input source.

        Concrete readers define whether recoverable parse failures are returned
        through a result object or recorded in :attr:`error`; callers should
        consult the reader-specific contract for that distinction.

        Parameters
        ----------
        filename : str or Path
            Path to the input file.

        Returns
        -------
        ReaderResult
            Parsed reader result defined by the concrete reader.

        Raises
        ------
        NotImplementedError
            Always in the abstract base implementation.
        """
        raise NotImplementedError
