# -*- coding: utf-8 -*-

"""Backend-neutral fitting solver protocol."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from .model import BaseFitModel
from .observations import FitObservations
from .options import FitOptions
from .parameters import ParameterMap
from .result import FitResult


@runtime_checkable
class FitSolver(Protocol):
    """Protocol implemented by general numerical fitting strategies."""

    def fit_problem(
        self,
        model: BaseFitModel,
        observations: FitObservations,
        parameters: ParameterMap,
        options: FitOptions,
    ) -> FitResult:
        """Fit one fully specified numerical problem and return diagnostics."""
        ...
