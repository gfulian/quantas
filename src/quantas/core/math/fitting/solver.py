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
        """Fit one fully specified numerical problem.

        Parameters
        ----------
        model : BaseFitModel
            Complete mathematical or physical model.
        observations : FitObservations
            Coordinates, responses, optional uncertainties, and selection mask.
        parameters : ParameterMap
            FREE/FIXED/IMPLIED/DERIVED parameter contract.
        options : FitOptions
            Method selection and backend-neutral numerical settings.

        Returns
        -------
        FitResult
            Complete physical parameters and method-neutral diagnostics. Expected
            input or numerical failures are represented by an unsuccessful result;
            unexpected programming errors may propagate from concrete solvers.
        """
        ...
