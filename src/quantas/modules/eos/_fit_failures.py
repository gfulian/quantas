# -*- coding: utf-8 -*-

"""Shared classification of documented EOS workflow request errors."""

from __future__ import annotations

from quantas.core.math.fitting import FitResult, FitStatus

from .models import EOSFitRequest, EOSFitResult

EXPECTED_EOS_FIT_EXCEPTIONS = (ValueError, TypeError, NotImplementedError)
"""Exceptions documented for invalid or unsupported EOS fit requests."""


def eos_invalid_request_result(
    request: EOSFitRequest,
    exc: Exception,
) -> EOSFitResult:
    """Convert one documented EOS request error to a persistent fit result.

    Parameters
    ----------
    request : EOSFitRequest
        Fit request being executed when validation failed.
    exc : Exception
        Documented request or validation exception.

    Returns
    -------
    EOSFitResult
        Failed result with ``INVALID_INPUT`` status and exception provenance.

    Notes
    -----
    Unexpected exceptions are deliberately not handled here. Numerical
    convergence and supported backend failures are already returned by the
    fitting services as structured :class:`FitResult` objects; an exception
    escaping those services indicates an unexpected workflow failure and must
    remain visible to the caller.
    """
    message = str(exc) or type(exc).__name__
    fit = FitResult.failed(
        message,
        status=FitStatus.INVALID_INPUT,
        method=(
            None
            if request.options.solver_options is None
            else request.options.solver_options.method
        ),
    )
    return EOSFitResult(
        request=request,
        fit=fit,
        warnings=[message],
        metadata={
            "workflow_exception": type(exc).__name__,
            "workflow_failure_classification": "request_error",
        },
    )
