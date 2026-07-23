# -*- coding: utf-8 -*-

"""Shared numerical policies used across Quantas workflows."""

from .rectilinear import (
    RectilinearFieldInterpolator,
    grid_step,
    regular_grid,
    validated_axis,
)

from .precision import (
    COMPLEX_DTYPE,
    FLOAT_DTYPE,
    NumericPrecisionPolicy,
    cast_floating_array,
    cast_floating_scalar,
    default_precision_policy,
    precision_policy_from_options,
)

__all__ = [
    "COMPLEX_DTYPE",
    "FLOAT_DTYPE",
    "NumericPrecisionPolicy",
    "RectilinearFieldInterpolator",
    "cast_floating_array",
    "cast_floating_scalar",
    "default_precision_policy",
    "grid_step",
    "precision_policy_from_options",
    "regular_grid",
    "validated_axis",
]
