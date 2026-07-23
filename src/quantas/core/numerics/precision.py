# -*- coding: utf-8 -*-

"""Central floating-point precision definition for Quantas.

Quantas scientific workflows use IEEE 754 double precision throughout the
validated numerical stack.  The same precision is used for native HDF5
floating-point datasets and scalar attributes.  Display formatting remains an
independent frontend concern and must never change the stored numerical values.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray


FLOAT_DTYPE = np.dtype(np.float64)
"""Validated real floating-point dtype used by Quantas."""

COMPLEX_DTYPE = np.dtype(np.complex128)
"""Validated complex floating-point dtype used by Quantas."""


@dataclass(frozen=True, slots=True)
class NumericPrecisionPolicy:
    """Fixed numerical precision metadata for Quantas workflows.

    Quantas currently exposes no user-selectable floating-point mode.  All
    scientific calculations and native HDF5 floating-point values use double
    precision.  This passive container exists only to make the policy explicit
    in metadata and reports.
    """

    working: str = field(default="double", init=False)
    storage: str = field(default="double", init=False)
    working_dtype: str = field(default="float64", init=False)
    storage_dtype: str = field(default="float64", init=False)
    working_bits: int = field(default=64, init=False)
    storage_bits: int = field(default=64, init=False)

    @classmethod
    def from_mapping(
        cls,
        values: Mapping[str, Any] | None,
    ) -> "NumericPrecisionPolicy":
        """Return the fixed Quantas precision policy.

        Historical mappings are accepted for compatibility but do not alter
        the validated double-precision policy.
        """
        return cls()

    def as_metadata(self) -> dict[str, str | int]:
        """Return serializable fixed precision metadata."""
        return {
            "working": self.working,
            "working_dtype": self.working_dtype,
            "working_bits": self.working_bits,
            "storage": self.storage,
            "storage_dtype": self.storage_dtype,
            "storage_bits": self.storage_bits,
        }


def default_precision_policy() -> NumericPrecisionPolicy:
    """Return the validated fixed precision policy."""
    return NumericPrecisionPolicy()


def precision_policy_from_options(
    options: Mapping[str, Any] | None,
) -> NumericPrecisionPolicy:
    """Return the fixed precision policy for any workflow options mapping."""
    return NumericPrecisionPolicy()


def cast_floating_array(
    value: ArrayLike,
    *,
    copy: bool = False,
) -> NDArray[Any]:
    """Cast a floating array to the validated Quantas dtype.

    Integer, boolean, string, and object arrays retain their original dtype.
    """
    array = np.asarray(value)
    if np.issubdtype(array.dtype, np.floating):
        return array.astype(FLOAT_DTYPE, copy=copy)
    if np.issubdtype(array.dtype, np.complexfloating):
        return array.astype(COMPLEX_DTYPE, copy=copy)
    return np.array(array, copy=copy) if copy else array


def cast_floating_scalar(value: Any) -> Any:
    """Cast one floating scalar to the validated Quantas dtype."""
    if isinstance(value, (float, np.floating)):
        return FLOAT_DTYPE.type(value)
    if isinstance(value, (complex, np.complexfloating)):
        return COMPLEX_DTYPE.type(value)
    return value


__all__ = [
    "COMPLEX_DTYPE",
    "FLOAT_DTYPE",
    "NumericPrecisionPolicy",
    "cast_floating_array",
    "cast_floating_scalar",
    "default_precision_policy",
    "precision_policy_from_options",
]
