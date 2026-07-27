# -*- coding: utf-8 -*-

"""Click helpers for optional elastic-tensor component transformations."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any, TypeVar

import click
import numpy as np

from quantas.cli.contracts import SCIENTIFIC_GROUP
from quantas.cli.grouped_options import grouped_option
from quantas.api.common import TensorRotation


F = TypeVar("F", bound=Callable[..., Any])


def tensor_rotation_options(function: F) -> F:
    """Add shared tensor-rotation options to a Click command.

    Parameters
    ----------
    function : callable
        Click command callback.

    Returns
    -------
    callable
        Decorated callback with ``--rotate-xyz`` and ``--rotation-matrix``.
    """
    function = grouped_option(
        "--rotation-matrix",
        group=SCIENTIFIC_GROUP,
        type=click.Path(exists=True, dir_okay=False, path_type=Path),
        default=None,
        help=(
            "Text file containing a proper 3 x 3 source-to-analysis "
            "component-transformation matrix."
        ),
    )(function)
    function = grouped_option(
        "--rotate-xyz",
        group=SCIENTIFIC_GROUP,
        type=float,
        nargs=3,
        default=None,
        metavar="X Y Z",
        help=(
            "Right-handed fixed-axis rotations in degrees, applied about "
            "source X, then Y, then Z."
        ),
    )(function)
    return function


def resolve_tensor_rotation(
    rotate_xyz: tuple[float, float, float] | None,
    rotation_matrix: Path | None,
) -> TensorRotation | None:
    """Build the requested tensor transformation from CLI values.

    Parameters
    ----------
    rotate_xyz : tuple of float or None
        Fixed-axis XYZ angles in degrees.
    rotation_matrix : Path or None
        Text file containing nine matrix elements.

    Returns
    -------
    TensorRotation or None
        Validated transformation, or ``None`` when no rotation was requested.

    Raises
    ------
    click.UsageError
        If both transformation modes are requested.
    click.BadParameter
        If the matrix file is malformed or the transformation is invalid.
    """
    if rotate_xyz is not None and rotation_matrix is not None:
        raise click.UsageError(
            "Use either --rotate-xyz or --rotation-matrix, not both."
        )
    try:
        if rotate_xyz is not None:
            return TensorRotation.from_xyz(*rotate_xyz, degrees=True)
        if rotation_matrix is not None:
            return TensorRotation.from_matrix(
                read_rotation_matrix(rotation_matrix),
                description=f"Loaded from {rotation_matrix}",
            )
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="tensor rotation") from exc
    return None


def read_rotation_matrix(filename: str | Path) -> np.ndarray:
    """Read a whitespace- or comma-separated ``3 x 3`` matrix.

    Parameters
    ----------
    filename : str or Path
        Matrix text file.

    Returns
    -------
    numpy.ndarray
        Floating-point matrix with shape ``(3, 3)``.

    Raises
    ------
    ValueError
        If the file does not contain exactly nine finite values.
    """
    path = Path(filename)
    text = path.read_text(encoding="utf-8").replace(",", " ")
    try:
        values = np.asarray([float(token) for token in text.split()], dtype=float)
    except ValueError as exc:
        raise ValueError(
            "The rotation-matrix file must contain numeric values."
        ) from exc
    if values.size != 9 or not np.all(np.isfinite(values)):
        raise ValueError(
            "The rotation-matrix file must contain exactly nine finite values."
        )
    return values.reshape(3, 3)
