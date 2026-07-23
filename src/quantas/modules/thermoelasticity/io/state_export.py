# -*- coding: utf-8 -*-

"""Portable single-state stiffness inputs for Elasticity and SEISMIC."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np

from quantas.modules.thermoelasticity.models import ThermoelasticResult


def write_thermoelastic_state_input(
    result: ThermoelasticResult,
    filename: str | Path,
    *,
    tensor_condition: Literal["isothermal", "adiabatic"] = "adiabatic",
    overwrite: bool = False,
) -> Path:
    r"""Write one reconstructed state in the shared elasticity text format.

    The output contains a job title, a complete symmetric ``6 x 6`` stiffness
    matrix in GPa, and density in kg m\ :sup:`-3`.  The same file can be read by
    both :mod:`quantas.modules.elasticity` and :mod:`quantas.modules.seismic`;
    Elasticity ignores the final density line, while SEISMIC requires it.

    Parameters
    ----------
    result : ThermoelasticResult
        Reconstructed result containing exactly one pressure-temperature state.
    filename : str or Path
        Destination text file.
    tensor_condition : {"isothermal", "adiabatic"}, optional
        Thermodynamic stiffness condition written to the file.
    overwrite : bool, optional
        Replace an existing destination.

    Returns
    -------
    Path
        Written path.

    Raises
    ------
    ValueError
        If the result does not contain exactly one valid state or the requested
        tensor condition is unavailable.
    FileExistsError
        If the destination exists and ``overwrite`` is false.
    """
    if result.temperature.size != 1 or result.pressure.size != 1:
        raise ValueError("single-state input requires exactly one P-T state")
    if tensor_condition == "isothermal":
        stiffness = result.stiffness_isothermal
    elif tensor_condition == "adiabatic":
        stiffness = result.stiffness_adiabatic
        if result.adiabatic_valid_mask is not None and not bool(
            result.adiabatic_valid_mask[0, 0]
        ):
            raise ValueError("adiabatic stiffness is invalid at the requested state")
    else:
        raise ValueError("tensor_condition must be isothermal or adiabatic")
    if stiffness is None:
        raise ValueError(f"{tensor_condition} stiffness is unavailable")
    density = float(result.density[0, 0])
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("density is unavailable at the requested state")
    matrix = np.asarray(stiffness[0, 0], dtype=np.float64)
    if matrix.shape != (6, 6) or not np.all(np.isfinite(matrix)):
        raise ValueError("stiffness must be a finite 6 x 6 matrix")
    path = Path(filename)
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    condition_symbol = "C^T" if tensor_condition == "isothermal" else "C^S"
    title = (
        f"{result.jobname} | P={float(result.pressure[0]):g} GPa, "
        f"T={float(result.temperature[0]):g} K, {condition_symbol}"
    )
    lines = [title]
    matrix = 0.5 * (matrix + matrix.T)
    for row in range(6):
        lines.append(" ".join(f"{value:14.6f}" for value in matrix[row]))
    lines.append(f"{density:.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


__all__ = ["write_thermoelastic_state_input"]
