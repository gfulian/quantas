# -*- coding: utf-8 -*-

"""Post-fit pressure-temperature reconstruction from thermoelastic archives.

This module evaluates the already fitted quasi-static elastic model without
repeating either the static EOS fit or the independent-component fits.  The
QHA equilibrium-volume surface archived by the fitting run supplies
``V(P, T)`` through rectilinear interpolation or controlled extrapolation.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from quantas.modules.thermoelasticity.models import (
    ThermoelasticResult,
)


def _normalized_cell_mass(source: ThermoelasticResult) -> float:
    value = source.metadata.get("normalized_cell_mass_kg")
    if value is not None and np.isfinite(float(value)) and float(value) > 0.0:
        return float(value)
    masses = source.density * source.equilibrium_volume * 1.0e-30
    mass = float(np.nanmedian(masses))
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("normalized cell mass is unavailable in the archive")
    return mass


def _elastic_volume_bounds(source: ThermoelasticResult) -> tuple[float, float]:
    try:
        lower = float(source.metadata["elastic_volume_min_A3"])
        upper = float(source.metadata["elastic_volume_max_A3"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            "elastic-volume bounds are unavailable in the archive"
        ) from exc
    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        raise ValueError("elastic-volume bounds in the archive are invalid")
    return lower, upper


def _required_metadata_text(source: ThermoelasticResult, key: str) -> str:
    value = source.metadata.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"required thermoelastic metadata '{key}' is unavailable")
    return value


def _enforce_extrapolation_policy(
    mask: NDArray[np.bool_],
    policy: str,
    context: str,
) -> None:
    count = int(np.count_nonzero(mask))
    if count == 0 or policy in {"warn", "allow"}:
        return
    if policy == "fail":
        raise ValueError(
            f"{count} of {mask.size} states in {context} require extrapolation"
        )
    raise ValueError(f"invalid extrapolation policy: {policy}")


def _enforce_stability_policy(
    mask: NDArray[np.bool_],
    policy: str,
    context: str,
) -> None:
    """Apply a fail-only mechanical-stability policy during post-fit analysis."""
    count = int(np.count_nonzero(mask))
    if count == 0 or policy in {"warn", "allow"}:
        return
    if policy == "fail":
        raise ValueError(
            f"{count} of {mask.size} states in {context} are mechanically "
            "unstable or indeterminate"
        )
    raise ValueError(f"invalid stability policy: {policy}")


def _stability_summary(stability: Any) -> dict[str, Any]:
    """Return compact serializable mechanical-stability diagnostics."""
    minimum = np.asarray(stability.minimum_eigenvalue, dtype=np.float64)
    finite = minimum[np.isfinite(minimum)]
    return {
        "criterion": stability.criterion,
        "stable_points": int(np.count_nonzero(stability.stable_mask)),
        "unstable_points": int(np.count_nonzero(stability.unstable_mask)),
        "indeterminate_points": int(np.count_nonzero(stability.indeterminate_mask)),
        "minimum_eigenvalue_GPa": float(np.min(finite)) if finite.size else None,
        "tolerance_GPa": float(stability.tolerance),
    }


def _append_stability_warning(
    warnings: list[str],
    stability: Any,
    context: str,
) -> None:
    """Append one readable warning for an unstable or indeterminate field."""
    if stability is None:
        return
    unstable = int(np.count_nonzero(stability.unstable_mask))
    indeterminate = int(np.count_nonzero(stability.indeterminate_mask))
    if unstable or indeterminate:
        warnings.append(
            f"{context} contains {unstable} mechanically unstable and "
            f"{indeterminate} indeterminate Wallace stiffness states"
        )


__all__ = []
