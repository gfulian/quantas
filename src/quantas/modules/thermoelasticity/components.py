# -*- coding: utf-8 -*-

"""Canonical Voigt-component semantics for thermoelastic analysis and I/O.

Plot-specific markers and line styles remain in
:mod:`quantas.modules.thermoelasticity.plot.components`.
"""

from __future__ import annotations

from typing import Literal, TypeAlias

import numpy as np

from quantas.modules.thermoelasticity.models import ThermoelasticResult

ThermoelasticComponentGroup: TypeAlias = Literal[
    "independent",
    "normal",
    "shear",
    "coupling",
    "offdiagonal",
    "all",
]

VOIGT_COMPONENTS: tuple[str, ...] = (
    "C11",
    "C12",
    "C13",
    "C14",
    "C15",
    "C16",
    "C22",
    "C23",
    "C24",
    "C25",
    "C26",
    "C33",
    "C34",
    "C35",
    "C36",
    "C44",
    "C45",
    "C46",
    "C55",
    "C56",
    "C66",
)


def normalize_component_label(label: str) -> str:
    """Return a canonical upper-triangular Voigt stiffness label.

    Parameters
    ----------
    label : str
        Component name, case-insensitive, optionally without the leading
        ``C``.

    Returns
    -------
    str
        Canonical label from :data:`VOIGT_COMPONENTS`.

    Raises
    ------
    ValueError
        If the label does not identify a symmetric ``6 x 6`` Voigt entry.
    """
    text = str(label).strip().upper().replace("_", "")
    if not text.startswith("C"):
        text = f"C{text}"
    if len(text) != 3 or not text[1:].isdigit():
        raise ValueError(f"invalid elastic component label: {label}")
    first = int(text[1])
    second = int(text[2])
    if first < 1 or first > 6 or second < 1 or second > 6:
        raise ValueError(f"invalid elastic component label: {label}")
    if second < first:
        first, second = second, first
    canonical = f"C{first}{second}"
    if canonical not in VOIGT_COMPONENTS:
        raise ValueError(f"invalid elastic component label: {label}")
    return canonical


def component_indices(label: str) -> tuple[int, int]:
    """Return zero-based Voigt indices for one component label.

    Parameters
    ----------
    label : str
        Elastic component label.

    Returns
    -------
    tuple of int
        Zero-based row and column indices.
    """
    canonical = normalize_component_label(label)
    return int(canonical[1]) - 1, int(canonical[2]) - 1


def component_symbol(label: str) -> str:
    """Return a math-text representation of one elastic component.

    Parameters
    ----------
    label : str
        Elastic component label.

    Returns
    -------
    str
        ``C_IJ`` math-text label.
    """
    canonical = normalize_component_label(label)
    return rf"$C_{{{canonical[1:]}}}$"


def resolve_components(
    result: ThermoelasticResult,
    components: tuple[str, ...] | list[str] | None = None,
    *,
    group: ThermoelasticComponentGroup = "independent",
    nonzero_tolerance: float = 1.0e-12,
) -> tuple[str, ...]:
    """Resolve explicit components or one semantic component group.

    Parameters
    ----------
    result : ThermoelasticResult
        Thermoelastic result containing component metadata and, optionally,
        reconstructed tensors.
    components : sequence of str or None, optional
        Explicit component labels.  When supplied, ``group`` is ignored.
    group : ThermoelasticComponentGroup, optional
        Semantic component group.
    nonzero_tolerance : float, optional
        Absolute GPa threshold used when full reconstructed tensors are
        available.

    Returns
    -------
    tuple of str
        Ordered canonical component labels.

    Raises
    ------
    ValueError
        If a requested component or group is unavailable.
    """
    available = _available_components(result, tolerance=nonzero_tolerance)
    if components:
        selected = tuple(normalize_component_label(item) for item in components)
        missing = [item for item in selected if item not in available]
        if missing:
            raise ValueError(
                "requested thermoelastic components are unavailable: "
                + ", ".join(missing)
            )
        return _deduplicate(selected)
    if group == "independent":
        return tuple(label for label in result.independent_labels if label in available)
    candidates: tuple[str, ...]
    if group == "normal":
        candidates = ("C11", "C22", "C33")
    elif group == "shear":
        candidates = ("C44", "C55", "C66")
    elif group == "coupling":
        candidates = ("C12", "C13", "C23")
    elif group == "offdiagonal":
        candidates = tuple(label for label in VOIGT_COMPONENTS if label[1] != label[2])
    elif group == "all":
        candidates = VOIGT_COMPONENTS
    else:
        raise ValueError(f"unknown thermoelastic component group: {group}")
    selected = tuple(label for label in candidates if label in available)
    if not selected:
        raise ValueError(f"component group '{group}' is empty for this result")
    return selected


def _available_components(result: ThermoelasticResult, *, tolerance: float) -> set[str]:
    available = set(result.independent_labels)
    tensor = result.stiffness_isothermal
    if tensor is None:
        return available
    for label in VOIGT_COMPONENTS:
        first, second = component_indices(label)
        values = np.asarray(tensor[..., first, second], dtype=np.float64)
        if np.any(np.isfinite(values) & (np.abs(values) > tolerance)):
            available.add(label)
    return available


def _deduplicate(labels: tuple[str, ...]) -> tuple[str, ...]:
    result: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label not in seen:
            result.append(label)
            seen.add(label)
    return tuple(result)


__all__ = [
    "ThermoelasticComponentGroup",
    "VOIGT_COMPONENTS",
    "component_indices",
    "component_symbol",
    "normalize_component_label",
    "resolve_components",
]
