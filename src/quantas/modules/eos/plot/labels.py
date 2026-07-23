# -*- coding: utf-8 -*-

"""Compatibility labels for EOS plot builders."""

from __future__ import annotations

from quantas.modules.eos.presentation import format_unit, property_label


def normalized_pressure_labels(
    metadata: dict[str, object], unit: str | None
) -> tuple[str, str]:
    """Return finite-strain and normalized-pressure axis labels."""
    strain_symbol = str(metadata.get("strain_symbol", "f"))
    pressure_symbol = str(metadata.get("normalized_pressure_symbol", "F"))
    x_label = rf"Finite strain, ${strain_symbol}$"
    unit_text = format_unit(unit)
    y_label = rf"Normalized pressure, ${pressure_symbol}$"
    if unit_text is not None:
        y_label += f" ({unit_text})"
    return x_label, y_label


__all__ = ["format_unit", "normalized_pressure_labels", "property_label"]
