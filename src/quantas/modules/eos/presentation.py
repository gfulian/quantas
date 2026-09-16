# -*- coding: utf-8 -*-

"""Human-readable names and units for EOS reports and plots.

Machine-readable identifiers remain unchanged in requests, HDF5, and CSV
metadata.  This module centralizes presentation-only transformations so every
frontend renders the same scientific terminology.
"""

from __future__ import annotations

import re
from typing import Any

from quantas.core.math.fitting import FitMethod
from quantas.core.physics.eos import (
    EOSModel,
    PVTModel,
    TemperatureEOSFamily,
    TemperatureEOSModel,
    TemperatureEOSVariant,
    ThermalPressureFamily,
)

from .models import EOSFitDomain


_SUPERSCRIPT_TRANSLATION = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")

_DOMAIN_LABELS = {
    EOSFitDomain.PRESSURE_VOLUME: "Pressure-volume",
    EOSFitDomain.ENERGY_VOLUME: "Energy-volume",
    EOSFitDomain.VOLUME_TEMPERATURE: "Volume-temperature",
    EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE: "Pressure-volume-temperature",
}

_TARGET_LABELS = {
    "pressure": "Pressure",
    "temperature": "Temperature",
    "volume": "Volume",
    "energy": "Energy",
    "a": "Cell parameter a",
    "b": "Cell parameter b",
    "c": "Cell parameter c",
    "residual": "Residual",
    "standardized_residual": "Standardized residual",
    "finite_strain": "Finite strain",
    "normalized_pressure": "Normalized pressure",
}

_PARAMETER_LABELS = {
    "E0": "E0",
    "K0": "K0",
    "KP": "K'",
    "KPP": "K''",
    "V0": "V0",
    "M0": "M0",
    "MP": "M'",
    "MPP": "M''",
    "L0": "L0",
    "temperature_ref": "Tref",
    "alpha0": "alpha0",
    "alpha1": "alpha1",
    "alpha2": "alpha2",
    "alpha_ref": "alpha_ref",
    "p1": "p1",
    "theta_sat": "theta_sat",
    "theta_e": "theta_E",
    "theta_d0": "theta_D0",
    "gamma0": "gamma0",
    "q": "q",
    "dK0_dT": "dK0/dT",
    "delta": "delta",
    "kp": "K'",
    "a0": "a0",
    "b0": "b0",
    "c0": "c0",
    "eta_a": "eta_a",
    "eta_b": "eta_b",
    "eta_c": "eta_c",
    "M_a": "M_a",
    "M_b": "M_b",
    "M_c": "M_c",
}

_SOLVER_LABELS = {
    FitMethod.OLS.value: "Ordinary least squares",
    FitMethod.WLS.value: "Weighted least squares",
    FitMethod.EFFECTIVE_VARIANCE.value: "Effective variance",
    FitMethod.ODR.value: "Orthogonal distance regression",
}

_TEMPERATURE_FAMILY_LABELS = {
    TemperatureEOSFamily.BERMAN: "Berman",
    TemperatureEOSFamily.FEI: "Fei",
    TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL: "Modified Holland-Powell",
    TemperatureEOSFamily.SALJE: "Salje",
    TemperatureEOSFamily.KROLL_HOLLAND_POWELL: "Kroll-Holland-Powell",
}

_VARIANT_LABELS = {
    TemperatureEOSVariant.LINEAR: "linear",
    TemperatureEOSVariant.QUADRATIC: "quadratic",
    TemperatureEOSVariant.INVERSE_SQUARE: "inverse-square",
    TemperatureEOSVariant.SIMPLIFIED: "simplified",
    TemperatureEOSVariant.GENERAL: "general",
    TemperatureEOSVariant.STANDARD: "standard",
}


def format_unit(unit: str | None) -> str | None:
    """Return a compact typographic representation of one EOS unit.

    This formatter is intended primarily for figure labels. Text reports and
    terminal output use :func:`format_text_unit` to avoid Unicode-dependent
    scientific notation.

    Parameters
    ----------
    unit : str | None
        Unit symbol or canonical unit name.

    Returns
    -------
    str | None
        A compact typographic representation of one EOS unit.
    """
    if unit is None:
        return None
    text = str(unit).strip()
    if text.lower() in {"", "1", "dimensionless", "none"}:
        return None
    text = text.replace("angstrom", "Å").replace("Angstrom", "Å")
    text = text.replace("cm^3/mol", "cm³ mol⁻¹")
    text = text.replace("J/mol", "J mol⁻¹")
    text = text.replace("kJ/mol", "kJ mol⁻¹")
    text = re.sub(
        r"\^\(?([+-]?\d+)\)?",
        lambda match: match.group(1).translate(_SUPERSCRIPT_TRANSLATION),
        text,
    )
    text = text.replace("/K", " K⁻¹")
    return text


def format_text_unit(unit: str | None) -> str | None:
    """Return a portable ASCII representation of one EOS unit.

    Machine-readable units are left unchanged elsewhere. This function is
    strictly a text-presentation boundary for terminal and plain-text reports.

    Parameters
    ----------
    unit : str | None
        Unit symbol or canonical unit name.

    Returns
    -------
    str | None
        A portable ASCII representation of one EOS unit.
    """
    if unit is None:
        return None
    text = str(unit).strip()
    if text.lower() in {"", "1", "dimensionless", "none"}:
        return None
    replacements = (
        ("Å", "angstrom"),
        ("Angstrom", "angstrom"),
        ("cm³ mol⁻¹", "cm^3/mol"),
        ("cm^3 mol^-1", "cm^3/mol"),
        ("J mol⁻¹", "J/mol"),
        ("J mol^-1", "J/mol"),
        ("kJ mol⁻¹", "kJ/mol"),
        ("kJ mol^-1", "kJ/mol"),
        ("⁻¹/²", "^-1/2"),
        ("⁻¹", "^-1"),
        ("²", "^2"),
        ("³", "^3"),
        (" K^-1", "/K"),
        (" K⁻¹", "/K"),
    )
    for source, target in replacements:
        text = text.replace(source, target)
    text = re.sub(r"^\(([^()]+)\)\^([+-]?\d+)$", r"\1^\2", text)
    return text


def domain_label(domain: EOSFitDomain | str) -> str:
    """Return a human-readable scientific-domain name.

    Parameters
    ----------
    domain : EOSFitDomain | str
        EOS domain identifier.

    Returns
    -------
    str
        A human-readable scientific-domain name.
    """
    return _DOMAIN_LABELS[EOSFitDomain(domain)]


def target_label(name: str) -> str:
    """Return a human-readable property or target name.

    Parameters
    ----------
    name : str
        Stable name or identifier.

    Returns
    -------
    str
        A human-readable property or target name.
    """
    normalized = str(name).strip()
    return _TARGET_LABELS.get(normalized, normalized.replace("_", " ").title())


def parameter_label(name: str) -> str:
    """Return a compact scientific symbol for one EOS parameter.

    Parameters
    ----------
    name : str
        Stable name or identifier.

    Returns
    -------
    str
        A compact scientific symbol for one EOS parameter.
    """
    normalized = str(name).strip()
    return _PARAMETER_LABELS.get(normalized, normalized.replace("_", " "))


def solver_label(method: Any) -> str:
    """Return the human-readable name of one numerical solver.

    Parameters
    ----------
    method : Any
        Stable numerical or scientific method identifier.

    Returns
    -------
    str
        The human-readable name of one numerical solver.
    """
    value = getattr(method, "value", method)
    normalized = str(value)
    return _SOLVER_LABELS.get(normalized, normalized.replace("_", " ").title())


def model_label(model: Any) -> str:
    """Return a human-readable label for an EOS model specification.

    Parameters
    ----------
    model : Any
        EOS or thermoelastic model used by the operation.

    Returns
    -------
    str
        A human-readable label for an EOS model specification.
    """
    if isinstance(model, EOSModel):
        return f"{model.name} [{model.tag}]"
    if isinstance(model, TemperatureEOSModel):
        temperature_variant = model.variant
        assert temperature_variant is not None
        return (
            f"{_TEMPERATURE_FAMILY_LABELS[model.family]}, "
            f"{_VARIANT_LABELS[temperature_variant]} [{model.tag}]"
        )
    if isinstance(model, PVTModel):
        pressure = model.pressure_spec.name
        coupling = model.coupling_family.value.replace("-", " ")
        if model.thermal_spec is not None:
            thermal = model_label(model.thermal_spec)
        else:
            thermal_pressure = model.thermal_pressure_spec
            assert thermal_pressure is not None
            if (
                thermal_pressure.family_name
                is ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN
            ):
                thermal = "Holland-Powell Einstein thermal pressure"
            else:
                mgd_variant = thermal_pressure.mgd_variant
                suffix = (
                    ""
                    if mgd_variant is None
                    else f", {mgd_variant.value.replace('-', ' ')}"
                )
                thermal = f"Mie-Gruneisen-Debye{suffix}"
        return f"{pressure} + {thermal} ({coupling}) [{model.tag}]"
    tag = getattr(model, "tag", None)
    if tag is not None:
        return str(tag)
    return str(model).replace("_", " ")


def property_label(name: str, unit: str | None = None) -> str:
    """Return a clean property label with an optional parenthesized unit.

    Parameters
    ----------
    name : str
        Stable name or identifier.
    unit : str | None
        Unit symbol or canonical unit name.

    Returns
    -------
    str
        A clean property label with an optional parenthesized unit.
    """
    label = target_label(name)
    rendered = format_unit(unit)
    return label if rendered is None else f"{label} ({rendered})"


__all__ = [
    "domain_label",
    "format_text_unit",
    "format_unit",
    "model_label",
    "parameter_label",
    "property_label",
    "solver_label",
    "target_label",
]
