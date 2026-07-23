# -*- coding: utf-8 -*-

"""Primitive value parsing and parameter-name normalization for EOS specs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .models import EOSFitDomain
from .spec import EOSSpecError, _Entry, _Section


def _merge_entries(*sections: _Section | None) -> dict[str, _Entry]:
    merged: dict[str, _Entry] = {}
    for section in sections:
        if section is not None:
            merged.update(section.mapping())
    return merged


def _parse_domain(entry: _Entry, source: Path | None, section: str) -> EOSFitDomain:
    aliases = {
        "pv": EOSFitDomain.PRESSURE_VOLUME,
        "p-v": EOSFitDomain.PRESSURE_VOLUME,
        "vt": EOSFitDomain.VOLUME_TEMPERATURE,
        "v-t": EOSFitDomain.VOLUME_TEMPERATURE,
        "pvt": EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE,
        "p-v-t": EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE,
    }
    value = aliases.get(entry.value.strip().lower())
    if value is None:
        raise EOSSpecError(
            "domain must be pv, vt, or pvt",
            source=source,
            line=entry.line,
            section=section,
        )
    return value


def _canonical_parameter_for_request(
    value: str,
    domain: EOSFitDomain,
    target: str,
) -> str:
    text = str(value).strip()
    key = text.lower().replace("-", "_").replace("'", "p")
    compact = key.replace("_", "")
    if domain is EOSFitDomain.VOLUME_TEMPERATURE:
        aliases = {
            "v0": "V0",
            "l0": "L0",
            "tref": "temperature_ref",
            "temperatureref": "temperature_ref",
            "alpharef": "alpha_ref",
            "thetae": "theta_e",
            "thetad": "theta_d0",
            "thetad0": "theta_d0",
            "debyetemperature": "theta_d0",
            "gamma": "gamma0",
            "gamma0": "gamma0",
            "gruneisen": "gamma0",
            "q": "q",
            "thetasat": "theta_sat",
            "kprime": "kp",
            "kp": "kp",
            "alpha0": "alpha0",
            "alpha1": "alpha1",
            "alpha2": "alpha2",
            "p1": "p1",
        }
        return aliases.get(compact, key)
    if domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
        aliases = {
            "b0": "K0",
            "k0": "K0",
            "kp": "KP",
            "kprime": "KP",
            "kpp": "KPP",
            "kdoubleprime": "KPP",
            "v0": "V0",
            "tref": "temperature_ref",
            "temperatureref": "temperature_ref",
            "alpharef": "alpha_ref",
            "thetae": "theta_e",
            "thetad": "theta_d0",
            "thetad0": "theta_d0",
            "debyetemperature": "theta_d0",
            "gamma": "gamma0",
            "gamma0": "gamma0",
            "gruneisen": "gamma0",
            "q": "q",
            "thetasat": "theta_sat",
            "dk0dt": "dK0_dT",
            "alpha0": "alpha0",
            "alpha1": "alpha1",
            "alpha2": "alpha2",
            "delta": "delta",
            "p1": "p1",
        }
        return aliases.get(compact, key)
    if target == "volume":
        return canonical_parameter_name(text)
    aliases = {
        "m0": "M0",
        "mp": "MP",
        "mpp": "MPP",
        "l0": "L0",
    }
    return aliases.get(compact, text)


def canonical_parameter_name(value: str) -> str:
    """Return the documented canonical EOS parameter name for one alias."""
    text = str(value).strip()
    key = text.lower().replace("_", "").replace("'", "p")
    aliases = {
        "b0": "K0",
        "k0": "K0",
        "kp": "KP",
        "kprime": "KP",
        "kpp": "KPP",
        "kdoubleprime": "KPP",
        "v0": "V0",
        "l0": "L0",
        "m0": "M0",
        "mp": "MP",
        "mpp": "MPP",
        "tref": "temperature_ref",
        "temperatureref": "temperature_ref",
        "thetae": "theta_e",
        "thetad": "theta_d0",
        "thetad0": "theta_d0",
        "debyetemperature": "theta_d0",
        "gamma": "gamma0",
        "gamma0": "gamma0",
        "gruneisen": "gamma0",
        "q": "q",
        "thetasat": "theta_sat",
        "dk0dt": "dK0_dT",
    }
    return aliases.get(key, text)


def _parse_bound(
    entry: _Entry, source: Path | None, section: str
) -> tuple[float | None, float | None]:
    parts = entry.value.split(":")
    if len(parts) != 2:
        raise EOSSpecError(
            "bounds must use LOW : HIGH",
            source=source,
            line=entry.line,
            section=section,
        )
    low = _nullable_float(parts[0], entry, source, section)
    high = _nullable_float(parts[1], entry, source, section)
    if low is not None and high is not None and low >= high:
        raise EOSSpecError(
            "bound lower value must be less than upper value",
            source=source,
            line=entry.line,
            section=section,
        )
    return low, high


def _nullable_float(
    text: str, entry: _Entry, source: Path | None, section: str
) -> float | None:
    stripped = text.strip().lower()
    if stripped in {"", "none", "null", "-inf", "+inf", "inf"}:
        return None
    try:
        value = float(stripped)
    except ValueError as exc:
        raise EOSSpecError(
            f"invalid bound value {text.strip()!r}",
            source=source,
            line=entry.line,
            section=section,
        ) from exc
    if not _is_finite(value):
        return None
    return value


def _finite_float(entry: _Entry, source: Path | None, section: str) -> float:
    try:
        value = float(entry.value)
    except ValueError as exc:
        raise EOSSpecError(
            f"{entry.key} must be a finite number",
            source=source,
            line=entry.line,
            section=section,
        ) from exc
    if not _is_finite(value):
        raise EOSSpecError(
            f"{entry.key} must be a finite number",
            source=source,
            line=entry.line,
            section=section,
        )
    return value


def _is_finite(value: float) -> bool:
    return value == value and value not in {float("inf"), float("-inf")}


def _positive_float(entry: _Entry, source: Path | None, section: str) -> float:
    """Return one strictly positive finite floating-point entry."""
    value = _optional_positive_float(entry, source, section)
    assert value is not None
    return value


def _optional_positive_float(
    entry: _Entry | None,
    source: Path | None,
    section: str,
) -> float | None:
    if entry is None:
        return None
    value = _finite_float(entry, source, section)
    if value <= 0.0:
        raise EOSSpecError(
            f"{entry.key} must be positive",
            source=source,
            line=entry.line,
            section=section,
        )
    return value


def _positive_int(entry: _Entry, source: Path | None, section: str) -> int:
    try:
        value = int(entry.value)
    except ValueError as exc:
        raise EOSSpecError(
            f"{entry.key} must be a positive integer",
            source=source,
            line=entry.line,
            section=section,
        ) from exc
    if value <= 0:
        raise EOSSpecError(
            f"{entry.key} must be a positive integer",
            source=source,
            line=entry.line,
            section=section,
        )
    return value


def _optional_positive_int(
    entry: _Entry | None,
    source: Path | None,
    section: str,
) -> int | None:
    return None if entry is None else _positive_int(entry, source, section)


def _parse_bool_entry(
    entry: _Entry | None,
    default: bool,
    source: Path | None,
    section: str,
) -> bool:
    if entry is None:
        return default
    values = {
        "yes": True,
        "true": True,
        "on": True,
        "1": True,
        "no": False,
        "false": False,
        "off": False,
        "0": False,
    }
    value = values.get(entry.value.strip().lower())
    if value is None:
        raise EOSSpecError(
            f"{entry.key} must be yes or no",
            source=source,
            line=entry.line,
            section=section,
        )
    return value


def _enum_value(
    entry: _Entry,
    values: dict[str, Any],
    source: Path | None,
    name: str,
) -> Any:
    value = values.get(entry.value.strip().lower())
    if value is None:
        raise EOSSpecError(
            f"unknown {name} {entry.value!r}",
            source=source,
            line=entry.line,
        )
    return value


def _entry_error(message: str, section: _Section, source: Path | None) -> EOSSpecError:
    return EOSSpecError(
        message,
        source=source,
        line=section.line,
        section=section.display_name,
    )
