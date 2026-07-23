# -*- coding: utf-8 -*-

"""Fit, normalization, solver, and constraint resolution for EOS specs."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

from quantas.core.math.fitting import (
    CovarianceScaling,
    EffectiveVarianceOptions,
    ODRDifferenceScheme,
    OLSOptions,
    OrthogonalDistanceOptions,
    WLSOptions,
)
from quantas.core.physics.eos import MGDNormalization, MGDVolumeBasis

from .domains.pvt import pvt_parameter_names
from .domains.vt import temperature_parameter_names
from .models import EOSFitDomain, ParameterConstraint
from .spec import EOSSpecError, _Entry
from ._spec_values import (
    _canonical_parameter_for_request,
    _finite_float,
    _optional_positive_float,
    _optional_positive_int,
    _parse_bound,
    _positive_float,
)

_PARAMETER_PREFIXES = ("fix.", "initial.", "bound.")
_MGD_NORMALIZATION_KEYS = {
    "volume_basis",
    "atoms_per_cell",
    "atoms_per_formula_unit",
    "formula",
    "formula_units_per_cell",
}


def _build_mgd_normalization(
    entries: dict[str, _Entry],
    *,
    source: Path | None,
    section: str,
) -> MGDNormalization:
    """Build MGD atom-count normalization from declarative entries."""
    basis_entry = entries.get("volume_basis")
    basis_text = "cell" if basis_entry is None else basis_entry.value.strip().lower()
    try:
        basis = MGDVolumeBasis(basis_text)
    except ValueError as exc:
        raise EOSSpecError(
            "volume_basis must be 'cell' or 'molar-formula-unit'",
            source=source,
            line=None if basis_entry is None else basis_entry.line,
            section=section,
        ) from exc

    formula_entry = entries.get("formula")
    formula = None if formula_entry is None else formula_entry.value.strip()
    try:
        if basis is MGDVolumeBasis.CELL:
            atoms_entry = entries.get("atoms_per_cell")
            z_entry = entries.get("formula_units_per_cell")
            if "atoms_per_formula_unit" in entries:
                entry = entries["atoms_per_formula_unit"]
                raise EOSSpecError(
                    "atoms_per_formula_unit is valid only for molar-formula-unit normalization",
                    source=source,
                    line=entry.line,
                    section=section,
                )
            atoms = (
                None
                if atoms_entry is None
                else _positive_float(atoms_entry, source, section)
            )
            z = None if z_entry is None else _positive_float(z_entry, source, section)
            return MGDNormalization.cell(
                atoms_per_cell=atoms,
                formula=formula,
                formula_units_per_cell=z,
            )
        if "atoms_per_cell" in entries or "formula_units_per_cell" in entries:
            name = (
                "atoms_per_cell"
                if "atoms_per_cell" in entries
                else "formula_units_per_cell"
            )
            entry = entries[name]
            raise EOSSpecError(
                f"{name} is valid only for cell normalization",
                source=source,
                line=entry.line,
                section=section,
            )
        atoms_entry = entries.get("atoms_per_formula_unit")
        atoms = (
            None
            if atoms_entry is None
            else _positive_float(atoms_entry, source, section)
        )
        return MGDNormalization.molar_formula_unit(
            atoms_per_formula_unit=atoms,
            formula=formula,
        )
    except EOSSpecError:
        raise
    except (TypeError, ValueError) as exc:
        candidates = [
            entries.get(name)
            for name in (
                "volume_basis",
                "atoms_per_cell",
                "atoms_per_formula_unit",
                "formula",
                "formula_units_per_cell",
            )
        ]
        location = next((entry for entry in candidates if entry is not None), None)
        raise EOSSpecError(
            str(exc),
            source=source,
            line=None if location is None else location.line,
            section=section,
        ) from exc


def _reject_mgd_normalization_keys(
    entries: dict[str, _Entry],
    *,
    source: Path | None,
    section: str,
) -> None:
    """Reject MGD-only normalization keys for a non-MGD model."""
    name = next((item for item in _MGD_NORMALIZATION_KEYS if item in entries), None)
    if name is None:
        return
    entry = entries[name]
    raise EOSSpecError(
        f"{name} is valid only for Mie-Gruneisen-Debye thermal pressure",
        source=source,
        line=entry.line,
        section=section,
    )


def _build_solver_options(
    entries: dict[str, _Entry],
    source: Path | None,
    section: str,
) -> Any:
    solver_entry = entries.get("solver")
    solver = (
        "ols"
        if solver_entry is None
        else solver_entry.value.strip().lower().replace("_", "-")
    )
    aliases = {
        "ols": "ols",
        "wls": "wls",
        "effective-variance": "effective-variance",
        "effectivevariance": "effective-variance",
        "ev": "effective-variance",
        "odr": "odr",
    }
    method = aliases.get(solver)
    if method is None:
        line = None if solver_entry is None else solver_entry.line
        raise EOSSpecError(
            f"unknown solver {solver!r}; expected ols, wls, effective-variance, or odr",
            source=source,
            line=line,
            section=section,
        )
    covariance = None
    covariance_entry = entries.get("covariance_scaling")
    if covariance_entry is not None:
        key = covariance_entry.value.strip().lower().replace("-", "_")
        try:
            covariance = CovarianceScaling(key)
        except ValueError as exc:
            raise EOSSpecError(
                "covariance_scaling must be absolute, reduced-chi-square, or inflate-only",
                source=source,
                line=covariance_entry.line,
                section=section,
            ) from exc
    max_iterations = _optional_positive_int(
        entries.get("max_iterations"), source, section
    )
    ftol = _optional_positive_float(entries.get("ftol"), source, section)
    xtol = _optional_positive_float(entries.get("xtol"), source, section)
    gtol = _optional_positive_float(entries.get("gtol"), source, section)
    if method != "effective-variance" and "inner_max_iterations" in entries:
        entry = entries["inner_max_iterations"]
        raise EOSSpecError(
            "inner_max_iterations is valid only for solver = effective-variance",
            source=source,
            line=entry.line,
            section=section,
        )
    if method != "odr":
        for name in ("odr_difference", "odr_ndigit"):
            if name in entries:
                entry = entries[name]
                raise EOSSpecError(
                    f"{name} is valid only for solver = odr",
                    source=source,
                    line=entry.line,
                    section=section,
                )
    try:
        if method == "ols":
            return OLSOptions(
                covariance_scaling=covariance,
                max_iterations=max_iterations,
                ftol=ftol,
                xtol=xtol,
                gtol=gtol,
            )
        if method == "wls":
            return WLSOptions(
                covariance_scaling=covariance,
                max_iterations=max_iterations,
                ftol=ftol,
                xtol=xtol,
                gtol=gtol,
            )
        if method == "effective-variance":
            return EffectiveVarianceOptions(
                covariance_scaling=covariance,
                max_iterations=max_iterations,
                ftol=ftol,
                xtol=xtol,
                gtol=gtol,
                inner_max_iterations=_optional_positive_int(
                    entries.get("inner_max_iterations"), source, section
                ),
            )
        difference_entry = entries.get("odr_difference")
        difference = ODRDifferenceScheme.CENTRAL
        if difference_entry is not None:
            try:
                difference = ODRDifferenceScheme(difference_entry.value.strip().lower())
            except ValueError as exc:
                raise EOSSpecError(
                    "odr_difference must be central or forward",
                    source=source,
                    line=difference_entry.line,
                    section=section,
                ) from exc
        return OrthogonalDistanceOptions(
            covariance_scaling=covariance,
            max_iterations=max_iterations,
            ftol=ftol,
            xtol=xtol,
            gtol=gtol,
            difference_scheme=difference,
            ndigit=_optional_positive_int(entries.get("odr_ndigit"), source, section),
        )
    except EOSSpecError:
        raise
    except (TypeError, ValueError) as exc:
        raise EOSSpecError(str(exc), source=source, section=section) from exc


def _build_constraints(
    entries: dict[str, _Entry],
    source: Path | None,
    section: str,
    *,
    domain: EOSFitDomain,
    target: str,
    model: Any,
) -> tuple[ParameterConstraint, ...]:
    allowed = _allowed_parameter_names(domain, target, model)
    parameters: dict[str, dict[str, _Entry]] = {}
    for key, entry in entries.items():
        lower_key = key.lower()
        if not lower_key.startswith(_PARAMETER_PREFIXES):
            continue
        kind, raw_name = key.split(".", 1)
        name = _canonical_parameter_for_request(raw_name, domain, target)
        if name not in allowed:
            raise EOSSpecError(
                f"parameter {raw_name!r} is not available for {domain.value}/{target}; "
                f"available parameters are {', '.join(allowed)}",
                source=source,
                line=entry.line,
                section=section,
            )
        state = parameters.setdefault(name, {})
        if kind in state:
            raise EOSSpecError(
                f"duplicate {kind} declaration for parameter {name!r}",
                source=source,
                line=entry.line,
                section=section,
            )
        state[kind] = entry
    constraints: list[ParameterConstraint] = []
    for name, state in parameters.items():
        if "fix" in state and ("initial" in state or "bound" in state):
            conflict = state.get("initial") or state.get("bound")
            assert conflict is not None
            raise EOSSpecError(
                f"fixed parameter {name!r} cannot also declare initial or bound",
                source=source,
                line=conflict.line,
                section=section,
            )
        if "fix" in state:
            entry = state["fix"]
            constraints.append(
                cast(
                    ParameterConstraint,
                    ParameterConstraint.fixed(
                        name, _finite_float(entry, source, section)
                    ),
                )
            )
            continue
        if "bound" in state and "initial" not in state:
            entry = state["bound"]
            raise EOSSpecError(
                f"bound for free parameter {name!r} requires initial.{name} = VALUE",
                source=source,
                line=entry.line,
                section=section,
            )
        if "initial" in state:
            initial = _finite_float(state["initial"], source, section)
            lower_bound, upper_bound = float("-inf"), float("inf")
            if "bound" in state:
                lower_value, upper_value = _parse_bound(state["bound"], source, section)
                if lower_value is not None:
                    lower_bound = lower_value
                if upper_value is not None:
                    upper_bound = upper_value
            try:
                constraints.append(
                    cast(
                        ParameterConstraint,
                        ParameterConstraint.free(
                            name,
                            initial,
                            lower_bound=lower_bound,
                            upper_bound=upper_bound,
                        ),
                    )
                )
            except ValueError as exc:
                raise EOSSpecError(
                    str(exc), source=source, line=state["initial"].line, section=section
                ) from exc
    return tuple(constraints)


def _allowed_parameter_names(
    domain: EOSFitDomain, target: str, model: Any
) -> tuple[str, ...]:
    if domain is EOSFitDomain.PRESSURE_VOLUME:
        return (
            ("K0", "KP", "KPP", "V0")
            if target == "volume"
            else ("M0", "MP", "MPP", "L0")
        )
    if domain is EOSFitDomain.VOLUME_TEMPERATURE:
        return temperature_parameter_names(
            model,
            axial=target in {"a", "b", "c"},
        )
    return pvt_parameter_names(model)
