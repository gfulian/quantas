# -*- coding: utf-8 -*-

"""Chemical-formula parsing and molar-mass calculations."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
import re

from .masses import atomic_mass
from .symbols import symbol2number

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?|\d+(?:\.\d+)?|\.\d+|[()\[\]{}])")
_OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}"}
_CLOSE_TO_OPEN = {value: key for key, value in _OPEN_TO_CLOSE.items()}


def _merge_counts(
    target: dict[str, float], source: Mapping[str, float], multiplier: float
) -> None:
    for symbol, count in source.items():
        target[symbol] = target.get(symbol, 0.0) + float(count) * multiplier


class _FormulaParser:
    """Recursive parser for one non-hydrate formula fragment."""

    def __init__(self, formula: str) -> None:
        compact = formula.replace(" ", "")
        if not compact:
            raise ValueError("chemical formula must be a non-empty string")
        tokens = _FORMULA_TOKEN.findall(compact)
        if "".join(tokens) != compact:
            raise ValueError(f"unsupported token in chemical formula: {formula!r}")
        self.tokens = tokens
        self.index = 0

    def parse(self) -> dict[str, float]:
        """Parse all tokens and return element counts."""
        counts = self._parse_group(stop_token=None)
        if self.index != len(self.tokens):
            token = self.tokens[self.index]
            raise ValueError(f"unexpected token in chemical formula: {token!r}")
        if not counts:
            raise ValueError("chemical formula does not contain any element")
        return counts

    def _parse_group(self, stop_token: str | None) -> dict[str, float]:
        counts: dict[str, float] = {}
        while self.index < len(self.tokens):
            token = self.tokens[self.index]
            if stop_token is not None and token == stop_token:
                self.index += 1
                return counts
            if token in _CLOSE_TO_OPEN:
                raise ValueError(f"unmatched closing bracket: {token!r}")
            if token in _OPEN_TO_CLOSE:
                self.index += 1
                subgroup = self._parse_group(_OPEN_TO_CLOSE[token])
                multiplier = self._parse_multiplier()
                _merge_counts(counts, subgroup, multiplier)
                continue
            if token[0].isalpha():
                symbol2number(token)
                self.index += 1
                counts[token] = counts.get(token, 0.0) + self._parse_multiplier()
                continue
            raise ValueError(f"unexpected multiplier in chemical formula: {token!r}")
        if stop_token is not None:
            raise ValueError(f"missing closing bracket: {stop_token!r}")
        return counts

    def _parse_multiplier(self) -> float:
        if self.index >= len(self.tokens):
            return 1.0
        token = self.tokens[self.index]
        try:
            multiplier = float(token)
        except ValueError:
            return 1.0
        if multiplier <= 0.0:
            raise ValueError("formula multipliers must be positive")
        self.index += 1
        return multiplier


def parse_formula(formula: str) -> dict[str, float]:
    """Parse a chemical formula into elemental stoichiometric counts.

    Parameters
    ----------
    formula : str
        Chemical formula. Parentheses, square brackets, braces and middle-dot
        hydrate notation are supported, for example ``"CaMg(CO3)2"`` and
        ``"CuSO4·5H2O"``.

    Returns
    -------
    dict
        Mapping from chemical symbol to stoichiometric count.

    Raises
    ------
    ValueError
        If the formula is empty, malformed, or contains an unknown element.
    """
    if not isinstance(formula, str) or not formula.strip():
        raise ValueError("chemical formula must be a non-empty string")
    totals: defaultdict[str, float] = defaultdict(float)
    for fragment in formula.replace(" ", "").split("·"):
        if not fragment:
            raise ValueError(f"empty formula fragment in {formula!r}")
        match = re.match(r"^(\d+(?:\.\d+)?|\.\d+)(?=[A-Z(\[{])", fragment)
        multiplier = float(match.group(1)) if match else 1.0
        if match:
            fragment = fragment[match.end() :]
        if multiplier <= 0.0:
            raise ValueError("formula multipliers must be positive")
        _merge_counts(totals, _FormulaParser(fragment).parse(), multiplier)
    return dict(totals)


def molar_mass(formula: str | Mapping[str, float]) -> float:
    """Return the molar mass of a chemical formula.

    Parameters
    ----------
    formula : str or mapping
        Formula string or pre-parsed element-count mapping.

    Returns
    -------
    float
        Molar mass in grams per mole.

    Raises
    ------
    ValueError
        If an element is unknown or a count is not positive and finite.
    """
    composition = parse_formula(formula) if isinstance(formula, str) else formula
    total = 0.0
    for symbol, count in composition.items():
        value = float(count)
        if value <= 0.0:
            raise ValueError("formula counts must be positive")
        total += atomic_mass(symbol) * value
    return float(total)
