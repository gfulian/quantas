"""Tests for chemical symbols and atomic-number conversion."""

from __future__ import annotations

import pytest

from quantas.core.chemistry.symbols import number2symbol, symbol2number


def test_number2symbol_returns_known_element_symbol():
    assert number2symbol(8) == "O"


def test_symbol2number_accepts_standard_symbol():
    assert symbol2number("Mg") == 12


def test_symbol_conversion_rejects_unknown_values():
    with pytest.raises(ValueError):
        number2symbol(0)
    with pytest.raises(ValueError):
        symbol2number("Xx")
