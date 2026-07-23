# -*- coding: utf-8 -*-

"""Chemical-symbol and atomic-number conversion utilities."""

from __future__ import annotations

CHEMICAL_SYMBOLS: tuple[str, ...] = (
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)

ATOMIC_NUMBERS: dict[str, int] = {
    symbol: number for number, symbol in enumerate(CHEMICAL_SYMBOLS) if symbol
}


def number2symbol(atomic_number: int) -> str:
    """Return the chemical symbol associated with an atomic number.

    Parameters
    ----------
    atomic_number : int
        Atomic number in the range 1--118.

    Returns
    -------
    str
        Chemical symbol.

    Raises
    ------
    ValueError
        If ``atomic_number`` is not an integer in the supported range.
    """
    if isinstance(atomic_number, bool) or not isinstance(atomic_number, int):
        raise ValueError("atomic_number must be an integer")
    if atomic_number <= 0 or atomic_number >= len(CHEMICAL_SYMBOLS):
        raise ValueError("atomic_number must be between 1 and 118")
    return CHEMICAL_SYMBOLS[atomic_number]


def symbol2number(symbol: str) -> int:
    """Return the atomic number associated with a chemical symbol.

    Parameters
    ----------
    symbol : str
        Element symbol. Leading and trailing whitespace are ignored and
        capitalization is normalized.

    Returns
    -------
    int
        Atomic number.

    Raises
    ------
    ValueError
        If ``symbol`` is empty or not a recognized element symbol.
    """
    if not isinstance(symbol, str) or not symbol.strip():
        raise ValueError("symbol must be a non-empty string")
    normalized = symbol.strip().capitalize()
    try:
        return ATOMIC_NUMBERS[normalized]
    except KeyError as exc:
        raise ValueError(f"unknown chemical symbol: {symbol!r}") from exc
