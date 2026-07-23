# -*- coding: utf-8 -*-

"""Stable public Python API for Quantas.

The package exposes one namespace per scientific domain.  Numerical and
workflow implementations remain under :mod:`quantas.core` and
:mod:`quantas.modules`; applications, notebooks, the CLI, and future graphical
frontends should depend on this package instead of those implementation paths.

Examples
--------
>>> from quantas.api import qha
>>> options = qha.Options()
>>> result = qha.run("input.yaml", options=options)
"""

from __future__ import annotations

from importlib import import_module
from types import ModuleType
from typing import TYPE_CHECKING

_PUBLIC_NAMESPACES = frozenset(
    {
        "common",
        "elasticity",
        "eos",
        "ha",
        "interop",
        "profiles",
        "qha",
        "rendering",
        "registry",
        "seismic",
        "thermoelasticity",
    }
)

if TYPE_CHECKING:
    from . import (
        common as common,
        elasticity as elasticity,
        eos as eos,
        ha as ha,
        interop as interop,
        profiles as profiles,
        qha as qha,
        rendering as rendering,
        registry as registry,
        seismic as seismic,
        thermoelasticity as thermoelasticity,
    )


def __getattr__(name: str) -> ModuleType:
    """Load a public API namespace on first access.

    Parameters
    ----------
    name : str
        Public namespace name.

    Returns
    -------
    ModuleType
        Imported public API module.

    Raises
    ------
    AttributeError
        If ``name`` is not a supported public namespace.
    """
    if name not in _PUBLIC_NAMESPACES:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(f"{__name__}.{name}")
    globals()[name] = module
    return module


def __dir__() -> list[str]:
    """Return only the supported public API namespace names."""
    return sorted(_PUBLIC_NAMESPACES)


__all__ = sorted(_PUBLIC_NAMESPACES)
