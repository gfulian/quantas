# -*- coding: utf-8 -*-

"""Post-fit pressure-temperature reconstruction from thermoelastic archives.

This module evaluates the already fitted quasi-static elastic model without
repeating either the static EOS fit or the independent-component fits.  The
QHA equilibrium-volume surface archived by the fitting run supplies
``V(P, T)`` through rectilinear interpolation or controlled extrapolation.
"""

from __future__ import annotations

from dataclasses import fields
from typing import Any, Mapping


from quantas.core.math.fitting import FitMethod
from quantas.modules.thermoelasticity.models import (
    ThermoelasticOptions,
)


def thermoelastic_options_from_mapping(
    values: Mapping[str, Any],
) -> ThermoelasticOptions:
    """Rebuild :class:`ThermoelasticOptions` from an HDF5 option mapping.

    Unknown keys are ignored so archives remain readable after new frontend
    controls are introduced.
    """
    names = {item.name for item in fields(ThermoelasticOptions)}
    kwargs = {name: values[name] for name in names if name in values}
    if "fit_method" in kwargs:
        kwargs["fit_method"] = FitMethod(kwargs["fit_method"])
    return ThermoelasticOptions(**kwargs)


__all__ = ["thermoelastic_options_from_mapping"]
