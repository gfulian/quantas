# -*- coding: utf-8 -*-

"""Definitions and parsing rules for equations of state.

An equation of state is identified by a physical family and, when meaningful,
by the order of its strain expansion.  The model specification is shared by
energy-volume and pressure-volume implementations so that a fitted integrated
energy EOS is always evaluated with the matching pressure EOS.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import re


class EOSFamily(str, Enum):
    """Supported isothermal equation-of-state families."""

    MURNAGHAN = "murnaghan"
    BIRCH_MURNAGHAN = "birchmurnaghan"
    NATURAL_STRAIN = "poirier-tarantola"
    VINET = "vinet"
    TAIT = "tait"


_FAMILY_ALIASES: dict[str, EOSFamily] = {
    "m": EOSFamily.MURNAGHAN,
    "murnaghan": EOSFamily.MURNAGHAN,
    "bm": EOSFamily.BIRCH_MURNAGHAN,
    "birchmurnaghan": EOSFamily.BIRCH_MURNAGHAN,
    "birch-murnaghan": EOSFamily.BIRCH_MURNAGHAN,
    "birch_murnaghan": EOSFamily.BIRCH_MURNAGHAN,
    "ns": EOSFamily.NATURAL_STRAIN,
    "pt": EOSFamily.NATURAL_STRAIN,
    "naturalstrain": EOSFamily.NATURAL_STRAIN,
    "natural-strain": EOSFamily.NATURAL_STRAIN,
    "natural_strain": EOSFamily.NATURAL_STRAIN,
    "poirier-tarantola": EOSFamily.NATURAL_STRAIN,
    "poiriertarantola": EOSFamily.NATURAL_STRAIN,
    "pouriertarantola": EOSFamily.NATURAL_STRAIN,
    "v": EOSFamily.VINET,
    "vinet": EOSFamily.VINET,
    "t": EOSFamily.TAIT,
    "tait": EOSFamily.TAIT,
}

_ALLOWED_ORDERS: dict[EOSFamily, tuple[int, ...]] = {
    EOSFamily.MURNAGHAN: (),
    EOSFamily.BIRCH_MURNAGHAN: (2, 3, 4),
    EOSFamily.NATURAL_STRAIN: (2, 3, 4),
    EOSFamily.VINET: (2, 3),
    EOSFamily.TAIT: (2, 3, 4),
}

_DEFAULT_ORDERS: dict[EOSFamily, int | None] = {
    EOSFamily.MURNAGHAN: None,
    EOSFamily.BIRCH_MURNAGHAN: 3,
    EOSFamily.NATURAL_STRAIN: 3,
    EOSFamily.VINET: 3,
    EOSFamily.TAIT: 3,
}

_TAG_PREFIX: dict[EOSFamily, str] = {
    EOSFamily.MURNAGHAN: "M",
    EOSFamily.BIRCH_MURNAGHAN: "BM",
    EOSFamily.NATURAL_STRAIN: "PT",
    EOSFamily.VINET: "V",
    EOSFamily.TAIT: "T",
}

_DISPLAY_NAME: dict[EOSFamily, str] = {
    EOSFamily.MURNAGHAN: "Murnaghan",
    EOSFamily.BIRCH_MURNAGHAN: "Birch-Murnaghan",
    EOSFamily.NATURAL_STRAIN: "Natural strain (Poirier-Tarantola)",
    EOSFamily.VINET: "Vinet",
    EOSFamily.TAIT: "Tait",
}

_ENERGY_ORDERS: dict[EOSFamily, tuple[int | None, ...]] = {
    EOSFamily.MURNAGHAN: (None,),
    EOSFamily.BIRCH_MURNAGHAN: (2, 3, 4),
    EOSFamily.NATURAL_STRAIN: (2, 3, 4),
    EOSFamily.VINET: (2, 3),
    EOSFamily.TAIT: (),
}


@dataclass(frozen=True, slots=True)
class EOSModel:
    """Describe one equation-of-state family and order.

    Parameters
    ----------
    family : EOSFamily
        Physical equation-of-state family.
    order : int or None
        Truncation order.  Murnaghan has no strain-expansion order.
    """

    family: EOSFamily
    order: int | None = None

    def __post_init__(self) -> None:
        """Validate the family-order combination."""
        allowed = _ALLOWED_ORDERS[self.family]
        if not allowed:
            if self.order is not None:
                raise ValueError(f"{self.family.value} does not define an EOS order")
            return
        if self.order not in allowed:
            choices = ", ".join(str(value) for value in allowed)
            raise ValueError(
                f"unsupported order {self.order!r} for {self.family.value}; "
                f"available orders are {choices}"
            )

    @property
    def tag(self) -> str:
        """Return the compact canonical tag, such as ``BM3``."""
        prefix = _TAG_PREFIX[self.family]
        return prefix if self.order is None else f"{prefix}{self.order}"

    @property
    def name(self) -> str:
        """Return a human-readable model name."""
        base = _DISPLAY_NAME[self.family]
        return base if self.order is None else f"{base}, order {self.order}"

    @property
    def supports_energy(self) -> bool:
        """Return whether a volume-integrated form is implemented."""
        return self.order in _ENERGY_ORDERS[self.family]

    @property
    def energy_parameter_names(self) -> tuple[str, ...]:
        """Return free parameters for an energy-volume fit."""
        if not self.supports_energy:
            raise ValueError(f"{self.tag} has no implemented energy-volume form")
        if self.family is EOSFamily.MURNAGHAN:
            return ("E0", "K0", "KP", "V0")
        if self.order == 2:
            return ("E0", "K0", "V0")
        if self.order == 3:
            return ("E0", "K0", "KP", "V0")
        return ("E0", "K0", "KP", "KPP", "V0")

    @property
    def pressure_parameter_names(self) -> tuple[str, ...]:
        """Return free parameters for a pressure-volume fit."""
        if self.family is EOSFamily.MURNAGHAN:
            return ("K0", "KP", "V0")
        if self.order == 2:
            return ("K0", "V0")
        if self.order == 3:
            return ("K0", "KP", "V0")
        return ("K0", "KP", "KPP", "V0")

    @property
    def parameter_sources(self) -> dict[str, str]:
        """Return whether physical parameters are fitted or implied."""
        sources = {"K0": "fitted", "V0": "fitted"}
        if self.family is EOSFamily.MURNAGHAN:
            sources.update({"KP": "fitted", "KPP": "implied"})
        elif self.order == 2:
            sources.update({"KP": "implied", "KPP": "implied"})
        elif self.order == 3:
            sources.update({"KP": "fitted", "KPP": "implied"})
        else:
            sources.update({"KP": "fitted", "KPP": "fitted"})
        return sources

    def as_dict(self) -> dict[str, object]:
        """Return a serializable representation."""
        return {
            "family": self.family.value,
            "order": self.order,
            "tag": self.tag,
            "name": self.name,
            "supports_energy": self.supports_energy,
        }


def available_eos_models(*, require_energy: bool = False) -> tuple[EOSModel, ...]:
    """Return all supported canonical EOS specifications.

    Parameters
    ----------
    require_energy : bool, optional
        If ``True``, return only models with a volume-integrated energy form.

    Returns
    -------
    tuple of EOSModel
        Canonical models in a stable presentation order.
    """
    models: list[EOSModel] = []
    for family in EOSFamily:
        orders = _ALLOWED_ORDERS[family]
        candidates = (
            (EOSModel(family),)
            if not orders
            else tuple(EOSModel(family, order) for order in orders)
        )
        for model in candidates:
            if not require_energy or model.supports_energy:
                models.append(model)
    return tuple(models)


def available_eos_tags(
    *,
    require_energy: bool = False,
    include_default_aliases: bool = False,
) -> tuple[str, ...]:
    """Return compact tags for the supported EOS specifications.

    Parameters
    ----------
    require_energy : bool, optional
        If ``True``, return only models with a volume-integrated energy form.
    include_default_aliases : bool, optional
        Include unqualified family tags such as ``BM`` and ``PT`` before the
        corresponding canonical ordered tags.

    Returns
    -------
    tuple of str
        Stable sequence of compact EOS tags.
    """
    tags: list[str] = []
    seen: set[str] = set()
    models = available_eos_models(require_energy=require_energy)
    for family in EOSFamily:
        family_models = tuple(model for model in models if model.family is family)
        if not family_models:
            continue
        if include_default_aliases:
            alias = _TAG_PREFIX[family]
            if alias not in seen:
                tags.append(alias)
                seen.add(alias)
        for model in family_models:
            if model.tag not in seen:
                tags.append(model.tag)
                seen.add(model.tag)
    return tuple(tags)


def parse_eos_model(
    eos: str | EOSFamily | EOSModel,
    order: int | None = None,
) -> EOSModel:
    """Return a validated equation-of-state model specification.

    Parameters
    ----------
    eos : str, EOSFamily or EOSModel
        Family name, compact tag, or an existing model specification.
        Compact tags may include the order, for example ``BM2`` or ``PT4``.
    order : int or None, optional
        Explicit order.  It must agree with an order embedded in ``eos``.

    Returns
    -------
    EOSModel
        Canonical family-order specification.

    Raises
    ------
    ValueError
        If the family is unknown, the order is unsupported, or two supplied
        orders conflict.
    """
    if isinstance(eos, EOSModel):
        if order is not None and order != eos.order:
            raise ValueError(
                f"conflicting EOS orders: model uses {eos.order}, explicit order is {order}"
            )
        return eos

    if isinstance(eos, EOSFamily):
        family = eos
        embedded_order = None
    else:
        raw = str(eos).strip().lower()
        match = re.fullmatch(r"(.+?)([234])?", raw)
        if match is None:
            raise ValueError(f"unknown equation of state: {eos!r}")
        base, suffix = match.groups()
        try:
            family = _FAMILY_ALIASES[base]
        except KeyError as exc:
            raise ValueError(f"unknown equation of state: {eos!r}") from exc
        embedded_order = None if suffix is None else int(suffix)

    if embedded_order is not None and order is not None and embedded_order != order:
        raise ValueError(
            f"conflicting EOS orders: tag uses {embedded_order}, explicit order is {order}"
        )
    selected_order = embedded_order if embedded_order is not None else order
    if selected_order is None:
        selected_order = _DEFAULT_ORDERS[family]
    return EOSModel(family=family, order=selected_order)


__all__ = [
    "EOSFamily",
    "EOSModel",
    "available_eos_models",
    "available_eos_tags",
    "parse_eos_model",
]
