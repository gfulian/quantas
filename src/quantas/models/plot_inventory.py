# -*- coding: utf-8 -*-

"""Frontend-neutral discovery contracts for scientific plot builders.

The descriptors in this module expose stable scientific metadata before a
plot is built. They deliberately describe quantities, representations, and
selection contexts without defining renderer controls or replacing the typed
module-specific build options.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal, TypeAlias, get_args


PlotKind: TypeAlias = Literal[
    "line",
    "contour",
    "polar",
    "spherical_map",
    "spherical_summary",
    "surface",
    "panel",
]
PlotContextValue: TypeAlias = str | int | float | bool

_KEY_PATTERN = re.compile(r"^[a-z][a-z0-9_]*$")
_PLOT_KINDS: frozenset[str] = frozenset(get_args(PlotKind))


def _validate_key(value: str, *, field_name: str) -> None:
    """Validate one stable machine-readable descriptor key."""
    if not isinstance(value, str) or not _KEY_PATTERN.fullmatch(value):
        raise ValueError(
            f"{field_name} must match {_KEY_PATTERN.pattern!r}; received {value!r}"
        )


def _validate_text(value: str, *, field_name: str) -> None:
    """Validate one required human-readable text field."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")


def _validate_unique_keys(items: tuple[object, ...], *, field_name: str) -> None:
    """Validate stable ``key`` attributes in one descriptor sequence."""
    keys = [getattr(item, "key") for item in items]
    if len(set(keys)) != len(keys):
        raise ValueError(f"{field_name} contains duplicate keys")


@dataclass(frozen=True, slots=True)
class PlotPropertyDescriptor:
    """Describe one scientifically meaningful plottable quantity.

    Parameters
    ----------
    key : str
        Stable key accepted by the relevant public builder when the builder
        selects properties directly.
    name : str
        Extended human-readable quantity name.
    symbol_math : str
        Mathematical symbol source without renderer delimiters such as ``$``.
    symbol_plain : str
        Plain-text or Unicode symbol suitable for non-mathematical frontends.
    unit : str or None
        Physical unit of the represented values. Dimensionless quantities use
        ``None``.
    description : str, optional
        Short scientific description.
    category : str, optional
        Stable module-local grouping key.
    components : tuple of str, optional
        Stable branches, modes, or components associated with the quantity.
    representations : tuple of str, optional
        Representation keys in the containing :class:`PlotInventory`.
    """

    key: str
    name: str
    symbol_math: str
    symbol_plain: str
    unit: str | None
    description: str = ""
    category: str = ""
    components: tuple[str, ...] = ()
    representations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate required scientific metadata and stable references."""
        _validate_key(self.key, field_name="property key")
        _validate_text(self.name, field_name="property name")
        _validate_text(self.symbol_math, field_name="mathematical symbol")
        _validate_text(self.symbol_plain, field_name="plain symbol")
        if self.unit is not None and not self.unit.strip():
            raise ValueError("property unit must be None or a non-empty string")
        if self.category:
            _validate_key(self.category, field_name="property category")
        for component in self.components:
            _validate_key(component, field_name="property component")
        if len(set(self.components)) != len(self.components):
            raise ValueError("property components must be unique")
        for representation in self.representations:
            _validate_key(representation, field_name="property representation")
        if len(set(self.representations)) != len(self.representations):
            raise ValueError("property representations must be unique")


@dataclass(frozen=True, slots=True)
class PlotContextDescriptor:
    """Describe one scientific selection or result context.

    Parameters
    ----------
    key : str
        Stable context identifier.
    name : str
        Human-readable context name.
    description : str, optional
        Scientific meaning of the context.
    values : tuple, optional
        Exact supported values. Numeric grid values are expressed in ``unit``.
    unit : str or None, optional
        Unit associated with numeric values.
    default : str, int, float, bool, or None, optional
        Public default when the context is selectable.
    required : bool, optional
        Whether a caller must choose a value explicitly.
    selectable : bool, optional
        ``False`` for informative result context such as a sampled grid or the
        calculation level already fixed in the stored result.
    """

    key: str
    name: str
    description: str = ""
    values: tuple[PlotContextValue, ...] = ()
    unit: str | None = None
    default: PlotContextValue | None = None
    required: bool = False
    selectable: bool = True

    def __post_init__(self) -> None:
        """Validate values and defaults without coercing scientific data."""
        _validate_key(self.key, field_name="context key")
        _validate_text(self.name, field_name="context name")
        if self.unit is not None and not self.unit.strip():
            raise ValueError("context unit must be None or a non-empty string")
        fingerprints = [(type(value), repr(value)) for value in self.values]
        if len(set(fingerprints)) != len(fingerprints):
            raise ValueError("context values must be unique")
        if self.default is not None and self.default not in self.values:
            raise ValueError("context default must be present in context values")
        if self.required and not self.values:
            raise ValueError("required contexts must expose at least one value")


@dataclass(frozen=True, slots=True)
class PlotRepresentationDescriptor:
    """Describe one scientific representation supported by a module.

    Parameters
    ----------
    key : str
        Stable representation identifier.
    name : str
        Human-readable representation name.
    plot_kind : PlotKind
        Frontend-neutral structural kind returned by the builder.
    description : str, optional
        Scientific meaning of the representation.
    property_keys : tuple of str, optional
        Compatible property keys. An empty tuple is valid for a representation
        whose quantity is selected entirely through contexts.
    supported_contexts : tuple of str, optional
        Context keys required or understood by the representation.
    constraints : tuple of str, optional
        Human-readable scientific limitations that cannot be expressed by the
        simple compatibility fields.
    """

    key: str
    name: str
    plot_kind: PlotKind
    description: str = ""
    property_keys: tuple[str, ...] = ()
    supported_contexts: tuple[str, ...] = ()
    constraints: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate stable property and context references."""
        _validate_key(self.key, field_name="representation key")
        _validate_text(self.name, field_name="representation name")
        if not isinstance(self.plot_kind, str) or self.plot_kind not in _PLOT_KINDS:
            raise ValueError(
                f"plot kind must be one of {sorted(_PLOT_KINDS)!r}; "
                f"received {self.plot_kind!r}"
            )
        for property_key in self.property_keys:
            _validate_key(property_key, field_name="representation property key")
        if len(set(self.property_keys)) != len(self.property_keys):
            raise ValueError("representation property keys must be unique")
        for context in self.supported_contexts:
            _validate_key(context, field_name="representation context")
        if len(set(self.supported_contexts)) != len(self.supported_contexts):
            raise ValueError("representation contexts must be unique")
        if any(not constraint.strip() for constraint in self.constraints):
            raise ValueError("representation constraints must be non-empty strings")


@dataclass(frozen=True, slots=True)
class PlotInventory:
    """Complete result-aware plot discovery response for one module.

    Parameters
    ----------
    module : str
        Stable Quantas module identifier.
    properties : tuple of PlotPropertyDescriptor
        Quantities that can be built from the supplied result.
    representations : tuple of PlotRepresentationDescriptor
        Available scientific representation families.
    contexts : tuple of PlotContextDescriptor, optional
        Scientific selections and informative result context.
    warnings : tuple of str, optional
        Non-fatal discovery limitations.
    """

    module: str
    properties: tuple[PlotPropertyDescriptor, ...]
    representations: tuple[PlotRepresentationDescriptor, ...]
    contexts: tuple[PlotContextDescriptor, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate uniqueness and cross-references across the inventory."""
        _validate_key(self.module, field_name="inventory module")
        _validate_unique_keys(self.properties, field_name="inventory properties")
        _validate_unique_keys(
            self.representations,
            field_name="inventory representations",
        )
        _validate_unique_keys(self.contexts, field_name="inventory contexts")
        if any(not warning.strip() for warning in self.warnings):
            raise ValueError("inventory warnings must be non-empty strings")

        property_keys = {item.key for item in self.properties}
        representation_keys = {item.key for item in self.representations}
        context_keys = {item.key for item in self.contexts}
        for item in self.properties:
            unknown = set(item.representations) - representation_keys
            if unknown:
                raise ValueError(
                    f"property {item.key!r} references unknown representations: "
                    f"{sorted(unknown)}"
                )
        for item in self.representations:
            unknown_properties = set(item.property_keys) - property_keys
            if unknown_properties:
                raise ValueError(
                    f"representation {item.key!r} references unknown properties: "
                    f"{sorted(unknown_properties)}"
                )
            unknown_contexts = set(item.supported_contexts) - context_keys
            if unknown_contexts:
                raise ValueError(
                    f"representation {item.key!r} references unknown contexts: "
                    f"{sorted(unknown_contexts)}"
                )

    def property_by_key(self, key: str) -> PlotPropertyDescriptor:
        """Return one property descriptor by stable key.

        Raises
        ------
        KeyError
            If the key is not present in this inventory.
        """
        for item in self.properties:
            if item.key == key:
                return item
        raise KeyError(f"unknown plot property {key!r} for module {self.module!r}")

    def representation_by_key(self, key: str) -> PlotRepresentationDescriptor:
        """Return one representation descriptor by stable key."""
        for item in self.representations:
            if item.key == key:
                return item
        raise KeyError(
            f"unknown plot representation {key!r} for module {self.module!r}"
        )

    def context_by_key(self, key: str) -> PlotContextDescriptor:
        """Return one context descriptor by stable key."""
        for item in self.contexts:
            if item.key == key:
                return item
        raise KeyError(f"unknown plot context {key!r} for module {self.module!r}")


__all__ = [
    "PlotContextDescriptor",
    "PlotContextValue",
    "PlotInventory",
    "PlotKind",
    "PlotPropertyDescriptor",
    "PlotRepresentationDescriptor",
]
