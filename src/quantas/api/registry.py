# -*- coding: utf-8 -*-

"""Capability-based discovery of Quantas public scientific modules."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from enum import Enum
from importlib import import_module
from pathlib import Path
from types import ModuleType
from typing import Any

import h5py

from quantas.io.hdf5.attrs import decode_text

from .common import _public_dir


class Capability(str, Enum):
    """Frontend-neutral operation supported by a public API namespace.

    The values describe workflow capabilities rather than implementation
    classes. Frontends use them to discover operations without importing Click,
    Rich, Dash, Matplotlib, or module-internal calculators.

    Notes
    -----
    A capability is declared only when the corresponding operation is part of
    the supported :mod:`quantas.api` contract. Scientific modules are not
    required to implement the same capability set.
    """

    READ_INPUT = "read_input"
    NORMALIZE_INPUT = "normalize_input"
    CREATE_INPUT = "create_input"
    INSPECT = "inspect"
    RUN = "run"
    RUN_CONTEXT = "run_context"
    FIT = "fit"
    BATCH = "batch"
    READ_RESULT = "read_result"
    WRITE_RESULT = "write_result"
    REPORT = "report"
    PLOT = "plot"
    PLOT_INVENTORY = "plot_inventory"
    EXPORT = "export"
    ARCHIVE = "archive"
    CALCULATE = "calculate"
    DIAGNOSE = "diagnose"
    PROFILE = "profile"
    INTEROP = "interop"


@dataclass(frozen=True, slots=True)
class ModuleDescriptor:
    """Describe one public scientific namespace without importing it eagerly.

    Parameters
    ----------
    name : str
        Stable module identifier used in result metadata.
    title : str
        Human-readable scientific module title.
    api_module : str
        Import path of the public namespace.
    result_key : str or None
        Key used in :class:`quantas.models.ResultData`, when applicable.
    capabilities : frozenset of Capability
        Operations supported by the namespace.
    operations : tuple of tuple
        Mapping from capabilities to public function names.
    input_type_name, options_type_name, result_type_name : str or None
        Public type aliases resolved lazily from the namespace.
    """

    name: str
    title: str
    api_module: str
    result_key: str | None
    capabilities: frozenset[Capability]
    operations: tuple[tuple[Capability, str], ...]
    input_type_name: str | None = None
    options_type_name: str | None = None
    result_type_name: str | None = None

    def load(self) -> ModuleType:
        """Import and return the public namespace lazily.

        Returns
        -------
        ModuleType
            Imported module identified by :attr:`api_module`.

        Raises
        ------
        ImportError
            If the declared namespace cannot be imported.
        """
        return import_module(self.api_module)

    def has(self, capability: Capability | str) -> bool:
        """Return whether one capability is declared.

        Parameters
        ----------
        capability : Capability or str
            Capability enum member or its stable string value.

        Returns
        -------
        bool
            ``True`` when the module advertises the capability.

        Raises
        ------
        ValueError
            If a string is not a valid :class:`Capability` value.
        """
        return Capability(capability) in self.capabilities

    def operation(self, capability: Capability | str) -> Callable[..., Any]:
        """Resolve the public function implementing one capability.

        Raises
        ------
        ValueError
            If the capability is not declared for this module.
        AttributeError
            If the declared public operation is absent.
        """
        resolved = Capability(capability)
        if resolved not in self.capabilities:
            raise ValueError(
                f"module {self.name!r} does not support {resolved.value!r}"
            )
        operation_name = dict(self.operations).get(resolved)
        if operation_name is None:
            raise ValueError(
                f"module {self.name!r} does not expose an operation for "
                f"{resolved.value!r}"
            )
        operation = getattr(self.load(), operation_name)
        if not callable(operation):
            raise AttributeError(
                f"{self.api_module}.{operation_name} is not callable"
            )
        return operation

    def resolve_type(self, name: str | None) -> type[Any] | None:
        """Resolve a declared public type alias lazily.

        Parameters
        ----------
        name : str or None
            Public type alias in the module namespace. ``None`` represents an
            unsupported or inapplicable contract type.

        Returns
        -------
        type or None
            Resolved public class, or ``None`` when no alias was declared.

        Raises
        ------
        AttributeError
            If the declared alias is missing from the namespace.
        TypeError
            If the resolved object is not a class.
        """
        if name is None:
            return None
        value = getattr(self.load(), name)
        if not isinstance(value, type):
            raise TypeError(f"{self.api_module}.{name} is not a type")
        return value

    @property
    def input_type(self) -> type[Any] | None:
        """Return the public input type, when applicable.

        Returns
        -------
        type or None
            Namespace ``Input``/dataset class, or ``None`` when the workflow
            has no single input contract.
        """
        return self.resolve_type(self.input_type_name)

    @property
    def options_type(self) -> type[Any] | None:
        """Return the public options type, when applicable.

        Returns
        -------
        type or None
            Namespace options class, or ``None`` when configuration is
            operation-specific.
        """
        return self.resolve_type(self.options_type_name)

    @property
    def result_type(self) -> type[Any] | None:
        """Return the public scientific result type, when applicable.

        Returns
        -------
        type or None
            Typed scientific payload or fit-result class exposed by the
            namespace, or ``None`` when no single result type applies.
        """
        return self.resolve_type(self.result_type_name)


_STANDARD_RESULT_OPERATIONS = (
    (Capability.READ_INPUT, "read_input"),
    (Capability.NORMALIZE_INPUT, "normalize_input"),
    (Capability.RUN, "run"),
    (Capability.READ_RESULT, "read_result"),
    (Capability.WRITE_RESULT, "write_result"),
    (Capability.REPORT, "build_report"),
    (Capability.PLOT, "build_plots"),
)

_PLOT_INVENTORY_OPERATION = (Capability.PLOT_INVENTORY, "describe_plots")


_DESCRIPTORS: tuple[ModuleDescriptor, ...] = (
    ModuleDescriptor(
        name="elasticity",
        title="Second-order elasticity",
        api_module="quantas.api.elasticity",
        result_key="elasticity",
        capabilities=frozenset(
            {
                *(capability for capability, _ in _STANDARD_RESULT_OPERATIONS),
                Capability.PLOT_INVENTORY,
            }
        ),
        operations=(*_STANDARD_RESULT_OPERATIONS, _PLOT_INVENTORY_OPERATION),
        input_type_name="Input",
        options_type_name="Options",
        result_type_name="Result",
    ),
    ModuleDescriptor(
        name="seismic",
        title="Directional seismic-wave analysis",
        api_module="quantas.api.seismic",
        result_key="seismic",
        capabilities=frozenset(
            {
                *(capability for capability, _ in _STANDARD_RESULT_OPERATIONS),
                Capability.EXPORT,
                Capability.PLOT_INVENTORY,
            }
        ),
        operations=(
            *_STANDARD_RESULT_OPERATIONS,
            _PLOT_INVENTORY_OPERATION,
            (Capability.EXPORT, "write_csv"),
        ),
        input_type_name="Input",
        options_type_name="Options",
        result_type_name="Result",
    ),
    ModuleDescriptor(
        name="ha",
        title="Harmonic thermodynamics",
        api_module="quantas.api.ha",
        result_key="ha",
        capabilities=frozenset(
            {
                *(capability for capability, _ in _STANDARD_RESULT_OPERATIONS),
                Capability.CREATE_INPUT,
                Capability.PLOT_INVENTORY,
            }
        ),
        operations=(
            *_STANDARD_RESULT_OPERATIONS,
            _PLOT_INVENTORY_OPERATION,
            (Capability.CREATE_INPUT, "create_input"),
        ),
        input_type_name="Input",
        options_type_name="Options",
        result_type_name="Result",
    ),
    ModuleDescriptor(
        name="qha",
        title="Quasi-harmonic thermodynamics",
        api_module="quantas.api.qha",
        result_key="qha",
        capabilities=frozenset(
            {
                *(capability for capability, _ in _STANDARD_RESULT_OPERATIONS),
                Capability.INSPECT,
            }
        ),
        operations=(*_STANDARD_RESULT_OPERATIONS, (Capability.INSPECT, "inspect")),
        input_type_name="Input",
        options_type_name="Options",
        result_type_name="Result",
    ),
    ModuleDescriptor(
        name="eos",
        title="Equation-of-state fitting and analysis",
        api_module="quantas.api.eos",
        result_key=None,
        capabilities=frozenset(
            {
                Capability.READ_INPUT,
                Capability.NORMALIZE_INPUT,
                Capability.FIT,
                Capability.BATCH,
                Capability.ARCHIVE,
                Capability.CALCULATE,
                Capability.DIAGNOSE,
                Capability.PLOT,
            }
        ),
        operations=(
            (Capability.READ_INPUT, "read_input"),
            (Capability.NORMALIZE_INPUT, "normalize_input"),
            (Capability.FIT, "fit"),
            (Capability.BATCH, "run_batch"),
            (Capability.ARCHIVE, "open_archive"),
            (Capability.CALCULATE, "calculate"),
            (Capability.DIAGNOSE, "diagnose"),
            (Capability.PLOT, "build_plots"),
        ),
        input_type_name="Dataset",
        options_type_name="FitOptions",
        result_type_name="FitResult",
    ),
    ModuleDescriptor(
        name="thermoelasticity",
        title="Quasi-static thermoelasticity",
        api_module="quantas.api.thermoelasticity",
        result_key="thermoelasticity",
        capabilities=frozenset(
            {
                *(capability for capability, _ in _STANDARD_RESULT_OPERATIONS),
                Capability.CREATE_INPUT,
                Capability.RUN_CONTEXT,
                Capability.PROFILE,
                Capability.EXPORT,
                Capability.INTEROP,
            }
        ),
        operations=(
            *_STANDARD_RESULT_OPERATIONS,
            (Capability.CREATE_INPUT, "create_input"),
            (Capability.RUN_CONTEXT, "run_context"),
            (Capability.PROFILE, "analyze_profiles"),
            (Capability.EXPORT, "write_tensor_export"),
            (Capability.INTEROP, "prepare_context"),
        ),
        input_type_name="Input",
        options_type_name="Options",
        result_type_name="Result",
    ),
)

_BY_NAME = {descriptor.name: descriptor for descriptor in _DESCRIPTORS}


def list_modules() -> tuple[ModuleDescriptor, ...]:
    """Return all public scientific module descriptors.

    Returns
    -------
    tuple of ModuleDescriptor
        Immutable descriptors in stable frontend display order.
    """
    return _DESCRIPTORS


def iter_modules() -> Iterator[ModuleDescriptor]:
    """Iterate over public scientific module descriptors.

    Returns
    -------
    iterator of ModuleDescriptor
        Lazy iterator over the stable descriptor sequence.
    """
    return iter(_DESCRIPTORS)


def get(name: str) -> ModuleDescriptor:
    """Return one public module descriptor by stable identifier.

    Parameters
    ----------
    name : str
        Stable module identifier such as ``qha`` or ``eos``.

    Returns
    -------
    ModuleDescriptor
        Public API and capability descriptor.

    Raises
    ------
    KeyError
        If the module identifier is unknown.
    """
    try:
        return _BY_NAME[name]
    except KeyError as exc:
        available = ", ".join(sorted(_BY_NAME))
        raise KeyError(f"unknown Quantas module {name!r}; available: {available}") from exc


def module_from_result(path: str | Path) -> ModuleDescriptor:
    """Inspect native HDF5 metadata and return the responsible module.

    Parameters
    ----------
    path : str or Path
        Native Quantas HDF5 result or archive.

    Returns
    -------
    ModuleDescriptor
        Descriptor selected from the persisted ``metadata/module`` value.

    Raises
    ------
    ValueError
        If Quantas metadata are absent.
    KeyError
        If the persisted module identifier is unsupported.
    """
    with h5py.File(path, "r") as handle:
        metadata = handle.get("metadata")
        if metadata is None:
            raise ValueError("HDF5 file does not contain Quantas metadata")
        module = decode_text(metadata.attrs.get("module", "unknown"))
    return get(module)


def open_result(path: str | Path, *, writable: bool = False) -> Any:
    """Open a native Quantas result using metadata-driven dispatch.

    Parameters
    ----------
    path : str or Path
        Native Quantas HDF5 result or archive.
    writable : bool, optional
        Request writable access. This is supported only by archive-style
        modules such as EOS.

    Returns
    -------
    Any
        :class:`quantas.models.ResultData` for single-shot modules or an active
        module-specific archive for archive-style workflows.

    Raises
    ------
    ValueError
        If writable access is requested for a single-shot result or the module
        has no registered native reader.
    KeyError
        If the persisted module identifier is unsupported.

    Notes
    -----
    EOS returns an active archive because it stores datasets and immutable fit
    records rather than one result envelope. Callers must close that archive.
    """
    descriptor = module_from_result(path)
    if descriptor.has(Capability.READ_RESULT):
        if writable:
            raise ValueError(
                f"writable mode is not supported for {descriptor.name!r} results"
            )
        return descriptor.operation(Capability.READ_RESULT)(path)
    if descriptor.has(Capability.ARCHIVE):
        return descriptor.operation(Capability.ARCHIVE)(path, writable=writable)
    raise ValueError(f"module {descriptor.name!r} has no native result reader")


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "Capability",
    "ModuleDescriptor",
    "get",
    "iter_modules",
    "list_modules",
    "module_from_result",
    "open_result",
]
