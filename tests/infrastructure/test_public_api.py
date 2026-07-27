"""Contracts for the organized Quantas public Python API."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from quantas import api
from quantas.api import (
    elasticity,
    eos,
    ha,
    plotting,
    qha,
    registry,
    rendering,
    seismic,
    thermoelasticity,
)
from quantas.api.registry import Capability


EXPECTED_NAMESPACES = {
    "common",
    "elasticity",
    "eos",
    "ha",
    "interop",
    "profiles",
    "plotting",
    "qha",
    "rendering",
    "registry",
    "seismic",
    "thermoelasticity",
}


def test_public_api_is_organized_by_scientific_namespace() -> None:
    """The facade exposes domain namespaces instead of a flat symbol list."""
    assert set(api.__all__) == EXPECTED_NAMESPACES
    assert api.elasticity is elasticity
    assert api.eos is eos
    assert api.ha is ha
    assert api.plotting is plotting
    assert api.qha is qha
    assert api.rendering is rendering
    assert api.seismic is seismic
    assert api.thermoelasticity is thermoelasticity
    assert not hasattr(api, "run_qha")
    assert not hasattr(api, "fit_eos")


def test_public_domain_aliases_avoid_repeated_domain_names() -> None:
    """Types and operations use their namespace as scientific context."""
    assert elasticity.Input.__name__ == "ElasticityInput"
    assert elasticity.Options.__name__ == "ElasticityOptions"
    assert seismic.Result.__name__ == "SeismicResult"
    assert ha.run.__name__ == "run"
    assert qha.inspect.__name__ == "inspect"
    assert eos.FitRequest.__name__ == "EOSFitRequest"
    assert thermoelasticity.Context.__name__ == "ThermoelasticContext"


def test_root_eos_exception_was_removed() -> None:
    """EOS is available through the same public hierarchy as every module."""
    assert importlib.util.find_spec("quantas.eos") is None


def test_registry_declares_all_scientific_modules_and_types() -> None:
    """The registry supports lazy, capability-based frontend discovery."""
    descriptors = registry.list_modules()
    assert tuple(item.name for item in descriptors) == (
        "elasticity",
        "seismic",
        "ha",
        "qha",
        "eos",
        "thermoelasticity",
    )
    for descriptor in descriptors:
        assert descriptor.load().__name__ == f"quantas.api.{descriptor.name}"
        assert descriptor.input_type is not None
        assert descriptor.options_type is not None
        assert descriptor.result_type is not None

    assert registry.get("qha").has(Capability.INSPECT)
    assert registry.get("qha").has(Capability.CREATE_INPUT)
    assert registry.get("qha").has(Capability.EXPORT)
    assert registry.get("eos").has(Capability.FIT)
    assert registry.get("eos").has(Capability.PLOT_INVENTORY)
    assert registry.get("eos").operation(Capability.PLOT_INVENTORY) is eos.describe_plots
    assert not registry.get("eos").has(Capability.RUN)
    assert registry.get("thermoelasticity").has(Capability.INTEROP)
    assert registry.get("eos").has(Capability.TEMPLATE)
    operation = registry.get("elasticity").operation(Capability.RUN)
    assert operation is elasticity.run


def test_registry_describes_multiple_named_operations() -> None:
    """Frontends can enumerate operation families without guessing CLI commands."""
    elasticity_descriptor = registry.get("elasticity")
    assert elasticity_descriptor.named_operation("create_input") is elasticity.create_input
    assert (
        elasticity_descriptor.named_operation("export_2d_table")
        is elasticity.write_table
    )

    eos_exports = registry.get("eos").operations_for(Capability.EXPORT)
    assert eos_exports == (eos.write_diagnostic_csv, eos.write_calculation_csv)

    thermo_exports = registry.get("thermoelasticity").list_operations(
        Capability.EXPORT
    )
    assert tuple(item.key for item in thermo_exports) == (
        "export_tensor",
        "export_grid_table",
        "export_profile_table",
        "export_seismic_input",
    )


def test_registry_named_operations_are_unique_and_public() -> None:
    """Every discovered lifecycle operation resolves through its public facade."""
    for descriptor in registry.list_modules():
        namespace = descriptor.load()
        operations = descriptor.list_operations()
        keys = tuple(item.key for item in operations)
        assert len(keys) == len(set(keys))
        for item in operations:
            assert item.capability in descriptor.capabilities
            assert item.function_name in namespace.__all__
            assert item.resolve(namespace) is getattr(namespace, item.function_name)


def test_registry_rejects_unknown_modules_and_capabilities() -> None:
    with pytest.raises(KeyError, match="unknown Quantas module"):
        registry.get("unknown")
    with pytest.raises(ValueError, match="does not support"):
        registry.get("ha").operation(Capability.FIT)


def test_import_quantas_remains_lightweight(tmp_path: Path) -> None:
    """Importing package metadata does not eagerly import numerical stacks."""
    script = tmp_path / "check_import.py"
    script.write_text(
        """
import json
import sys
import quantas
names = ("scipy", "h5py", "matplotlib")
print(json.dumps({name: name in sys.modules for name in names}))
""",
        encoding="utf-8",
    )
    env = os.environ.copy()
    source_root = str(Path(__file__).resolve().parents[2] / "src")
    env["PYTHONPATH"] = os.pathsep.join(
        item for item in (source_root, env.get("PYTHONPATH", "")) if item
    )
    completed = subprocess.run(
        [sys.executable, str(script)],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    loaded = json.loads(completed.stdout)
    assert loaded == {
        "scipy": False,
        "h5py": False,
        "matplotlib": False,
    }


def test_public_input_and_options_annotations_are_closed() -> None:
    """Public workflow contracts require no imports from implementation packages."""
    import dataclasses
    import typing

    public_namespaces = (
        api.common,
        elasticity,
        eos,
        ha,
        plotting,
        qha,
        registry,
        seismic,
        thermoelasticity,
    )
    public_identities = {
        id(getattr(namespace, name))
        for namespace in public_namespaces
        for name in namespace.__all__
    }
    contracts = (
        elasticity.Input,
        elasticity.Options,
        elasticity.SurfaceOptions,
        seismic.Input,
        seismic.Options,
        seismic.SurfaceOptions,
        ha.Input,
        ha.Options,
        ha.PlotOptions,
        qha.Input,
        qha.Options,
        qha.PlotOptions,
        thermoelasticity.Input,
        thermoelasticity.Options,
        thermoelasticity.Context,
        eos.Dataset,
        eos.FitOptions,
        eos.FitRequest,
        eos.PVTModel,
        eos.MGDNormalization,
    )

    def iter_types(annotation):
        if isinstance(annotation, type):
            yield annotation
        for argument in typing.get_args(annotation):
            yield from iter_types(argument)

    missing: list[str] = []
    for contract in contracts:
        assert dataclasses.is_dataclass(contract)
        for field_name, annotation in typing.get_type_hints(contract).items():
            for referenced in iter_types(annotation):
                if not referenced.__module__.startswith("quantas."):
                    continue
                if id(referenced) not in public_identities:
                    missing.append(
                        f"{contract.__name__}.{field_name}: "
                        f"{referenced.__module__}.{referenced.__name__}"
                    )

    assert missing == []


def test_public_function_annotations_are_closed() -> None:
    """Public callables resolve hints using only public Quantas identities."""
    import inspect
    import typing

    public_namespaces = (
        api.common,
        elasticity,
        eos,
        ha,
        plotting,
        qha,
        registry,
        seismic,
        thermoelasticity,
    )
    public_identities = {
        id(getattr(namespace, name))
        for namespace in public_namespaces
        for name in namespace.__all__
    }

    def iter_types(annotation):
        if isinstance(annotation, type):
            yield annotation
        for argument in typing.get_args(annotation):
            yield from iter_types(argument)

    missing: list[str] = []
    unresolved: list[str] = []
    for namespace in public_namespaces:
        for name in namespace.__all__:
            value = getattr(namespace, name)
            if not inspect.isfunction(value):
                continue
            try:
                annotations = typing.get_type_hints(value)
            except Exception as exc:
                unresolved.append(
                    f"{namespace.__name__}.{name}: {type(exc).__name__}: {exc}"
                )
                continue
            for parameter, annotation in annotations.items():
                for referenced in iter_types(annotation):
                    if not referenced.__module__.startswith("quantas."):
                        continue
                    if id(referenced) not in public_identities:
                        missing.append(
                            f"{namespace.__name__}.{name}.{parameter}: "
                            f"{referenced.__module__}.{referenced.__name__}"
                        )

    assert unresolved == []
    assert missing == []
