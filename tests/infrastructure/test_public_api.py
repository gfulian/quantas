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
    assert registry.get("eos").has(Capability.FIT)
    assert registry.get("eos").has(Capability.PLOT_INVENTORY)
    assert registry.get("eos").operation(Capability.PLOT_INVENTORY) is eos.describe_plots
    assert not registry.get("eos").has(Capability.RUN)
    assert registry.get("thermoelasticity").has(Capability.INTEROP)
    operation = registry.get("elasticity").operation(Capability.RUN)
    assert operation is elasticity.run


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
