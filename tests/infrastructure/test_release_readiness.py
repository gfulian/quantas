"""Release-readiness contracts for publication tooling."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

import quantas


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _release_identity_module():
    """Load the release-identity tool as a Python module."""
    script = PROJECT_ROOT / "tools" / "check_release_identity.py"
    specification = importlib.util.spec_from_file_location(
        "quantas_check_release_identity",
        script,
    )
    assert specification is not None
    assert specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_release_identity_accepts_current_source_metadata() -> None:
    """Release-facing metadata should agree with the authoritative version."""
    module = _release_identity_module()
    assert module.validate_release_identity() == quantas.__version__


def test_release_identity_requires_exact_version_tag() -> None:
    """Publication tags must use the exact ``v<version>`` form."""
    module = _release_identity_module()
    expected = f"v{quantas.__version__}"
    assert module.validate_release_identity(tag=expected) == quantas.__version__

    with pytest.raises(RuntimeError, match="does not match package version"):
        module.validate_release_identity(tag=quantas.__version__)


def test_publication_workflow_cannot_manually_dispatch_to_pypi() -> None:
    """Production publication should require a tagged GitHub release."""
    workflow = (PROJECT_ROOT / ".github/workflows/release.yml").read_text(
        encoding="utf-8"
    )

    assert "- testpypi" in workflow
    assert "- pypi" not in workflow
    assert "if: github.event_name == 'release'" in workflow
    assert "github.event.release.tag_name" in workflow
    assert "python tools/check_release_identity.py --tag" in workflow


def test_publication_workflow_repeats_the_release_quality_gates() -> None:
    """The publishing workflow should not be weaker than source validation."""
    workflow = (PROJECT_ROOT / ".github/workflows/release.yml").read_text(
        encoding="utf-8"
    )
    required_commands = (
        "python tools/update_examples_manifest.py --check",
        "ruff check src tests tools docs/tools",
        "mypy",
        "python -m compileall -q src tests tools docs/tools",
        "python tools/run_tests.py all -- -q",
        "python tools/check_architecture.py --root .",
        "python -m sphinx -E -a -b html -W --keep-going docs/source docs/_build/html",
        "python -m build",
        "python -m twine check dist/*",
        "python tools/check_distribution.py dist",
        "python tools/check_repository.py",
    )
    for command in required_commands:
        assert command in workflow
