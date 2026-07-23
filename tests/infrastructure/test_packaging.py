"""Tests for distribution metadata and GUI/package boundaries."""

from __future__ import annotations

from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as tomllib


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _pyproject() -> dict[str, object]:
    """Return parsed project metadata."""
    with (PROJECT_ROOT / "pyproject.toml").open("rb") as stream:
        return tomllib.load(stream)


def test_license_and_typed_marker_are_packaged() -> None:
    """Distribution metadata must expose license and PEP 561 contracts."""
    metadata = _pyproject()
    project = metadata["project"]
    assert isinstance(project, dict)
    assert project["license"] == "BSD-3-Clause"
    assert project["license-files"] == ["LICENSE"]
    assert (PROJECT_ROOT / "LICENSE").is_file()
    assert (PROJECT_ROOT / "src/quantas/py.typed").is_file()

    setuptools = metadata["tool"]["setuptools"]
    assert setuptools["package-data"]["quantas"] == ["py.typed"]


def test_matplotlib_is_the_only_supported_plotting_backend() -> None:
    """The scientific package exposes one optional static renderer."""
    metadata = _pyproject()
    extras = metadata["project"]["optional-dependencies"]
    assert "gui" not in extras
    assert extras["plot"] == ["matplotlib>=3.7"]
    flattened = " ".join(item for values in extras.values() for item in values)
    assert "dash" not in flattened.lower()


def test_mypy_checks_the_complete_quantas_package() -> None:
    """Static checking must cover every package module."""
    metadata = _pyproject()
    mypy = metadata["tool"]["mypy"]
    assert mypy["files"] == ["src/quantas"]


def test_typecheck_extra_pins_numpy_stub_generation() -> None:
    """The Python 3.10 type-check baseline uses compatible NumPy stubs."""
    metadata = _pyproject()
    extras = metadata["project"]["optional-dependencies"]
    assert "numpy>=1.23.5,<2.3" in extras["typecheck"]
    assert (PROJECT_ROOT / "requirements" / "typecheck.txt").is_file()
