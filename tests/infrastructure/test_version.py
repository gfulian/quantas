"""Tests for the authoritative Quantas package version."""

from __future__ import annotations

from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as tomllib

from click.testing import CliRunner

import quantas
from quantas._version import __version__ as authoritative_version
from quantas.cli.main import main
from quantas.models import ResultMetadata


EXPECTED_VERSION = "2.0.0b13"


def test_public_version_uses_authoritative_version():
    """The public package version should expose the authoritative value."""
    assert authoritative_version == EXPECTED_VERSION
    assert quantas.__version__ == authoritative_version


def test_result_metadata_uses_package_version_by_default():
    """New result metadata should record the current Quantas version."""
    assert ResultMetadata().version == authoritative_version


def test_cli_version_reports_package_version():
    """The top-level CLI should print the current Quantas version."""
    result = CliRunner().invoke(main, ["-v"])

    assert result.exit_code == 0
    assert f"v{authoritative_version}" in result.output


def test_build_metadata_reads_authoritative_version_module():
    """Project metadata should obtain its version from ``quantas._version``."""
    root = Path(__file__).resolve().parents[2]
    configuration = tomllib.loads((root / "pyproject.toml").read_text())

    assert "version" not in configuration["project"]
    assert "version" in configuration["project"]["dynamic"]
    assert configuration["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "quantas._version.__version__"
    }


def test_project_readme_is_declared_and_present() -> None:
    """The source distribution must expose an explicit project README."""
    root = Path(__file__).parents[2]
    metadata = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    assert metadata["project"]["readme"] == "README.md"
    readme = root / metadata["project"]["readme"]
    assert readme.is_file()
    assert readme.read_text(encoding="utf-8").startswith("# Quantas")


def test_root_namespace_exposes_only_the_version() -> None:
    """Package metadata should not become an accidental public API."""
    assert quantas.__all__ == ["__version__"]
    assert not hasattr(quantas, "__author__")
    assert not hasattr(quantas, "__citation_key__")


def test_release_metadata_matches_authoritative_version() -> None:
    """Current release-facing metadata must track the source version."""
    root = Path(__file__).resolve().parents[2]
    version = authoritative_version

    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    project_state = (root / "PROJECT_STATE.md").read_text(encoding="utf-8")
    roadmap = (root / "ROADMAP.md").read_text(encoding="utf-8")
    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    lifecycle_tool = (
        root / "tools" / "validate_public_lifecycle_api.ps1"
    ).read_text(encoding="utf-8")
    tools_readme = (root / "tools" / "README.md").read_text(encoding="utf-8")

    assert f"version: {version}" in citation
    assert f"| Current development version | `{version}` |" in project_state
    assert f"Development is now on ``{version}``" in roadmap
    assert f"## [{version}] - Unreleased" in changelog
    assert version in lifecycle_tool
    assert version in tools_readme
    assert not (root / "tools" / "validate_public_lifecycle_api_step6.ps1").exists()
