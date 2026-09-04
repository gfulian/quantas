"""Source-tree hygiene contracts for the frozen Quantas test suite."""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest


pytestmark = pytest.mark.infrastructure
_TEST_ROOT = Path(__file__).resolve().parents[1]
_MAX_STATIC_NODE_LENGTH = 110
_FORBIDDEN_NAME_TERMS = re.compile(r"(?:^|_)(?:phase\d+|legacy)(?:_|$)")
_NUMBERED_TEST_FILE = re.compile(r"^test_\d+(?:[-_]|$)")


def _test_functions(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return [
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name.startswith("test_")
    ]


def test_test_paths_describe_domains_not_refactoring_history() -> None:
    """Test paths must not encode implementation phases or numeric ordering."""
    violations: list[str] = []
    for path in sorted(_TEST_ROOT.rglob("test_*.py")):
        relative = path.relative_to(_TEST_ROOT).as_posix()
        if _NUMBERED_TEST_FILE.match(path.name):
            violations.append(relative)
        if any(_FORBIDDEN_NAME_TERMS.search(part.lower()) for part in path.parts):
            violations.append(relative)
    assert not violations, "historical test paths: " + ", ".join(violations)


def test_test_names_are_concise_and_descriptive() -> None:
    """Static pytest node identifiers must remain readable in console output."""
    violations: list[str] = []
    for path in sorted(_TEST_ROOT.rglob("test_*.py")):
        relative = path.relative_to(_TEST_ROOT.parent).as_posix()
        for function_name in _test_functions(path):
            nodeid = f"{relative}::{function_name}"
            if _FORBIDDEN_NAME_TERMS.search(function_name.lower()):
                violations.append(f"historical name: {nodeid}")
            if len(nodeid) > _MAX_STATIC_NODE_LENGTH:
                violations.append(
                    f"{len(nodeid)} characters (max {_MAX_STATIC_NODE_LENGTH}): {nodeid}"
                )
    assert not violations, "unreadable test identifiers:\n" + "\n".join(violations)


def test_pytest_uses_importlib_mode_for_readable_module_names() -> None:
    """Repeated semantic filenames are safe only with importlib collection."""
    pyproject = (_TEST_ROOT.parent / "pyproject.toml").read_text(encoding="utf-8")
    assert '"--import-mode=importlib"' in pyproject


def test_developer_tools_use_public_human_readable_names() -> None:
    """The public tools directory must expose only current maintenance names."""
    tools = _TEST_ROOT.parent / "tools"
    names = {path.name for path in tools.glob("*.py")}
    expected = {
        "benchmark_elasticity_sampling.py",
        "benchmark_seismic_vectorization.py",
        "check_architecture.py",
        "check_distribution.py",
        "check_repository.py",
        "compare_qha_workflows.py",
        "project_stats.py",
        "run_tests.py",
        "update_examples_manifest.py",
        "update_scientific_reference.py",
        "update_seismic_reference.py",
    }
    assert names == expected


def test_source_files_use_the_current_header_policy() -> None:
    """First-party modules follow the current encoding and banner policy."""
    project_root = _TEST_ROOT.parent
    source_root = project_root / "src" / "quantas"
    first_party_roots = (
        source_root,
        project_root / "tests",
        project_root / "tools",
        project_root / "examples",
        project_root / "docs" / "tools",
    )
    banner_phrase = "This file is part of the Quantas code."
    paths = sorted(
        {
            path
            for search_root in first_party_roots
            if search_root.is_dir()
            for path in search_root.rglob("*.py")
            if "__pycache__" not in path.parts
        }
    )
    violations: list[str] = []
    for path in paths:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(project_root).as_posix()
        if any(line.strip(" #") == banner_phrase for line in text.splitlines()):
            violations.append(f"historical banner: {relative}")
        in_source = path.is_relative_to(source_root)
        first_line = text.splitlines()[0] if text.splitlines() else ""
        if in_source and first_line != "# -*- coding: utf-8 -*-":
            violations.append(f"missing source encoding declaration: {relative}")
        if not in_source and first_line == "# -*- coding: utf-8 -*-":
            violations.append(f"unnecessary encoding declaration: {relative}")
    assert not violations, "invalid Python headers:\n" + "\n".join(violations)


def test_repository_configuration_files_are_present() -> None:
    """Public source snapshots include portable repository configuration."""
    project_root = _TEST_ROOT.parent
    for name in (".gitignore", ".gitattributes", ".editorconfig"):
        assert (project_root / name).is_file(), f"missing repository file: {name}"
    assert not (project_root / "uv.lock").exists()


def test_windows_docs_build_starts_from_repository_root() -> None:
    """The Windows Sphinx entry point keeps Git discovery at repository root."""
    project_root = _TEST_ROOT.parent
    batch = (project_root / "docs" / "make.bat").read_text(encoding="utf-8")
    assert 'set "ROOT_DIR=%~dp0.."' in batch
    assert 'pushd "%ROOT_DIR%"' in batch
    assert '"docs\\source" "docs\\_build"' in batch
    assert 'pushd "%~dp0"' not in batch


def test_public_source_objects_have_docstrings() -> None:
    """Every public source module, class, function, and method is documented."""
    source_root = _TEST_ROOT.parent / "src" / "quantas"
    violations: list[str] = []
    for path in sorted(source_root.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        relative = path.relative_to(source_root).as_posix()
        if ast.get_docstring(tree) is None:
            violations.append(f"module: {relative}")
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                if not node.name.startswith("_") and ast.get_docstring(node) is None:
                    violations.append(f"function: {relative}::{node.name}")
            elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
                if ast.get_docstring(node) is None:
                    violations.append(f"class: {relative}::{node.name}")
                for member in node.body:
                    if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        if (
                            not member.name.startswith("_")
                            and ast.get_docstring(member) is None
                        ):
                            violations.append(
                                f"method: {relative}::{node.name}.{member.name}"
                            )
    assert not violations, "missing public docstrings:\n" + "\n".join(violations)
