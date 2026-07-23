#!/usr/bin/env python3
"""Write a JSON inventory of the Quantas source and test tree."""

from __future__ import annotations

import argparse
import ast
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def imports_from(path: Path) -> set[str]:
    """Return absolute imports contained in a Python source file.

    Parameters
    ----------
    path : Path
        Python source file.

    Returns
    -------
    set of str
        Imported module names.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module)
        elif isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
    return imports


def build_snapshot(root: Path) -> dict[str, Any]:
    """Build a project inventory for a Quantas source tree.

    Parameters
    ----------
    root : Path
        Project root containing ``src`` and ``tests``.

    Returns
    -------
    dict
        Serializable architecture inventory.
    """
    package_root = root / "src" / "quantas"
    source_files = sorted(package_root.rglob("*.py"))
    test_files = sorted((root / "tests").rglob("test_*.py"))
    cross_dependencies: set[tuple[str, str]] = set()
    for path in source_files:
        relative = path.relative_to(package_root)
        if len(relative.parts) < 3 or relative.parts[0] != "modules":
            continue
        source_module = relative.parts[1]
        for imported in imports_from(path):
            if not imported.startswith("quantas.modules."):
                continue
            parts = imported.split(".")
            if len(parts) >= 3 and parts[2] != source_module:
                cross_dependencies.add((source_module, parts[2]))

    domains: dict[str, dict[str, int]] = {}
    total_tests = 0
    for path in test_files:
        relative = path.relative_to(root / "tests")
        domain = relative.parts[0]
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        tests_in_file = sum(
            1
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name.startswith("test_")
        )
        record = domains.setdefault(domain, {"files": 0, "tests": 0})
        record["files"] += 1
        record["tests"] += tests_in_file
        total_tests += tests_in_file

    return {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_file_count": len(source_files),
        "test_file_count": len(test_files),
        "static_test_function_count": total_tests,
        "source_packages": sorted(
            str(path.parent.relative_to(package_root))
            for path in source_files
            if path.name == "__init__.py"
        ),
        "test_domains": domains,
        "cross_module_dependencies": [
            {"source": source, "target": target}
            for source, target in sorted(cross_dependencies)
        ],
    }


def main() -> None:
    """Write the requested JSON project inventory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "project_internal/checks/project_stats.json"
        ),
    )
    args = parser.parse_args()
    snapshot = build_snapshot(args.root.resolve())
    output = args.output
    if not output.is_absolute():
        output = args.root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
