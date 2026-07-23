#!/usr/bin/env python3
"""Check the Quantas package architecture without changing source behavior.

The check is intentionally static: it inspects imports, file sizes, class and
function sizes, cross-module dependencies, direct UI dependencies in non-UI
layers, print usage, HDF5 concentration, duplicated helper names, and test
coverage distribution.  The output is meant to guide refactoring discussions;
it is not a strict quality gate by default.
"""

from __future__ import annotations

import argparse
import ast
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

UI_IMPORTS = ("click", "dash", "matplotlib")
TABLE_RENDER_IMPORTS = ("rich", "prettytable", "tabulate", "pandas")
HDF5_IMPORTS = ("h5py",)
CORE_ALLOWED_UI_IMPORTS: tuple[str, ...] = ()
LONG_FILE_LINES = 800
VERY_LONG_FILE_LINES = 1200
LONG_FUNCTION_LINES = 120
LONG_CLASS_LINES = 300


@dataclass(slots=True)
class PythonObjectMetric:
    """Static metric for a function or class.

    Parameters
    ----------
    path : str
        Project-relative source path.
    name : str
        Object name.
    kind : str
        Object kind, usually ``function`` or ``class``.
    start_line : int
        First source line.
    end_line : int
        Last source line.
    line_count : int
        Object size in lines.
    """

    path: str
    name: str
    kind: str
    start_line: int
    end_line: int
    line_count: int


@dataclass(slots=True)
class FileMetric:
    """Static metric for one Python source file.

    Parameters
    ----------
    path : str
        Project-relative source path.
    line_count : int
        Source size in lines.
    functions : int
        Number of top-level and nested functions.
    classes : int
        Number of classes.
    imports : list[str]
        Imported top-level or absolute module names.
    """

    path: str
    line_count: int
    functions: int
    classes: int
    imports: list[str]


@dataclass(slots=True)
class AuditIssue:
    """Architecture check issue or observation.

    Parameters
    ----------
    severity : str
        One of ``OK``, ``INFO``, ``WARN`` or ``FAIL``.
    category : str
        Short category name.
    message : str
        Human-readable description.
    path : str | None, optional
        Project-relative path, when relevant.
    line : int | None, optional
        Source line, when relevant.
    """

    severity: str
    category: str
    message: str
    path: str | None = None
    line: int | None = None


def relpath(path: Path, root: Path) -> str:
    """Return a POSIX project-relative path.

    Parameters
    ----------
    path : Path
        Path to normalize.
    root : Path
        Project root.

    Returns
    -------
    str
        POSIX relative path.
    """

    return path.relative_to(root).as_posix()


def read_text(path: Path) -> str:
    """Read source text using UTF-8.

    Parameters
    ----------
    path : Path
        File to read.

    Returns
    -------
    str
        File contents.
    """

    return path.read_text(encoding="utf-8")


def parse_source(path: Path) -> ast.Module:
    """Parse a Python source file.

    Parameters
    ----------
    path : Path
        Source file to parse.

    Returns
    -------
    ast.Module
        Parsed syntax tree.
    """

    return ast.parse(read_text(path), filename=str(path))


def iter_imports(tree: ast.AST) -> set[str]:
    """Collect imported module names from an AST.

    Parameters
    ----------
    tree : ast.AST
        Parsed source tree.

    Returns
    -------
    set[str]
        Imported module names as written in import statements.
    """

    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module)
    return imports


def has_import(imports: Iterable[str], prefixes: Iterable[str]) -> bool:
    """Check whether imports contain any requested prefix.

    Parameters
    ----------
    imports : Iterable[str]
        Imported module names.
    prefixes : Iterable[str]
        Prefixes to match.

    Returns
    -------
    bool
        ``True`` if at least one import starts with a prefix.
    """

    for imported in imports:
        for prefix in prefixes:
            if imported == prefix or imported.startswith(prefix + "."):
                return True
    return False


def object_metric(
    path: Path, root: Path, node: ast.AST, kind: str
) -> PythonObjectMetric | None:
    """Build a metric record for a class or function AST node.

    Parameters
    ----------
    path : Path
        Source file path.
    root : Path
        Project root.
    node : ast.AST
        Class or function node.
    kind : str
        Object kind.

    Returns
    -------
    PythonObjectMetric | None
        Metric record, or ``None`` when line numbers are unavailable.
    """

    start = getattr(node, "lineno", None)
    end = getattr(node, "end_lineno", None)
    name = getattr(node, "name", None)
    if start is None or end is None or name is None:
        return None
    return PythonObjectMetric(
        path=relpath(path, root),
        name=name,
        kind=kind,
        start_line=start,
        end_line=end,
        line_count=end - start + 1,
    )


def file_layer(path: str) -> str:
    """Classify a source path into a high-level package layer.

    Parameters
    ----------
    path : str
        Project-relative path.

    Returns
    -------
    str
        Layer name.
    """

    parts = Path(path).parts
    if len(parts) < 3 or parts[0] != "src" or parts[1] != "quantas":
        return "other"
    if len(parts) == 3:
        return "root"
    return parts[2]


def module_name_for_path(path: str) -> str | None:
    """Return the Quantas workflow module name for a source path.

    Parameters
    ----------
    path : str
        Project-relative path.

    Returns
    -------
    str | None
        Workflow module name, or ``None`` for non-module paths.
    """

    parts = Path(path).parts
    if len(parts) >= 4 and parts[0:3] == ("src", "quantas", "modules"):
        return parts[3]
    return None


def summarize_tests(root: Path) -> dict[str, dict[str, int]]:
    """Summarize static test counts by top-level test domain.

    Parameters
    ----------
    root : Path
        Project root.

    Returns
    -------
    dict[str, dict[str, int]]
        Domain summaries with file and test counts.
    """

    tests_root = root / "tests"
    domains: dict[str, dict[str, int]] = {}
    if not tests_root.exists():
        return domains
    for path in sorted(tests_root.rglob("test_*.py")):
        relative = path.relative_to(tests_root)
        domain = relative.parts[0]
        tree = parse_source(path)
        count = sum(
            1
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name.startswith("test_")
        )
        record = domains.setdefault(domain, {"files": 0, "tests": 0})
        record["files"] += 1
        record["tests"] += count
    return domains


def analyze_source(root: Path) -> dict[str, Any]:
    """Run the static architecture check.

    Parameters
    ----------
    root : Path
        Project root.

    Returns
    -------
    dict[str, Any]
        Serializable architecture-check result.
    """

    package_root = root / "src" / "quantas"
    source_files = sorted(
        path for path in package_root.rglob("*.py") if "__pycache__" not in path.parts
    )
    metrics: list[FileMetric] = []
    functions: list[PythonObjectMetric] = []
    classes: list[PythonObjectMetric] = []
    issues: list[AuditIssue] = []
    imports_by_file: dict[str, set[str]] = {}
    helper_names: dict[str, list[str]] = defaultdict(list)
    calls_by_file: dict[str, list[tuple[str, int]]] = defaultdict(list)

    for path in source_files:
        text = read_text(path)
        tree = parse_source(path)
        imports = iter_imports(tree)
        relative = relpath(path, root)
        imports_by_file[relative] = imports
        line_count = len(text.splitlines())
        function_count = 0
        class_count = 0
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                function_count += 1
                metric = object_metric(path, root, node, "function")
                if metric:
                    functions.append(metric)
                    helper_names[node.name].append(relative)
            elif isinstance(node, ast.ClassDef):
                class_count += 1
                metric = object_metric(path, root, node, "class")
                if metric:
                    classes.append(metric)
                    helper_names[node.name].append(relative)
            elif isinstance(node, ast.Call):
                func = node.func
                name: str | None = None
                if isinstance(func, ast.Name):
                    name = func.id
                elif isinstance(func, ast.Attribute):
                    name = func.attr
                if name:
                    calls_by_file[relative].append((name, getattr(node, "lineno", 0)))
        metrics.append(
            FileMetric(
                path=relative,
                line_count=line_count,
                functions=function_count,
                classes=class_count,
                imports=sorted(imports),
            )
        )

        if line_count >= VERY_LONG_FILE_LINES:
            issues.append(
                AuditIssue(
                    "WARN",
                    "large-file",
                    f"Very large source file ({line_count} lines); strong candidate for decomposition.",
                    relative,
                )
            )
        elif line_count >= LONG_FILE_LINES:
            issues.append(
                AuditIssue(
                    "INFO",
                    "large-file",
                    f"Large source file ({line_count} lines); monitor or split when touched.",
                    relative,
                )
            )

        layer = file_layer(relative)
        if layer == "core" and has_import(imports, UI_IMPORTS):
            issues.append(
                AuditIssue(
                    "FAIL",
                    "ui-import-in-core",
                    "Core layer imports a UI/rendering package.",
                    relative,
                )
            )
        if layer in {"core", "models", "modules", "io", "interfaces"} and has_import(
            imports, ("click", "dash")
        ):
            issues.append(
                AuditIssue(
                    "FAIL",
                    "frontend-import",
                    "Non-frontend layer imports Click or Dash.",
                    relative,
                )
            )
        if layer in {
            "core",
            "models",
            "modules",
            "io",
            "interfaces",
            "cli",
        } and has_import(imports, ("matplotlib",)):
            issues.append(
                AuditIssue(
                    "FAIL",
                    "renderer-import",
                    "Non-renderer layer imports Matplotlib.",
                    relative,
                )
            )
        if layer == "core" and has_import(imports, HDF5_IMPORTS):
            issues.append(
                AuditIssue(
                    "FAIL",
                    "hdf5-import-in-core",
                    "Core physics/math layer imports HDF5 directly.",
                    relative,
                )
            )

        for name, lineno in calls_by_file[relative]:
            if name == "print" and layer not in {"cli"}:
                issues.append(
                    AuditIssue(
                        "WARN",
                        "print-call",
                        "Direct print() outside cli layer; prefer events, logging, or explicit renderer output.",
                        relative,
                        lineno,
                    )
                )

    for func in functions:
        if func.line_count >= LONG_FUNCTION_LINES:
            issues.append(
                AuditIssue(
                    "WARN",
                    "long-function",
                    f"Long function/method {func.name!r} ({func.line_count} lines).",
                    func.path,
                    func.start_line,
                )
            )
    for cls in classes:
        if cls.line_count >= LONG_CLASS_LINES:
            issues.append(
                AuditIssue(
                    "WARN",
                    "long-class",
                    f"Large class {cls.name!r} ({cls.line_count} lines).",
                    cls.path,
                    cls.start_line,
                )
            )

    cross_module_dependencies: list[dict[str, str]] = []
    for relative, imports in imports_by_file.items():
        source_module = module_name_for_path(relative)
        if source_module is None:
            continue
        for imported in imports:
            if not imported.startswith("quantas.modules."):
                continue
            parts = imported.split(".")
            if len(parts) >= 3:
                target_module = parts[2]
                if target_module != source_module:
                    cross_module_dependencies.append(
                        {
                            "source": source_module,
                            "target": target_module,
                            "path": relative,
                        }
                    )
                    issues.append(
                        AuditIssue(
                            "WARN",
                            "module-cross-dependency",
                            f"Workflow module imports another workflow module: {source_module} -> {target_module}.",
                            relative,
                        )
                    )

    hdf5_files = [
        metric.path
        for metric in metrics
        if has_import(imports_by_file[metric.path], HDF5_IMPORTS)
    ]
    duplicated_helpers = [
        {"name": name, "count": len(paths), "paths": sorted(set(paths))}
        for name, paths in helper_names.items()
        if len(set(paths)) >= 3
        and name.startswith("_")
        and not (name.startswith("__") and name.endswith("__"))
        and name not in {"__post_init__"}
    ]
    duplicated_helpers.sort(key=lambda item: (-item["count"], item["name"]))

    module_contracts: dict[str, bool] = {}
    modules_root = package_root / "modules"
    if modules_root.exists():
        for module_dir in sorted(
            p for p in modules_root.iterdir() if p.is_dir() and p.name != "__pycache__"
        ):
            init_file = module_dir / "__init__.py"
            module_contracts[module_dir.name] = (
                init_file.exists() and "MODULE_CONTRACT" in read_text(init_file)
            )
            if not module_contracts[module_dir.name]:
                issues.append(
                    AuditIssue(
                        "INFO",
                        "module-contract",
                        "Workflow module has no MODULE_CONTRACT declaration in __init__.py.",
                        relpath(init_file if init_file.exists() else module_dir, root),
                    )
                )

    severity_counts = Counter(issue.severity for issue in issues)
    layer_counts = Counter(file_layer(metric.path) for metric in metrics)
    source_lines_by_layer: dict[str, int] = defaultdict(int)
    for metric in metrics:
        source_lines_by_layer[file_layer(metric.path)] += metric.line_count

    return {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "summary": {
            "source_files": len(source_files),
            "source_lines": sum(metric.line_count for metric in metrics),
            "functions": len(functions),
            "classes": len(classes),
            "issues": dict(severity_counts),
            "layers": dict(layer_counts),
            "source_lines_by_layer": dict(sorted(source_lines_by_layer.items())),
        },
        "files": [asdict(metric) for metric in metrics],
        "largest_files": [
            asdict(metric)
            for metric in sorted(metrics, key=lambda m: m.line_count, reverse=True)[:25]
        ],
        "longest_functions": [
            asdict(metric)
            for metric in sorted(functions, key=lambda m: m.line_count, reverse=True)[
                :25
            ]
        ],
        "largest_classes": [
            asdict(metric)
            for metric in sorted(classes, key=lambda m: m.line_count, reverse=True)[:25]
        ],
        "hdf5_files": sorted(hdf5_files),
        "duplicated_private_helpers": duplicated_helpers[:50],
        "cross_module_dependencies": cross_module_dependencies,
        "module_contracts": module_contracts,
        "test_domains": summarize_tests(root),
        "issues": [asdict(issue) for issue in issues],
    }


def severity_rank(issue: dict[str, Any]) -> int:
    """Return sort rank for issue severity.

    Parameters
    ----------
    issue : dict[str, Any]
        Serialized issue record.

    Returns
    -------
    int
        Numeric rank.
    """

    return {"FAIL": 0, "WARN": 1, "INFO": 2, "OK": 3}.get(issue["severity"], 4)


def write_markdown(audit: dict[str, Any], output: Path) -> None:
    """Write a Markdown architecture check report.

    Parameters
    ----------
    audit : dict[str, Any]
        Audit result.
    output : Path
        Destination Markdown path.
    """

    summary = audit["summary"]
    lines: list[str] = []
    lines.append("# Quantas architecture check\n")
    lines.append(f"Generated at: `{audit['created_at']}`\n")
    lines.append("## Executive summary\n")
    lines.append(f"- Source files: **{summary['source_files']}**")
    lines.append(f"- Source lines: **{summary['source_lines']}**")
    lines.append(f"- Functions/methods: **{summary['functions']}**")
    lines.append(f"- Classes: **{summary['classes']}**")
    lines.append(f"- Issue counts: `{summary['issues']}`")
    lines.append("")
    lines.append("## Layer sizes\n")
    lines.append("| Layer | Files | Lines |")
    lines.append("|---|---:|---:|")
    for layer, files in sorted(summary["layers"].items()):
        lines.append(
            f"| `{layer}` | {files} | {summary['source_lines_by_layer'].get(layer, 0)} |"
        )
    lines.append("")

    lines.append("## Largest files\n")
    lines.append("| Lines | File |")
    lines.append("|---:|---|")
    for item in audit["largest_files"][:15]:
        lines.append(f"| {item['line_count']} | `{item['path']}` |")
    lines.append("")

    lines.append("## Longest functions and methods\n")
    lines.append("| Lines | Function | File | Line |")
    lines.append("|---:|---|---|---:|")
    for item in audit["longest_functions"][:15]:
        lines.append(
            f"| {item['line_count']} | `{item['name']}` | `{item['path']}` | {item['start_line']} |"
        )
    lines.append("")

    lines.append("## Largest classes\n")
    lines.append("| Lines | Class | File | Line |")
    lines.append("|---:|---|---|---:|")
    for item in audit["largest_classes"][:15]:
        lines.append(
            f"| {item['line_count']} | `{item['name']}` | `{item['path']}` | {item['start_line']} |"
        )
    lines.append("")

    lines.append("## HDF5 concentration\n")
    lines.append("Files importing `h5py`:")
    for path in audit["hdf5_files"]:
        lines.append(f"- `{path}`")
    lines.append("")

    lines.append("## Cross-module workflow dependencies\n")
    if audit["cross_module_dependencies"]:
        lines.append("| Source | Target | File |")
        lines.append("|---|---|---|")
        for item in audit["cross_module_dependencies"]:
            lines.append(
                f"| `{item['source']}` | `{item['target']}` | `{item['path']}` |"
            )
    else:
        lines.append("No direct `quantas.modules.<other>` dependencies found.")
    lines.append("")

    lines.append("## Module contracts\n")
    lines.append("| Module | MODULE_CONTRACT present |")
    lines.append("|---|---:|")
    for name, present in sorted(audit["module_contracts"].items()):
        lines.append(f"| `{name}` | {'yes' if present else 'no'} |")
    lines.append("")

    lines.append("## Duplicated private helper names\n")
    lines.append(
        "These are not automatically wrong, but repeated names often indicate reusable mechanics."
    )
    lines.append("| Helper | Count | Files |")
    lines.append("|---|---:|---|")
    for item in audit["duplicated_private_helpers"][:20]:
        paths = "<br>".join(f"`{path}`" for path in item["paths"][:8])
        extra = len(item["paths"]) - 8
        if extra > 0:
            paths += f"<br>... +{extra} more"
        lines.append(f"| `{item['name']}` | {item['count']} | {paths} |")
    lines.append("")

    lines.append("## Test inventory\n")
    lines.append("| Domain | Files | Static test functions |")
    lines.append("|---|---:|---:|")
    for domain, record in sorted(audit["test_domains"].items()):
        lines.append(f"| `{domain}` | {record['files']} | {record['tests']} |")
    lines.append("")

    lines.append("## Issues\n")
    lines.append("| Severity | Category | Location | Message |")
    lines.append("|---|---|---|---|")
    for issue in sorted(
        audit["issues"],
        key=lambda i: (severity_rank(i), i["category"], i.get("path") or ""),
    ):
        location = issue.get("path") or ""
        if issue.get("line"):
            location += f":{issue['line']}"
        lines.append(
            f"| {issue['severity']} | `{issue['category']}` | `{location}` | {issue['message']} |"
        )
    lines.append("")

    lines.append("## Recommended interpretation\n")
    lines.append(
        "This check is a map, not a refactoring plan. FAIL entries identify architectural boundaries that should be corrected or explicitly justified. WARN entries identify good candidates for decomposition or shared helpers. INFO entries are mostly planning signals."
    )
    lines.append("")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run the architecture check command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--json",
        type=Path,
        default=Path(
            "project_internal/checks/architecture.json"
        ),
    )
    parser.add_argument(
        "--markdown",
        type=Path,
        default=Path(
            "project_internal/checks/architecture.md"
        ),
    )
    args = parser.parse_args()
    root = args.root.resolve()
    audit = analyze_source(root)

    json_output = args.json if args.json.is_absolute() else root / args.json
    markdown_output = (
        args.markdown if args.markdown.is_absolute() else root / args.markdown
    )
    json_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    write_markdown(audit, markdown_output)

    print(f"JSON report: {json_output}")
    print(f"Markdown report: {markdown_output}")
    print(f"Issue counts: {audit['summary']['issues']}")


if __name__ == "__main__":
    main()
