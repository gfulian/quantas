"""Dependency-boundary tests for the Quantas package architecture."""

from __future__ import annotations

import ast
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[2] / "src" / "quantas"


def _absolute_imports(path: Path) -> set[str]:
    """Return absolute modules imported by one Python source file."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module)
        elif isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
    return imports


def _contains_print_call(path: Path) -> bool:
    """Return whether a source file calls the built-in ``print`` function."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "print"
        for node in ast.walk(tree)
    )


def test_core_does_not_depend_on_modules_interfaces_or_user_interfaces() -> None:
    """The shared numerical core remains independent from higher layers."""
    forbidden = (
        "quantas.modules",
        "quantas.interfaces",
        "quantas.cli",
        "quantas.gui",
    )
    violations: list[str] = []
    for path in sorted((PACKAGE_ROOT / "core").rglob("*.py")):
        for imported in _absolute_imports(path):
            if imported.startswith(forbidden):
                violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> {imported}")
    assert not violations, "core dependency violations:\n" + "\n".join(violations)


def test_models_do_not_depend_on_scientific_modules() -> None:
    """Passive shared models must not depend on scientific workflows."""
    violations: list[str] = []
    for path in sorted((PACKAGE_ROOT / "models").rglob("*.py")):
        for imported in _absolute_imports(path):
            if imported.startswith("quantas.modules"):
                violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> {imported}")
    assert not violations, "model dependency violations:\n" + "\n".join(violations)


def test_interfaces_do_not_depend_on_scientific_modules() -> None:
    """Code-specific parsers remain independent from module workflows."""
    violations: list[str] = []
    for path in sorted((PACKAGE_ROOT / "interfaces").rglob("*.py")):
        for imported in _absolute_imports(path):
            if imported.startswith("quantas.modules"):
                violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> {imported}")
    assert not violations, "interface dependency violations:\n" + "\n".join(violations)


def test_current_cross_module_dependencies_are_acyclic_and_explicit() -> None:
    """Scientific workflows do not import one another."""
    allowed: set[tuple[str, str]] = set()
    dependencies: set[tuple[str, str]] = set()
    module_root = PACKAGE_ROOT / "modules"
    for path in sorted(module_root.rglob("*.py")):
        relative = path.relative_to(module_root)
        if len(relative.parts) < 2:
            continue
        source = relative.parts[0]
        for imported in _absolute_imports(path):
            prefix = "quantas.modules."
            if not imported.startswith(prefix):
                continue
            parts = imported.split(".")
            if len(parts) >= 3 and parts[2] != source:
                dependencies.add((source, parts[2]))

    assert dependencies <= allowed
    for source, target in dependencies:
        assert (target, source) not in dependencies


def test_source_tree_does_not_depend_on_numba() -> None:
    """Numba and its removed backend options must not return."""
    violations: list[str] = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        for imported in _absolute_imports(path):
            if imported == "numba" or imported.startswith("numba."):
                violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> {imported}")
        source = path.read_text(encoding="utf-8")
        if "use_numba" in source or "use-numba" in source:
            violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> numba option")
    assert not violations, "Numba dependency violations:\n" + "\n".join(violations)


def test_report_table_has_one_authoritative_definition() -> None:
    """Only the frontend-neutral models package defines ReportTable."""
    definitions: list[Path] = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if any(
            isinstance(node, ast.ClassDef) and node.name == "ReportTable"
            for node in ast.walk(tree)
        ):
            definitions.append(path.relative_to(PACKAGE_ROOT))

    assert definitions == [Path("models/report.py")]


def test_plot_specs_have_one_authoritative_model_module() -> None:
    """Neutral plot specification classes are defined only in models.plot."""
    names = {
        "PlotAxis",
        "PlotSeriesStyle",
        "PlotSeries",
        "LinePlotSpec",
        "ContourPlotSpec",
        "PolarPlotPanel",
        "PolarPlotSpec",
        "SurfaceStyle",
        "SurfaceLayer",
        "SurfacePlotSpec",
        "PlotCollection",
        "SphericalMarker",
        "AxisFieldLayer",
        "SphericalMapSpec",
    }
    definitions: dict[str, list[Path]] = {name: [] for name in names}
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name in definitions:
                definitions[node.name].append(path.relative_to(PACKAGE_ROOT))

    assert all(paths == [Path("models/plot.py")] for paths in definitions.values())


def test_matplotlib_is_confined_to_the_plotting_backend() -> None:
    """Scientific modules and neutral models must not import Matplotlib."""
    violations: list[str] = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        relative = path.relative_to(PACKAGE_ROOT)
        is_backend = relative == Path("renderers/plots/__init__.py") or relative.parts[
            :3
        ] == (
            "renderers",
            "plots",
            "matplotlib",
        )
        for imported in _absolute_imports(path):
            if imported == "matplotlib" or imported.startswith("matplotlib."):
                if not is_backend:
                    violations.append(f"{relative} -> {imported}")
    assert not violations, "Matplotlib dependency violations:\n" + "\n".join(violations)


def test_module_plot_builders_do_not_depend_on_concrete_renderers() -> None:
    """Module plot packages prepare specs without importing renderer packages."""
    violations: list[str] = []
    for path in sorted((PACKAGE_ROOT / "modules").glob("*/plot/*.py")):
        for imported in _absolute_imports(path):
            if imported.startswith("quantas.renderers.plots") or imported.startswith(
                "matplotlib"
            ):
                violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> {imported}")
    assert not violations, "plot-builder dependency violations:\n" + "\n".join(
        violations
    )


def test_obsolete_module_plotting_paths_are_removed() -> None:
    """Historical module-owned Matplotlib implementations must not return."""
    obsolete = [
        PACKAGE_ROOT / "modules" / module / "plot" / "plotting.py"
        for module in ("ha", "qha", "elasticity")
    ]
    existing = [path.relative_to(PACKAGE_ROOT) for path in obsolete if path.exists()]
    assert not existing, "obsolete plotting modules found:\n" + "\n".join(
        str(path) for path in existing
    )


def test_elasticity_core_is_a_frontend_neutral_scientific_domain() -> None:
    """Elasticity physics must remain independent from workflows and frontends."""
    root = PACKAGE_ROOT / "core" / "physics" / "elasticity"
    forbidden = (
        "click",
        "matplotlib",
        "dash",
        "quantas.cli",
        "quantas.renderers.plots",
        "quantas.modules",
    )
    violations: list[str] = []
    for path in sorted(root.rglob("*.py")):
        for imported in _absolute_imports(path):
            if imported.startswith(forbidden):
                violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> {imported}")
        if _contains_print_call(path):
            violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> print")
    assert not violations, "elasticity core boundary violations:\n" + "\n".join(
        violations
    )


def test_elasticity_workflow_analysis_does_not_depend_on_frontends() -> None:
    """High-level elasticity analysis remains usable without CLI or plotting."""
    path = PACKAGE_ROOT / "modules" / "elasticity" / "analysis.py"
    forbidden = (
        "click",
        "matplotlib",
        "dash",
        "quantas.cli",
        "quantas.renderers.plots",
        "quantas.modules.elasticity.cli",
        "quantas.modules.elasticity.plot",
    )
    violations = [
        f"{path.relative_to(PACKAGE_ROOT)} -> {imported}"
        for imported in _absolute_imports(path)
        if imported.startswith(forbidden)
    ]
    if _contains_print_call(path):
        violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> print")
    assert not violations, "elasticity workflow boundary violations:\n" + "\n".join(
        violations
    )


def test_elasticity_core_has_separated_scientific_responsibilities() -> None:
    """The extracted elasticity domain exposes the agreed module layout."""
    root = PACKAGE_ROOT / "core" / "physics" / "elasticity"
    expected = {
        "__init__.py",
        "adiabatic.py",
        "tensor.py",
        "symmetry.py",
        "stability.py",
        "surface.py",
        "averages.py",
        "conventions.py",
        "directional.py",
        "extrema.py",
        "quasistatic.py",
        "sampling.py",
        "validation.py",
    }
    assert {path.name for path in root.glob("*.py")} == expected


def test_elasticity_contains_no_seismic_physics() -> None:
    """All acoustic-wave physics is reserved for the seismic domain."""
    root = PACKAGE_ROOT / "core" / "physics" / "elasticity"
    forbidden_terms = {
        "christoffel",
        "phase_velocity",
        "group_velocity",
        "polarization",
        "longitudinal_wave",
        "shear_wave_1",
        "shear_wave_2",
        "isotropic_seismic",
        "seismic_velocity",
    }
    violations: list[str] = []
    for path in sorted(root.rglob("*.py")):
        source = path.read_text(encoding="utf-8").lower()
        for term in forbidden_terms:
            if term in source:
                violations.append(f"{path.relative_to(PACKAGE_ROOT)} -> {term}")
    assert not violations, "seismic code found in elasticity:\n" + "\n".join(violations)


def test_quasistatic_elasticity_exposes_frontend_neutral_relations() -> None:
    """The activated quasi-static core exposes only numerical relations."""
    path = PACKAGE_ROOT / "core" / "physics" / "elasticity" / "quasistatic.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    functions = {
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assert functions == {
        "cold_finite_strain_component",
        "cold_finite_strain_component_jacobian",
        "cold_finite_strain_stiffness",
        "eulerian_finite_strain",
        "wallace_hydrostatic_delta_voigt",
    }
    forbidden = (
        "click",
        "rich",
        "matplotlib",
        "quantas.modules",
        "quantas.interfaces",
    )
    imported = _absolute_imports(path)
    assert not {name for name in imported if name.startswith(forbidden)}
    assert not _contains_print_call(path)


def test_obsolete_soec_namespace_and_module_owned_core_are_removed() -> None:
    """Quantas 2.0 does not retain old SOEC packages or numerical shims."""
    obsolete = [
        PACKAGE_ROOT / "modules" / "soec",
        PACKAGE_ROOT / "modules" / "elasticity" / "core",
        PACKAGE_ROOT / "modules" / "elasticity" / "optimization.py",
        PACKAGE_ROOT / "interfaces" / "crystal" / "soec.py",
        PACKAGE_ROOT / "interfaces" / "vasp" / "soec.py",
    ]
    existing = [path.relative_to(PACKAGE_ROOT) for path in obsolete if path.exists()]
    assert not existing, "obsolete SOEC paths found:\n" + "\n".join(
        str(path) for path in existing
    )


def test_old_soec_public_class_names_are_absent_from_source() -> None:
    """No backward-compatibility aliases preserve the old public class names."""
    forbidden_names = {
        "SOEC",
        "SOECOrtho",
        "SOECInput",
        "SOECOptions",
        "SOECResult",
        "SOECCalculator",
        "SOECInputFileReader",
    }
    definitions: list[str] = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if (
                isinstance(node, (ast.ClassDef, ast.FunctionDef))
                and node.name in forbidden_names
            ):
                definitions.append(f"{path.relative_to(PACKAGE_ROOT)} -> {node.name}")
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id in forbidden_names:
                        definitions.append(
                            f"{path.relative_to(PACKAGE_ROOT)} -> {target.id}"
                        )
    assert not definitions, "old SOEC API aliases found:\n" + "\n".join(definitions)
