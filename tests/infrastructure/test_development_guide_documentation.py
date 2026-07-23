"""Structural coverage tests for the Quantas Development Guide."""

from __future__ import annotations

from pathlib import Path


ROOT = Path("docs/source/developer")


EXPECTED_PAGES = {
    "getting_started.rst",
    "architecture.rst",
    "module_anatomy.rst",
    "numerical_methods.rst",
    "public_api.rst",
    "change_recipes.rst",
    "extending.rst",
    "interfaces.rst",
    "events.rst",
    "persistence.rst",
    "rendering_frontends.rst",
    "citation_registry.rst",
    "testing.rst",
    "documentation.rst",
    "thermoelastic_architecture.rst",
    "packaging.rst",
    "review_checklist.rst",
}


def _text(name: str) -> str:
    return (ROOT / name).read_text(encoding="utf-8")


def test_development_guide_has_no_placeholder_pages() -> None:
    """All approved developer topics are implemented rather than placeholders."""
    assert EXPECTED_PAGES <= {path.name for path in ROOT.glob("*.rst")}
    for name in EXPECTED_PAGES:
        content = _text(name)
        assert "Work in progress" not in content
        assert len(content.splitlines()) >= 100


def test_architecture_topics_are_covered() -> None:
    """The guide preserves the central layer and frontend rules."""
    content = _text("architecture.rst")
    for phrase in (
        "Scientific results must not depend on the frontend",
        "Allowed dependency direction",
        "Active objects and passive data",
        "Scientific module isolation",
        "Precision and units",
        "Frontend ownership",
    ):
        assert phrase in content


def test_module_and_change_guides_are_practical() -> None:
    """Contributors receive complete module and cross-cutting procedures."""
    module = _text("module_anatomy.rst")
    for phrase in (
        "Passive contracts",
        "Analysis functions",
        "Calculator",
        "Result envelope",
        "Persistence adapter",
        "Internal and public facades",
        "CLI adapter",
    ):
        assert phrase in module

    recipes = _text("change_recipes.rst")
    for phrase in (
        "Add a scientific option",
        "Add a result field",
        "Add a CLI option",
        "Add a plot",
        "Add an HDF5 field or change a schema",
        "Add a public API symbol",
        "Add an external-code parser",
        "Change a numerical default",
    ):
        assert phrase in recipes


def test_extension_starts_from_science() -> None:
    """A new module is developed from physics outward to frontends."""
    content = _text("extending.rst")
    assert "Do not begin by writing a Click command" in content
    for phase in range(1, 13):
        assert f"Phase {phase}:" in content
    assert "CLI/API equivalence" in content
    assert "scientific regression" in content


def test_cross_cutting_contracts_are_documented() -> None:
    """Persistence, events, rendering, interfaces, and citations are explicit."""
    expectations = {
        "events.rst": (
            "Persistent and operational events",
            "PROGRESS",
            "EventRecord",
        ),
        "persistence.rst": (
            "Shared envelope",
            "Package version versus schema version",
            "Migration policy",
            "Round-trip testing",
        ),
        "rendering_frontends.rst": (
            "Report pipeline",
            "Plot pipeline",
            "Matplotlib backend",
            "Future GUI requirements",
        ),
        "interfaces.rst": (
            "File recognition and completion",
            "Units and conventions",
            "Fixture strategy",
        ),
        "citation_registry.rst": (
            "Module and method sets",
            "Scientific-background pages",
            "DOI",
        ),
    }
    for name, phrases in expectations.items():
        content = _text(name)
        for phrase in phrases:
            assert phrase in content


def test_project_tools_are_documented() -> None:
    """The guide points contributors to the deterministic project workflow."""
    testing = _text("testing.rst")
    assert "python tools/run_tests.py all -- -q" in testing
    assert "Characterization versus validation" in testing
    assert "python tools/check_distribution.py dist" in testing

    packaging = _text("packaging.rst")
    assert "python -m build" in packaging
    assert "python -m twine check dist/*" in packaging
    assert "py.typed" in packaging


def test_review_checklist_covers_lifecycle() -> None:
    """The final checklist spans science, architecture, docs, tests, and packages."""
    content = _text("review_checklist.rst")
    for heading in (
        "Scientific contract",
        "Architecture",
        "Results and persistence",
        "Public API and CLI",
        "Tests",
        "Documentation and citations",
        "Packaging",
    ):
        assert heading in content
