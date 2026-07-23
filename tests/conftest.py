"""Shared pytest configuration for Quantas tests."""

from __future__ import annotations

import os
from collections import Counter

import pytest

from tests.suite_registry import (
    STAGED_ALL_SUITES,
    SUITE_BY_NAME,
    classify_nodeid,
    domain_for_nodeid,
    list_suite_rows,
)

os.environ.setdefault("MPLBACKEND", "Agg")

_DOMAIN_ORDER = (
    "infrastructure",
    "physics",
    "math",
    "chemistry",
    "interfaces",
    "elasticity",
    "seismic",
    "ha",
    "qha",
    "eos",
    "thermoelasticity",
    "examples",
)
_DOMAIN_LABELS = {
    "infrastructure": "Infrastructure",
    "physics": "Shared physics",
    "math": "Mathematics",
    "chemistry": "Chemistry",
    "interfaces": "QM interfaces",
    "elasticity": "Elasticity",
    "seismic": "SEISMIC",
    "ha": "HA",
    "qha": "QHA",
    "eos": "EOS",
    "thermoelasticity": "Thermoelasticity",
    "examples": "Real examples",
}
_LAYER_ORDER = {
    "library": 0,
    "cli": 1,
    "plotting": 2,
}

_verbose_sections = False
_active_suite = "all"


def pytest_addoption(parser: pytest.Parser) -> None:
    """Register Quantas-specific suite selection options."""
    group = parser.getgroup("quantas test suites")
    group.addoption(
        "--suite",
        action="store",
        default="all",
        choices=tuple(SUITE_BY_NAME),
        help="Run one named Quantas test suite (default: all).",
    )
    group.addoption(
        "--list-suites",
        action="store_true",
        default=False,
        help="List available Quantas test suites and exit.",
    )


def pytest_configure(config: pytest.Config) -> None:
    """Initialize suite selection and optional verbose reporting."""
    global _verbose_sections, _active_suite
    _verbose_sections = bool(config.getoption("verbose") > 0)
    _active_suite = str(config.getoption("--suite"))

    if config.getoption("--list-suites"):
        width = max(len(name) for name, _ in list_suite_rows())
        print("Quantas test suites:")
        for name, description in list_suite_rows():
            print(f"  {name:<{width}}  {description}")
        pytest.exit("Suite list requested", returncode=0)


def _item_markers(item: pytest.Item) -> frozenset[str]:
    explicit = (marker.name for marker in item.iter_markers())
    return classify_nodeid(item.nodeid, explicit)


def _layer_priority(markers: frozenset[str]) -> int:
    if "plotting" in markers:
        return _LAYER_ORDER["plotting"]
    if "cli" in markers:
        return _LAYER_ORDER["cli"]
    return _LAYER_ORDER["library"]


def _validate_staged_partition(
    items: list[pytest.Item],
    classified: dict[str, frozenset[str]],
) -> None:
    """Ensure the complete logical partition covers every test exactly once."""
    invalid: list[str] = []
    for item in items:
        markers = classified[item.nodeid]
        memberships = [
            name for name in STAGED_ALL_SUITES if SUITE_BY_NAME[name].predicate(markers)
        ]
        if len(memberships) != 1:
            invalid.append(f"{item.nodeid}: {memberships or ['unclassified']}")
    if invalid:
        details = "\n".join(f"  {entry}" for entry in invalid[:20])
        raise pytest.UsageError(
            "Quantas staged suites must partition all tests exactly once:\n" + details
        )


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(
    config: pytest.Config,
    items: list[pytest.Item],
) -> None:
    """Assign orthogonal markers, select a suite and order collected tests."""
    domain_order = {name: index for index, name in enumerate(_DOMAIN_ORDER)}
    classified: dict[str, frozenset[str]] = {}

    for item in items:
        markers = _item_markers(item)
        classified[item.nodeid] = markers
        for marker in markers:
            item.add_marker(getattr(pytest.mark, marker))

    _validate_staged_partition(items, classified)

    suite = SUITE_BY_NAME[_active_suite]
    selected: list[pytest.Item] = []
    deselected: list[pytest.Item] = []
    for item in items:
        if suite.predicate(classified[item.nodeid]):
            selected.append(item)
        else:
            deselected.append(item)

    selected.sort(
        key=lambda item: (
            _layer_priority(classified[item.nodeid]),
            domain_order.get(domain_for_nodeid(item.nodeid) or "", len(domain_order)),
            item.nodeid,
        )
    )
    items[:] = selected
    if deselected:
        config.hook.pytest_deselected(items=deselected)


def pytest_collection_finish(session: pytest.Session) -> None:
    """Print selected test counts by domain and execution concern."""
    if not _verbose_sections:
        return
    reporter = session.config.pluginmanager.getplugin("terminalreporter")
    if reporter is None:
        return

    domain_counts: Counter[str] = Counter()
    concern_counts: Counter[str] = Counter()
    for item in session.items:
        markers = _item_markers(item)
        domain = domain_for_nodeid(item.nodeid)
        if domain is not None:
            domain_counts[domain] += 1
        for concern in ("library", "cli", "plotting", "hdf5", "export", "io"):
            if concern in markers:
                concern_counts[concern] += 1

    reporter.write_sep("-", f"Quantas suite: {_active_suite}")
    for domain in _DOMAIN_ORDER:
        if domain_counts[domain]:
            reporter.write_line(
                f"  {_DOMAIN_LABELS[domain]:<18} {domain_counts[domain]:>4} tests"
            )
    reporter.write_line(
        "  Concerns: "
        + ", ".join(
            f"{name}={concern_counts[name]}"
            for name in ("library", "cli", "plotting", "hdf5", "export", "io")
            if concern_counts[name]
        )
    )


def pytest_report_header(config: pytest.Config) -> str:
    """Describe the selected suite and deterministic execution order."""
    suite = config.getoption("--suite")
    return (
        f"Quantas suite: {suite}; order: numbered domains within "
        "library -> CLI -> plotting"
    )
