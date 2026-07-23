"""Tests for the logical Quantas pytest-suite registry."""

from __future__ import annotations

import pytest

from tests.suite_registry import SUITE_BY_NAME, classify_nodeid


@pytest.mark.parametrize(
    ("nodeid", "expected"),
    [
        (
            "tests/modules/qha/test_analysis.py::test_analysis",
            {"qha", "module", "library", "scientific", "fast"},
        ),
        (
            "tests/modules/qha/test_plotting.py::test_renderer",
            {"qha", "module", "plotting", "matplotlib", "slow"},
        ),
        (
            "tests/modules/ha/test_cli.py::test_ha_run_help",
            {"ha", "module", "cli", "scientific", "fast"},
        ),
        (
            "tests/physics/elasticity/test_batch_sampling.py::test_batch",
            {"elasticity", "library", "scientific", "fast"},
        ),
        (
            "tests/physics/earth_profiles/test_pressure.py::test_prem",
            {"physics", "library", "scientific", "fast"},
        ),
    ],
)
def test_classifier_assigns_orthogonal_domain_and_concern_markers(
    nodeid: str,
    expected: set[str],
) -> None:
    markers = classify_nodeid(nodeid)
    assert expected.issubset(markers)


def test_qha_suite_excludes_plotting_but_qha_all_includes_it() -> None:
    plotting_markers = classify_nodeid(
        "tests/modules/qha/test_plotting.py::test_renderer"
    )
    library_markers = classify_nodeid(
        "tests/modules/qha/test_analysis.py::test_analysis"
    )

    assert SUITE_BY_NAME["qha"].predicate(library_markers)
    assert not SUITE_BY_NAME["qha"].predicate(plotting_markers)
    assert SUITE_BY_NAME["qha-all"].predicate(plotting_markers)
    assert SUITE_BY_NAME["qha-plotting"].predicate(plotting_markers)


def test_scientific_suite_excludes_cli_and_plotting_surfaces() -> None:
    library = classify_nodeid(
        "tests/modules/elasticity/test_analysis.py::test_analysis"
    )
    cli = classify_nodeid(
        "tests/modules/elasticity/test_cli.py::test_run_command"
    )
    plotting = classify_nodeid(
        "tests/modules/elasticity/test_plotting.py::test_rendered_plot"
    )

    predicate = SUITE_BY_NAME["scientific"].predicate
    assert predicate(library)
    assert not predicate(cli)
    assert not predicate(plotting)


def test_explicit_concern_marker_overrides_filename_inference() -> None:
    markers = classify_nodeid(
        "tests/modules/qha/test_plotting.py::test_data_contract",
        explicit_markers=("library",),
    )

    assert "library" in markers
    assert "plotting" not in markers


def test_eos_plotting_suite_is_part_of_the_plotting_partition() -> None:
    markers = classify_nodeid(
        "tests/modules/eos/test_plotting.py::test_rendered_plot"
    )

    assert SUITE_BY_NAME["eos-all"].predicate(markers)
    assert SUITE_BY_NAME["eos-plotting"].predicate(markers)
    assert not SUITE_BY_NAME["eos"].predicate(markers)


def test_thermoelasticity_plotting_suite_is_staged_separately() -> None:
    """Thermoelastic plots remain outside library/CLI partitions."""
    markers = classify_nodeid(
        "tests/modules/thermoelasticity/test_plotting.py::test_rendered_plot"
    )

    assert SUITE_BY_NAME["thermoelasticity-all"].predicate(markers)
    assert SUITE_BY_NAME["thermoelasticity-plotting"].predicate(markers)
    assert not SUITE_BY_NAME["thermoelasticity"].predicate(markers)
