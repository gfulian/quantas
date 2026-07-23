"""Logical test-suite registry used by the Quantas pytest configuration.

The registry keeps scientific domains and test concerns orthogonal.  A test can
therefore belong to a scientific domain such as ``qha`` and, independently, to
an execution concern such as ``plotting`` or ``cli``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import Callable, Iterable


MarkerSet = frozenset[str]
SuitePredicate = Callable[[MarkerSet], bool]


@dataclass(frozen=True)
class SuiteDefinition:
    """Describe one selectable pytest suite.

    Parameters
    ----------
    name
        Public suite name accepted by ``pytest --suite``.
    description
        Human-readable explanation shown by ``--list-suites``.
    predicate
        Function deciding whether a classified test belongs to the suite.
    """

    name: str
    description: str
    predicate: SuitePredicate


_DOMAIN_PATHS: tuple[tuple[str, str], ...] = (
    ("tests/infrastructure/", "infrastructure"),
    ("tests/physics/units/", "physics"),
    ("tests/physics/eos/", "physics"),
    ("tests/physics/thermodynamics/", "physics"),
    ("tests/physics/elasticity/", "elasticity"),
    ("tests/physics/seismic/", "seismic"),
    ("tests/physics/earth_profiles/", "physics"),
    ("tests/math/", "math"),
    ("tests/chemistry/", "chemistry"),
    ("tests/interfaces/", "interfaces"),
    ("tests/modules/elasticity/", "elasticity"),
    ("tests/modules/seismic/", "seismic"),
    ("tests/modules/ha/", "ha"),
    ("tests/modules/qha/", "qha"),
    ("tests/modules/eos/", "eos"),
    ("tests/modules/thermoelasticity/", "thermoelasticity"),
    ("tests/examples/", "examples"),
)

_DOMAIN_MARKERS = frozenset(
    {
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
    }
)

_CORE_DOMAINS = frozenset(
    {"infrastructure", "physics", "math", "chemistry", "interfaces"}
)
_MODULE_DOMAINS = frozenset(
    {"elasticity", "seismic", "ha", "qha", "eos", "thermoelasticity"}
)

_PURE_PLOTTING_FILES = frozenset({"test_plotting.py"})

_PURE_CLI_FILES = frozenset(
    {"test_cli.py", "test_cli_output.py", "test_cli_run.py"}
)


def _normalized_nodeid(nodeid: str) -> str:
    """Return a platform-independent, lower-case pytest node identifier."""
    return nodeid.replace("\\", "/").lower()


def _path_part(nodeid: str) -> str:
    """Return only the filesystem portion of a pytest node identifier."""
    return _normalized_nodeid(nodeid).split("::", maxsplit=1)[0]


def _file_name(nodeid: str) -> str:
    """Return the test module file name from a pytest node identifier."""
    return PurePosixPath(_path_part(nodeid)).name


def _test_name(nodeid: str) -> str:
    """Return the unparameterized test function name from a node identifier."""
    parts = _normalized_nodeid(nodeid).split("::")
    for part in parts[1:]:
        candidate = part.split("[", maxsplit=1)[0]
        if candidate.startswith("test_"):
            return candidate
    return ""


def domain_for_nodeid(nodeid: str) -> str | None:
    """Return the scientific domain associated with a pytest node identifier.

    Parameters
    ----------
    nodeid
        Pytest node identifier.

    Returns
    -------
    str or None
        Domain marker, or ``None`` when the path is not classified.
    """
    normalized = _path_part(nodeid)
    for path_fragment, marker in _DOMAIN_PATHS:
        if path_fragment in normalized or normalized.startswith(
            path_fragment.removeprefix("tests/")
        ):
            return marker
    return None


def _is_cli_test(nodeid: str) -> bool:
    """Return whether a test exercises the Click/text command surface."""
    if "tests/modules/" not in _path_part(nodeid):
        return False
    filename = _file_name(nodeid)
    test_name = _test_name(nodeid)
    if filename in _PURE_CLI_FILES:
        return True
    if filename == "test_renderer_cli.py":
        return (
            "_command_" in test_name
            or "root_help" in test_name
            or test_name.startswith("test_cli_")
            or "plot_defaults" in test_name
            or "plot_3d" in test_name
        )
    return "_cli" in filename or "cli_" in filename


def _is_plotting_test(nodeid: str) -> bool:
    """Return whether a test builds or renders graphical output."""
    if "tests/modules/" not in _path_part(nodeid):
        return False
    filename = _file_name(nodeid)
    test_name = _test_name(nodeid)
    if filename in _PURE_PLOTTING_FILES:
        return True
    if filename == "test_renderer_cli.py" and test_name in {
        "test_seismic_run_command_writes_report_and_hdf5_without_figures",
        "test_root_help_exposes_seismic_run_plot_and_export",
    }:
        return False
    plotting_terms = (
        "plot",
        "render",
        "figure",
        "colormap",
        "cmap",
        "matplotlib",
        "surface_mesh",
        "surface_renderer",
        "three_dimensional",
        "unit_sphere",
        "physical_3d",
        "polarization_overlay",
        "polarization_variants",
        "full_sphere_map",
        "spherical_renderer",
    )
    return any(term in test_name for term in plotting_terms)


def _is_hdf5_test(nodeid: str) -> bool:
    """Return whether a test exercises native HDF5 persistence."""
    normalized = _path_part(nodeid)
    return "hdf5" in normalized or "round_trip" in _test_name(nodeid)


def _is_export_test(nodeid: str) -> bool:
    """Return whether a test exercises a public exporter."""
    normalized = _path_part(nodeid) + "::" + _test_name(nodeid)
    return any(term in normalized for term in ("export", "csv", "table_export"))


def _is_io_test(nodeid: str) -> bool:
    """Return whether a test exercises readers or input generation."""
    normalized = _path_part(nodeid) + "::" + _test_name(nodeid)
    return any(
        term in normalized
        for term in ("reader", "inpgen", "_io.py", "interfaces/test_elasticity")
    )


def _is_architecture_test(nodeid: str) -> bool:
    """Return whether a test protects package architecture or public contracts."""
    normalized = _path_part(nodeid) + "::" + _test_name(nodeid)
    return any(
        term in normalized
        for term in (
            "architecture",
            "package_layout",
            "public_api",
            "module_contracts",
            "result_schema",
        )
    )


def classify_nodeid(nodeid: str, explicit_markers: Iterable[str] = ()) -> MarkerSet:
    """Classify one pytest item along domain and execution-concern axes.

    Parameters
    ----------
    nodeid
        Pytest node identifier.
    explicit_markers
        Marker names already attached directly to the test.

    Returns
    -------
    frozenset of str
        Complete marker set used for suite selection and ordering.
    """
    markers = set(explicit_markers) - _DOMAIN_MARKERS - {"module"}
    domain = domain_for_nodeid(nodeid)
    if domain is not None:
        markers.add(domain)
        if domain in _MODULE_DOMAINS:
            markers.add("module")

    explicit_concern = bool({"library", "cli", "plotting"}.intersection(markers))
    if not explicit_concern:
        if _is_cli_test(nodeid):
            markers.add("cli")
        if _is_plotting_test(nodeid):
            markers.add("plotting")
    if _is_hdf5_test(nodeid):
        markers.add("hdf5")
    if _is_export_test(nodeid):
        markers.add("export")
    if _is_io_test(nodeid):
        markers.add("io")
    if _is_architecture_test(nodeid):
        markers.add("architecture")

    if "baseline" in markers or "baseline" in _path_part(nodeid):
        markers.add("baseline")

    if (
        "plotting" not in markers
        and "cli" not in markers
        and domain != "examples"
    ):
        markers.add("library")
    if domain in _DOMAIN_MARKERS and domain not in {"infrastructure", "interfaces"}:
        markers.add("scientific")
    if domain == "examples":
        markers.add("scientific_regression")

    if "plotting" not in markers:
        markers.add("fast")
    else:
        markers.add("slow")

    if "plotting" in markers:
        markers.add("matplotlib")

    return frozenset(markers)


def _has_domain(markers: MarkerSet, domain: str) -> bool:
    return domain in markers


def _domain_library(domain: str) -> SuitePredicate:
    return lambda markers: _has_domain(markers, domain) and "library" in markers


def _domain_all(domain: str) -> SuitePredicate:
    return lambda markers: _has_domain(markers, domain)


def _domain_concern(domain: str, concern: str) -> SuitePredicate:
    return lambda markers: _has_domain(markers, domain) and concern in markers


def _domain_cli(domain: str) -> SuitePredicate:
    return lambda markers: (
        _has_domain(markers, domain) and "cli" in markers and "plotting" not in markers
    )


SUITES: tuple[SuiteDefinition, ...] = (
    SuiteDefinition("all", "Every collected test.", lambda markers: True),
    SuiteDefinition(
        "core",
        "Infrastructure, shared physics, mathematics, chemistry and interfaces.",
        lambda markers: (
            bool(_CORE_DOMAINS.intersection(markers)) and "library" in markers
        ),
    ),
    SuiteDefinition(
        "library",
        "All Python-library tests, including infrastructure and interfaces.",
        lambda markers: "library" in markers,
    ),
    SuiteDefinition(
        "scientific",
        "Scientific calculations excluding infrastructure, CLI and plotting.",
        lambda markers: "scientific" in markers and "library" in markers,
    ),
    SuiteDefinition(
        "architecture",
        "Package layout, public API and architectural contract tests.",
        lambda markers: "architecture" in markers,
    ),
    SuiteDefinition(
        "interfaces",
        "External-code parsers and adapters.",
        _domain_all("interfaces"),
    ),
    SuiteDefinition(
        "physics",
        "Shared physical-library tests.",
        _domain_all("physics"),
    ),
    SuiteDefinition(
        "math",
        "Numerical mathematics and fitting tests.",
        _domain_all("math"),
    ),
    SuiteDefinition(
        "chemistry",
        "Chemical utility tests.",
        _domain_all("chemistry"),
    ),
    SuiteDefinition(
        "elasticity",
        "Elasticity library/workflow tests, excluding CLI and plots.",
        _domain_library("elasticity"),
    ),
    SuiteDefinition(
        "elasticity-all",
        "Every Elasticity test.",
        _domain_all("elasticity"),
    ),
    SuiteDefinition(
        "elasticity-cli",
        "Elasticity command-line tests excluding graphical plot commands.",
        _domain_cli("elasticity"),
    ),
    SuiteDefinition(
        "elasticity-plotting",
        "Elasticity plotting and rendering tests.",
        _domain_concern("elasticity", "plotting"),
    ),
    SuiteDefinition(
        "seismic",
        "SEISMIC library/workflow tests, excluding CLI and plots.",
        _domain_library("seismic"),
    ),
    SuiteDefinition(
        "seismic-all",
        "Every SEISMIC test.",
        _domain_all("seismic"),
    ),
    SuiteDefinition(
        "seismic-cli",
        "SEISMIC command-line tests excluding graphical plot commands.",
        _domain_cli("seismic"),
    ),
    SuiteDefinition(
        "seismic-plotting",
        "SEISMIC plotting and rendering tests.",
        _domain_concern("seismic", "plotting"),
    ),
    SuiteDefinition(
        "ha",
        "HA library/workflow tests, excluding CLI and plots.",
        _domain_library("ha"),
    ),
    SuiteDefinition("ha-all", "Every HA test.", _domain_all("ha")),
    SuiteDefinition(
        "ha-cli",
        "HA command-line tests excluding graphical plot commands.",
        _domain_cli("ha"),
    ),
    SuiteDefinition(
        "ha-plotting",
        "HA plotting and rendering tests.",
        _domain_concern("ha", "plotting"),
    ),
    SuiteDefinition(
        "qha",
        "QHA library/workflow tests, excluding CLI and plots.",
        _domain_library("qha"),
    ),
    SuiteDefinition("qha-all", "Every QHA test.", _domain_all("qha")),
    SuiteDefinition(
        "qha-cli",
        "QHA command-line tests excluding graphical plot commands.",
        _domain_cli("qha"),
    ),
    SuiteDefinition(
        "qha-plotting",
        "QHA plotting and rendering tests.",
        _domain_concern("qha", "plotting"),
    ),
    SuiteDefinition(
        "eos",
        "EOS library/workflow tests, excluding CLI and plots.",
        _domain_library("eos"),
    ),
    SuiteDefinition(
        "eos-all",
        "Every EOS test.",
        _domain_all("eos"),
    ),
    SuiteDefinition(
        "eos-plotting",
        "EOS plotting and rendering tests.",
        _domain_concern("eos", "plotting"),
    ),
    SuiteDefinition(
        "thermoelasticity",
        "Thermoelasticity library/workflow tests, excluding CLI and plots.",
        _domain_library("thermoelasticity"),
    ),
    SuiteDefinition(
        "thermoelasticity-all",
        "Every Thermoelasticity test.",
        _domain_all("thermoelasticity"),
    ),
    SuiteDefinition(
        "thermoelasticity-cli",
        "Thermoelasticity command-line tests excluding graphical plot commands.",
        _domain_cli("thermoelasticity"),
    ),
    SuiteDefinition(
        "thermoelasticity-plotting",
        "Thermoelasticity plotting and rendering tests.",
        _domain_concern("thermoelasticity", "plotting"),
    ),
    SuiteDefinition(
        "cli",
        "Command-line and text-rendering tests excluding graphical plot tests.",
        lambda markers: "cli" in markers and "plotting" not in markers,
    ),
    SuiteDefinition(
        "cli-all",
        "Every command-line test, including graphical plot commands.",
        lambda markers: "cli" in markers,
    ),
    SuiteDefinition(
        "plotting",
        "All plotting/specification/rendering tests.",
        lambda markers: "plotting" in markers,
    ),
    SuiteDefinition(
        "matplotlib",
        "Static Matplotlib plotting tests.",
        lambda markers: "matplotlib" in markers,
    ),
    SuiteDefinition(
        "hdf5",
        "Native HDF5 persistence and round-trip tests.",
        lambda markers: "hdf5" in markers,
    ),
    SuiteDefinition(
        "export",
        "Table/CSV/exporter tests.",
        lambda markers: "export" in markers,
    ),
    SuiteDefinition(
        "io",
        "Reader, parser and input-generator tests.",
        lambda markers: "io" in markers,
    ),
    SuiteDefinition(
        "baseline",
        "Scientific regression baselines.",
        lambda markers: "baseline" in markers,
    ),
    SuiteDefinition(
        "examples",
        "Curated real-data examples and end-to-end scientific regressions.",
        _domain_all("examples"),
    ),
    SuiteDefinition(
        "scientific-regression",
        "Real-data and frozen scientific-regression tests.",
        lambda markers: (
            "scientific_regression" in markers or "baseline" in markers
        ),
    ),
    SuiteDefinition(
        "fast",
        "All tests except plotting/rendering tests.",
        lambda markers: "fast" in markers,
    ),
)

SUITE_BY_NAME = {suite.name: suite for suite in SUITES}

STAGED_ALL_SUITES: tuple[str, ...] = (
    "library",
    "examples",
    "cli",
    "elasticity-plotting",
    "seismic-plotting",
    "ha-plotting",
    "qha-plotting",
    "eos-plotting",
    "thermoelasticity-plotting",
)


def get_suite(name: str) -> SuiteDefinition:
    """Return a suite definition by name.

    Parameters
    ----------
    name
        Registered suite name.

    Returns
    -------
    SuiteDefinition
        Matching suite definition.

    Raises
    ------
    KeyError
        If the suite name is unknown.
    """
    return SUITE_BY_NAME[name]


def list_suite_rows() -> tuple[tuple[str, str], ...]:
    """Return suite names and descriptions for terminal output."""
    return tuple((suite.name, suite.description) for suite in SUITES)
