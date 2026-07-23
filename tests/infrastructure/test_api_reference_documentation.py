"""Coverage contracts for the curated public API reference."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re

from quantas import api


DOCS_ROOT = Path("docs/source/api")
_DIRECTIVE = re.compile(
    r"^\.\. auto(?:class|function|data):: "
    r"quantas\.api\.([a-z_]+)\.([A-Za-z0-9_]+)\s*$",
    flags=re.MULTILINE,
)


def _documented_symbols() -> dict[str, dict[str, list[Path]]]:
    """Return public names and the RST pages documenting them."""
    found: dict[str, dict[str, list[Path]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for path in DOCS_ROOT.rglob("*.rst"):
        text = path.read_text(encoding="utf-8")
        for namespace, name in _DIRECTIVE.findall(text):
            found[namespace][name].append(path)
    return found


def test_public_symbols_are_documented_once() -> None:
    """The reference follows each namespace's explicit ``__all__`` contract."""
    documented = _documented_symbols()
    assert set(documented) == set(api.__all__)

    for namespace_name in api.__all__:
        namespace = getattr(api, namespace_name)
        expected = set(namespace.__all__)
        actual = set(documented[namespace_name])
        assert actual == expected
        duplicates = {
            name: paths
            for name, paths in documented[namespace_name].items()
            if len(paths) != 1
        }
        assert duplicates == {}


def test_reference_excludes_internal_members() -> None:
    """Autodoc lists only explicit public aliases, never whole internal modules."""
    forbidden = (
        ":imported-members:",
        ".. automodule::",
        "quantas.core",
        "quantas.modules",
        "quantas.models",
        "quantas.renderers",
        "quantas.io",
    )
    for path in DOCS_ROOT.rglob("*.rst"):
        text = path.read_text(encoding="utf-8")
        for token in forbidden:
            assert token not in text, f"{path}: exposes {token}"


def test_large_namespaces_use_subpages() -> None:
    """EOS and Thermoelasticity remain navigable despite broad public surfaces."""
    eos = (DOCS_ROOT / "eos.rst").read_text(encoding="utf-8")
    thermo = (DOCS_ROOT / "thermoelasticity.rst").read_text(encoding="utf-8")
    for target in (
        "eos/contracts",
        "eos/fitting",
        "eos/batch_archive",
        "eos/postfit",
    ):
        assert target in eos
    for target in (
        "thermoelasticity/contracts",
        "thermoelasticity/calibration",
        "thermoelasticity/analysis",
        "thermoelasticity/plotting",
    ):
        assert target in thermo


def test_domain_pages_link_related_guides() -> None:
    """A reader can move between contract, rationale, example, and file schema."""
    expectations = {
        "ha.rst": ("../workflows/ha", "../tutorials/ha", "../formats/phonon_yaml"),
        "qha.rst": ("../workflows/qha", "../tutorials/qha", "../formats/ha_qha_hdf5"),
        "elasticity.rst": (
            "../workflows/elasticity",
            "../tutorials/elasticity",
            "../formats/elasticity_input",
        ),
        "seismic.rst": (
            "../workflows/seismic",
            "../tutorials/seismic",
            "../formats/elasticity_seismic_hdf5",
        ),
        "eos.rst": ("../workflows/eos", "../tutorials/eos", "../formats/eos_hdf5"),
        "thermoelasticity.rst": (
            "../workflows/thermoelasticity",
            "../tutorials/thermoelasticity",
            "../formats/thermoelastic_hdf5",
        ),
    }
    for filename, targets in expectations.items():
        text = (DOCS_ROOT / filename).read_text(encoding="utf-8")
        for target in targets:
            assert target in text
