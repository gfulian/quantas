#!/usr/bin/env python3
"""Build-independent smoke validation for Quantas wheel and sdist artifacts.

The validator creates a fresh virtual environment for each selected artifact,
installs the distribution with its declared runtime dependencies, and verifies
that the command-line entry point, organized public API, package metadata, and
PEP 561 marker are available from the installed package rather than the source
tree.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import os
import subprocess
import tarfile
import tempfile
from typing import Literal, Sequence
import venv
import zipfile


ArtifactKind = Literal["all", "wheel", "sdist"]
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _parser() -> argparse.ArgumentParser:
    """Return the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "dist_dir",
        nargs="?",
        type=Path,
        default=PROJECT_ROOT / "dist",
        help="Directory containing built distributions (default: ./dist).",
    )
    parser.add_argument(
        "--kind",
        choices=("all", "wheel", "sdist"),
        default="all",
        help="Artifact kind to validate (default: all).",
    )
    parser.add_argument(
        "--keep-environments",
        action="store_true",
        help="Keep temporary virtual environments for manual inspection.",
    )
    return parser


def _artifacts(dist_dir: Path, kind: ArtifactKind) -> tuple[Path, ...]:
    """Return the unique artifacts selected for validation."""
    candidates: list[Path] = []
    if kind in {"all", "wheel"}:
        candidates.extend(sorted(dist_dir.glob("quantas-*.whl")))
    if kind in {"all", "sdist"}:
        candidates.extend(sorted(dist_dir.glob("quantas-*.tar.gz")))
    selected = tuple(dict.fromkeys(path.resolve() for path in candidates))
    if not selected:
        raise FileNotFoundError(
            f"No Quantas {kind} distribution artifacts found in {dist_dir}."
        )
    return selected




def _archive_names(artifact: Path) -> tuple[str, ...]:
    """Return normalized member names from one wheel or source archive."""
    if artifact.suffix == ".whl":
        with zipfile.ZipFile(artifact) as archive:
            return tuple(archive.namelist())
    with tarfile.open(artifact, mode="r:gz") as archive:
        return tuple(member.name for member in archive.getmembers() if member.isfile())


def _validate_archive_contents(artifacts: Sequence[Path]) -> None:
    """Validate public-file boundaries before installing distributions."""
    required_sdist_suffixes = (
        "/CHANGELOG.md",
        "/CONTRIBUTING.md",
        "/SECURITY.md",
        "/ROADMAP.md",
        "/PROJECT_STATE.md",
        "/RELEASE.md",
        "/CITATION.cff",
        "/requirements/minimum.txt",
        "/requirements/typecheck.txt",
        "/tools/validate_release.sh",
        "/examples/manifest.json",
        "/examples/MANIFEST.sha256",
        "/examples/qha/crystal-qha/mgo-b3lyp-crystal-qha.out",
    )
    forbidden_fragments = (
        "project_internal",
        "validation/phase",
        "PROJECT_VALIDATION",
        "PACKAGE_MANIFEST_SHA256",
        "__pycache__",
        ".pyc",
        ".pyo",
        "/src/quantas/api.py",
        "/src/quantas/eos.py",
    )
    for artifact in artifacts:
        names = _archive_names(artifact)
        lowered = tuple(name.lower() for name in names)
        for fragment in forbidden_fragments:
            if any(fragment.lower() in name for name in lowered):
                raise RuntimeError(
                    f"{artifact.name} contains forbidden internal material: {fragment}"
                )
        if artifact.suffix == ".whl":
            if "quantas/py.typed" not in names:
                raise RuntimeError(f"{artifact.name} is missing quantas/py.typed")
            if any(name.startswith("examples/") for name in names):
                raise RuntimeError(
                    f"{artifact.name} must not embed repository examples"
                )
            continue
        missing = [
            suffix
            for suffix in required_sdist_suffixes
            if not any(name.endswith(suffix) for name in names)
        ]
        if missing:
            raise RuntimeError(
                f"{artifact.name} is missing source-distribution files: "
                + ", ".join(missing)
            )

def _venv_command(environment: Path) -> Path:
    """Return the Quantas console entry point inside a virtual environment."""
    if os.name == "nt":
        return environment / "Scripts" / "quantas.exe"
    return environment / "bin" / "quantas"


def _venv_python(environment: Path) -> Path:
    """Return the Python executable inside a virtual environment."""
    if os.name == "nt":
        return environment / "Scripts" / "python.exe"
    return environment / "bin" / "python"


def _run(command: Sequence[str], *, cwd: Path | None = None) -> None:
    """Run one subprocess and raise on failure."""
    print("$ " + " ".join(str(part) for part in command), flush=True)
    subprocess.run(list(command), cwd=cwd, check=True)


def _smoke_script() -> str:
    """Return the installed-package smoke-test program."""
    return r"""
from importlib import metadata, resources
from pathlib import Path

import quantas
from quantas.api import elasticity, eos, ha, plotting, qha, registry, seismic, thermoelasticity
from quantas.api.registry import Capability
from quantas.core.physics.elasticity import cold_finite_strain_component

version = metadata.version("quantas")
assert version == quantas.__version__
assert registry.get("eos").has(Capability.FIT)
assert registry.get("thermoelasticity").has(Capability.RUN_CONTEXT)
plot_namespaces = {
    "elasticity": elasticity,
    "seismic": seismic,
    "ha": ha,
    "qha": qha,
    "thermoelasticity": thermoelasticity,
    "eos": eos,
}
for name, namespace in plot_namespaces.items():
    descriptor = registry.get(name)
    assert descriptor.has(Capability.PLOT_INVENTORY)
    assert callable(descriptor.operation(Capability.PLOT_INVENTORY))
    assert callable(namespace.describe_plots)
assert plotting.PlotInventory.__module__ == "quantas.models.plot_inventory"
assert {item.name for item in registry.list_modules()} == {
    "elasticity", "seismic", "ha", "qha", "eos", "thermoelasticity"
}
capabilities = {
    item.domain.value: item.status.value for item in eos.DOMAIN_CAPABILITIES
}
assert capabilities["pv"] == "public"
assert capabilities["ev"] == "core_only"
assert callable(thermoelasticity.run_context)
assert callable(thermoelasticity.write_result)
assert thermoelasticity.Options(report_level="debug").solver_debug
reference_component = cold_finite_strain_component(
    100.0,
    reference_volume=100.0,
    bulk_modulus=100.0,
    bulk_modulus_derivative=4.0,
    reference_component=200.0,
    component_pressure_derivative=5.0,
    wallace_delta=-3.0,
)
assert float(reference_component) == 200.0

requirements = metadata.requires("quantas") or []
assert not any("dash" in item.lower() for item in requirements)

marker = resources.files("quantas").joinpath("py.typed")
assert marker.is_file(), "Installed wheel is missing quantas/py.typed"

package_path = Path(quantas.__file__).resolve()
print(f"Quantas {version} installed from {package_path}")
print("Organized public API:", ", ".join(item.name for item in registry.list_modules()))
print("Thermoelastic QSA API: available")
"""


def _validate_artifacts(artifacts: Sequence[Path], root: Path) -> None:
    """Install and smoke-test artifacts sequentially in one fresh environment."""
    environment = root / "environment"
    venv.EnvBuilder(with_pip=True, clear=True).create(environment)
    python = _venv_python(environment)
    _run(
        [
            str(python),
            "-m",
            "pip",
            "install",
            "--no-compile",
            "setuptools>=77",
            "wheel",
        ]
    )

    for index, artifact in enumerate(artifacts):
        if index:
            _run([str(python), "-m", "pip", "uninstall", "-y", "quantas"])
        command = [
            str(python),
            "-m",
            "pip",
            "install",
            "--no-compile",
            "--timeout",
            "120",
            "--retries",
            "8",
        ]
        if index:
            command.extend(["--no-deps", "--no-build-isolation"])
        command.append(str(artifact))
        _run(command)
        _run([str(_venv_command(environment)), "--version"], cwd=root)
        _run([str(python), "-c", _smoke_script()], cwd=root)


def main(argv: Sequence[str] | None = None) -> int:
    """Validate selected Quantas distribution artifacts."""
    args = _parser().parse_args(argv)
    dist_dir = args.dist_dir.resolve()
    artifacts = _artifacts(dist_dir, args.kind)
    _validate_archive_contents(artifacts)

    if args.keep_environments:
        root = Path(tempfile.mkdtemp(prefix="quantas-dist-check-"))
        print(f"Keeping validation environments in {root}")
        _validate_artifacts(artifacts, root)
    else:
        with tempfile.TemporaryDirectory(prefix="quantas-dist-check-") as temporary:
            root = Path(temporary)
            _validate_artifacts(artifacts, root)

    print(f"Validated {len(artifacts)} Quantas distribution artifact(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
