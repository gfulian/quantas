"""Integration and scientific-regression checks for curated real examples."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest

from quantas.interfaces.crystal.phonons import CrystalPhononReader
from quantas.api import elasticity, eos, ha, qha, thermoelasticity
from quantas.io.hdf5 import read_input_data, write_input_data
from quantas.modules.qha.calculator import QHACalculator


pytestmark = [pytest.mark.examples, pytest.mark.scientific_regression]

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = PROJECT_ROOT / "examples"
_IGNORED_EXAMPLE_DIRECTORIES = {
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
}
_IGNORED_EXAMPLE_SUFFIXES = {".pyc", ".pyo"}


def _is_curated_example_file(path: Path) -> bool:
    """Return whether *path* belongs to the distributed example corpus."""
    if not path.is_file():
        return False
    relative = path.relative_to(EXAMPLES)
    if path.name in {"manifest.json", "MANIFEST.sha256"}:
        return False
    if any(part in _IGNORED_EXAMPLE_DIRECTORIES for part in relative.parts):
        return False
    return path.suffix not in _IGNORED_EXAMPLE_SUFFIXES


def _assert_examples_manifest_current() -> None:
    """Assert that distributed examples match the checked-in manifest."""
    completed = subprocess.run(
        [sys.executable, "tools/update_examples_manifest.py", "--check"],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr

    manifest_path = EXAMPLES / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    recorded = {entry["path"] for entry in manifest["files"]}
    actual = {
        path.relative_to(EXAMPLES).as_posix()
        for path in EXAMPLES.rglob("*")
        if _is_curated_example_file(path)
    }
    assert recorded == actual
    digest = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    assert (EXAMPLES / "MANIFEST.sha256").read_text(encoding="utf-8") == (
        f"{digest}  manifest.json\n"
    )


def test_examples_manifest_is_current_and_complete() -> None:
    """The checked-in manifest must cover every distributed example exactly."""
    _assert_examples_manifest_current()


def test_python_examples_use_only_the_supported_quantas_api() -> None:
    """Distributed scripts must not promote implementation paths as public API."""
    forbidden = ("quantas.core", "quantas.modules", "quantas.renderers")
    violations: list[str] = []
    for path in sorted(EXAMPLES.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if module.startswith(forbidden):
                    violations.append(f"{path.relative_to(EXAMPLES)}: {module}")
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith(forbidden):
                        violations.append(
                            f"{path.relative_to(EXAMPLES)}: {alias.name}"
                        )
    assert violations == []


@pytest.mark.parametrize(
    ("filename", "jobname"),
    [
        ("calcite.dat", "Calcite (B3LYP)"),
        ("hydroxylapatite.dat", "Hydroxylapatite"),
    ],
)
def test_real_elasticity_inputs_run_and_persist(
    filename: str,
    jobname: str,
    tmp_path: Path,
) -> None:
    """Real elasticity inputs must survive reader, calculation, and HDF5 output."""
    input_data = elasticity.read_input(EXAMPLES / "elasticity" / filename)
    assert input_data.jobname == jobname
    assert input_data.stiffness.shape == (6, 6)
    assert input_data.stiffness.dtype == np.dtype("float64")
    result = elasticity.run(input_data)
    output = elasticity.write_result(result, tmp_path / filename)
    with h5py.File(output, "r") as h5:
        assert h5["results/stiffness"].dtype == np.dtype("float64")
        assert np.all(np.isfinite(h5["results/stiffness"][...]))


@pytest.mark.parametrize(
    ("filename", "required_columns"),
    [
        ("PV_quartz.dat", {"pressure", "volume"}),
        ("PV_topaz.dat", {"pressure", "volume", "a", "b", "c"}),
        ("T_triclinic.dat", {"temperature", "volume"}),
        ("VT_rutile.dat", {"temperature", "volume", "a", "c"}),
        ("PVT_NaF.dat", {"temperature", "pressure", "volume"}),
    ],
)
def test_real_eos_inputs_are_normalized(
    filename: str,
    required_columns: set[str],
) -> None:
    """Every public EOS dataset must parse into a finite float64 table."""
    dataset = eos.read_input(EXAMPLES / "eos" / filename)
    assert dataset.npoints > 10
    assert required_columns.issubset(dataset.columns)
    for values in dataset.columns.values():
        assert values.dtype == np.dtype("float64")
        assert np.all(np.isfinite(values))


def test_real_molar_volume_convention_is_preserved() -> None:
    """Absolute molar volume declared through VSCALE remains molar volume."""
    dataset = eos.read_input(EXAMPLES / "eos" / "PVT_NaF.dat")
    assert dataset.units["volume"] == "cm^3/mol"
    assert dataset.raw_units["volume"] == "cm^3/mol"
    assert dataset.column("volume")[0] == pytest.approx(14.75734)


def test_real_quartz_dataset_recovers_reference_bm3_fit() -> None:
    """A complete real P-V path must produce the characterized quartz fit."""
    dataset = eos.read_input(EXAMPLES / "eos" / "PV_quartz.dat")
    result = eos.fit(dataset, eos.FitRequest(model="BM3"))
    assert result.fit.success
    assert result.parameter_values["V0"] == pytest.approx(112.96752, rel=2.0e-5)
    assert result.parameter_values["K0"] == pytest.approx(37.28543, rel=2.0e-4)
    assert result.parameter_values["KP"] == pytest.approx(5.93351, rel=2.0e-4)


@pytest.mark.parametrize(
    ("relative", "nvolumes", "nqpoints", "nmodes"),
    [
        ("crystal-qha/mgo_b3lyp.yaml", 11, 32, 6),
        ("crystal-phonons/dol_pbe0.yaml", 7, 27, 30),
    ],
)
def test_real_qha_yaml_contracts(
    relative: str,
    nvolumes: int,
    nqpoints: int,
    nmodes: int,
) -> None:
    """Normalized QHA YAML must preserve real volume, q-point, and mode axes."""
    data = qha.read_input(EXAMPLES / "qha" / relative)
    assert data.volume.shape == (nvolumes,)
    assert data.frequencies.shape == (nqpoints, nmodes, nvolumes)
    assert data.volume.dtype == np.dtype("float64")
    assert data.frequencies.dtype == np.dtype("float64")


@pytest.fixture(scope="module")
def generated_mgo_qha_yaml(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Generate one normalized MgO QHA input from the complete CRYSTAL output."""
    destination = tmp_path_factory.mktemp("mgo-qha-example") / "mgo_generated.yaml"
    ha.create_input(
        EXAMPLES / "qha/crystal-qha/mgo-b3lyp-crystal-qha.out",
        destination,
        interface="crystal-qha",
        jobname="MgO Periclase (CRYSTAL17 QHA)",
    )
    return destination


def test_native_crystal_qha_generation_matches_curated_yaml(
    generated_mgo_qha_yaml: Path,
) -> None:
    """Native CRYSTAL ingestion must reproduce the curated numerical contract."""
    generated = qha.read_input(generated_mgo_qha_yaml)
    expected = qha.read_input(EXAMPLES / "qha/crystal-qha/mgo_b3lyp.yaml")
    assert generated.jobname == expected.jobname
    assert generated.nvol == expected.nvol == 11
    assert generated.qpoints == expected.qpoints == 32
    assert generated.nmodes == expected.nmodes == 6
    np.testing.assert_allclose(generated.volume, expected.volume, atol=2.0e-8)
    np.testing.assert_allclose(generated.energy, expected.energy, atol=1.0e-10)
    np.testing.assert_allclose(generated.weights, expected.weights, atol=0.0)
    np.testing.assert_allclose(generated.frequencies, expected.frequencies, atol=1.0e-8)
    assert generated.qcoords is None
    assert generated.structure is not None
    assert generated.structure.natoms == 2
    assert generated.structure.nvol == 11
    assert generated.structure.normalization.repetitions == 32
    assert generated.structure.symmetry is not None
    assert generated.structure.symmetry.space_group_number == 225


def test_generated_mgo_structure_survives_hdf5_input_envelope(
    generated_mgo_qha_yaml: Path,
    tmp_path: Path,
) -> None:
    """The real reconstructed structure must survive shared HDF5 persistence."""
    qha_input = qha.read_input(generated_mgo_qha_yaml)
    calculator = QHACalculator(qha_input)
    destination = tmp_path / "mgo_input_structure.hdf5"
    with h5py.File(destination, "w") as handle:
        write_input_data(handle, calculator.input_data)
    with h5py.File(destination, "r") as handle:
        restored = read_input_data(handle)
    assert restored is not None
    assert restored.data["has_structure"] is True
    structure = restored.data["structure"]
    assert structure["representation"] == "primitive"
    assert structure["normalization"]["repetitions"] == 32
    np.testing.assert_allclose(
        structure["volume_series"]["volume"],
        qha_input.volume,
        atol=2.0e-8,
    )


def test_native_crystal_phonon_series_is_consistent() -> None:
    """All seven dolomite CRYSTAL phonon calculations must parse consistently."""
    directory = EXAMPLES / "qha/crystal-phonons"
    names = [
        line.strip()
        for line in (directory / "files.txt").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert len(names) == 7
    energies: list[float] = []
    for name in names:
        reader = CrystalPhononReader(directory / name)
        assert reader.completed, f"{name}: {reader.error}"
        assert reader.qpoints == 27
        assert reader.nphonon == 30
        energies.append(float(reader.energy))
    assert np.all(np.isfinite(energies))
    assert len(set(energies)) == 7


def test_real_thermoelastic_input_preserves_elastic_volume_series() -> None:
    """The dolomite QSA input must preserve all real elastic-volume points."""
    data = thermoelasticity.read_input(
        EXAMPLES / "thermoelasticity/dol_pbe0_thermoelastic.yaml"
    )
    assert data.jobname == "Dolomite PBE0 QSA thermoelasticity"
    assert data.method == "quasistatic"
    assert len(data.elastic_series.points) == 11
    volumes = np.asarray([point.volume for point in data.elastic_series.points])
    assert np.all(np.isfinite(volumes))
    assert np.all(np.diff(np.sort(volumes)) > 0.0)

    list_file = EXAMPLES / "thermoelasticity/crystal_outputs/dol_pbe0_soec.txt"
    listed = {
        line.strip()
        for line in list_file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    archived_sources = {point.source for point in data.elastic_series.points}
    assert listed == archived_sources


def test_examples_remain_immutable_after_regression_suite() -> None:
    """End-to-end checks must not modify their distributed source examples."""
    _assert_examples_manifest_current()


@pytest.mark.parametrize(
    ("data_name", "spec_name", "job_count", "accepted_request"),
    [
        ("PV_quartz.dat", "quartz_pv_tutorial.spec", 6, "bm3-effective-variance"),
        ("VT_rutile.dat", "rutile_vt_tutorial.spec", 7, "berman-effective-variance"),
        ("PVT_NaF.dat", "naf_pvt_tutorial.spec", 4, "mgd-full"),
    ],
)
def test_eos_tutorial_batch_specifications_resolve(
    data_name: str,
    spec_name: str,
    job_count: int,
    accepted_request: str,
) -> None:
    """Distributed EOS tutorial plans resolve into valid public requests."""
    dataset = eos.read_input(EXAMPLES / "eos" / data_name)
    document = eos.read_spec(EXAMPLES / "eos/specs" / spec_name)
    resolved = eos.resolve_spec(document, dataset)
    assert len(resolved.plan.jobs) == job_count
    accepted = [
        job.request.request_id for job in resolved.plan.jobs if job.accept
    ]
    assert accepted == [accepted_request]
    for job in resolved.plan.jobs:
        eos.validate_request(dataset, job.request)


def test_naf_tutorial_mgd_uses_molar_formula_unit_normalization() -> None:
    """The MGD tutorial must preserve the molar-volume convention of NaF."""
    dataset = eos.read_input(EXAMPLES / "eos/PVT_NaF.dat")
    document = eos.read_spec(EXAMPLES / "eos/specs/naf_pvt_tutorial.spec")
    resolved = eos.resolve_spec(document, dataset)
    job = next(
        item
        for item in resolved.plan.jobs
        if item.request.request_id == "mgd-full"
    )
    normalization = job.request.model.mgd_normalization
    assert normalization is not None
    assert normalization.volume_basis.value == "molar-formula-unit"
    assert normalization.formula == "NaF"
    assert normalization.formula_units_per_cell is None
    assert normalization.atoms_per_unit == pytest.approx(2.0)
