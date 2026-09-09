"""Tests for thermoelastic input generation and QHA context preparation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from quantas.api.thermoelasticity import prepare_context
from quantas.models import ResultData, ResultMetadata
from quantas.modules.qha.models import QHAResult
from quantas.modules.thermoelasticity.api import (
    create_thermoelastic_input,
    read_thermoelastic_input,
)


def _write_crystal_soec(
    path: Path,
    *,
    pressure: float,
    volume: float,
    density: float,
    energy: float,
) -> None:
    """Write a compact cubic CRYSTAL-like SOEC output fixture."""
    a = volume ** (1.0 / 3.0)
    c11 = 200.0 + 4.0 * pressure
    c12 = 100.0 + 2.0 * pressure
    c44 = 80.0 + pressure
    path.write_text(
        f"""ELASTCON OPTION
PRESSURE
{pressure}
FINAL OPTIMIZED GEOMETRY - DIMENSIONALITY OF THE SYSTEM      3
LATTICE PARAMETERS (ANGSTROMS AND DEGREES) - BOHR = 0.5291772083 ANGSTROM
PRIMITIVE CELL - CENTRING CODE 1/0 VOLUME= {volume:.10f} - DENSITY {density:.6f} g/cm^3
        A              B              C           ALPHA      BETA       GAMMA
 {a:.12f} {a:.12f} {a:.12f} 90.000000 90.000000 90.000000
ATOMS IN THE ASYMMETRIC UNIT    1 - ATOMS IN THE UNIT CELL:    1
     ATOM                 X/A                 Y/B                 Z/C
      1 T  12 MG    0.000000000000E+00 0.000000000000E+00 0.000000000000E+00
DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)
          X                    Y                    Z
 {a:.12f} 0.000000000000 0.000000000000
 0.000000000000 {a:.12f} 0.000000000000
 0.000000000000 0.000000000000 {a:.12f}
TOTAL ENERGY(DFT)(AU)(  2) {energy:.12E} DE 1.0E-12 tester 1.0E-12
PRESSURE IN GIGAPASCAL: {pressure:.8E}
VOLUME OF THE CELL: {volume:.10f}
DENSITY OF THE CRYSTAL = {density:.8f} g/cm^3
FINAL RESULTS START
ELASTIC PROPERTIES AT PRESSURE (GPa) = {pressure:.8f}
SYMMETRIZED ELASTIC CONSTANTS FOR CUBIC CASE, IN GPa

| {c11:.6f} {c12:.6f} {c12:.6f} 0.000000 0.000000 0.000000 |
|          {c11:.6f} {c12:.6f} 0.000000 0.000000 0.000000 |
|                   {c11:.6f} 0.000000 0.000000 0.000000 |
|                            {c44:.6f} 0.000000 0.000000 |
|                                     {c44:.6f} 0.000000 |
|                                              {c44:.6f} |
""",
        encoding="utf-8",
    )


def test_generator_sorts_points_and_round_trips(tmp_path: Path) -> None:
    """Unordered CRYSTAL files produce a readable volume-sorted YAML file."""
    files = []
    for name, pressure, volume in (
        ("expanded.out", -2.0, 110.0),
        ("compressed.out", 2.0, 90.0),
        ("reference.out", 0.0, 100.0),
    ):
        path = tmp_path / name
        _write_crystal_soec(
            path,
            pressure=pressure,
            volume=volume,
            density=3.0,
            energy=-100.0 + pressure,
        )
        files.append(path)

    output = create_thermoelastic_input(
        [files[0], files[1], files[2]],
        tmp_path / "thermoelastic.yaml",
        jobname="Synthetic cubic",
    )
    text = output.read_text(encoding="utf-8")
    assert "fractional_positions:\n    - [" in text
    assert "stiffness:\n  - [" in text

    parsed = read_thermoelastic_input(output)
    assert parsed.jobname == "Synthetic cubic"
    assert parsed.elastic_series.npoints == 3
    assert np.allclose(parsed.elastic_series.volumes, [90.0, 100.0, 110.0])
    assert parsed.elastic_series.reference_index == 1
    assert parsed.elastic_series.elastic_symmetry == "cubic"
    assert parsed.elastic_series.symmetry.space_group_number == 221
    frame = parsed.elastic_series.metadata["frame_normalization"]
    assert frame["method"] == "right_polar_decomposition_corotation"
    assert frame["maximum_removed_rotation_degrees"] == 0.0



def test_generator_reuses_kieffer_qha_static_energy(
    tmp_path: Path,
) -> None:
    """QSA pressure fitting consumes E(V) independently of a Kieffer block."""
    volumes = [90.0, 95.0, 100.0, 105.0, 110.0]
    energies = [-100.0 + 1.0e-4 * (volume - 100.0) ** 2 for volume in volumes]
    qha_input = tmp_path / "qha-kieffer.yaml"
    qha_input.write_text(
        yaml.safe_dump(
            {
                "job": "QSA static E(V) reuse",
                "natom": 1,
                "formula_units": 1,
                "units": {
                    "energy": "Ha",
                    "volume": "angstrom^3",
                    "frequency": "cm^-1",
                    "length": "angstrom",
                },
                "supercell": np.eye(3, dtype=int).tolist(),
                "q_position_source": "crystal-output",
                "q_position_convention": "fractional-reciprocal",
                "qpoints": 1,
                "volume": volumes,
                "energy": energies,
                "phonon": [
                    {
                        "q-position": [0.0, 0.0, 0.0],
                        "weight": 1.0,
                        "band": [
                            {"frequency": [0.0] * len(volumes)} for _ in range(3)
                        ],
                    }
                ],
                "kieffer": {"model": "sine-wave", "provenance": {"test": True}},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    files: list[Path] = []
    for index, volume in enumerate(reversed(volumes)):
        path = tmp_path / f"elastic-energy-{index}.out"
        _write_crystal_soec(
            path,
            pressure=0.0,
            volume=volume,
            density=3.0,
            energy=-1.0,
        )
        text = path.read_text(encoding="utf-8")
        text = text.replace("PRESSURE\n0.0\n", "")
        text = text.replace("PRESSURE IN GIGAPASCAL: 0.00000000E+00\n", "")
        path.write_text(text, encoding="utf-8")
        files.append(path)

    output = create_thermoelastic_input(
        files,
        tmp_path / "energy-input.yaml",
        pressure_source="energy_polynomial",
        polynomial_degree=2,
        energy_input=qha_input,
    )
    parsed = read_thermoelastic_input(output)
    model = parsed.elastic_series.metadata["pressure_resolution"]["energy_model"]
    assert model["source_dataset"] == str(qha_input)
    assert len(model["volume_matches"]) == len(volumes)
    assert model["fit"]["success"] is True

def test_prepare_context_reports_elastic_extrapolation(tmp_path: Path) -> None:
    """The API marks QHA volumes outside the sampled elastic interval."""
    files = []
    for index, (pressure, volume) in enumerate(
        ((2.0, 90.0), (0.0, 100.0), (-2.0, 110.0))
    ):
        path = tmp_path / f"point_{index}.out"
        _write_crystal_soec(
            path,
            pressure=pressure,
            volume=volume,
            density=3.0,
            energy=-100.0 + pressure,
        )
        files.append(path)
    yaml_path = create_thermoelastic_input(files, tmp_path / "thermoelastic.yaml")

    qha = QHAResult(
        temperature=np.asarray([300.0, 400.0]),
        pressure=np.asarray([0.0, 1.0]),
        volume=np.asarray([90.0, 100.0, 110.0]),
        static_energy=np.asarray([-9.8, -10.0, -9.7]),
        equilibrium_volume=np.asarray([[95.0, 105.0], [111.0, 89.0]]),
        isochoric_heat_capacity=np.ones((2, 2)),
        isothermal_bulk_modulus=np.ones((2, 2)) * 100.0,
        bulk_modulus_derivative=np.ones((2, 2)) * 4.0,
        thermal_expansion=np.ones((2, 2)),
        axial_thermal_expansion=np.ones((2, 2, 3)),
        thermal_expansion_tensor=np.ones((2, 2, 3, 3)),
        equilibrium_lattice=np.ones((2, 2, 3, 3)),
        lattice_parameters=np.ones((2, 2, 6)),
    )
    result = ResultData(
        metadata=ResultMetadata(module="qha", method="quasi-harmonic"),
        results={"qha": qha},
    )
    context = prepare_context(yaml_path, result)
    assert context.extrapolation_mask.tolist() == [[False, False], [True, True]]
    assert context.has_complete_quasistatic_inputs
    assert context.has_complete_adiabatic_inputs
    assert context.metadata["extrapolated_points"] == 2


def test_list_file_paths_are_local_and_sorted_by_volume(tmp_path: Path) -> None:
    """Relative list entries are resolved locally and input order is ignored."""
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()
    files: list[Path] = []
    for name, pressure, volume in (
        ("expanded.out", -2.0, 110.0),
        ("compressed.out", 2.0, 90.0),
        ("reference.out", 0.0, 100.0),
    ):
        path = output_dir / name
        _write_crystal_soec(
            path,
            pressure=pressure,
            volume=volume,
            density=3.0,
            energy=-100.0 + pressure,
        )
        files.append(path)
    list_file = tmp_path / "soec.txt"
    list_file.write_text(
        "# deliberately unordered\n"
        "outputs/reference.out\n"
        "outputs/expanded.out\n"
        "outputs/compressed.out\n",
        encoding="utf-8",
    )

    yaml_path = create_thermoelastic_input(
        list_file,
        tmp_path / "from_list.yaml",
        is_list=True,
    )
    parsed = read_thermoelastic_input(yaml_path)
    assert np.allclose(parsed.elastic_series.volumes, [90.0, 100.0, 110.0])
    assert [point.source for point in parsed.elastic_series.points] == [
        files[1].name,
        files[2].name,
        files[0].name,
    ]


def test_generator_preserves_backend_pressure_keyword(
    tmp_path: Path,
) -> None:
    """A backend PRESSURE value remains authoritative if the elastic line is absent."""
    path = tmp_path / "corrected-no-elastic-pressure.out"
    _write_crystal_soec(
        path,
        pressure=2.0,
        volume=100.0,
        density=3.0,
        energy=-100.0,
    )
    text = path.read_text(encoding="utf-8")
    text = text.replace("ELASTIC PROPERTIES AT PRESSURE (GPa) = 2.00000000\n", "")
    path.write_text(text, encoding="utf-8")

    output = create_thermoelastic_input(path, tmp_path / "corrected.yaml")
    parsed = read_thermoelastic_input(output)
    resolution = parsed.elastic_series.metadata["pressure_resolution"]
    assert parsed.elastic_series.points[0].pressure == 2.0
    assert resolution["states"][0]["pressure_source"] == "applied_prestress"
    assert resolution["states"][0]["correction_applied_by"] == "crystal"


def test_generator_auto_corrects_raw_tensor_from_output_stress(tmp_path: Path) -> None:
    """Raw CRYSTAL tensors use output stress and record one Wallace correction."""
    path = tmp_path / "raw.out"
    _write_crystal_soec(
        path,
        pressure=2.0,
        volume=100.0,
        density=3.0,
        energy=-100.0,
    )
    text = path.read_text(encoding="utf-8")
    text = text.replace("PRESSURE\n2.0\n", "")
    path.write_text(text, encoding="utf-8")

    output = create_thermoelastic_input(path, tmp_path / "auto.yaml")
    parsed = read_thermoelastic_input(output)
    resolution = parsed.elastic_series.metadata["pressure_resolution"]
    assert resolution["requested_source"] == "auto"
    assert resolution["states"][0]["pressure_source"] == "output_stress"
    assert resolution["states"][0]["correction_method"] == (
        "crystal-erba-2014-hydrostatic"
    )
    assert resolution["states"][0]["correction_applied_by"] == (
        "quantas-thermoelastic-inpgen"
    )
    point = parsed.elastic_series.points[0]
    assert point.pressure == 2.0
    assert point.stiffness[0, 0] == pytest.approx(208.0)
    assert point.stiffness[0, 1] == pytest.approx(106.0)
    assert point.stiffness[3, 3] == pytest.approx(81.0)


def test_generator_rejects_raw_tensor_without_pressure(
    tmp_path: Path,
) -> None:
    """A raw tensor without pressure information cannot enter QSA silently."""
    path = tmp_path / "raw-no-pressure.out"
    _write_crystal_soec(
        path,
        pressure=0.0,
        volume=100.0,
        density=3.0,
        energy=-100.0,
    )
    text = path.read_text(encoding="utf-8")
    text = text.replace("PRESSURE\n0.0\n", "")
    text = text.replace("PRESSURE IN GIGAPASCAL: 0.00000000E+00\n", "")
    path.write_text(text, encoding="utf-8")

    with np.testing.assert_raises_regex(ValueError, "Select an explicit pressure source"):
        create_thermoelastic_input(path, tmp_path / "invalid.yaml")


def test_generator_energy_polynomial_corrects_raw_series(tmp_path: Path) -> None:
    """Static elastic-output E(V) can supply pressure for raw QSA tensors."""
    files: list[Path] = []
    volumes = np.asarray([90.0, 95.0, 100.0, 105.0, 110.0])
    for index, volume in enumerate(volumes):
        path = tmp_path / f"raw-{index}.out"
        energy = -100.0 + 1.0e-4 * (volume - 100.0) ** 2
        _write_crystal_soec(
            path,
            pressure=0.0,
            volume=float(volume),
            density=3.0,
            energy=float(energy),
        )
        text = path.read_text(encoding="utf-8")
        text = text.replace("PRESSURE\n0.0\n", "")
        text = text.replace("PRESSURE IN GIGAPASCAL: 0.00000000E+00\n", "")
        path.write_text(text, encoding="utf-8")
        files.append(path)

    output = create_thermoelastic_input(
        files,
        tmp_path / "energy.yaml",
        pressure_source="energy_polynomial",
        polynomial_degree=2,
    )
    parsed = read_thermoelastic_input(output)
    resolution = parsed.elastic_series.metadata["pressure_resolution"]
    assert resolution["requested_source"] == "energy_polynomial"
    assert resolution["energy_model"]["relation"] == "P(V) = -dE/dV"
    assert resolution["energy_model"]["fit"]["success"] is True
    assert resolution["correction_formulation"] == "crystal-erba-2014-hydrostatic"
    assert all(
        state["pressure_source"] == "energy_polynomial"
        for state in resolution["states"]
    )
    for point in parsed.elastic_series.points:
        pressure = point.pressure
        assert point.stiffness[0, 0] == pytest.approx(200.0)
        assert point.stiffness[0, 1] == pytest.approx(100.0 + pressure)
        assert point.stiffness[3, 3] == pytest.approx(80.0 - 0.5 * pressure)
