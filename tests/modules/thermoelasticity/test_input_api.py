"""Tests for thermoelastic input generation and QHA context preparation."""

from __future__ import annotations

from pathlib import Path

import numpy as np

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


def test_generator_requires_crystal_pressure_keyword(tmp_path: Path) -> None:
    """A corrected tensor without explicit PRESSURE provenance is rejected."""
    path = tmp_path / "uncorrected.out"
    _write_crystal_soec(
        path,
        pressure=0.0,
        volume=100.0,
        density=3.0,
        energy=-100.0,
    )
    text = path.read_text(encoding="utf-8")
    text = text.replace("PRESSURE\n0.0\n", "")
    path.write_text(text, encoding="utf-8")

    with np.testing.assert_raises_regex(ValueError, "PRESSURE keyword is required"):
        create_thermoelastic_input(path, tmp_path / "invalid.yaml")
