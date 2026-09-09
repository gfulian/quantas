"""Tests for CRYSTAL multi-volume elastic-state import."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.interfaces.crystal import (
    CrystalPressurePolicy,
    crystal_hydrostatic_stiffness,
    read_crystal_elastic_series,
)
from quantas.models.elastic_states import ElasticTensorKind, PressureSource


def _write_output(
    path: Path,
    *,
    volume: float,
    density_g_cm3: float,
    energy_hartree: float,
    stress_pressure_gpa: float | None,
    crystal_pressure_gpa: float | None = None,
    corrected_total_hartree: float | None = None,
) -> Path:
    """Write the minimal completed CRYSTAL sections required by the reader."""
    lines = [
        "ELAPIEZO OPTION",
        f"VOLUME OF THE CELL: {volume:.8f}",
        f"DENSITY OF THE CRYSTAL = {density_g_cm3:.8f}",
        f"TOTAL ENERGY(DFT)(AU)( 12) {energy_hartree:.12f}",
    ]
    if corrected_total_hartree is not None:
        lines.append(
            "TOTAL ENERGY + DISP (AU) " f"{corrected_total_hartree:.12f}"
        )
    if stress_pressure_gpa is not None:
        lines.append(f"PRESSURE IN GIGAPASCAL: {stress_pressure_gpa:.8f}")
    if crystal_pressure_gpa is not None:
        lines.extend(
            [
                "PRESSURE",
                f"{crystal_pressure_gpa:.8f}",
                f"ELASTIC PROPERTIES AT PRESSURE {crystal_pressure_gpa:.8f}",
            ]
        )
    lines.extend(
        [
            "FINAL RESULTS START",
            "SYMMETRIZED ELASTIC CONSTANTS",
            "header",
            "200 80 70 0 0 0",
            "190 65 0 0 0",
            "180 0 0 0",
            "60 0 0",
            "55 0",
            "50",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_import_uses_corrected_crystal_total_energy(tmp_path) -> None:
    """Elastic states must use CRYSTAL's corrected total rather than raw SCF energy."""
    output = _write_output(
        tmp_path / "d3.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.0,
        corrected_total_hartree=-10.025,
        stress_pressure_gpa=0.0,
    )

    state = read_crystal_elastic_series([output]).states[0]

    assert state.energy == pytest.approx(-10.025)
    energy = state.metadata["energy"]
    assert energy["selected_quantity"] == "total_energy"
    assert energy["scf_energy_hartree"] == pytest.approx(-10.0)
    assert energy["total_energy_hartree"] == pytest.approx(-10.025)
    assert energy["source_marker"] == "TOTAL ENERGY + DISP (AU)"
    assert energy["corrections"] == ["DISP"]


def test_crystal_pressure_correction_matches_erba_matrix() -> None:
    """CRYSTAL raw coefficients follow Erba et al. Eq. (6)--(7)."""
    raw = np.diag([200.0, 210.0, 220.0, 70.0, 75.0, 80.0])

    compressed = crystal_hydrostatic_stiffness(raw, 2.0)
    expected_compressed = raw.copy()
    expected_compressed[0, 1] = expected_compressed[1, 0] = 2.0
    expected_compressed[0, 2] = expected_compressed[2, 0] = 2.0
    expected_compressed[1, 2] = expected_compressed[2, 1] = 2.0
    expected_compressed[3, 3] = 69.0
    expected_compressed[4, 4] = 74.0
    expected_compressed[5, 5] = 79.0
    np.testing.assert_allclose(compressed, expected_compressed)

    tensile = crystal_hydrostatic_stiffness(raw, -2.0)
    expected_tensile = raw.copy()
    expected_tensile[0, 1] = expected_tensile[1, 0] = -2.0
    expected_tensile[0, 2] = expected_tensile[2, 0] = -2.0
    expected_tensile[1, 2] = expected_tensile[2, 1] = -2.0
    expected_tensile[3, 3] = 71.0
    expected_tensile[4, 4] = 76.0
    expected_tensile[5, 5] = 81.0
    np.testing.assert_allclose(tensile, expected_tensile)
    np.testing.assert_allclose(crystal_hydrostatic_stiffness(raw, 0.0), raw)


def test_import_sorts_states_selects_reference_and_corrects(tmp_path) -> None:
    """Raw outputs become an increasing, Wallace-corrected volume series."""
    large = _write_output(
        tmp_path / "large.out",
        volume=110.0,
        density_g_cm3=3.0,
        energy_hartree=-100.2,
        stress_pressure_gpa=-1.5,
    )
    small = _write_output(
        tmp_path / "small.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-100.0,
        stress_pressure_gpa=2.0,
    )

    series = read_crystal_elastic_series([large, small])

    np.testing.assert_allclose(series.volumes, [100.0, 110.0])
    assert series.reference_index == 1
    assert series.orientation == "crystal-cartesian"
    assert series.states[0].density == pytest.approx(3300.0)
    assert (
        series.states[0].prestress.tensor_kind is ElasticTensorKind.WALLACE_HYDROSTATIC
    )
    assert series.states[0].prestress.pressure_source is PressureSource.OUTPUT_STRESS
    assert series.states[0].prestress.correction_applied_by == "quantas-crystal-import"
    assert series.states[0].stiffness[0, 0] == pytest.approx(200.0)
    assert series.states[0].stiffness[0, 1] == pytest.approx(82.0)
    assert series.states[0].stiffness[3, 3] == pytest.approx(59.0)
    assert series.states[0].prestress.correction_method == (
        "crystal-erba-2014-hydrostatic"
    )


def test_manual_pressure_follows_input_order_before_volume_sort(tmp_path) -> None:
    """Manual pressures remain associated with their source output."""
    large = _write_output(
        tmp_path / "large.out",
        volume=110.0,
        density_g_cm3=3.0,
        energy_hartree=-10.1,
        stress_pressure_gpa=None,
    )
    small = _write_output(
        tmp_path / "small.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.0,
        stress_pressure_gpa=None,
    )

    series = read_crystal_elastic_series(
        [large, small],
        pressure_policy=CrystalPressurePolicy.MANUAL,
        manual_pressures_gpa=[-2.0, 3.0],
        apply_prestress_correction=False,
    )

    assert [state.prestress.pressure_gpa for state in series.states] == [3.0, -2.0]
    assert all(
        state.prestress.pressure_source is PressureSource.MANUAL
        for state in series.states
    )
    assert all(
        state.prestress.tensor_kind is ElasticTensorKind.RAW_ENERGY_STRAIN
        for state in series.states
    )


def test_auto_preserves_crystal_pressure_correction(tmp_path) -> None:
    """A CRYSTAL-corrected tensor is marked incremental without recorrection."""
    output = _write_output(
        tmp_path / "corrected.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.0,
        stress_pressure_gpa=2.1,
        crystal_pressure_gpa=2.0,
    )

    state = read_crystal_elastic_series([output]).states[0]

    assert state.prestress.tensor_kind is ElasticTensorKind.WALLACE_HYDROSTATIC
    assert state.prestress.pressure_gpa == pytest.approx(2.0)
    assert state.prestress.pressure_source is PressureSource.APPLIED_PRESTRESS
    assert state.prestress.correction_applied_by == "crystal"
    assert state.stiffness[0, 0] == pytest.approx(200.0)



def test_auto_preserves_crystal_presseos_correction(tmp_path) -> None:
    """PRESSEOS is recognized as backend-applied Wallace provenance."""
    output = _write_output(
        tmp_path / "presseos.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.0,
        stress_pressure_gpa=2.1,
        crystal_pressure_gpa=2.0,
    )
    text = output.read_text(encoding="utf-8").replace("\nPRESSURE\n", "\nPRESSEOS\n")
    output.write_text(text, encoding="utf-8")

    state = read_crystal_elastic_series([output]).states[0]

    assert state.prestress.tensor_kind is ElasticTensorKind.WALLACE_HYDROSTATIC
    assert state.prestress.pressure_gpa == pytest.approx(2.0)
    assert state.prestress.correction_method == "crystal-presseos-keyword"
    assert state.metadata["prestress_keyword"] == "PRESSEOS"

def test_missing_stress_requests_manual_pressure(tmp_path) -> None:
    """Raw tensors without output stress do not acquire an implicit pressure."""
    output = _write_output(
        tmp_path / "missing-pressure.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.0,
        stress_pressure_gpa=None,
    )

    with pytest.raises(ValueError, match="pressure_policy='manual'"):
        read_crystal_elastic_series([output])


def test_deferred_pressure_retains_an_explicitly_raw_tensor(tmp_path) -> None:
    """External E(V) workflows can parse raw tensors before assigning pressure."""
    output = _write_output(
        tmp_path / "deferred.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.0,
        stress_pressure_gpa=None,
    )

    state = read_crystal_elastic_series(
        [output],
        pressure_policy=CrystalPressurePolicy.DEFERRED,
        apply_prestress_correction=False,
    ).states[0]

    assert state.prestress.tensor_kind is ElasticTensorKind.RAW_ENERGY_STRAIN
    assert state.prestress.pressure_gpa is None
    assert state.prestress.pressure_source is PressureSource.UNAVAILABLE

    with pytest.raises(ValueError, match="apply_prestress_correction=False"):
        read_crystal_elastic_series(
            [output],
            pressure_policy=CrystalPressurePolicy.DEFERRED,
        )


def test_manual_pressure_count_and_duplicate_volumes_are_rejected(tmp_path) -> None:
    """Ambiguous pressure and volume associations fail before use."""
    first = _write_output(
        tmp_path / "first.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.0,
        stress_pressure_gpa=1.0,
    )
    second = _write_output(
        tmp_path / "second.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.1,
        stress_pressure_gpa=1.0,
    )

    with pytest.raises(ValueError, match="one pressure per"):
        read_crystal_elastic_series(
            [first, second],
            pressure_policy="manual",
            manual_pressures_gpa=[1.0],
        )
    with pytest.raises(ValueError, match="duplicate volumes"):
        read_crystal_elastic_series([first, second])


def test_inconsistent_parsed_lattice_is_recorded_but_not_attached(
    tmp_path,
    monkeypatch,
) -> None:
    """The final elastic volume remains authoritative over a stale cell block."""
    output = _write_output(
        tmp_path / "state.out",
        volume=100.0,
        density_g_cm3=3.3,
        energy_hartree=-10.0,
        stress_pressure_gpa=1.0,
    )
    from quantas.models.structures import CrystalStructure

    structure = CrystalStructure(
        lattice=np.diag([5.0, 5.0, 5.0]),
        fractional_positions=np.zeros((1, 3)),
        atomic_numbers=np.array([1]),
    )

    monkeypatch.setattr(
        "quantas.interfaces.crystal.elasticity."
        "CrystalElasticityReader._read_reference_structure",
        staticmethod(lambda _lines, *, volume: (structure, 0)),
    )

    state = read_crystal_elastic_series([output]).states[0]

    assert state.lattice is None
    assert state.metadata["parsed_structure_volume_angstrom3"] == pytest.approx(125.0)
    assert state.metadata["parsed_structure_matches_elastic_volume"] is False
