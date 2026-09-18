"""Backend-neutral VASP elastic-state series characterization tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.interfaces.vasp import (
    VASP_HYDROSTATIC_PRESTRESS_METHOD,
    assign_vasp_manual_pressures,
    convert_vasp_hydrostatic_elastic_series,
    convert_vasp_hydrostatic_elastic_state,
    read_vasp_elastic_series,
    resolve_vasp_energy_derived_pressures,
    vasp_hydrostatic_incremental_stiffness,
)
from quantas.models.elastic_states import ElasticTensorKind, PressureSource


DATA = Path(__file__).parent / "data"


def _write_volume_variant(source: Path, destination: Path, volume: float) -> None:
    """Write one OUTCAR fixture variant with a different reference volume."""
    text = source.read_text(encoding="utf-8")
    text = text.replace("volume of cell :       19.28", f"volume of cell :       {volume:.2f}")
    destination.write_text(text, encoding="utf-8")


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_raw_elastic_series_preserves_reader_facts(tmp_path: Path) -> None:
    """The VASP series adapter preserves raw tensor and pressure provenance."""
    run = tmp_path / "state00"
    run.mkdir()
    (run / "OUTCAR").write_bytes((DATA / "vasp_mgo_soec_00_v544.OUTCAR").read_bytes())

    series = read_vasp_elastic_series([run])

    assert series.nstates == 1
    assert series.reference_index == 0
    assert series.orientation == "vasp-cartesian"
    assert series.metadata["backend"] == "vasp"
    assert series.metadata["prestress_correction_applied"] is False
    state = series.states[0]
    assert state.volume == pytest.approx(19.28)
    assert state.density > 0.0
    assert state.stiffness[0, 0] == pytest.approx(271.25302)
    assert state.stiffness[0, 1] == pytest.approx(87.63414)
    assert state.prestress.tensor_kind is ElasticTensorKind.RAW_STRESS_STRAIN
    assert state.prestress.pressure_gpa == pytest.approx(-0.781768)
    assert state.prestress.pressure_source is PressureSource.OUTPUT_STRESS
    assert state.metadata["prestress_correction_applied"] is False
    assert state.metadata["source_index"] == 0

    with pytest.raises(ValueError, match="incremental"):
        series.require_incremental()


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_raw_series_preserves_reference_after_sorting(
    tmp_path: Path,
) -> None:
    """Volume sorting must not change the reference state chosen in source order."""
    source = DATA / "vasp_mgo_soec_00_v544.OUTCAR"
    high = tmp_path / "high.OUTCAR"
    low = tmp_path / "low.OUTCAR"
    _write_volume_variant(source, high, 21.00)
    _write_volume_variant(source, low, 18.00)

    series = read_vasp_elastic_series(
        [high, low],
        reference_source_index=0,
    )

    np.testing.assert_allclose(series.volumes, [18.0, 21.0])
    assert series.reference_index == 1
    assert series.states[0].metadata["source_index"] == 1
    assert series.states[1].metadata["source_index"] == 0


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_raw_elastic_series_rejects_duplicate_sources() -> None:
    """One calculation source may not appear twice in an elastic series."""
    source = DATA / "vasp_mgo_soec_00_v544.OUTCAR"
    with pytest.raises(ValueError, match="must be unique"):
        read_vasp_elastic_series([source, source])


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_pressure_adjustment_matches_documented_matrix() -> None:
    """The VASP pressure adjustment follows the documented Voigt relation."""
    raw = np.asarray(
        [
            [100.0, 10.0, 20.0, 1.0, 2.0, 3.0],
            [10.0, 110.0, 25.0, 4.0, 5.0, 6.0],
            [20.0, 25.0, 120.0, 7.0, 8.0, 9.0],
            [1.0, 4.0, 7.0, 30.0, 0.5, 0.6],
            [2.0, 5.0, 8.0, 0.5, 40.0, 0.7],
            [3.0, 6.0, 9.0, 0.6, 0.7, 50.0],
        ],
        dtype=np.float64,
    )

    corrected = vasp_hydrostatic_incremental_stiffness(raw, 5.0)

    expected = raw.copy()
    expected[np.diag_indices(6)] -= 5.0
    for first, second in ((0, 1), (0, 2), (1, 2)):
        expected[first, second] += 5.0
        expected[second, first] += 5.0
    np.testing.assert_allclose(corrected, expected, rtol=0.0, atol=0.0)
    assert corrected[3, 4] == raw[3, 4]
    assert corrected[0, 3] == raw[0, 3]


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_output_stress_conversion_matches_mgo_fixture(tmp_path: Path) -> None:
    """The MgO fixture converts once from raw VASP to incremental stiffness."""
    run = tmp_path / "state00"
    run.mkdir()
    (run / "OUTCAR").write_bytes((DATA / "vasp_mgo_soec_00_v544.OUTCAR").read_bytes())
    raw_series = read_vasp_elastic_series([run])
    raw_state = raw_series.states[0]

    corrected = convert_vasp_hydrostatic_elastic_state(raw_state)

    assert raw_state.prestress.tensor_kind is ElasticTensorKind.RAW_STRESS_STRAIN
    assert corrected.prestress.tensor_kind is ElasticTensorKind.WALLACE_HYDROSTATIC
    assert corrected.prestress.source_tensor_kind is ElasticTensorKind.RAW_STRESS_STRAIN
    assert corrected.prestress.correction_method == VASP_HYDROSTATIC_PRESTRESS_METHOD
    assert corrected.prestress.pressure_gpa == pytest.approx(-0.781768)
    assert corrected.stiffness[0, 0] == pytest.approx(272.034788)
    assert corrected.stiffness[0, 1] == pytest.approx(86.852372)
    assert corrected.stiffness[3, 3] == pytest.approx(141.729318)
    assert corrected.metadata["prestress_correction_applied"] is True
    assert corrected.metadata["prestress_correction"][
        "hydrostatic_stress_residual_gpa"
    ] == pytest.approx(0.0)

    with pytest.raises(ValueError, match="raw stress-strain"):
        convert_vasp_hydrostatic_elastic_state(corrected)


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_hydrostatic_series_becomes_incremental(tmp_path: Path) -> None:
    """Only the explicitly converted VASP series may pass the incremental gate."""
    run = tmp_path / "state00"
    run.mkdir()
    (run / "OUTCAR").write_bytes((DATA / "vasp_mgo_soec_00_v544.OUTCAR").read_bytes())
    raw_series = read_vasp_elastic_series([run])

    with pytest.raises(ValueError, match="incremental"):
        raw_series.require_incremental()

    corrected = convert_vasp_hydrostatic_elastic_series(raw_series)
    corrected.require_incremental()
    assert corrected.metadata["prestress_correction_applied"] is True
    assert corrected.states[0].prestress.tensor_kind is ElasticTensorKind.WALLACE_HYDROSTATIC


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_conversion_rejects_deviatoric_stress(
    tmp_path: Path,
) -> None:
    """The scalar hydrostatic correction must not hide deviatoric pre-stress."""
    run = tmp_path / "state00"
    run.mkdir()
    (run / "OUTCAR").write_bytes((DATA / "vasp_mgo_soec_00_v544.OUTCAR").read_bytes())
    state = read_vasp_elastic_series([run]).states[0]
    stress = np.asarray(state.metadata["reference_stress_gpa"], dtype=np.float64)
    stress[0, 0] += 0.02
    state.metadata["reference_stress_gpa"] = stress.tolist()

    with pytest.raises(ValueError, match="not hydrostatic"):
        convert_vasp_hydrostatic_elastic_state(state, hydrostatic_atol_gpa=0.01)


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_manual_pressure_override_preserves_raw_tensor(
    tmp_path: Path,
) -> None:
    """Manual pressure replacement is separate from VASP stiffness conversion."""
    run = tmp_path / "state00"
    run.mkdir()
    (run / "OUTCAR").write_bytes(
        (DATA / "vasp_mgo_soec_00_v544.OUTCAR").read_bytes()
    )
    raw = read_vasp_elastic_series([run])

    assigned = assign_vasp_manual_pressures(
        raw,
        [0.5],
        assignment_method="manual-test",
    )

    np.testing.assert_allclose(assigned.states[0].stiffness, raw.states[0].stiffness)
    assert assigned.states[0].prestress.tensor_kind is ElasticTensorKind.RAW_STRESS_STRAIN
    assert assigned.states[0].prestress.pressure_gpa == pytest.approx(0.5)
    assert assigned.states[0].prestress.pressure_source is PressureSource.MANUAL
    assignment = assigned.states[0].metadata["pressure_assignment"]
    assert assignment["replaced_pressure_gpa"] == pytest.approx(-0.781768)
    assert assignment["replaced_pressure_source"] == "output_stress"
    with pytest.raises(ValueError, match="incremental"):
        assigned.require_incremental()

    corrected = convert_vasp_hydrostatic_elastic_series(assigned)
    corrected.require_incremental()
    assert corrected.states[0].stiffness[0, 0] == pytest.approx(270.75302)
    correction = corrected.states[0].metadata["prestress_correction"]
    assert correction["pressure_source"] == "manual"
    assert correction["output_stress_pressure_gpa"] == pytest.approx(-0.781768)
    assert correction["pressure_minus_output_stress_gpa"] == pytest.approx(1.281768)


@pytest.mark.interfaces
@pytest.mark.elasticity
def test_vasp_energy_pressure_resolution_replaces_output_pressure(
    tmp_path: Path,
) -> None:
    """Backend-neutral E(V) pressure fitting may feed raw VASP conversion."""
    source = DATA / "vasp_mgo_soec_00_v544.OUTCAR"
    volumes = np.asarray([18.0, 19.0, 20.0, 21.0], dtype=np.float64)
    runs: list[Path] = []
    for index, volume in enumerate(volumes):
        path = tmp_path / f"state{index}.OUTCAR"
        _write_volume_variant(source, path, float(volume))
        runs.append(path)
    raw = read_vasp_elastic_series(runs, reference_source_index=1)
    energies = -10.0 + 0.01 * (volumes - 19.5) ** 2

    resolution = resolve_vasp_energy_derived_pressures(
        raw,
        volumes,
        energies,
        pressure_source=PressureSource.ENERGY_POLYNOMIAL,
        source_dataset="synthetic-v-e",
        energy_unit="eV",
        volume_length_unit="angstrom",
        volume_unit="angstrom^3",
        polynomial_degree=2,
    )

    assert resolution.series.reference_index == raw.reference_index
    np.testing.assert_allclose(resolution.series.stiffness, raw.stiffness)
    assert all(
        state.prestress.pressure_source is PressureSource.ENERGY_POLYNOMIAL
        for state in resolution.series.states
    )
    assignment = resolution.series.states[0].metadata["pressure_assignment"]
    assert assignment["replaced_pressure_source"] == "output_stress"
    assert resolution.provenance["source_dataset"] == "synthetic-v-e"

    corrected = convert_vasp_hydrostatic_elastic_series(resolution.series)
    corrected.require_incremental()
    assert all(
        state.prestress.pressure_source is PressureSource.ENERGY_POLYNOMIAL
        for state in corrected.states
    )
