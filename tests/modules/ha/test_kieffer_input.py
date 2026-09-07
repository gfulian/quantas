"""Tests for Kieffer enrichment of phonon YAML inputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml
from click.testing import CliRunner

from quantas.api import ha, qha
from quantas.cli.ha import ha as ha_command
from quantas.cli.qha import qha as qha_command
from quantas.core.physics.eos import EnergyEOS
from quantas.core.physics.units import energy_to_pressure
from quantas.io.kieffer import (
    kieffer_series_from_mapping,
    kieffer_series_to_mapping,
)
from quantas.models.kieffer import KiefferCutoffState, KiefferVolumeSeries


def _phonon_input(path: Path, volumes: list[float]) -> Path:
    """Write a primitive Gamma-only phonon input."""
    data = {
        "job": "Kieffer enrichment test",
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
        "energy": [-10.0 - 0.01 * index for index in range(len(volumes))],
        "phonon": [
            {
                "q-position": [0.0, 0.0, 0.0],
                "weight": 1.0,
                "band": [
                    {"frequency": [value] * len(volumes)} for value in (0.0, 0.0, 0.0)
                ],
            }
        ],
    }
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return path


def _elastic_output(path: Path, volume: float, energy: float, pressure: float) -> Path:
    """Write compact completed ELAPIEZO output data."""
    rows = (
        "240 80 80 0 0 0",
        "240 80 0 0 0",
        "240 0 0 0",
        "80 0 0",
        "80 0",
        "80",
    )
    lines = [
        "ELAPIEZO OPTION",
        f"VOLUME OF THE CELL: {volume}",
        "DENSITY OF THE CRYSTAL = 3.2",
        f"TOTAL ENERGY(DFT)(AU)( 10) {energy}",
        f"PRESSURE IN GIGAPASCAL: {pressure}",
        "FINAL RESULTS START",
        "SYMMETRIZED ELASTIC CONSTANTS",
        "header",
        *rows,
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _direct_cutoffs(volumes: list[float]) -> KiefferVolumeSeries:
    """Return increasing direct cutoff states for selected volumes."""
    states = []
    for index, volume in enumerate(sorted(volumes)):
        scale = (10.0 / volume) ** (1.0 / 3.0)
        states.append(
            KiefferCutoffState(
                volume=volume,
                frequencies_hz=scale * np.array([3.0e12, 4.0e12, 5.0e12]),
                effective_velocities_km_s=np.array([3.0, 4.0, 7.0]),
                source_elastic_indices=(index,),
            )
        )
    return KiefferVolumeSeries(states=tuple(states))


def test_kieffer_mapping_round_trip_preserves_direct_states() -> None:
    """The YAML contract restores cutoff and velocity arrays exactly."""
    series = KiefferVolumeSeries(
        states=(
            KiefferCutoffState(
                volume=80.0,
                frequencies_hz=[3.0e12, 4.0e12, 5.0e12],
                effective_velocities_km_s=[3.0, 4.0, 7.0],
                source_elastic_indices=(0,),
                metadata={"pressure_gpa": 2.0},
            ),
        ),
        metadata={"source": "test"},
    )

    restored = kieffer_series_from_mapping(
        kieffer_series_to_mapping(series, provenance={"backend": "crystal"})
    )

    np.testing.assert_array_equal(restored.volumes, series.volumes)
    np.testing.assert_array_equal(restored.frequencies_hz, series.frequencies_hz)
    np.testing.assert_array_equal(
        restored.effective_velocities_km_s,
        series.effective_velocities_km_s,
    )
    assert restored.metadata["input_provenance"]["backend"] == "crystal"


def test_kieffer_mapping_canonicalizes_flexible_constructor_types() -> None:
    """Serialization narrows ArrayLike and string fields to canonical values."""
    series = KiefferVolumeSeries(
        states=(
            KiefferCutoffState(
                volume=np.float32(80.0),
                frequencies_hz=(3.0e12, 4.0e12, 5.0e12),
                effective_velocities_km_s=(3, 4, 7),
                source="direct",
                source_elastic_indices=(0,),
            ),
        )
    )

    state = kieffer_series_to_mapping(series)["states"][0]

    assert state["cutoff_frequency"] == [3.0e12, 4.0e12, 5.0e12]
    assert state["effective_velocity"] == [3.0, 4.0, 7.0]
    assert state["source"] == "direct"


def test_public_ha_api_creates_new_enriched_input(tmp_path) -> None:
    """HA enrichment leaves its source untouched and embeds one cutoff state."""
    source = _phonon_input(tmp_path / "ha.yaml", [80.0])
    original = source.read_bytes()
    elastic = _elastic_output(tmp_path / "ha-elastic.out", 80.0, -10.0, 2.0)
    destination = tmp_path / "ha-kieffer.yaml"

    output = ha.add_kieffer_input(
        source,
        destination,
        [elastic],
        mu_order=2,
        phi_order=4,
        refinement_factor=2,
    )

    assert output == destination
    assert source.read_bytes() == original
    cutoffs = ha.read_kieffer_input(destination)
    assert len(cutoffs.states) == 1
    assert np.all(cutoffs.frequencies_hz > 0.0)
    raw = yaml.safe_load(destination.read_text(encoding="utf-8"))
    assert raw["kieffer"]["composition"] == "additional-acoustic-branches"
    assert raw["kieffer"]["states"][0]["metadata"]["pressure_gpa"] == 2.0


def test_public_qha_api_matches_outputs_independently_of_file_order(tmp_path) -> None:
    """QHA enrichment writes increasing cutoff states and validates all volumes."""
    source = _phonon_input(tmp_path / "qha.yaml", [80.0, 90.0])
    large = _elastic_output(tmp_path / "large.out", 90.0, -10.1, -1.0)
    small = _elastic_output(tmp_path / "small.out", 80.0, -10.0, 2.0)
    destination = tmp_path / "qha-kieffer.yaml"

    qha.add_kieffer_input(
        source,
        destination,
        [large, small],
        mu_order=2,
        phi_order=4,
        refinement_factor=2,
    )

    cutoffs = qha.read_kieffer_input(destination)
    np.testing.assert_array_equal(cutoffs.volumes, [80.0, 90.0])
    assert all(np.all(state.frequencies_hz > 0.0) for state in cutoffs.states)


def test_qha_enrichment_derives_pressure_from_energy_polynomial(tmp_path) -> None:
    """Polynomial P(V) values and complete fit provenance are stored in YAML."""
    volumes = [80.0, 85.0, 90.0, 95.0, 100.0]
    source = _phonon_input(tmp_path / "qha.yaml", volumes)
    raw = yaml.safe_load(source.read_text(encoding="utf-8"))
    raw["energy"] = [-10.0 + 1.0e-5 * (volume - 90.0) ** 2 for volume in volumes]
    source.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    outputs = [
        _elastic_output(
            tmp_path / f"elastic-{index}.out",
            volume,
            -10.0,
            99.0,
        )
        for index, volume in enumerate(reversed(volumes))
    ]
    destination = tmp_path / "qha-kieffer.yaml"

    qha.add_kieffer_input(
        source,
        destination,
        outputs,
        interface="crystal",
        pressure_policy="energy_polynomial",
        polynomial_degree=2,
        mu_order=2,
        phi_order=4,
        refinement_factor=2,
    )

    enriched = yaml.safe_load(destination.read_text(encoding="utf-8"))
    provenance = enriched["kieffer"]["provenance"]
    expected = energy_to_pressure(
        -2.0e-5 * (np.asarray(volumes) - 90.0),
        "Ha",
        "angstrom",
        "GPa",
    )
    actual = [
        state["metadata"]["pressure_gpa"] for state in enriched["kieffer"]["states"]
    ]
    np.testing.assert_allclose(actual, expected, rtol=2.0e-11, atol=1.0e-10)
    assert provenance["elastic_interface"] == "crystal"
    assert provenance["pressure_source"] == "energy_polynomial"
    assert provenance["pressure_model"]["settings"] == {"degree": 2}
    assert provenance["pressure_model"]["relation"] == "P(V) = -dE/dV"
    assert provenance["pressure_model"]["fit"]["success"] is True
    assert len(provenance["pressure_model"]["volume_matches"]) == len(volumes)
    assert all(
        state["metadata"]["pressure_source"] == "energy_polynomial"
        for state in enriched["kieffer"]["states"]
    )


def test_qha_enrichment_derives_pressure_from_energy_eos(tmp_path) -> None:
    """EOS selection is canonicalized and retained beside evaluated pressures."""
    volumes = np.linspace(80.0, 100.0, 5)
    source = _phonon_input(tmp_path / "qha.yaml", volumes.tolist())
    parameters = np.array([-10.0, 0.005, 4.2, 90.0])
    model = EnergyEOS()
    raw = yaml.safe_load(source.read_text(encoding="utf-8"))
    raw["energy"] = model.evaluate("BM3", volumes, parameters).tolist()
    source.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    outputs = [
        _elastic_output(
            tmp_path / f"elastic-{index}.out",
            float(volume),
            -10.0,
            99.0,
        )
        for index, volume in enumerate(volumes)
    ]
    destination = tmp_path / "qha-kieffer.yaml"

    qha.add_kieffer_input(
        source,
        destination,
        outputs,
        pressure_policy="energy-eos",
        eos="BM3",
        mu_order=2,
        phi_order=4,
        refinement_factor=2,
    )

    enriched = yaml.safe_load(destination.read_text(encoding="utf-8"))
    pressure_model = enriched["kieffer"]["provenance"]["pressure_model"]
    expected = energy_to_pressure(
        model.pressure("BM3", parameters, volumes),
        "Ha",
        "angstrom",
        "GPa",
    )
    np.testing.assert_allclose(
        pressure_model["evaluated_pressures_gpa"],
        expected,
        rtol=2.0e-8,
        atol=1.0e-7,
    )
    assert pressure_model["method"] == "energy_eos"
    assert pressure_model["settings"] == {
        "eos": "BM3",
        "eos_family": "birchmurnaghan",
        "eos_order": 3,
    }


def test_ha_rejects_energy_derived_pressure(tmp_path) -> None:
    """A single-volume HA dataset cannot define an E(V) derivative."""
    source = _phonon_input(tmp_path / "ha.yaml", [80.0])
    elastic = _elastic_output(tmp_path / "elastic.out", 80.0, -10.0, 1.0)

    with pytest.raises(ValueError, match="multi-volume QHA"):
        ha.add_kieffer_input(
            source,
            tmp_path / "ha-kieffer.yaml",
            [elastic],
            pressure_policy="energy_eos",
        )


def test_enrichment_rejects_source_overwrite_and_existing_block(tmp_path) -> None:
    """Neither implicit in-place writes nor silent Kieffer replacement occur."""
    source = _phonon_input(tmp_path / "ha.yaml", [80.0])
    elastic = _elastic_output(tmp_path / "elastic.out", 80.0, -10.0, 2.0)
    with pytest.raises(ValueError, match="new output path"):
        ha.add_kieffer_input(source, source, [elastic])

    destination = tmp_path / "ha-kieffer.yaml"
    ha.add_kieffer_input(
        source,
        destination,
        [elastic],
        mu_order=2,
        phi_order=4,
        refinement_factor=2,
    )
    with pytest.raises(ValueError, match="already contains"):
        ha.add_kieffer_input(
            destination,
            tmp_path / "second.yaml",
            [elastic],
            mu_order=2,
            phi_order=4,
            refinement_factor=2,
        )


def test_ha_cli_runs_embedded_kieffer_contribution_end_to_end(tmp_path) -> None:
    """HA uses embedded cutoffs only when the execution flag is present."""
    filename = _phonon_input(tmp_path / "ha.yaml", [80.0])
    raw = yaml.safe_load(filename.read_text(encoding="utf-8"))
    raw["kieffer"] = kieffer_series_to_mapping(_direct_cutoffs([80.0]))
    filename.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    output = tmp_path / "ha-result.hdf5"

    response = CliRunner().invoke(
        ha_command,
        [
            "run",
            str(filename),
            "--kieffer",
            "--temperature",
            "0",
            "300",
            "300",
            "--output",
            str(output),
            "--force",
            "--quiet",
            "--no-progress",
        ],
    )

    assert response.exit_code == 0, response.output
    result = ha.get_result(ha.read_result(output))
    assert result.kieffer_contribution is not None
    assert result.metadata["kieffer"]["composition"] == ("additional-acoustic-branches")

    baseline_output = tmp_path / "ha-baseline.hdf5"
    baseline_response = CliRunner().invoke(
        ha_command,
        [
            "run",
            str(filename),
            "--temperature",
            "0",
            "300",
            "300",
            "--output",
            str(baseline_output),
            "--force",
            "--quiet",
            "--no-progress",
        ],
    )
    assert baseline_response.exit_code == 0, baseline_response.output
    baseline = ha.get_result(ha.read_result(baseline_output))
    assert baseline.kieffer_contribution is None


@pytest.mark.parametrize("scheme", ["freq", "td"])
def test_qha_cli_runs_kieffer_end_to_end(tmp_path, scheme: str) -> None:
    """Both QHA schemes pass embedded cutoffs through minimization and HDF5."""
    volumes = [10.0, 9.0, 11.0, 9.5, 10.5]
    filename = _phonon_input(tmp_path / "qha.yaml", volumes)
    raw = yaml.safe_load(filename.read_text(encoding="utf-8"))
    raw["energy"] = [0.02 * (volume - 10.0) ** 2 - 10.0 for volume in volumes]
    raw["mode_continuity"] = "assumed"
    raw["kieffer"] = kieffer_series_to_mapping(_direct_cutoffs(volumes))
    filename.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    output = tmp_path / "qha-result.hdf5"

    response = CliRunner().invoke(
        qha_command,
        [
            "run",
            str(filename),
            "--kieffer",
            "--scheme",
            scheme,
            "--energy-degree",
            "2",
            "--frequency-degree",
            "2",
            "--temperature",
            "0",
            "300",
            "300",
            "--pressure",
            "0",
            "0",
            "1",
            "--output",
            str(output),
            "--force",
            "--quiet",
            "--no-progress",
        ],
    )

    assert response.exit_code == 0, response.output
    result = qha.get_result(qha.read_result(output))
    assert result.completed is True
    assert result.kieffer_sampled_contribution is not None
    assert result.metadata["kieffer"]["composition"] == ("additional-acoustic-branches")
