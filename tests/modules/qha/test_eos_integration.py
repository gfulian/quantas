from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.units import pressure_to_energy
from quantas.core.physics.eos import FittedEnergyEOS
from quantas.core.physics.eos import EnergyEOS
from quantas.modules.qha.models import QHAInput, QHAOptions
from quantas.modules.qha.workflow import run_volume_minimization_workflow


def _eos_input(
    *, temperatures: int = 1, noisy: bool = False
) -> tuple[QHAInput, np.ndarray]:
    volume = np.linspace(68.0, 76.0, 13)
    k0_gpa = 160.0
    k0_native = float(pressure_to_energy(k0_gpa, "eV", "A", "GPa"))
    parameters = np.array([-100.0, k0_native, 4.2, 72.0], dtype=np.float64)
    energy = EnergyEOS().birchmurnaghan(volume, *parameters)
    if noisy:
        energy = (
            energy
            + np.array(
                [0.0, 1.0, -0.5, 0.8, -0.7, 0.2, 0.0, -0.2, 0.6, -0.8, 0.5, -1.0, 0.0]
            )
            * 1.0e-5
        )
    if temperatures == 1:
        free_energy = energy[np.newaxis, :]
    else:
        free_energy = np.vstack(
            [energy + index * 1.0e-4 for index in range(temperatures)]
        )
    input_data = QHAInput(jobname="eos-workflow", volume=volume, energy=energy)
    return input_data, free_energy


def _options(**kwargs: object) -> QHAOptions:
    values = dict(
        temperature_min=300.0,
        temperature_max=300.0,
        temperature_step=100.0,
        pressure_min=0.0,
        pressure_max=20.0,
        pressure_step=10.0,
        minimization="eos",
        eos="BM",
        energy_unit="eV",
        volume_unit="A",
        pressure_unit="GPa",
        debug=True,
    )
    values.update(kwargs)
    return QHAOptions(**values)


def test_eos_workflow_fits_once_per_temperature(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_data, free_energy = _eos_input(temperatures=2)
    options = _options(
        temperature_max=400.0,
        estimate_uncertainties=False,
    )
    calls = 0
    original = EnergyEOS.fit

    def counted_fit(self, *args, **kwargs):
        nonlocal calls
        calls += 1
        return original(self, *args, **kwargs)

    monkeypatch.setattr(EnergyEOS, "fit", counted_fit)

    result = run_volume_minimization_workflow(
        input_data,
        options,
        free_energy=free_energy,
    )

    assert result.completed is True
    assert calls == 2
    assert result.valid_mask.all()
    assert len(result.fit_records) == 2
    assert all(record.quantity == "F" for record in result.fit_records)
    assert all(record.pressure is None for record in result.fit_records)
    assert result.metadata["eos_workflow"]["fit_count"] == 2
    assert result.metadata["eos_workflow"]["state_count"] == 6


def test_eos_workflow_evaluates_all_pressures_from_one_fitted_model() -> None:
    input_data, free_energy = _eos_input()
    options = _options(estimate_uncertainties=False)

    result = run_volume_minimization_workflow(
        input_data,
        options,
        free_energy=free_energy,
    )

    k0_gpa = 160.0
    reference = FittedEnergyEOS(
        "BM",
        [-100.0, k0_gpa, 4.2, 72.0],
        sampled_volumes=input_data.volume,
    )
    expected = reference.states_at_pressures(result.pressure)

    np.testing.assert_allclose(
        result.equilibrium_volume[0],
        [state.volume for state in expected],
        rtol=1.0e-8,
        atol=1.0e-8,
    )
    np.testing.assert_allclose(
        result.isothermal_bulk_modulus[0],
        [state.bulk_modulus for state in expected],
        rtol=1.0e-8,
        atol=1.0e-8,
    )
    np.testing.assert_allclose(
        result.bulk_modulus_derivative[0],
        [state.bulk_modulus_derivative for state in expected],
        rtol=1.0e-8,
        atol=1.0e-8,
    )


def test_eos_workflow_populates_fit_uncertainties() -> None:
    input_data, free_energy = _eos_input(noisy=True)
    options = _options(
        pressure_max=10.0,
        estimate_uncertainties=True,
        uncertainty_method="covariance",
    )

    result = run_volume_minimization_workflow(
        input_data,
        options,
        free_energy=free_energy,
    )

    assert set(result.uncertainties) >= {"sigma_VT", "sigma_KT", "sigma_Kp"}
    for key in ("sigma_VT", "sigma_KT", "sigma_Kp"):
        values = result.uncertainties[key]
        assert values.shape == result.valid_mask.shape
        assert np.all(np.isfinite(values[result.valid_mask]))
        assert np.all(values[result.valid_mask] >= 0.0)
    assert np.any(result.uncertainties["sigma_VT"][result.valid_mask] > 0.0)


def test_eos_rejects_out_of_sample_volume() -> None:
    input_data, free_energy = _eos_input()
    options = _options(
        pressure_min=100.0,
        pressure_max=100.0,
        pressure_step=1.0,
        estimate_uncertainties=False,
        extrapolation_policy="fail",
        fit_failure_policy="continue",
    )

    result = run_volume_minimization_workflow(
        input_data,
        options,
        free_energy=free_energy,
    )

    assert not result.valid_mask[0, 0]
    assert len(result.failed_points) == 1
    assert result.failed_points[0].stage == "volume_extrapolation"


def test_analysis_eos_path_also_fits_once_per_temperature(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from quantas.modules.qha.analysis import analyze_volume_minimization

    input_data, free_energy = _eos_input(temperatures=2)
    options = _options(
        temperature_max=400.0,
        estimate_uncertainties=False,
    )
    calls = 0
    original = EnergyEOS.fit

    def counted_fit(self, *args, **kwargs):
        nonlocal calls
        calls += 1
        return original(self, *args, **kwargs)

    monkeypatch.setattr(EnergyEOS, "fit", counted_fit)

    result = analyze_volume_minimization(
        input_data,
        options,
        free_energy=free_energy,
    )

    assert calls == 2
    assert result.valid_mask.all()
    assert result.metadata["eos_workflow"]["fit_count"] == 2


@pytest.mark.parametrize(
    ("eos_tag", "parameters"),
    [
        ("BM2", np.array([-100.0, 160.0, 72.0], dtype=np.float64)),
        ("BM3", np.array([-100.0, 160.0, 4.2, 72.0], dtype=np.float64)),
        ("BM4", np.array([-100.0, 160.0, 4.2, -0.02, 72.0], dtype=np.float64)),
    ],
)
def test_qha_eos_workflow_preserves_selected_birch_murnaghan_order(
    eos_tag: str,
    parameters: np.ndarray,
) -> None:
    volume = np.linspace(68.0, 76.0, 17)
    native = parameters.copy()
    native[1] = pressure_to_energy(parameters[1], "eV", "A", "GPa")
    if eos_tag == "BM4":
        native[3] = parameters[3] / pressure_to_energy(1.0, "eV", "A", "GPa")
    energy = EnergyEOS().evaluate(eos_tag, volume, native)
    input_data = QHAInput(jobname=eos_tag, volume=volume, energy=energy)
    options = _options(eos=eos_tag, estimate_uncertainties=False)

    result = run_volume_minimization_workflow(
        input_data,
        options,
        free_energy=energy[np.newaxis, :],
    )

    assert result.completed
    assert result.valid_mask.all()
    assert result.metadata["eos_workflow"]["eos"] == eos_tag
    assert result.fit_records[0].fit.metadata["eos_tag"] == eos_tag
    np.testing.assert_allclose(result.equilibrium_volume[0, 0], 72.0, atol=1.0e-7)
    np.testing.assert_allclose(result.isothermal_bulk_modulus[0, 0], 160.0, rtol=1.0e-7)
    expected_kp = 4.0 if eos_tag == "BM2" else 4.2
    np.testing.assert_allclose(
        result.bulk_modulus_derivative[0, 0], expected_kp, rtol=1.0e-7
    )


def test_bm_alias_is_identical_to_bm3_in_qha() -> None:
    input_data, free_energy = _eos_input(noisy=True)
    common = dict(estimate_uncertainties=False)

    alias = run_volume_minimization_workflow(
        input_data,
        _options(eos="BM", **common),
        free_energy=free_energy,
    )
    explicit = run_volume_minimization_workflow(
        input_data,
        _options(eos="BM3", **common),
        free_energy=free_energy,
    )

    np.testing.assert_array_equal(alias.equilibrium_volume, explicit.equilibrium_volume)
    np.testing.assert_array_equal(
        alias.isothermal_bulk_modulus, explicit.isothermal_bulk_modulus
    )
    np.testing.assert_array_equal(
        alias.bulk_modulus_derivative, explicit.bulk_modulus_derivative
    )
