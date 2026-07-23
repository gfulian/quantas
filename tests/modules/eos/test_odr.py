"""P-V integration tests for the ODRPACK95 adapter."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.math.fitting import FitMethod, ODROptions, ParameterState
from quantas.core.physics.eos import (
    PressureEOS,
    available_eos_models,
    implied_kp,
    implied_kpp,
)
from quantas.modules.eos import (
    EOSDataset,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    ParameterConstraint,
    read_eos_input,
)

_DATA = Path(__file__).with_name("data")


@pytest.mark.parametrize(
    ("filename", "expected"),
    (
        (
            "PV_quartz.dat",
            {"K0": 37.12518, "KP": 5.98867, "V0": 112.980885},
        ),
        (
            "PV_topaz.dat",
            {"K0": 161.71745, "KP": 3.00184, "V0": 345.50950},
        ),
    ),
)
def test_bm3_odr_is_stable_on_experimental_pv_datasets(
    filename: str,
    expected: dict[str, float],
) -> None:
    dataset = read_eos_input(_DATA / filename)
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            options=EOSFitOptions(solver_options=ODROptions()),
        ),
    )

    assert result.fit.success
    assert result.fit.method is FitMethod.ODR
    assert result.fit.diagnostics is not None
    for name, value in expected.items():
        tolerance = 3.0e-3 if name != "V0" else 3.0e-4
        assert result.parameter_values[name] == pytest.approx(value, abs=tolerance)
    diagnostics = result.fit.diagnostics
    assert diagnostics.x_corrections is not None
    assert diagnostics.y_corrections is not None
    assert diagnostics.metadata["backend"] == "odrpack95"
    assert diagnostics.metadata["sum_square_x"] > 0.0
    assert diagnostics.metadata["sum_square_y"] > 0.0
    adjusted = np.asarray(result.fit.metadata["adjusted_x"], dtype=np.float64)
    assert adjusted.shape == result.fit.fitted.shape


def test_bm3_odr_respects_fixed_parameters_and_propagates_implied_covariance() -> None:
    dataset = read_eos_input(_DATA / "PV_quartz.dat")
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            constraints=(ParameterConstraint.fixed("KP", 6.0),),
            options=EOSFitOptions(solver_options=ODROptions()),
        ),
    )

    assert result.fit.success
    assert result.parameter_values["KP"] == 6.0
    assert result.fit.parameter_states == (
        ParameterState.FREE,
        ParameterState.FIXED,
        ParameterState.IMPLIED,
        ParameterState.FREE,
    )
    assert result.fit.covariance is not None
    assert result.fit.covariance[1, 1] == 0.0
    assert result.fit.covariance[2, 2] > 0.0


@pytest.mark.parametrize(
    "model",
    available_eos_models(),
    ids=lambda item: item.tag,
)
def test_odr_fits_every_supported_pressure_family_and_order(model) -> None:
    """Verify that the backend-neutral adapter is not specialized to BM3."""
    kp = implied_kp(model) if model.order == 2 else 4.5
    if kp is None:
        kp = 4.5
    kpp = implied_kpp(model, 150.0, kp) if model.order != 4 else -0.02
    physical = {"K0": 150.0, "KP": kp, "KPP": kpp, "V0": 100.0}
    volume = np.linspace(82.0, 100.0, 24, dtype=np.float64)
    pressure = PressureEOS().pressure(model, physical, volume)
    pressure = pressure + 5.0e-4 * np.sin(np.arange(volume.size))
    dataset = EOSDataset(
        columns={
            "pressure": pressure,
            "sigma_pressure": np.full(volume.shape, 0.01),
            "volume": volume,
            "sigma_volume": np.full(volume.shape, 0.02),
        }
    )

    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model=model,
            options=EOSFitOptions(solver_options=ODROptions()),
        ),
    )

    assert result.fit.success
    assert result.parameter_values["K0"] == pytest.approx(150.0, rel=1.0e-3)
    assert result.parameter_values["V0"] == pytest.approx(100.0, rel=1.0e-5)
