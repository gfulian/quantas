"""Scientific regression checks for the Quantas 2.0 baseline."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from quantas.api import elasticity, ha, qha

PROJECT_ROOT = Path(__file__).parents[2]
BASELINE_DIR = PROJECT_ROOT / "tests" / "baselines"


@pytest.fixture(scope="module")
def baseline() -> tuple[dict[str, np.ndarray], dict]:
    """Load baseline arrays and comparison tolerances."""
    arrays_file = np.load(BASELINE_DIR / "scientific_reference.npz")
    arrays = {key: arrays_file[key] for key in arrays_file.files}
    metadata = json.loads(
        (BASELINE_DIR / "scientific_reference.json").read_text(encoding="utf-8")
    )
    return arrays, metadata


def _qha_tolerance(metadata: dict, property_name: str) -> dict[str, float]:
    """Return the baseline tolerance for one QHA property.

    Primary state variables use the strict default QHA tolerance. Quantities
    obtained from fitted derivatives use a slightly wider, but still stringent,
    cross-platform tolerance because nonlinear solvers and BLAS/LAPACK builds
    may differ in their final few digits.

    Parameters
    ----------
    metadata
        Baseline metadata loaded from the JSON manifest.
    property_name
        Name of the QHA result property being compared.

    Returns
    -------
    dict
        Keyword arguments suitable for :func:`numpy.testing.assert_allclose`.
    """
    tolerances = metadata["tolerances"]
    return tolerances.get("qha_derived", {}).get(
        property_name,
        tolerances["qha"],
    )


@pytest.mark.baseline
def test_ha_scientific_baseline(baseline) -> None:
    """HA reference arrays remain unchanged within strict numerical tolerance."""
    expected, metadata = baseline
    envelope = ha.run(
        ha.read_input(PROJECT_ROOT / "tests/modules/ha/data/mgo_b3lyp_qha.yaml"),
        options=ha.Options(
            temperature_min=0.0,
            temperature_max=200.0,
            temperature_step=100.0,
        ),
    )
    result = ha.get_result(envelope)
    tolerance = metadata["tolerances"]["ha"]
    actual = {
        "ha_temperature": result.temperature,
        "ha_volume": result.volume,
        "ha_zero_point_energy": result.zero_point_energy,
        "ha_free_energy": result.free_energy,
        "ha_isochoric_heat_capacity": result.isochoric_heat_capacity,
    }
    for key, value in actual.items():
        np.testing.assert_allclose(value, expected[key], **tolerance)


@pytest.mark.baseline
def test_qha_four_workflow_scientific_baseline(baseline) -> None:
    """All four QHA workflow combinations preserve their reference arrays."""
    expected, metadata = baseline
    input_data = qha.read_input(
        PROJECT_ROOT / "tests/modules/qha/data/mgo_b3lyp_qha.yaml"
    )
    for scheme in ("freq", "td"):
        for minimization in ("poly", "eos"):
            name = f"{scheme}_{minimization}"
            envelope = qha.run(
                input_data,
                options=qha.Options(
                    temperature_min=0.0,
                    temperature_max=200.0,
                    temperature_step=200.0,
                    pressure_min=0.0,
                    pressure_max=5.0,
                    pressure_step=5.0,
                    scheme=scheme,
                    minimization=minimization,
                    eos="BM",
                    energy_degree=3,
                    free_energy_degree=3,
                    frequency_degree=3,
                    polynomial_derivative_method="local_grid",
                    estimate_uncertainties=minimization == "eos",
                    uncertainty_method="covariance",
                ),
            )
            result = qha.get_result(envelope)
            actual = {
                "equilibrium_volume": result.equilibrium_volume,
                "isothermal_bulk_modulus": result.isothermal_bulk_modulus,
                "free_energy": result.free_energy,
                "isobaric_heat_capacity": result.isobaric_heat_capacity,
            }
            for property_name, value in actual.items():
                key = f"qha_{name}_{property_name}"
                np.testing.assert_allclose(
                    value,
                    expected[key],
                    err_msg=f"QHA baseline mismatch for {name}.{property_name}",
                    **_qha_tolerance(metadata, property_name),
                )
            np.testing.assert_array_equal(
                result.valid_mask,
                expected[f"qha_{name}_valid_mask"],
            )


@pytest.mark.baseline
def test_elasticity_scientific_baseline(baseline) -> None:
    """Elasticity tensors, averages, stability and polar curves remain stable."""
    expected, metadata = baseline
    envelope = elasticity.run(
        elasticity.read_input(
            PROJECT_ROOT / "tests/modules/elasticity/data/hydroxylapatite.dat"
        ),
        options=elasticity.Options(calculate_2d=True, ntheta=9),
    )
    result = elasticity.get_result(envelope)
    tolerance = metadata["tolerances"]["elasticity"]
    actual = {
        "elasticity_stiffness": result.stiffness,
        "elasticity_compliance": result.compliance,
        "elasticity_averages": result.averages.as_array(),
        "elasticity_eigenvalues": result.stability.eigenvalues,
        "elasticity_xy_young_modulus": result.properties_2d["xy"]["young_modulus"],
        "elasticity_xy_shear_modulus": result.properties_2d["xy"]["shear_modulus"],
        "elasticity_xy_poisson_ratio": result.properties_2d["xy"]["poisson_ratio"],
    }
    for key, value in actual.items():
        np.testing.assert_allclose(value, expected[key], **tolerance)
