#!/usr/bin/env python3
"""Generate the Quantas 2.0 scientific regression baseline."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy

from quantas import __version__
from quantas.api import elasticity, ha, qha

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BASELINE_DIR = PROJECT_ROOT / "tests" / "baselines"
BASELINE_ARRAYS = BASELINE_DIR / "scientific_reference.npz"
BASELINE_METADATA = BASELINE_DIR / "scientific_reference.json"


def build_baseline() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Calculate representative HA, QHA, and elasticity reference arrays.

    Returns
    -------
    tuple
        Numerical arrays and serializable baseline metadata.
    """
    ha_input = ha.read_input(
        PROJECT_ROOT / "tests/modules/ha/data/mgo_b3lyp_qha.yaml"
    )
    ha_options = ha.Options(
        temperature_min=0.0,
        temperature_max=200.0,
        temperature_step=100.0,
    )
    ha_result = ha.run(ha_input, options=ha_options)
    ha_payload = ha.get_result(ha_result)

    arrays: dict[str, np.ndarray] = {
        "ha_temperature": np.asarray(ha_payload.temperature),
        "ha_volume": np.asarray(ha_payload.volume),
        "ha_zero_point_energy": np.asarray(ha_payload.zero_point_energy),
        "ha_free_energy": np.asarray(ha_payload.free_energy),
        "ha_isochoric_heat_capacity": np.asarray(ha_payload.isochoric_heat_capacity),
    }

    qha_input = qha.read_input(
        PROJECT_ROOT / "tests/modules/qha/data/mgo_b3lyp_qha.yaml"
    )
    qha_options: dict[str, dict[str, Any]] = {}
    for scheme in ("freq", "td"):
        for minimization in ("poly", "eos"):
            name = f"{scheme}_{minimization}"
            options = qha.Options(
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
            )
            qha_result = qha.run(qha_input, options=options)
            qha_payload = qha.get_result(qha_result)
            arrays[f"qha_{name}_equilibrium_volume"] = np.asarray(
                qha_payload.equilibrium_volume
            )
            arrays[f"qha_{name}_isothermal_bulk_modulus"] = np.asarray(
                qha_payload.isothermal_bulk_modulus
            )
            arrays[f"qha_{name}_free_energy"] = np.asarray(qha_payload.free_energy)
            arrays[f"qha_{name}_isobaric_heat_capacity"] = np.asarray(
                qha_payload.isobaric_heat_capacity
            )
            arrays[f"qha_{name}_valid_mask"] = np.asarray(qha_payload.valid_mask)
            qha_options[name] = {
                "scheme": scheme,
                "minimization": minimization,
                "eos": "BM",
                "temperature": [0.0, 200.0, 200.0],
                "pressure": [0.0, 5.0, 5.0],
            }

    elasticity_input = elasticity.read_input(
        PROJECT_ROOT / "tests/modules/elasticity/data/hydroxylapatite.dat"
    )
    elasticity_options = elasticity.Options(
        calculate_2d=True,
        ntheta=9,
    )
    elasticity_result = elasticity.run(
        elasticity_input,
        options=elasticity_options,
    )
    elasticity_payload = elasticity.get_result(elasticity_result)
    arrays.update(
        {
            "elasticity_stiffness": np.asarray(elasticity_payload.stiffness),
            "elasticity_compliance": np.asarray(elasticity_payload.compliance),
            "elasticity_averages": np.asarray(elasticity_payload.averages.as_array()),
            "elasticity_eigenvalues": np.asarray(elasticity_payload.stability.eigenvalues),
            "elasticity_xy_young_modulus": np.asarray(
                elasticity_payload.properties_2d["xy"]["young_modulus"]
            ),
            "elasticity_xy_shear_modulus": np.asarray(
                elasticity_payload.properties_2d["xy"]["shear_modulus"]
            ),
            "elasticity_xy_poisson_ratio": np.asarray(
                elasticity_payload.properties_2d["xy"]["poisson_ratio"]
            ),
        }
    )

    metadata = {
        "baseline": "Quantas 2.0 scientific checkpoint",
        "quantas_version": __version__,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "tolerances": {
            "ha": {"rtol": 1.0e-12, "atol": 1.0e-12},
            "qha": {"rtol": 1.0e-7, "atol": 1.0e-8},
            "qha_derived": {
                "isothermal_bulk_modulus": {"rtol": 1.0e-6, "atol": 1.0e-6},
                "isobaric_heat_capacity": {"rtol": 1.0e-6, "atol": 1.0e-6},
            },
            "elasticity": {"rtol": 1.0e-10, "atol": 1.0e-10},
        },
        "ha_options": {
            "temperature_min": 0.0,
            "temperature_max": 200.0,
            "temperature_step": 100.0,
        },
        "qha_options": qha_options,
        "elasticity_options": {
            "calculate_2d": True,
            "ntheta": 9,
        },
        "arrays": sorted(arrays),
    }
    return arrays, metadata


def main() -> None:
    """Generate the compressed baseline arrays and metadata manifest."""
    BASELINE_DIR.mkdir(parents=True, exist_ok=True)
    arrays, metadata = build_baseline()
    np.savez_compressed(BASELINE_ARRAYS, **arrays)
    metadata["npz_sha256"] = hashlib.sha256(BASELINE_ARRAYS.read_bytes()).hexdigest()
    BASELINE_METADATA.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(BASELINE_ARRAYS)
    print(BASELINE_METADATA)


if __name__ == "__main__":
    main()
