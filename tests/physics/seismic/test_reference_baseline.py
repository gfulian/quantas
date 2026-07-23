"""Stored Quantas 0.9 SEISMIC directional regression baseline."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from tests.reference.seismic_reference import SeismicFormulaReference


BASELINE_DIR = Path(__file__).parents[2] / "baselines"


@pytest.mark.physics
@pytest.mark.seismic
@pytest.mark.baseline
def test_directional_reference_baseline(
    hydroxylapatite_reference: SeismicFormulaReference,
) -> None:
    """All primary reference arrays remain fixed at selected generic directions."""
    arrays_file = np.load(BASELINE_DIR / "seismic_reference.npz")
    expected = {key: arrays_file[key] for key in arrays_file.files}
    metadata = json.loads(
        (BASELINE_DIR / "seismic_reference.json").read_text(encoding="utf-8")
    )
    directions = expected["directions"]
    actual = [hydroxylapatite_reference.solve(direction) for direction in directions]
    np.testing.assert_allclose(
        hydroxylapatite_reference.christoffel_hessian,
        expected["christoffel_hessian"],
        **metadata["tolerance"],
    )
    fields = {
        "christoffel": np.stack([result.christoffel for result in actual]),
        "christoffel_gradient": np.stack(
            [result.christoffel_gradient for result in actual]
        ),
        "eigenvalues": np.stack([result.eigenvalues for result in actual]),
        "polarizations": np.stack([result.polarizations for result in actual]),
        "phase_speeds": np.stack([result.phase_speeds for result in actual]),
        "eigenvalue_gradients": np.stack(
            [result.eigenvalue_gradients for result in actual]
        ),
        "eigenvalue_hessians": np.stack(
            [result.eigenvalue_hessians for result in actual]
        ),
        "group_velocities": np.stack([result.group_velocities for result in actual]),
        "group_speeds": np.stack([result.group_speeds for result in actual]),
        "group_directions": np.stack([result.group_directions for result in actual]),
        "power_flow_angles": np.stack([result.power_flow_angles for result in actual]),
        "enhancement": np.stack([result.enhancement for result in actual]),
        "log10_enhancement": np.log10(
            np.stack([result.enhancement for result in actual])
        ),
    }
    tolerance = metadata["tolerance"]
    for name, values in fields.items():
        if name == "polarizations":
            overlap = np.abs(np.einsum("nmi,nmi->nm", values, expected[name]))
            np.testing.assert_allclose(overlap, np.ones_like(overlap), **tolerance)
        else:
            np.testing.assert_allclose(values, expected[name], **tolerance)


@pytest.mark.physics
@pytest.mark.seismic
@pytest.mark.baseline
def test_logarithmic_enhancement_baseline_is_base_ten() -> None:
    """The persisted logarithmic field is exactly log10 of raw enhancement."""
    arrays_file = np.load(BASELINE_DIR / "seismic_reference.npz")
    np.testing.assert_allclose(
        arrays_file["log10_enhancement"],
        np.log10(arrays_file["enhancement"]),
        rtol=0.0,
        atol=np.finfo(float).eps,
    )
    np.testing.assert_allclose(
        10.0 ** arrays_file["log10_enhancement"],
        arrays_file["enhancement"],
        rtol=1.0e-14,
        atol=1.0e-14,
    )
