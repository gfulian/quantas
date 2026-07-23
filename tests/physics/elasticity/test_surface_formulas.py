"""Frozen anisotropic surface baseline for a future batched sampler."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    ElasticTensor,
    linear_compressibility,
    sample_elasticity_surfaces,
    transverse_extrema,
    young_modulus,
)

_BASELINE_DIR = Path(__file__).parents[2] / "baselines"


def _signed_field(
    surfaces: dict[str, object],
    positive_key: str,
    negative_key: str,
    shape: tuple[int, int],
) -> np.ndarray:
    """Reconstruct a signed field from positive and negative plot branches."""
    result = np.zeros(shape, dtype=float)
    for key in (positive_key, negative_key):
        surface = surfaces.get(key)
        if surface is not None:
            result += np.nan_to_num(surface.values, nan=0.0)  # type: ignore[attr-defined]
    return result


def _poisson_minimum_field(
    surfaces: dict[str, object],
    shape: tuple[int, int],
) -> np.ndarray:
    """Reconstruct the signed minimum Poisson field from plot branches."""
    result = np.zeros(shape, dtype=float)
    for key in ("poisson_negative", "poisson_minimum"):
        surface = surfaces.get(key)
        if surface is None:
            continue
        values = surface.values  # type: ignore[attr-defined]
        finite = np.isfinite(values)
        result[finite] = values[finite]
    return result


@pytest.mark.physics
@pytest.mark.elasticity
@pytest.mark.baseline
def test_anisotropic_surfaces_match_reference(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """All signed physical fields remain fixed before batch optimization."""
    arrays_file = np.load(_BASELINE_DIR / "elasticity_reference.npz")
    expected = {key: arrays_file[key] for key in arrays_file.files}
    metadata = json.loads(
        (_BASELINE_DIR / "elasticity_reference.json").read_text(encoding="utf-8")
    )
    grid = metadata["grid"]
    sampled = sample_elasticity_surfaces(
        hydroxylapatite_tensor,
        ntheta=grid["ntheta"],
        nphi=grid["nphi"],
    )
    surfaces = sampled.surfaces
    shape = surfaces["young"].values.shape
    actual = {
        "theta": surfaces["young"].theta,
        "phi": surfaces["young"].phi,
        "young": surfaces["young"].values,
        "compressibility": _signed_field(
            surfaces,
            "compressibility_positive",
            "compressibility_negative",
            shape,
        ),
        "shear_minimum": surfaces["shear_minimum"].values,
        "shear_maximum": surfaces["shear_maximum"].values,
        "poisson_minimum": _poisson_minimum_field(surfaces, shape),
        "poisson_maximum": surfaces["poisson_maximum"].values,
    }
    tolerance = metadata["tolerance"]
    for name, values in actual.items():
        np.testing.assert_allclose(values, expected[name], **tolerance)


@pytest.mark.physics
@pytest.mark.elasticity
def test_surface_sampler_matches_pointwise_formulas(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """The sampled surface is constrained by pointwise physical operations."""
    sampled = sample_elasticity_surfaces(
        hydroxylapatite_tensor,
        ntheta=4,
        nphi=7,
    )
    surfaces = sampled.surfaces
    shape = surfaces["young"].values.shape
    compressibility = _signed_field(
        surfaces,
        "compressibility_positive",
        "compressibility_negative",
        shape,
    )
    poisson_minimum = _poisson_minimum_field(surfaces, shape)
    indexes = ((0, 0), (1, 1), (1, 5), (2, 3), (3, 6))

    for row, column in indexes:
        theta = float(surfaces["young"].theta[row, column])
        phi = float(surfaces["young"].phi[row, column])
        assert surfaces["young"].values[row, column] == pytest.approx(
            young_modulus(hydroxylapatite_tensor, (theta, phi)),
            rel=2.0e-13,
            abs=2.0e-13,
        )
        assert compressibility[row, column] == pytest.approx(
            linear_compressibility(hydroxylapatite_tensor, (theta, phi)),
            rel=2.0e-13,
            abs=2.0e-13,
        )
        shear = transverse_extrema(
            hydroxylapatite_tensor,
            theta,
            phi,
            kind="shear",
        )
        poisson = transverse_extrema(
            hydroxylapatite_tensor,
            theta,
            phi,
            kind="poisson",
        )
        assert surfaces["shear_minimum"].values[row, column] == pytest.approx(
            shear.minimum,
            rel=2.0e-12,
            abs=2.0e-12,
        )
        assert surfaces["shear_maximum"].values[row, column] == pytest.approx(
            shear.maximum,
            rel=2.0e-12,
            abs=2.0e-12,
        )
        assert poisson_minimum[row, column] == pytest.approx(
            poisson.minimum,
            rel=2.0e-12,
            abs=2.0e-12,
        )
        assert surfaces["poisson_maximum"].values[row, column] == pytest.approx(
            poisson.maximum,
            rel=2.0e-12,
            abs=2.0e-12,
        )


@pytest.mark.physics
@pytest.mark.elasticity
def test_surface_progress_is_monotonic_and_finishes_at_total(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """Progress reports completed grid-property work in a stable unit."""
    calls: list[tuple[int, int]] = []
    sample_elasticity_surfaces(
        hydroxylapatite_tensor,
        ntheta=3,
        nphi=5,
        progress_callback=lambda current, total: calls.append((current, total)),
    )
    assert calls
    totals = {total for _, total in calls}
    assert totals == {60}
    completed = [current for current, _ in calls]
    assert completed == sorted(completed)
    assert completed[-1] == 60
