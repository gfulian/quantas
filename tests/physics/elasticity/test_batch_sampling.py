"""Tests for exact batched directional elasticity sampling."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    ElasticTensor,
    exact_transverse_extrema,
    linear_compressibility,
    poisson_ratio,
    sample_elastic_directional_field,
    shear_modulus,
    young_modulus,
)


def _isotropic_tensor() -> ElasticTensor:
    """Return a stable isotropic elastic tensor."""
    lam = 50.0
    mu = 30.0
    stiffness = np.array(
        [
            [lam + 2.0 * mu, lam, lam, 0.0, 0.0, 0.0],
            [lam, lam + 2.0 * mu, lam, 0.0, 0.0, 0.0],
            [lam, lam, lam + 2.0 * mu, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, mu, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, mu, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, mu],
        ],
        dtype=float,
    )
    return ElasticTensor(stiffness)


@pytest.mark.physics
@pytest.mark.elasticity
def test_batch_longitudinal_fields_match_pointwise_operations(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """Batched Young and compressibility fields reproduce scalar formulas."""
    theta = np.array([0.0, 0.37, 1.123, 2.41, np.pi])
    phi = np.array([0.0, 1.91, 0.789, 5.77, 3.12])
    field = sample_elastic_directional_field(
        hydroxylapatite_tensor,
        theta,
        phi,
        properties=("young", "compressibility"),
    )
    expected_young = np.array(
        [
            young_modulus(hydroxylapatite_tensor, (polar, azimuth))
            for polar, azimuth in zip(theta, phi, strict=True)
        ]
    )
    expected_compressibility = np.array(
        [
            linear_compressibility(
                hydroxylapatite_tensor,
                (polar, azimuth),
            )
            for polar, azimuth in zip(theta, phi, strict=True)
        ]
    )
    np.testing.assert_allclose(field.young_modulus, expected_young, rtol=2.0e-13)
    np.testing.assert_allclose(
        field.linear_compressibility,
        expected_compressibility,
        rtol=2.0e-13,
        atol=2.0e-13,
    )


@pytest.mark.physics
@pytest.mark.elasticity
def test_transverse_extrema_match_returned_angles(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """The algebraic extrema are verified against independent scalar formulas."""
    directions = ((0.31, 0.12), (1.123, 0.789), (2.47, 5.13))
    for theta, phi in directions:
        shear = exact_transverse_extrema(
            hydroxylapatite_tensor,
            theta,
            phi,
            kind="shear",
        )
        poisson = exact_transverse_extrema(
            hydroxylapatite_tensor,
            theta,
            phi,
            kind="poisson",
        )
        assert shear_modulus(
            hydroxylapatite_tensor,
            (theta, phi, shear.minimum_angle),
        ) == pytest.approx(shear.minimum, rel=2.0e-13, abs=2.0e-13)
        assert shear_modulus(
            hydroxylapatite_tensor,
            (theta, phi, shear.maximum_angle),
        ) == pytest.approx(shear.maximum, rel=2.0e-13, abs=2.0e-13)
        assert poisson_ratio(
            hydroxylapatite_tensor,
            (theta, phi, poisson.minimum_angle),
        ) == pytest.approx(poisson.minimum, rel=2.0e-13, abs=2.0e-13)
        assert poisson_ratio(
            hydroxylapatite_tensor,
            (theta, phi, poisson.maximum_angle),
        ) == pytest.approx(poisson.maximum, rel=2.0e-13, abs=2.0e-13)


@pytest.mark.physics
@pytest.mark.elasticity
def test_batch_transverse_fields_match_exact_pointwise_solution(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """Batched transverse values reproduce exact pointwise eigensolutions."""
    theta = np.array([[0.23, 0.91, 1.71], [2.37, 2.82, np.pi]])
    phi = np.array([[0.0, 0.73, 2.11], [3.89, 5.77, 1.2]])
    field = sample_elastic_directional_field(
        hydroxylapatite_tensor,
        theta,
        phi,
        properties=("shear", "poisson"),
    )
    expected = {
        name: np.empty(theta.shape, dtype=float)
        for name in (
            "shear_minimum",
            "shear_maximum",
            "poisson_minimum",
            "poisson_maximum",
        )
    }
    for index in np.ndindex(theta.shape):
        shear = exact_transverse_extrema(
            hydroxylapatite_tensor,
            float(theta[index]),
            float(phi[index]),
            kind="shear",
        )
        poisson = exact_transverse_extrema(
            hydroxylapatite_tensor,
            float(theta[index]),
            float(phi[index]),
            kind="poisson",
        )
        expected["shear_minimum"][index] = shear.minimum
        expected["shear_maximum"][index] = shear.maximum
        expected["poisson_minimum"][index] = poisson.minimum
        expected["poisson_maximum"][index] = poisson.maximum

    for name, values in expected.items():
        np.testing.assert_allclose(
            getattr(field, name),
            values,
            rtol=3.0e-13,
            atol=3.0e-13,
        )


@pytest.mark.physics
@pytest.mark.elasticity
def test_isotropic_degeneracy_is_constant_without_branch_spikes() -> None:
    """Exact transverse eigenvalues stay constant at complete degeneracy."""
    theta = np.linspace(0.0, np.pi, 41)[:, np.newaxis]
    phi = np.linspace(0.0, 2.0 * np.pi, 83, endpoint=False)[np.newaxis, :]
    field = sample_elastic_directional_field(_isotropic_tensor(), theta, phi)

    np.testing.assert_allclose(field.young_modulus, 78.75, rtol=2.0e-14)
    np.testing.assert_allclose(field.linear_compressibility, 1000.0 / 210.0)
    np.testing.assert_allclose(field.shear_minimum, 30.0, rtol=2.0e-14)
    np.testing.assert_allclose(field.shear_maximum, 30.0, rtol=2.0e-14)
    np.testing.assert_allclose(field.poisson_minimum, 0.3125, rtol=2.0e-14)
    np.testing.assert_allclose(field.poisson_maximum, 0.3125, rtol=2.0e-14)
    assert field.diagnostics.shear_degeneracy_count == theta.size * phi.size
    assert field.diagnostics.poisson_degeneracy_count == theta.size * phi.size


@pytest.mark.physics
@pytest.mark.elasticity
def test_closed_plane_field_is_periodic_and_locally_smooth(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """Closed polar fields contain no isolated optimizer-induced outliers."""
    phi = np.linspace(0.0, 2.0 * np.pi, 721, endpoint=True)
    theta = np.full_like(phi, 1.017)
    field = sample_elastic_directional_field(
        hydroxylapatite_tensor,
        theta,
        phi,
    )

    for name in (
        "young_modulus",
        "linear_compressibility",
        "shear_minimum",
        "shear_maximum",
        "poisson_minimum",
        "poisson_maximum",
    ):
        values = getattr(field, name)
        assert values is not None
        assert values[0] == pytest.approx(values[-1], rel=2.0e-13, abs=2.0e-13)

        # A smooth tensor contraction has a locally second-order residual on a
        # uniform angular grid.  A single failed local optimizer would exceed
        # this scale by orders of magnitude.
        residual = np.abs(values[1:-1] - 0.5 * (values[:-2] + values[2:]))
        scale = max(float(np.ptp(values)), float(np.max(np.abs(values))), 1.0)
        assert float(np.max(residual)) < 2.0e-4 * scale


@pytest.mark.physics
@pytest.mark.elasticity
def test_field_arrays_and_progress_contract(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """The neutral field is immutable and reports monotonic batched progress."""
    calls: list[tuple[int, int]] = []
    field = sample_elastic_directional_field(
        hydroxylapatite_tensor,
        np.array([0.1, 0.2, 0.3]),
        np.array([0.4, 0.5, 0.6]),
        properties=("young", "shear", "poisson"),
        batch_size=1,
        progress_callback=lambda current, total: calls.append((current, total)),
    )
    assert calls == [(3, 9), (6, 9), (9, 9)]
    assert field.batch_size == 1
    assert field.young_modulus is not None
    with pytest.raises(ValueError, match="read-only"):
        field.young_modulus[0] = 0.0
