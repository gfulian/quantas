"""Characterization tests for directional elastic properties."""

from __future__ import annotations

import numpy as np

from quantas.core.physics.elasticity import (
    ElasticTensor,
    OrthorhombicElasticTensor,
    direction_vector,
    linear_compressibility,
    poisson_ratio,
    sample_linear_compressibility,
    sample_poisson_ratio,
    sample_shear_modulus,
    sample_young_modulus,
    shear_modulus,
    transverse_poisson_extrema,
    transverse_shear_extrema,
    young_modulus,
)


def test_direction_vectors_follow_quantas_angular_convention() -> None:
    np.testing.assert_allclose(direction_vector(0.0, 0.0), [0.0, 0.0, 1.0])
    np.testing.assert_allclose(
        direction_vector(np.pi / 2.0, 0.0),
        [1.0, 0.0, 0.0],
        atol=1.0e-15,
    )


def test_reference_directional_values(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    angles_2d = [1.123, 0.789]
    angles_3d = [1.123, 0.789, 0.456]
    np.testing.assert_allclose(
        young_modulus(hydroxylapatite_tensor, angles_2d), 126.65842872115124, rtol=1e-12
    )
    np.testing.assert_allclose(
        linear_compressibility(hydroxylapatite_tensor, angles_2d),
        3.0502230801977754,
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        shear_modulus(hydroxylapatite_tensor, angles_3d), 50.751953617499105, rtol=1e-12
    )
    np.testing.assert_allclose(
        poisson_ratio(hydroxylapatite_tensor, angles_3d),
        0.37466871471294666,
        rtol=1e-12,
    )


def test_orthorhombic_formulas_match_general_tensor(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    specialized = OrthorhombicElasticTensor(hydroxylapatite_tensor)
    for angles in ([1.123, 0.789], [np.pi / 2.0, 0.0]):
        np.testing.assert_allclose(
            young_modulus(specialized, angles),
            young_modulus(hydroxylapatite_tensor, angles),
            rtol=1e-10,
        )
        np.testing.assert_allclose(
            linear_compressibility(specialized, angles),
            linear_compressibility(hydroxylapatite_tensor, angles),
            rtol=1e-10,
        )


def test_reference_transverse_extrema(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    angles = [1.123, 0.789]
    np.testing.assert_allclose(
        transverse_shear_extrema(hydroxylapatite_tensor, angles),
        [49.74293181127166, 55.42544215914348],
        rtol=1e-11,
    )
    np.testing.assert_allclose(
        transverse_poisson_extrema(hydroxylapatite_tensor, angles),
        [0.0, 0.19601706948447487, 0.417646467928685],
        rtol=1e-11,
    )


def test_reference_xy_polar_values_and_progress(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    phi = np.linspace(0.0, 2.0 * np.pi, 5)
    theta = np.full(5, np.pi / 2.0)
    calls: list[tuple[int, int]] = []
    young = sample_young_modulus(
        hydroxylapatite_tensor,
        theta,
        phi,
        progress_callback=lambda current, total: calls.append((current, total)),
    )
    np.testing.assert_allclose(young, np.full(5, 147.9715908465), rtol=1e-10)
    np.testing.assert_allclose(
        sample_linear_compressibility(hydroxylapatite_tensor, theta, phi),
        np.tile([3.29642322, 0.0], (5, 1)),
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        sample_shear_modulus(hydroxylapatite_tensor, theta, phi),
        np.tile([39.687, 61.007], (5, 1)),
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        sample_poisson_ratio(hydroxylapatite_tensor, theta, phi),
        np.tile([0.0, 0.21273278569441567, 0.29949022650578877], (5, 1)),
        rtol=1e-10,
    )
    assert calls == [(1, 5), (2, 5), (3, 5), (4, 5), (5, 5)]
