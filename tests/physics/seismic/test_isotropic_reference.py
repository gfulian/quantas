"""Elastic-medium and isotropic-reference characterization."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor, compute_elastic_averages
from quantas.core.physics.seismic import (
    ElasticMedium,
    IsotropicSeismicVelocities,
    isotropic_seismic_velocities,
)


@pytest.mark.physics
@pytest.mark.seismic
def test_isotropic_velocities_match_reference_values(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    """Hill-average seismic references remain numerically unchanged."""
    stiffness, density = hydroxylapatite_data
    medium = ElasticMedium(ElasticTensor(stiffness), density)
    averages = compute_elastic_averages(medium.elastic_tensor)
    velocities = isotropic_seismic_velocities(averages, medium.density)

    assert isinstance(velocities, IsotropicSeismicVelocities)
    np.testing.assert_allclose(
        velocities.as_array(),
        [4.00911178, 4.00911178, 7.64303712],
        rtol=1.0e-8,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_isotropic_reference_velocities_follow_density_scaling(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    """Both isotropic wave speeds scale as the inverse square root of density."""
    stiffness, density = hydroxylapatite_data
    averages = compute_elastic_averages(ElasticTensor(stiffness))
    reference = isotropic_seismic_velocities(averages, density)
    scaled = isotropic_seismic_velocities(averages, 4.0 * density)

    np.testing.assert_allclose(scaled.as_array(), reference.as_array() / 2.0)


@pytest.mark.physics
@pytest.mark.seismic
@pytest.mark.parametrize("density", [0.0, -1.0, np.nan, np.inf])
def test_isotropic_reference_rejects_invalid_density(
    hydroxylapatite_data: tuple[np.ndarray, float],
    density: float,
) -> None:
    """Isotropic references require a finite positive density."""
    stiffness, _ = hydroxylapatite_data
    averages = compute_elastic_averages(ElasticTensor(stiffness))
    with pytest.raises(ValueError, match="finite positive density"):
        isotropic_seismic_velocities(averages, density)
