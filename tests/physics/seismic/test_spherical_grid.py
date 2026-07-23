"""Shared spherical grids used by acoustic field sampling."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.geometry import (
    Hemisphere,
    cartesian_to_spherical,
    close_periodic_seam,
    create_spherical_grid,
    spherical_direction,
)


@pytest.mark.physics
@pytest.mark.seismic
def test_spherical_grid_domains_and_directions_are_explicit() -> None:
    upper = create_spherical_grid(5, 8, hemisphere=Hemisphere.UPPER)
    lower = create_spherical_grid(5, 8, hemisphere="lower")
    full = create_spherical_grid(7, 8)

    assert upper.hemisphere is Hemisphere.UPPER
    assert lower.hemisphere is Hemisphere.LOWER
    assert full.hemisphere is Hemisphere.FULL
    assert upper.shape == (5, 8)
    assert upper.size == 40
    assert upper.flat_directions.shape == (40, 3)
    assert upper.theta[0] == pytest.approx(0.0)
    assert upper.theta[-1] == pytest.approx(np.pi / 2.0)
    assert lower.theta[0] == pytest.approx(np.pi / 2.0)
    assert lower.theta[-1] == pytest.approx(np.pi)
    assert full.theta[-1] == pytest.approx(np.pi)
    assert upper.phi[0] == pytest.approx(0.0)
    assert upper.phi[-1] < 2.0 * np.pi
    np.testing.assert_allclose(
        np.linalg.norm(upper.directions, axis=2),
        1.0,
        atol=1.0e-15,
    )
    np.testing.assert_array_equal(full.directions[0, :, :2], 0.0)
    np.testing.assert_array_equal(full.directions[-1, :, :2], 0.0)
    np.testing.assert_array_equal(upper.directions[-1, :, 2], 0.0)
    assert all(
        not array.flags.writeable
        for array in (
            upper.theta,
            upper.phi,
            upper.theta_grid,
            upper.phi_grid,
            upper.directions,
            upper.flat_directions,
        )
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_spherical_cartesian_round_trip_uses_standard_azimuth() -> None:
    theta = np.asarray([0.2, 0.7, 1.4, 2.5])
    phi = np.asarray([0.1, 1.2, 3.4, 5.8])
    directions = spherical_direction(theta, phi)
    recovered_theta, recovered_phi = cartesian_to_spherical(directions)

    np.testing.assert_allclose(recovered_theta, theta, atol=1.0e-14)
    np.testing.assert_allclose(recovered_phi, phi, atol=1.0e-14)


@pytest.mark.physics
@pytest.mark.seismic
def test_periodic_seam_is_added_only_for_presentation() -> None:
    values = np.arange(12.0).reshape(3, 4)
    closed = close_periodic_seam(values)

    assert closed.shape == (3, 5)
    np.testing.assert_array_equal(closed[:, :-1], values)
    np.testing.assert_array_equal(closed[:, -1], values[:, 0])
    assert not closed.flags.writeable


@pytest.mark.physics
@pytest.mark.seismic
@pytest.mark.parametrize(
    ("ntheta", "nphi"),
    [(1, 8), (5, 2), (True, 8), (5, False)],
)
def test_spherical_grid_rejects_invalid_dimensions(ntheta: int, nphi: int) -> None:
    with pytest.raises(ValueError):
        create_spherical_grid(ntheta, nphi)
