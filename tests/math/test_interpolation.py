"""Tests for the shared rectilinear field interpolator."""

from __future__ import annotations

import numpy as np

from quantas.core.numerics import RectilinearFieldInterpolator


def test_rectilinear_interpolator_handles_nonuniform_axes_and_trailing_fields() -> None:
    """Bilinear fields and trailing tensor dimensions share one exact engine."""
    x = np.asarray([0.0, 1.0, 3.0])
    y = np.asarray([0.0, 2.0, 5.0])
    xx, yy = np.meshgrid(x, y, indexing="ij")
    scalar = 2.0 * xx + 3.0 * yy
    field = np.stack((scalar, -scalar), axis=-1)
    interpolator = RectilinearFieldInterpolator(x, y, field)
    values, extrapolated = interpolator.evaluate_points([0.5, 2.0], [1.0, 4.0])
    np.testing.assert_allclose(values, [[4.0, -4.0], [16.0, -16.0]])
    assert not np.any(extrapolated)


def test_rectilinear_interpolator_handles_singletons_and_marks_extrapolation() -> None:
    """Singleton axes remain deterministic and extrapolation is explicit."""
    interpolator = RectilinearFieldInterpolator(
        [300.0],
        [0.0, 5.0],
        np.asarray([[10.0, 20.0]]),
    )
    values, extrapolated = interpolator.evaluate_grid([300.0, 400.0], [-1.0, 2.5])
    np.testing.assert_allclose(values, [[8.0, 15.0], [8.0, 15.0]])
    assert extrapolated.tolist() == [[True, False], [True, True]]
