"""Units and validation of the production Christoffel solver."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ChristoffelSolver, ElasticMedium


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    return stiffness


@pytest.mark.physics
@pytest.mark.seismic
def test_reduced_stiffness_and_phase_speeds_match_full_si_units() -> None:
    stiffness = _isotropic_stiffness(lame=70.0, shear=50.0)
    density = 3000.0
    tensor = ElasticTensor(stiffness)
    medium = ElasticMedium(tensor, density)

    for value in (-1.0, 1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match="pseudoinverse_rcond"):
            ChristoffelSolver(medium, pseudoinverse_rcond=value)

    solver = ChristoffelSolver(medium)
    direction = np.asarray([1.0, 2.0, 3.0], dtype=float)
    direction /= np.linalg.norm(direction)

    np.testing.assert_allclose(
        solver.reduced_stiffness,
        tensor.stiffness_tensor * 1000.0 / density,
        rtol=1.0e-15,
        atol=1.0e-15,
    )
    result = solver.solve_direction(direction)

    stiffness_pa = tensor.stiffness_tensor * 1.0e9
    gamma_si = (
        np.einsum(
            "j,ijkl,l->ik",
            direction,
            stiffness_pa,
            direction,
            optimize=True,
        )
        / density
    )
    expected_km_s = np.sqrt(np.linalg.eigvalsh(gamma_si)) / 1000.0
    np.testing.assert_allclose(
        result.phase_speeds,
        expected_km_s,
        rtol=1.0e-13,
        atol=1.0e-13,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_phase_speeds_follow_inverse_square_root_density_scaling(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    direction = [1.0, 2.0, 3.0]
    first = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), density)
    ).solve_direction(direction)
    second = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), 4.0 * density)
    ).solve_direction(direction)

    np.testing.assert_allclose(second.eigenvalues, first.eigenvalues / 4.0)
    np.testing.assert_allclose(second.phase_speeds, first.phase_speeds / 2.0)
    np.testing.assert_allclose(second.polarizations, first.polarizations)
    np.testing.assert_array_equal(second.valid_mask, first.valid_mask)
    np.testing.assert_array_equal(second.degeneracy_mask, first.degeneracy_mask)


@pytest.mark.physics
@pytest.mark.seismic
def test_elastic_medium_requires_a_tensor_and_positive_density() -> None:
    stiffness = np.eye(6)
    tensor = ElasticTensor(stiffness)

    with pytest.raises(TypeError, match="ElasticTensor"):
        ElasticMedium(stiffness, 3000.0)  # type: ignore[arg-type]

    for density in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match="finite and positive"):
            ElasticMedium(tensor, density)

    medium = ElasticMedium(tensor, 3000.0)
    assert medium.elastic_tensor is tensor
    assert medium.density == pytest.approx(3000.0)
    assert ChristoffelSolver(medium).medium is medium
    assert ChristoffelSolver(medium).density == pytest.approx(3000.0)

    with pytest.raises(TypeError, match="ElasticMedium"):
        ChristoffelSolver(tensor)  # type: ignore[arg-type]


@pytest.mark.physics
@pytest.mark.seismic
def test_solver_rejects_invalid_tensor_tolerances_and_directions() -> None:
    with pytest.raises(TypeError, match="ElasticMedium"):
        ChristoffelSolver(np.eye(6))  # type: ignore[arg-type]

    medium = ElasticMedium(ElasticTensor(np.eye(6)), 3000.0)
    for name in (
        "eigenvalue_rtol",
        "eigenvalue_atol",
        "degeneracy_rtol",
        "degeneracy_atol",
        "caustic_rtol",
        "caustic_atol",
    ):
        for value in (-1.0, np.nan, np.inf):
            kwargs = {name: value}
            with pytest.raises(ValueError, match=name):
                ChristoffelSolver(medium, **kwargs)

    for value in (-1.0, 1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match="pseudoinverse_rcond"):
            ChristoffelSolver(medium, pseudoinverse_rcond=value)

    solver = ChristoffelSolver(medium)
    invalid_directions = (
        [1.0, 2.0],
        [1.0, 2.0, 3.0, 4.0],
        [0.0, 0.0, 0.0],
        [1.0, np.nan, 0.0],
        [1.0, np.inf, 0.0],
    )
    for direction in invalid_directions:
        with pytest.raises(ValueError):
            solver.solve_direction(direction)
