"""Reference tests for elastic averages and mechanical stability."""

from __future__ import annotations

import numpy as np

from quantas.core.physics.elasticity import (
    ElasticTensor,
    check_positive_definiteness,
    compute_elastic_averages,
)


def test_vrh_averages_match_reference_values(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    expected = np.array(
        [
            [118.4746666667, 136.6398917700, 52.2412000000, 0.3077790304],
            [116.6044090776, 131.0543620594, 49.9186434654, 0.3126795217],
            [117.5395378721, 133.8503581846, 51.0799217327, 0.3102052004],
        ]
    )
    averages = compute_elastic_averages(hydroxylapatite_tensor)
    np.testing.assert_allclose(
        averages.as_array(), expected, rtol=1.0e-10, atol=1.0e-10
    )


def test_positive_definiteness_returns_structured_result(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    stability = check_positive_definiteness(hydroxylapatite_tensor)
    assert stability.is_stable is True
    assert stability.failed_eigenvalues.size == 0
    unstable = np.diag([1.0, 2.0, 3.0, -1.0, 5.0, 6.0])
    result = check_positive_definiteness(unstable)
    assert result.is_stable is False
    np.testing.assert_allclose(result.failed_eigenvalues, [-1.0])
