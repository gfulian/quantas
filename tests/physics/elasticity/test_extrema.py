"""Tests for directional optimization and elasticity namespaces."""

from __future__ import annotations

import importlib

import numpy as np

from quantas.core.physics.elasticity import find_directional_extrema


def test_directional_extrema_preserve_search_contract() -> None:
    result = find_directional_extrema(
        lambda angles: 2.0 + np.sin(angles[0]) ** 2,
        dimensions=2,
    )
    np.testing.assert_allclose(result.minimum, 2.0, atol=1.0e-8)
    np.testing.assert_allclose(result.maximum, 3.0, atol=1.0e-8)
    np.testing.assert_allclose(result.anisotropy, 1.5, atol=1.0e-8)


def test_quasistatic_namespace_exposes_cold_finite_strain_relations() -> None:
    module = importlib.import_module("quantas.core.physics.elasticity.quasistatic")
    assert module.__all__ == [
        "FiniteStrainOrder",
        "cold_finite_strain_component",
        "cold_finite_strain_component_jacobian",
        "cold_finite_strain_stiffness",
        "eulerian_finite_strain",
        "wallace_hydrostatic_delta_voigt",
    ]


def test_elasticity_namespace_has_no_seismic_models() -> None:
    module = importlib.import_module("quantas.core.physics.elasticity")
    assert not hasattr(module, "IsotropicSeismicVelocities")
    assert not hasattr(module, "isotropic_seismic_velocities")


def test_pure_elastic_tensor_has_no_density_attribute() -> None:
    from quantas.core.physics.elasticity import ElasticTensor

    tensor = ElasticTensor(np.eye(6))
    assert not hasattr(tensor, "density")
