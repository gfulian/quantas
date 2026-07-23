"""End-to-end library workflow tests for elasticity."""

from __future__ import annotations

from quantas.modules.elasticity.calculator import ElasticityCalculator
from quantas.modules.elasticity.models import ElasticityOptions


def test_elasticity_calculator_basic_workflow(hydroxylapatite_input) -> None:
    calculator = ElasticityCalculator(
        hydroxylapatite_input,
        options=ElasticityOptions(calculate_2d=False),
    )
    result = calculator.execute()

    assert calculator.completed is True
    assert result.metadata.module == "elasticity"
    assert result.metadata.method == "second_order"
    assert result.metadata.schema_version == "2.0"
    elasticity = result.results["elasticity"]
    assert elasticity.crystal_system == "hexagonal"
    assert elasticity.averages.as_array().shape == (3, 4)
    assert elasticity.stability.is_stable is True
    assert not hasattr(elasticity, "isotropic_velocities")
    assert not hasattr(elasticity, "density")
    assert set(elasticity.variations) == {
        "young_modulus",
        "linear_compressibility",
        "shear_modulus",
        "poisson_ratio",
    }
