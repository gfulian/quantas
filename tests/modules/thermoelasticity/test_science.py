"""Scientific support and stability tests for thermoelastic QSA."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    evaluate_stability_field,
    rotate_voigt_stiffness,
)
from quantas.modules.thermoelasticity.frames import normalize_elastic_frame
from quantas.modules.thermoelasticity.fitting import ColdFiniteStrainComponentModel
from quantas.modules.thermoelasticity.quality import assess_component_fit_quality
from quantas.modules.thermoelasticity.io.reader import read_thermoelastic_input


def test_frame_corotation_recovers_reference_components() -> None:
    """Rigidly rotated lattice and stiffness return to one reference frame."""
    reference_lattice = np.asarray(
        [[4.2, 0.0, 0.0], [0.3, 5.1, 0.0], [0.1, 0.4, 6.0]],
        dtype=np.float64,
    )
    stiffness = np.asarray(
        [
            [240.0, 82.0, 71.0, 2.0, 3.0, 4.0],
            [82.0, 210.0, 68.0, 5.0, 1.0, 2.0],
            [71.0, 68.0, 190.0, 4.0, 3.0, 1.0],
            [2.0, 5.0, 4.0, 65.0, 2.0, 1.0],
            [3.0, 1.0, 3.0, 2.0, 72.0, 2.5],
            [4.0, 2.0, 1.0, 1.0, 2.5, 78.0],
        ],
        dtype=np.float64,
    )
    angle = np.deg2rad(31.0)
    rotation = np.asarray(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    sampled_lattice = reference_lattice @ rotation.T
    sampled_stiffness = rotate_voigt_stiffness(stiffness, rotation)

    normalized = normalize_elastic_frame(
        sampled_lattice,
        sampled_stiffness,
        reference_lattice,
    )

    np.testing.assert_allclose(normalized.lattice, reference_lattice, atol=1.0e-11)
    np.testing.assert_allclose(normalized.stiffness, stiffness, atol=1.0e-9)
    np.testing.assert_allclose(normalized.removed_rotation_degrees, 31.0, atol=1.0e-10)


def test_quality_classifies_well_supported_and_narrow_paths() -> None:
    """Coverage diagnostics distinguish robust and weakly supported paths."""
    model = ColdFiniteStrainComponentModel(
        reference_volume=100.0,
        bulk_modulus=120.0,
        bulk_modulus_derivative=4.0,
        wallace_delta=-3.0,
        order=3,
        label="C11",
    )
    volumes = np.linspace(90.0, 110.0, 9)
    design = np.column_stack(
        (
            model.evaluate(volumes, [1.0, 0.0]) - model.evaluate(volumes, [0.0, 0.0]),
            model.evaluate(volumes, [0.0, 1.0]) - model.evaluate(volumes, [0.0, 0.0]),
        )
    )
    observed = model.evaluate(volumes, [220.0, 6.0])
    quality = assess_component_fit_quality(
        volumes=volumes,
        observed=observed,
        symmetry_spread=np.zeros_like(volumes),
        reference_volume=100.0,
        design_matrix=design,
        fitted_parameters=np.asarray([220.0, 6.0]),
        leave_one_out_parameters=np.tile([220.0, 6.0], (volumes.size, 1)),
        alternate_order_parameters=np.asarray([220.0, 6.0]),
        minimum_points=4,
        minimum_strain_span=5.0e-3,
        maximum_design_condition=1.0e6,
        maximum_relative_symmetry_spread=1.0e-2,
        maximum_leave_one_out_change=0.5,
        maximum_order_sensitivity=0.5,
        relative_floor=1.0e-8,
    )
    assert quality.level == "supported"
    assert quality.reference_volume_bracketed

    narrow = volumes[:4] * 0.0 + np.asarray([99.99, 99.995, 100.005, 100.01])
    narrow_design = np.column_stack(
        (
            model.evaluate(narrow, [1.0, 0.0]) - model.evaluate(narrow, [0.0, 0.0]),
            model.evaluate(narrow, [0.0, 1.0]) - model.evaluate(narrow, [0.0, 0.0]),
        )
    )
    narrow_quality = assess_component_fit_quality(
        volumes=narrow,
        observed=model.evaluate(narrow, [220.0, 6.0]),
        symmetry_spread=np.zeros_like(narrow),
        reference_volume=100.0,
        design_matrix=narrow_design,
        fitted_parameters=np.asarray([220.0, 6.0]),
        leave_one_out_parameters=np.tile([220.0, 6.0], (narrow.size, 1)),
        alternate_order_parameters=np.asarray([220.0, 6.0]),
        minimum_points=4,
        minimum_strain_span=5.0e-3,
        maximum_design_condition=1.0e6,
        maximum_relative_symmetry_spread=1.0e-2,
        maximum_leave_one_out_change=0.5,
        maximum_order_sensitivity=0.5,
        relative_floor=1.0e-8,
    )
    assert narrow_quality.level == "caution"
    assert "narrow_eulerian_strain_span" in narrow_quality.issues


def test_stability_preserves_nonstable_states() -> None:
    """Stability assessment reports states without altering their matrices."""
    stable = np.diag([200.0, 180.0, 160.0, 70.0, 65.0, 60.0])
    unstable = stable.copy()
    unstable[-1, -1] = -2.0
    indeterminate = stable.copy()
    indeterminate[0, 0] = np.nan
    field = np.stack((stable, unstable, indeterminate))
    original = field.copy()

    result = evaluate_stability_field(field)

    assert result.stable_mask.tolist() == [True, False, False]
    assert result.unstable_mask.tolist() == [False, True, False]
    assert result.indeterminate_mask.tolist() == [False, False, True]
    np.testing.assert_equal(field, original)
    assert result.minimum_eigenvalue[1] == -2.0


def test_schema_10_requires_frame_normalization(tmp_path: Path) -> None:
    """Pre-release schema 1.0 rejects inputs without frame provenance."""
    text = """\
schema:
  name: quantas-thermoelastic-input
  version: '1.0'
job: minimal
method: quasistatic
interface: crystal
conventions:
  tensor_orientation: crystal
reference:
  index: 0
  source: point.out
  structure:
    natom: 1
    lattice:
    - [ 4.0, 0.0, 0.0 ]
    - [ 0.0, 4.0, 0.0 ]
    - [ 0.0, 0.0, 4.0 ]
    atomic_numbers: [ 12 ]
    fractional_positions:
    - [ 0.0, 0.0, 0.0 ]
  symmetry:
    space_group_number: 221
    international_symbol: Pm-3m
    elastic_system: cubic
elastic_data:
- source: point.out
  pressure: 0.0
  stress_pressure: 0.0
  volume: 64.0
  density: 3000.0
  energy: -100.0
  lattice:
  - [ 4.0, 0.0, 0.0 ]
  - [ 0.0, 4.0, 0.0 ]
  - [ 0.0, 0.0, 4.0 ]
  stiffness:
  - [ 200.0, 80.0, 80.0, 0.0, 0.0, 0.0 ]
  - [ 80.0, 200.0, 80.0, 0.0, 0.0, 0.0 ]
  - [ 80.0, 80.0, 200.0, 0.0, 0.0, 0.0 ]
  - [ 0.0, 0.0, 0.0, 60.0, 0.0, 0.0 ]
  - [ 0.0, 0.0, 0.0, 0.0, 60.0, 0.0 ]
  - [ 0.0, 0.0, 0.0, 0.0, 0.0, 60.0 ]
"""
    path = tmp_path / "minimal.yaml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="frame is required"):
        read_thermoelastic_input(path)
