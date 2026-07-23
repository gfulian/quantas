# -*- coding: utf-8 -*-

"""High-level numerical operations for the elasticity workflow."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict

import numpy as np

from quantas.core.geometry import TensorRotation
from quantas.core.physics.elasticity import (
    ElasticTensor,
    check_positive_definiteness,
    compute_elastic_averages,
    detect_elastic_symmetry,
    find_directional_extrema,
    linear_compressibility,
    poisson_ratio,
    sample_elastic_directional_field,
    shear_modulus,
    specialize_elastic_tensor,
    validate_stiffness_matrix,
    young_modulus,
)
from quantas.modules.elasticity.models import (
    ElasticityInput,
    ElasticityOptions,
    ElasticityResult,
)


ProgressCallback = Callable[[str, int, int], None]
StepCallback = Callable[[str, str], None]


def validate_input(input_data: ElasticityInput, options: ElasticityOptions) -> None:
    """Validate elasticity input and scientific options.

    Parameters
    ----------
    input_data : ElasticityInput
        Input data to validate.
    options : ElasticityOptions
        Scientific options to validate.

    Raises
    ------
    ValueError
        If the stiffness matrix or options are invalid.
    """
    if input_data.stiffness is None:
        raise ValueError("The elastic stiffness matrix is missing.")

    validate_stiffness_matrix(input_data.stiffness, copy=False)
    if options.pressure_unit != "GPa":
        raise ValueError(
            "Elastic stiffness values are currently supported only in GPa."
        )
    if options.ntheta < 2:
        raise ValueError("ntheta must be at least 2.")


def create_elastic_tensor(
    input_data: ElasticityInput,
    rotation: TensorRotation | None = None,
) -> ElasticTensor:
    """Create the tensor used by the elasticity analysis.

    Parameters
    ----------
    input_data : ElasticityInput
        Source stiffness matrix.
    rotation : TensorRotation or None, optional
        Optional source-to-analysis component transformation.

    Returns
    -------
    ElasticTensor
        Source tensor or transformed tensor used by the workflow.
    """
    stiffness = input_data.stiffness
    if stiffness is None:
        raise ValueError("The elastic stiffness matrix is missing.")
    tensor = ElasticTensor(stiffness)
    if rotation is None:
        return tensor
    return tensor.rotate(rotation.matrix)


def calculate_basic_properties(
    tensor: ElasticTensor,
    input_data: ElasticityInput,
) -> ElasticityResult:
    """Calculate symmetry, averages, compliance, and stability."""
    result = ElasticityResult(
        jobname=input_data.jobname,
        stiffness=tensor.stiffness.copy(),
        compliance=tensor.compliance.copy(),
        averages=compute_elastic_averages(tensor),
        stability=check_positive_definiteness(tensor),
    )
    result.crystal_system = detect_elastic_symmetry(tensor.stiffness)
    return result


def specialize_tensor(
    tensor: ElasticTensor,
    result: ElasticityResult,
) -> ElasticTensor:
    """Select an available symmetry-specific tensor representation."""
    return specialize_elastic_tensor(tensor, result.crystal_system or "triclinic")


def calculate_directional_variations(
    tensor: ElasticTensor,
    result: ElasticityResult,
) -> None:
    """Calculate global extrema of the four directional elastic properties."""
    result.add_variation(
        "young_modulus",
        find_directional_extrema(
            lambda angles: young_modulus(tensor, angles.tolist()), 2
        ),
    )
    result.add_variation(
        "linear_compressibility",
        find_directional_extrema(
            lambda angles: linear_compressibility(tensor, angles.tolist()), 2
        ),
    )
    result.add_variation(
        "shear_modulus",
        find_directional_extrema(
            lambda angles: shear_modulus(tensor, angles.tolist()), 3
        ),
    )
    result.add_variation(
        "poisson_ratio",
        find_directional_extrema(
            lambda angles: poisson_ratio(tensor, angles.tolist()), 3
        ),
    )


def calculate_2d_properties(
    tensor: ElasticTensor,
    result: ElasticityResult,
    options: ElasticityOptions,
    progress_callback: ProgressCallback | None = None,
    step_callback: StepCallback | None = None,
) -> None:
    """Calculate elastic properties on the principal Cartesian planes.

    All four property families are evaluated in one batched contraction per
    plane.  Transverse shear and Poisson extrema are solved algebraically and
    therefore do not depend on local optimizer convergence.
    """
    if not options.calculate_2d:
        return

    planes = create_principal_plane_grids(options.ntheta)
    property_names = (
        "young_modulus",
        "linear_compressibility",
        "shear_modulus",
        "poisson_ratio",
    )

    for plane, angles in planes.items():
        theta = angles["theta"]
        phi = angles["phi"]
        result.add_2d_data(plane, "theta", theta)
        result.add_2d_data(plane, "phi", phi)

        for property_name in property_names:
            _notify_step(plane, property_name, step_callback)

        field = sample_elastic_directional_field(
            tensor,
            theta,
            phi,
            progress_callback=_wrap_progress(
                f"{plane}: directional elastic field",
                progress_callback,
            ),
        )
        assert field.young_modulus is not None
        assert field.linear_compressibility is not None
        assert field.shear_minimum is not None
        assert field.shear_maximum is not None
        assert field.poisson_minimum is not None
        assert field.poisson_maximum is not None

        compressibility = field.linear_compressibility
        poisson_minimum = field.poisson_minimum
        result.add_2d_data(
            plane,
            "young_modulus",
            np.array(field.young_modulus, dtype=float, copy=True),
        )
        result.add_2d_data(
            plane,
            "linear_compressibility",
            np.column_stack(
                (
                    np.maximum(compressibility, 0.0),
                    np.maximum(-compressibility, 0.0),
                )
            ),
        )
        result.add_2d_data(
            plane,
            "shear_modulus",
            np.column_stack((field.shear_minimum, field.shear_maximum)),
        )
        result.add_2d_data(
            plane,
            "poisson_ratio",
            np.column_stack(
                (
                    np.minimum(poisson_minimum, 0.0),
                    np.maximum(poisson_minimum, 0.0),
                    field.poisson_maximum,
                )
            ),
        )
        result.metadata.setdefault("sampling_diagnostics_2d", {})[plane] = asdict(
            field.diagnostics
        )


def create_principal_plane_grids(points: int) -> dict[str, dict[str, np.ndarray]]:
    """Create closed angular grids for the ``xy``, ``xz``, and ``yz`` planes."""
    angles = np.linspace(0.0, 2.0 * np.pi, points, endpoint=True)
    return {
        "xy": {
            "theta": np.full(points, np.pi / 2.0, dtype=float),
            "phi": angles.copy(),
        },
        "xz": {
            "theta": angles.copy(),
            "phi": np.zeros(points, dtype=float),
        },
        "yz": {
            "theta": angles.copy(),
            "phi": np.full(points, np.pi / 2.0, dtype=float),
        },
    }


def _wrap_progress(
    label: str,
    callback: ProgressCallback | None,
) -> Callable[[int, int], None] | None:
    """Adapt a workflow progress callback to the numerical callback contract."""
    if callback is None:
        return None

    def wrapped(current: int, total: int) -> None:
        callback(label, current, total)

    return wrapped


def _notify_step(
    plane: str,
    property_name: str,
    callback: StepCallback | None,
) -> None:
    """Notify the workflow that one two-dimensional property calculation is starting."""
    if callback is not None:
        callback(plane, property_name)
