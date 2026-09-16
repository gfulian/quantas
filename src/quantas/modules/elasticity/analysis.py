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
    """Create the elastic tensor used by the analysis workflow.

    Parameters
    ----------
    input_data : ElasticityInput
        Source stiffness matrix in Voigt notation and GPa.
    rotation : TensorRotation or None, optional
        Optional source-to-analysis component transformation. The transformation
        changes tensor components, not the underlying physical tensor.

    Returns
    -------
    ElasticTensor
        Elastic tensor expressed in the analysis Cartesian frame.

    Raises
    ------
    ValueError
        If the input does not contain a stiffness matrix.
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
    """Calculate the frame-dependent tensor representation and bulk summaries.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor in the analysis Cartesian frame. Stiffness is expressed in
        GPa and compliance in GPa^-1.
    input_data : ElasticityInput
        Workflow input providing the job description.

    Returns
    -------
    ElasticityResult
        Result containing copies of stiffness and compliance, the detected elastic
        crystal system, Voigt-Reuss-Hill averages, and the positive-definiteness
        diagnostic.
    """
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
    """Return the symmetry-specialized representation used for directional analysis.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor in the analysis Cartesian frame.
    result : ElasticityResult
        Result whose ``crystal_system`` selects the specialization. Missing
        symmetry information falls back to triclinic behavior.

    Returns
    -------
    ElasticTensor
        Symmetry-specialized tensor representation with stiffness in GPa.
    """
    return specialize_elastic_tensor(tensor, result.crystal_system or "triclinic")


def calculate_directional_variations(
    tensor: ElasticTensor,
    result: ElasticityResult,
) -> None:
    """Calculate exact global extrema of the directional elastic properties.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor in the analysis Cartesian frame.
    result : ElasticityResult
        Result object updated in place under ``variations``.

    Notes
    -----
    Young's modulus and linear compressibility depend on one direction. Shear
    modulus and Poisson's ratio additionally optimize over an orthogonal transverse
    direction. The routine stores the extrema returned by the shared elasticity
    core and performs no display rounding.
    """
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
    """Calculate directional elastic properties on the principal Cartesian planes.

    All property families are evaluated through the shared vectorized directional
    field sampler. Transverse shear and Poisson extrema are solved algebraically
    and therefore do not depend on local optimizer convergence.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor in the analysis Cartesian frame.
    result : ElasticityResult
        Result updated in place under ``properties_2d`` and sampling diagnostics.
    options : ElasticityOptions
        Workflow options. No work is performed when ``calculate_2d`` is false.
    progress_callback : callable or None, optional
        Workflow callback receiving ``(label, current, total)`` for numerical
        sampling progress.
    step_callback : callable or None, optional
        Callback receiving ``(plane, property_name)`` before each reported
        property family.

    Notes
    -----
    Angles are stored in radians. Young's modulus and shear modulus are stored in
    GPa, linear compressibility is stored as separated positive and negative
    branches in TPa^-1, and Poisson's ratio is dimensionless.
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
    """Create closed angular grids for the Cartesian principal planes.

    Parameters
    ----------
    points : int
        Number of angular samples per plane, including both equivalent endpoints
        at 0 and 2*pi.

    Returns
    -------
    dict
        Mapping for ``xy``, ``xz`` and ``yz``. Each entry contains ``theta`` and
        ``phi`` arrays of shape ``(points,)`` in radians.
    """
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
