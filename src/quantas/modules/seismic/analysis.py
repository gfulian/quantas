# -*- coding: utf-8 -*-

"""High-level numerical operations for seismic workflows."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from quantas.core.geometry import (
    Hemisphere,
    create_spherical_grid,
    tensor_frame_mapping,
)
from quantas.core.physics.elasticity import (
    ElasticTensor,
    StiffnessSymmetryCriterion,
    check_positive_definiteness,
    compute_elastic_averages,
    validate_stiffness_matrix,
)
from quantas.core.physics.seismic import (
    ChristoffelSolver,
    ElasticMedium,
    SamplingLevel,
    isotropic_seismic_velocities,
    sample_seismic_field,
)
from quantas.modules.seismic.models import (
    SeismicInput,
    SeismicOptions,
    SeismicResult,
)


ProgressCallback = Callable[[int, int], None]


def validate_input(input_data: SeismicInput, options: SeismicOptions) -> None:
    """Validate seismic input data and scientific options.

    Parameters
    ----------
    input_data : SeismicInput
        Stiffness matrix and density.
    options : SeismicOptions
        Sampling and numerical options.

    Raises
    ------
    ValueError
        If the input matrix, density or options are invalid.
    """
    if input_data.stiffness is None:
        raise ValueError("The elastic stiffness matrix is missing.")
    validate_stiffness_matrix(
        input_data.stiffness,
        symmetry_tolerance=1.0e-8,
        symmetry_criterion=StiffnessSymmetryCriterion.ELEMENTWISE,
        copy=False,
    )
    if not np.isfinite(input_data.density) or input_data.density <= 0.0:
        raise ValueError("Density must be finite and positive.")

    create_spherical_grid(
        options.ntheta,
        options.nphi,
        hemisphere=Hemisphere(options.hemisphere),
    )
    SamplingLevel(options.level)
    if (
        isinstance(options.batch_size, bool)
        or int(options.batch_size) != options.batch_size
        or options.batch_size < 1
    ):
        raise ValueError("batch_size must be a positive integer.")

    tolerance_values = {
        "eigenvalue_rtol": options.eigenvalue_rtol,
        "eigenvalue_atol": options.eigenvalue_atol,
        "degeneracy_rtol": options.degeneracy_rtol,
        "degeneracy_atol": options.degeneracy_atol,
        "caustic_rtol": options.caustic_rtol,
        "caustic_atol": options.caustic_atol,
    }
    for name, value in tolerance_values.items():
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative.")
    if (
        not np.isfinite(options.pseudoinverse_rcond)
        or options.pseudoinverse_rcond < 0.0
        or options.pseudoinverse_rcond >= 1.0
    ):
        raise ValueError("pseudoinverse_rcond must be in the interval [0, 1).")


def run_analysis(
    input_data: SeismicInput,
    options: SeismicOptions,
    *,
    progress_callback: ProgressCallback | None = None,
) -> SeismicResult:
    """Calculate a complete sampled seismic field.

    Parameters
    ----------
    input_data : SeismicInput
        Validated stiffness matrix and density.
    options : SeismicOptions
        Sampling and numerical options.
    progress_callback : callable or None, optional
        Callback receiving ``(current, total)`` after each batch.

    Returns
    -------
    SeismicResult
        Sampled acoustic fields and isotropic reference quantities.

    Raises
    ------
    ValueError
        If the stiffness matrix is not positive definite.
    """
    validate_input(input_data, options)
    assert input_data.stiffness is not None
    source_stiffness = np.asarray(input_data.stiffness, dtype=float)
    source_tensor = ElasticTensor(source_stiffness)
    tensor = (
        source_tensor
        if options.rotation is None
        else source_tensor.rotate(options.rotation.matrix)
    )
    stiffness = tensor.stiffness
    stability = check_positive_definiteness(tensor)
    if not stability.is_stable:
        raise ValueError(
            "Seismic propagation requires a positive-definite stiffness matrix."
        )

    averages = compute_elastic_averages(tensor)
    medium = ElasticMedium(tensor, input_data.density)
    isotropic = isotropic_seismic_velocities(averages, medium.density)
    solver = ChristoffelSolver(
        medium,
        eigenvalue_rtol=options.eigenvalue_rtol,
        eigenvalue_atol=options.eigenvalue_atol,
        degeneracy_rtol=options.degeneracy_rtol,
        degeneracy_atol=options.degeneracy_atol,
        pseudoinverse_rcond=options.pseudoinverse_rcond,
        caustic_rtol=options.caustic_rtol,
        caustic_atol=options.caustic_atol,
    )
    grid = create_spherical_grid(
        options.ntheta,
        options.nphi,
        hemisphere=options.hemisphere,
    )
    field = sample_seismic_field(
        solver,
        grid,
        level=options.level,
        batch_size=options.batch_size,
        track_polarization_axes=options.track_polarization_axes,
        progress_callback=progress_callback,
    )

    metadata = {
        "tensor_frame": tensor_frame_mapping(options.rotation),
        "n_points": field.n_points,
        "invalid_phase_points": int(np.count_nonzero(~field.phase.valid_mask)),
        "degenerate_mode_points": int(np.count_nonzero(field.phase.degeneracy_mask)),
        "shear_acoustic_axis_candidates": int(
            np.count_nonzero(field.phase.shear_axis_candidate_mask)
        ),
    }
    if field.enhancement is not None:
        metadata["caustic_candidates"] = int(
            np.count_nonzero(field.enhancement.caustic_candidate_mask)
        )
        metadata["non_finite_enhancement_points"] = int(
            np.count_nonzero(
                field.enhancement.valid_mask & ~field.enhancement.finite_mask
            )
        )

    return SeismicResult(
        jobname=input_data.jobname,
        density=float(input_data.density),
        stiffness=np.array(stiffness, dtype=float, copy=True),
        stability=stability,
        averages=averages,
        isotropic_velocities=isotropic,
        grid=grid,
        field=field,
        metadata=metadata,
    )
