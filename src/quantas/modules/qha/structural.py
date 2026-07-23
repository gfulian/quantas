# -*- coding: utf-8 -*-

"""Structural properties derived from the QHA equilibrium-volume surface.

The functions in this module connect the frontend-neutral structural-path model
with QHA result containers.  The equilibrium volume and selected volumetric
thermal-expansion coefficient remain authoritative QHA quantities.  Cell shape
is reconstructed from the static volume-constrained structural path stored in
the input.
"""

from __future__ import annotations

import numpy as np

from quantas.core.geometry.thermal_expansion import StructuralPathModel
from quantas.models.structures import StructureVolumeSeries
from quantas.modules.qha.models import QHAOptions, QHAResult


def calculate_structural_thermal_expansion(
    result: QHAResult,
    structure: StructureVolumeSeries,
    options: QHAOptions,
) -> None:
    r"""Calculate temperature-dependent cell and expansion tensor data.

    The QHA equilibrium volume :math:`V(P,T)` is obtained from Helmholtz/Gibbs
    free-energy minimization.  The crystallographic cell at that volume is
    reconstructed from the stored volume-constrained structural path.  Axial
    and tensorial expansion coefficients are then evaluated by the chain rule,
    using the selected QHA volumetric expansion coefficient:

    .. math::

        \alpha_a = \frac{\partial \ln a}{\partial \ln V}\,\alpha_V.

    The same construction is applied to the complete deformation gradient, so
    that the tensor trace reproduces :math:`\alpha_V` by construction.

    Parameters
    ----------
    result : QHAResult
        QHA result containing equilibrium volumes and, when available, the
        selected volumetric thermal-expansion coefficient.
    structure : StructureVolumeSeries
        Compact primitive structural path sampled on the QHA volume grid.
    options : QHAOptions
        QHA options controlling the structural-path polynomial degree.

    Returns
    -------
    None
        ``result`` is updated in place.

    Raises
    ------
    ValueError
        If equilibrium volumes and the structural path use inconsistent sampled
        volumes, or if the structural model cannot be built.
    """
    if result.equilibrium_volume is None:
        return
    if result.volume is not None:
        sampled = np.asarray(result.volume, dtype=np.float64)
        structural = np.asarray(structure.volumes, dtype=np.float64)
        if sampled.shape != structural.shape or not np.allclose(
            sampled,
            structural,
            rtol=1.0e-8,
            atol=1.0e-8,
        ):
            raise ValueError(
                "QHA sampled volumes and structural-path volumes are inconsistent"
            )

    model = StructuralPathModel(structure, degree=options.structural_degree)
    sigma_volume = result.uncertainties.get("sigma_VT")
    evaluated = model.evaluate(
        np.asarray(result.equilibrium_volume, dtype=np.float64),
        None
        if result.thermal_expansion is None
        else np.asarray(result.thermal_expansion, dtype=np.float64),
        None if sigma_volume is None else np.asarray(sigma_volume, dtype=np.float64),
    )
    result.equilibrium_lattice = evaluated.lattice
    result.lattice_parameters = evaluated.lattice_parameters
    result.lattice_parameter_derivatives = evaluated.lattice_parameter_derivatives
    result.axial_thermal_expansion = evaluated.axial_thermal_expansion
    result.thermal_expansion_tensor = evaluated.thermal_expansion_tensor
    result.structural_extrapolation_mask = evaluated.extrapolation_mask
    if evaluated.lattice_parameter_uncertainties is not None:
        result.uncertainties["sigma_cell_parameters"] = (
            evaluated.lattice_parameter_uncertainties
        )

    metadata = result.metadata.setdefault("structural_thermal_expansion", {})
    metadata.update(evaluated.metadata)
    metadata.update(
        {
            "automatic": True,
            "available": True,
            "equilibrium_volume_source": "qha_free_energy_minimization",
            "volumetric_expansion_source": result.metadata.get(
                "thermal_expansion",
                {},
            ).get("selected_method", "selected_qha_alphaV"),
            "cell_shape_source": "volume_constrained_static_relaxation_path",
            "approximation": "one_dimensional_volume_structural_path",
            "full_anisotropic_qha": False,
            "lattice_parameter_order": ["a", "b", "c", "alpha", "beta", "gamma"],
            "lattice_parameter_units": [
                "angstrom",
                "angstrom",
                "angstrom",
                "degree",
                "degree",
                "degree",
            ],
            "axial_expansion_order": ["alpha_a", "alpha_b", "alpha_c"],
            "thermal_expansion_unit": "K^-1",
            "cell_parameter_uncertainty_key": (
                "sigma_cell_parameters"
                if evaluated.lattice_parameter_uncertainties is not None
                else None
            ),
            "cell_parameter_uncertainty_order": [
                "sigma_a",
                "sigma_b",
                "sigma_c",
                "sigma_alpha",
                "sigma_beta",
                "sigma_gamma",
            ],
            "cell_parameter_uncertainty_scope": (
                "equilibrium-volume covariance plus structural-path fit covariance; "
                "source DFT structure uncertainties are not included"
            ),
            "extrapolated_points": int(np.count_nonzero(evaluated.extrapolation_mask)),
        }
    )
