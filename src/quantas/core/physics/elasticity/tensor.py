# -*- coding: utf-8 -*-

"""Elastic-tensor representations and low-level tensor transformations."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .conventions import (
    VOIGT_INDEX_MAP,
    rotate_voigt_stiffness,
    voigt_compliance_to_cartesian,
    voigt_stiffness_to_cartesian,
)
from .validation import invert_stiffness_matrix, validate_stiffness_matrix


class ElasticTensor:
    """Represent a second-order elastic stiffness tensor.

    Parameters
    ----------
    stiffness : array_like
        Elastic stiffness matrix in Voigt notation, with shape ``(6, 6)``.
        Values are expressed in GPa.

    Raises
    ------
    ValueError
        If the stiffness matrix is invalid or singular.
    """

    def __init__(self, stiffness: ArrayLike) -> None:
        self._stiffness = validate_stiffness_matrix(stiffness, copy=True)
        self._compliance = invert_stiffness_matrix(self._stiffness)
        self._stiffness_tensor = voigt_stiffness_to_cartesian(self._stiffness)
        self._compliance_tensor = voigt_compliance_to_cartesian(self._compliance)

        for array in (
            self._stiffness,
            self._compliance,
            self._stiffness_tensor,
            self._compliance_tensor,
        ):
            array.setflags(write=False)

    @property
    def stiffness(self) -> NDArray[np.float64]:
        """Return the ``6 x 6`` stiffness matrix in Voigt notation."""
        return self._stiffness

    @property
    def compliance(self) -> NDArray[np.float64]:
        """Return the ``6 x 6`` compliance matrix in Voigt notation."""
        return self._compliance

    @property
    def stiffness_tensor(self) -> NDArray[np.float64]:
        """Return the full Cartesian stiffness tensor ``C_ijkl``."""
        return self._stiffness_tensor

    @property
    def compliance_tensor(self) -> NDArray[np.float64]:
        """Return the full Cartesian compliance tensor ``S_ijkl``."""
        return self._compliance_tensor

    @property
    def eigenvalues(self) -> NDArray[np.float64]:
        """Return the eigenvalues of the stiffness matrix."""
        return np.asarray(np.linalg.eigvalsh(self._stiffness), dtype=np.float64)

    def stiffness_component(self, indexes: tuple[int, int, int, int]) -> float:
        """Return one Cartesian stiffness-tensor component.

        Parameters
        ----------
        indexes : tuple of int
            Four Cartesian indexes ``(i, j, k, l)``.

        Returns
        -------
        float
            Selected stiffness component.
        """
        i, j, k, ell = indexes
        return float(self._stiffness[VOIGT_INDEX_MAP[i, j], VOIGT_INDEX_MAP[k, ell]])

    def compliance_component(self, indexes: tuple[int, int, int, int]) -> float:
        """Return one Cartesian compliance-tensor component.

        Parameters
        ----------
        indexes : tuple of int
            Four Cartesian indexes ``(i, j, k, l)``.

        Returns
        -------
        float
            Selected compliance component, including engineering-shear factors.
        """
        i, j, k, ell = indexes
        return float(self._compliance_tensor[i, j, k, ell])

    def rotate(self, rotation: ArrayLike) -> "ElasticTensor":
        """Return tensor components transformed to a rotated Cartesian frame.

        The component transformation follows

        ``C'_{ijkl} = R_ia R_jb R_kc R_ld C_abcd``.

        Parameters
        ----------
        rotation : array_like
            Proper orthogonal ``3 x 3`` transformation matrix. When its rows
            are target basis vectors expressed in the current frame, the result
            describes the same physical tensor in that target frame.

        Returns
        -------
        ElasticTensor
            General elastic tensor in the transformed frame. A symmetry-tagged
            subclass is deliberately not preserved because an arbitrary frame
            need not retain its canonical sparse matrix form.

        Raises
        ------
        ValueError
            If the rotation matrix is invalid.
        """
        stiffness = rotate_voigt_stiffness(self._stiffness, rotation)
        return ElasticTensor(stiffness)


class OrthorhombicElasticTensor(ElasticTensor):
    """Elastic tensor tagged for orthorhombic directional specializations.

    Parameters
    ----------
    tensor : ElasticTensor
        General tensor to specialize.

    Raises
    ------
    TypeError
        If ``tensor`` is not an :class:`ElasticTensor` instance.
    """

    def __init__(self, tensor: ElasticTensor) -> None:
        if not isinstance(tensor, ElasticTensor):
            raise TypeError(
                "OrthorhombicElasticTensor requires an ElasticTensor instance."
            )
        super().__init__(tensor.stiffness)
