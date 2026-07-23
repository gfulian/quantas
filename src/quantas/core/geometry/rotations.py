# -*- coding: utf-8 -*-

"""Proper rotations and Cartesian tensor transformations."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray


DEFAULT_ROTATION_TOLERANCE = 1.0e-12


class TensorRotationKind(str, Enum):
    """Supported user descriptions of a tensor-component transformation."""

    MATRIX = "matrix"
    XYZ = "xyz"


@dataclass(frozen=True, slots=True)
class TensorRotation:
    """User-defined transformation from the source to the analysis frame.

    The stored matrix follows the Quantas component convention

    ``T'_{ij...} = R_ia R_jb ... T_ab...``.

    Parameters
    ----------
    matrix : array_like
        Proper orthogonal matrix with shape ``(3, 3)``. Its rows are the
        analysis-frame basis vectors expressed in the source frame.
    kind : TensorRotationKind or str, optional
        Description used to construct the matrix.
    angles : tuple of float or None, optional
        Input angles for an ``xyz`` construction.
    angle_unit : str or None, optional
        Unit of ``angles``. Currently ``"degree"`` or ``"radian"``.
    description : str or None, optional
        Optional user-facing description saved with result provenance.

    Raises
    ------
    ValueError
        If the matrix or angular provenance is invalid.
    """

    matrix: NDArray[np.float64]
    kind: TensorRotationKind = TensorRotationKind.MATRIX
    angles: tuple[float, float, float] | None = None
    angle_unit: str | None = None
    description: str | None = None

    def __post_init__(self) -> None:
        matrix = validate_rotation_matrix(self.matrix, copy=True)
        matrix.setflags(write=False)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "kind", TensorRotationKind(self.kind))

        if self.angles is None:
            if self.kind is TensorRotationKind.XYZ:
                raise ValueError("XYZ tensor rotations require three angles.")
            if self.angle_unit is not None:
                raise ValueError("angle_unit requires explicit rotation angles.")
            return

        angles = tuple(float(value) for value in self.angles)
        if len(angles) != 3 or not np.all(np.isfinite(angles)):
            raise ValueError("Rotation angles must contain three finite values.")
        if self.kind is not TensorRotationKind.XYZ:
            raise ValueError("Rotation angles are valid only for kind='xyz'.")
        if self.angle_unit not in {"degree", "radian"}:
            raise ValueError("angle_unit must be 'degree' or 'radian'.")
        object.__setattr__(self, "angles", angles)

    @classmethod
    def from_matrix(
        cls,
        matrix: ArrayLike,
        *,
        description: str | None = None,
    ) -> "TensorRotation":
        """Build a transformation from an explicit rotation matrix.

        Parameters
        ----------
        matrix : array_like
            Proper orthogonal ``3 x 3`` component-transformation matrix.
        description : str or None, optional
            Optional provenance note.

        Returns
        -------
        TensorRotation
            Validated matrix transformation.
        """
        return cls(
            matrix=np.asarray(matrix, dtype=float),
            kind=TensorRotationKind.MATRIX,
            description=description,
        )

    @classmethod
    def from_xyz(
        cls,
        x: float,
        y: float,
        z: float,
        *,
        degrees: bool = True,
        description: str | None = None,
    ) -> "TensorRotation":
        """Build a transformation from ordered right-handed XYZ rotations.

        The rotations are applied about the fixed source axes in the order
        ``x``, then ``y``, then ``z``. For column-vector matrices this gives
        ``R = Rz(z) @ Ry(y) @ Rx(x)``. The resulting matrix is used directly
        in the Quantas tensor-component transformation.

        Parameters
        ----------
        x, y, z : float
            Rotation angles about the fixed source Cartesian axes.
        degrees : bool, optional
            Interpret the supplied angles as degrees when ``True``.
        description : str or None, optional
            Optional provenance note.

        Returns
        -------
        TensorRotation
            Validated XYZ transformation.
        """
        angles = (float(x), float(y), float(z))
        matrix = xyz_rotation_matrix(angles, degrees=degrees)
        return cls(
            matrix=matrix,
            kind=TensorRotationKind.XYZ,
            angles=angles,
            angle_unit="degree" if degrees else "radian",
            description=description,
        )

    def as_mapping(self) -> dict[str, Any]:
        """Return a serialization-friendly provenance mapping.

        Returns
        -------
        dict
            Rotation kind, convention, matrix and optional angular input.
        """
        data: dict[str, Any] = {
            "kind": self.kind.value,
            "component_transform": self.matrix.copy(),
            "convention": "C'_ijkl = R_ia R_jb R_kc R_ld C_abcd",
            "source_frame": "source",
            "analysis_frame": "rotated",
        }
        if self.angles is not None:
            data["angles"] = self.angles
            data["angle_unit"] = self.angle_unit
            data["sequence"] = "fixed_xyz"
        if self.description is not None:
            data["description"] = self.description
        return data


def tensor_frame_mapping(
    rotation: TensorRotation | None,
) -> dict[str, Any]:
    """Return source and analysis frame provenance for a workflow run.

    Parameters
    ----------
    rotation : TensorRotation or None
        Optional transformation applied before numerical analysis.

    Returns
    -------
    dict
        Serialization-friendly source/analysis frame metadata.
    """
    if rotation is None:
        return {
            "source_frame": "source",
            "analysis_frame": "source",
            "transformed": False,
        }
    data = rotation.as_mapping()
    data["transformed"] = True
    return data


def validate_rotation_matrix(
    rotation: ArrayLike,
    *,
    tolerance: float = DEFAULT_ROTATION_TOLERANCE,
    copy: bool = True,
) -> NDArray[np.float64]:
    """Validate a proper orthogonal rotation matrix.

    Parameters
    ----------
    rotation : array_like
        Candidate matrix with shape ``(3, 3)``.
    tolerance : float, optional
        Absolute tolerance used for orthogonality and determinant checks.
    copy : bool, optional
        Whether to return an independent copy.

    Returns
    -------
    numpy.ndarray
        Validated rotation matrix.

    Raises
    ------
    ValueError
        If the matrix is malformed, non-finite, non-orthogonal or left-handed.
    """
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("tolerance must be finite and non-negative.")
    if copy:
        matrix = np.array(rotation, dtype=float, copy=True)
    else:
        matrix = np.asarray(rotation, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError("The rotation matrix must have shape (3, 3).")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("The rotation matrix must contain finite values.")
    if not np.allclose(
        matrix @ matrix.T,
        np.eye(3),
        atol=tolerance,
        rtol=0.0,
    ):
        raise ValueError("The rotation matrix must be orthogonal.")
    if not np.isclose(
        np.linalg.det(matrix),
        1.0,
        atol=tolerance,
        rtol=0.0,
    ):
        raise ValueError("The rotation matrix must be right-handed.")
    return matrix


def axis_angle_rotation(
    axis: ArrayLike,
    angle: float,
) -> NDArray[np.float64]:
    """Return a right-handed active rotation from an axis and angle.

    Parameters
    ----------
    axis : array_like
        Non-zero Cartesian rotation axis.
    angle : float
        Rotation angle in radians.

    Returns
    -------
    numpy.ndarray
        Proper orthogonal rotation matrix.

    Raises
    ------
    ValueError
        If the axis or angle is invalid.
    """
    vector = np.asarray(axis, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError("The rotation axis must be a finite vector of shape (3,).")
    if not np.isfinite(angle):
        raise ValueError("The rotation angle must be finite.")
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        raise ValueError("The rotation axis must have non-zero length.")
    x, y, z = vector / norm
    cosine = float(np.cos(angle))
    sine = float(np.sin(angle))
    complement = 1.0 - cosine
    matrix = np.asarray(
        [
            [
                cosine + x * x * complement,
                x * y * complement - z * sine,
                x * z * complement + y * sine,
            ],
            [
                y * x * complement + z * sine,
                cosine + y * y * complement,
                y * z * complement - x * sine,
            ],
            [
                z * x * complement - y * sine,
                z * y * complement + x * sine,
                cosine + z * z * complement,
            ],
        ],
        dtype=float,
    )
    matrix.setflags(write=False)
    return matrix


def xyz_rotation_matrix(
    angles: ArrayLike,
    *,
    degrees: bool = True,
) -> NDArray[np.float64]:
    """Return the ordered fixed-axis XYZ rotation matrix.

    Parameters
    ----------
    angles : array_like
        Three angles ``(x, y, z)``.
    degrees : bool, optional
        Interpret angles as degrees when ``True``.

    Returns
    -------
    numpy.ndarray
        Matrix ``Rz(z) @ Ry(y) @ Rx(x)``.

    Raises
    ------
    ValueError
        If the angles are malformed or non-finite.
    """
    values = np.asarray(angles, dtype=float)
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise ValueError("XYZ rotation angles must be a finite vector of shape (3,).")
    if degrees:
        values = np.radians(values)
    rx = axis_angle_rotation((1.0, 0.0, 0.0), float(values[0]))
    ry = axis_angle_rotation((0.0, 1.0, 0.0), float(values[1]))
    rz = axis_angle_rotation((0.0, 0.0, 1.0), float(values[2]))
    matrix = np.asarray(rz @ ry @ rx, dtype=float)
    matrix.setflags(write=False)
    return matrix


def rotate_cartesian_rank4(
    tensor: ArrayLike,
    rotation: ArrayLike,
) -> NDArray[np.float64]:
    """Transform a fourth-rank Cartesian tensor with a proper rotation.

    The implemented convention is

    ``T'_{ijkl} = R_ia R_jb R_kc R_ld T_abcd``.

    Parameters
    ----------
    tensor : array_like
        Cartesian fourth-rank tensor with shape ``(3, 3, 3, 3)``.
    rotation : array_like
        Proper orthogonal rotation matrix.

    Returns
    -------
    numpy.ndarray
        Transformed fourth-rank tensor.

    Raises
    ------
    ValueError
        If the tensor or rotation is invalid.
    """
    cartesian = np.asarray(tensor, dtype=float)
    if cartesian.shape != (3, 3, 3, 3):
        raise ValueError("The Cartesian tensor must have shape (3, 3, 3, 3).")
    if not np.all(np.isfinite(cartesian)):
        raise ValueError("The Cartesian tensor must contain finite values.")
    matrix = validate_rotation_matrix(rotation, copy=False)
    return np.einsum(
        "ia,jb,kc,ld,abcd->ijkl",
        matrix,
        matrix,
        matrix,
        matrix,
        cartesian,
        optimize=True,
    )
