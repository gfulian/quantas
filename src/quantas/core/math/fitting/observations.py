# -*- coding: utf-8 -*-

"""Passive observation contracts shared by numerical fitting algorithms."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(slots=True)
class FitObservations:
    """Store coordinates, observations, uncertainties, and selection mask.

    Parameters
    ----------
    x, y : array-like
        Independent coordinates and observed dependent values. ``x`` may be
        a vector for one explanatory coordinate or a matrix with shape
        ``(n_coordinates, n_observations)``. ``y`` is always a vector.
    sigma_x, sigma_y : array-like or None, optional
        Standard uncertainties matching the corresponding data shape. Zero
        values are retained in the passive data
        contract but methods that use a sigma as a divisor must require it to
        be strictly positive for selected observations.
    mask : array-like of bool or None, optional
        Selection mask. Omitted masks select every observation.
    x_name, y_name : str, optional
        Stable quantity names.
    x_unit, y_unit : str or None, optional
        Unit labels.
    metadata : dict, optional
        Additional passive provenance and method-independent information.

    Raises
    ------
    ValueError
        If arrays are non-finite, inconsistent, or the mask is empty.
    """

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    sigma_x: NDArray[np.float64] | None = None
    sigma_y: NDArray[np.float64] | None = None
    mask: NDArray[np.bool_] | None = None
    x_name: str = "x"
    y_name: str = "y"
    x_unit: str | None = None
    y_unit: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays to stable NumPy dtypes and validate the contract."""
        self.x = _coordinate_array(self.x, "x")
        self.y = _vector(self.y, "y")
        if self.x.shape[-1] != self.y.size:
            raise ValueError(
                "fit observation x must contain one coordinate column per y value"
            )
        self.sigma_x = _sigma(self.sigma_x, self.x.shape, "sigma_x", coordinate=True)
        self.sigma_y = _sigma(self.sigma_y, self.y.shape, "sigma_y")
        if self.mask is None:
            self.mask = np.ones(self.y.shape, dtype=np.bool_)
        else:
            mask = np.asarray(self.mask, dtype=np.bool_)
            if mask.ndim != 1 or mask.shape != self.y.shape:
                raise ValueError(
                    "fit observation mask must match the observation shape"
                )
            self.mask = mask.copy()
        if not np.any(self.mask):
            raise ValueError("fit observation mask must select at least one point")
        self.x_name = str(self.x_name)
        self.y_name = str(self.y_name)
        self.x_unit = None if self.x_unit is None else str(self.x_unit)
        self.y_unit = None if self.y_unit is None else str(self.y_unit)
        self.metadata = dict(self.metadata)

    @property
    def size(self) -> int:
        """Return the total number of observations."""
        return int(self.y.size)

    @property
    def selected_size(self) -> int:
        """Return the number of selected observations."""
        assert self.mask is not None
        return int(np.count_nonzero(self.mask))

    def selected(self) -> FitObservations:
        """Return a compact observation object containing selected points only."""
        assert self.mask is not None
        selected = self.mask
        return FitObservations(
            x=(self.x[selected] if self.x.ndim == 1 else self.x[:, selected]),
            y=self.y[selected],
            sigma_x=(
                None
                if self.sigma_x is None
                else (
                    self.sigma_x[selected]
                    if self.sigma_x.ndim == 1
                    else self.sigma_x[:, selected]
                )
            ),
            sigma_y=None if self.sigma_y is None else self.sigma_y[selected],
            x_name=self.x_name,
            y_name=self.y_name,
            x_unit=self.x_unit,
            y_unit=self.y_unit,
            metadata={**self.metadata, "source_mask": self.mask.copy()},
        )

    def require_positive_sigma(self, axis: str) -> NDArray[np.float64]:
        """Return a selected uncertainty vector suitable for weighting.

        Parameters
        ----------
        axis : {"x", "y"}
            Requested uncertainty axis.

        Returns
        -------
        ndarray
            Strictly positive selected standard uncertainties.

        Raises
        ------
        ValueError
            If the uncertainty is absent or contains selected zero values.
        """
        if axis not in {"x", "y"}:
            raise ValueError("sigma axis must be 'x' or 'y'")
        sigma = self.sigma_x if axis == "x" else self.sigma_y
        if sigma is None:
            raise ValueError(f"sigma_{axis} is required for this fit method")
        assert self.mask is not None
        selected = sigma[self.mask] if sigma.ndim == 1 else sigma[:, self.mask]
        if np.any(selected <= 0.0):
            raise ValueError(f"selected sigma_{axis} values must be strictly positive")
        return selected.copy()

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable observation representation."""
        return {
            "x": self.x.tolist(),
            "y": self.y.tolist(),
            "sigma_x": None if self.sigma_x is None else self.sigma_x.tolist(),
            "sigma_y": None if self.sigma_y is None else self.sigma_y.tolist(),
            "mask": None if self.mask is None else self.mask.tolist(),
            "x_name": self.x_name,
            "y_name": self.y_name,
            "x_unit": self.x_unit,
            "y_unit": self.y_unit,
            "metadata": dict(self.metadata),
        }


def _vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    """Return a finite one-dimensional ``float64`` array."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"fit observation {name} must be a non-empty vector")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"fit observation {name} must contain finite values")
    return array.copy()


def _coordinate_array(values: ArrayLike, name: str) -> NDArray[np.float64]:
    """Return finite one- or two-dimensional explanatory coordinates."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim not in {1, 2} or array.size == 0 or array.shape[-1] == 0:
        raise ValueError(
            f"fit observation {name} must be a non-empty vector or coordinate matrix"
        )
    if not np.all(np.isfinite(array)):
        raise ValueError(f"fit observation {name} must contain finite values")
    return array.copy()


def _sigma(
    values: ArrayLike | None,
    shape: tuple[int, ...],
    name: str,
    *,
    coordinate: bool = False,
) -> NDArray[np.float64] | None:
    """Normalize optional standard uncertainties."""
    if values is None:
        return None
    array = _coordinate_array(values, name) if coordinate else _vector(values, name)
    if array.shape != shape:
        raise ValueError(f"fit observation {name} must match the data shape")
    if np.any(array < 0.0):
        raise ValueError(f"fit observation {name} cannot contain negative values")
    return array
