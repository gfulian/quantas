# -*- coding: utf-8 -*-

"""Dataset and coordinate contracts for EOS workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np

from quantas.core.math.fitting import FitObservations

from ._model_utils import as_float64_vector, optional_sigma

EOS_COLUMN_NAMES = frozenset(
    {
        "pressure",
        "sigma_pressure",
        "temperature",
        "sigma_temperature",
        "volume",
        "sigma_volume",
        "energy",
        "sigma_energy",
        "a",
        "sigma_a",
        "b",
        "sigma_b",
        "c",
        "sigma_c",
        "alpha",
        "sigma_alpha",
        "beta",
        "sigma_beta",
        "gamma",
        "sigma_gamma",
    }
)

EOS_TARGET_NAMES = frozenset(
    {"volume", "energy", "a", "b", "c", "alpha", "beta", "gamma"}
)

_SIGMA_TO_VALUE = {
    "sigma_pressure": "pressure",
    "sigma_temperature": "temperature",
    "sigma_volume": "volume",
    "sigma_energy": "energy",
    "sigma_a": "a",
    "sigma_b": "b",
    "sigma_c": "c",
    "sigma_alpha": "alpha",
    "sigma_beta": "beta",
    "sigma_gamma": "gamma",
}
class EOSCrystalSystem(str, Enum):
    """Canonical crystal systems accepted by EOS input files.

    Notes
    -----
    The crystal system is metadata in the current EOS workflows.  It is
    normalized now so future directional ``[UVW]`` and ``(hkl)`` targets can
    rely on stable names without changing the public input contract.
    """

    TRICLINIC = "triclinic"
    MONOCLINIC = "monoclinic"
    ORTHORHOMBIC = "orthorhombic"
    TETRAGONAL = "tetragonal"
    TRIGONAL = "trigonal"
    HEXAGONAL = "hexagonal"
    CUBIC = "cubic"

    @property
    def independent_axes(self) -> tuple[str, ...]:
        """Return independent conventional cell-axis labels."""
        if self is EOSCrystalSystem.CUBIC:
            return ("a",)
        if self in {
            EOSCrystalSystem.TETRAGONAL,
            EOSCrystalSystem.TRIGONAL,
            EOSCrystalSystem.HEXAGONAL,
        }:
            return ("a", "c")
        return ("a", "b", "c")


def parse_eos_crystal_system(value: str | EOSCrystalSystem) -> EOSCrystalSystem:
    """Return a canonical EOS crystal system.

    Parameters
    ----------
    value : str or EOSCrystalSystem
        Crystal-system label. ``rhombohedral`` is accepted as an alias of
        ``trigonal``.

    Returns
    -------
    EOSCrystalSystem
        Canonical crystal system.

    Raises
    ------
    ValueError
        If ``value`` is not one of the seven crystal systems.
    """
    if isinstance(value, EOSCrystalSystem):
        return value
    normalized = str(value).strip().lower().replace("_", "-")
    aliases = {
        "rhombohedral": EOSCrystalSystem.TRIGONAL,
        "rhomboedral": EOSCrystalSystem.TRIGONAL,
    }
    if normalized in aliases:
        return aliases[normalized]
    try:
        return EOSCrystalSystem(normalized)
    except ValueError as exc:
        allowed = ", ".join(item.value for item in EOSCrystalSystem)
        raise ValueError(
            f"Unsupported EOS crystal system {value!r}; expected one of {allowed}."
        ) from exc


class EOSCoordinateVariation(str, Enum):
    """Numerical variation state of one tabulated EOS coordinate."""

    CONSTANT = "constant"
    VARIABLE = "variable"


@dataclass(frozen=True, slots=True)
class EOSCoordinateProfile:
    """Numerical range and variation state of one dataset coordinate.

    Parameters
    ----------
    name : str
        Canonical EOS column name.
    variation : EOSCoordinateVariation
        Whether the selected values are constant or variable within the
        declared numerical tolerance.
    minimum, maximum, span : float
        Range statistics in the stored normalized unit.
    reference_value : float or None
        Mean selected value when the coordinate is constant, otherwise
        ``None``.
    unit : str or None
        Normalized unit label.
    npoints : int
        Number of selected observations.
    absolute_tolerance, relative_tolerance : float
        Tolerances used for the constancy decision.
    """

    name: str
    variation: EOSCoordinateVariation
    minimum: float
    maximum: float
    span: float
    reference_value: float | None
    unit: str | None
    npoints: int
    absolute_tolerance: float
    relative_tolerance: float

    @property
    def is_constant(self) -> bool:
        """Return whether the selected coordinate is constant."""
        return self.variation is EOSCoordinateVariation.CONSTANT

    @property
    def is_variable(self) -> bool:
        """Return whether the selected coordinate contains usable variation."""
        return self.variation is EOSCoordinateVariation.VARIABLE

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready representation of the profile."""
        return {
            "name": self.name,
            "variation": self.variation.value,
            "minimum": self.minimum,
            "maximum": self.maximum,
            "span": self.span,
            "reference_value": self.reference_value,
            "unit": self.unit,
            "npoints": self.npoints,
            "absolute_tolerance": self.absolute_tolerance,
            "relative_tolerance": self.relative_tolerance,
        }


@dataclass(frozen=True, slots=True)
class EOSDatasetClassification:
    """Scientific classification of the coordinates available in a dataset.

    Parameters
    ----------
    profiles : dict
        Coordinate profiles keyed by canonical data-column name.
    variable_coordinates, constant_coordinates : tuple
        Coordinates classified from the selected observations.
    is_isobaric, is_isothermal : bool
        Whether pressure or temperature, respectively, is present and constant.
    reference_pressure, reference_temperature : float or None
        Constant reference conditions when available.

    Notes
    -----
    A constant pressure remains scientifically useful metadata for an isobaric
    V-T dataset, but it cannot constrain a P-V equation of state. Likewise, a
    constant temperature describes an isothermal condition without supplying
    information for a V-T fit.
    """

    profiles: dict[str, EOSCoordinateProfile]
    variable_coordinates: tuple[str, ...]
    constant_coordinates: tuple[str, ...]
    is_isobaric: bool
    is_isothermal: bool
    reference_pressure: float | None
    reference_temperature: float | None

    def profile(self, name: str) -> EOSCoordinateProfile:
        """Return the profile for one available canonical coordinate."""
        return self.profiles[name]

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready dataset classification."""
        return {
            "profiles": {
                name: profile.as_dict() for name, profile in self.profiles.items()
            },
            "variable_coordinates": list(self.variable_coordinates),
            "constant_coordinates": list(self.constant_coordinates),
            "is_isobaric": self.is_isobaric,
            "is_isothermal": self.is_isothermal,
            "reference_pressure": self.reference_pressure,
            "reference_temperature": self.reference_temperature,
        }


@dataclass(slots=True)
class EOSSeries:
    """One normalized dependent-variable series selected from an EOS dataset.

    Parameters
    ----------
    independent : str
        Canonical name of the independent variable.
    target : str
        Canonical name of the fitted dependent variable.
    x, y : ndarray
        Independent and dependent values as one-dimensional ``float64`` arrays.
    sigma_x, sigma_y : ndarray or None, optional
        Standard uncertainties associated with ``x`` and ``y``.
    mask : ndarray or None, optional
        Boolean selection mask. If omitted, every observation is selected.
    units : dict, optional
        Unit labels for the selected quantities.
    metadata : dict, optional
        Additional passive series metadata.

    Raises
    ------
    ValueError
        If arrays have incompatible shape, non-finite values, or an empty
        selection.
    """

    independent: str
    target: str
    x: np.ndarray
    y: np.ndarray
    sigma_x: np.ndarray | None = None
    sigma_y: np.ndarray | None = None
    mask: np.ndarray | None = None
    units: dict[str, str] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays and validate the selected series."""
        self.x = as_float64_vector(self.x, name="x")
        self.y = as_float64_vector(self.y, name="y")
        if self.x.shape != self.y.shape:
            raise ValueError("EOS series x and y arrays must have equal length.")
        self.sigma_x = optional_sigma(self.sigma_x, self.x.shape, "sigma_x")
        self.sigma_y = optional_sigma(self.sigma_y, self.y.shape, "sigma_y")
        if self.mask is None:
            self.mask = np.ones(self.x.shape, dtype=bool)
        else:
            mask = np.asarray(self.mask, dtype=bool)
            if mask.ndim != 1 or mask.shape != self.x.shape:
                raise ValueError("EOS series mask must match the data shape.")
            self.mask = mask.copy()
        if not np.any(self.mask):
            raise ValueError("EOS series mask must select at least one observation.")
        self.units = dict(self.units)
        self.metadata = dict(self.metadata)

    @property
    def size(self) -> int:
        """Return the total number of observations in the series."""
        return int(self.x.size)

    @property
    def selected_size(self) -> int:
        """Return the number of selected observations."""
        if self.mask is None:  # Defensive guard after post-init normalization.
            raise RuntimeError("EOS series mask was not initialized.")
        return int(np.count_nonzero(self.mask))

    def observations(self) -> FitObservations:
        """Return observations in the stored independent/target orientation."""
        return FitObservations(
            x=self.x,
            y=self.y,
            sigma_x=self.sigma_x,
            sigma_y=self.sigma_y,
            mask=self.mask,
            x_name=self.independent,
            y_name=self.target,
            x_unit=self.units.get(self.independent),
            y_unit=self.units.get(self.target),
            metadata=self.metadata,
        )

    def pressure_observations(self) -> FitObservations:
        """Return experimental data in the ``P(structural quantity)`` form.

        The input table naturally describes pressure as the independent
        experimental control and volume or a length as the selected target.
        Isothermal EOS fitting instead minimizes pressure residuals and thus
        evaluates pressure as a function of the structural quantity.

        Raises
        ------
        ValueError
            If the series is not pressure against a volume or linear target.
        """
        if self.independent != "pressure" or self.target not in {
            "volume",
            "a",
            "b",
            "c",
        }:
            raise ValueError(
                "pressure observations require pressure with volume or a linear target"
            )
        return FitObservations(
            x=self.y,
            y=self.x,
            sigma_x=self.sigma_y,
            sigma_y=self.sigma_x,
            mask=self.mask,
            x_name=self.target,
            y_name="pressure",
            x_unit=self.units.get(self.target),
            y_unit=self.units.get("pressure"),
            metadata={**self.metadata, "source_orientation": "pressure-control"},
        )

    def axial_pressure_observations(self) -> FitObservations:
        r"""Return a linear target as :math:`P(q)` with :math:`q=x^3`.

        The finite-strain linear-EOS convention cubes a cell-axis length and
        treats the transformed quantity as an auxiliary volume. Standard
        uncertainties are propagated to first order as

        .. math::

            \sigma_q = 3x^2\sigma_x.

        Returns
        -------
        FitObservations
            Pressure observations against the cubed length, suitable for OLS,
            WLS, effective variance, and ODR through the same solver contract.

        Raises
        ------
        ValueError
            If the series is not pressure against ``a``, ``b``, or ``c``, or
            if a length is not strictly positive.
        """
        pressure = self.pressure_observations()
        if self.target not in {"a", "b", "c"}:
            raise ValueError("axial observations require target 'a', 'b', or 'c'")
        length = np.asarray(pressure.x, dtype=np.float64)
        if np.any(length <= 0.0):
            raise ValueError("cell-axis lengths must be strictly positive")
        sigma_q = None
        if pressure.sigma_x is not None:
            sigma_q = 3.0 * length**2 * np.asarray(pressure.sigma_x, dtype=np.float64)
        unit = self.units.get(self.target)
        transformed_unit = None if unit is None else f"({unit})^3"
        return FitObservations(
            x=length**3,
            y=pressure.y,
            sigma_x=sigma_q,
            sigma_y=pressure.sigma_y,
            mask=pressure.mask,
            x_name=f"{self.target}^3",
            y_name="pressure",
            x_unit=transformed_unit,
            y_unit=self.units.get("pressure"),
            metadata={
                **pressure.metadata,
                "linear_target": self.target,
                "coordinate_transform": "q=x^3",
                "uncertainty_transform": "sigma_q=3*x^2*sigma_x",
                "original_length_unit": unit,
            },
        )

    def temperature_observations(self) -> FitObservations:
        """Return V--T or L--T data in the stored ``X(T)`` orientation.

        Returns
        -------
        FitObservations
            Temperature as the independent coordinate and the selected
            structural quantity as the dependent observation.

        Raises
        ------
        ValueError
            If the series is not temperature against volume or a cell axis.
        """
        if self.independent != "temperature" or self.target not in {
            "volume",
            "a",
            "b",
            "c",
        }:
            raise ValueError(
                "temperature observations require temperature with volume "
                "or a linear target"
            )
        return FitObservations(
            x=self.x,
            y=self.y,
            sigma_x=self.sigma_x,
            sigma_y=self.sigma_y,
            mask=self.mask,
            x_name="temperature",
            y_name=self.target,
            x_unit=self.units.get("temperature"),
            y_unit=self.units.get(self.target),
            metadata={**self.metadata, "source_orientation": "temperature-control"},
        )

    def axial_temperature_observations(self) -> FitObservations:
        r"""Return an axial V--T series as :math:`q(T)=x(T)^3`.

        Standard uncertainties are propagated to first order as

        .. math::

            \sigma_q = 3x^2\sigma_x.

        Returns
        -------
        FitObservations
            Temperature and auxiliary cubed-length observations.

        Raises
        ------
        ValueError
            If the target is not ``a``, ``b``, or ``c``, or if any selected
            length is non-positive.
        """
        observations = self.temperature_observations()
        if self.target not in {"a", "b", "c"}:
            raise ValueError("axial V-T observations require target a, b, or c")
        length = np.asarray(observations.y, dtype=np.float64)
        if np.any(length <= 0.0):
            raise ValueError("cell-axis lengths must be strictly positive")
        sigma_q = None
        if observations.sigma_y is not None:
            sigma_q = (
                3.0 * length**2 * np.asarray(observations.sigma_y, dtype=np.float64)
            )
        unit = self.units.get(self.target)
        transformed_unit = (
            "dimensionless"
            if unit == "dimensionless"
            else (None if unit is None else f"({unit})^3")
        )
        return FitObservations(
            x=observations.x,
            y=length**3,
            sigma_x=observations.sigma_x,
            sigma_y=sigma_q,
            mask=observations.mask,
            x_name="temperature",
            y_name=f"{self.target}^3",
            x_unit=self.units.get("temperature"),
            y_unit=transformed_unit,
            metadata={
                **observations.metadata,
                "linear_target": self.target,
                "coordinate_transform": "q=x^3",
                "uncertainty_transform": "sigma_q=3*x^2*sigma_x",
                "original_length_unit": unit,
            },
        )


@dataclass(slots=True)
class EOSDataset:
    """Normalized tabular data available to EOS workflows.

    Parameters
    ----------
    jobname : str, optional
        Dataset title or short description.
    columns : dict, optional
        Mapping from canonical EOS column names to one-dimensional arrays.
    units : dict, optional
        Mapping from canonical column names to normalized unit labels.
    raw_columns : dict, optional
        Original numeric columns before unit or scale normalization.
    raw_units : dict, optional
        Unit labels associated with ``raw_columns``.
    provenance : dict, optional
        Provenance information associated with individual columns.
    source : str, Path or None, optional
        Source file used to construct the dataset.
    metadata : dict, optional
        Additional dataset-level metadata.
    default_mask : ndarray or None, optional
        Input-level non-destructive selection. When omitted, all observations
        are selected.
    groups : ndarray or None, optional
        Positive integer data-group identifiers aligned with the observations.

    Raises
    ------
    ValueError
        If no data are supplied, a column name is unsupported, arrays have
        inconsistent lengths, values are non-finite, or an uncertainty is
        negative or has no associated measured quantity.
    """

    jobname: str = "Unknown"
    columns: dict[str, np.ndarray] = field(default_factory=dict)
    units: dict[str, str] = field(default_factory=dict)
    raw_columns: dict[str, np.ndarray] = field(default_factory=dict)
    raw_units: dict[str, str] = field(default_factory=dict)
    provenance: dict[str, str] = field(default_factory=dict)
    source: str | Path | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    default_mask: np.ndarray | None = None
    groups: np.ndarray | None = None

    def __post_init__(self) -> None:
        """Normalize columns to ``float64`` and validate the table."""
        if not self.columns:
            raise ValueError("EOS dataset must contain at least one data column.")
        normalized: dict[str, np.ndarray] = {}
        expected_size: int | None = None
        for name, values in self.columns.items():
            if name not in EOS_COLUMN_NAMES:
                raise ValueError(f"Unsupported EOS data column: {name}")
            array = as_float64_vector(values, name=name)
            if expected_size is None:
                expected_size = int(array.size)
            elif array.size != expected_size:
                raise ValueError("All EOS data columns must have equal length.")
            if name.startswith("sigma_") and np.any(array < 0.0):
                raise ValueError(f"EOS uncertainty column '{name}' cannot be negative.")
            normalized[name] = array
        for sigma_name, value_name in _SIGMA_TO_VALUE.items():
            if sigma_name in normalized and value_name not in normalized:
                raise ValueError(
                    f"EOS uncertainty column '{sigma_name}' requires '{value_name}'."
                )
        self.columns = normalized
        raw_normalized: dict[str, np.ndarray] = {}
        for name, values in self.raw_columns.items():
            if name not in EOS_COLUMN_NAMES:
                raise ValueError(f"Unsupported raw EOS data column: {name}")
            array = as_float64_vector(values, name=f"raw_{name}")
            if array.size != expected_size:
                raise ValueError(
                    "Raw and normalized EOS data columns must have equal length."
                )
            raw_normalized[name] = array
        if raw_normalized and set(raw_normalized) != set(normalized):
            raise ValueError(
                "Raw and normalized EOS datasets must contain the same columns."
            )
        self.raw_columns = raw_normalized
        self.units = {
            str(name): str(unit)
            for name, unit in self.units.items()
            if name in self.columns
        }
        self.raw_units = {
            str(name): str(unit)
            for name, unit in self.raw_units.items()
            if name in self.raw_columns
        }
        self.provenance = {
            str(name): str(value)
            for name, value in self.provenance.items()
            if name in self.columns
        }
        assert expected_size is not None
        if self.default_mask is None:
            self.default_mask = np.ones(expected_size, dtype=np.bool_)
        else:
            default_mask = np.asarray(self.default_mask, dtype=np.bool_)
            if default_mask.ndim != 1 or default_mask.size != expected_size:
                raise ValueError(
                    "EOS dataset default_mask must match the observation count."
                )
            if not np.any(default_mask):
                raise ValueError(
                    "EOS dataset default_mask must select at least one observation."
                )
            self.default_mask = default_mask.copy()
        if self.groups is not None:
            raw_groups = np.asarray(self.groups)
            if raw_groups.ndim != 1 or raw_groups.size != expected_size:
                raise ValueError("EOS dataset groups must match the observation count.")
            if not np.all(np.isfinite(raw_groups.astype(np.float64))):
                raise ValueError("EOS dataset groups must be finite integers.")
            group_values = raw_groups.astype(np.int64)
            if not np.all(raw_groups.astype(np.float64) == group_values):
                raise ValueError("EOS dataset groups must contain integers.")
            if np.any(group_values < 1):
                raise ValueError("EOS dataset group identifiers must be positive.")
            self.groups = group_values.copy()
        self.metadata = dict(self.metadata)
        if isinstance(self.source, str):
            self.source = Path(self.source)

    @property
    def npoints(self) -> int:
        """Return the number of tabulated observations."""
        first = next(iter(self.columns.values()))
        return int(first.size)

    @property
    def selected_npoints(self) -> int:
        """Return observations selected by the dataset default mask."""
        if self.default_mask is None:  # Defensive guard after normalization.
            raise RuntimeError("EOS dataset default mask was not initialized.")
        return int(np.count_nonzero(self.default_mask))

    @property
    def excluded_npoints(self) -> int:
        """Return observations excluded by the dataset default mask."""
        return self.npoints - self.selected_npoints

    @property
    def group_ids(self) -> tuple[int, ...]:
        """Return sorted group identifiers available in the dataset."""
        if self.groups is None:
            return ()
        return tuple(int(value) for value in np.unique(self.groups))

    @property
    def crystal_system(self) -> EOSCrystalSystem | None:
        """Return the canonical crystal system stored in dataset metadata."""
        value = self.metadata.get("crystal_system")
        return None if value is None else parse_eos_crystal_system(str(value))

    def selection_mask(self, mask: np.ndarray | None = None) -> np.ndarray:
        """Return a validated effective observation mask.

        Parameters
        ----------
        mask : ndarray or None, optional
            Explicit final selection.  When omitted, the non-destructive
            dataset default mask produced from ``USE`` or trailing ``*``
            markers is returned.

        Returns
        -------
        ndarray
            Detached one-dimensional boolean selection.
        """
        if mask is None:
            if self.default_mask is None:  # Defensive guard.
                raise RuntimeError("EOS dataset default mask was not initialized.")
            return self.default_mask.copy()
        selected = np.asarray(mask, dtype=np.bool_)
        if selected.ndim != 1 or selected.size != self.npoints:
            raise ValueError("EOS selection mask must match the observation count.")
        if not np.any(selected):
            raise ValueError("EOS selection mask must select at least one observation.")
        return selected.copy()

    def group_summary(
        self, mask: np.ndarray | None = None
    ) -> tuple[dict[str, int], ...]:
        """Return total, selected, and excluded counts for each data group."""
        selected = self.selection_mask(mask)
        groups = self.groups
        if groups is None:
            return ()
        rows: list[dict[str, int]] = []
        for group_id in self.group_ids:
            member = groups == group_id
            count = int(np.count_nonzero(member))
            used = int(np.count_nonzero(member & selected))
            rows.append(
                {
                    "group": group_id,
                    "total": count,
                    "selected": used,
                    "excluded": count - used,
                }
            )
        return tuple(rows)

    @property
    def column_names(self) -> tuple[str, ...]:
        """Return canonical column names in input order."""
        return tuple(self.columns)

    def has(self, name: str) -> bool:
        """Return whether a canonical column is present."""
        return name in self.columns

    def column(self, name: str) -> np.ndarray:
        """Return one canonical data column.

        Parameters
        ----------
        name : str
            Canonical EOS column name.

        Returns
        -------
        ndarray
            Stored one-dimensional ``float64`` array.

        Raises
        ------
        KeyError
            If the requested column is unavailable.
        """
        return self.columns[name]

    def coordinate_profile(
        self,
        name: str,
        *,
        mask: np.ndarray | None = None,
        relative_tolerance: float = 1.0e-12,
        absolute_tolerance: float | None = None,
    ) -> EOSCoordinateProfile:
        """Classify one coordinate as constant or variable.

        Parameters
        ----------
        name : str
            Canonical non-uncertainty column name.
        mask : ndarray or None, optional
            Boolean observation selection.
        relative_tolerance, absolute_tolerance : float or None, optional
            Numerical constancy tolerances.

        Returns
        -------
        EOSCoordinateProfile
            Range, reference condition, and variation state.
        """
        from .classification import coordinate_profile

        return coordinate_profile(
            self,
            name,
            mask=self.selection_mask(mask),
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )

    def classify(
        self,
        *,
        mask: np.ndarray | None = None,
        relative_tolerance: float = 1.0e-12,
        absolute_tolerance: float | None = None,
    ) -> EOSDatasetClassification:
        """Classify all measured coordinates and reference conditions.

        Parameters
        ----------
        mask : ndarray or None, optional
            Boolean observation selection.
        relative_tolerance, absolute_tolerance : float or None, optional
            Numerical constancy tolerances.

        Returns
        -------
        EOSDatasetClassification
            Coordinate profiles and isobaric/isothermal conditions.
        """
        from .classification import classify_dataset

        return classify_dataset(
            self,
            mask=self.selection_mask(mask),
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )

    def require_variable_coordinate(
        self,
        name: str,
        *,
        purpose: str,
        mask: np.ndarray | None = None,
    ) -> EOSCoordinateProfile:
        """Return a profile or reject a constant fitting coordinate.

        Parameters
        ----------
        name : str
            Canonical measured coordinate.
        purpose : str
            Human-readable scientific operation used in the error message.
        mask : ndarray or None, optional
            Boolean observation selection.

        Returns
        -------
        EOSCoordinateProfile
            Variable coordinate profile.
        """
        from .classification import require_variable_coordinate

        return require_variable_coordinate(
            self, name, purpose=purpose, mask=self.selection_mask(mask)
        )

    def series(
        self,
        target: str,
        *,
        independent: str = "pressure",
        mask: np.ndarray | None = None,
    ) -> EOSSeries:
        """Select one fitting series from the dataset.

        Parameters
        ----------
        target : str
            Dependent quantity, such as ``"volume"`` or ``"a"``.
        independent : str, optional
            Independent quantity. The experimental P-V workflow uses
            ``"pressure"``; future V-T and E-V workflows can select another
            available column explicitly.
        mask : ndarray or None, optional
            Boolean observation mask.

        Returns
        -------
        EOSSeries
            Normalized selected series with available uncertainties.

        Raises
        ------
        ValueError
            If a name is invalid or both names identify the same quantity.
        KeyError
            If a required data column is unavailable.
        """
        if target not in EOS_TARGET_NAMES:
            raise ValueError(f"Unsupported EOS fitting target: {target}")
        if independent not in {"pressure", "temperature", "volume"}:
            raise ValueError(f"Unsupported EOS independent variable: {independent}")
        if target == independent:
            raise ValueError("EOS target and independent variable must differ.")
        sigma_x_name = f"sigma_{independent}"
        sigma_y_name = f"sigma_{target}"
        selected_units = {
            name: self.units[name]
            for name in (independent, target, sigma_x_name, sigma_y_name)
            if name in self.units
        }
        column_scales = self.metadata.get("column_scales", {})
        scale_metadata = dict(column_scales) if isinstance(column_scales, dict) else {}
        return EOSSeries(
            independent=independent,
            target=target,
            x=self.columns[independent],
            y=self.columns[target],
            sigma_x=self.columns.get(sigma_x_name),
            sigma_y=self.columns.get(sigma_y_name),
            mask=self.selection_mask(mask),
            units=selected_units,
            metadata={
                "jobname": self.jobname,
                "independent_scale": scale_metadata.get(independent, "absolute"),
                "target_scale": scale_metadata.get(target, "absolute"),
            },
        )

    def pvt_observations(
        self,
        *,
        mask: np.ndarray | None = None,
    ) -> FitObservations:
        r"""Return pressure observations at paired volume and temperature.

        The explanatory coordinate matrix has rows ``volume`` and
        ``temperature``.  Missing explanatory-coordinate uncertainties are
        represented by zeros when the other coordinate uncertainty is
        available.  This allows effective variance to include the available
        terms while ODR still rejects any selected zero uncertainty.

        Parameters
        ----------
        mask : ndarray or None, optional
            Boolean observation selection.

        Returns
        -------
        FitObservations
            Multicoordinate observations for :math:`P(V,T)`.

        Raises
        ------
        KeyError
            If pressure, volume, or temperature is absent.
        """
        pressure = self.columns["pressure"]
        volume = self.columns["volume"]
        temperature = self.columns["temperature"]
        sigma_volume = self.columns.get("sigma_volume")
        sigma_temperature = self.columns.get("sigma_temperature")
        sigma_x = None
        if sigma_volume is not None or sigma_temperature is not None:
            sigma_x = np.vstack(
                [
                    np.zeros_like(volume) if sigma_volume is None else sigma_volume,
                    np.zeros_like(temperature)
                    if sigma_temperature is None
                    else sigma_temperature,
                ]
            )
        return FitObservations(
            x=np.vstack([volume, temperature]),
            y=pressure,
            sigma_x=sigma_x,
            sigma_y=self.columns.get("sigma_pressure"),
            mask=self.selection_mask(mask),
            x_name="volume,temperature",
            y_name="pressure",
            x_unit=f"{self.units.get('volume', 'angstrom^3')},{self.units.get('temperature', 'K')}",
            y_unit=self.units.get("pressure", "GPa"),
            metadata={
                "jobname": self.jobname,
                "coordinate_order": ["volume", "temperature"],
                "coordinate_units": [
                    self.units.get("volume", "angstrom^3"),
                    self.units.get("temperature", "K"),
                ],
            },
        )

__all__ = [
    "EOS_COLUMN_NAMES",
    "EOS_TARGET_NAMES",
    "EOSCoordinateProfile",
    "EOSCoordinateVariation",
    "EOSCrystalSystem",
    "EOSDataset",
    "EOSDatasetClassification",
    "EOSSeries",
    "parse_eos_crystal_system",
]
