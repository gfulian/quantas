# -*- coding: utf-8 -*-

"""Build backend-neutral raw elastic volume series from VASP outputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.physics.elasticity.pressure_resolution import (
    EnergyPressureResolution,
    resolve_energy_derived_pressures,
)
from quantas.core.physics.elasticity.prestress import assign_hydrostatic_pressures
from quantas.models.elastic_states import (
    ElasticState,
    ElasticStateSeries,
    ElasticTensorKind,
    PressureSource,
    PrestressProvenance,
)

from .elasticity import VASPElasticityReader


def read_vasp_elastic_series(
    sources: Sequence[str | Path],
    *,
    reference_source_index: int = 0,
) -> ElasticStateSeries:
    """Read a raw volume-resolved elastic series from VASP calculations.

    This adapter performs only source normalization. It preserves the selected
    VASP elastic tensor, reference stress/pressure provenance, density, and
    primitive-cell volume reported by :class:`VASPElasticityReader`. No
    finite-prestress correction is applied here.

    Parameters
    ----------
    sources : sequence of str or pathlib.Path
        VASP run directories, ``OUTCAR`` files, or ``vasprun.xml`` files with
        sibling ``OUTCAR`` files.
    reference_source_index : int, optional
        Zero-based reference index in input-source order. The reference follows
        the same physical state after sorting by increasing volume.

    Returns
    -------
    ElasticStateSeries
        Increasing-volume backend-neutral series containing the raw VASP
        stiffness tensors and their source pressure provenance.

    Raises
    ------
    ValueError
        If the source list, reference index, parsed elasticity data, volume, or
        density is invalid or ambiguous.
    """
    paths = tuple(Path(source) for source in sources)
    if not paths:
        raise ValueError("at least one VASP elastic source is required")
    if len(set(paths)) != len(paths):
        raise ValueError("VASP elastic sources must be unique")
    reference_source_index = int(reference_source_index)
    if not 0 <= reference_source_index < len(paths):
        raise ValueError("reference_source_index is outside the VASP source list")

    states: list[ElasticState] = []
    for source_index, path in enumerate(paths):
        reader = VASPElasticityReader(path)
        if not reader.completed:
            detail = reader.error or "unknown reader error"
            raise ValueError(f"unable to read VASP elastic source {path}: {detail}")
        volume = reader.reference_volume_angstrom3
        density = reader.density
        if volume is None or not np.isfinite(volume) or volume <= 0.0:
            raise ValueError(f"VASP elastic source {path} lacks a valid reference volume")
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError(f"VASP elastic source {path} lacks a valid density")

        metadata = dict(reader.metadata)
        metadata.update(
            {
                "backend": "vasp",
                "source_index": source_index,
                "series_adapter": "read_vasp_elastic_series",
                "prestress_correction_applied": False,
            }
        )
        states.append(
            ElasticState(
                volume=volume,
                density=density,
                stiffness=reader.stiffness,
                prestress=reader.prestress,
                source=path,
                metadata=metadata,
            )
        )

    states.sort(key=lambda state: state.volume)
    volumes = np.asarray([state.volume for state in states], dtype=np.float64)
    if np.any(np.diff(volumes) <= 0.0):
        raise ValueError("VASP elastic sources contain duplicate volumes")

    reference_index = next(
        index
        for index, state in enumerate(states)
        if int(state.metadata["source_index"]) == reference_source_index
    )
    return ElasticStateSeries(
        states=tuple(states),
        reference_index=reference_index,
        orientation="vasp-cartesian",
        metadata={
            "backend": "vasp",
            "reader": "VASPElasticityReader",
            "reference_policy": "input_source_index",
            "reference_source_index": reference_source_index,
            "input_source_count": len(paths),
            "prestress_correction_applied": False,
        },
    )


def assign_vasp_manual_pressures(
    series: ElasticStateSeries,
    pressures_gpa: ArrayLike,
    *,
    assignment_method: str = "manual",
    metadata: Mapping[str, object] | None = None,
) -> ElasticStateSeries:
    """Replace raw VASP output pressures with explicit hydrostatic values.

    Parameters
    ----------
    series : ElasticStateSeries
        Raw VASP elastic series produced by :func:`read_vasp_elastic_series`.
    pressures_gpa : array_like
        One hydrostatic pressure per elastic state, in GPa and positive in
        compression.
    assignment_method : str, optional
        Stable provenance label for the pressure assignment.
    metadata : mapping, optional
        Additional series-level provenance.

    Returns
    -------
    ElasticStateSeries
        Independent raw VASP series with replacement pressure provenance.

    Raises
    ------
    TypeError
        If ``series`` is not an :class:`ElasticStateSeries`.
    ValueError
        If the series is not raw VASP stress--strain data or the pressure
        assignment is invalid.

    Notes
    -----
    This operation changes pressure provenance only. Stiffness coefficients
    remain raw until :func:`convert_vasp_hydrostatic_elastic_series` is called.
    The original VASP output pressure remains recorded in each state's
    ``pressure_assignment`` metadata.
    """
    _require_raw_vasp_series(series)
    return assign_hydrostatic_pressures(
        series,
        pressures_gpa,
        pressure_source=PressureSource.MANUAL,
        assignment_method=assignment_method,
        metadata=metadata,
        replace_existing=True,
    )


def resolve_vasp_energy_derived_pressures(
    series: ElasticStateSeries,
    energy_volumes: ArrayLike,
    energies: ArrayLike,
    *,
    pressure_source: PressureSource | str,
    source_dataset: str,
    energy_unit: str,
    volume_length_unit: str,
    volume_unit: str | None = None,
    eos: str = "BM3",
    polynomial_degree: int = 3,
    maxfev: int | None = None,
) -> EnergyPressureResolution:
    """Assign EOS- or polynomial-derived pressures to raw VASP elastic states.

    Parameters
    ----------
    series : ElasticStateSeries
        Raw VASP stress--strain series with output-stress provenance.
    energy_volumes, energies : array_like
        Static energy-volume samples normalized to the same primitive-cell
        basis as ``series``.
    pressure_source : PressureSource or str
        ``energy_eos`` or ``energy_polynomial``.
    source_dataset : str
        Provenance label identifying the static energy-volume dataset.
    energy_unit : str
        Unit of ``energies``.
    volume_length_unit : str
        Length unit whose cube defines ``energy_volumes``.
    volume_unit : str or None, optional
        Human-readable volume unit retained in provenance.
    eos : str, optional
        Integrated energy EOS used for ``energy_eos``.
    polynomial_degree : int, optional
        Polynomial degree used for ``energy_polynomial``.
    maxfev : int or None, optional
        Optional EOS fitter evaluation limit.

    Returns
    -------
    EnergyPressureResolution
        Raw VASP series carrying energy-derived pressure provenance together
        with fit diagnostics and explicit volume matches.

    Raises
    ------
    TypeError
        If ``series`` has an unsupported type.
    ValueError
        If the series is not raw VASP data or the energy-derived pressure
        resolution fails.

    Notes
    -----
    The VASP stiffness tensors remain unchanged. The output-stress pressure is
    retained as replaced provenance, while the selected E(V)-derived pressure
    becomes the value used by the subsequent VASP hydrostatic conversion.
    """
    _require_raw_vasp_series(series)
    return resolve_energy_derived_pressures(
        series,
        energy_volumes,
        energies,
        pressure_source=pressure_source,
        source_dataset=source_dataset,
        energy_unit=energy_unit,
        volume_length_unit=volume_length_unit,
        volume_unit=volume_unit,
        eos=eos,
        polynomial_degree=polynomial_degree,
        maxfev=maxfev,
        replace_existing_pressure=True,
    )


VASP_HYDROSTATIC_PRESTRESS_METHOD = "vasp-residual-pressure-hydrostatic"
VASP_HYDROSTATIC_PRESTRESS_REFERENCE_DOI = "10.1016/j.cpc.2021.108068"


def vasp_hydrostatic_incremental_stiffness(
    raw_stiffness: ArrayLike,
    pressure_gpa: float,
) -> NDArray[np.float64]:
    r"""Apply the hydrostatic pressure adjustment required for VASP stiffnesses.

    VASP ``TOTAL ELASTIC MODULI`` are treated as raw finite-difference
    stress--strain coefficients.  For hydrostatic pressure ``P`` (positive in
    compression), the VASP post-processing relation documented by Singh et al.
    (Computer Physics Communications 267, 108068, 2021, Appendix A) is

    .. math::

       B_{aa} = C_{aa} - P, \quad a=1,\ldots,6,

       B_{12} = C_{12} + P, \quad
       B_{13} = C_{13} + P, \quad
       B_{23} = C_{23} + P,

    with symmetric partners adjusted identically and all other off-diagonal
    coefficients unchanged.  The six-index diagonal includes the three shear
    entries ``44``, ``55`` and ``66``.

    Parameters
    ----------
    raw_stiffness : array_like
        Symmetric VASP stiffness matrix in GPa and Quantas Voigt order
        ``(11, 22, 33, 23, 13, 12)``.
    pressure_gpa : float
        Hydrostatic pressure in GPa, positive in compression.

    Returns
    -------
    ndarray
        Hydrostatic incremental stiffness matrix in GPa.

    Raises
    ------
    ValueError
        If the stiffness matrix or pressure is invalid.

    Notes
    -----
    This is a VASP ingestion rule, not the CRYSTAL/Erba energy--strain
    transformation and not Quantas' generic Eulerian finite-strain operator.
    """
    stiffness = np.asarray(raw_stiffness, dtype=np.float64)
    pressure = float(pressure_gpa)
    if stiffness.shape != (6, 6) or not np.all(np.isfinite(stiffness)):
        raise ValueError("raw_stiffness must be finite with shape (6, 6)")
    if not np.allclose(stiffness, stiffness.T, rtol=0.0, atol=1.0e-10):
        raise ValueError("raw_stiffness must be symmetric")
    if not np.isfinite(pressure):
        raise ValueError("pressure_gpa must be finite")

    corrected = stiffness.copy()
    diagonal = np.diag_indices(6)
    corrected[diagonal] -= pressure
    for first, second in ((0, 1), (0, 2), (1, 2)):
        corrected[first, second] += pressure
        corrected[second, first] += pressure
    return np.asarray(0.5 * (corrected + corrected.T), dtype=np.float64)


def convert_vasp_hydrostatic_elastic_state(
    state: ElasticState,
    *,
    correction_applied_by: str = "quantas",
    hydrostatic_atol_gpa: float = 1.0e-2,
) -> ElasticState:
    """Convert one raw VASP state to hydrostatic incremental stiffness.

    Parameters
    ----------
    state : ElasticState
        Raw VASP stress--strain state with explicit hydrostatic pressure
        provenance and the original reference stress tensor retained in metadata.
    correction_applied_by : str, optional
        Provenance label for the caller applying the conversion.
    hydrostatic_atol_gpa : float, optional
        Maximum absolute component residual, in GPa, allowed between the
        original VASP reference stress and its own hydrostatic projection.

    Returns
    -------
    ElasticState
        Independent state containing pressure-adjusted incremental stiffness.

    Raises
    ------
    TypeError
        If ``state`` is not an :class:`ElasticState`.
    ValueError
        If tensor kind, pressure provenance, reference stress, tolerance, or
        correction provenance is incompatible with this VASP-only conversion.
    """
    if not isinstance(state, ElasticState):
        raise TypeError("state must be an ElasticState")
    source_kind = ElasticTensorKind(state.prestress.tensor_kind)
    if source_kind is not ElasticTensorKind.RAW_STRESS_STRAIN:
        raise ValueError(
            "VASP hydrostatic conversion requires a raw stress-strain tensor"
        )
    pressure = state.prestress.pressure_gpa
    pressure_source = PressureSource(state.prestress.pressure_source)
    if pressure is None or pressure_source in {
        PressureSource.UNAVAILABLE,
        PressureSource.APPLIED_PRESTRESS,
    }:
        raise ValueError(
            "VASP hydrostatic conversion requires explicit hydrostatic pressure provenance"
        )
    tolerance = float(hydrostatic_atol_gpa)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("hydrostatic_atol_gpa must be finite and non-negative")
    applied_by = str(correction_applied_by).strip()
    if not applied_by:
        raise ValueError("correction_applied_by must be non-empty")

    output_stress_pressure, hydrostatic_residual = _reference_hydrostaticity(
        state,
        selected_pressure_gpa=float(pressure),
        pressure_source=pressure_source,
        tolerance_gpa=tolerance,
    )

    metadata = dict(state.metadata)
    metadata["prestress_correction_applied"] = True
    metadata["prestress_correction"] = {
        "method": VASP_HYDROSTATIC_PRESTRESS_METHOD,
        "reference_doi": VASP_HYDROSTATIC_PRESTRESS_REFERENCE_DOI,
        "pressure_gpa": float(pressure),
        "pressure_source": pressure_source.value,
        "output_stress_pressure_gpa": output_stress_pressure,
        "pressure_minus_output_stress_gpa": float(pressure) - output_stress_pressure,
        "hydrostatic_stress_residual_gpa": hydrostatic_residual,
        "hydrostatic_atol_gpa": tolerance,
        "applied_by": applied_by,
        "source_tensor_kind": source_kind.value,
        "target_tensor_kind": ElasticTensorKind.WALLACE_HYDROSTATIC.value,
    }
    return ElasticState(
        volume=state.volume,
        density=state.density,
        stiffness=vasp_hydrostatic_incremental_stiffness(
            state.stiffness,
            float(pressure),
        ),
        prestress=PrestressProvenance(
            tensor_kind=ElasticTensorKind.WALLACE_HYDROSTATIC,
            pressure_gpa=float(pressure),
            pressure_source=pressure_source,
            correction_method=VASP_HYDROSTATIC_PRESTRESS_METHOD,
            correction_applied_by=applied_by,
            source_tensor_kind=source_kind,
        ),
        energy=state.energy,
        energy_unit=state.energy_unit,
        lattice=state.lattice,
        source=state.source,
        metadata=metadata,
    )


def convert_vasp_hydrostatic_elastic_series(
    series: ElasticStateSeries,
    *,
    correction_applied_by: str = "quantas",
    hydrostatic_atol_gpa: float = 1.0e-2,
) -> ElasticStateSeries:
    """Convert a raw VASP elastic series using its selected pressures.

    Parameters
    ----------
    series : ElasticStateSeries
        Increasing series produced by :func:`read_vasp_elastic_series`.
    correction_applied_by : str, optional
        Provenance label recorded on every corrected state.
    hydrostatic_atol_gpa : float, optional
        Hydrostatic-stress residual tolerance passed to every state conversion.

    Returns
    -------
    ElasticStateSeries
        Independent series of hydrostatic incremental stiffness tensors.

    Raises
    ------
    TypeError
        If ``series`` is not an :class:`ElasticStateSeries`.
    ValueError
        If any state cannot be converted exactly once from raw VASP output.
    """
    if not isinstance(series, ElasticStateSeries):
        raise TypeError("series must be an ElasticStateSeries")
    corrected_states: list[ElasticState] = []
    for index, state in enumerate(series.states):
        try:
            corrected_states.append(
                convert_vasp_hydrostatic_elastic_state(
                    state,
                    correction_applied_by=correction_applied_by,
                    hydrostatic_atol_gpa=hydrostatic_atol_gpa,
                )
            )
        except ValueError as exc:
            raise ValueError(f"elastic state {index}: {exc}") from exc

    metadata = dict(series.metadata)
    metadata["prestress_correction_applied"] = True
    metadata["prestress_correction"] = {
        "method": VASP_HYDROSTATIC_PRESTRESS_METHOD,
        "reference_doi": VASP_HYDROSTATIC_PRESTRESS_REFERENCE_DOI,
        "applied_by": str(correction_applied_by).strip(),
        "hydrostatic_atol_gpa": float(hydrostatic_atol_gpa),
        "state_count": series.nstates,
    }
    return ElasticStateSeries(
        states=tuple(corrected_states),
        reference_index=series.reference_index,
        orientation=series.orientation,
        metadata=metadata,
    )


def _reference_hydrostaticity(
    state: ElasticState,
    *,
    selected_pressure_gpa: float,
    pressure_source: PressureSource,
    tolerance_gpa: float,
) -> tuple[float, float]:
    """Return output pressure and residual after validating VASP hydrostaticity."""
    stress_value = state.metadata.get("reference_stress_gpa")
    if stress_value is None:
        raise ValueError(
            "VASP hydrostatic conversion requires the full reference stress tensor"
        )
    stress = np.asarray(stress_value, dtype=np.float64)
    if stress.shape != (3, 3) or not np.all(np.isfinite(stress)):
        raise ValueError("reference_stress_gpa must be finite with shape (3, 3)")
    if not np.allclose(stress, stress.T, rtol=0.0, atol=1.0e-10):
        raise ValueError("reference_stress_gpa must be symmetric")

    output_pressure = float(np.trace(stress) / 3.0)
    if pressure_source is PressureSource.OUTPUT_STRESS and not np.isclose(
        selected_pressure_gpa,
        output_pressure,
        rtol=0.0,
        atol=tolerance_gpa,
    ):
        raise ValueError(
            "VASP output-stress pressure disagrees with the retained reference "
            "stress tensor"
        )
    target = output_pressure * np.eye(3, dtype=np.float64)
    residual = float(np.max(np.abs(stress - target)))
    if residual > tolerance_gpa:
        raise ValueError(
            "VASP reference stress is not hydrostatic within the requested "
            f"tolerance: residual={residual:.6g} GPa, "
            f"tolerance={tolerance_gpa:.6g} GPa"
        )
    return output_pressure, residual


def _require_raw_vasp_series(series: ElasticStateSeries) -> None:
    """Require a VASP raw stress--strain series before pressure replacement."""
    if not isinstance(series, ElasticStateSeries):
        raise TypeError("series must be an ElasticStateSeries")
    if str(series.metadata.get("backend", "")).lower() != "vasp":
        raise ValueError("pressure replacement requires backend='vasp'")
    for index, state in enumerate(series.states):
        if (
            ElasticTensorKind(state.prestress.tensor_kind)
            is not ElasticTensorKind.RAW_STRESS_STRAIN
        ):
            raise ValueError(
                f"elastic state {index}: pressure replacement requires a raw "
                "VASP stress-strain tensor"
            )


__all__ = [
    "VASP_HYDROSTATIC_PRESTRESS_METHOD",
    "VASP_HYDROSTATIC_PRESTRESS_REFERENCE_DOI",
    "convert_vasp_hydrostatic_elastic_series",
    "convert_vasp_hydrostatic_elastic_state",
    "assign_vasp_manual_pressures",
    "resolve_vasp_energy_derived_pressures",
    "read_vasp_elastic_series",
    "vasp_hydrostatic_incremental_stiffness",
]
