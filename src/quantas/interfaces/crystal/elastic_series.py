# -*- coding: utf-8 -*-

"""Build backend-neutral elastic volume series from CRYSTAL outputs."""

from __future__ import annotations

from collections.abc import Sequence
from enum import Enum
from pathlib import Path

import numpy as np

from quantas.core.physics.elasticity import correct_hydrostatic_elastic_state
from quantas.models.elastic_states import (
    ElasticState,
    ElasticStateSeries,
    ElasticTensorKind,
    PressureSource,
    PrestressProvenance,
)

from .elasticity import CrystalElasticityReader


class CrystalPressurePolicy(str, Enum):
    """Select the pressure attached to raw CRYSTAL elastic tensors."""

    AUTO = "auto"
    OUTPUT_STRESS = "output_stress"
    MANUAL = "manual"
    DEFERRED = "deferred"


def read_crystal_elastic_series(
    filenames: Sequence[str | Path],
    *,
    pressure_policy: CrystalPressurePolicy | str = CrystalPressurePolicy.AUTO,
    manual_pressures_gpa: Sequence[float] | None = None,
    apply_prestress_correction: bool = True,
    correction_applied_by: str = "quantas-crystal-import",
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
) -> ElasticStateSeries:
    """Read a volume-resolved elastic series from CRYSTAL output files.

    The returned states are sorted by increasing primitive-cell volume.  The
    reference state is the one with the lowest static energy.  Manual
    pressures, when requested, correspond to ``filenames`` before sorting.

    Parameters
    ----------
    filenames : sequence of str or Path
        Completed CRYSTAL ELASTCON or ELAPIEZO output files.
    pressure_policy : CrystalPressurePolicy or str, optional
        ``"auto"`` accepts a CRYSTAL-corrected tensor when ``PRESSURE`` or
        ``PRESSEOS`` is present and otherwise uses pressure from the unstrained
        output stress. ``"output_stress"`` requires raw tensors and uses that
        stress explicitly. ``"manual"`` requires ``manual_pressures_gpa``.
        ``"deferred"`` retains an explicitly raw tensor without assigning a
        pressure so that another workflow can attach independently fitted
        values before correction.
    manual_pressures_gpa : sequence of float or None, optional
        Hydrostatic pressures in GPa, positive in compression, in input-file
        order.  Accepted only with the manual policy.
    apply_prestress_correction : bool, optional
        Convert raw energy--strain tensors to Wallace hydrostatic tensors.
        Tensors already corrected by CRYSTAL are retained unchanged.
    correction_applied_by : str, optional
        Provenance label used when Quantas applies the Wallace correction.
    symprec, angle_tolerance : float, optional
        Symmetry tolerances forwarded to :class:`CrystalElasticityReader`.

    Returns
    -------
    ElasticStateSeries
        Increasing-volume series with complete source and pressure provenance.

    Raises
    ------
    ValueError
        If the file list, pressure policy, parsed data, volumes, or energies
        are invalid or ambiguous.
    """
    paths = tuple(Path(filename) for filename in filenames)
    if not paths:
        raise ValueError("at least one CRYSTAL elastic output is required")
    if len(set(paths)) != len(paths):
        raise ValueError("CRYSTAL elastic output paths must be unique")

    policy = CrystalPressurePolicy(pressure_policy)
    if policy is CrystalPressurePolicy.DEFERRED and apply_prestress_correction:
        raise ValueError(
            "pressure_policy='deferred' requires apply_prestress_correction=False"
        )
    manual = _manual_pressures(policy, manual_pressures_gpa, len(paths))
    states: list[ElasticState] = []
    for index, path in enumerate(paths):
        reader = CrystalElasticityReader(
            path,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
        )
        if not reader.completed:
            detail = reader.error or "unknown reader error"
            raise ValueError(f"unable to read CRYSTAL elastic output {path}: {detail}")
        _require_state_scalars(reader, path)
        state = _state_from_reader(reader, path, policy, manual[index])
        tensor_kind = ElasticTensorKind(state.prestress.tensor_kind)
        if apply_prestress_correction and not tensor_kind.is_incremental:
            state = correct_hydrostatic_elastic_state(
                state,
                correction_applied_by=correction_applied_by,
            )
        states.append(state)

    states.sort(key=lambda state: state.volume)
    volumes = np.asarray([state.volume for state in states], dtype=np.float64)
    if np.any(np.diff(volumes) <= 0.0):
        raise ValueError("CRYSTAL elastic outputs contain duplicate volumes")
    energies = np.asarray(
        [state.energy if state.energy is not None else np.nan for state in states],
        dtype=np.float64,
    )
    reference_index = int(np.argmin(energies))
    return ElasticStateSeries(
        states=tuple(states),
        reference_index=reference_index,
        orientation="crystal-cartesian",
        metadata={
            "backend": "crystal",
            "reader": "CrystalElasticityReader",
            "pressure_policy": policy.value,
            "prestress_correction_requested": bool(apply_prestress_correction),
            "reference_policy": "minimum_static_energy",
            "input_file_count": len(paths),
        },
    )


def _manual_pressures(
    policy: CrystalPressurePolicy,
    values: Sequence[float] | None,
    count: int,
) -> tuple[float | None, ...]:
    """Validate manual pressures and return one optional value per file."""
    if policy is not CrystalPressurePolicy.MANUAL:
        if values is not None:
            raise ValueError("manual_pressures_gpa requires pressure_policy='manual'")
        return (None,) * count
    if values is None or len(values) != count:
        raise ValueError(
            "manual pressure policy requires one pressure per CRYSTAL output"
        )
    pressures = tuple(float(value) for value in values)
    if not np.all(np.isfinite(pressures)):
        raise ValueError("manual pressures must be finite")
    return pressures


def _require_state_scalars(reader: CrystalElasticityReader, path: Path) -> None:
    """Require the volume, density, and energy needed by a volume series."""
    required = {
        "volume": reader.volume,
        "density": reader.density,
        "energy": reader.energy,
    }
    invalid = [name for name, value in required.items() if not np.isfinite(value)]
    if invalid:
        fields = ", ".join(invalid)
        raise ValueError(f"CRYSTAL elastic output {path} lacks finite {fields}")
    if reader.volume <= 0.0 or reader.density <= 0.0:
        raise ValueError(f"CRYSTAL elastic output {path} has invalid volume or density")


def _state_from_reader(
    reader: CrystalElasticityReader,
    path: Path,
    policy: CrystalPressurePolicy,
    manual_pressure: float | None,
) -> ElasticState:
    """Translate one completed reader into a backend-neutral state."""
    if reader.prestress_applied:
        if policy is not CrystalPressurePolicy.AUTO:
            raise ValueError(
                f"CRYSTAL output {path} already contains a {reader.prestress_keyword or 'pre-stress'} correction; "
                "use pressure_policy='auto' to preserve it"
            )
        pressure = (
            reader.pressure
            if np.isfinite(reader.pressure)
            else reader.pressure_keyword_value
        )
        if not np.isfinite(pressure):
            raise ValueError(
                f"CRYSTAL output {path} reports a pre-stress correction without a value"
            )
        prestress = PrestressProvenance(
            tensor_kind=ElasticTensorKind.WALLACE_HYDROSTATIC,
            pressure_gpa=pressure,
            pressure_source=PressureSource.APPLIED_PRESTRESS,
            correction_method=(
                f"crystal-{(reader.prestress_keyword or 'pressure').lower()}-keyword"
            ),
            correction_applied_by="crystal",
            source_tensor_kind=ElasticTensorKind.RAW_ENERGY_STRAIN,
        )
    else:
        raw_pressure, source = _raw_pressure(reader, path, policy, manual_pressure)
        prestress = PrestressProvenance(
            tensor_kind=ElasticTensorKind.RAW_ENERGY_STRAIN,
            pressure_gpa=raw_pressure,
            pressure_source=source,
        )

    structure = reader.structure
    lattice = None
    symmetry = reader.symmetry
    metadata: dict[str, object] = {
        "backend": "crystal",
        "calculation": "elastic_constants",
        "prestress_applied_by_backend": reader.prestress_applied,
        "prestress_keyword": reader.prestress_keyword,
    }
    if structure is not None:
        metadata["parsed_structure_volume_angstrom3"] = structure.volume
        lattice_matches_volume = bool(
            np.isclose(
                structure.volume,
                reader.volume,
                rtol=2.0e-6,
                atol=1.0e-6,
            )
        )
        metadata["parsed_structure_matches_elastic_volume"] = lattice_matches_volume
        if lattice_matches_volume:
            lattice = structure.lattice
    if symmetry is not None:
        metadata["space_group_number"] = symmetry.space_group_number
        metadata["space_group_symbol"] = symmetry.international_symbol
    return ElasticState(
        volume=reader.volume,
        density=reader.density,
        stiffness=reader.stiffness,
        prestress=prestress,
        energy=reader.energy,
        energy_unit="hartree",
        lattice=lattice,
        source=path,
        metadata=metadata,
    )


def _raw_pressure(
    reader: CrystalElasticityReader,
    path: Path,
    policy: CrystalPressurePolicy,
    manual_pressure: float | None,
) -> tuple[float | None, PressureSource]:
    """Select pressure and provenance for a raw CRYSTAL tensor."""
    if policy is CrystalPressurePolicy.MANUAL:
        assert manual_pressure is not None
        return manual_pressure, PressureSource.MANUAL
    if policy is CrystalPressurePolicy.DEFERRED:
        return None, PressureSource.UNAVAILABLE
    if not np.isfinite(reader.stress_pressure):
        raise ValueError(
            f"CRYSTAL output {path} lacks pressure from the unstrained stress; "
            "supply pressure_policy='manual'"
        )
    return reader.stress_pressure, PressureSource.OUTPUT_STRESS


__all__ = ["CrystalPressurePolicy", "read_crystal_elastic_series"]
