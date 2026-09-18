# -*- coding: utf-8 -*-

"""Enrich Quantas phonon YAML inputs with Kieffer acoustic cutoffs."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any, Literal

import numpy as np

from quantas.core.physics.elasticity import (
    EnergyPressureResolution,
    resolve_energy_derived_pressures,
)
from quantas.core.physics.kieffer import build_kieffer_volume_series
from quantas.interfaces.crystal import (
    CrystalPressurePolicy,
    correct_crystal_hydrostatic_elastic_series,
    read_crystal_elastic_series,
)
from quantas.interfaces.vasp import (
    VASP_HYDROSTATIC_PRESTRESS_METHOD,
    assign_vasp_manual_pressures,
    convert_vasp_hydrostatic_elastic_series,
    read_vasp_elastic_series,
    resolve_vasp_energy_derived_pressures,
)
from quantas.io.kieffer import (
    kieffer_series_from_mapping,
    kieffer_series_to_mapping,
)
from quantas.io.phonons import PhononInputFileReader
from quantas.models.elastic_states import ElasticStateSeries
from quantas.models.kieffer import KiefferVolumeSeries
from quantas.models.kieffer_application import (
    validate_kieffer_phonon_applicability,
)
from quantas.models.phonons import PhononInputData
from quantas.models.volume_matching import VolumeMatch
from quantas.modules.ha.io.inpgen import format_quantas_yaml


KiefferInputWorkflow = Literal["ha", "qha"]
KiefferElasticInterface = Literal["crystal", "vasp"]

_ENERGY_PRESSURE_SOURCES = {"energy_eos", "energy_polynomial"}
_PRESSURE_SOURCES = {
    "auto",
    "output_stress",
    "manual",
    *_ENERGY_PRESSURE_SOURCES,
}


def add_kieffer_to_phonon_input(
    source: str | Path,
    destination: str | Path,
    elastic_outputs: Sequence[str | Path],
    *,
    workflow: KiefferInputWorkflow,
    interface: KiefferElasticInterface | str = "crystal",
    pressure_policy: CrystalPressurePolicy | str = CrystalPressurePolicy.AUTO,
    manual_pressures_gpa: Sequence[float] | None = None,
    eos: str = "BM3",
    polynomial_degree: int = 3,
    maxfev: int | None = None,
    mu_order: int = 12,
    phi_order: int = 24,
    refinement_factor: int = 2,
    batch_size: int = 512,
) -> Path:
    """Create a new HA/QHA YAML input containing Kieffer acoustic cutoff data.

    The source YAML is never overwritten. Elastic tensors are read through the
    selected external-code interface, converted to the hydrostatic finite-pressure
    form when required, and used to construct the three Kieffer branches. For QHA,
    ``energy_eos`` and ``energy_polynomial`` obtain ``P_static(V)`` from the static
    energy-volume series before the selected backend-specific finite-prestress
    conversion is applied exactly once. Energy-derived pressure is intentionally
    unavailable for single-volume HA.

    Parameters
    ----------
    source : str or Path
        Existing Quantas phonon YAML input.
    destination : str or Path
        New YAML path that will receive the Kieffer data.
    elastic_outputs : sequence of str or Path
        Elastic sources associated with the HA/QHA volume states. CRYSTAL uses
        output files; VASP accepts calculation directories, ``OUTCAR`` files, or
        ``vasprun.xml`` files with sibling ``OUTCAR`` files.
    workflow : {"ha", "qha"}
        Applicability contract to enforce.
    interface : str, optional
        Elastic interface. Supported values are ``"crystal"`` and ``"vasp"``.
    pressure_policy : CrystalPressurePolicy or str, optional
        Pressure source: ``auto``, ``output_stress``, ``manual``, ``energy_eos``, or
        ``energy_polynomial``.
    manual_pressures_gpa : sequence of float or None, optional
        Hydrostatic pressures in GPa for ``manual`` pressure assignment.
    eos : str, optional
        Energy EOS used when ``pressure_policy="energy_eos"``.
    polynomial_degree : int, optional
        Degree used when ``pressure_policy="energy_polynomial"``.
    maxfev : int or None, optional
        Maximum EnergyEOS optimizer evaluations when applicable.
    mu_order, phi_order : int, optional
        Polar and azimuthal spherical quadrature orders used for acoustic averaging.
    refinement_factor : int, optional
        Angular refinement factor used by the Kieffer velocity integration.
    batch_size : int, optional
        Number of directions processed per numerical batch.

    Returns
    -------
    Path
        Path to the newly written enriched YAML file.

    Raises
    ------
    ValueError
        If source/destination are invalid, Kieffer data already exist, the elastic
        interface or pressure source is unsupported, volumes cannot be matched, the
        phonon calculation violates the primitive Gamma-only contract, or an
        energy-derived pressure source is requested for HA.
    """
    source_path = Path(source)
    destination_path = Path(destination)
    if source_path.resolve() == destination_path.resolve():
        raise ValueError("Kieffer enrichment requires a new output path")
    reader = PhononInputFileReader(source_path)
    if not reader.completed:
        raise ValueError(reader.error or "Unable to read phonon YAML input")
    raw = reader.data
    if raw is None:
        raise ValueError("Phonon YAML input does not contain a mapping")
    if "kieffer" in raw:
        raise ValueError("Phonon YAML input already contains Kieffer data")

    phonon_input = reader.to_input(source=source_path)
    selected_interface = _normalize_interface(interface)
    selected_pressure = _normalize_pressure_source(pressure_policy)
    elastic_series, pressure_model = _prepare_elastic_series(
        phonon_input,
        elastic_outputs,
        workflow=workflow,
        interface=selected_interface,
        pressure_source=selected_pressure,
        manual_pressures_gpa=manual_pressures_gpa,
        eos=eos,
        polynomial_degree=polynomial_degree,
        maxfev=maxfev,
    )
    cutoffs = build_kieffer_volume_series(
        elastic_series,
        mu_order=mu_order,
        phi_order=phi_order,
        refinement_factor=refinement_factor,
        batch_size=batch_size,
    )
    validate_kieffer_phonon_applicability(
        phonon_input,
        cutoffs,
        workflow=workflow,
    )

    provenance: dict[str, Any] = {
        "phonon_input": str(source_path),
        "elastic_outputs": [str(state.source) for state in elastic_series.states],
        "elastic_interface": selected_interface,
        "pressure_source": selected_pressure,
        "prestress_correction": _prestress_correction_name(selected_interface),
        "reference_elastic_index": elastic_series.reference_index,
    }
    if pressure_model is not None:
        provenance["pressure_model"] = pressure_model
    raw["kieffer"] = kieffer_series_to_mapping(cutoffs, provenance=provenance)
    destination_path.write_text(format_quantas_yaml(raw), encoding="utf-8")
    return destination_path


def _prepare_elastic_series(
    phonon_input: PhononInputData,
    elastic_outputs: Sequence[str | Path],
    *,
    workflow: KiefferInputWorkflow,
    interface: KiefferElasticInterface,
    pressure_source: str,
    manual_pressures_gpa: Sequence[float] | None,
    eos: str,
    polynomial_degree: int,
    maxfev: int | None,
) -> tuple[ElasticStateSeries, dict[str, Any] | None]:
    """Read elastic outputs and resolve their hydrostatic pressure source."""
    if interface == "vasp":
        return _prepare_vasp_elastic_series(
            phonon_input,
            elastic_outputs,
            workflow=workflow,
            pressure_source=pressure_source,
            manual_pressures_gpa=manual_pressures_gpa,
            eos=eos,
            polynomial_degree=polynomial_degree,
            maxfev=maxfev,
        )

    if pressure_source not in _ENERGY_PRESSURE_SOURCES:
        series = read_crystal_elastic_series(
            elastic_outputs,
            pressure_policy=CrystalPressurePolicy(pressure_source),
            manual_pressures_gpa=manual_pressures_gpa,
            apply_prestress_correction=True,
        )
        return series, None
    if manual_pressures_gpa is not None:
        raise ValueError(
            "manual_pressures_gpa cannot be combined with an energy pressure source"
        )
    if workflow != "qha":
        raise ValueError(
            "energy-derived pressure requires a multi-volume QHA input; "
            "use output stress or manual pressure for HA"
        )
    if phonon_input.volume is None or phonon_input.energy is None:
        raise ValueError("energy-derived pressure requires QHA volume and energy data")

    volume = np.asarray(phonon_input.volume, dtype=np.float64)
    energy = np.asarray(phonon_input.energy, dtype=np.float64)
    if volume.size < 3:
        raise ValueError(
            "energy-derived pressure requires at least three volume-energy "
            "points; use --pressure-source manual"
        )
    energy_unit = str(phonon_input.units.get("energy", "Ha"))
    volume_unit = str(phonon_input.units.get("volume", "angstrom^3"))
    length_unit = str(phonon_input.units.get("length", "angstrom"))
    raw_series = read_crystal_elastic_series(
        elastic_outputs,
        pressure_policy=CrystalPressurePolicy.DEFERRED,
        apply_prestress_correction=False,
    )
    resolution = resolve_energy_derived_pressures(
        raw_series,
        volume,
        energy,
        pressure_source=pressure_source,
        source_dataset="phonon_input_static_energy",
        energy_unit=energy_unit,
        volume_length_unit=length_unit,
        volume_unit=volume_unit,
        eos=eos,
        polynomial_degree=polynomial_degree,
        maxfev=maxfev,
    )
    shared_provenance = resolution.provenance
    fit_provenance = {
        "method": shared_provenance["method"],
        "relation": shared_provenance["relation"],
        "source_dataset": shared_provenance["source_dataset"],
        "energy_unit": shared_provenance["energy_unit"],
        "volume_unit": volume_unit,
        "volume_length_unit": shared_provenance["volume_length_unit"],
        "pressure_unit": shared_provenance["pressure_unit"],
        "settings": shared_provenance["settings"],
        "evaluated_pressures_gpa": resolution.estimate.pressure.tolist(),
        "fit": shared_provenance["fit"],
        "warnings": shared_provenance["warnings"],
        "volume_matches": _volume_match_provenance(resolution.matches),
    }
    corrected = correct_crystal_hydrostatic_elastic_series(
        resolution.series,
        correction_applied_by="quantas-kieffer-enrichment",
    )
    return corrected, fit_provenance


def _prepare_vasp_elastic_series(
    phonon_input: PhononInputData,
    elastic_outputs: Sequence[str | Path],
    *,
    workflow: KiefferInputWorkflow,
    pressure_source: str,
    manual_pressures_gpa: Sequence[float] | None,
    eos: str,
    polynomial_degree: int,
    maxfev: int | None,
) -> tuple[ElasticStateSeries, dict[str, Any] | None]:
    """Prepare incremental VASP elastic states for Kieffer acoustics.

    VASP elastic moduli enter this adapter as raw stress--strain coefficients.
    Pressure selection is kept separate from tensor conversion: output stress,
    manual values, or static-energy-derived pressures are attached to the raw
    series first, then the VASP hydrostatic conversion is applied exactly once.

    Parameters
    ----------
    phonon_input : PhononInputData
        HA/QHA phonon input supplying static E(V) data when required.
    elastic_outputs : sequence of str or Path
        VASP run directories, ``OUTCAR`` files, or resolvable ``vasprun.xml``
        files.
    workflow : {"ha", "qha"}
        Thermodynamic workflow requesting the enrichment.
    pressure_source : str
        Normalized pressure source.
    manual_pressures_gpa : sequence of float or None
        Manual hydrostatic pressures in input-source order.
    eos : str
        Integrated energy EOS name used for ``energy_eos``.
    polynomial_degree : int
        Polynomial degree used for ``energy_polynomial``.
    maxfev : int or None
        Optional EOS fitter evaluation limit.

    Returns
    -------
    tuple
        Incremental elastic series and optional energy-pressure fit provenance.

    Raises
    ------
    ValueError
        If pressure selection is inconsistent with the VASP source data or the
        requested workflow.
    """
    raw_series = read_vasp_elastic_series(elastic_outputs)

    if pressure_source not in _ENERGY_PRESSURE_SOURCES:
        if pressure_source == "manual":
            if manual_pressures_gpa is None:
                raise ValueError(
                    "manual pressure source requires one pressure per VASP "
                    "elastic source"
                )
            selected = assign_vasp_manual_pressures(
                raw_series,
                manual_pressures_gpa,
                assignment_method="quantas-kieffer-manual-pressure",
            )
        else:
            if manual_pressures_gpa is not None:
                raise ValueError(
                    "manual_pressures_gpa requires pressure source 'manual'"
                )
            # VASP has no backend-corrected elastic tensor analogous to CRYSTAL
            # PRESSURE/PRESSEOS.  Both "auto" and "output_stress" therefore use
            # the unstrained VASP reference stress already attached by the reader.
            selected = raw_series

        corrected = convert_vasp_hydrostatic_elastic_series(
            selected,
            correction_applied_by="quantas-kieffer-enrichment",
        )
        return corrected, None

    if manual_pressures_gpa is not None:
        raise ValueError(
            "manual_pressures_gpa cannot be combined with an energy pressure source"
        )
    volume, energy, energy_unit, volume_unit, length_unit = _energy_pressure_inputs(
        phonon_input,
        workflow=workflow,
    )
    resolution = resolve_vasp_energy_derived_pressures(
        raw_series,
        volume,
        energy,
        pressure_source=pressure_source,
        source_dataset="phonon_input_static_energy",
        energy_unit=energy_unit,
        volume_length_unit=length_unit,
        volume_unit=volume_unit,
        eos=eos,
        polynomial_degree=polynomial_degree,
        maxfev=maxfev,
    )
    fit_provenance = _energy_pressure_provenance(
        resolution,
        volume_unit=volume_unit,
    )
    corrected = convert_vasp_hydrostatic_elastic_series(
        resolution.series,
        correction_applied_by="quantas-kieffer-enrichment",
    )
    return corrected, fit_provenance


def _energy_pressure_inputs(
    phonon_input: PhononInputData,
    *,
    workflow: KiefferInputWorkflow,
) -> tuple[np.ndarray, np.ndarray, str, str, str]:
    """Return validated static E(V) arrays and their declared units."""
    if workflow != "qha":
        raise ValueError(
            "energy-derived pressure requires a multi-volume QHA input; "
            "use output stress or manual pressure for HA"
        )
    if phonon_input.volume is None or phonon_input.energy is None:
        raise ValueError("energy-derived pressure requires QHA volume and energy data")
    volume = np.asarray(phonon_input.volume, dtype=np.float64)
    energy = np.asarray(phonon_input.energy, dtype=np.float64)
    if volume.size < 3:
        raise ValueError(
            "energy-derived pressure requires at least three volume-energy "
            "points; use --pressure-source manual"
        )
    energy_unit = str(phonon_input.units.get("energy", "Ha"))
    volume_unit = str(phonon_input.units.get("volume", "angstrom^3"))
    length_unit = str(phonon_input.units.get("length", "angstrom"))
    return volume, energy, energy_unit, volume_unit, length_unit


def _energy_pressure_provenance(
    resolution: EnergyPressureResolution,
    *,
    volume_unit: str,
) -> dict[str, Any]:
    """Return serialization-ready provenance for one E(V)-pressure resolution."""
    shared_provenance = resolution.provenance
    return {
        "method": shared_provenance["method"],
        "relation": shared_provenance["relation"],
        "source_dataset": shared_provenance["source_dataset"],
        "energy_unit": shared_provenance["energy_unit"],
        "volume_unit": volume_unit,
        "volume_length_unit": shared_provenance["volume_length_unit"],
        "pressure_unit": shared_provenance["pressure_unit"],
        "settings": shared_provenance["settings"],
        "evaluated_pressures_gpa": resolution.estimate.pressure.tolist(),
        "fit": shared_provenance["fit"],
        "warnings": shared_provenance["warnings"],
        "volume_matches": _volume_match_provenance(resolution.matches),
    }


def _prestress_correction_name(interface: KiefferElasticInterface) -> str:
    """Return the backend-specific correction identifier stored in YAML."""
    if interface == "vasp":
        return VASP_HYDROSTATIC_PRESTRESS_METHOD
    return "crystal-erba-2014-hydrostatic"


def _volume_match_provenance(matches: Sequence[VolumeMatch]) -> list[dict[str, Any]]:
    """Return serialization-ready elastic-to-phonon volume associations."""
    return [
        {
            "elastic_index": match.target_index,
            "phonon_index": match.source_index,
            "elastic_volume": match.target_volume,
            "phonon_volume": match.source_volume,
            "absolute_difference": match.absolute_difference,
            "relative_difference": match.relative_difference,
        }
        for match in matches
    ]


def _normalize_interface(
    value: KiefferElasticInterface | str,
) -> KiefferElasticInterface:
    """Return the supported elastic-output interface identifier."""
    normalized = str(value).strip().lower()
    if normalized not in {"crystal", "vasp"}:
        raise ValueError(f"unsupported Kieffer elastic interface: {value!r}")
    if normalized == "crystal":
        return "crystal"
    if normalized == "vasp":
        return "vasp"
    raise AssertionError("unreachable Kieffer interface normalization")


def _normalize_pressure_source(value: CrystalPressurePolicy | str) -> str:
    """Normalize CLI/API spelling of a Kieffer pressure source."""
    if isinstance(value, CrystalPressurePolicy):
        normalized = value.value
    else:
        normalized = str(value).strip().lower().replace("-", "_")
    if normalized not in _PRESSURE_SOURCES:
        choices = ", ".join(sorted(_PRESSURE_SOURCES))
        raise ValueError(
            f"unsupported pressure source {value!r}; choose from {choices}"
        )
    return normalized


def read_kieffer_from_phonon_input(source: str | Path) -> KiefferVolumeSeries:
    """Read embedded Kieffer cutoff data from a Quantas phonon YAML file.

    Parameters
    ----------
    source : str or Path
        Quantas HA/QHA phonon YAML file.

    Returns
    -------
    KiefferVolumeSeries
        Validated cutoff states, effective acoustic velocities, and provenance.

    Raises
    ------
    ValueError
        If the YAML cannot be read or does not contain a valid ``kieffer`` mapping.
    """
    reader = PhononInputFileReader(source)
    if not reader.completed:
        raise ValueError(reader.error or "Unable to read phonon YAML input")
    raw = reader.data or {}
    if "kieffer" not in raw:
        raise ValueError("Phonon YAML input does not contain Kieffer data")
    return kieffer_series_from_mapping(raw["kieffer"])


__all__ = [
    "KiefferElasticInterface",
    "KiefferInputWorkflow",
    "add_kieffer_to_phonon_input",
    "read_kieffer_from_phonon_input",
]
