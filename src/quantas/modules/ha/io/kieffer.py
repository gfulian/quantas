# -*- coding: utf-8 -*-

"""Enrich Quantas phonon YAML inputs with Kieffer acoustic cutoffs."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any, Literal

import numpy as np

from quantas.core.physics.elasticity import assign_hydrostatic_pressures
from quantas.core.physics.eos import (
    PressureEstimate,
    pressure_from_energy_eos,
    pressure_from_energy_polynomial,
)
from quantas.core.physics.kieffer import build_kieffer_volume_series
from quantas.interfaces.crystal import (
    CrystalPressurePolicy,
    correct_crystal_hydrostatic_elastic_series,
    read_crystal_elastic_series,
)
from quantas.io.kieffer import (
    kieffer_series_from_mapping,
    kieffer_series_to_mapping,
)
from quantas.io.phonons import PhononInputFileReader
from quantas.models.elastic_states import ElasticStateSeries, PressureSource
from quantas.models.kieffer import KiefferVolumeSeries
from quantas.models.kieffer_application import (
    validate_kieffer_phonon_applicability,
)
from quantas.models.phonons import PhononInputData
from quantas.models.volume_matching import VolumeMatch, match_sampled_volumes
from quantas.modules.ha.io.inpgen import format_quantas_yaml


KiefferInputWorkflow = Literal["ha", "qha"]
KiefferElasticInterface = Literal["crystal"]

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
    """Create a new HA/QHA YAML input containing Kieffer cutoff data.

    The source input is never overwritten implicitly.  Existing Kieffer data
    are rejected so replacement remains an explicit caller decision. For QHA,
    ``energy_eos`` and ``energy_polynomial`` derive hydrostatic pressures from
    the static energy-volume arrays in the phonon input, attach them to the raw
    elastic tensors, and then apply the Wallace correction exactly once.
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
        "prestress_correction": "crystal-erba-2014-hydrostatic",
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
    if interface != "crystal":
        raise ValueError(f"unsupported Kieffer elastic interface: {interface!r}")
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

    estimate, fit_provenance = _fit_phonon_input_pressure(
        phonon_input,
        pressure_source=pressure_source,
        eos=eos,
        polynomial_degree=polynomial_degree,
        maxfev=maxfev,
    )
    raw_series = read_crystal_elastic_series(
        elastic_outputs,
        pressure_policy=CrystalPressurePolicy.DEFERRED,
        apply_prestress_correction=False,
    )
    if phonon_input.volume is None:
        raise ValueError("QHA input does not contain sampled volumes")
    matches = match_sampled_volumes(raw_series.volumes, phonon_input.volume)
    pressures = np.asarray(
        [estimate.pressure[match.source_index] for match in matches],
        dtype=np.float64,
    )
    pressure_enum = (
        PressureSource.ENERGY_EOS
        if pressure_source == "energy_eos"
        else PressureSource.ENERGY_POLYNOMIAL
    )
    match_provenance = _volume_match_provenance(matches)
    fit_provenance["volume_matches"] = match_provenance
    assigned = assign_hydrostatic_pressures(
        raw_series,
        pressures,
        pressure_source=pressure_enum,
        assignment_method=pressure_source,
        metadata=fit_provenance,
    )
    corrected = correct_crystal_hydrostatic_elastic_series(
        assigned,
        correction_applied_by="quantas-kieffer-enrichment",
    )
    return corrected, fit_provenance


def _fit_phonon_input_pressure(
    phonon_input: PhononInputData,
    *,
    pressure_source: str,
    eos: str,
    polynomial_degree: int,
    maxfev: int | None,
) -> tuple[PressureEstimate, dict[str, Any]]:
    """Fit QHA static energy data and return pressures with provenance."""
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
    if pressure_source == "energy_eos":
        estimate = pressure_from_energy_eos(
            volume,
            energy,
            eos=eos,
            energy_unit=energy_unit,
            volume_unit=length_unit,
            pressure_unit="GPa",
            maxfev=maxfev,
        )
    else:
        estimate = pressure_from_energy_polynomial(
            volume,
            energy,
            degree=polynomial_degree,
            energy_unit=energy_unit,
            volume_unit=length_unit,
            pressure_unit="GPa",
        )
    if not estimate.success:
        detail = estimate.fit.message or "fit did not return finite pressures"
        raise ValueError(f"{pressure_source} pressure fit failed: {detail}")
    provenance = {
        "method": pressure_source,
        "relation": "P(V) = -dE/dV",
        "source_dataset": "phonon_input_static_energy",
        "energy_unit": energy_unit,
        "volume_unit": volume_unit,
        "volume_length_unit": length_unit,
        "pressure_unit": estimate.unit,
        "settings": dict(estimate.metadata),
        "evaluated_pressures_gpa": estimate.pressure.tolist(),
        "fit": estimate.fit.as_dict(),
        "warnings": list(estimate.warnings),
    }
    return estimate, provenance


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
    if normalized != "crystal":
        raise ValueError(f"unsupported Kieffer elastic interface: {value!r}")
    return "crystal"


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
    """Read embedded Kieffer cutoff data from a Quantas phonon YAML file."""
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
