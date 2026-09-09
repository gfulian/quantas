#!/usr/bin/env python3
"""Run the complete Kieffer QHA validation on the OHAp source dataset.

The utility deliberately keeps two comparisons separate:

* direct acoustic/optical composition at one sampled volume;
* the net effect of Kieffer on independently minimized QHA states.

It starts from the CRYSTAL output lists, exercises every supported pressure
route, refines the spherical quadrature, runs the four QHA workflow
combinations, checks HDF5 persistence, and records expected failure modes.
The raw source files are never modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable, Iterable

import h5py
import numpy as np
import scipy
import yaml

from quantas import __version__
from quantas.api import qha
from quantas.core.physics.thermodynamics import (
    entropy,
    isochoric_heat_capacity,
    kieffer_entropy,
    kieffer_isochoric_heat_capacity,
    kieffer_thermal_energy,
    kieffer_vibrational_free_energy,
    kieffer_zero_point_energy,
    thermal_energy,
    vibrational_free_energy,
    zero_point_energy,
)
from quantas.core.physics.units import convert_frequency
from quantas.interfaces.crystal import read_crystal_elastic_series
from quantas.io.kieffer import (
    kieffer_series_from_mapping,
    kieffer_series_to_mapping,
)
from quantas.models.kieffer import KiefferVolumeSeries
from quantas.modules.qha.kieffer import (
    matched_kieffer_arrays,
    validate_kieffer_qha_applicability,
)
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult


PRESSURE_METHODS = (
    "auto",
    "output_stress",
    "manual",
    "energy_eos",
    "energy_polynomial",
)
WORKFLOWS = (
    ("freq", "poly"),
    ("freq", "eos"),
    ("td", "poly"),
    ("td", "eos"),
)
PROBE_TEMPERATURES = np.asarray(
    [
        1.0,
        2.0,
        3.0,
        5.0,
        7.5,
        10.0,
        15.0,
        20.0,
        25.0,
        30.0,
        50.0,
        100.0,
        300.0,
        500.0,
        1000.0,
        1500.0,
        1.0e7,
    ],
    dtype=np.float64,
)
QHA_PROPERTIES: tuple[tuple[str, str], ...] = (
    ("equilibrium_volume", "angstrom^3"),
    ("zero_point_energy", "Ha cell^-1"),
    ("thermal_energy", "Ha cell^-1"),
    ("entropy", "Ha cell^-1 K^-1"),
    ("vibrational_free_energy", "Ha cell^-1"),
    ("isochoric_heat_capacity", "Ha cell^-1 K^-1"),
    ("isobaric_heat_capacity", "Ha cell^-1 K^-1"),
    ("heat_capacity_difference", "Ha cell^-1 K^-1"),
    ("isothermal_bulk_modulus", "GPa"),
    ("adiabatic_bulk_modulus", "GPa"),
    ("bulk_modulus_derivative", "1"),
    ("thermal_expansion", "K^-1"),
    ("gruneisen", "1"),
)


def parse_arguments() -> argparse.Namespace:
    """Return command-line settings for the OHAp validation."""
    parser = argparse.ArgumentParser(
        description=(
            "Validate Kieffer input generation and complete QHA execution "
            "from the OHAp CRYSTAL dataset."
        )
    )
    parser.add_argument(
        "dataset",
        type=Path,
        help="Directory containing qha-files.txt and elastic-files.txt",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("ohap_kieffer_validation"),
    )
    parser.add_argument("--source-archive", type=Path, default=None)
    parser.add_argument("--reference", type=int, default=4)
    parser.add_argument("--formula-units", type=int, default=2)
    parser.add_argument("--eos", default="BM3")
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument(
        "--temperature",
        nargs=3,
        type=float,
        metavar=("MIN", "MAX", "STEP"),
        default=(0.0, 1500.0, 10.0),
    )
    parser.add_argument(
        "--pressure",
        nargs=3,
        type=float,
        metavar=("MIN", "MAX", "STEP"),
        default=(0.0, 0.0, 1.0),
    )
    parser.add_argument("--mu-order", type=int, default=12)
    parser.add_argument("--phi-order", type=int, default=24)
    parser.add_argument("--refinement-factor", type=int, default=2)
    parser.add_argument(
        "--pressure-agreement-limit",
        type=float,
        default=1.0,
        help="Maximum accepted fitted-vs-output pressure difference in GPa",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace files generated previously in the selected output directory",
    )
    return parser.parse_args()


def _listed_paths(filename: Path) -> tuple[Path, ...]:
    """Resolve non-comment entries relative to a portable list file."""
    paths: list[Path] = []
    for line in filename.read_text(encoding="utf-8").splitlines():
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        path = Path(value)
        if not path.is_absolute():
            path = filename.parent / path
        if not path.is_file():
            raise FileNotFoundError(f"listed source does not exist: {path}")
        paths.append(path.resolve())
    if not paths:
        raise ValueError(f"source list is empty: {filename}")
    return tuple(paths)


def _sha256(path: Path) -> str:
    """Return the SHA-256 digest of one file without loading it all at once."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_manifest(paths: Iterable[Path], root: Path) -> dict[str, Any]:
    """Build file-level and aggregate provenance for the QM sources."""
    records: list[dict[str, Any]] = []
    aggregate = hashlib.sha256()
    for path in sorted(paths, key=lambda item: item.as_posix()):
        relative = path.relative_to(root).as_posix()
        checksum = _sha256(path)
        records.append(
            {"path": relative, "size_bytes": path.stat().st_size, "sha256": checksum}
        )
        aggregate.update(relative.encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(bytes.fromhex(checksum))
    return {"aggregate_sha256": aggregate.hexdigest(), "files": records}


def _prepare_output_directory(path: Path, force: bool) -> Path:
    """Create an output directory and protect previous validation results."""
    path = path.resolve()
    path.mkdir(parents=True, exist_ok=True)
    generated = tuple(path.glob("ohap-*"))
    if generated and not force:
        raise FileExistsError(
            f"{path} already contains OHAp validation files; use --force"
        )
    return path


def _input_order_pressures(elastic_paths: tuple[Path, ...]) -> tuple[float, ...]:
    """Read output-stress pressures and restore the list-file order."""
    series = read_crystal_elastic_series(
        elastic_paths,
        pressure_policy="output_stress",
        apply_prestress_correction=False,
    )
    by_source = {
        Path(state.source).resolve(): float(state.prestress.pressure_gpa)
        for state in series.states
    }
    return tuple(by_source[path.resolve()] for path in elastic_paths)


def _create_inputs(
    dataset: Path,
    output: Path,
    arguments: argparse.Namespace,
    elastic_paths: tuple[Path, ...],
) -> tuple[Path, dict[str, Path], tuple[float, ...]]:
    """Generate the phonon input and all pressure/quadrature variants."""
    qha_list = dataset / "qha-files.txt"
    base = output / "ohap-qha.yaml"
    qha.create_input(
        qha_list,
        base,
        interface="crystal",
        is_list=True,
        reference=arguments.reference,
        jobname="OHAp Kieffer end-to-end validation",
        formula_units=arguments.formula_units,
    )
    manual_pressures = _input_order_pressures(elastic_paths)
    generated: dict[str, Path] = {}
    for pressure_method in PRESSURE_METHODS:
        destination = output / f"ohap-{pressure_method.replace('_', '-')}.yaml"
        kwargs: dict[str, Any] = {}
        if pressure_method == "manual":
            kwargs["manual_pressures_gpa"] = manual_pressures
        qha.add_kieffer_input(
            base,
            destination,
            elastic_paths,
            interface="crystal",
            pressure_policy=pressure_method,
            eos=arguments.eos,
            polynomial_degree=arguments.degree,
            mu_order=arguments.mu_order,
            phi_order=arguments.phi_order,
            refinement_factor=arguments.refinement_factor,
            **kwargs,
        )
        generated[pressure_method] = destination

    reversed_path = output / "ohap-output-stress-reversed.yaml"
    qha.add_kieffer_input(
        base,
        reversed_path,
        tuple(reversed(elastic_paths)),
        pressure_policy="output_stress",
        mu_order=arguments.mu_order,
        phi_order=arguments.phi_order,
        refinement_factor=arguments.refinement_factor,
    )
    generated["output_stress_reversed"] = reversed_path

    for multiplier, label in ((2, "dense"), (4, "fine")):
        destination = output / f"ohap-energy-eos-{label}.yaml"
        qha.add_kieffer_input(
            base,
            destination,
            elastic_paths,
            pressure_policy="energy_eos",
            eos=arguments.eos,
            mu_order=arguments.mu_order * multiplier,
            phi_order=arguments.phi_order * multiplier,
            refinement_factor=arguments.refinement_factor,
        )
        generated[f"energy_eos_{label}"] = destination
    return base, generated, manual_pressures


def _state_pressures(series: KiefferVolumeSeries) -> np.ndarray:
    """Return the pressure recorded on every direct cutoff state."""
    return np.asarray(
        [float(state.metadata["pressure_gpa"]) for state in series.states],
        dtype=np.float64,
    )


def _pressure_analysis(
    inputs: dict[str, Path],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Compare pressure routes, cutoff changes, and quadrature refinement."""
    series = {name: qha.read_kieffer_input(path) for name, path in inputs.items()}
    reference = series["output_stress"]
    rows: list[dict[str, Any]] = []
    pressure_arrays = {
        name: _state_pressures(series[name]) for name in PRESSURE_METHODS
    }
    for index, volume in enumerate(reference.volumes):
        row: dict[str, Any] = {"volume_angstrom3": float(volume)}
        for name in PRESSURE_METHODS:
            row[f"pressure_{name}_gpa"] = float(pressure_arrays[name][index])
        row["energy_eos_minus_output_gpa"] = float(
            pressure_arrays["energy_eos"][index]
            - pressure_arrays["output_stress"][index]
        )
        row["energy_polynomial_minus_output_gpa"] = float(
            pressure_arrays["energy_polynomial"][index]
            - pressure_arrays["output_stress"][index]
        )
        rows.append(row)

    exact_names = ("auto", "manual", "output_stress_reversed")
    exact_equivalence = {
        name: bool(
            np.array_equal(series[name].frequencies_hz, reference.frequencies_hz)
            and np.array_equal(
                series[name].effective_velocities_km_s,
                reference.effective_velocities_km_s,
            )
        )
        for name in exact_names
    }
    fitted_sensitivity: dict[str, Any] = {}
    for name in ("energy_eos", "energy_polynomial"):
        relative = (
            series[name].frequencies_hz - reference.frequencies_hz
        ) / reference.frequencies_hz
        fitted_sensitivity[name] = {
            "maximum_pressure_difference_gpa": float(
                np.max(np.abs(pressure_arrays[name] - pressure_arrays["output_stress"]))
            ),
            "maximum_cutoff_relative_difference": float(np.max(np.abs(relative))),
        }

    default = series["energy_eos"]
    dense = series["energy_eos_dense"]
    fine = series["energy_eos_fine"]
    quadrature = {
        "default_to_dense_maximum_relative_change": float(
            np.max(
                np.abs(
                    (dense.frequencies_hz - default.frequencies_hz)
                    / dense.frequencies_hz
                )
            )
        ),
        "dense_to_fine_maximum_relative_change": float(
            np.max(
                np.abs(
                    (fine.frequencies_hz - dense.frequencies_hz) / fine.frequencies_hz
                )
            )
        ),
        "fine_reported_maximum_relative_error": float(
            max(
                max(state.metadata["quadrature"]["relative_errors"])
                for state in fine.states
            )
        ),
    }
    return rows, {
        "exact_equivalence": exact_equivalence,
        "fitted_sensitivity": fitted_sensitivity,
        "quadrature": quadrature,
    }


def _qha_options(
    arguments: argparse.Namespace,
    scheme: str,
    minimization: str,
) -> QHAOptions:
    """Build one controlled QHA option set for paired calculations."""
    tmin, tmax, tstep = arguments.temperature
    pmin, pmax, pstep = arguments.pressure
    return qha.Options(
        temperature_min=tmin,
        temperature_max=tmax,
        temperature_step=tstep,
        pressure_min=pmin,
        pressure_max=pmax,
        pressure_step=pstep,
        scheme=scheme,
        minimization=minimization,
        eos=arguments.eos,
        energy_degree=arguments.degree,
        free_energy_degree=arguments.degree,
        frequency_degree=arguments.degree,
        thermal_expansion_method="mixed_derivative",
        calculate_mode_gruneisen=False,
        estimate_uncertainties=False,
        fit_failure_policy="raise",
        extrapolation_policy="fail",
    )


def _run_workflows(
    input_data: QHAInput,
    cutoffs: KiefferVolumeSeries,
    arguments: argparse.Namespace,
) -> tuple[
    dict[str, tuple[QHAResult, QHAResult]],
    dict[str, Any],
    dict[str, tuple[Any, Any]],
]:
    """Execute normal and Kieffer QHA for every supported workflow."""
    results: dict[str, tuple[QHAResult, QHAResult]] = {}
    summaries: dict[str, Any] = {}
    envelopes: dict[str, tuple[Any, Any]] = {}
    for scheme, minimization in WORKFLOWS:
        name = f"{scheme}_{minimization}"
        options = _qha_options(arguments, scheme, minimization)
        normal_envelope = qha.run(input_data, options=options)
        kieffer_envelope = qha.run(
            input_data,
            options=options,
            kieffer_cutoffs=cutoffs,
        )
        normal = qha.get_result(normal_envelope)
        enriched = qha.get_result(kieffer_envelope)
        results[name] = (normal, enriched)
        envelopes[name] = (normal_envelope, kieffer_envelope)
        summaries[name] = {
            "phonons_only": qha.validate_result(normal, input_data).as_dict(),
            "with_kieffer": qha.validate_result(enriched, input_data).as_dict(),
            "differences": [
                {
                    "property": item.property_name,
                    "maximum_absolute": item.maximum_absolute,
                    "maximum_relative": item.maximum_relative,
                    "root_mean_square": item.root_mean_square,
                    "maximum_absolute_temperature": item.maximum_absolute_temperature,
                    "maximum_absolute_pressure": item.maximum_absolute_pressure,
                }
                for item in qha.compare_results(normal, enriched)
            ],
        }
    return results, summaries, envelopes


def _direct_composition(
    input_data: QHAInput,
    cutoffs: KiefferVolumeSeries,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Evaluate exact acoustic and optical shares at the static minimum volume."""
    if input_data.frequencies is None or input_data.energy is None:
        raise ValueError("OHAp validation requires energies and frequencies")
    reference_index = int(np.argmin(np.asarray(input_data.energy, dtype=np.float64)))
    matches = validate_kieffer_qha_applicability(input_data, cutoffs)
    matched_cutoffs, _ = matched_kieffer_arrays(cutoffs, matches)
    frequencies = np.asarray(input_data.frequencies, dtype=np.float64)[
        :, :, reference_index : reference_index + 1
    ]
    frequency_hz = np.asarray(
        convert_frequency(
            frequencies,
            str(input_data.units.get("frequency", "cm^-1")),
            "Hz",
        ),
        dtype=np.float64,
    )
    cutoff_hz = matched_cutoffs[:, reference_index : reference_index + 1]
    weights = input_data.normalized_weights()
    temperature = PROBE_TEMPERATURES

    optical_zpe = np.broadcast_to(
        zero_point_energy(np.zeros(1), frequency_hz, weights),
        (temperature.size, 1),
    )
    acoustic_zpe = np.broadcast_to(
        kieffer_zero_point_energy(np.zeros(1), cutoff_hz),
        (temperature.size, 1),
    )
    values = {
        "zero_point_energy": (optical_zpe, acoustic_zpe, "kJ mol^-1"),
        "thermal_energy": (
            thermal_energy(temperature, frequency_hz, weights),
            kieffer_thermal_energy(temperature, cutoff_hz),
            "kJ mol^-1",
        ),
        "entropy": (
            entropy(temperature, frequency_hz, weights),
            kieffer_entropy(temperature, cutoff_hz),
            "J mol^-1 K^-1",
        ),
        "isochoric_heat_capacity": (
            isochoric_heat_capacity(temperature, frequency_hz, weights),
            kieffer_isochoric_heat_capacity(temperature, cutoff_hz),
            "J mol^-1 K^-1",
        ),
        "vibrational_free_energy": (
            vibrational_free_energy(temperature, frequency_hz, weights),
            kieffer_vibrational_free_energy(temperature, cutoff_hz),
            "kJ mol^-1",
        ),
    }
    rows: list[dict[str, Any]] = []
    for index, temperature_value in enumerate(temperature):
        for name, (optical, acoustic, unit) in values.items():
            optical_value = float(optical[index, 0])
            acoustic_value = float(acoustic[index, 0])
            total = optical_value + acoustic_value
            positive_fraction = name in {
                "zero_point_energy",
                "thermal_energy",
                "entropy",
                "isochoric_heat_capacity",
            }
            rows.append(
                {
                    "temperature_k": float(temperature_value),
                    "property": name,
                    "unit": unit,
                    "optical": optical_value,
                    "acoustic": acoustic_value,
                    "total": total,
                    "acoustic_percent": (
                        100.0 * acoustic_value / total
                        if positive_fraction and total > 0.0
                        else None
                    ),
                    "optical_percent": (
                        100.0 * optical_value / total
                        if positive_fraction and total > 0.0
                        else None
                    ),
                }
            )

    positive_modes = int(np.sum(frequencies > 0.0))
    translational_modes = int(np.sum(frequencies <= 0.0))
    return rows, {
        "reference_index": reference_index,
        "reference_volume_angstrom3": float(input_data.volume[reference_index]),
        "positive_gamma_modes": positive_modes,
        "nonpositive_gamma_modes": translational_modes,
        "classical_acoustic_cv_percent": 100.0 * 3.0 / (positive_modes + 3.0),
    }


def _qha_effect_rows(
    results: dict[str, tuple[QHAResult, QHAResult]],
) -> list[dict[str, Any]]:
    """Return long-form paired differences on common pressure-temperature grids."""
    rows: list[dict[str, Any]] = []
    for workflow, (normal, enriched) in results.items():
        temperature = np.asarray(normal.temperature, dtype=np.float64)
        pressure = np.asarray(normal.pressure, dtype=np.float64)
        for attribute, unit in QHA_PROPERTIES:
            normal_value = getattr(normal, attribute)
            enriched_value = getattr(enriched, attribute)
            if normal_value is None or enriched_value is None:
                continue
            first = np.asarray(normal_value, dtype=np.float64)
            second = np.asarray(enriched_value, dtype=np.float64)
            for it, temperature_value in enumerate(temperature):
                for ip, pressure_value in enumerate(pressure):
                    baseline = float(first[it, ip])
                    corrected = float(second[it, ip])
                    delta = corrected - baseline
                    rows.append(
                        {
                            "workflow": workflow,
                            "temperature_k": float(temperature_value),
                            "pressure_gpa": float(pressure_value),
                            "property": attribute,
                            "unit": unit,
                            "phonons_only": baseline,
                            "with_kieffer": corrected,
                            "delta": delta,
                            "delta_over_phonons_percent": (
                                100.0 * delta / baseline if baseline != 0.0 else None
                            ),
                            "delta_over_enriched_percent": (
                                100.0 * delta / corrected if corrected != 0.0 else None
                            ),
                        }
                    )
    return rows


def _expect_error(
    label: str,
    operation: Callable[[], Any],
    expected_text: str,
) -> dict[str, Any]:
    """Run one destructive-input test and retain its expected diagnostic."""
    try:
        operation()
    except Exception as exc:  # the exception family is part of the record
        message = str(exc)
        return {
            "name": label,
            "passed": expected_text.lower() in message.lower(),
            "exception": type(exc).__name__,
            "message": message,
            "expected_text": expected_text,
        }
    return {
        "name": label,
        "passed": False,
        "exception": None,
        "message": "operation unexpectedly succeeded",
        "expected_text": expected_text,
    }


def _adversarial_checks(
    base: Path,
    enriched: Path,
    input_data: QHAInput,
    cutoffs: KiefferVolumeSeries,
    elastic_paths: tuple[Path, ...],
    output: Path,
    arguments: argparse.Namespace,
) -> list[dict[str, Any]]:
    """Exercise scientific and I/O failure boundaries without editing sources."""
    checks: list[dict[str, Any]] = []
    checks.append(
        _expect_error(
            "base input has no embedded Kieffer block",
            lambda: qha.read_kieffer_input(base),
            "does not contain Kieffer data",
        )
    )

    non_gamma = deepcopy(input_data)
    non_gamma.qcoords = np.asarray([[0.1, 0.0, 0.0]], dtype=np.float64)
    checks.append(
        _expect_error(
            "non-Gamma q-point",
            lambda: validate_kieffer_qha_applicability(non_gamma, cutoffs),
            "Gamma q-point",
        )
    )
    supercell = deepcopy(input_data)
    supercell.supercell = np.diag([2.0, 1.0, 1.0])
    checks.append(
        _expect_error(
            "phonon supercell",
            lambda: validate_kieffer_qha_applicability(supercell, cutoffs),
            "identity supercell",
        )
    )
    incomplete = KiefferVolumeSeries(states=cutoffs.states[:-1])
    checks.append(
        _expect_error(
            "missing cutoff volume",
            lambda: validate_kieffer_qha_applicability(input_data, incomplete),
            "one cutoff state per sampled volume",
        )
    )
    shifted = deepcopy(cutoffs)
    shifted.states[0].volume += 0.1
    checks.append(
        _expect_error(
            "mismatched cutoff volume",
            lambda: validate_kieffer_qha_applicability(input_data, shifted),
            "no source match",
        )
    )

    invalid_mapping = kieffer_series_to_mapping(cutoffs)
    invalid_mapping["states"][0]["cutoff_frequency"][0] = -1.0
    checks.append(
        _expect_error(
            "negative cutoff frequency",
            lambda: kieffer_series_from_mapping(invalid_mapping),
            "three finite positive values",
        )
    )

    modal_options = replace(
        _qha_options(arguments, "freq", "poly"),
        calculate_mode_gruneisen=True,
    )
    checks.append(
        _expect_error(
            "incomplete mode-Gruneisen coupling",
            lambda: qha.run(
                input_data,
                options=modal_options,
                kieffer_cutoffs=cutoffs,
            ),
            "does not yet support mode-Gruneisen",
        )
    )
    checks.append(
        _expect_error(
            "duplicate elastic source",
            lambda: read_crystal_elastic_series(
                (elastic_paths[0], elastic_paths[0]),
                pressure_policy="output_stress",
            ),
            "paths must be unique",
        )
    )
    checks.append(
        _expect_error(
            "in-place enrichment",
            lambda: qha.add_kieffer_input(base, base, elastic_paths),
            "new output path",
        )
    )
    checks.append(
        _expect_error(
            "silent replacement of existing block",
            lambda: qha.add_kieffer_input(
                enriched,
                output / "ohap-invalid-replacement.yaml",
                elastic_paths,
            ),
            "already contains Kieffer data",
        )
    )
    checks.append(
        _expect_error(
            "incomplete elastic volume series",
            lambda: qha.add_kieffer_input(
                base,
                output / "ohap-invalid-incomplete.yaml",
                elastic_paths[:-1],
                pressure_policy="output_stress",
                mu_order=2,
                phi_order=4,
                refinement_factor=2,
            ),
            "one cutoff state per sampled volume",
        )
    )
    return checks


def _round_trip_check(
    output: Path,
    envelopes: tuple[Any, Any],
    results: tuple[QHAResult, QHAResult],
) -> dict[str, Any]:
    """Persist and reload the primary pair, checking every public property."""
    files = (
        output / "ohap-phonons-only.hdf5",
        output / "ohap-with-kieffer.hdf5",
    )
    property_checks: dict[str, bool] = {}
    contribution_checks: dict[str, bool] = {}
    for label, envelope, original, filename in zip(
        ("phonons_only", "with_kieffer"),
        envelopes,
        results,
        files,
    ):
        qha.write_result(envelope, filename)
        restored = qha.get_result(qha.read_result(filename))
        original_properties = original.as_property_dict()
        restored_properties = restored.as_property_dict()
        property_equal = all(
            np.array_equal(
                np.asarray(original_value),
                np.asarray(restored_properties[property_name]),
                equal_nan=True,
            )
            for property_name, original_value in original_properties.items()
        )
        first = original.kieffer_sampled_contribution
        second = restored.kieffer_sampled_contribution
        if first is None or second is None:
            contribution_equal = first is None and second is None
        else:
            contribution_equal = all(
                np.array_equal(
                    np.asarray(getattr(first, attribute)),
                    np.asarray(getattr(second, attribute)),
                    equal_nan=True,
                )
                for attribute in (
                    "cutoff_frequencies_hz",
                    "effective_velocities_km_s",
                    "zero_point_energy",
                    "thermal_energy",
                    "entropy",
                    "vibrational_free_energy",
                    "isochoric_heat_capacity",
                )
            ) and json.dumps(
                first.metadata,
                sort_keys=True,
                default=_json_default,
            ) == json.dumps(
                second.metadata,
                sort_keys=True,
                default=_json_default,
            )
        property_checks[label] = bool(property_equal)
        contribution_checks[label] = bool(contribution_equal)
    return {
        "files": [path.name for path in files],
        "exact_public_property_round_trip": property_checks,
        "exact_kieffer_contribution_round_trip": contribution_checks,
    }


def _fraction(
    rows: list[dict[str, Any]],
    property_name: str,
    temperature: float,
) -> float:
    """Select one direct acoustic percentage from the long table."""
    for row in rows:
        if row["property"] == property_name and row["temperature_k"] == temperature:
            return float(row["acoustic_percent"])
    raise KeyError((property_name, temperature))


def _acceptance_checks(
    pressure: dict[str, Any],
    workflow: dict[str, Any],
    composition_rows: list[dict[str, Any]],
    composition: dict[str, Any],
    adversarial: list[dict[str, Any]],
    round_trip: dict[str, Any],
    pressure_limit: float,
) -> dict[str, bool]:
    """Evaluate declared numerical and scientific acceptance criteria."""
    classical = float(composition["classical_acoustic_cv_percent"])
    cv_low = _fraction(composition_rows, "isochoric_heat_capacity", 1.0)
    cv_room = _fraction(composition_rows, "isochoric_heat_capacity", 300.0)
    cv_high = _fraction(composition_rows, "isochoric_heat_capacity", 1500.0)
    cv_classical_probe = _fraction(
        composition_rows,
        "isochoric_heat_capacity",
        1.0e7,
    )
    all_workflows = [
        summary[variant]
        for summary in workflow.values()
        for variant in ("phonons_only", "with_kieffer")
    ]
    sensitivity = pressure["fitted_sensitivity"]
    return {
        "auto_output_manual_and_reordered_are_exact": all(
            pressure["exact_equivalence"].values()
        ),
        "fitted_pressures_agree_with_output_within_limit": all(
            sensitivity[name]["maximum_pressure_difference_gpa"] <= pressure_limit
            for name in ("energy_eos", "energy_polynomial")
        ),
        "fine_quadrature_refinement_below_2e-4": pressure["quadrature"][
            "fine_reported_maximum_relative_error"
        ]
        < 2.0e-4,
        "all_qha_workflows_complete_and_finite": all(
            summary["completed"]
            and summary["finite_properties"]
            and summary["valid_points"] == summary["total_points"]
            for summary in all_workflows
        ),
        "all_equilibrium_volumes_inside_sampled_range": all(
            summary["volumes_below_sampled_range"] == 0
            and summary["volumes_above_sampled_range"] == 0
            for summary in all_workflows
        ),
        "low_to_high_temperature_acoustic_cv_share_decreases": (
            cv_low > cv_room > cv_high > classical
        ),
        "high_temperature_cv_reaches_branch_count_limit": abs(
            cv_classical_probe - classical
        )
        < 1.0e-6,
        "all_expected_failures_are_rejected": all(
            check["passed"] for check in adversarial
        ),
        "primary_hdf5_round_trip_is_exact": all(
            round_trip["exact_public_property_round_trip"].values()
        )
        and all(round_trip["exact_kieffer_contribution_round_trip"].values()),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write deterministic long-form CSV data."""
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _git_revision() -> str | None:
    """Return the current Git revision when the tool runs in a work tree."""
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def _selected_rows(
    rows: list[dict[str, Any]],
    *,
    workflow: str,
    property_name: str,
    pressure: float,
    temperatures: tuple[float, ...],
) -> list[dict[str, Any]]:
    """Select exact report rows without silently snapping grid values."""
    return [
        row
        for row in rows
        if row.get("workflow") == workflow
        and row["property"] == property_name
        and row.get("pressure_gpa") == pressure
        and row["temperature_k"] in temperatures
    ]


def _render_markdown(data: dict[str, Any]) -> str:
    """Render the compact human-readable OHAp validation report."""
    lines = [
        "# OHAp Kieffer end-to-end validation",
        "",
        f"- Quantas: `{data['software']['quantas_version']}`",
        f"- Git revision: `{data['software']['git_revision']}`",
        f"- QM-source aggregate SHA-256: `{data['dataset']['sources']['aggregate_sha256']}`",
        f"- Source archive SHA-256: `{data['dataset'].get('source_archive_sha256')}`",
        f"- Volumes: {data['input']['nvolumes']}",
        f"- Atoms / Gamma modes: {data['input']['natoms']} / {data['input']['nmodes']}",
        "",
        "## Pressure and quadrature checks",
        "",
        "| Check | Result |",
        "|---|---:|",
    ]
    for name, value in data["pressure_analysis"]["exact_equivalence"].items():
        lines.append(f"| `{name}` equals explicit output stress | {value} |")
    for name, values in data["pressure_analysis"]["fitted_sensitivity"].items():
        lines.append(
            f"| `{name}` max |P-P_output| (GPa) | "
            f"{values['maximum_pressure_difference_gpa']:.6f} |"
        )
        lines.append(
            f"| `{name}` max cutoff change (%) | "
            f"{100.0 * values['maximum_cutoff_relative_difference']:.6f} |"
        )
    for name, value in data["pressure_analysis"]["quadrature"].items():
        lines.append(f"| {name} | {value:.8e} |")

    lines.extend(
        [
            "",
            "## Direct acoustic share at the static minimum volume",
            "",
            "Percentages use acoustic / (acoustic + optical) at one common volume.",
            "",
            "| T (K) | Uth acoustic (%) | S acoustic (%) | Cv acoustic (%) |",
            "|---:|---:|---:|---:|",
        ]
    )
    composition_rows = data["direct_composition"]["rows"]
    for temperature in (
        1.0,
        5.0,
        10.0,
        20.0,
        50.0,
        100.0,
        300.0,
        1000.0,
        1500.0,
        1.0e7,
    ):
        values = {
            row["property"]: row["acoustic_percent"]
            for row in composition_rows
            if row["temperature_k"] == temperature
        }
        lines.append(
            f"| {temperature:g} | {values['thermal_energy']:.6f} | "
            f"{values['entropy']:.6f} | "
            f"{values['isochoric_heat_capacity']:.6f} |"
        )

    lines.extend(
        [
            "",
            "## Net frequency-polynomial QHA effect at zero pressure",
            "",
            "These are differences between independent equilibrium calculations, not pure branch fractions.",
            "",
            "| T (K) | Delta V (A^3) | Delta V (%) | Delta Cv (%) | Delta S (%) | Delta KT (%) |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )
    effect_rows = data["qha_effect_rows"]
    selected_temperatures = (0.0, 10.0, 50.0, 100.0, 300.0, 1000.0, 1500.0)
    by_property = {
        name: {
            row["temperature_k"]: row
            for row in _selected_rows(
                effect_rows,
                workflow="freq_poly",
                property_name=name,
                pressure=0.0,
                temperatures=selected_temperatures,
            )
        }
        for name in (
            "equilibrium_volume",
            "isochoric_heat_capacity",
            "entropy",
            "isothermal_bulk_modulus",
        )
    }
    for temperature in selected_temperatures:
        volume = by_property["equilibrium_volume"].get(temperature)
        if volume is None:
            continue
        cv = by_property["isochoric_heat_capacity"][temperature]
        entropy_row = by_property["entropy"][temperature]
        bulk = by_property["isothermal_bulk_modulus"][temperature]
        lines.append(
            f"| {temperature:g} | {volume['delta']:.8f} | "
            f"{_optional_number(volume['delta_over_phonons_percent'])} | "
            f"{_optional_number(cv['delta_over_phonons_percent'])} | "
            f"{_optional_number(entropy_row['delta_over_phonons_percent'])} | "
            f"{_optional_number(bulk['delta_over_phonons_percent'])} |"
        )

    lines.extend(
        [
            "",
            "## Expected failures",
            "",
            "| Case | Rejected | Diagnostic |",
            "|---|---:|---|",
        ]
    )
    for check in data["adversarial_checks"]:
        message = check["message"].replace("|", "\\|")
        lines.append(f"| {check['name']} | {check['passed']} | `{message}` |")

    lines.extend(
        [
            "",
            "## Acceptance",
            "",
            "| Criterion | Passed |",
            "|---|---:|",
        ]
    )
    for name, passed in data["acceptance"].items():
        lines.append(f"| {name} | {passed} |")
    lines.append("")
    return "\n".join(lines)


def _optional_number(value: float | None) -> str:
    """Format one optional percentage for Markdown."""
    return "n/a" if value is None else f"{value:.6f}"


def _json_default(value: Any) -> Any:
    """Convert NumPy scalars and arrays for deterministic JSON output."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"cannot serialize {type(value).__name__}")


def main() -> None:
    """Run the validation and write machine- and human-readable evidence."""
    arguments = parse_arguments()
    dataset = arguments.dataset.resolve()
    if not dataset.is_dir():
        raise NotADirectoryError(dataset)
    output = _prepare_output_directory(arguments.output_dir, arguments.force)
    qha_list = dataset / "qha-files.txt"
    elastic_list = dataset / "elastic-files.txt"
    qha_paths = _listed_paths(qha_list)
    elastic_paths = _listed_paths(elastic_list)
    all_sources = (*qha_paths, *elastic_paths, qha_list, elastic_list)

    base, inputs, manual_pressures = _create_inputs(
        dataset,
        output,
        arguments,
        elastic_paths,
    )
    primary = inputs["energy_eos_fine"]
    input_data = qha.read_input(primary)
    cutoffs = qha.read_kieffer_input(primary)
    pressure_rows, pressure_analysis = _pressure_analysis(inputs)
    results, workflow_summaries, envelopes = _run_workflows(
        input_data,
        cutoffs,
        arguments,
    )
    composition_rows, composition_summary = _direct_composition(
        input_data,
        cutoffs,
    )
    effect_rows = _qha_effect_rows(results)
    adversarial = _adversarial_checks(
        base,
        primary,
        input_data,
        cutoffs,
        elastic_paths,
        output,
        arguments,
    )
    round_trip = _round_trip_check(
        output,
        envelopes["freq_poly"],
        results["freq_poly"],
    )
    acceptance = _acceptance_checks(
        pressure_analysis,
        workflow_summaries,
        composition_rows,
        composition_summary,
        adversarial,
        round_trip,
        arguments.pressure_agreement_limit,
    )

    source_archive_sha256 = None
    if arguments.source_archive is not None:
        source_archive_sha256 = _sha256(arguments.source_archive.resolve())
    raw_mapping = yaml.safe_load(primary.read_text(encoding="utf-8"))
    data = {
        "title": "OHAp Kieffer end-to-end validation",
        "software": {
            "quantas_version": __version__,
            "git_revision": _git_revision(),
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "h5py_version": h5py.__version__,
        },
        "dataset": {
            "root": str(dataset),
            "source_archive_sha256": source_archive_sha256,
            "sources": _source_manifest(all_sources, dataset),
        },
        "input": {
            "nvolumes": input_data.nvol,
            "natoms": input_data.natoms,
            "nmodes": int(np.asarray(input_data.frequencies).shape[1]),
            "formula_units": input_data.formula_units,
            "mode_continuity": input_data.mode_continuity,
            "volume_range_angstrom3": [
                float(np.min(input_data.volume)),
                float(np.max(input_data.volume)),
            ],
            "manual_pressures_input_order_gpa": manual_pressures,
            "primary_pressure_provenance": raw_mapping["kieffer"]["provenance"],
        },
        "options": {
            "temperature": arguments.temperature,
            "pressure": arguments.pressure,
            "workflows": [f"{scheme}_{method}" for scheme, method in WORKFLOWS],
            "eos": arguments.eos,
            "polynomial_degree": arguments.degree,
            "quadrature_initial": [arguments.mu_order, arguments.phi_order],
            "quadrature_fine": [arguments.mu_order * 4, arguments.phi_order * 4],
            "refinement_factor": arguments.refinement_factor,
        },
        "pressure_rows": pressure_rows,
        "pressure_analysis": pressure_analysis,
        "workflow_summaries": workflow_summaries,
        "direct_composition": {
            **composition_summary,
            "rows": composition_rows,
        },
        "qha_effect_rows": effect_rows,
        "hdf5_round_trip": round_trip,
        "adversarial_checks": adversarial,
        "acceptance": acceptance,
    }

    _write_csv(output / "ohap-kieffer-pressure-sources.csv", pressure_rows)
    _write_csv(output / "ohap-kieffer-acoustic-fractions.csv", composition_rows)
    _write_csv(output / "ohap-kieffer-qha-effect.csv", effect_rows)
    json_path = output / "ohap-kieffer-validation.json"
    json_path.write_text(
        json.dumps(data, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    report_path = output / "ohap-kieffer-validation.md"
    report_path.write_text(_render_markdown(data), encoding="utf-8")
    print(report_path)
    print(json_path)
    if not all(acceptance.values()):
        failed = ", ".join(name for name, passed in acceptance.items() if not passed)
        raise SystemExit(f"OHAp Kieffer validation failed: {failed}")


if __name__ == "__main__":
    main()
