# -*- coding: utf-8 -*-

"""Generate readable quasi-static thermoelastic YAML input files."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np
import yaml
from numpy.typing import NDArray

from quantas.core.physics.elasticity import (
    assign_hydrostatic_pressures,
    detect_elastic_symmetry,
)
from quantas.core.physics.eos import (
    PressureEstimate,
    pressure_from_energy_eos,
    pressure_from_energy_polynomial,
)
from quantas.interfaces.crystal import (
    CrystalPressurePolicy,
    correct_crystal_hydrostatic_elastic_series,
    read_crystal_elastic_series,
)
from quantas.interfaces.crystal.elasticity import CrystalElasticityReader
from quantas.io.phonons import PhononInputFileReader
from quantas.models.elastic_states import (
    ElasticState,
    ElasticStateSeries,
    ElasticTensorKind,
    PressureSource,
)
from quantas.models.volume_matching import match_sampled_volumes
from quantas.references import method_citation_keys, render_citation_inline
from quantas.models.structures import CrystalStructure, SymmetryMetadata
from quantas.modules.thermoelasticity.models import (
    ElasticVolumePoint,
    ElasticVolumeSeries,
    ThermoelasticInput,
)
from quantas.modules.thermoelasticity.frames import (
    maximum_ordered_fractional_displacement,
    normalize_elastic_frame,
)


THERMOELASTIC_INPUT_SCHEMA = "1.0"


class ThermoelasticInputCreator:
    """Create a quasi-static thermoelastic input from elastic output files.

    Parameters
    ----------
    interface : str, optional
        Electronic-structure interface.  The first implementation supports
        only ``"crystal"``.
    symprec : float, optional
        Cartesian tolerance in angstrom used for structural symmetry analysis.
    angle_tolerance : float, optional
        Angular tolerance in degrees used for structural symmetry analysis.
    elastic_tolerance : float, optional
        Absolute tolerance in GPa used to detect the elastic crystal system.
    pressure_tolerance : float, optional
        Maximum permitted difference in GPa between an explicit CRYSTAL
        ``PRESSURE``/``PRESSEOS`` value and the pressure printed for the
        corrected elastic properties.
    structure_correspondence_tolerance : float, optional
        Maximum ordered-atom displacement in angstrom, after wrapping
        fractional differences into the nearest periodic image.  This guards
        against fitting unrelated structural paths while allowing internal
        relaxation along the same phase branch.
    """

    def __init__(
        self,
        *,
        interface: str = "crystal",
        symprec: float = 1.0e-5,
        angle_tolerance: float = -1.0,
        elastic_tolerance: float = 1.0e-3,
        pressure_tolerance: float = 5.0e-2,
        structure_correspondence_tolerance: float = 5.0e-1,
    ) -> None:
        if interface.strip().lower() != "crystal":
            raise ValueError("only the CRYSTAL thermoelastic interface is available")
        self.interface = "crystal"
        self.symprec = float(symprec)
        self.angle_tolerance = float(angle_tolerance)
        self.elastic_tolerance = float(elastic_tolerance)
        self.pressure_tolerance = float(pressure_tolerance)
        self.structure_correspondence_tolerance = float(
            structure_correspondence_tolerance
        )
        if (
            not np.isfinite(self.structure_correspondence_tolerance)
            or self.structure_correspondence_tolerance < 0.0
        ):
            raise ValueError(
                "structure_correspondence_tolerance must be finite and non-negative"
            )

    def create(
        self,
        sources: str | Path | Sequence[str | Path],
        *,
        jobname: str = "Quantas quasi-static thermoelastic input",
        is_list: bool = False,
        reference: int | None = None,
        pressure_source: str = "auto",
        manual_pressures_gpa: Sequence[float] | None = None,
        eos: str = "BM3",
        polynomial_degree: int = 3,
        maxfev: int | None = None,
        energy_input: str | Path | None = None,
    ) -> ThermoelasticInput:
        """Read CRYSTAL outputs and return normalized thermoelastic input.

        CRYSTAL ``PRESSURE`` and ``PRESSEOS`` outputs already contain the
        hydrostatic Barron--Klein/Wallace correction and are preserved.  Raw
        energy--strain tensors are corrected exactly once after resolving the
        hydrostatic pressure from output stress, explicit values, or an
        energy-volume relation.

        Parameters
        ----------
        sources : str, Path, or sequence
            Output path, list-file path, or explicit sequence of output paths.
        jobname : str, optional
            Human-readable description stored in the YAML file.
        is_list : bool, optional
            Interpret a scalar ``sources`` value as a text file containing one
            output path per line.
        reference : int or None, optional
            Reference index after sorting by increasing volume.  If omitted,
            the point with the smallest absolute pressure is selected.
        pressure_source : str, optional
            Pressure resolution for raw tensors: ``auto``, ``output_stress``,
            ``manual``, ``energy_eos``, or ``energy_polynomial``.  ``auto``
            preserves backend-corrected tensors and otherwise uses the final
            unstrained stress pressure.
        manual_pressures_gpa : sequence of float or None, optional
            One pressure per source output in source order for ``manual``.
        eos : str, optional
            Integrated energy EOS used by ``energy_eos``.
        polynomial_degree : int, optional
            Polynomial degree used by ``energy_polynomial``.
        maxfev : int or None, optional
            Optional nonlinear-fit iteration limit for the energy EOS.
        energy_input : str, Path, or None, optional
            Optional HA/QHA YAML whose sampled static ``E(V)`` dataset is used
            for energy-derived pressure.  Extra blocks such as Kieffer data do
            not alter the static energy-volume relation.

        Returns
        -------
        ThermoelasticInput
            Validated input contract containing only Wallace/incremental
            stiffness tensors.

        Raises
        ------
        ValueError
            If a file is invalid, a raw tensor has no defensible selected
            pressure source, or the elastic series is inconsistent.
        """
        files = self._normalize_sources(sources, is_list=is_list)
        if not files:
            raise ValueError("no CRYSTAL elastic output files were provided")

        selected_pressure = _normalize_pressure_source(pressure_source)
        elastic_states, pressure_resolution = _resolve_pressure_series(
            files,
            pressure_source=selected_pressure,
            manual_pressures_gpa=manual_pressures_gpa,
            eos=eos,
            polynomial_degree=polynomial_degree,
            maxfev=maxfev,
            energy_input=energy_input,
            symprec=self.symprec,
            angle_tolerance=self.angle_tolerance,
        )
        elastic_states.require_incremental()
        state_by_source = {
            str(Path(state.source).resolve()): state
            for state in elastic_states.states
            if state.source is not None
        }

        records: list[
            tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]
        ] = []
        for path in files:
            reader = CrystalElasticityReader(
                path,
                symprec=self.symprec,
                angle_tolerance=self.angle_tolerance,
            )
            if not reader.completed:
                raise ValueError(reader.error or f"unable to read {path}")
            state = state_by_source.get(str(path.resolve()))
            if state is None:
                raise ValueError(f"{path}: corrected elastic-state provenance is missing")
            if reader.prestress_applied and np.isfinite(reader.pressure):
                keyword = reader.prestress_keyword or "PRESSURE"
                if not np.isclose(
                    reader.pressure_keyword_value,
                    reader.pressure,
                    rtol=0.0,
                    atol=self.pressure_tolerance,
                ):
                    raise ValueError(
                        f"{path}: {keyword} value {reader.pressure_keyword_value:.6g} "
                        f"GPa differs from elastic pressure {reader.pressure:.6g} GPa"
                    )
            structure = reader.structure
            symmetry = reader.symmetry
            if structure is None or symmetry is None:
                raise ValueError(
                    f"{path}: elastic-reference structure or symmetry is unavailable"
                )
            lattice = state.lattice
            if lattice is None:
                raise ValueError(
                    f"{path}: elastic-reference lattice is incompatible with the "
                    "reported primitive-cell volume"
                )
            pressure = state.prestress.pressure_gpa
            if pressure is None or not np.isfinite(pressure):
                raise ValueError(f"{path}: corrected tensor lacks finite pressure provenance")
            stiffness = np.asarray(state.stiffness, dtype=np.float64)
            elastic_symmetry = detect_elastic_symmetry(
                stiffness,
                tolerance=self.elastic_tolerance,
            )
            point = ElasticVolumePoint(
                source=path.name,
                pressure=pressure,
                stress_pressure=reader.stress_pressure,
                volume=reader.volume,
                density=reader.density,
                energy=reader.energy,
                stiffness=stiffness,
                lattice=lattice,
                prestress_applied=True,
                metadata={
                    "prestress": _prestress_mapping(state),
                    "backend_prestress_keyword": reader.prestress_keyword,
                },
            )
            records.append((point, structure, symmetry, elastic_symmetry))

        records.sort(key=lambda item: item[0].volume)
        self._validate_records(records)
        points = tuple(item[0] for item in records)
        if reference is None:
            reference_index = int(
                np.argmin(np.abs([point.pressure for point in points]))
            )
        else:
            reference_index = int(reference)
            if reference_index < 0 or reference_index >= len(points):
                raise ValueError("reference index is outside the sorted elastic series")

        reference_structure = records[reference_index][1]
        reference_symmetry = records[reference_index][2]
        elastic_symmetry = records[reference_index][3]
        records, frame_metadata = self._normalize_frames(
            records,
            reference_index=reference_index,
        )
        points = tuple(item[0] for item in records)
        series = ElasticVolumeSeries(
            points=points,
            reference_structure=reference_structure,
            symmetry=reference_symmetry,
            elastic_symmetry=elastic_symmetry,
            reference_index=reference_index,
            orientation="crystal",
            metadata={
                "interface": self.interface,
                "ordering": "volume-ascending",
                "pressure_resolution": pressure_resolution,
                "symprec": self.symprec,
                "angle_tolerance": self.angle_tolerance,
                "elastic_tolerance_gpa": self.elastic_tolerance,
                "frame_normalization": frame_metadata,
                "structure_correspondence_tolerance_A": (
                    self.structure_correspondence_tolerance
                ),
            },
        )
        return ThermoelasticInput(
            jobname=str(jobname),
            elastic_series=series,
            method="quasistatic",
            metadata={
                "schema_name": "quantas-thermoelastic-input",
                "schema_version": THERMOELASTIC_INPUT_SCHEMA,
                "pressure_resolution": pressure_resolution,
            },
        )

    def write(
        self,
        input_data: ThermoelasticInput,
        outfile: str | Path,
    ) -> Path:
        """Write a readable thermoelastic YAML file.

        Parameters
        ----------
        input_data : ThermoelasticInput
            Input contract to serialize.
        outfile : str or Path
            Destination YAML path.

        Returns
        -------
        Path
            Written path.
        """
        path = Path(outfile)
        path.write_text(format_thermoelastic_yaml(input_data), encoding="utf-8")
        return path

    @staticmethod
    def _normalize_sources(
        sources: str | Path | Sequence[str | Path],
        *,
        is_list: bool,
    ) -> list[Path]:
        """Return explicit source paths, resolving list entries locally."""
        if is_list:
            if isinstance(sources, Sequence) and not isinstance(sources, (str, Path)):
                raise ValueError("is_list requires one list-file path")
            list_path = Path(sources)
            files: list[Path] = []
            for line in list_path.read_text(encoding="utf-8").splitlines():
                value = line.strip()
                if not value or value.startswith("#"):
                    continue
                item = Path(value)
                if not item.is_absolute():
                    item = list_path.parent / item
                files.append(item)
        elif isinstance(sources, (str, Path)):
            files = [Path(sources)]
        else:
            files = [Path(value) for value in sources]

        resolved: list[Path] = []
        seen: set[Path] = set()
        for path in files:
            normalized = path.expanduser().resolve()
            if normalized in seen:
                continue
            if not normalized.is_file():
                raise ValueError(f"elastic output file does not exist: {path}")
            seen.add(normalized)
            resolved.append(normalized)
        return resolved

    @staticmethod
    def _validate_records(
        records: list[
            tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]
        ],
    ) -> None:
        """Validate structural, crystallographic, and elastic consistency."""
        reference_point, reference_structure, reference_symmetry, elastic_symmetry = (
            records[0]
        )
        del reference_point
        for point, structure, symmetry, current_elastic_symmetry in records[1:]:
            if not np.array_equal(
                structure.atomic_numbers,
                reference_structure.atomic_numbers,
            ):
                raise ValueError(
                    f"{point.source}: primitive atomic species/order differs across inputs"
                )
            if symmetry.space_group_number != reference_symmetry.space_group_number:
                raise ValueError(
                    f"{point.source}: space group {symmetry.space_group_number} differs "
                    f"from reference {reference_symmetry.space_group_number}"
                )
            if (
                reference_symmetry.hall_number > 0
                and symmetry.hall_number > 0
                and symmetry.hall_number != reference_symmetry.hall_number
            ):
                raise ValueError(
                    f"{point.source}: Hall number {symmetry.hall_number} differs "
                    f"from reference {reference_symmetry.hall_number}"
                )
            if (
                reference_symmetry.choice
                and symmetry.choice
                and symmetry.choice != reference_symmetry.choice
            ):
                raise ValueError(
                    f"{point.source}: crystallographic setting choice "
                    f"{symmetry.choice!r} differs from reference "
                    f"{reference_symmetry.choice!r}"
                )
            if current_elastic_symmetry != elastic_symmetry:
                raise ValueError(
                    f"{point.source}: elastic symmetry {current_elastic_symmetry} differs "
                    f"from reference {elastic_symmetry}"
                )
        volumes = np.asarray([item[0].volume for item in records], dtype=np.float64)
        if np.any(np.diff(volumes) <= 1.0e-10):
            raise ValueError("elastic outputs contain duplicate or unresolved volumes")

    def _normalize_frames(
        self,
        records: list[
            tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]
        ],
        *,
        reference_index: int,
    ) -> tuple[
        list[tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]],
        dict[str, Any],
    ]:
        """Co-rotate every tensor into the selected reference Cartesian frame.

        Parameters
        ----------
        records : list
            Parsed and volume-sorted elastic records.
        reference_index : int
            Index defining the common Cartesian frame.

        Returns
        -------
        tuple
            Normalized records and series-level diagnostics.

        Raises
        ------
        ValueError
            If atom correspondence is inconsistent or a frame transformation
            is improper.
        """
        reference_structure = records[reference_index][1]
        normalized: list[
            tuple[ElasticVolumePoint, CrystalStructure, SymmetryMetadata, str]
        ] = []
        angles: list[float] = []
        displacements: list[float] = []
        for point, structure, symmetry, elastic_symmetry in records:
            displacement = maximum_ordered_fractional_displacement(
                reference_structure.fractional_positions,
                structure.fractional_positions,
                reference_structure.lattice,
            )
            if displacement > self.structure_correspondence_tolerance:
                raise ValueError(
                    f"{point.source}: maximum ordered-atom displacement "
                    f"{displacement:.6g} A exceeds the configured structural "
                    f"path tolerance {self.structure_correspondence_tolerance:.6g} A"
                )
            frame = normalize_elastic_frame(
                point.lattice,
                point.stiffness,
                reference_structure.lattice,
            )
            metadata = {
                **point.metadata,
                "frame_normalization": {
                    **frame.metadata,
                    "status": "normalized",
                    "rotation_to_reference": frame.rotation_to_reference.copy(),
                    "removed_rotation_degrees": frame.removed_rotation_degrees,
                    "principal_logarithmic_strain": (
                        frame.principal_logarithmic_strain.copy()
                    ),
                    "source_lattice": point.lattice.copy(),
                    "maximum_ordered_atom_displacement_A": displacement,
                },
            }
            normalized_point = ElasticVolumePoint(
                source=point.source,
                pressure=point.pressure,
                stress_pressure=point.stress_pressure,
                volume=point.volume,
                density=point.density,
                energy=point.energy,
                stiffness=frame.stiffness,
                lattice=frame.lattice,
                prestress_applied=point.prestress_applied,
                metadata=metadata,
            )
            normalized.append((normalized_point, structure, symmetry, elastic_symmetry))
            angles.append(frame.removed_rotation_degrees)
            displacements.append(displacement)
        return normalized, {
            "status": "normalized",
            "method": "right_polar_decomposition_corotation",
            "reference_index": int(reference_index),
            "maximum_removed_rotation_degrees": float(max(angles, default=0.0)),
            "maximum_ordered_atom_displacement_A": float(
                max(displacements, default=0.0)
            ),
            "reference": "; ".join(
                render_citation_inline(key)
                for key in method_citation_keys("wallace_stress_strain")
            ),
            "citation_keys": list(method_citation_keys("wallace_stress_strain")),
        }



_PRESSURE_SOURCES = {
    "auto",
    "output_stress",
    "manual",
    "energy_eos",
    "energy_polynomial",
}
_ENERGY_PRESSURE_SOURCES = {"energy_eos", "energy_polynomial"}


def _normalize_pressure_source(value: str) -> str:
    """Return one canonical QSA pressure-resolution identifier."""
    normalized = str(value).strip().lower().replace("-", "_")
    if normalized not in _PRESSURE_SOURCES:
        allowed = ", ".join(sorted(_PRESSURE_SOURCES))
        raise ValueError(f"unsupported pressure source {value!r}; choose one of {allowed}")
    return normalized


def _resolve_pressure_series(
    files: Sequence[Path],
    *,
    pressure_source: str,
    manual_pressures_gpa: Sequence[float] | None,
    eos: str,
    polynomial_degree: int,
    maxfev: int | None,
    energy_input: str | Path | None,
    symprec: float,
    angle_tolerance: float,
) -> tuple[ElasticStateSeries, dict[str, Any]]:
    """Return Wallace tensors and complete pressure-resolution provenance."""
    if pressure_source not in _ENERGY_PRESSURE_SOURCES:
        if energy_input is not None:
            raise ValueError(
                "energy_input is only valid with pressure_source='energy_eos' "
                "or 'energy_polynomial'"
            )
        if pressure_source != "manual" and manual_pressures_gpa is not None:
            raise ValueError(
                "manual_pressures_gpa requires pressure_source='manual'"
            )
        try:
            series = read_crystal_elastic_series(
                files,
                pressure_policy=CrystalPressurePolicy(pressure_source),
                manual_pressures_gpa=manual_pressures_gpa,
                apply_prestress_correction=True,
                correction_applied_by="quantas-thermoelastic-inpgen",
                symprec=symprec,
                angle_tolerance=angle_tolerance,
            )
        except ValueError as exc:
            if pressure_source == "auto" and "lacks pressure" in str(exc):
                raise ValueError(
                    f"{exc}. The tensor is raw and cannot enter QSA without a "
                    "hydrostatic Barron-Klein/Wallace correction. Select an "
                    "explicit pressure source: output-stress, manual, energy-eos, "
                    "or energy-polynomial."
                ) from exc
            raise
        return series, _series_pressure_resolution(
            series,
            requested_source=pressure_source,
            energy_model=None,
        )

    if manual_pressures_gpa is not None:
        raise ValueError(
            "manual_pressures_gpa cannot be combined with an energy pressure source"
        )
    raw_series = read_crystal_elastic_series(
        files,
        pressure_policy=CrystalPressurePolicy.DEFERRED,
        apply_prestress_correction=False,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    )
    estimate, model, pressures = _energy_pressure_for_elastic_series(
        raw_series,
        pressure_source=pressure_source,
        eos=eos,
        polynomial_degree=polynomial_degree,
        maxfev=maxfev,
        energy_input=energy_input,
    )
    pressure_enum = (
        PressureSource.ENERGY_EOS
        if pressure_source == "energy_eos"
        else PressureSource.ENERGY_POLYNOMIAL
    )
    assigned = assign_hydrostatic_pressures(
        raw_series,
        pressures,
        pressure_source=pressure_enum,
        assignment_method=pressure_source,
        metadata=model,
    )
    corrected = correct_crystal_hydrostatic_elastic_series(
        assigned,
        correction_applied_by="quantas-thermoelastic-inpgen",
    )
    model["evaluated_pressures_gpa"] = np.asarray(
        estimate.pressure, dtype=np.float64
    ).tolist()
    model["elastic_pressures_gpa"] = np.asarray(
        pressures, dtype=np.float64
    ).tolist()
    return corrected, _series_pressure_resolution(
        corrected,
        requested_source=pressure_source,
        energy_model=model,
    )


def _energy_pressure_for_elastic_series(
    raw_series: ElasticStateSeries,
    *,
    pressure_source: str,
    eos: str,
    polynomial_degree: int,
    maxfev: int | None,
    energy_input: str | Path | None,
) -> tuple[PressureEstimate, dict[str, Any], NDArray[np.float64]]:
    """Fit an energy-volume relation and evaluate it at elastic volumes."""
    source_dataset = "elastic_outputs_static_energy"
    volume = raw_series.volumes
    energy = np.asarray(
        [state.energy if state.energy is not None else np.nan for state in raw_series.states],
        dtype=np.float64,
    )
    energy_unit = "hartree"
    length_unit = "angstrom"
    matches = None
    if energy_input is not None:
        input_path = Path(energy_input)
        reader = PhononInputFileReader(input_path)
        if not reader.completed:
            raise ValueError(reader.error or f"unable to read energy input {input_path}")
        phonon_input = reader.to_input(source=input_path)
        if phonon_input.volume is None or phonon_input.energy is None:
            raise ValueError("energy input does not contain sampled static volume-energy data")
        volume = np.asarray(phonon_input.volume, dtype=np.float64)
        energy = np.asarray(phonon_input.energy, dtype=np.float64)
        energy_unit = str(phonon_input.units.get("energy", "Ha"))
        length_unit = str(phonon_input.units.get("length", "angstrom"))
        source_dataset = str(input_path)
        matches = match_sampled_volumes(raw_series.volumes, volume)
    if volume.size < 3 or energy.shape != volume.shape or not np.all(np.isfinite(energy)):
        raise ValueError(
            "energy-derived pressure requires at least three finite aligned "
            "volume-energy points; select output-stress or manual pressure instead"
        )
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
    if matches is None:
        pressures = np.asarray(estimate.pressure, dtype=np.float64)
        if pressures.shape != (raw_series.nstates,):
            raise ValueError("energy-pressure fit is not aligned with elastic volumes")
        volume_matches: list[dict[str, Any]] = []
    else:
        pressures = np.asarray(
            [estimate.pressure[match.source_index] for match in matches],
            dtype=np.float64,
        )
        volume_matches = [
            {
                "elastic_index": match.target_index,
                "energy_index": match.source_index,
                "elastic_volume": match.target_volume,
                "energy_volume": match.source_volume,
                "absolute_difference": match.absolute_difference,
                "relative_difference": match.relative_difference,
            }
            for match in matches
        ]
    model = _plain_data(
        {
            "method": pressure_source,
            "relation": "P(V) = -dE/dV",
            "source_dataset": source_dataset,
            "energy_unit": energy_unit,
            "volume_length_unit": length_unit,
            "pressure_unit": estimate.unit,
            "settings": dict(estimate.metadata),
            "fit": estimate.fit.as_dict(),
            "warnings": list(estimate.warnings),
            "volume_matches": volume_matches,
        }
    )
    return estimate, model, pressures


def _prestress_mapping(state: ElasticState) -> dict[str, Any]:
    """Return serialization-ready provenance for one corrected elastic state."""
    prestress = state.prestress
    tensor_kind = ElasticTensorKind(prestress.tensor_kind)
    pressure_source = PressureSource(prestress.pressure_source)
    source_tensor_kind = (
        None
        if prestress.source_tensor_kind is None
        else ElasticTensorKind(prestress.source_tensor_kind)
    )
    return {
        "tensor_kind": tensor_kind.value,
        "pressure_gpa": prestress.pressure_gpa,
        "pressure_source": pressure_source.value,
        "correction_method": prestress.correction_method,
        "correction_applied_by": prestress.correction_applied_by,
        "source_tensor_kind": (
            None if source_tensor_kind is None else source_tensor_kind.value
        ),
    }


def _series_pressure_resolution(
    series: ElasticStateSeries,
    *,
    requested_source: str,
    energy_model: dict[str, Any] | None,
) -> dict[str, Any]:
    """Summarize how every source tensor became QSA-ready."""
    states = []
    for state in series.states:
        entry = _prestress_mapping(state)
        entry.update({"source": state.source, "volume": state.volume})
        keyword = state.metadata.get("prestress_keyword")
        if keyword:
            entry["backend_keyword"] = keyword
        states.append(entry)
    payload: dict[str, Any] = {
        "interface": "crystal",
        "requested_source": requested_source,
        "target_tensor_kind": "wallace_hydrostatic",
        "correction_policy": "preserve-backend-or-apply-once",
        "correction_formulation": "crystal-erba-2014-hydrostatic",
        "correction_reference_doi": "10.1063/1.4869144",
        "states": states,
    }
    if energy_model is not None:
        payload["energy_model"] = energy_model
    return _plain_data(payload)


def _plain_data(value: Any) -> Any:
    """Recursively convert NumPy scalars and arrays to YAML/HDF5-safe values."""
    if isinstance(value, np.ndarray):
        return [_plain_data(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _plain_data(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain_data(item) for item in value]
    return value

def create_thermoelastic_input(
    sources: str | Path | Sequence[str | Path],
    outfile: str | Path,
    *,
    interface: str = "crystal",
    is_list: bool = False,
    jobname: str = "Quantas quasi-static thermoelastic input",
    reference: int | None = None,
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
    elastic_tolerance: float = 1.0e-3,
    pressure_tolerance: float = 5.0e-2,
    structure_correspondence_tolerance: float = 5.0e-1,
    pressure_source: str = "auto",
    manual_pressures_gpa: Sequence[float] | None = None,
    eos: str = "BM3",
    polynomial_degree: int = 3,
    maxfev: int | None = None,
    energy_input: str | Path | None = None,
) -> Path:
    """Create a quasi-static thermoelastic YAML input through the Python API.

    Parameters
    ----------
    sources : str, Path, or sequence
        CRYSTAL output path, list-file path, or explicit output sequence.
    outfile : str or Path
        Destination YAML file.
    interface : str, optional
        Electronic-structure interface.  Currently ``"crystal"`` only.
    is_list : bool, optional
        Interpret a scalar source as a list file.
    jobname : str, optional
        Input description.
    reference : int or None, optional
        Reference index after volume sorting.
    symprec : float, optional
        spglib Cartesian tolerance in angstrom.
    angle_tolerance : float, optional
        spglib angular tolerance in degrees.
    elastic_tolerance : float, optional
        Elastic symmetry tolerance in GPa.
    pressure_tolerance : float, optional
        Maximum permitted difference, in GPa, between CRYSTAL
        ``PRESSURE``/``PRESSEOS`` and the pressure reported for the corrected
        elastic tensor.
    structure_correspondence_tolerance : float, optional
        Maximum ordered-atom displacement in angstrom along the structural
        path.
    pressure_source : str, optional
        Pressure source for raw elastic tensors.
    manual_pressures_gpa : sequence of float or None, optional
        Explicit pressures in input-file order for the manual policy.
    eos : str, optional
        Integrated energy EOS used by the energy-EOS policy.
    polynomial_degree : int, optional
        Polynomial degree used by the energy-polynomial policy.
    maxfev : int or None, optional
        Optional EOS fitting iteration limit.
    energy_input : str, Path, or None, optional
        Optional HA/QHA YAML providing the static energy-volume series.

    Returns
    -------
    Path
        Written YAML path.
    """
    creator = ThermoelasticInputCreator(
        interface=interface,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        elastic_tolerance=elastic_tolerance,
        pressure_tolerance=pressure_tolerance,
        structure_correspondence_tolerance=structure_correspondence_tolerance,
    )
    input_data = creator.create(
        sources,
        jobname=jobname,
        is_list=is_list,
        reference=reference,
        pressure_source=pressure_source,
        manual_pressures_gpa=manual_pressures_gpa,
        eos=eos,
        polynomial_degree=polynomial_degree,
        maxfev=maxfev,
        energy_input=energy_input,
    )
    return creator.write(input_data, outfile)


def format_thermoelastic_yaml(input_data: ThermoelasticInput) -> str:
    """Serialize thermoelastic input with compact row-wise vectors and matrices.

    Parameters
    ----------
    input_data : ThermoelasticInput
        Validated thermoelastic input contract.

    Returns
    -------
    str
        Readable YAML text.
    """
    series = input_data.elastic_series
    structure = series.reference_structure
    symmetry = series.symmetry
    lines: list[str] = [
        "schema:",
        "  name: quantas-thermoelastic-input",
        f"  version: '{THERMOELASTIC_INPUT_SCHEMA}'",
        f"job: {_yaml_string(input_data.jobname)}",
        f"method: {input_data.method}",
        "interface: crystal",
        "conventions:",
        "  strain: eulerian-finite-strain",
        "  stiffness: wallace-hydrostatic-stress-strain",
        "  prestress: wallace-hydrostatic-resolved-by-interface",
        "  tensor_orientation: crystal",
        "  voigt_order: [ 11, 22, 33, 23, 13, 12 ]",
        "units:",
        "  pressure: GPa",
        "  stiffness: GPa",
        "  volume: angstrom^3",
        "  density: kg/m^3",
        "  energy: hartree",
        "reference:",
        f"  index: {series.reference_index}",
        f"  source: {_yaml_string(str(series.points[series.reference_index].source))}",
        "  structure:",
        f"    natom: {structure.natoms}",
        "    lattice:",
    ]
    lines.extend(_matrix_lines(structure.lattice, indent="    ", precision=12))
    atomic_numbers = ", ".join(str(int(value)) for value in structure.atomic_numbers)
    lines.append(f"    atomic_numbers: [ {atomic_numbers} ]")
    lines.append("    fractional_positions:")
    lines.extend(
        _matrix_lines(
            structure.fractional_positions,
            indent="    ",
            precision=12,
        )
    )
    pressure_resolution = series.metadata.get("pressure_resolution")
    if isinstance(pressure_resolution, dict):
        reference_line = lines.index("reference:")
        resolution_lines = yaml.safe_dump(
            {"pressure_resolution": _plain_data(pressure_resolution)},
            sort_keys=False,
            default_flow_style=False,
        ).rstrip().splitlines()
        lines[reference_line:reference_line] = resolution_lines

    lines.extend(
        [
            "  symmetry:",
            f"    space_group_number: {symmetry.space_group_number}",
            f"    international_symbol: {_yaml_string(symmetry.international_symbol)}",
            f"    hall_number: {symmetry.hall_number}",
            f"    hall_symbol: {_yaml_string(symmetry.hall_symbol)}",
            f"    choice: {_yaml_string(symmetry.choice)}",
            f"    point_group: {_yaml_string(symmetry.point_group)}",
            f"    elastic_system: {series.elastic_symmetry}",
            f"    symprec: {symmetry.symprec:.6E}",
            f"    angle_tolerance: {symmetry.angle_tolerance:.6f}",
            "elastic_data:",
        ]
    )
    for point in series.points:
        lines.extend(
            [
                f"- source: {_yaml_string(str(point.source))}",
                f"  pressure: {point.pressure: .8f}",
                (
                    "  stress_pressure: null"
                    if not np.isfinite(point.stress_pressure)
                    else f"  stress_pressure: {point.stress_pressure: .8f}"
                ),
                f"  volume: {point.volume: .12f}",
                f"  density: {point.density: .8f}",
                f"  energy: {point.energy: .12E}",
                "  lattice:",
            ]
        )
        lines.extend(_matrix_lines(point.lattice, indent="  ", precision=12))
        lines.append("  stiffness:")
        lines.extend(_matrix_lines(point.stiffness, indent="  ", precision=8))
        frame = point.metadata.get("frame_normalization")
        if isinstance(frame, dict):
            lines.extend(
                [
                    "  frame:",
                    f"    status: {_yaml_string(str(frame.get('status', 'unknown')))}",
                    f"    method: {_yaml_string(str(frame.get('method', 'unknown')))}",
                    "    removed_rotation_degrees: "
                    f"{float(frame.get('removed_rotation_degrees', 0.0)): .12E}",
                    "    maximum_ordered_atom_displacement_A: "
                    f"{float(frame.get('maximum_ordered_atom_displacement_A', 0.0)): .12E}",
                    "    rotation_to_reference:",
                ]
            )
            lines.extend(
                _matrix_lines(
                    frame.get("rotation_to_reference", np.eye(3)),
                    indent="    ",
                    precision=12,
                )
            )
            principal = np.asarray(
                frame.get("principal_logarithmic_strain", np.zeros(3)),
                dtype=np.float64,
            )
            principal_text = ", ".join(f"{value: .12E}" for value in principal)
            lines.append(f"    principal_logarithmic_strain: [ {principal_text} ]")
            lines.append("    source_lattice:")
            lines.extend(
                _matrix_lines(
                    frame.get("source_lattice", point.lattice),
                    indent="    ",
                    precision=12,
                )
            )
    return "\n".join(lines) + "\n"


def _matrix_lines(
    matrix: Any,
    *,
    indent: str,
    precision: int,
) -> list[str]:
    """Format a two-dimensional array with one compact row per YAML line."""
    array = np.asarray(matrix, dtype=np.float64)
    rows: list[str] = []
    threshold = 0.5 * 10.0 ** (-precision)
    for row in array:
        cleaned = [
            0.0 if abs(float(value)) < threshold else float(value) for value in row
        ]
        values = ", ".join(f"{value: .{precision}f}" for value in cleaned)
        rows.append(f"{indent}- [ {values} ]")
    return rows


def _yaml_string(value: str) -> str:
    """Return a single-quoted YAML scalar with escaped apostrophes."""
    return "'" + str(value).replace("'", "''") + "'"


__all__ = [
    "THERMOELASTIC_INPUT_SCHEMA",
    "ThermoelasticInputCreator",
    "create_thermoelastic_input",
    "format_thermoelastic_yaml",
]
