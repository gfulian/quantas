# -*- coding: utf-8 -*-

"""CRYSTAL structure--energy parsing for static Energy EOS workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from quantas.interfaces.crystal.geometry import CrystalGeometryParser
from quantas.interfaces.crystal.output import CrystalOutputDocument, CrystalOutputParser
from quantas.models.computation import (
    EnergyKind,
    EnergyRecord,
    OptimizationStatus,
    RunTerminationStatus,
    SourceProvenance,
    StructureEnergyPoint,
    StructureEnergySeries,
)
from quantas.models.structures import CrystalStructure


_EOS_SORTING_MARKER = "SORTING VOLUMES/ENERGIES"
_EOS_FITTING_MARKER = "FITTING USING"
_NATIVE_VOLUME_ATOL = 5.0e-5
_NATIVE_ENERGY_ATOL = 5.0e-9


@dataclass(frozen=True, slots=True)
class CrystalEnergyVolumeParseResult:
    """Parsed structure--energy states from one CRYSTAL output.

    Parameters
    ----------
    series : StructureEnergySeries
        Backend-neutral structure--energy observations.
    run_kind : str
        Classified CRYSTAL run kind: ``"static"``, ``"optimization"``, or
        ``"native_eos"``.
    corrections : tuple of str
        Canonical correction signature shared by all returned states.
    """

    series: StructureEnergySeries
    run_kind: str
    corrections: tuple[str, ...]


class CrystalEnergyVolumeReader:
    """Read one CRYSTAL output as one or more Energy EOS observations.

    The reader classifies the output from its semantic markers. A static SCF
    contributes one point, a normal geometry optimization contributes its one
    final optimized state, and a native CRYSTAL ``EOS`` calculation contributes
    every state in the final ``SORTING VOLUMES/ENERGIES`` table.

    Parameters
    ----------
    source : str, Path, or iterable of str
        CRYSTAL output path or pre-read output lines.

    Raises
    ------
    ValueError
        If the output is incomplete, ambiguous, structurally inconsistent, or
        cannot associate one authoritative total energy with each state.
    """

    def __init__(self, source: str | Path | Iterable[str]) -> None:
        """Initialize the reader and retain one shared in-memory document."""
        self.document = CrystalOutputDocument(source)
        self.output = CrystalOutputParser(self.document)
        self.geometry = CrystalGeometryParser(self.document.lines)
        self.source = self.document.source

    def read(self) -> CrystalEnergyVolumeParseResult:
        """Return normalized Energy EOS states from the CRYSTAL output.

        Returns
        -------
        CrystalEnergyVolumeParseResult
            Parsed series, run classification, and energy-correction signature.

        Raises
        ------
        ValueError
            If CRYSTAL did not terminate normally or the scientific state cannot
            be reconstructed unambiguously.
        """
        termination = self.output.termination()
        if termination.status is not RunTerminationStatus.NORMAL:
            raise ValueError("CRYSTAL output did not terminate normally")

        summary = self._native_eos_summary()
        optimized = self.geometry.optimized_cells()
        if summary:
            points = self._read_native_eos(summary, optimized)
            run_kind = "native_eos"
        elif optimized:
            points = self._read_single_optimization(optimized)
            run_kind = "optimization"
        else:
            points = self._read_static()
            run_kind = "static"

        corrections = self._validate_correction_signature(points)
        ordered = tuple(sorted(points, key=_point_volume))
        series = StructureEnergySeries(
            points=ordered,
            reference_index=min(
                range(len(ordered)),
                key=lambda index: ordered[index].energy.value,
            ),
            metadata={
                "interface": "crystal",
                "run_kind": run_kind,
                "corrections": corrections,
                "source": None if self.source is None else str(self.source),
            },
        )
        return CrystalEnergyVolumeParseResult(
            series=series,
            run_kind=run_kind,
            corrections=corrections,
        )

    def _read_static(self) -> list[StructureEnergyPoint]:
        """Return one static SCF structure--energy state."""
        totals = self.output.total_energies()
        if len(totals) != 1:
            raise ValueError(
                "static CRYSTAL Energy EOS input requires exactly one resolved "
                f"SCF/total-energy state; found {len(totals)}"
            )
        structure = self._static_structure()
        return [self._point(structure, totals[0], record_index=0, state_kind="static")]

    def _read_single_optimization(
        self,
        optimized: list[CrystalStructure],
    ) -> list[StructureEnergyPoint]:
        """Return the final state of one ordinary geometry optimization."""
        if len(optimized) != 1:
            raise ValueError(
                "multiple FINAL OPTIMIZED GEOMETRY blocks require a completed "
                "native CRYSTAL EOS summary"
            )
        optimizations = self.output.optimizations()
        if not optimizations:
            raise ValueError("optimized CRYSTAL geometry has no optimization record")
        if any(item.status is not OptimizationStatus.CONVERGED for item in optimizations):
            raise ValueError("CRYSTAL geometry optimization did not converge")

        marker_lines = self._final_geometry_marker_lines()
        marker_line = marker_lines[0]
        energy = self._first_total_after(marker_line, stop_line=len(self.document.lines))
        if energy is None:
            corrected = self.output.corrected_total_energies()
            if corrected:
                raise ValueError(
                    "corrected CRYSTAL energies are present but none can be "
                    "associated with the final optimized geometry"
                )
            final = optimizations[-1].final_energy
            if final is None:
                raise ValueError("final optimized CRYSTAL energy is unavailable")
            energy = EnergyRecord(
                value=final.value,
                unit=final.unit,
                kind=EnergyKind.TOTAL,
                metadata={
                    **final.metadata,
                    "corrections": (),
                    "energy_resolution": "optimization_end_fallback",
                },
            )

        return [
            self._point(
                optimized[0],
                energy,
                record_index=0,
                state_kind="optimization",
            )
        ]

    def _read_native_eos(
        self,
        summary: list[tuple[float, float]],
        optimized: list[CrystalStructure],
    ) -> list[StructureEnergyPoint]:
        """Return all states from one completed native CRYSTAL EOS run."""
        if len(optimized) != len(summary):
            raise ValueError(
                "CRYSTAL native EOS state count is inconsistent: "
                f"{len(summary)} summary points but {len(optimized)} final geometries"
            )
        marker_lines = self._final_geometry_marker_lines()
        if len(marker_lines) != len(optimized):
            raise ValueError("unable to locate every CRYSTAL final-geometry marker")
        summary_line = self._native_eos_summary_line()
        if summary_line is None:
            raise ValueError("CRYSTAL native EOS summary marker is unavailable")

        source_states: list[tuple[CrystalStructure, EnergyRecord]] = []
        for index, (structure, marker_line) in enumerate(zip(optimized, marker_lines, strict=True)):
            stop_line = (
                marker_lines[index + 1]
                if index + 1 < len(marker_lines)
                else summary_line
            )
            energy = self._first_total_after(marker_line, stop_line=stop_line)
            if energy is None:
                raise ValueError(
                    "unable to associate a resolved total energy with CRYSTAL "
                    f"native EOS optimized state {index}"
                )
            source_states.append((structure, energy))

        points: list[StructureEnergyPoint] = []
        unused = set(range(len(source_states)))
        for record_index, (summary_volume, summary_energy) in enumerate(summary):
            matches = [
                index
                for index in unused
                if np.isclose(
                    _structure_reported_volume(source_states[index][0]),
                    summary_volume,
                    rtol=0.0,
                    atol=_NATIVE_VOLUME_ATOL,
                )
            ]
            if len(matches) != 1:
                raise ValueError(
                    "CRYSTAL native EOS summary volume cannot be matched uniquely "
                    f"to a final optimized geometry: {summary_volume:.8g} Å^3"
                )
            source_index = matches[0]
            unused.remove(source_index)
            structure, energy = source_states[source_index]
            if not np.isclose(
                energy.value,
                summary_energy,
                rtol=1.0e-12,
                atol=_NATIVE_ENERGY_ATOL,
            ):
                raise ValueError(
                    "CRYSTAL native EOS summary energy disagrees with the "
                    "state-resolved total energy at "
                    f"V={summary_volume:.8g} Å^3"
                )
            point = self._point(
                structure,
                energy,
                record_index=record_index,
                state_kind="native_eos",
                reported_volume=summary_volume,
            )
            point.metadata.update(
                {
                    "native_eos_summary_energy": float(summary_energy),
                    "native_eos_source_state_index": int(source_index),
                }
            )
            points.append(point)

        if unused:
            raise ValueError("unmatched CRYSTAL native EOS optimized geometries remain")
        return points

    def _native_eos_summary(self) -> list[tuple[float, float]]:
        """Parse the final CRYSTAL sorted volume--energy table when present."""
        marker = self._native_eos_summary_line()
        if marker is None:
            return []
        rows: list[tuple[float, float]] = []
        for line in self.document.lines[marker + 1 :]:
            stripped = line.strip()
            if rows and _EOS_FITTING_MARKER in stripped.upper():
                break
            tokens = stripped.split()
            if len(tokens) != 2:
                continue
            try:
                volume = float(tokens[0].replace("D", "E").replace("d", "e"))
                energy = float(tokens[1].replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            if volume > 0.0 and np.isfinite(volume) and np.isfinite(energy):
                rows.append((volume, energy))
        if not rows:
            raise ValueError("CRYSTAL EOS sorting marker has no volume--energy rows")
        return rows

    def _native_eos_summary_line(self) -> int | None:
        """Return the final native-EOS sorting-table marker line."""
        for index in range(len(self.document.lines) - 1, -1, -1):
            if _EOS_SORTING_MARKER in self.document.lines[index].upper():
                return index
        return None

    def _final_geometry_marker_lines(self) -> list[int]:
        """Return source lines of all final optimized geometry markers."""
        return [
            index
            for index, line in enumerate(self.document.lines)
            if "FINAL OPTIMIZED GEOMETRY" in line.upper()
        ]

    def _first_total_after(self, marker_line: int, *, stop_line: int) -> EnergyRecord | None:
        """Return the first resolved total energy after a final geometry marker."""
        candidates = [
            record
            for record in self.output.total_energies()
            if marker_line < int(record.metadata.get("line_index", -1)) < stop_line
        ]
        if not candidates:
            return None
        return min(candidates, key=lambda item: int(item.metadata.get("line_index", -1)))

    def _static_structure(self) -> CrystalStructure:
        """Return the wave-function cell for a static CRYSTAL calculation."""
        try:
            return self.geometry.wavefunction_cell()
        except ValueError:
            return self.geometry.initial_primitive_cell()

    def _point(
        self,
        structure: CrystalStructure,
        energy: EnergyRecord,
        *,
        record_index: int,
        state_kind: str,
        reported_volume: float | None = None,
    ) -> StructureEnergyPoint:
        """Build one backend-neutral point with CRYSTAL provenance."""
        volume = (
            _structure_reported_volume(structure)
            if reported_volume is None
            else float(reported_volume)
        )
        source = None if self.source is None else str(self.source)
        return StructureEnergyPoint(
            structure=structure,
            energy=energy,
            provenance=SourceProvenance(
                interface="crystal",
                source=source,
                record_index=record_index,
                metadata={"run_kind": state_kind},
            ),
            metadata={
                "run_kind": state_kind,
                "reported_volume_angstrom3": volume,
                "corrections": _correction_signature(energy),
            },
        )

    @staticmethod
    def _validate_correction_signature(
        points: list[StructureEnergyPoint],
    ) -> tuple[str, ...]:
        """Return the common correction signature or reject mixed states."""
        signatures = {_correction_signature(point.energy) for point in points}
        if len(signatures) != 1:
            rendered = ", ".join("+".join(item) if item else "DFT" for item in sorted(signatures))
            raise ValueError(
                "CRYSTAL Energy EOS states use inconsistent energy corrections: "
                f"{rendered}"
            )
        return next(iter(signatures))


def read_crystal_energy_volume(
    source: str | Path | Iterable[str],
) -> CrystalEnergyVolumeParseResult:
    """Read one CRYSTAL output into backend-neutral Energy EOS states.

    Parameters
    ----------
    source : str, Path, or iterable of str
        CRYSTAL output source.

    Returns
    -------
    CrystalEnergyVolumeParseResult
        Parsed state series and CRYSTAL run classification.
    """
    return CrystalEnergyVolumeReader(source).read()


def _correction_signature(energy: EnergyRecord) -> tuple[str, ...]:
    """Return a canonical order-insensitive correction signature."""
    return tuple(sorted(str(item).strip().upper() for item in energy.metadata.get("corrections", ())))


def _structure_reported_volume(structure: CrystalStructure) -> float:
    """Return the backend-reported volume when available."""
    value = structure.metadata.get("reported_volume_angstrom3")
    return structure.volume if value is None else float(value)


def _point_volume(point: StructureEnergyPoint) -> float:
    """Return the volume used to order one parsed Energy EOS point."""
    value = point.metadata.get("reported_volume_angstrom3")
    return point.volume if value is None else float(value)


__all__ = [
    "CrystalEnergyVolumeParseResult",
    "CrystalEnergyVolumeReader",
    "read_crystal_energy_volume",
]
