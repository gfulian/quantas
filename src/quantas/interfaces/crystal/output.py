# -*- coding: utf-8 -*-

"""Backend-specific parsing of generic CRYSTAL run and convergence records."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

from quantas.interfaces.crystal import patterns
from quantas.models.computation import (
    EnergyKind,
    EnergyRecord,
    OptimizationResult,
    OptimizationStatus,
    OptimizationStep,
    RunTermination,
    RunTerminationStatus,
    SCFResult,
    SCFStatus,
)


def _as_float(value: str) -> float:
    """Convert a CRYSTAL/Fortran numeric literal to Python ``float``."""
    return float(value.replace("D", "E").replace("d", "e"))


class CrystalOutputDocument:
    """In-memory text representation of one CRYSTAL output.

    Parameters
    ----------
    source : str or Path or iterable of str
        CRYSTAL output path or pre-read output lines.
    """

    def __init__(self, source: str | Path | Iterable[str]) -> None:
        """Read and normalize CRYSTAL output lines exactly once."""
        self.source: Path | None
        if isinstance(source, (str, Path)):
            self.source = Path(source)
            self.lines = self.source.read_text(
                encoding="utf-8",
                errors="strict",
            ).splitlines()
        else:
            self.source = None
            self.lines = [str(line).rstrip("\n") for line in source]

    def __len__(self) -> int:
        """Return the number of source lines.

        Returns
        -------
        int
            Number of lines in the output document.
        """
        return len(self.lines)

    def __iter__(self) -> Iterator[str]:
        """Iterate over normalized output lines.

        Returns
        -------
        iterator of str
            Iterator over source lines in original order.
        """
        return iter(self.lines)


class CrystalOutputParser:
    """Parse backend-generic run information from a CRYSTAL output.

    The parser deliberately does not decide whether a parsed result is
    scientifically acceptable for HA, QHA, elasticity, or any other Quantas
    workflow. It reports facts printed by CRYSTAL and leaves workflow policy to
    the specialized reader or calculator.

    Parameters
    ----------
    source : CrystalOutputDocument or str or Path or iterable of str
        Output document or source used to construct one.
    """

    def __init__(
        self,
        source: CrystalOutputDocument | str | Path | Iterable[str],
    ) -> None:
        """Initialize a generic parser from one CRYSTAL output document."""
        if isinstance(source, CrystalOutputDocument):
            self.document = source
        else:
            self.document = CrystalOutputDocument(source)

    @property
    def lines(self) -> list[str]:
        """Return the normalized output lines.

        Returns
        -------
        list of str
            Source lines in original order.
        """
        return self.document.lines

    def termination(self) -> RunTermination:
        """Return the observed global CRYSTAL termination state.

        Returns
        -------
        RunTermination
            ``NORMAL`` when the standard CRYSTAL termination banner is found;
            otherwise ``INCOMPLETE`` for a non-empty output and ``UNKNOWN`` for
            an empty document.
        """
        for index in range(len(self.lines) - 1, -1, -1):
            if patterns.NORMAL_TERMINATION_RE.search(self.lines[index]):
                return RunTermination(
                    status=RunTerminationStatus.NORMAL,
                    line_index=index,
                    metadata={"source_marker": "EEEEEEEEEE TERMINATION"},
                )
        status = (
            RunTerminationStatus.INCOMPLETE
            if self.lines
            else RunTerminationStatus.UNKNOWN
        )
        return RunTermination(status=status)

    def scf_energies(self) -> tuple[EnergyRecord, ...]:
        """Return energies printed on CRYSTAL ``SCF ENDED`` records.

        These values represent the converged electronic SCF energy before any
        optional a-posteriori correction such as DFT-D or gCP.  The parser
        reports every record in source order and leaves convergence policy to
        the consuming workflow.

        Returns
        -------
        tuple of EnergyRecord
            SCF energies in hartree.  The record kind is
            :attr:`~quantas.models.computation.EnergyKind.DFT` for historical
            compatibility with Quantas' backend-neutral energy taxonomy.
        """
        records: list[EnergyRecord] = []
        for index, line in enumerate(self.lines):
            end_match = patterns.SCF_END_RE.search(line)
            energy_match = patterns.SCF_END_ENERGY_RE.search(line)
            if end_match is None or energy_match is None:
                continue
            reason = end_match.group("reason").strip()
            records.append(
                EnergyRecord(
                    value=_as_float(energy_match.group("energy")),
                    unit="Ha",
                    kind=EnergyKind.DFT,
                    metadata={
                        "source_marker": "SCF ENDED",
                        "line_index": index,
                        "reported_cycles": int(energy_match.group("cycles")),
                        "end_reason": reason,
                        "converged": "CONVERGENCE" in reason.upper(),
                    },
                )
            )
        return tuple(records)

    def dft_energies(self) -> tuple[EnergyRecord, ...]:
        """Return all ``TOTAL ENERGY(DFT)(AU)`` records in source order.

        Returns
        -------
        tuple of EnergyRecord
            DFT energies in hartree. No convergence policy is applied.
        """
        records: list[EnergyRecord] = []
        for index, line in enumerate(self.lines):
            match = patterns.TOTAL_DFT_ENERGY_RE.search(line)
            if match is None:
                continue
            metadata: dict[str, object] = {
                "source_marker": "TOTAL ENERGY(DFT)(AU)",
                "line_index": index,
                "cycles": int(match.group("cycles")),
            }
            delta_match = patterns.TOTAL_DFT_DELTA_RE.search(line)
            if delta_match is not None:
                metadata["delta_energy"] = _as_float(delta_match.group("delta"))
            records.append(
                EnergyRecord(
                    value=_as_float(match.group("energy")),
                    unit="Ha",
                    kind=EnergyKind.DFT,
                    metadata=metadata,
                )
            )
        return tuple(records)

    def corrected_total_energies(self) -> tuple[EnergyRecord, ...]:
        """Return explicitly corrected ``TOTAL ENERGY + ...`` records.

        Returns
        -------
        tuple of EnergyRecord
            Corrected total energies in hartree. An empty tuple is returned
            when CRYSTAL did not print this optional form.
        """
        records: list[EnergyRecord] = []
        for index, line in enumerate(self.lines):
            match = patterns.CORRECTED_TOTAL_ENERGY_RE.search(line)
            if match is None:
                continue
            marker = " ".join(match.group("marker").split())
            corrections = tuple(
                term.strip().upper()
                for term in match.group("corrections").split("+")
                if term.strip()
            )
            records.append(
                EnergyRecord(
                    value=_as_float(match.group("energy")),
                    unit="Ha",
                    kind=EnergyKind.TOTAL,
                    metadata={
                        "source_marker": marker,
                        "line_index": index,
                        "corrections": corrections,
                    },
                )
            )
        return tuple(records)

    def total_energies(self) -> tuple[EnergyRecord, ...]:
        """Resolve the physical total energy for each CRYSTAL SCF state.

        CRYSTAL always reports the electronic SCF energy and may subsequently
        print a corrected total such as ``TOTAL ENERGY + DISP (AU)``,
        ``TOTAL ENERGY + GCP (AU)``, or
        ``TOTAL ENERGY + DISP + GCP (AU)``.  Quantas treats the most complete
        corrected total printed for the same SCF state as authoritative.  If
        no corrected total is present, the SCF/DFT energy is also the total
        energy.

        State association is performed in source order.  ``SCF ENDED`` is the
        preferred state anchor because it establishes the electronic energy of
        one completed SCF calculation.  ``TOTAL ENERGY(DFT)(AU)`` and optional
        corrected totals printed after that marker are considered only before
        the next SCF state, preventing a correction from being attached to a
        later calculation.  Older outputs without an ``SCF ENDED`` energy fall
        back to ``TOTAL ENERGY(DFT)(AU)`` as the state anchor.

        Returns
        -------
        tuple of EnergyRecord
            One total-energy record per resolved SCF state, in source order.
            Metadata retains the corresponding SCF energy, correction labels,
            and source markers.
        """
        dft_records = list(self.dft_energies())
        scf_records = list(self.scf_energies())
        corrected_records = list(self.corrected_total_energies())

        anchors = scf_records if scf_records else dft_records
        if not anchors:
            return ()

        totals: list[EnergyRecord] = []
        for position, anchor in enumerate(anchors):
            anchor_line = int(anchor.metadata.get("line_index", -1))
            next_anchor_line = (
                int(anchors[position + 1].metadata.get("line_index", len(self.lines)))
                if position + 1 < len(anchors)
                else len(self.lines)
            )

            associated_scf = anchor
            dft_candidates = [
                record
                for record in dft_records
                if anchor_line <= int(record.metadata.get("line_index", -1)) < next_anchor_line
            ]
            dft = max(
                dft_candidates,
                key=lambda record: int(record.metadata.get("line_index", -1)),
                default=None,
            )
            base_energy = dft if dft is not None else associated_scf

            candidates = [
                record
                for record in corrected_records
                if anchor_line < int(record.metadata.get("line_index", -1)) < next_anchor_line
            ]
            corrected = max(
                candidates,
                key=lambda record: (
                    len(tuple(record.metadata.get("corrections", ()))),
                    int(record.metadata.get("line_index", -1)),
                ),
                default=None,
            )

            if corrected is None:
                value = float(base_energy.value)
                source_marker = str(
                    base_energy.metadata.get("source_marker", "SCF ENDED")
                )
                source_line = int(base_energy.metadata.get("line_index", anchor_line))
                corrections: tuple[str, ...] = ()
            else:
                value = float(corrected.value)
                source_marker = str(corrected.metadata.get("source_marker", "TOTAL ENERGY +"))
                source_line = int(corrected.metadata.get("line_index", anchor_line))
                corrections = tuple(corrected.metadata.get("corrections", ()))

            scf_value = float(associated_scf.value)
            totals.append(
                EnergyRecord(
                    value=value,
                    unit="Ha",
                    kind=EnergyKind.TOTAL,
                    metadata={
                        "source_marker": source_marker,
                        "line_index": source_line,
                        "corrections": corrections,
                        "scf_energy": scf_value,
                        "scf_energy_unit": "Ha",
                        "scf_source_marker": associated_scf.metadata.get(
                            "source_marker", "SCF ENDED"
                        ),
                        "scf_line_index": associated_scf.metadata.get("line_index"),
                        "total_correction_energy": value - scf_value,
                    },
                )
            )
        return tuple(totals)

    def reference_energies(self) -> tuple[EnergyRecord, ...]:
        """Return CRYSTAL ``CENTRAL POINT`` reference energies.

        Returns
        -------
        tuple of EnergyRecord
            Reference energies in hartree, in source order.
        """
        records: list[EnergyRecord] = []
        for index, line in enumerate(self.lines):
            match = patterns.CENTRAL_POINT_RE.search(line)
            if match is None:
                continue
            records.append(
                EnergyRecord(
                    value=_as_float(match.group("energy")),
                    unit="Ha",
                    kind=EnergyKind.REFERENCE,
                    metadata={
                        "source_marker": "CENTRAL POINT",
                        "line_index": index,
                    },
                )
            )
        return tuple(records)

    def scf_results(self) -> tuple[SCFResult, ...]:
        """Return every SCF convergence block printed in the output.

        Returns
        -------
        tuple of SCFResult
            SCF blocks in source order. A trailing block without an ``SCF
            ENDED`` marker is retained with status ``INCOMPLETE``.
        """
        results: list[SCFResult] = []
        cycle_numbers: list[int] = []
        energies: list[float] = []
        delta_energies: list[float] = []
        block_start: int | None = None

        for index, line in enumerate(self.lines):
            cycle_match = patterns.SCF_CYCLE_RE.search(line)
            if cycle_match is not None:
                if block_start is None:
                    block_start = index
                cycle_numbers.append(int(cycle_match.group("cycle")))
                energies.append(_as_float(cycle_match.group("energy")))
                delta_energies.append(_as_float(cycle_match.group("delta")))
                continue

            end_match = patterns.SCF_END_RE.search(line)
            if end_match is None or not cycle_numbers:
                continue

            reason = end_match.group("reason").upper()
            if "CONVERGENCE" in reason:
                status = SCFStatus.CONVERGED
            elif "TOO MANY CYCLES" in reason:
                status = SCFStatus.TOO_MANY_CYCLES
            else:
                status = SCFStatus.UNKNOWN
            final_energy: EnergyRecord | None = None
            end_energy_match = patterns.SCF_END_ENERGY_RE.search(line)
            if end_energy_match is not None:
                final_energy = EnergyRecord(
                    value=_as_float(end_energy_match.group("energy")),
                    unit="Ha",
                    kind=EnergyKind.DFT,
                    metadata={
                        "source_marker": "SCF ENDED",
                        "line_index": index,
                        "reported_cycles": int(end_energy_match.group("cycles")),
                    },
                )
            results.append(
                self._make_scf_result(
                    status=status,
                    cycle_numbers=cycle_numbers,
                    energies=energies,
                    delta_energies=delta_energies,
                    start_line=block_start,
                    end_line=index,
                    end_reason=end_match.group("reason").strip(),
                    final_energy=final_energy,
                )
            )
            cycle_numbers = []
            energies = []
            delta_energies = []
            block_start = None

        if cycle_numbers:
            results.append(
                self._make_scf_result(
                    status=SCFStatus.INCOMPLETE,
                    cycle_numbers=cycle_numbers,
                    energies=energies,
                    delta_energies=delta_energies,
                    start_line=block_start,
                    end_line=len(self.lines) - 1,
                    end_reason="missing SCF end marker",
                    final_energy=None,
                )
            )
        return tuple(results)

    def optimizations(self) -> tuple[OptimizationResult, ...]:
        """Return geometry-optimization runs in source order.

        Returns
        -------
        tuple of OptimizationResult
            Each run ends at one ``OPT END`` marker. A trailing optimization
            without that marker is retained as ``INCOMPLETE``.
        """
        results: list[OptimizationResult] = []
        points: list[tuple[int, int]] = []

        for index, line in enumerate(self.lines):
            point_match = patterns.OPT_POINT_RE.search(line)
            if point_match is not None:
                points.append((index, int(point_match.group("point"))))
                continue

            end_match = patterns.OPT_END_RE.search(line)
            if end_match is None or not points:
                continue

            reason = end_match.group("reason").upper()
            if "CONVERGED" in reason:
                status = OptimizationStatus.CONVERGED
            elif "FAILED" in reason:
                status = OptimizationStatus.FAILED
            else:
                status = OptimizationStatus.UNKNOWN
            results.append(
                self._make_optimization_result(
                    points=points,
                    end_line=index,
                    status=status,
                    final_energy=EnergyRecord(
                        value=_as_float(end_match.group("energy")),
                        unit="Ha",
                        kind=EnergyKind.TOTAL,
                        metadata={
                            "source_marker": "OPT END",
                            "line_index": index,
                        },
                    ),
                    metadata={
                        "reported_points": int(end_match.group("points")),
                        "end_reason": end_match.group("reason").strip(),
                    },
                )
            )
            points = []

        if points:
            results.append(
                self._make_optimization_result(
                    points=points,
                    end_line=len(self.lines) - 1,
                    status=OptimizationStatus.INCOMPLETE,
                    final_energy=None,
                    metadata={"end_reason": "missing OPT END marker"},
                )
            )
        return tuple(results)

    def _make_scf_result(
        self,
        *,
        status: SCFStatus,
        cycle_numbers: list[int],
        energies: list[float],
        delta_energies: list[float],
        start_line: int | None,
        end_line: int,
        end_reason: str,
        final_energy: EnergyRecord | None,
    ) -> SCFResult:
        """Build one normalized SCF result from parsed source values."""
        if final_energy is None and energies:
            final_energy = EnergyRecord(
                value=energies[-1],
                unit="Ha",
                kind=EnergyKind.DFT,
                metadata={"source_marker": "CYC ETOT(AU)"},
            )
        return SCFResult(
            status=status,
            cycle_numbers=np.asarray(cycle_numbers, dtype=np.int64),
            energies=np.asarray(energies, dtype=np.float64),
            delta_energies=np.asarray(delta_energies, dtype=np.float64),
            energy_unit="Ha",
            final_energy=final_energy,
            metadata={
                "start_line": start_line,
                "end_line": end_line,
                "end_reason": end_reason,
            },
        )

    def _make_optimization_result(
        self,
        *,
        points: list[tuple[int, int]],
        end_line: int,
        status: OptimizationStatus,
        final_energy: EnergyRecord | None,
        metadata: dict[str, object],
    ) -> OptimizationResult:
        """Build one optimization result from point and end markers."""
        steps: list[OptimizationStep] = []
        for offset, (start_line, point_index) in enumerate(points):
            next_line = points[offset + 1][0] if offset + 1 < len(points) else end_line
            steps.append(
                self._parse_optimization_step(
                    point_index=point_index,
                    start_line=start_line,
                    end_line=next_line,
                )
            )
        return OptimizationResult(
            status=status,
            steps=tuple(steps),
            final_energy=final_energy,
            metadata={
                "start_line": points[0][0],
                "end_line": end_line,
                **metadata,
            },
        )

    def _parse_optimization_step(
        self,
        *,
        point_index: int,
        start_line: int,
        end_line: int,
    ) -> OptimizationStep:
        """Parse one optimization point without applying convergence policy."""
        energy: EnergyRecord | None = None
        delta_energy: float | None = None
        metrics: dict[str, float | None] = {
            "max_gradient": None,
            "rms_gradient": None,
            "max_displacement": None,
            "rms_displacement": None,
        }
        metric_patterns = (
            ("max_gradient", patterns.MAX_GRADIENT_RE),
            ("rms_gradient", patterns.RMS_GRADIENT_RE),
            ("max_displacement", patterns.MAX_DISPLACEMENT_RE),
            ("rms_displacement", patterns.RMS_DISPLACEMENT_RE),
        )

        for index in range(start_line, min(end_line + 1, len(self.lines))):
            line = self.lines[index]
            energy_match = patterns.TOTAL_DFT_ENERGY_RE.search(line)
            if energy_match is not None:
                energy = EnergyRecord(
                    value=_as_float(energy_match.group("energy")),
                    unit="Ha",
                    kind=EnergyKind.DFT,
                    metadata={
                        "source_marker": "TOTAL ENERGY(DFT)(AU)",
                        "line_index": index,
                    },
                )
                delta_match = patterns.TOTAL_DFT_DELTA_RE.search(line)
                if delta_match is not None:
                    delta_energy = _as_float(delta_match.group("delta"))
            for name, regex in metric_patterns:
                match = regex.search(line)
                if match is not None:
                    metrics[name] = _as_float(match.group("value"))

        return OptimizationStep(
            index=point_index,
            energy=energy,
            delta_energy=delta_energy,
            max_gradient=metrics["max_gradient"],
            rms_gradient=metrics["rms_gradient"],
            max_displacement=metrics["max_displacement"],
            rms_displacement=metrics["rms_displacement"],
            metadata={"start_line": start_line, "end_line": end_line},
        )
