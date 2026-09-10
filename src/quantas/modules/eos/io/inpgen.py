# -*- coding: utf-8 -*-

"""Frontend-neutral Energy EOS input generation from external-code outputs."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any, ClassVar

import numpy as np

from quantas.core.events import Event, EventLevel, NullObserver, Observer
from quantas.core.geometry.cells import lattice_parameters
from quantas.interfaces.crystal.energy_volume import (
    CrystalEnergyVolumeParseResult,
    read_crystal_energy_volume,
)
from quantas.models.computation import StructureEnergyPoint, StructureEnergySeries


_DUPLICATE_VOLUME_ATOL = 5.0e-6


class EOSEnergyInputCreator:
    """Create Quantas Energy EOS text inputs from backend output files.

    One source file may contribute one or several structure--energy states. In
    particular, the CRYSTAL interface accepts static calculations, ordinary
    geometry optimizations, and native CRYSTAL ``EOS`` runs. When ``is_list``
    is true, every listed source is parsed independently and the resulting
    points are flattened into one compatible Energy EOS series.

    Parameters
    ----------
    interface : str, optional
        Backend interface identifier. ``"crystal"`` is currently supported.
    observer : Observer or None, optional
        Frontend-neutral observer receiving input-generation events.
    """

    interface_filter: ClassVar[dict[str, Any]] = {
        "crystal": read_crystal_energy_volume,
    }

    def __init__(
        self,
        interface: str = "crystal",
        observer: Observer | None = None,
    ) -> None:
        """Initialize one Energy EOS input creator."""
        self.interface = str(interface).strip().lower()
        self.observer = observer if observer is not None else NullObserver()
        self.files: list[Path] = []
        self.results: list[CrystalEnergyVolumeParseResult] = []
        self.series: StructureEnergySeries | None = None

    def emit(
        self,
        message: str,
        *,
        level: EventLevel = EventLevel.INFO,
        data: dict[str, Any] | None = None,
    ) -> None:
        """Emit one frontend-neutral input-generation event.

        Parameters
        ----------
        message : str
            Human-readable event text.
        level : EventLevel, optional
            Event severity.
        data : dict or None, optional
            Structured event payload.
        """
        self.observer(Event(message=message, level=level, data=data or {}))

    def read(
        self,
        sources: str | Path | Sequence[str | Path],
        *,
        is_list: bool = False,
    ) -> StructureEnergySeries:
        """Read and merge one or more backend output files.

        Parameters
        ----------
        sources : str, Path, or sequence of path-like
            One backend output, a direct sequence of outputs, or a list-file
            path when ``is_list`` is true.
        is_list : bool, optional
            Interpret a scalar source as a text file containing one output path
            per non-comment line.

        Returns
        -------
        StructureEnergySeries
            Compatible structure--energy states sorted by increasing volume.

        Raises
        ------
        ValueError
            If the interface is unsupported, a source is invalid, states use
            incompatible chemistry/energy semantics, or duplicate volumes are
            present.
        OSError
            If an input file cannot be read.
        """
        if self.interface not in self.interface_filter:
            raise ValueError(f"Unsupported EOS input interface: {self.interface}")
        files = self._normalize_sources(sources, is_list=is_list)
        if not files:
            raise ValueError("No Energy EOS source files were provided")

        reader = self.interface_filter[self.interface]
        self.emit(
            f"Reading {len(files)} Energy EOS source file(s) with the {self.interface} interface",
            data={
                "kind": "eos_input_sources",
                "interface": self.interface,
                "source_count": len(files),
            },
        )

        results: list[CrystalEnergyVolumeParseResult] = []
        points: list[StructureEnergyPoint] = []
        for source_index, path in enumerate(files):
            result = reader(path)
            results.append(result)
            for point in result.series.points:
                point.metadata.setdefault("source_file_index", source_index)
                points.append(point)
            self.emit(
                f"Parsed Energy EOS source {path}",
                level=EventLevel.DEBUG,
                data={
                    "kind": "eos_input_source_parsed",
                    "source": str(path),
                    "source_index": source_index,
                    "state_count": result.series.npoints,
                    "run_kind": result.run_kind,
                },
            )

        self._validate_structure_compatibility(points)
        self._validate_correction_signatures(results)
        ordered = sorted(points, key=_point_volume)
        self._validate_duplicate_volumes(ordered)
        series = StructureEnergySeries(
            points=tuple(ordered),
            reference_index=min(
                range(len(ordered)),
                key=lambda index: ordered[index].energy.value,
            ),
            metadata={
                "interface": self.interface,
                "source_files": tuple(str(path) for path in files),
                "source_run_kinds": tuple(result.run_kind for result in results),
                "corrections": results[0].corrections,
            },
        )
        self.files = files
        self.results = results
        self.series = series
        self.emit(
            "Energy EOS source parsing completed",
            level=EventLevel.RESULT,
            data={
                "kind": "eos_input_series",
                "source_count": len(files),
                "state_count": series.npoints,
                "volume_min_angstrom3": float(series.volumes.min()),
                "volume_max_angstrom3": float(series.volumes.max()),
            },
        )
        return series

    def write(
        self,
        destination: str | Path,
        *,
        jobname: str = "Quantas Energy EOS input",
    ) -> Path:
        """Write the parsed series in the canonical Quantas EOS text format.

        Parameters
        ----------
        destination : str or Path
            Output text path.
        jobname : str, optional
            Human-readable dataset title.

        Returns
        -------
        Path
            Written path.

        Raises
        ------
        RuntimeError
            If :meth:`read` has not completed successfully.
        OSError
            If the destination cannot be written.
        """
        if self.series is None:
            raise RuntimeError("Energy EOS sources must be read before writing")
        path = Path(destination)
        rows = [
            _format_point(point)
            for point in self.series.points
        ]
        provenance = f"Generated by Quantas from {self.interface.upper()} structure-energy output"
        comments = [
            f"COMMENT source {index + 1}: {source}"
            for index, source in enumerate(self.series.metadata.get("source_files", ()))
        ]
        text = "\n".join(
            [
                f"JOB {jobname}",
                f"PROVENANCE {provenance}",
                *comments,
                (
                    "UNITS V=angstrom^3 A=angstrom B=angstrom C=angstrom "
                    "ALPHA=degree BETA=degree GAMMA=degree "
                    f"E={self.series.energy_unit}"
                ),
                "FORMAT V A B C ALPHA BETA GAMMA E",
                "DATA",
                *rows,
                "",
            ]
        )
        path.write_text(text, encoding="utf-8")
        return path

    @staticmethod
    def _normalize_sources(
        sources: str | Path | Sequence[str | Path],
        *,
        is_list: bool,
    ) -> list[Path]:
        """Return explicit unique source paths with list entries resolved locally."""
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
                raise ValueError(f"Energy EOS source file does not exist: {path}")
            seen.add(normalized)
            resolved.append(normalized)
        return resolved

    @staticmethod
    def _validate_structure_compatibility(
        points: Sequence[StructureEnergyPoint],
    ) -> None:
        """Reject source files that do not share one cell composition."""
        reference = points[0].structure
        numbers = np.sort(reference.atomic_numbers)
        for point in points[1:]:
            structure = point.structure
            if structure.natoms != reference.natoms:
                raise ValueError(
                    "Energy EOS source files use inconsistent atom counts"
                )
            if not np.array_equal(np.sort(structure.atomic_numbers), numbers):
                raise ValueError(
                    "Energy EOS source files use inconsistent chemical compositions"
                )

    @staticmethod
    def _validate_correction_signatures(
        results: Sequence[CrystalEnergyVolumeParseResult],
    ) -> None:
        """Reject source files that describe different total-energy surfaces."""
        signatures = {tuple(result.corrections) for result in results}
        if len(signatures) == 1:
            return
        rendered = ", ".join(
            "+".join(signature) if signature else "DFT"
            for signature in sorted(signatures)
        )
        raise ValueError(
            "Energy EOS source files use incompatible energy corrections: "
            f"{rendered}"
        )

    @staticmethod
    def _validate_duplicate_volumes(points: Sequence[StructureEnergyPoint]) -> None:
        """Reject numerically duplicate volumes without merging nearby states."""
        for left, right in zip(points, points[1:]):
            if np.isclose(
                _point_volume(left),
                _point_volume(right),
                rtol=1.0e-12,
                atol=_DUPLICATE_VOLUME_ATOL,
            ):
                raise ValueError(
                    "Energy EOS source series contains duplicate volumes near "
                    f"{_point_volume(left):.8g} Å^3"
                )


def create_eos_energy_input(
    sources: str | Path | Sequence[str | Path],
    destination: str | Path,
    *,
    interface: str = "crystal",
    is_list: bool = False,
    jobname: str = "Quantas Energy EOS input",
    observer: Observer | None = None,
) -> Path:
    """Create one Quantas Energy EOS input from external-code outputs.

    Parameters
    ----------
    sources : str, Path, or sequence of path-like
        Backend output source(s), or a list-file path when ``is_list`` is true.
    destination : str or Path
        Destination EOS text file.
    interface : str, optional
        Backend interface identifier. Currently ``"crystal"``.
    is_list : bool, optional
        Interpret a scalar source as a text file listing output paths.
    jobname : str, optional
        Human-readable dataset title.
    observer : Observer or None, optional
        Frontend-neutral input-generation observer.

    Returns
    -------
    Path
        Written EOS input path.
    """
    creator = EOSEnergyInputCreator(interface=interface, observer=observer)
    creator.read(sources, is_list=is_list)
    return creator.write(destination, jobname=jobname)


def _point_volume(point: StructureEnergyPoint) -> float:
    """Return the backend-reported point volume when available."""
    value = point.metadata.get("reported_volume_angstrom3")
    return point.volume if value is None else float(value)


def _cell_parameters(point: StructureEnergyPoint) -> np.ndarray:
    """Return source cell parameters or reconstruct them from the lattice."""
    values = point.structure.metadata.get("cell_parameters")
    if values is None:
        return lattice_parameters(point.structure.lattice)
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (6,):
        raise ValueError("structure cell_parameters metadata must contain six values")
    return array


def _format_point(point: StructureEnergyPoint) -> str:
    """Serialize one Energy EOS state without display-precision rounding."""
    parameters = _cell_parameters(point)
    values = [
        _point_volume(point),
        *parameters.tolist(),
        float(point.energy.value),
    ]
    return " ".join(format(float(value), ".15g") for value in values)


__all__ = ["EOSEnergyInputCreator", "create_eos_energy_input"]
