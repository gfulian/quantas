# -*- coding: utf-8 -*-

"""VASP structure--energy adaptation for static Energy EOS workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from quantas.core.geometry.symmetry import analyze_symmetry, reduce_to_primitive_cell
from quantas.interfaces.vasp.output import VaspOutputParser
from quantas.models.computation import (
    EnergyKind,
    EnergyRecord,
    RunTerminationStatus,
    SourceProvenance,
    StructureEnergyPoint,
    StructureEnergySeries,
)


_SIGNATURE_PARAMETERS = (
    "PREC",
    "ENCUT",
    "ICHARG",
    "LREAL",
    "ISYM",
    "ISPIN",
    "NELECT",
    "GGA",
    "METAGGA",
    "LHFCALC",
    "AEXX",
    "HFSCREEN",
    "LDAU",
    "LDAUTYPE",
    "IVDW",
    "LUSE_VDW",
    "LASPH",
    "LSORBIT",
    "LNONCOLLINEAR",
    "LMAXMIX",
)
_SIGNATURE_VECTOR_PARAMETERS = (
    "LDAUL",
    "LDAUU",
    "LDAUJ",
    "MAGMOM",
)
_SIGMA_INACTIVE_ISMEAR = {-5, -4}


@dataclass(frozen=True, slots=True)
class VaspEnergyVolumeParseResult:
    """One primitive-normalized VASP state for an Energy EOS series.

    Parameters
    ----------
    series : StructureEnergySeries
        One-point backend-neutral structure--energy series.
    run_kind : str
        VASP run classification.  The b13 adapter currently accepts only a
        single resolved ionic state and reports ``"single_state"``.
    energy_signature : tuple of str
        Canonical signature of the selected VASP energy quantity and the
        electronic settings that define its energy surface.
    source_to_primitive : ndarray
        Row-basis transformation from the VASP source cell to the normalized
        primitive representation.  It is the identity when the source is
        already primitive and its backend basis is preserved.
    source_repetitions : int
        Number of primitive cells represented by the VASP source cell.
    primitive_to_crystallographic : ndarray or None, optional
        Explicit primitive-to-crystallographic transformation.  VASP does not
        print a CRYSTAL-style transformation, so the current adapter leaves
        this as ``None`` and lets the shared EOS layer use spglib.
    space_group_number : int or None, optional
        Space group determined from the normalized primitive structure.
    metadata : dict, optional
        Additional VASP normalization and energy-selection provenance.
    """

    series: StructureEnergySeries
    run_kind: str
    energy_signature: tuple[str, ...]
    source_to_primitive: np.ndarray
    source_repetitions: int
    primitive_to_crystallographic: np.ndarray | None = None
    space_group_number: int | None = None
    metadata: dict[str, Any] | None = None


class VaspEnergyVolumeReader:
    """Adapt one completed VASP run to a primitive Energy EOS observation.

    The generic VASP output parser remains responsible for reconstructing the
    run.  This adapter imposes the first EOS-specific scientific policy:

    * one source directory contributes exactly one ionic state;
    * ``e_0_energy`` / ``energy(sigma->0)`` is the static energy quantity;
    * source cells are reduced with the shared spglib primitive-cell helper;
    * source energies are divided by the corresponding cell multiplicity;
    * the selected energy semantics and relevant electronic settings are
      recorded in an explicit compatibility signature.

    Parameters
    ----------
    source : str or pathlib.Path
        VASP calculation directory, ``vasprun.xml``, or sibling ``OUTCAR``.
    symprec : float, optional
        Cartesian spglib tolerance in angstrom used for primitive reduction.

    Raises
    ------
    ValueError
        If the VASP run is incomplete, contains more than one ionic state, uses
        an energy integration mode outside the current static-EOS policy, or
        cannot be normalized to a primitive structure consistently.
    ImportError
        If spglib is unavailable.
    """

    def __init__(self, source: str | Path, *, symprec: float = 1.0e-5) -> None:
        self.output = VaspOutputParser(source)
        self.symprec = float(symprec)

    def read(self) -> VaspEnergyVolumeParseResult:
        """Return the primitive-normalized VASP Energy EOS state.

        Returns
        -------
        VaspEnergyVolumeParseResult
            One-point series with VASP and primitive-normalization provenance.

        Raises
        ------
        ValueError
            If the run cannot define one unambiguous static E(V) observation.
        ImportError
            If spglib is unavailable.
        """
        termination = self.output.termination()
        if termination.status is not RunTerminationStatus.NORMAL:
            raise ValueError("VASP Energy EOS source did not terminate normally")

        steps = self.output.ionic_steps()
        if len(steps) != 1:
            raise ValueError(
                "VASP Energy EOS source must contain exactly one ionic state; "
                f"found {len(steps)}. Use a dedicated single-state calculation "
                "for every E(V) point rather than mixing an optimization history "
                "into the EOS dataset."
            )
        step = steps[0]
        signature = _energy_signature(self.output)
        reduction = reduce_to_primitive_cell(
            step.structure,
            symprec=self.symprec,
            no_idealize=True,
        )
        symmetry = analyze_symmetry(reduction.structure, symprec=self.symprec)

        source_energy = step.energies.sigma_zero_energy
        primitive_energy = EnergyRecord(
            value=source_energy.value / reduction.repetitions,
            unit=source_energy.unit,
            kind=EnergyKind.DFT,
            metadata={
                **source_energy.metadata,
                "energy_selection": "sigma_zero_energy",
                "vasp_tag": "e_0_energy",
                "source_energy_eV": float(source_energy.value),
                "primitive_repetitions": int(reduction.repetitions),
                "normalization": "per_primitive_cell",
                "energy_signature": signature,
            },
        )
        source = self.output.source
        point = StructureEnergyPoint(
            structure=reduction.structure,
            energy=primitive_energy,
            provenance=SourceProvenance(
                interface="vasp",
                source=str(source.directory),
                record_index=0,
                metadata={
                    "vasprun_xml": str(source.vasprun_xml),
                    "outcar": None if source.outcar is None else str(source.outcar),
                    "vasp_version": self.output.version,
                    "source_ionic_step": int(step.index),
                },
            ),
            metadata={
                "run_kind": "single_state",
                "reported_volume_angstrom3": float(reduction.structure.volume),
                "source_volume_angstrom3": float(reduction.source_volume),
                "source_atoms": int(reduction.source_atoms),
                "primitive_repetitions": int(reduction.repetitions),
                "energy_signature": signature,
            },
        )
        series = StructureEnergySeries(
            points=(point,),
            reference_index=0,
            metadata={
                "interface": "vasp",
                "run_kind": "single_state",
                "energy_signature": signature,
                "source": str(source.directory),
                "vasp_version": self.output.version,
            },
        )
        return VaspEnergyVolumeParseResult(
            series=series,
            run_kind="single_state",
            energy_signature=signature,
            source_to_primitive=reduction.source_to_primitive.copy(),
            source_repetitions=reduction.repetitions,
            space_group_number=symmetry.space_group_number,
            metadata={
                "energy_selection": "sigma_zero_energy",
                "symprec": self.symprec,
                "pseudopotentials": self.output.document.pseudopotential_labels(),
            },
        )


def read_vasp_energy_volume(
    source: str | Path,
    *,
    symprec: float = 1.0e-5,
) -> VaspEnergyVolumeParseResult:
    """Read one VASP run into one primitive-normalized Energy EOS state.

    Parameters
    ----------
    source : str or pathlib.Path
        VASP calculation directory, ``vasprun.xml``, or sibling ``OUTCAR``.
    symprec : float, optional
        Cartesian spglib tolerance in angstrom.

    Returns
    -------
    VaspEnergyVolumeParseResult
        One-point Energy EOS parse result.
    """
    return VaspEnergyVolumeReader(source, symprec=symprec).read()


def _energy_signature(parser: VaspOutputParser) -> tuple[str, ...]:
    """Return the canonical VASP Energy EOS compatibility signature."""
    document = parser.document
    raw_ismear = document.parameter("ISMEAR")
    if raw_ismear is None:
        raise ValueError("VASP Energy EOS source does not report ISMEAR")
    ismear = int(raw_ismear)
    if ismear == -3:
        raise ValueError(
            "VASP ISMEAR=-3 performs multiple smearing evaluations and is not "
            "an unambiguous single Energy EOS energy surface"
        )
    if ismear == -2:
        raise ValueError(
            "VASP ISMEAR=-2 uses fixed occupancies and is outside the current "
            "ground-state Energy EOS policy"
        )

    version = parser.version
    if not version:
        raise ValueError("VASP Energy EOS source does not report the VASP version")
    values = [
        "quantity=e_0_energy",
        f"VASP_VERSION={version}",
        f"ISMEAR={ismear}",
    ]
    if ismear not in _SIGMA_INACTIVE_ISMEAR:
        sigma = document.parameter("SIGMA")
        if sigma is None:
            raise ValueError(
                f"VASP Energy EOS source with ISMEAR={ismear} does not report SIGMA"
            )
        values.append(f"SIGMA={_format_signature_value(sigma)}")

    for name in _SIGNATURE_PARAMETERS:
        value = document.parameter(name)
        if value is not None:
            values.append(f"{name}={_format_signature_value(value)}")
    for name in _SIGNATURE_VECTOR_PARAMETERS:
        vector = document.parameter_vector(name)
        if vector is not None:
            rendered = ",".join(_format_signature_value(value) for value in vector)
            values.append(f"{name}={rendered}")
    kpoints = document.kpoint_signature()
    if not kpoints:
        raise ValueError(
            "VASP Energy EOS source does not report Brillouin-zone sampling"
        )
    values.extend(f"KPOINTS:{item}" for item in kpoints)
    labels = document.pseudopotential_labels()
    if not labels:
        raise ValueError(
            "VASP Energy EOS source does not report pseudopotential labels"
        )
    values.append("POTENTIALS=" + "|".join(labels))
    return tuple(values)


def _format_signature_value(value: Any) -> str:
    """Return a stable scalar representation for compatibility signatures."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".15g")
    return str(value).strip()


__all__ = [
    "VaspEnergyVolumeParseResult",
    "VaspEnergyVolumeReader",
    "read_vasp_energy_volume",
]
